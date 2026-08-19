#=
################################################################################
Effect of bodysize specification rules
################################################################################


################################################################################
=#

# --- 1. General Set-up ---

using CSV
using DataFrames
using DifferentialEquations
using Distributions
using EcologicalNetworksDynamics
using Extinctions
using FoodWebTools
using JLD2
using pfim
using ProgressMeter
using SpeciesInteractionNetworks
using Statistics

include("libs/internals.jl")

import Random
Random.seed!(66)

# --- storage ---

network_archive = Dict()

species_metadata_store = DataFrame()

adjacency_store = DataFrame(
    t=Int[],
    net_id=Int[],
    community=String[],
    net_type=String[],
    stage=String[],
    S=Int[],
    C=Float64[],
    adjacency=Any[]
)

summary_rows = Dict[]

# --- data ---
communities = Dict("turkiye" => Dict(
    "traits" => CSV.read("data/turkiye_community.csv", DataFrame),
    "plankton" => CSV.read("data/turkiye_plankton.csv", DataFrame)
))
feeding_rules = CSV.read("data/feeding_rules.csv", DataFrame)

# --- global params ---
n_networks = 30  # number of runs per network
survival_threshold = 1e-12
C_min = 0.05
C_max = 0.15

# Define the t-values we want to test
t_values = [5000]

# --- Distributions ---
C_dist = truncated(Normal(0.15, 0.05), C_min, C_max)
global_dist = LogNormal(log(30), 1.5)

size_bounds = Dict(
    "primary" => (0.01, 0.1),
    "tiny" => (0.1, 10.0),
    "small" => (10.0, 50.0),
    "medium" => (50.0, 100.0),
    "large" => (100.0, 300.0),
    "very_large" => (300.0, 500.0),
    "gigantic" => (500.0, Inf)
)

# --- 2. Build and realise networks ---

for t in t_values

    # 2. Inner Loop: Process networks
    for i in 1:n_networks

        C_targ = rand(C_dist)

        for (community_name, community_data) in communities

            traits = deepcopy(community_data["traits"])
            plankton = deepcopy(community_data["plankton"])

            println(">>> Processing run $i of $n_networks (t = $t) for $community_name...")

            # --- 1. Sample parameters & Body sizes ---
            y = collect(String, traits.size)

            bodysize = [
                begin
                    lo, hi = size_bounds[s]
                    rand(truncated(global_dist, lo, hi))
                end
                for s in y
            ]

            traits[!, :bodymass] = bodysize
            traits.biomass = fill(missing, nrow(traits))
            df = vcat(traits, plankton)

            # --- 2. Biomass estimates ---
            # create initial biomass from uniform rand dist
            B_init = rand(nrow(df))

            # --- 3. Base Networks Generation (Creation) ---
            mass_rule = (res, con) -> con >= 0.5 * res ? 1 : 0

            pfim_meta = PFIM(df, feeding_rules; return_type=:matrix)

            # Create the Niche-Downsampled Network from pfim_cont
            pfim_niche_down = Foodweb(Matrix(downsample(pfim_meta, :niche; sigma_scale=1.0, target_co=C_targ, max_iter=300)))

            niche_fw = Foodweb(:niche; S=size(pfim_meta, 1), C=C_targ)

            # For ATN we will vary thresholding to find best fit to Co
            prods = map(==("primary"), string.(df.tiering))

            # build range
            atn_range = atn_build_range(df, prods)
            # find best fit
            best = argmin(abs.(getfield.(atn_range, :connectance) .- C_targ))
            # get Foofweb object
            atn_fw = Foodweb(atn_range[best].adjacency)

            # Consolidate all initial networks into one dictionary
            initial_networks = Dict(
                "niche" => niche_fw,
                "down_niche" => pfim_niche_down,
                "atn" => atn_fw,
                "metaweb" => Foodweb(pfim_meta),
                "down_niche_empiricalbm" => pfim_niche_down,
                "atn_empiricalbm" => atn_fw,
                "metaweb_empiricalbm" => Foodweb(pfim_meta)
            )

            # Dictionary to keep initial S and C values for summary reference
            creation_metrics = Dict()

            for (net_name, fw) in initial_networks

                # summary metrics
                S_init = size(fw.A, 1)
                C_init = sum(fw.A) / (S_init^2)
                creation_metrics[net_name] = (S=S_init, C=C_init)

                # Extract Species Metadata at creation
                # Create a temporary container for this specific record
                temp_metadata = DataFrame()

                # Let the function record metadata into our temporary DataFrame
                # we can set default params
                params = default_model(fw)

                record_species_stage!(temp_metadata, i, net_name, "creation", fw.A, params.M, params, nothing, B_init)

                # Inject the current t-value into the temporary DataFrame
                temp_metadata[!, :t_val] .= t
                temp_metadata[!, :community] .= community_name

                # Append the completed DataFrame to our global store
                append!(species_metadata_store, temp_metadata, cols=:union)

                # Save the adjacency matrix as well
                push!(adjacency_store, (
                    t=t,
                    net_id=i,
                    community=community_name,
                    net_type=net_name,
                    stage="creation",
                    S=size(fw.A, 1),
                    C=sum(fw.A) / size(fw.A, 1)^2,
                    adjacency=copy(fw.A)
                ))

            end

            # --- 4. Burn-In & Realisation ---
            realised_networks = Dict()

            for (net_name, fw) in initial_networks

                # for everything BUT the niche we do an empirical bodymass run
                if net_name ∈ ["metaweb_empiricalbm", "down_niche_empiricalbm", "atn_empiricalbm"]
                    realised = realise_network(
                        fw;
                        bodymasses=df.bodymass,
                        t=t,
                        threshold=survival_threshold,
                        B0=B_init
                    )
                else
                    # and then we do the normal (internal bodysize specification)
                    realised = realise_network(
                        fw;
                        t=t,
                        threshold=survival_threshold,
                        B0=B_init
                    )

                end

                if realised !== nothing
                    realised_networks[net_name] = realised
                end
            end

            # Skip this replicate iteration if no networks successfully survived burn-in
            if isempty(realised_networks)
                @warn "Iteration $i (t = $t): no realised networks. Skipping."
                continue
            end

            # Record realised networks and associated metadata

            for (net_name, realised) in realised_networks

                A_realised = realised.A
                params = realised.params
                final_biomasses = realised.biomasses
                survivors = realised.survivors

                S_realised = length(survivors)
                C_realised = sum(A_realised) / (S_realised^2)

                # Save realised adjacency matrix

                push!(adjacency_store, (
                    t=t,
                    net_id=i,
                    community=community_name,
                    net_type=net_name,
                    stage="burnin",
                    S=S_realised,
                    C=C_realised,
                    adjacency=copy(A_realised)
                ))

                # Record species metadata after burn-in

                temp_metadata = DataFrame()

                record_species_stage!(
                    temp_metadata,
                    i,
                    net_name,
                    "post_burn_in",
                    A_realised,
                    params.M,
                    params,
                    survivors,
                    final_biomasses
                )

                temp_metadata.t_val .= t
                temp_metadata.community .= community_name

                append!(species_metadata_store, temp_metadata, cols=:union)

                # Network summary

                row = Dict(
                    :t_val => t,
                    :net_id => i,
                    :community => community_name,
                    :net_type => net_name,
                    :C_target => C_targ,
                    :S_creation => creation_metrics[net_name].S,
                    :C_creation => creation_metrics[net_name].C,
                    :C_realised => C_realised
                )

                push!(summary_rows, row)

                # Archive realised network

                network_archive[(community_name, i, net_name)] = Dict("generation" => Dict(
                        "community" => community_name,
                        "burnin_time" => t,
                        "C_target" => C_targ,
                        "bodymass" => copy(df.bodymass),
                        "B0" => copy(B_init)
                    ), "creation" => Dict(
                        "A" => copy(initial_networks[net_name].A),
                        "S" => creation_metrics[net_name].S,
                        "C" => creation_metrics[net_name].C
                    ), "realised" => Dict(
                        "A" => copy(A_realised),
                        "params" => params,
                        "biomasses" => final_biomasses,
                        "survivors" => survivors,
                        "S" => S_realised,
                        "C" => C_realised
                    )
                )
            end
        end
    end
end


# --- 3. Collate and export relevant data ---

archive = Dict(
    # Global metadata
    "metadata" => Dict(
        "seed" => 66,
        "n_networks" => n_networks,
        "survival_threshold" => survival_threshold,
        "burnin_times" => t_values
    ),
    # Realised network collection
    "networks" => network_archive,
    # Species metadata
    "species_metadata" => species_metadata_store,
    # Stored adjacency matrices
    "adjacency_store" => adjacency_store,
    # Network summaries
    "summary" => DataFrame(summary_rows)
)

JLD2.save_object(
    "outputs/S1_bodysize/network_archive.jld2",
    archive
)

CSV.write("outputs/S1_bodysize/species_metadata.csv", archive["species_metadata"])
JLD2.save_object("outputs/S1_bodysize/adjacency_matrices.jld2", archive["adjacency_store"])

# get topology

include("libs/topology.jl")
include("libs/internals.jl")

# import adjacency matrices

networks = load_object("outputs/S1_bodysize/adjacency_matrices.jld2")

# data frame

rows = Dict[]

for i in 1:nrow(networks)

    N = build_network(Matrix(networks.adjacency[i]))

    # push network summary
    d = _network_summary(N)
    d[:net_id] = networks.net_id[i]
    d[:net_type] = networks.net_type[i]
    d[:community] = networks.community[i]
    d[:stage] = networks.stage[i]

    # send to df
    push!(rows, d)

end

topology_df = DataFrame(rows)

CSV.write("outputs/S1_bodysize/topology.csv", topology_df)

# get gen and vul of each species from each network

deg_rows = DataFrame(
    spp_id=Any[],
    gen=Any[],
    vul=Any[],
    can=Any[],
    net_id=Any[],
    net_type=Any[],
    community=Any[],
    stage=Any[],
    S4_consumer=Any[],
    S4_resource=Any[],
    S5_consumer=Any[],
    S5_resource=Any[])

for i in 1:nrow(networks)

    N = build_network(Matrix(networks.adjacency[i]))

    gen = SpeciesInteractionNetworks.generality(N)
    vul = SpeciesInteractionNetworks.vulnerability(N)

    # also get motif membership
    spp = SpeciesInteractionNetworks.species(N)
    S4 = findmotif(motifs(Unipartite, 3)[4], N)
    S5 = findmotif(motifs(Unipartite, 3)[5], N)

    S4_counts = [count(t -> t[i] == s, S4)
              for s in spp, i in 1:3]

    
    S5_counts = [count(t -> t[i] == s, S5)
              for s in spp, i in 1:3]

    d = DataFrame(
        spp_id=collect(keys(gen)),
        gen=collect(values(gen)),
        vul=collect(values(vul)),
        can=string.(diag(Matrix(N.edges.edges))),
        net_id=fill(networks.net_id[i], length(gen)),
        net_type=fill(networks.net_type[i], length(gen)),
        stage=fill(networks.stage[i], length(gen)),
        community=fill(networks.community[i], length(gen)),
        S4_consumer=S4_counts[:, 1] .+ S4_counts[:, 3],
        S4_resource=S4_counts[:, 2],
        S5_consumer=S5_counts[:, 1],
        S5_resource=S5_counts[:, 2] .+ S5_counts[:, 3]
    )

    # send to df
    append!(deg_rows, d)

end

CSV.write("outputs/S1_bodysize/spp_degrees.csv", deg_rows)