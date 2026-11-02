################################################################################
# BUILD REALISED PRE AND POST NETWORKS
################################################################################

# Purpose
# -------
# Generate paired pre- and post-community ecological networks for each run.
# Post communities are restricted to species present in the corresponding
# pre community. Shared species retain the same body masses and initial
# biomasses in both states. Plankton data are shared between states.
#
# Outputs
# -------
# outputs/pre/realised_network_archive.jld2
# outputs/post/realised_network_archive.jld2
#
# Each archive contains realised networks, creation networks, species
# metadata, adjacency matrices, generation parameters and summaries.

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

function main()
    Random.seed!(66)

    # --- Storage ---

    states = ["pre", "post"]

    network_archive = Dict(state => Dict() for state in states)

    species_metadata_store = DataFrame()

    adjacency_store = DataFrame(
        t=Int[],
        net_id=Int[],
        community=String[],
        state=String[],
        net_type=String[],
        stage=String[],
        S=Int[],
        C=Float64[],
        adjacency=Any[]
    )

    summary_rows = Dict[]

    # --- Data ---

    communities = Dict(
        "dolomites" => Dict(
            "pre" => CSV.read("data/pre/dolomites_community.csv", DataFrame),
            "post" => CSV.read("data/post/dolomites_community.csv", DataFrame),
            "plankton" => CSV.read("data/plankton/dolomites_plankton.csv", DataFrame)
        ),
        "greenland" => Dict(
            "pre" => CSV.read("data/pre/greenland_community.csv", DataFrame),
            "post" => CSV.read("data/post/greenland_community.csv", DataFrame),
            "plankton" => CSV.read("data/plankton/greenland_plankton.csv", DataFrame)
        ),
        "tibet" => Dict(
            "pre" => CSV.read("data/pre/tibet_community.csv", DataFrame),
            "post" => CSV.read("data/post/tibet_community.csv", DataFrame),
            "plankton" => CSV.read("data/plankton/tibet_plankton.csv", DataFrame)
        ),
        "turkiye" => Dict(
            "pre" => CSV.read("data/pre/turkiye_community.csv", DataFrame),
            "post" => CSV.read("data/post/turkiye_community.csv", DataFrame),
            "plankton" => CSV.read("data/plankton/turkiye_plankton.csv", DataFrame)
        ),
        "meishan" => Dict(
            "pre" => CSV.read("data/pre/meishan_community.csv", DataFrame),
            "post" => CSV.read("data/post/meishan_community.csv", DataFrame),
            "plankton" => CSV.read("data/plankton/meishan_plankton.csv", DataFrame)
        ),
        "kashmir" => Dict(
            "pre" => CSV.read("data/pre/kashmir_community.csv", DataFrame),
            "post" => CSV.read("data/post/kashmir_community.csv", DataFrame),
            "plankton" => CSV.read("data/plankton/kashmir_plankton.csv", DataFrame)
        )
    )

    feeding_rules = CSV.read("data/feeding_rules.csv", DataFrame)

    # --- Global parameters ---

    n_networks = 30                # number of paired runs per community
    survival_threshold = 1e-12
    C_min = 0.05
    C_max = 0.15

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

    ################################################################################
    # 2. Build and realise paired networks
    ################################################################################

    for t in t_values

        for i in 1:n_networks

            # One target connectance shared by pre and post for this replicate.
            C_targ = rand(C_dist)

            for (community_name, community_data) in communities

                println(">>> Preparing run $i/$n_networks, t=$t, community=$community_name")

                # Work on copies so the original loaded data are unchanged.
                pre_traits = deepcopy(community_data["pre"])
                post_traits = deepcopy(community_data["post"])
                plankton = deepcopy(community_data["plankton"])

                # Check the species identifier is available.
                for df in (pre_traits, post_traits, plankton)
                    @assert :species in propertynames(df) "Missing :species column in $community_name data."
                end

                # Remove post species that were not present in the pre community.
                pre_species = Set(String.(pre_traits.species))
                filter!(row -> String(row.species) in pre_species, post_traits)

                # Sample body masses once from the pre community.
                y = String.(pre_traits.size)

                bodysize = [
                    begin
                        lo, hi = size_bounds[s]
                        rand(truncated(global_dist, lo, hi))
                    end
                    for s in y
                ]

                bodymass_lookup = Dict(zip(String.(pre_traits.species), bodysize))

                pre_traits[!, :bodymass] = bodysize
                post_traits[!, :bodymass] = [
                    bodymass_lookup[String(s)] for s in post_traits.species
                ]

                pre_traits[!, :biomass] = fill(missing, nrow(pre_traits))
                post_traits[!, :biomass] = fill(missing, nrow(post_traits))

                # Add the same plankton community to both states.
                pre_df = vcat(pre_traits, deepcopy(plankton); cols=:union)
                post_df = vcat(post_traits, deepcopy(plankton); cols=:union)

                # Sample initial biomasses once, then match by species.
                # This ensures shared species (including plankton) start identically.
                @assert length(unique(String.(pre_df.species))) == nrow(pre_df) "Species identifiers must be unique in the combined pre community and plankton data for $community_name."

                B_init_pre = rand(nrow(pre_df))
                B0_lookup = Dict(zip(String.(pre_df.species), B_init_pre))

                @assert all(
                    s -> haskey(B0_lookup, String(s)),
                    post_df.species
                ) "Some post species are missing from the pre biomass lookup for $community_name."

                B_init_post = [
                    B0_lookup[String(s)] for s in post_df.species
                ]

                # Process pre and post independently, using paired inputs.
                for (state, df, B_init) in (
                    ("pre", pre_df, B_init_pre),
                    ("post", post_df, B_init_post)
                )

                    println("    >>> Building $state networks")

                    # --- 1. Base network generation ---

                    mass_rule = (res, con) -> con >= 0.5 * res ? 1 : 0

                    pfim_meta = PFIM(df, feeding_rules; return_type=:matrix)

                    pfim_power_down = Foodweb(Matrix(
                        downsample(pfim_meta, :powerlaw;
                            y=2.5, target_co=C_targ, max_iter=300)
                    ))

                    pfim_niche_down = Foodweb(Matrix(
                        downsample(pfim_meta, :niche;
                            sigma_scale=1.0, target_co=C_targ, max_iter=300)
                    ))

                    pfim_link_down = Foodweb(Matrix(
                        downsample(pfim_meta, :degree;
                            target_co=C_targ, max_iter=300)
                    ))

                    pfim_rand_down = Foodweb(Matrix(
                        downsample(pfim_meta, :random;
                            target_co=C_targ, max_iter=300)
                    ))

                    niche_fw = Foodweb(:niche; S=size(pfim_meta, 1), C=C_targ)

                    # ATN: select the candidate with connectance closest to C_targ.
                    prods = map(==("primary"), string.(df.tiering))
                    atn_range = atn_build_range(df, prods)
                    best = argmin(abs.(getfield.(atn_range, :connectance) .- C_targ))
                    atn_fw = Foodweb(atn_range[best].adjacency)

                    initial_networks = Dict(
                        "down_link" => pfim_link_down,
                        "down_power" => pfim_power_down,
                        "down_niche" => pfim_niche_down,
                        "down_rand" => pfim_rand_down,
                        "niche" => niche_fw,
                        "atn" => atn_fw,
                        "metaweb" => Foodweb(pfim_meta)
                    )

                    # --- 2. Record creation networks ---

                    creation_metrics = Dict()

                    for (net_name, fw) in initial_networks

                        S_init = size(fw.A, 1)
                        C_init = sum(fw.A) / S_init^2
                        creation_metrics[net_name] = (S=S_init, C=C_init)

                        temp_metadata = DataFrame()
                        params = default_model(fw)

                        record_species_stage!(
                            temp_metadata, i, net_name, "creation",
                            fw.A, params.M, params, nothing, B_init
                        )

                        temp_metadata[!, :t_val] .= t
                        temp_metadata[!, :community] .= community_name
                        temp_metadata[!, :state] .= state

                        append!(species_metadata_store, temp_metadata; cols=:union)

                        push!(adjacency_store, (
                            t=t,
                            net_id=i,
                            community=community_name,
                            state=state,
                            net_type=net_name,
                            stage="creation",
                            S=S_init,
                            C=C_init,
                            adjacency=copy(fw.A)
                        ))
                    end

                    # --- 3. Burn-in and realisation ---

                    realised_networks = Dict()

                    for (net_name, fw) in initial_networks

                        realised = realise_network(
                            fw;
                            t=t,
                            threshold=survival_threshold,
                            B0=B_init
                        )

                        if realised !== nothing
                            realised_networks[net_name] = realised
                        end
                    end

                    if isempty(realised_networks)
                        @warn "No networks survived burn-in for $state, run $i, t=$t, community=$community_name. Skipping this state."
                        continue
                    end

                    # --- 4. Record realised networks ---

                    for (net_name, realised) in realised_networks

                        A_realised = realised.A
                        params = realised.params
                        final_biomasses = realised.biomasses
                        survivors = realised.survivors

                        S_realised = length(survivors)
                        C_realised = sum(A_realised) / S_realised^2

                        push!(adjacency_store, (
                            t=t,
                            net_id=i,
                            community=community_name,
                            state=state,
                            net_type=net_name,
                            stage="burnin",
                            S=S_realised,
                            C=C_realised,
                            adjacency=copy(A_realised)
                        ))

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

                        temp_metadata[!, :t_val] .= t
                        temp_metadata[!, :community] .= community_name
                        temp_metadata[!, :state] .= state

                        append!(species_metadata_store, temp_metadata; cols=:union)

                        push!(summary_rows, Dict(
                            :t_val => t,
                            :net_id => i,
                            :community => community_name,
                            :state => state,
                            :net_type => net_name,
                            :C_target => C_targ,
                            :S_creation => creation_metrics[net_name].S,
                            :C_creation => creation_metrics[net_name].C,
                            :C_realised => C_realised
                        ))

                        # Store each state's networks in its own archive.
                        network_archive[state][(community_name, i, net_name)] = Dict(
                            "generation" => Dict(
                                "community" => community_name,
                                "state" => state,
                                "burnin_time" => t,
                                "net_id" => i,
                                "C_target" => C_targ,
                                "bodymass" => copy(df.bodymass),
                                "B0" => copy(B_init),
                                "species" => String.(df.species),
                            ),
                            "creation" => Dict(
                                "A" => copy(initial_networks[net_name].A),
                                "S" => creation_metrics[net_name].S,
                                "C" => creation_metrics[net_name].C
                            ),
                            "realised" => Dict(
                                "A" => copy(A_realised),
                                "params" => params,
                                "biomasses" => final_biomasses,
                                "survivors" => survivors,
                                "S" => S_realised,
                                "C" => C_realised
                            )
                        )
                    end
                end # state
            end # community
        end # replicate
    end # burn-in time

    ################################################################################
    # 3. Export separate pre and post archives
    ################################################################################

    for state in states

        output_dir = "outputs/$state"
        mkpath(output_dir)

        state_species_metadata = filter(
            row -> row.state == state,
            species_metadata_store
        )

        state_adjacency_store = filter(
            row -> row.state == state,
            adjacency_store
        )

        state_summary_rows = filter(
            row -> row[:state] == state,
            summary_rows
        )

        archive = Dict(
            "metadata" => Dict(
                "seed" => 66,
                "n_networks" => n_networks,
                "survival_threshold" => survival_threshold,
                "burnin_times" => t_values,
                "state" => state
            ),
            "networks" => network_archive[state],
            "species_metadata" => state_species_metadata,
            "adjacency_store" => state_adjacency_store,
            "summary" => DataFrame(state_summary_rows)
        )

        JLD2.save_object(
            joinpath(output_dir, "realised_network_archive.jld2"),
            archive
        )

        CSV.write(
            joinpath(output_dir, "species_metadata.csv"),
            state_species_metadata
        )

        JLD2.save_object(
            joinpath(output_dir, "adjacency_matrices.jld2"),
            state_adjacency_store
        )

        # Keep the original downsample metadata export, now separated by state.
        metadata_rows = Dict[]

        for ((community, net_id, net_name), network) in network_archive[state]

            push!(metadata_rows, Dict(
                :state => state,
                :targ_Co => network["generation"]["C_target"],
                :net_id => net_id,
                :community => community,
                :net_name => net_name,
                :burnin_time => network["generation"]["burnin_time"]
            ))
        end

        CSV.write(
            joinpath(output_dir, "downsample_metadata.csv"),
            DataFrame(metadata_rows)
        )

        println("Saved $state archive with $(length(network_archive[state])) realised networks.")
    end

end # main

main()