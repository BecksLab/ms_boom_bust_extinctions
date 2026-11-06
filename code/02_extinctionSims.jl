#=
################################################################################
RUN EXTINCTION SIMULATIONS
################################################################################

Purpose
-------
Load previously generated realised ecological networks and perform extinction
experiments.

Network generation and burn-in are deliberately separated from this workflow.
This allows extinction scenarios to be modified, expanded, or rerun without
regenerating ecological networks.

Pipeline
--------
1. Load realised network archive
2. Select realised networks
3. Run topological extinction experiments
4. Run dynamic extinction experiments
5. Calculate robustness metrics
6. Export extinction curves and summaries

Inputs
------
realised_network_archive.jld2

Outputs
-------
outputs/extinction_summary.csv       # one row per network, wide robustness metrics
outputs/extinction_curves.csv        # one row per extinction step
outputs/tss_networks/species.csv     # species identity and matrix position

For threaded dynamic scenarios, start Julia with JULIA_NUM_THREADS set, e.g.
JULIA_NUM_THREADS=auto julia --project=. 02_extinctionSims.jl


################################################################################
=#

# --- 1. General setup ---------------------------------------------------------

using CSV
using DataFrames
using DifferentialEquations
using EcologicalNetworksDynamics
using Extinctions
using FoodWebTools
using JLD2
using ProgressMeter
using SpeciesInteractionNetworks
using Statistics

include("libs/internals.jl")

import Random
Random.seed!(66)


# ============================================================
# 1. LOAD PRE AND POST NETWORK ARCHIVES
# ============================================================

pre_archive = JLD2.load_object(
    "outputs/pre/realised_network_archive.jld2"
)

post_archive = JLD2.load_object(
    "outputs/post/realised_network_archive.jld2"
)

pre_networks = pre_archive["networks"]
post_networks = post_archive["networks"]

survival_threshold = pre_archive["metadata"]["survival_threshold"]


# ============================================================
# 2. STORAGE
# ============================================================

summary_rows = Dict[]

topo_curve_store = DataFrame()
dyn_curve_store = DataFrame()

# Extracted networks for TSS analyses
tss_pre_networks = Dict()

extinction_times = [5000]


# ============================================================
# 3. HELPER FUNCTIONS
# ============================================================

# Find the closest richness in an extinction trajectory.
# If two networks are equally close, the first is selected.

function extract_at_richness(result, target_S)

    richness_values = SpeciesInteractionNetworks.richness.(
        result.networks
    )

    isempty(richness_values) && return nothing

    idx = argmin(abs.(richness_values .- target_S))

    return (
        network=result.networks[idx],
        richness=richness_values[idx],
        richness_difference=richness_values[idx] - target_S,
        step=idx
    )
end


# Add network metadata to exported extinction curves.

function add_curve_metadata!(
    df,
    community,
    net_name,
    extinction_t
)
    df[!, :community] = fill(community, nrow(df))
    df[!, :net_name] = fill(net_name, nrow(df))
    df[!, :extinction_t] = fill(extinction_t, nrow(df))

    return df
end


# Extract and store closest-richness networks from a set
# of extinction scenarios.

function extract_scenario_networks!(
    storage,
    results,
    target_S,
    community,
    net_id,
    net_name,
    extinction_t,
    extinction_type
)

    for (scenario, result) in results

        extracted = extract_at_richness(result, target_S)

        if extracted === nothing
            continue
        end

        key = (
            community,
            net_id,
            net_name,
            extinction_t,
            extinction_type,
            scenario
        )

        storage[key] = Dict(
            "target_S" => target_S,
            "actual_S" => extracted.richness,
            "richness_difference" => extracted.richness_difference,
            "step" => extracted.step,
            "species" => String.(
                SpeciesInteractionNetworks.species(extracted.network)
            ),
            "A" => Matrix(extracted.network.edges.edges)
        )
    end
end


# ============================================================
# 4. EXTINCTION SIMULATIONS
# ============================================================

network_keys = sort!(collect(keys(pre_networks)))
println("Dynamic scenarios will use $(Threads.nthreads()) Julia thread(s).")

for (time_index, extinction_t) in enumerate(extinction_times)

    for (network_index, key) in enumerate(network_keys)
        community, net_id, net_name = key
        network = pre_networks[key]
        network_seed = 66 + (time_index - 1) * length(network_keys) + network_index
        Random.seed!(network_seed)

        # ----------------------------------------------------
        # Match corresponding POST network
        # ----------------------------------------------------

        key = (community, net_id, net_name)

        if !haskey(post_networks, key)

            @warn "No matching POST network found: $key"

            continue
        end

        post_network = post_networks[key]

        target_S = post_network["realised"]["S"]

        # ----------------------------------------------------
        # Extract PRE network information
        # ----------------------------------------------------

        creation = network["creation"]
        realised = network["realised"]
        generation_species = String.(network["generation"]["species"])
        realised_species = generation_species[realised["survivors"]]

        params = realised["params"]
        biomasses = realised["biomasses"]


        # ====================================================
        # 4A. TOPOLOGICAL EXTINCTIONS: CREATION NETWORK
        # ====================================================

        A_creation = Matrix(creation["A"])

        N_creation = build_network(A_creation, generation_species)

        topo_creation_results = run_topological_extinctions(
            N_creation,
            network["generation"]["bodymass"]
        )

        topo_creation_robustness = compute_robustness(
            topo_creation_results
        )

        topo_creation_curves = Dict(
            k => extinction_breakdown(v)
            for (k, v) in topo_creation_results
        )

        topo_creation_df = export_curves(
            topo_creation_curves,
            "topological_creation",
            net_id
        )

        add_curve_metadata!(
            topo_creation_df,
            community,
            net_name,
            extinction_t
        )

        append!(
            topo_curve_store,
            topo_creation_df;
            cols=:union
        )


        # ====================================================
        # 4B. TOPOLOGICAL EXTINCTIONS: REALISED NETWORK
        # ====================================================

        A_realised = Matrix(realised["A"])

        N_realised = build_network(A_realised, realised_species)

        topo_realised_results = run_topological_extinctions(
            N_realised,
            params.M[realised["survivors"]]
        )

        topo_realised_robustness = compute_robustness(
            topo_realised_results
        )

        topo_realised_curves = Dict(
            k => extinction_breakdown(v)
            for (k, v) in topo_realised_results
        )

        topo_realised_df = export_curves(
            topo_realised_curves,
            "topological_realised",
            net_id
        )

        add_curve_metadata!(
            topo_realised_df,
            community,
            net_name,
            extinction_t
        )

        append!(
            topo_curve_store,
            topo_realised_df;
            cols=:union
        )


        # ----------------------------------------------------
        # Extract PRE topological networks at closest POST S
        # ----------------------------------------------------

        extract_scenario_networks!(
            tss_pre_networks,
            topo_realised_results,
            target_S,
            community,
            net_id,
            net_name,
            extinction_t,
            "topological"
        )


        # ====================================================
        # 4C. DYNAMIC EXTINCTIONS
        # ====================================================

        realised_params = default_model(
            Foodweb(A_realised),
            BodyMass(params.M[realised["survivors"]]),
            ClassicResponse(; h=params.h)
        )

        dyn_results = run_dynamic_extinctions(
            realised_params,
            biomasses[realised["survivors"]];
            t=extinction_t,
            species_ids=realised_species,
            seed=network_seed * 10
        )

        dyn_robustness = compute_robustness(
            dyn_results
        )

        dyn_curves = Dict(
            k => extinction_breakdown(v)
            for (k, v) in dyn_results
        )

        dyn_df = export_curves(
            dyn_curves,
            "dynamic",
            net_id
        )

        add_curve_metadata!(
            dyn_df,
            community,
            net_name,
            extinction_t
        )

        append!(
            dyn_curve_store,
            dyn_df;
            cols=:union
        )


        # ----------------------------------------------------
        # Extract PRE dynamic networks at closest POST S
        # ----------------------------------------------------

        extract_scenario_networks!(
            tss_pre_networks,
            dyn_results,
            target_S,
            community,
            net_id,
            net_name,
            extinction_t,
            "dynamic"
        )


        # ====================================================
        # 4D. STORE SUMMARY
        # ====================================================

        summary = Dict{Symbol, Any}(
            :community => community,
            :net_id => net_id,
            :net_type => net_name,
            :extinction_time => extinction_t,
            :C_target => network["generation"]["C_target"],
            :S_post_target => target_S,
            :S_creation => creation["S"],
            :C_creation => creation["C"],
            :S_realised => realised["S"],
            :C_realised => realised["C"]
        )

        for (scenario, value) in topo_creation_robustness
            summary[Symbol("topo_creation_", scenario)] = value
        end
        for (scenario, value) in topo_realised_robustness
            summary[Symbol("topo_realised_", scenario)] = value
        end
        for (scenario, value) in dyn_robustness
            summary[Symbol("dyn_realised_", scenario)] = value
        end
        push!(summary_rows, summary)

    end
end


# ============================================================
# 5. EXPORT EXTINCTION RESULTS
# ============================================================

curve_df = vcat(
    topo_curve_store,
    dyn_curve_store;
    cols=:union
)

CSV.write(
    "outputs/extinction_summary.csv",
    DataFrame(summary_rows)
)

CSV.write(
    "outputs/extinction_curves.csv",
    curve_df
)


# ============================================================
# 6. EXPORT TSS NETWORKS
# ============================================================

JLD2.save_object(
    "outputs/tss_networks/pre_at_post_richness.jld2",
    tss_pre_networks
)

species_rows = Dict[]
for (key, network) in tss_pre_networks
    community, net_id, net_name, extinction_time, extinction_type, scenario = key
    for (position, species) in enumerate(network["species"])
        push!(species_rows, Dict(
            :community => community,
            :net_id => net_id,
            :net_name => net_name,
            :extinction_time => extinction_time,
            :extinction_type => extinction_type,
            :scenario => scenario,
            :network_role => "simulated_pre",
            :position => position,
            :species_id => species
        ))
    end
end

CSV.write(
    "outputs/tss_networks/species.csv",
    DataFrame(species_rows)
)