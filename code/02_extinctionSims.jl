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
extinction_summary.csv
extinction_curves.csv


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
tss_post_networks = Dict()

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
            "network" => extracted.network
        )
    end
end


# ============================================================
# 4. EXTINCTION SIMULATIONS
# ============================================================

for extinction_t in extinction_times

    for ((community, net_id, net_name), network) in pre_networks

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
        # Store POST realised network
        # ----------------------------------------------------

        tss_post_networks[key] = Dict(
            "target_S" => target_S,
            "A" => copy(post_network["realised"]["A"]),
            "network" => build_network(
                Matrix(post_network["realised"]["A"])
            )
        )

        # ----------------------------------------------------
        # Extract PRE network information
        # ----------------------------------------------------

        creation = network["creation"]
        realised = network["realised"]

        params = realised["params"]
        biomasses = realised["biomasses"]


        # ====================================================
        # 4A. TOPOLOGICAL EXTINCTIONS: CREATION NETWORK
        # ====================================================

        A_creation = Matrix(creation["A"])

        N_creation = build_network(A_creation)

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

        N_realised = build_network(A_realised)

        topo_realised_results = run_topological_extinctions(
            N_realised,
            params.body_mass
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

        dyn_results = run_dynamic_extinctions(
            params,
            biomasses;
            t=extinction_t
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

        push!(
            summary_rows,
            Dict(
                "community" => community,
                "net_id" => net_id,
                "net_name" => net_name,
                "extinction_t" => extinction_t,
                "post_target_S" => target_S,
                "topological_creation_robustness" => topo_creation_robustness,
                "topological_realised_robustness" => topo_realised_robustness,
                "dynamic_robustness" => dyn_robustness
            )
        )

    end
end


# ============================================================
# 5. EXPORT EXTINCTION RESULTS
# ============================================================

summary_df = DataFrame(summary_rows)

curve_df = vcat(
    topo_curve_store,
    dyn_curve_store;
    cols=:union
)

CSV.write(
    "outputs/extinction_summary.csv",
    summary_df
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

JLD2.save_object(
    "outputs/tss_networks/post_realised.jld2",
    tss_post_networks
)