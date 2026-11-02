
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

import Random

include("libs/internals.jl")

Random.seed!(66)


# ============================================================
# 1. LOAD NETWORK ARCHIVES
# ============================================================

pre_archive = JLD2.load_object(
    "outputs/pre/realised_network_archive.jld2"
)

post_archive = JLD2.load_object(
    "outputs/post/realised_network_archive.jld2"
)

pre_networks = pre_archive["networks"]
post_networks = post_archive["networks"]


# ============================================================
# 2. SETTINGS AND STORAGE
# ============================================================

extinction_times = [5000]

tss_rows = NamedTuple[]

# ============================================================
# 3. HELPER FUNCTIONS
# ============================================================

# Construct a network while preserving species identifiers.

function labelled_network(A, species_ids)

    nodes = Unipartite(Symbol.(species_ids))
    edges = Binary(.!iszero.(A))

    return SpeciesInteractionNetwork(nodes, edges)
end

# Find the network with richness closest to the POST target.

function closest_richness_index(networks, target_S)

    isempty(networks) && return nothing

    richness_values = SpeciesInteractionNetworks.richness.(
        networks
    )

    return argmin(abs.(richness_values .- target_S))
end

# Calculate and record TSS for topological extinction scenarios.

function record_topological_tss!(
    rows,
    results,
    N_real,
    N_original,
    target_S,
    community,
    net_id,
    net_name,
    extinction_t
)

    for (scenario, result) in results

        networks = result.networks

        idx = closest_richness_index(
            networks,
            target_S
        )

        isnothing(idx) && continue

        N_sim = networks[idx]

        sim_S = SpeciesInteractionNetworks.richness(N_sim)

        # Calculate TSS
        tss_link, tss_node = TSS(
            N_real,
            N_sim,
            N_original
        )

        push!(
            rows,
            (
                community=community,
                net_id=net_id,
                net_name=net_name,
                extinction_t=extinction_t,
                extinction_type="topological",
                scenario=scenario,
                target_richness=target_S,
                sim_richness=sim_S,
                richness_difference=sim_S - target_S,
                extinction_step=idx,
                tss_link=tss_link,
                tss_node=tss_node
            )
        )
    end
end

# Calculate and record TSS for dynamic extinction scenarios.
# Dynamic extinction networks are reconstructed using their
# retained original species indices.

function record_dynamic_tss!(
    rows,
    results,
    A_original,
    original_species,
    N_real,
    N_original,
    target_S,
    community,
    net_id,
    net_name,
    extinction_t
)

    for (scenario, result) in results

        networks = result.networks

        idx = closest_richness_index(
            networks,
            target_S
        )

        isnothing(idx) && continue

        # Original PRE species indices retained at this step
        alive = result.survivors[idx]

        isempty(alive) && continue

        # Reconstruct network with original species identities
        A_sim = A_original[alive, alive]

        sim_species = original_species[alive]

        N_sim = labelled_network(
            A_sim,
            sim_species
        )

        sim_S = SpeciesInteractionNetworks.richness(N_sim)

        # Calculate TSS
        tss_link, tss_node = TSS(
            N_real,
            N_sim,
            N_original
        )

        push!(
            rows,
            (community=community,
                net_id=net_id,
                net_name=net_name,
                extinction_t=extinction_t,
                extinction_type="dynamic",
                scenario=scenario,
                target_richness=target_S,
                sim_richness=sim_S,
                richness_difference=sim_S - target_S,
                extinction_step=idx,
                tss_link=tss_link,
                tss_node=tss_node
            )
        )
    end
end


# ============================================================
# 4. MAIN TSS WORKFLOW
# ============================================================

for extinction_t in extinction_times

    for ((community, net_id, net_name), pre_network) in pre_networks

        # ----------------------------------------------------
        # 4A. MATCH CORRESPONDING POST NETWORK
        # ----------------------------------------------------

        key = (community, net_id, net_name)

        if !haskey(post_networks, key)

            @warn "No matching POST network found: $key"

            continue
        end

        post_network = post_networks[key]


        # ----------------------------------------------------
        # 4B. GET PRE AND POST SPECIES IDENTITIES
        # ----------------------------------------------------

        pre_generation_species = Symbol.(
            pre_network["generation"]["species"]
        )

        post_generation_species = Symbol.(
            post_network["generation"]["species"]
        )


        # Realised PRE network
        pre_realised = pre_network["realised"]

        pre_survivors = pre_realised["survivors"]

        pre_species = pre_generation_species[pre_survivors]

        A_pre = Matrix(pre_realised["A"])

        N_original = labelled_network(
            A_pre,
            pre_species
        )


        # Realised POST network
        post_realised = post_network["realised"]

        post_survivors = post_realised["survivors"]

        post_species = post_generation_species[post_survivors]

        A_post = Matrix(post_realised["A"])

        N_real = labelled_network(
            A_post,
            post_species
        )

        # POST realised richness is the target
        target_S = SpeciesInteractionNetworks.richness(N_real)


        # ----------------------------------------------------
        # 4C. BUILD REALISED PRE MODEL PARAMETERS
        # ----------------------------------------------------

        # The stored params refer to the original creation network.
        # Subset body masses and biomasses to the realised PRE species.

        params_full = pre_realised["params"]

        M_pre = params_full.M[pre_survivors]

        B_pre = pre_realised["biomasses"][pre_survivors]

        fw_pre = Foodweb(A_pre)

        params_pre = default_model(
            fw_pre,
            BodyMass(M_pre),
            ClassicResponse(; h=params_full.h)
        )


        # ----------------------------------------------------
        # 4D. TOPOLOGICAL EXTINCTIONS
        # ----------------------------------------------------

        topo_results = run_topological_extinctions(
            N_original,
            M_pre
        )

        record_topological_tss!(
            tss_rows,
            topo_results,
            N_real,
            N_original,
            target_S,
            community,
            net_id,
            net_name,
            extinction_t
        )


        # ----------------------------------------------------
        # 4E. DYNAMIC EXTINCTIONS
        # ----------------------------------------------------

        dyn_results = run_dynamic_extinctions(
            params_pre,
            B_pre;
            t=extinction_t
        )

        record_dynamic_tss!(
            tss_rows,
            dyn_results,
            A_pre,
            pre_species,
            N_real,
            N_original,
            target_S,
            community,
            net_id,
            net_name,
            extinction_t
        )

    end
end


# ============================================================
# 5. EXPORT TSS RESULTS
# ============================================================

tss_df = DataFrame(tss_rows)

CSV.write(
    "outputs/tss_results.csv",
    tss_df
)