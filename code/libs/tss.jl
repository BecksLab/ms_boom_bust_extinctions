function realised_species(network_data)
    generation_species = String.(network_data["generation"]["species"])
    survivors = network_data["realised"]["survivors"]
    species = generation_species[survivors]
    adjacency = network_data["realised"]["A"]

    size(adjacency) == (length(species), length(species)) ||
        error("Species IDs do not match realised adjacency dimensions.")
    length(unique(species)) == length(species) ||
        error("Species IDs must be unique.")

    return species
end

function validate_network(adjacency, species_ids, label)
    size(adjacency) == (length(species_ids), length(species_ids)) ||
        error("$label: species IDs do not match adjacency dimensions.")
    length(unique(species_ids)) == length(species_ids) ||
        error("$label: species IDs are not unique.")
    all(value -> value in (0, 1, false, true), adjacency) ||
        error("$label: adjacency matrix must be binary.")
end

function ordered_adjacency(adjacency, species_ids, requested_species)
    positions = Dict(id => i for (i, id) in enumerate(species_ids))
    all(id -> haskey(positions, id), requested_species) ||
        error("Requested species are missing from the network.")
    indices = [positions[id] for id in requested_species]
    return adjacency[indices, indices]
end

function tss_component(tp, fn, tn, fp)
    sensitivity_denominator = tp + fn
    specificity_denominator = tn + fp
    (sensitivity_denominator == 0 || specificity_denominator == 0) &&
        return missing
    return tp / sensitivity_denominator + tn / specificity_denominator - 1
end

function calculate_tss(
    original_adjacency, original_species,
    realised_adjacency, realised_species,
    simulated_adjacency, simulated_species
)
    validate_network(original_adjacency, original_species, "Original PRE")
    validate_network(realised_adjacency, realised_species, "Realised POST")
    validate_network(simulated_adjacency, simulated_species, "Simulated PRE")

    original_set = Set(original_species)
    realised_set = Set(realised_species)
    simulated_set = Set(simulated_species)
    issubset(realised_set, original_set) ||
        error("POST species are not a subset of PRE species.")
    issubset(simulated_set, original_set) ||
        error("Simulated species are not a subset of PRE species.")

    node_tp = length(intersect(simulated_set, realised_set))
    node_fn = length(setdiff(realised_set, simulated_set))
    extinct_species = setdiff(original_set, realised_set)
    node_tn = length(setdiff(extinct_species, simulated_set))
    node_fp = length(setdiff(simulated_set, realised_set))
    node_tss = tss_component(node_tp, node_fn, node_tn, node_fp)

    shared = sort!(collect(intersect(simulated_set, realised_set)))
    link_counts = (tp=0, tn=0, fp=0, fn=0)
    link_tss = missing

    if length(shared) > 1
        simulated = ordered_adjacency(
            simulated_adjacency, simulated_species, shared
        )
        realised = ordered_adjacency(
            realised_adjacency, realised_species, shared
        )
        link_counts = (
            tp=sum((simulated .== 1) .& (realised .== 1)),
            tn=sum((simulated .== 0) .& (realised .== 0)),
            fp=sum((simulated .== 1) .& (realised .== 0)),
            fn=sum((simulated .== 0) .& (realised .== 1))
        )
        link_tss = tss_component(
            link_counts.tp, link_counts.fn,
            link_counts.tn, link_counts.fp
        )
    end

    return (
        tss_node=node_tss,
        tss_link=link_tss,
        tp_node=node_tp,
        tn_node=node_tn,
        fp_node=node_fp,
        fn_node=node_fn,
        tp_link=link_counts.tp,
        tn_link=link_counts.tn,
        fp_link=link_counts.fp,
        fn_link=link_counts.fn,
        n_original=length(original_species),
        n_simulated=length(simulated_species),
        n_realised=length(realised_species),
        n_shared=length(shared)
    )
end

function calculate_tss(pre_data::AbstractDict, post_data::AbstractDict, simulation::AbstractDict)
    original_species = String.(pre_data["generation"]["species"])
    post_species = realised_species(post_data)
    simulated_species = String.(simulation["species"])
    length(simulated_species) == simulation["actual_S"] ||
        error("Simulated richness does not match stored species IDs.")

    return calculate_tss(
        Matrix(pre_data["creation"]["A"]), original_species,
        Matrix(post_data["realised"]["A"]), post_species,
        Matrix(simulation["A"]), simulated_species
    )
end
