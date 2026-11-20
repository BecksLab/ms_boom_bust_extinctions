using EcologicalNetworksDynamics
using FoodWebTools
using SpeciesInteractionNetworks
using Test

include("../libs/internals.jl")
include("../libs/tss.jl")

species_ids = ["species_A", "species_B", "species_C", "species_D"]
creation_adjacency = Bool[
    false true  false false;
    false false true  false;
    false false false true;
    false false false false
]
creation_network = build_network(creation_adjacency, species_ids)

post_species = ["species_A", "species_C", "species_D"]
post_adjacency = Bool[
    false true  false;
    false false true;
    true  false false
]
post_network = build_network(post_adjacency, post_species)

simulated_species = ["species_D", "species_A", "species_C"]
simulated_adjacency_canonical = Bool[
    false true  false;
    false false true;
    false false false
]
simulated_order = [3, 1, 2]
simulated_adjacency =
    simulated_adjacency_canonical[simulated_order, simulated_order]
simulated_network = build_network(simulated_adjacency, simulated_species)

pre_data = Dict(
    "generation" => Dict("species" => species_ids),
    "creation" => Dict("A" => Matrix(creation_network.edges.edges)),
    "realised" => Dict(
        "survivors" => [1, 2],
        "A" => creation_adjacency[[1, 2], [1, 2]]
    )
)
post_data = Dict(
    "generation" => Dict("species" => species_ids),
    "realised" => Dict(
        "survivors" => [1, 3, 4],
        "A" => Matrix(post_network.edges.edges)
    )
)
simulation_data = Dict(
    "species" => String.(SpeciesInteractionNetworks.species(simulated_network)),
    "A" => Matrix(simulated_network.edges.edges),
    "actual_S" => length(simulated_species)
)

@testset "species identity from network creation through TSS" begin
    @test String.(SpeciesInteractionNetworks.species(creation_network)) == species_ids
    @test Matrix(creation_network.edges.edges) == creation_adjacency
    @test realised_species(post_data) == post_species
    @test String.(SpeciesInteractionNetworks.species(simulated_network)) == simulated_species

    scores = calculate_tss(pre_data, post_data, simulation_data)
    @test scores.n_original == 4
    @test scores.n_realised == 3
    @test scores.n_simulated == 3
    @test scores.n_shared == 3
    @test scores.tp_node == 3
    @test scores.tn_node == 1
    @test scores.fp_node == 0
    @test scores.fn_node == 0
    @test scores.tss_node == 1.0
    @test scores.tp_link == 2
    @test scores.tn_link == 6
    @test scores.fp_link == 0
    @test scores.fn_link == 1
    @test scores.tss_link ≈ 2 / 3

    unknown_species = merge(
        simulation_data,
        Dict(
            "species" => ["not_in_PRE"],
            "A" => falses(1, 1),
            "actual_S" => 1
        )
    )
    @test_throws ErrorException calculate_tss(pre_data, post_data, unknown_species)
end
