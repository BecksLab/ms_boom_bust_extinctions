using DifferentialEquations
using EcologicalNetworksDynamics
using FoodWebTools
using SpeciesInteractionNetworks
using Test

include("../libs/internals.jl")

const SCENARIOS = Dict(
    "degree_high" => (:degree, true),
    "degree_low" => (:degree, false),
    "vul_high" => (:vulnerability, true),
    "vul_low" => (:vulnerability, false),
    "gen_high" => (:generality, true),
    "gen_low" => (:generality, false),
    "bm_high" => (:bodymass, true),
    "bm_low" => (:bodymass, false),
    "rand_basal" => (:random_basal, true),
    "rand_consumer" => (:random_consumer, true)
)

A = falses(7, 7) |> Matrix
A[4, 1] = true
A[4, 2] = true
A[4, 5] = true
A[4, 7] = true
A[5, 1] = true
A[5, 2] = true
A[5, 4] = true
A[6, 1] = true
A[6, 2] = true
A[6, 3] = true
A[7, 1] = true

species_ids = [
    "basal_p", "basal_q", "basal_r",
    "consumer_a", "consumer_b", "consumer_c", "consumer_d"
]
body_masses = collect(1.0:7.0)
biomasses = ones(length(species_ids))
consumer_indices = findall(vec(sum(A; dims=2)) .> 0)
basal_indices = setdiff(collect(eachindex(species_ids)), consumer_indices)
params = default_model(
    Foodweb(A),
    BodyMass(body_masses),
    ClassicResponse(; h=2.0)
)

function assert_dynamic_order(result, scenario)
    criterion, descending = SCENARIOS[scenario]
    @test result.primary == collect(0:length(result.primary_removed))
    @test length(result.survivors) == length(result.primary_removed) + 1
    @test length(result.networks) == length(result.survivors)
    @test length(unique(result.primary_removed)) == length(result.primary_removed)
    @test length(result.secondary_removed) == length(result.primary_removed)

    for (web, survivors) in zip(result.networks, result.survivors)
        ids = String.(SpeciesInteractionNetworks.species(web))
        @test ids == species_ids[survivors]
        @test length(ids) == SpeciesInteractionNetworks.richness(web)
    end

    for (step, removed_id) in enumerate(result.primary_removed)
        current_indices = result.survivors[step]
        current_A = A[current_indices, current_indices]
        current_consumers = intersect(current_indices, consumer_indices)
        current_basals = intersect(current_indices, basal_indices)
        removed_index = findfirst(==(removed_id), species_ids)
        @test removed_index in current_indices

        if criterion == :random_basal
            @test removed_index in current_basals
        elseif criterion == :random_consumer
            @test removed_index in current_consumers
        else
            @test removed_index in current_consumers
            candidates = current_consumers
            local_candidates = [
                findfirst(==(index), current_indices) for index in candidates
            ]
            generality = vec(sum(current_A; dims=2))
            vulnerability = vec(sum(current_A; dims=1))
            degree = generality .+ vulnerability
            scores = criterion == :generality ? generality :
                     criterion == :vulnerability ? vulnerability :
                     criterion == :degree ? degree : body_masses[current_indices]
            candidate_scores = scores[local_candidates]
            target = descending ? maximum(candidate_scores) : minimum(candidate_scores)
            @test scores[findfirst(==(removed_index), current_indices)] == target
        end

        next_indices = result.survivors[step + 1]
        @test !(removed_index in next_indices)
        @test issubset(Set(next_indices), setdiff(Set(current_indices), Set([removed_index])))
        @test !(removed_id in result.secondary_removed[step])
        @test all(id -> id in species_ids, result.secondary_removed[step])
        @test isempty(intersect(
            Set(result.secondary_removed[step]), Set(next_indices)
        ))
    end
end

@testset "dynamic extinction order and identity" begin
    @test_throws ArgumentError dynamic_extinction_adaptive(
        params, biomasses;
        criterion=:degree,
        t=1,
        show_progress=false,
        species_ids=species_ids[1:end-1]
    )
    @test_throws ArgumentError dynamic_extinction_adaptive(
        params, biomasses;
        criterion=:degree,
        t=1,
        show_progress=false,
        species_ids=["duplicate", "duplicate", species_ids[3:end]...]
    )

    results = run_dynamic_extinctions(
        params,
        biomasses;
        t=1,
        show_progress=false,
        species_ids=species_ids,
        seed=7
    )
    @test Set(keys(results)) == Set(keys(SCENARIOS))

    for (scenario, result) in results
        assert_dynamic_order(result, scenario)
    end

    # This fixture makes the degree ranking change after the first removal:
    # initially B ranks above C, but removing A reduces B's degree and makes C
    # the next target under adaptive re-ranking.
    initial_scores = vec(sum(A; dims=2)) .+ vec(sum(A; dims=1))
    static_order = sort(consumer_indices; by=index -> -initial_scores[index])
    degree_high_order = results["degree_high"].primary_removed
    @test species_ids[static_order[1:2]] == ["consumer_a", "consumer_b"]
    @test degree_high_order[1:2] == ["consumer_a", "consumer_c"]
    @test results["bm_high"].primary_removed ==
          ["consumer_d", "consumer_c", "consumer_b", "consumer_a"]
    @test results["bm_low"].primary_removed ==
          ["consumer_a", "consumer_b", "consumer_c", "consumer_d"]
end

@testset "dynamic threshold-driven secondary extinction" begin
    chain_A = Bool[
        false false false;
        true  false false;
        false true  false
    ]
    chain_ids = ["basal_alpha", "consumer_beta", "consumer_gamma"]
    chain_params = default_model(
        Foodweb(chain_A),
        BodyMass([1.0, 10.0, 100.0]),
        ClassicResponse(; h=2.0)
    )
    threshold_result = dynamic_extinction_adaptive(
        chain_params,
        ones(3);
        criterion=:bodymass,
        descending=true,
        t=1,
        survival_threshold=1e100,
        show_progress=false,
        species_ids=chain_ids
    )

    @test threshold_result.primary_removed == ["consumer_gamma"]
    @test threshold_result.secondary_removed == [
        ["basal_alpha", "consumer_beta"]
    ]
    @test threshold_result.survivors == [[1, 2, 3]]

end

@testset "consumers starve after prey loss" begin
    chain_A = Bool[
        false false false;
        true  false false;
        false true  false
    ]
    chain_ids = ["basal_alpha", "consumer_beta", "consumer_gamma"]
    chain_params = default_model(
        Foodweb(chain_A),
        BodyMass([1.0, 10.0, 100.0]),
        ClassicResponse(; h=2.0)
    )
    starvation_result = dynamic_extinction_adaptive(
        chain_params,
        ones(3);
        criterion=:random_basal,
        t=5000,
        survival_threshold=1e-12,
        show_progress=false,
        species_ids=chain_ids
    )
    @test starvation_result.primary_removed == ["basal_alpha"]
    @test starvation_result.secondary_removed == [
        ["consumer_beta", "consumer_gamma"]
    ]
    @test starvation_result.survivors == [[1, 2, 3]]
end

println(
    "Validated adaptive extinction order, prey-loss starvation, and species identity."
)
