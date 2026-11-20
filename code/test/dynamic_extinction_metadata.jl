using DifferentialEquations
using DataFrames
using EcologicalNetworksDynamics
using Extinctions
using FoodWebTools
using ProgressMeter
using SpeciesInteractionNetworks
using Test

include(joinpath(@__DIR__, "..", "libs", "internals.jl"))

const EXPECTED_SCENARIOS = Set([
    "degree_high", "degree_low",
    "vul_high", "vul_low",
    "gen_high", "gen_low",
    "bm_high", "bm_low",
    "rand_basal", "rand_consumer"
])

A = Bool[
    false false false;
    true  false false;
    false true  false
]
species_ids = ["basal_alpha", "consumer_beta", "consumer_gamma"]
params = default_model(
    Foodweb(A),
    BodyMass([1.0, 10.0, 100.0]),
    ClassicResponse(; h=2.0)
)

results = run_dynamic_extinctions(
    params,
    [1.0, 1.0, 1.0];
    t=1,
    show_progress=false,
    species_ids=species_ids,
    seed=7
)

curve_data = Dict(
    scenario => extinction_breakdown(result)
    for (scenario, result) in results
)
curves = export_curves(curve_data, "dynamic", 42)
add_curve_metadata!(curves, "test_community", "test_network", 1)

@testset "dynamic extinction function metadata" begin
    @test Set(keys(results)) == EXPECTED_SCENARIOS
    @test Set(curves.type) == Set(["dynamic"])
    @test Set(curves.scenario) == EXPECTED_SCENARIOS
    @test all(==(42), curves.net_id)
    @test Set(curves.community) == Set(["test_community"])
    @test Set(curves.net_name) == Set(["test_network"])
    @test Set(curves.extinction_t) == Set([1])

    for result in values(results), (network, survivors) in
        zip(result.networks, result.survivors)
        ids = String.(SpeciesInteractionNetworks.species(network))
        @test ids == species_ids[survivors]
        @test length(ids) == SpeciesInteractionNetworks.richness(network)
        @test length(unique(ids)) == length(ids)
        @test all(id -> id in species_ids, ids)
    end

    for curve in groupby(curves, :scenario)
        @test collect(curve.step) == collect(1:nrow(curve))
        @test all(diff(curve.total) .>= -1e-12)
        @test all(
            isapprox.(
                curve.total,
                curve.primary .+ curve.secondary;
                atol=1e-12,
                rtol=0
            )
        )
    end
end

println(
    "Exercised run_dynamic_extinctions with $(Threads.nthreads()) thread(s) " *
    "and validated $(nrow(curves)) exported curve rows."
)
