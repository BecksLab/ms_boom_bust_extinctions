using DifferentialEquations
using EcologicalNetworksDynamics
using FoodWebTools
using SpeciesInteractionNetworks
using Test

include("../libs/internals.jl")

function run_reproducibility_fixture()
    adjacency = falses(5, 5)
    adjacency[4, 1:3] .= true
    adjacency[5, 1:3] .= true
    species_ids = [
        "basal_1", "basal_2", "basal_3", "consumer_1", "consumer_2"
    ]
    params = default_model(
        Foodweb(adjacency),
        BodyMass([1.0, 2.0, 3.0, 10.0, 20.0]),
        ClassicResponse(; h=2.0)
    )
    run_dynamic_extinctions(
        params,
        ones(length(species_ids));
        t=1,
        show_progress=false,
        species_ids,
        seed=712
    )
end

function result_trace(results, species_ids)
    lines = String[]
    for scenario in sort(collect(keys(results)))
        result = results[scenario]
        secondary = [
            join(removed, ",") for removed in result.secondary_removed
        ]
        survivors = [
            join((species_ids[index] for index in group), ",")
            for group in result.survivors
        ]
        push!(
            lines,
            join((
                scenario,
                join(result.primary_removed, ","),
                join(secondary, ";"),
                join(survivors, ";")
            ), "\t")
        )
    end
    join(lines, "\n")
end

if "--worker" in ARGS
    species_ids = [
        "basal_1", "basal_2", "basal_3", "consumer_1", "consumer_2"
    ]
    first_trace = result_trace(run_reproducibility_fixture(), species_ids)
    repeated_trace = result_trace(run_reproducibility_fixture(), species_ids)
    first_trace == repeated_trace ||
        error("Dynamic extinction results changed between runs with the same seed")
    println(first_trace)
else
    @testset "dynamic extinction reproducibility across random runs and threads" begin
        project_dir = normpath(joinpath(@__DIR__, ".."))
        test_file = abspath(@__FILE__)

        function run_worker(thread_count)
            command = `$(Base.julia_cmd()) --startup-file=no --project=$project_dir --threads=$thread_count $test_file --worker`
            read(command, String)
        end

        one_thread_trace = run_worker(1)
        two_thread_trace = run_worker(2)

        @test occursin("rand_basal", one_thread_trace)
        @test occursin("rand_consumer", one_thread_trace)
        @test one_thread_trace == two_thread_trace
    end
end
