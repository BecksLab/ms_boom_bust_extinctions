using Extinctions
using Random
using SpeciesInteractionNetworks
using Test

include("../libs/internals.jl")

const EXPECTED_SCENARIOS = Set([
    "degree_high", "degree_low",
    "vul_high", "vul_low",
    "gen_high", "gen_low",
    "bm_high", "bm_low",
    "rand_basal", "rand_consumer"
])

function ranked_fixture()
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
    return A, species_ids
end

web_species(web) = String.(SpeciesInteractionNetworks.species(web))

function consumer_scores(A, original_ids, current_ids, consumer_ids, metric)
    original_positions = [findfirst(==(id), original_ids) for id in current_ids]
    current_A = A[original_positions, original_positions]
    generality = vec(sum(current_A; dims=2))
    vulnerability = vec(sum(current_A; dims=1))
    degree = generality .+ vulnerability
    values = metric == :generality ? generality :
             metric == :vulnerability ? vulnerability : degree

    return Dict(
        id => values[findfirst(==(id), current_ids)]
        for id in current_ids if id in consumer_ids
    )
end

function consumer_removal_order(result, consumer_ids)
    order = String[]
    for (before, after) in zip(result.networks, result.networks[2:end])
        before_ids = Set(web_species(before))
        after_ids = Set(web_species(after))
        append!(order, [id for id in consumer_ids if id in before_ids && !(id in after_ids)])
    end
    return order
end

function check_ranked_trajectory(result, A, original_ids, consumer_ids, metric, descending)
    @test result.primary == collect(0:(length(result.networks) - 1))

    for (before, after) in zip(result.networks, result.networks[2:end])
        current_ids = web_species(before)
        next_ids = web_species(after)
        current_set = Set(current_ids)
        next_set = Set(next_ids)
        removed_consumers = [
            id for id in consumer_ids if id in current_set && !(id in next_set)
        ]
        @test length(removed_consumers) == 1
        isempty(removed_consumers) && continue
        selected = only(removed_consumers)

        scores = consumer_scores(A, original_ids, current_ids, consumer_ids, metric)
        target = descending ? maximum(values(scores)) : minimum(values(scores))
        @test scores[selected] == target

        # Apart from the selected consumer, only species made completely
        # disconnected by that removal should disappear in this fixture.
        after_primary = [id for id in current_ids if id != selected]
        positions = [findfirst(==(id), original_ids) for id in after_primary]
        remaining_A = A[positions, positions]
        linked = (vec(sum(remaining_A; dims=1)) .+
                  vec(sum(remaining_A; dims=2))) .> 0
        expected_after = after_primary[linked]
        @test Set(next_ids) == Set(expected_after)
    end
end

function random_fixture()
    A = falses(5, 5) |> Matrix
    A[4, 1:3] .= true
    A[5, 1:3] .= true
    ids = ["basal_1", "basal_2", "basal_3", "consumer_1", "consumer_2"]
    return A, ids
end

A, species_ids = ranked_fixture()
consumer_ids = species_ids[4:7]
body_masses = [100.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
network = build_network(A, species_ids)
Random.seed!(12)
results = run_topological_extinctions(network, body_masses)

@testset "topological ranked extinction order and cascade" begin
    @test Set(keys(results)) == EXPECTED_SCENARIOS
    @test_throws ArgumentError run_topological_extinctions(network, body_masses[1:6])

    for (scenario, metric, descending) in (
        ("degree_high", :degree, true),
        ("degree_low", :degree, false),
        ("vul_high", :vulnerability, true),
        ("vul_low", :vulnerability, false),
        ("gen_high", :generality, true),
        ("gen_low", :generality, false)
    )
        check_ranked_trajectory(
            results[scenario], A, species_ids, consumer_ids, metric, descending
        )
    end

    # The degree ranking is recalculated after each removal: after A is removed,
    # C becomes highest-degree, ahead of B's original second-place score.
    initial_scores = consumer_scores(
        A, species_ids, species_ids, consumer_ids, :degree
    )
    initial_order = sort(consumer_ids; by=id -> -initial_scores[id])
    high_degree_order = consumer_removal_order(results["degree_high"], consumer_ids)
    @test initial_order[1:2] == ["consumer_a", "consumer_b"]
    @test high_degree_order[1:2] == ["consumer_a", "consumer_c"]

    basal_ids = species_ids[1:3]
    for descending in (true, false)
        result = results[descending ? "bm_high" : "bm_low"]
        current_ids = species_ids
        ranked_consumers = sort(
            consumer_ids;
            by=id -> body_masses[findfirst(==(id), species_ids)],
            rev=descending
        )
        for after in result.networks[2:end]
            after_ids = web_species(after)
            lost = setdiff(Set(current_ids), Set(after_ids))
            expected_primary = first(
                id for id in ranked_consumers if id in current_ids
            )
            @test expected_primary in lost
            current_ids = after_ids
        end
        @test issubset(Set(basal_ids), Set(web_species(result.networks[2])))
    end

    @test consumer_removal_order(results["bm_high"], consumer_ids) ==
          ["consumer_d", "consumer_c", "consumer_b", "consumer_a"]
    @test [web_species(web) for web in results["bm_low"].networks] == [
        species_ids,
        ["basal_p", "basal_q", "basal_r", "consumer_b", "consumer_c", "consumer_d"],
        ["basal_p", "basal_q", "basal_r", "consumer_c", "consumer_d"],
        ["basal_p", "consumer_d"],
        String[]
    ]
end

@testset "topological random roles and disconnected removal" begin
    random_A, random_ids = random_fixture()
    Random.seed!(19)
    random_results = run_topological_extinctions(
        build_network(random_A, random_ids),
        ones(length(random_ids))
    )
    random_basals = Set(random_ids[1:3])
    random_consumers = Set(random_ids[4:5])

    basal_series = random_results["rand_basal"].networks
    @test length(basal_series) == 4
    for (step, (before, after)) in enumerate(zip(basal_series[1:3], basal_series[2:4]))
        lost = setdiff(Set(web_species(before)), Set(web_species(after)))
        @test length(intersect(lost, random_basals)) == 1
        if step < 3
            @test isempty(intersect(lost, random_consumers))
        else
            @test issubset(random_consumers, lost)
        end
    end
    @test isempty(web_species(last(basal_series)))

    consumer_series = random_results["rand_consumer"].networks
    @test length(consumer_series) == 3
    first_consumer_loss = setdiff(
        Set(web_species(consumer_series[1])), Set(web_species(consumer_series[2]))
    )
    final_loss = setdiff(
        Set(web_species(consumer_series[2])), Set(web_species(consumer_series[3]))
    )
    @test length(intersect(first_consumer_loss, random_consumers)) == 1
    @test isempty(intersect(first_consumer_loss, random_basals))
    @test length(intersect(final_loss, random_consumers)) == 1
    @test issubset(Set(random_ids[1:3]), final_loss)
    @test isempty(web_species(last(consumer_series)))

    isolated_A = falses(8, 8) |> Matrix
    isolated_A[1:7, 1:7] .= A
    isolated_ids = [species_ids; "isolated_species"]
    isolated_results = run_topological_extinctions(
        build_network(isolated_A, isolated_ids),
        [body_masses; 0.5]
    )
    @test "isolated_species" in web_species(first(isolated_results["degree_high"].networks))
    @test !("isolated_species" in web_species(isolated_results["degree_high"].networks[2]))
end

println("Validated ranked order, cascade behaviour, random roles, and isolated-node removal.")
