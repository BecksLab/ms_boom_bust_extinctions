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
using Test

outputs_dir = normpath(joinpath(@__DIR__, "..", "outputs"))
pre_archive_path = joinpath(outputs_dir, "pre", "realised_network_archive.jld2")
post_archive_path = joinpath(outputs_dir, "post", "realised_network_archive.jld2")
summary_path = joinpath(outputs_dir, "extinction_summary.csv")
curves_path = joinpath(outputs_dir, "extinction_curves.csv")
tss_networks_path = joinpath(outputs_dir, "tss_networks", "pre_at_post_richness.jld2")
tss_species_path = joinpath(outputs_dir, "tss_networks", "species.csv")

for path in (
    pre_archive_path, post_archive_path, summary_path, curves_path,
    tss_networks_path, tss_species_path
)
    isfile(path) || error(
        "Missing $(basename(path)); run 01_buildNetworks.jl and " *
        "02_extinctionSims.jl before this test."
    )
end

pre_networks = JLD2.load_object(pre_archive_path)["networks"]
post_networks = JLD2.load_object(post_archive_path)["networks"]
summary = CSV.read(summary_path, DataFrame)
curves = CSV.read(curves_path, DataFrame)
tss_networks = JLD2.load_object(tss_networks_path)
tss_species = CSV.read(tss_species_path, DataFrame)

network_key(key) = (String(key[1]), Int(key[2]), String(key[3]))
pre_keys = Set(network_key(key) for key in keys(pre_networks))
post_keys = Set(network_key(key) for key in keys(post_networks))
summary_keys = [
    (String(row.community), Int(row.net_id), String(row.net_type))
    for row in eachrow(summary)
]
summary_by_key = Dict(
    key => row
    for (key, row) in zip(summary_keys, eachrow(summary))
)

expected_scenarios = Set([
    "degree_high", "degree_low",
    "vul_high", "vul_low",
    "gen_high", "gen_low",
    "bm_high", "bm_low",
    "rand_basal", "rand_consumer"
])
curve_prefix = Dict(
    "topological_creation" => "topo_creation_",
    "topological_realised" => "topo_realised_",
    "dynamic" => "dyn_realised_"
)
expected_curve_groups = Set(
    (curve_type, scenario)
    for curve_type in keys(curve_prefix)
    for scenario in expected_scenarios
)
expected_tss_groups = Set(
    (extinction_type, scenario)
    for extinction_type in ("topological", "dynamic")
    for scenario in expected_scenarios
)
species_key_columns = [
    :community, :net_id, :net_name, :extinction_time, :extinction_type, :scenario
]
species_groups = Dict(
    (
        String(group.community[1]),
        Int(group.net_id[1]),
        String(group.net_name[1]),
        Int(group.extinction_time[1]),
        String(group.extinction_type[1]),
        String(group.scenario[1])
    ) => group
    for group in groupby(tss_species, species_key_columns)
)
pre_networks_by_key = Dict(network_key(key) => value for (key, value) in pre_networks)
post_networks_by_key = Dict(network_key(key) => value for (key, value) in post_networks)

function check_archive_state(networks, state)
    for (raw_key, network) in networks
        community, net_id, net_type = network_key(raw_key)
        generation = network["generation"]
        species = String.(generation["species"])
        creation = Matrix(network["creation"]["A"])
        realised = network["realised"]
        survivors = Int.(realised["survivors"])
        realised_species = species[survivors]
        realised_adjacency = Matrix(realised["A"])

        @test String(generation["community"]) == community
        @test Int(generation["net_id"]) == net_id
        @test String(generation["state"]) == state
        @test length(unique(species)) == length(species)
        @test size(creation) == (length(species), length(species))
        @test all(i -> 1 <= i <= length(species), survivors)
        @test length(unique(survivors)) == length(survivors)
        @test size(realised_adjacency) ==
              (length(realised_species), length(realised_species))

        if state == "post"
            pre_species = String.(pre_networks[raw_key]["generation"]["species"])
            @test issubset(Set(species), Set(pre_species))
        end
    end
end

function robustness_from_curve(curve, initial_richness)
    ordered = sort(curve, :step)
    richness = round.(Int, initial_richness .* (1 .- ordered.total))
    target = initial_richness / 2

    for i in 2:length(richness)
        if richness[i] <= target
            previous, current = richness[i - 1], richness[i]
            x1 = (i - 2) / initial_richness
            x2 = (i - 1) / initial_richness
            previous == current && return x2
            fraction = (target - previous) / (current - previous)
            return x1 + fraction * (x2 - x1)
        end
    end

    return (length(richness) - 1) / initial_richness
end

@testset "paired network and extinction metadata" begin
    @test pre_keys == post_keys
    @test length(summary_keys) == length(unique(summary_keys))
    @test Set(summary_keys) == pre_keys
    @test Set(String.(curves.type)) == Set(keys(curve_prefix))

    check_archive_state(pre_networks, "pre")
    check_archive_state(post_networks, "post")

    curve_networks = groupby(curves, [:community, :net_id, :net_name])
    curve_network_keys = Set(
        (
            String(group.community[1]),
            Int(group.net_id[1]),
            String(group.net_name[1])
        )
        for group in curve_networks
    )
    @test curve_network_keys == Set(summary_keys)

    for network_curves in curve_networks
        observed_groups = Set(
            (String(row.type), String(row.scenario))
            for row in eachrow(network_curves)
        )
        @test observed_groups == expected_curve_groups
    end

    scenario_curves = groupby(
        curves,
        [:community, :net_id, :net_name, :extinction_t, :type, :scenario]
    )
    for curve in scenario_curves
        key = (
            String(curve.community[1]),
            Int(curve.net_id[1]),
            String(curve.net_name[1])
        )
        summary_row = summary_by_key[key]
        curve_type = String(curve.type[1])
        scenario = String(curve.scenario[1])
        prefix = curve_prefix[curve_type]
        metric = Symbol(prefix * scenario)
        ordered_curve = sort(curve, :step)
        initial_richness = curve_type == "topological_creation" ?
                           Int(summary_row.S_creation) :
                           Int(summary_row.S_realised)

        @test all(!ismissing, curve.step)
        @test Int(curve.extinction_t[1]) == Int(summary_row.extinction_time)
        @test collect(ordered_curve.step) == collect(1:nrow(curve))
        @test first(ordered_curve.total) == 0
        @test all(diff(ordered_curve.total) .>= -1e-12)
        @test metric in propertynames(summary)
        @test isapprox(
            robustness_from_curve(curve, initial_richness),
            Float64(summary_row[metric]);
            atol=1e-10,
            rtol=0
        )
    end

    @test Set(network_key(key) for key in keys(tss_networks)) ⊆ pre_keys
    @test length(tss_networks) == length(pre_keys) * length(expected_tss_groups)

    nonempty_tss_species_keys = Set{Tuple{String,Int,String,Int,String,String}}()
    tss_network_groups = Dict{Tuple{String,Int,String},Set{Tuple{String,String}}}()
    empty_tss_networks = 0

    for (raw_key, simulation) in tss_networks
        key = (
            String(raw_key[1]),
            Int(raw_key[2]),
            String(raw_key[3]),
            Int(raw_key[4]),
            String(raw_key[5]),
            String(raw_key[6])
        )
        network = pre_networks_by_key[key[1:3]]
        post_network = post_networks_by_key[key[1:3]]
        generation_species = String.(network["generation"]["species"])
        simulation_species = String.(simulation["species"])
        adjacency = Matrix(simulation["A"])

        @test key[4] == Int(network["generation"]["burnin_time"])
        @test simulation["target_S"] == post_network["realised"]["S"]
        @test Int(simulation["actual_S"]) == length(simulation_species)
        @test size(adjacency) == (length(simulation_species), length(simulation_species))
        @test length(unique(simulation_species)) == length(simulation_species)
        @test issubset(Set(simulation_species), Set(generation_species))
        @test key[5] in ("topological", "dynamic")

        pair_key = key[1:3]
        push!(get!(tss_network_groups, pair_key, Set{Tuple{String,String}}()), (key[5], key[6]))

        if isempty(simulation_species)
            empty_tss_networks += 1
        else
            push!(nonempty_tss_species_keys, key)
            @test haskey(species_groups, key)
            if haskey(species_groups, key)
                rows = sort(species_groups[key], :position)
                @test all(rows.network_role .== "simulated_pre")
                @test Int.(rows.position) == collect(1:length(simulation_species))
                @test String.(rows.species_id) == simulation_species
            end
        end
    end

    @test Set(keys(species_groups)) == nonempty_tss_species_keys
    @test all(groups == expected_tss_groups for groups in values(tss_network_groups))
    println("Validated $empty_tss_networks empty selected TSS networks.")
end

println(
    "Validated $(length(pre_keys)) paired networks and " *
    "$(nrow(curves)) extinction-curve rows, plus $(length(tss_networks)) " *
    "extinction-to-TSS handoffs."
)
