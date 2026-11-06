################################################################################
# EXPORT NETWORK AND SPECIES-LEVEL TOPOLOGY
#
# Species IDs in every row follow the adjacency matrix ordering.
################################################################################

using CSV
using DataFrames
using JLD2
using SpeciesInteractionNetworks

include("libs/topology.jl")
include("libs/internals.jl")

networks = JLD2.load_object("outputs/pre/adjacency_matrices.jld2")

topology_rows = Dict[]
degree_rows = DataFrame(
    spp_id=String[],
    gen=Int[],
    vul=Int[],
    net_id=Int[],
    net_type=String[],
    community=String[],
    stage=String[]
)

for row in eachrow(networks)
    species_ids = String.(row.species)
    network = build_network(Matrix(row.adjacency), species_ids)

    summary = _network_summary(network)
    merge!(summary, Dict(
        :net_id => row.net_id,
        :net_type => row.net_type,
        :community => row.community,
        :stage => row.stage
    ))
    push!(topology_rows, summary)

    generality = SpeciesInteractionNetworks.generality(network)
    vulnerability = SpeciesInteractionNetworks.vulnerability(network)
    append!(degree_rows, DataFrame(
        spp_id=species_ids,
        gen=[generality[Symbol(id)] for id in species_ids],
        vul=[vulnerability[Symbol(id)] for id in species_ids],
        net_id=fill(row.net_id, length(species_ids)),
        net_type=fill(row.net_type, length(species_ids)),
        community=fill(row.community, length(species_ids)),
        stage=fill(row.stage, length(species_ids))
    ))
end

CSV.write("outputs/paleo_topology.csv", DataFrame(topology_rows))
CSV.write("outputs/paleo_spp_degrees.csv", degree_rows)