################################################################################
# SCORE PRE EXTINCTION NETWORKS AGAINST THE REALISED POST COMMUNITY
#
# Reads networks and their explicit species-ID ordering from outputs/tss_networks.
# CSV scores are saved under outputs/tss/ and to outputs/tss_results.csv for R.
################################################################################

using CSV
using DataFrames
using EcologicalNetworksDynamics
using JLD2
using SpeciesInteractionNetworks
using Statistics

include("libs/tss.jl")

function mean_or_missing(values)
    observed = collect(skipmissing(values))
    return isempty(observed) ? missing : mean(observed)
end

function sd_or_missing(values)
    observed = collect(skipmissing(values))
    return length(observed) < 2 ? missing : std(observed)
end

pre_networks = JLD2.load_object("outputs/pre/realised_network_archive.jld2")["networks"]
post_networks = JLD2.load_object("outputs/post/realised_network_archive.jld2")["networks"]
simulations = JLD2.load_object("outputs/tss_networks/pre_at_post_richness.jld2")

rows = Dict[]
for (key, simulation) in simulations
    community, net_id, net_name, extinction_time, extinction_type, scenario = key
    network_key = (community, net_id, net_name)
    haskey(pre_networks, network_key) ||
        error("Missing matching PRE network: $network_key")
    haskey(post_networks, network_key) ||
        error("Missing matching POST network: $network_key")

    pre = pre_networks[network_key]
    post = post_networks[network_key]
    scores = calculate_tss(pre, post, simulation)

    push!(rows, merge(Dict(
        :community => community,
        :net_id => net_id,
        :net_name => net_name,
        :extinction_time => extinction_time,
        :extinction_type => extinction_type,
        :scenario => scenario,
        :target_S => simulation["target_S"],
        :actual_S => simulation["actual_S"]
    ), Dict(pairs(scores))))
end

scores = DataFrame(rows)
CSV.write("outputs/tss_scores.csv", scores)
JLD2.save_object("outputs/tss_networks/tss_scores.jld2", scores)
CSV.write("outputs/tss_results.csv", scores)