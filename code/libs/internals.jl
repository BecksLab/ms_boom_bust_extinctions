using Random

"""
    build_network(adj_mat::Matrix{Bool}, species_ids)

Construct a `SpeciesInteractionNetwork` from a Boolean adjacency matrix.

Converts a square adjacency matrix into a unipartite, binary food web. When
`species_ids` are omitted, nodes are labelled `:1`, ..., `:S`; otherwise the
supplied identifiers are kept in matrix order.

# Arguments
- `adj_mat::Matrix{Bool}`: Square adjacency matrix (S × S), where:
    - `adj_mat[i, j] = true` indicates a trophic interaction (i consumes j)
    - Rows and columns refer to the same ordered species set
- `species_ids`: Optional identifiers in row/column order

# Returns
- `SpeciesInteractionNetwork{Unipartite{Symbol}, Binary{Bool}}`

# Notes
- Species identifiers must be unique and match matrix dimensions
- Interaction direction is preserved
- Assumes adjacency matrix is already ecologically valid

# Errors
- Throws `ArgumentError` if the matrix is not square
"""
function build_network(adj_mat::Matrix{Bool}, species_ids=Symbol.(1:size(adj_mat, 1)))

    if size(adj_mat, 1) != size(adj_mat, 2)
        throw(ArgumentError("Adjacency matrix must be square"))
    end

    spp_rich = size(adj_mat, 1)
    length(species_ids) == spp_rich ||
        throw(ArgumentError("Species identifiers must match adjacency matrix dimensions"))
    spp_list = Symbol.(string.(species_ids))
    length(unique(spp_list)) == spp_rich ||
        throw(ArgumentError("Species identifiers must be unique"))

    edges = Binary(.!iszero.(adj_mat))
    nodes = Unipartite(spp_list)

    return SpeciesInteractionNetwork(nodes, edges)
end

"""
    extract_valid_network(A, B, threshold)

Extract an ecologically valid subnetwork based on biomass and connectivity.

Species are retained if they:
1. Have biomass above a survival threshold
2. Participate in at least one trophic interaction

# Arguments
- `A::Matrix`: Full adjacency matrix
- `B::Vector`: Species biomasses
- `threshold::Float64`: Extinction threshold

# Returns
- `(survivors, A_sub)`:
    - `survivors::Vector{Int}`: Indices of retained species
    - `A_sub::Matrix`: Induced subnetwork adjacency matrix
- `(nothing, nothing)` if no valid species remain

# Ecological Interpretation
- Species with `B ≤ threshold` are considered extinct
- Species with no prey AND no consumers are removed as non-functional
- Ensures resulting network is trophically connected

# Notes
- Mirrors assumptions used in topological extinction workflows
- Does NOT distinguish basal vs consumer roles explicitly
"""
function extract_valid_network(A, B, threshold)

    survivors = findall(x -> x > threshold, B)

    if isempty(survivors)
        return nothing, nothing
    end

    A_sub = A[survivors, survivors]

    has_prey = vec(sum(A_sub; dims=2)) .> 0
    has_consumer = vec(sum(A_sub; dims=1)) .> 0

    keep_mask = has_prey .| has_consumer
    survivors = survivors[keep_mask]

    if isempty(survivors)
        return nothing, nothing
    end

    return survivors, A[survivors, survivors]
end

"""
    run_topological_extinctions(N, params)

Run topological (structural) extinction sequences on a food web.

Species are removed sequentially according to structural traits, and
secondary extinctions occur via trophic disconnection (cascade mechanism).

# Arguments
- `N`: `SpeciesInteractionNetwork`
- `body_mass`: Body mass values in network species order

# Returns
- `Dict{String, NamedTuple}`:
    - `networks`: sequence of networks after each extinction step
    - `primary`: cumulative number of primary deletions

# Extinction Scenarios
- `"degree_high"` / `"degree_low"`: connectivity-based removal
- `"vul_high"` / `"vul_low"`: vulnerability (in-degree)
- `"gen_high"` / `"gen_low"`: generality (out-degree)
- `"bm_high"` / `"bm_low"`: body-mass ordering among consumers; basal species
  are protected from primary removal
- `"rand_basal"`: random removal of basal species only
- `"rand_consumer"`: random removal of consumers only

# Ecological Assumptions
- Secondary extinctions occur when species lose all prey
- Basal species are defined by generality = 0 (in initial network)
- Random removal respects initial trophic roles (fixed classification)

# Important Implementation Detail
- Disconnected species are removed after each primary extinction.
- Consumers that lose all prey are removed through the cascade mechanism.

# Notes
- No population dynamics are simulated
- Fully comparable to classical topological robustness studies
"""
function run_topological_extinctions(N, body_mass)
    length(body_mass) == SpeciesInteractionNetworks.richness(N) ||
        throw(ArgumentError("Body masses must match network richness"))

    results = Dict()

    scenarios = Dict(
        "degree_high" => extinction(
            N, "degree", true; remove_disconnected=true
        ),
        "degree_low" => extinction(
            N, "degree", false; remove_disconnected=true
        ),
        "vul_high" => extinction(
            N, "vulnerability", true; remove_disconnected=true
        ),
        "vul_low" => extinction(
            N, "vulnerability", false; remove_disconnected=true
        ),
        "gen_high" => extinction(
            N, "generality", true; remove_disconnected=true
        ),
        "gen_low" => extinction(
            N, "generality", false; remove_disconnected=true
        ),
        "bm_high" => extinction(
            N, species(N)[sortperm(body_mass, rev=true)];
            protect=:basal,
            remove_disconnected=true
        ),
        "bm_low" => extinction(
            N, species(N)[sortperm(body_mass, rev=false)];
            protect=:basal,
            remove_disconnected=true
        ),
        "rand_basal" => extinction(
            N; protect=:consumer, remove_disconnected=true
        ),
        "rand_consumer" => extinction(
            N; protect=:basal, remove_disconnected=true
        )
    )

    for (name, Ns) in scenarios

        # primary = 1 species removed per step
        primary_counts = collect(0:(length(Ns)-1))

        results[name] = (
            networks=Ns,
            primary=primary_counts
        )
    end

    return results
end

"""
    dynamic_extinction_adaptive(params, B0; kwargs...)

Run adaptive dynamical extinction simulations with sequential species removal.

At each step:
1. Select a target species (based on a criterion)
2. Force its extinction (primary removal)
3. Simulate system dynamics for `t` time steps
4. Record secondary extinctions (biomass collapse)
5. Repeat until system collapse

# Arguments
- `params`: Model parameters (includes adjacency matrix and body masses)
- `B0::Vector`: Initial biomasses after burn-in

# Keyword Arguments
- `criterion::Symbol`: Removal rule
    - `:degree`, `:generality`, `:vulnerability`, `:bodymass`
    - `:random_basal`, `:random_consumer`
- `descending::Bool`: Remove highest (true) or lowest (false)
- `t::Int`: Simulation duration after each deletion
- `survival_threshold::Float64`: Biomass extinction threshold
- `show_progress::Bool`
- `debug::Bool`

# Returns
- NamedTuple:
    - `networks`: sequence of surviving subnetworks
    - `primary`: cumulative number of primary deletions
    - `primary_removed`: species IDs removed as primary targets, in order
    - `secondary_removed`: species IDs falling below threshold after each
      primary deletion, grouped by removal step

# Ecological Interpretation
- Combines press perturbations with relaxation dynamics
- Secondary extinctions emerge from population dynamics
- Captures indirect effects and trophic cascades

# Key Assumptions
- Body masses are FIXED from initial network (no re-sampling)
- Basal/consumer roles for random deletions are FIXED from initial network
- Species go extinct if `B < survival_threshold`

# Notes
- Multiple species may go extinct per step
- Extinction order is adaptive (recomputed each step)
- This corresponds to the dynamical approach in the Curtsdottir 2011
"""
function dynamic_extinction_adaptive(params, B0;
    criterion=:degree,
    descending=true,
    t=5000,
    survival_threshold=1e-12,
    show_progress=true,
    debug=false,
    species_ids=Symbol.(1:length(B0))
)

    S0 = length(B0)
    length(species_ids) == S0 ||
        throw(ArgumentError("Species identifiers must match initial biomass dimensions"))
    length(unique(species_ids)) == S0 ||
        throw(ArgumentError("Species identifiers must be unique"))
    isfinite(survival_threshold) && survival_threshold > 0 ||
        throw(ArgumentError("survival_threshold must be a finite positive number"))

    A = params.A
    B = copy(B0)

    # --- species roles ---
    gen0 = vec(sum(A, dims=2))
    basal0 = gen0 .== 0
    consumer0 = gen0 .> 0

    alive = trues(S0)

    # --- store networks ---
    Nseq = SpeciesInteractionNetwork[]
    push!(Nseq, build_network(Matrix(A), species_ids))

    # Original species indices retained at each extinction step
    survivor_seq = Vector{Vector{Int}}()
    push!(survivor_seq, collect(1:S0))

    primary_counts = Int[0]
    primary_removed = eltype(species_ids)[]
    secondary_removed = Vector{Vector{eltype(species_ids)}}()

    prog = show_progress ? Progress(S0, "Dynamic extinctions") : nothing
    step = 0

    # --- stopping condition ---
    function continue_condition()
        if criterion == :random_basal
            return any(alive .& basal0)
        elseif criterion == :random_consumer
            return any(alive .& consumer0)
        else
            # all trait-based = consumers only
            return any(alive .& consumer0)
        end
    end

    while continue_condition()
        step += 1

        alive_idx = findall(alive)
        isempty(alive_idx) && break

        # --- build subnetwork for ranking ---
        A_sub = A[alive_idx, alive_idx]

        vuln = vec(sum(A_sub, dims=1))
        gen = vec(sum(A_sub, dims=2))
        deg = vuln .+ gen

        # --- SELECT SPECIES ---
        if criterion == :random_basal
            candidates = filter(i -> alive[i] && basal0[i], eachindex(alive))
            isempty(candidates) && break
            sp = rand(candidates)

        elseif criterion == :random_consumer
            candidates = filter(i -> alive[i] && consumer0[i], eachindex(alive))
            isempty(candidates) && break
            sp = rand(candidates)

        else
            # trait-based → consumers only
            candidates = filter(i -> alive[i] && consumer0[i], eachindex(alive))

            # metric evaluated on alive_idx
            metric = if criterion == :degree
                deg
            elseif criterion == :generality
                gen
            elseif criterion == :vulnerability
                vuln
            elseif criterion == :bodymass
                params.body_mass[alive_idx]
            elseif criterion == :random
                nothing
            else
                error("Unknown criterion: $criterion")
            end

            if !isempty(candidates)
                # map candidates to local indices
                local_inds = findall(i -> alive_idx[i] in candidates, eachindex(alive_idx))

                if metric === nothing
                    local_choice = rand(local_inds)
                else
                    vals = metric[local_inds]
                    target = descending ? maximum(vals) : minimum(vals)
                    tied = local_inds[findall(x -> x == target, vals)]
                    local_choice = rand(tied)
                end

                sp = alive_idx[local_choice]
            else
                break
            end
        end

        isempty(candidates) && break

        # --- primary extinction ---
        alive[sp] = false
        push!(primary_removed, species_ids[sp])
        push!(secondary_removed, eltype(species_ids)[])
        push!(primary_counts, primary_counts[end] + 1)

        # --- simulate ALL remaining species ---
        alive_idx = findall(alive)
        isempty(alive_idx) && break

        A_sim = A[alive_idx, alive_idx]
        B_sim = B[alive_idx]
        M_sim = params.M[alive_idx]

        fw = Foodweb(A_sim)

        params_sim = default_model(
            fw,
            BodyMass(M_sim),
            ClassicResponse(;
                h=params.h,
                #c = params.intraspecific_interference[1],
            ),
        )
        disable_producer_growth_for_former_consumers!(
            params_sim,
            params,
            alive_idx,
            basal0
        )

        try
            sol = simulate(params_sim, B_sim, t;
                show_degenerated=false,
                abstol=min(1e-12, survival_threshold / 10),
                callback=CallbackSet(
                    extinction_callback(params_sim, survival_threshold)
                )
            )
            B_sim = sol.u[end]

        catch e
            @warn "Simulation failed" exception=e
            break
        end

        # --- write back ---
        B .= 0.0
        B[alive_idx] = B_sim

        # --- apply extinction threshold ---
        secondary_indices = findall(
            i -> alive[i] && B[i] < survival_threshold,
            eachindex(alive)
        )
        secondary_removed[end] = species_ids[secondary_indices]
        alive[secondary_indices] .= false
        B[secondary_indices] .= 0.0

        # --- record network ---
        survivors = findall(alive)
        isempty(survivors) && break

        A_valid = A[survivors, survivors]
        push!(Nseq, build_network(A_valid, species_ids[survivors]))
        push!(survivor_seq, copy(survivors))

        if show_progress
            next!(prog)
        end

        if debug
            println("Step $step | remaining: $(length(survivors))")
        end
    end

    return (
        networks=Nseq,
        primary=primary_counts,
        survivors=survivor_seq,
        primary_removed=primary_removed,
        secondary_removed=secondary_removed
    )
end

function disable_producer_growth_for_former_consumers!(
    submodel,
    original_model,
    original_indices,
    originally_producers
)
    for local_index in submodel.producers.indices
        original_index = original_indices[local_index]
        if !originally_producers[original_index]
            submodel.growth_rate[local_index] = 0.0
            submodel.metabolism[local_index] = original_model.metabolism[original_index]
        end
    end
end

"""
    run_dynamic_extinctions(params, B; t=5000, show_progress=true)

Run all adaptive dynamic extinction scenarios.

# Arguments
- `params`: Model parameters
- `B`: Initial biomasses after burn-in

# Keyword Arguments
- `t::Int`: Simulation duration after each primary deletion (passed to dynamic simulation)
- `show_progress::Bool`

# Returns
- `Dict{String, NamedTuple}`: Mapping scenario name to its extinction result,
  including the network trajectory and removal traces
"""
function run_dynamic_extinctions(
    params, B;
    t=5000,
    show_progress=true,
    species_ids=Symbol.(1:length(B)),
    seed=66
)

    # EcologicalNetworksDynamics lazily caches model property accessors.
    # Initialize them before worker threads first access the shared model.
    _ = (params.A, params.M, params.body_mass, params.h)

    scenarios = [
        ("degree_high", :degree, true),
        ("degree_low", :degree, false),
        ("vul_high", :vulnerability, true),
        ("vul_low", :vulnerability, false),
        ("gen_high", :generality, true),
        ("gen_low", :generality, false),
        ("bm_high", :bodymass, true),
        ("bm_low", :bodymass, false),
        ("rand_basal", :random_basal, true),
        ("rand_consumer", :random_consumer, true)
    ]

    prog = show_progress ? Progress(length(scenarios), "Dynamic scenarios") : nothing
    progress_lock = ReentrantLock()
    results = Vector{Any}(undef, length(scenarios))

    Threads.@threads for i in eachindex(scenarios)
        name, criterion, descending = scenarios[i]
        Random.seed!(seed + i)

        try
            results[i] = dynamic_extinction_adaptive(
                params, B;
                criterion=criterion,
                descending=descending,
                t=t,
                show_progress=false,
                species_ids=species_ids
            )
        finally
            if show_progress
                lock(progress_lock) do
                    next!(prog)
                end
            end
        end
    end

    return Dict(scenarios[i][1] => results[i] for i in eachindex(scenarios))
end

"""
    compute_robustness(results_dict)

Compute robustness metrics for a set of extinction sequences.

# Arguments
- `results_dict::Dict`: Mapping scenario name → network sequence

# Returns
- `Dict{String, Float64}`: Mapping scenario name → robustness value (e.g. R50)

# Notes
- Assumes all sequences are compatible with `robustness()`
- Works for both topological and dynamic extinction outputs
"""
function compute_robustness(results_dict)

    R = Dict()

    for (k, v) in results_dict
        # extract only networks
        R[k] = robustness_integral(v.networks)
    end

    return R
end

"""
    extinction_breakdown(result)

Decompose extinction trajectories into primary and secondary components.

# Arguments
- `result`: NamedTuple with:
    - `networks`: sequence of networks
    - `primary`: cumulative primary deletions

# Returns
- NamedTuple:
    - `primary`: proportion of primary removals (PDEL)
    - `secondary`: proportion of secondary extinctions
    - `total`: total extinction proportion

# Definitions
- Primary extinction: directly removed species
- Secondary extinction: species lost due to cascading effects
- Total extinction = primary + secondary

# Notes
- All values are normalised by initial richness (S0)
- Primary is constrained to not exceed total extinction
- Ensures physically consistent decomposition

# Warnings
- Emits warnings if extinction curve is non-monotonic
"""
function extinction_breakdown(result)

    Ns = result.networks
    primary_counts = result.primary

    if isempty(Ns)
        return NaN
    end

    n = min(length(Ns), length(primary_counts))

    Ns = Ns[1:n]
    primary_counts = primary_counts[1:n]

    # --- richness ---
    richness_vals = SpeciesInteractionNetworks.richness.(Ns)

    if any(diff(richness_vals) .> 0)
        @warn "Non-monotonic richness detected (INVALID CURVE)"
    end

    S0 = richness_vals[1]

    if S0 <= 0
        error("Invalid S0 (zero or negative richness baseline)")
    end

    if any(richness_vals .> S0)
        @warn "Richness exceeds baseline S0 → normalization broken"
    end

    total_loss = (S0 .- richness_vals) ./ S0

    # --- primary loss (event-based) ---
    primary_prop = primary_counts ./ S0

    # enforce physical constraint:
    # primary cannot exceed observed loss
    primary_prop = min.(primary_prop, total_loss)

    # secondary becomes true residual
    secondary_prop = total_loss .- primary_prop

    if any(diff(total_loss) .< 0)
        @warn "Non-monotonic extinction curve detected"
    end

    return (
        primary=primary_prop,
        secondary=secondary_prop,
        total=total_loss
    )
end

function export_curves(curves_dict, type_label, net_id)

    rows = DataFrame()

    for (scenario, c) in curves_dict

        if !hasproperty(c, :total)

            rows = DataFrame(
                net_id=Int[],
                type=String[],
                scenario=String[],
                step=Union{Missing,Int}[],
                primary=Union{Missing,Float64}[],
                secondary=Union{Missing,Float64}[],
                total=Union{Missing,Float64}[]
            )

        else
            n = length(c.total)

            append!(rows, DataFrame(
                net_id=fill(net_id, n),
                type=fill(type_label, n),
                scenario=fill(scenario, n),
                step=1:n,
                primary=c.primary,
                secondary=c.secondary,
                total=c.total
            ))
        end
    end

    return rows
end

function add_curve_metadata!(df, community, net_name, extinction_t)
    df[!, :community] = fill(community, nrow(df))
    df[!, :net_name] = fill(net_name, nrow(df))
    df[!, :extinction_t] = fill(extinction_t, nrow(df))

    return df
end

"""
    robustness_integral(network_sequence)

Compute robustness (R50) from an extinction sequence.

R50 is defined as the proportion of species that must be removed
(primary deletions) to cause 50% total species loss.

# Arguments
- `network_sequence`: Vector of `SpeciesInteractionNetwork`

# Returns
- `Float64`: R50 value (∈ [0, 0.5])

# Method
- Tracks species richness decline across extinction steps
- Identifies when richness falls below 50% of initial richness (S0)
- Uses linear interpolation between steps to estimate R50

# Interpretation
- R50 = 0.5 → no secondary extinctions (max robustness)
- R50 < 0.5 → cascading extinctions reduce robustness

# Notes
- Assumes 1 primary deletion per step
- Uses step index as proxy for PDEL (as in topological approaches)
"""
function robustness_integral(network_sequence::Vector{<:SpeciesInteractionNetwork})

    if isempty(network_sequence)
        return NaN
    end

    S = SpeciesInteractionNetworks.richness.(network_sequence)
    S0 = S[1]

    target = 0.5 * S0

    for i in 2:length(S)

        if S[i] <= target

            S1, S2 = S[i-1], S[i]

            # PDEL values
            x1 = (i - 2) / S0
            x2 = (i - 1) / S0

            if S1 == S2
                return x2
            end

            # linear interpolation
            frac = (target - S1) / (S2 - S1)

            R50 = x1 + frac * (x2 - x1)

            return R50
        end
    end

    return (length(S) - 1) / S0

end

function robustness_auc(network_sequence)
    S = SpeciesInteractionNetworks.richness.(network_sequence)
    S0 = S[1]

    frac = S ./ S0
    x = (0:(length(S)-1)) ./ S0

    return sum((frac[1:(end-1)] .+ frac[2:end]) ./ 2 .* diff(x))
end

function compute_trophic_levels(A)

    S = size(A, 1) # Species richness.
    out_degree = sum(A; dims=2)
    D = -(A ./ out_degree) # Diet matrix.
    D[isnan.(D)] .= 0.0
    D[diagind(D)] .= 1.0 .- D[diagind(D)]
    # Solve with the inverse matrix.
    inverse = iszero(det(D)) ? pinv : inv
    tls = inverse(D) * ones(S)

    return tls
end

function realise_network(
    fw;
    bodymasses=nothing,
    t=5000,
    threshold=1e-12,
    B0
)

    if isnothing(bodymasses)

        params = default_model(
            fw,
            BodyMass(; Z=10),
            ClassicResponse(; h=2.0),
        )

    else

        params = default_model(
            fw,
            BodyMass(bodymasses),
            ClassicResponse(; h=2.0),
        )

    end

    S = size(fw.A, 1)

    sol = simulate(
        params,
        B0,
        t;
        show_degenerated=false,
        callback=CallbackSet(
            extinction_callback(params, threshold),
            TerminateSteadyState(1e-14, 1e-12, DiffEqCallbacks.allDerivPass),
        )
    )

    final_biomasses = sol.u[end]

    survivors, final_adj =
        extract_valid_network(
            params.A,
            final_biomasses,
            threshold
        )

    survivors === nothing && return nothing

    S_real = length(survivors)
    L_real = sum(final_adj)
    C_real = L_real / S_real^2

    return (
        A=final_adj,
        params=params,
        biomasses=final_biomasses,
        survivors=survivors,
        S=S_real,
        L=L_real,
        C=C_real
    )
end

function get_trophic_class(A)
    # vul: number of predators (summing columns, i.e., across rows for each column)
    vul = vec(sum(A, dims=1))
    # gen: number of prey (summing rows, i.e., across columns for each row)
    gen = vec(sum(A, dims=2))

    # Define the classification logic for a single species
    classify(g, v) =
        if g == 0 && v > 0
            :basal
        elseif v == 0 && g > 0
            :top
        elseif v > 0 && g > 0
            :intermediate
        else
            :isolated
        end

    # Broadcast the classify function over the gen and vul vectors
    return classify.(gen, vul)
end

# --- Helper Functions for Clean Extraction ---

"""
Helper to extract, package, and append species level metadata (TL & TC) 
at different lifecycle stages.
"""
function record_species_stage!(
    store, run_id, net_name, stage, A, bodysizes,
    params=nothing, survivors=nothing, biomasses=nothing;
    species_ids=1:size(A, 1)
)
    S = size(A, 1)
    indices = something(survivors, 1:S)
    length(species_ids) >= maximum(indices; init=0) ||
        throw(ArgumentError("Species identifiers do not cover the recorded indices"))

    # Safely get Trophic Levels if parameters are available
    TL = (params !== nothing && hasproperty(params, :trophic)) ? params.trophic.levels[indices] : trophic_level(Matrix(A))

    # Biomasses - if provided. if both biomass and survivors then subset by indices else full 'list'
    BM = (biomasses !== nothing && survivors !== nothing) ? biomasses[indices] : biomasses

    #bodysizes - if we have survivors then we need to subset by indices else full 'list'
    BS = (survivors !== nothing) ? params.M[indices] : params.M

    #metabolism - subset by survivors but then also return nothing if not assigned
    MB = (params !== nothing && hasproperty(params, :metabolism)) ? params.metabolism[indices] : fill(:nothing, length(indices))

    # Get Trophic Classes for the active matrix block
    TC = get_trophic_class(A)

    # Package into a temporary DataFrame and append to store
    df = DataFrame(
        net_id=fill(run_id, length(indices)),
        net_type=fill(net_name, length(indices)),
        stage=fill(stage, length(indices)),
        species_id=1:length(indices),
        original_id=indices,
        species=String.(species_ids[indices]),
        trophic_level=TL,
        trophic_class=TC,
        biomass=BM,
        bodysize=BS,
        metabolism=MB
    )
    append!(store, df)
end

"""
Helper to find best threshold value that returns an lmatrix network matching 
the target connectance
"""
function atn_build_range(df, prods;
    thresholds=0.01:0.01:0.12)

    results = NamedTuple[]

    for thresh in thresholds
        adj = lmatrix(df.species, df.bodymass, prods;
            threshold=thresh)

        S = size(adj, 1)
        L = sum(adj)
        Co = L / S^2

        push!(results, (
            threshold=thresh,
            connectance=Co,
            adjacency=adj
        ))
    end

    return results
end