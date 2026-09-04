# for dimulating bodysizes following the appraoch of Stouffer and Bascompte 2010
# doi: /10.1111/j.1461-0248.2009.01407.x

using Distributions

function initial_bodymasses(traits, size_bounds, global_dist)
    [
        begin
            lo, hi = size_bounds[String(row.size)]
            rand(truncated(global_dist, lo, hi))
        end
        for row in eachrow(traits)
    ]
end

function loglikelihood(logM, A; μ=6.1, σ=5.75)

    ll = 0.0

    for predator in axes(A,1), prey in axes(A,2)

        if A[predator,prey] == 1
            ratio = logM[predator] - logM[prey]
            ll += logpdf(Normal(μ,σ), ratio)
        end

    end

    return ll
end

function propose!(logM, idx, bounds)

    proposal = copy(logM)

    proposal[idx] += rand(Normal(0,0.25))

    lo, hi = log.(bounds)

    proposal[idx] = clamp(proposal[idx], lo, hi)

    return proposal
end

function metropolis_step(logM, A, species_bounds)

    current = loglikelihood(logM, A)

    idx = rand(eachindex(logM))

    proposal = propose!(logM, idx, species_bounds[idx])

    proposed = loglikelihood(proposal, A)

    if log(rand()) < proposed - current
        return proposal
    else
        return logM
    end
end

function nested_bodymasses(traits, A, size_bounds;
                           μ=6.1,
                           σ=5.75,
                           iterations=20000)

    M = initial_bodymasses(traits, size_bounds, global_dist)

    logM = log.(M)

    species_bounds = [
        size_bounds[String(s)]
        for s in traits.size
    ]

    for _ in 1:iterations
        logM = metropolis_step(logM, A, species_bounds)
    end

    return exp.(logM)
end