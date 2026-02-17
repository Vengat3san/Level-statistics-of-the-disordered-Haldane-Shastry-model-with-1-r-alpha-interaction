
using Dierckx,Statistics,StatsBase


function level_stat_ratio_test(E::Vector{<:Real})
    tol = 1e-10
    dE = diff(sort(E))
    Sn = dE[2:end]
    Sn1 = dE[1:end-1]

    r_tilde = [min(s, s1) / max(s, s1) for (s, s1) in zip(Sn, Sn1) if max(abs(s), abs(s1)) > tol]
    return r_tilde, Statistics.mean(r_tilde), Statistics.var(r_tilde)
end


function cumulative_density(E::Vector{<:Real}; tol=1e-10)
    unique_E = Float64[]
    nE = Int[]
    count = 0
    i = 1
    N = length(E)

        while i <= N
            current_val = E[i]
            degenerate_count = 1

            while i < N && abs(E[i+1] - current_val) < tol
                degenerate_count += 1
                i += 1
            end

            count += degenerate_count
            push!(unique_E, current_val)
            push!(nE, count)
            i += 1
        end

        return unique_E, nE
end

function level_stat_distribution_test(E)

    unique_E, Cumulative_count = cumulative_density(E; tol=1e-10)
    if length(unique_E) < 35
        @warn "Not enough unique energy levels for distribution test."
        return [], 0.0, 0.0
    end
    n_sparse = min(round(length(unique_E)/10), length(unique_E))
    indices = round.(Int, range(1, length(unique_E), length=Int(n_sparse)))
    E_sparse = unique_E[indices]
    nE_sparse = Float64.(Cumulative_count[indices])

    spline = Spline1D(E_sparse, nE_sparse, k=3)  

    x_unfolded = spline.(E)

    x_unfolded_sorted = sort(x_unfolded)
    s = diff(x_unfolded_sorted)

    return s, Statistics.mean(s), Statistics.var(s)
end

function level_stat_distribution_test(E,bin_width::Float64)

    unique_E, Cumulative_count = cumulative_density(E; tol=1e-10)
    if length(unique_E) < 35
        @warn "Not enough unique energy levels for distribution test."
        return [], 0.0, 0.0
    end
    n_sparse = min(round(length(unique_E)/10), length(unique_E))
    indices = round.(Int, range(1, length(unique_E), length=Int(n_sparse)))
    E_sparse = unique_E[indices]
    nE_sparse = Float64.(Cumulative_count[indices])

    spline = Spline1D(E_sparse, nE_sparse, k=3)  

    x_unfolded = spline.(E)

    x_unfolded_sorted = sort(x_unfolded)
    s = diff(x_unfolded_sorted)


    edges = 0:bin_width:5
    s_hist = fit(StatsBase.Histogram, s, edges)
    s_pdf = normalize(s_hist, mode=:pdf)

    return s_pdf, Statistics.mean(s), Statistics.var(s)
end
