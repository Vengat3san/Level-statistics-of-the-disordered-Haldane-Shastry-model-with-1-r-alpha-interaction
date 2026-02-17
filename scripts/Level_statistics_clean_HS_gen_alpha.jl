using Base.Threads,JLD2

include("../src/Hamiltonians.jl")
include("../src/Level_statistics.jl")
include("../src/helper_functions.jl")


function compute_and_save_level_stat_alpha_range(L::Int, m::Float64,k::Float64, alpha_list)

    k_basis_rep, k_basis_period = m_k_states_period(L, k, m)
    m_k_states = (k_basis_rep, k_basis_period)
    Level_stat_result = Dict{Tuple{Int64, Float64, Float64, Float64}, NamedTuple}()

    Threads.@threads for alpha in alpha_list
        println("Computing for alpha = $alpha")
        dij = dij_matrix(L, alpha)
    
        H = H_S_Hamiltonian_m_k_block(L, m_k_states, dij, k)
        E = eigvals!(H)
        r_mean = level_stat_ratio_test(E)[2]
        s_pdf, s_mean, s_var  = level_stat_distribution_test(E,0.05)
        avg_s_pdf = s_pdf
        r_mean_avg = r_mean
        s_mean_avg = s_mean
        s_var_avg  = s_var
        edges = 0:0.05:5
        bin_list = 0.5 .* (edges[1:end-1] .+ edges[2:end])

        Level_stat_result[(L, alpha, m, k)] = (
            r_mean_avg  = r_mean_avg,
            s_mean_avg  = s_mean_avg,
            s_var_avg   = s_var_avg,
            avg_s_pdf = avg_s_pdf,
            bin_list = bin_list
            )
    end


    results_dir = joinpath("results","Level_statistics_clean_HS_gen_alpha")
    isdir(results_dir) || mkpath(results_dir)
    filename = joinpath(results_dir,"Level_Statistics_L$(L)_m$(m)_k$(k)_alpha0.0_0.2_3.0_all_alphas.jld2")
    @save filename Level_stat_result
    println("Saved all results to $filename")
end



L=length(ARGS) >= 1 ? parse(Int64, ARGS[1])    : 10
m=0.0
k=1.0
alpha_list= 0.0:0.2:3.0 

compute_and_save_level_stat_alpha_range(L, m,k, alpha_list)
