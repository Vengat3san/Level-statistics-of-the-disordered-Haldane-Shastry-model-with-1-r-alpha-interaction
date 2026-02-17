using Base.Threads,JLD2

include("../src/Hamiltonians.jl")
include("../src/Level_statistics.jl")
include("../src/helper_functions.jl")


function compute_and_save_level_stat_alpha_range(L::Int64, m::Float64, alpha_list, hj::Float64,delta::Float64, N::Int64)

    m_states = m_block(L, m)

    Level_stat_result = Dict{Tuple{Int64, Float64, Float64, Float64, Float64}, NamedTuple}()

    for alpha in alpha_list
        println("Computing for alpha = $alpha")
        r_mean_list = zeros(N)
        s_mean_list = zeros(N)
        s_var_list  = zeros(N)
        s_pdf_list = Vector{StatsBase.Histogram}(undef, N)


        Threads.@threads for i in 1:N
            H = z_h_disordered_H_S_Hamiltonian_m_block(L,m_states,alpha,delta,hj)
            E = eigvals!(H)

            r_mean = level_stat_ratio_test(E)[2]
            s_pdf, s_mean, s_var  = level_stat_distribution_test(E,0.05)

            r_mean_list[i] = r_mean
            s_mean_list[i] = s_mean
            s_var_list[i]  = s_var
            s_pdf_list[i]  = s_pdf 
        end

        bin_list,avg_s_pdf = average_histograms(s_pdf_list)
        r_mean_avg = Statistics.mean(r_mean_list)
        s_mean_avg = Statistics.mean(s_mean_list)
        s_var_avg  = Statistics.mean(s_var_list)

        Level_stat_result[(L, alpha, m, hj,delta)] = (
            r_mean_list = r_mean_list,
            s_mean_list = s_mean_list,
            s_var_list  = s_var_list,
            r_mean_avg  = r_mean_avg,
            s_mean_avg  = s_mean_avg,
            s_var_avg   = s_var_avg,
            avg_s_pdf = avg_s_pdf,
            bin_list = bin_list
            )

             results_dir = joinpath("results","Level_statistics_z_h_dis_HS_gen_alpha")
             isdir(results_dir) || mkpath(results_dir)
             filename = joinpath(results_dir,"Level_Statistics_L$(L)_m$(m)_alpha$(alpha)_delta_$(delta)_h_by_J$(hj)_ensemble_$(N).jld2")
             @save filename Level_stat_result
             println("Saved all results to $filename")

    end


 L=16
 m=0.0 
 alpha_start=length(ARGS) >= 1 ? parse(Float64, ARGS[1])    : 0.2
 alpha_end=length(ARGS) >= 2 ? parse(Float64, ARGS[2])    : 1.0

 alpha_list=alpha_start:0.2:alpha_end

 delta=length(ARGS) >= 3 ? parse(Float64, ARGS[3])    : 0.1
 hj=length(ARGS) >= 4 ? parse(Float64, ARGS[4])    : 0.1 

 N=1024

 println("Using $(nthreads()) threads ")
 println("Using $(MKL.get_num_threads()) MKL threads")
 println(" ")
 println("Computing for L=$L, m=$m, N=$N disorder realizations for alpha in $alpha_list and delta in $delta and h by J $hj")
 println(" ")

 compute_and_save_level_stat_alpha_range(L, m, alpha_list, hj, delta, N)
