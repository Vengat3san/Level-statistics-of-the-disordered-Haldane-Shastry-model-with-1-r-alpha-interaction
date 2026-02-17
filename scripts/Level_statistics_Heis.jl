# using Base.Threads,JLD2

# include("../src/Hamiltonians.jl")
# include("../src/Level_statistics.jl")
# include("../src/helper_functions.jl")

# function compute_and_save_level_stat_alpha_range(L::Int, m::Float64, alpha_list, delta::Float64, N::Int)

#     m_states = m_block(L, m)

#     Level_stat_result = Dict{Tuple{Int64, Float64, Float64, Float64}, NamedTuple}()

#     for alpha in alpha_list
#         println("Computing for alpha = $alpha")

#         r_mean_list = zeros(N)
#         # s_mean_list = zeros(N)
#         # s_var_list  = zeros(N)
#         # s_pdf_list = Vector{StatsBase.Histogram}(undef, N)

#         Threads.@threads for i in 1:N
#             H0 = Heis_XXX_rand_Hamiltonian_m_block(L, m_states, alpha, delta)
#             H  = h_disordered_Hamiltonian_m_block(L,H0,m_states,h)
#             E = eigvals!(H)

#             r_mean = level_stat_ratio_test(E)[2]
#             # s_pdf, s_mean, s_var  = level_stat_distribution_test(E,0.05)

#              r_mean_list[i] = r_mean
#             # s_mean_list[i] = s_mean
#             # s_var_list[i]  = s_var
#             # s_pdf_list[i]  = s_pdf 
#         end

#         # bin_list,avg_s_pdf = average_histograms(s_pdf_list)
#         r_mean_avg = Statistics.mean(r_mean_list)
#         # s_mean_avg = Statistics.mean(s_mean_list)
#         # s_var_avg  = Statistics.mean(s_var_list)

#         Level_stat_result[(L, alpha, m, delta)] = (
#             r_mean_list = r_mean_list,
#             # s_mean_list = s_mean_list,
#             # s_var_list  = s_var_list,
#             r_mean_avg  = r_mean_avg,
#             # s_mean_avg  = s_mean_avg,
#             # s_var_avg   = s_var_avg,
#             # avg_s_pdf = avg_s_pdf,
#             # bin_list = bin_list
#             )
#     end


#     results_dir = joinpath("results","Heis_Level_statistics_alpha")
#     isdir(results_dir) || mkpath(results_dir)
#     filename = joinpath(results_dir,"Heis_Level_Statistics_L$(L)_m$(m)_alpha0.2_0.2_3.0_delta$(delta)_ensemble_$(N)_all_alphas.jld2")
#     @save filename Level_stat_result
#     println("Saved all results to $filename")
# end


# L=12
# m=0.0
# alpha_list=0.2:0.2:3.0
# delta_list=[0.1,0.5,1.0]
# N=100
# for delta in delta_list
#     compute_and_save_level_stat_alpha_range(L, m, alpha_list, delta, N)
# end

using Base.Threads,JLD2

include("../src/Hamiltonians.jl")
include("../src/Level_statistics.jl")
include("../src/helper_functions.jl")


function compute_and_save_level_stat_alpha_range(L::Int64, m::Float64, alpha_list, hj::Float64,delta::Float64, N::Int64)

    m_states = m_block(L, m)

    Level_stat_result = Dict{Tuple{Int64, Float64, Float64, Float64, Float64}, NamedTuple}()

    for alpha in alpha_list
        println("Computing for alpha = $alpha")
        h=hj*((pi/L)/sin(pi/L))^alpha
        r_mean_list = zeros(N)
        # s_mean_list = zeros(N)
        # s_var_list  = zeros(N)
        # s_pdf_list = Vector{StatsBase.Histogram}(undef, N)


        Threads.@threads for i in 1:N
            H0 = Heis_XXX_rand_Hamiltonian_m_block(L, m_states, alpha, delta)
            H  = h_disordered_Hamiltonian_m_block(L,H0,m_states,h)
            E = eigvals!(H)

            r_mean = level_stat_ratio_test(E)[2]
            # s_pdf, s_mean, s_var  = level_stat_distribution_test(E,0.05)

            r_mean_list[i] = r_mean
        #     s_mean_list[i] = s_mean
        #     s_var_list[i]  = s_var
        #     s_pdf_list[i]  = s_pdf 
        end

        # bin_list,avg_s_pdf = average_histograms(s_pdf_list)
        r_mean_avg = Statistics.mean(r_mean_list)
        # s_mean_avg = Statistics.mean(s_mean_list)
        # s_var_avg  = Statistics.mean(s_var_list)

        Level_stat_result[(L, alpha, m, hj,delta)] = (
            r_mean_list = r_mean_list,
            # s_mean_list = s_mean_list,
            # s_var_list  = s_var_list,
            r_mean_avg  = r_mean_avg,
            # s_mean_avg  = s_mean_avg,
            # s_var_avg   = s_var_avg,
            # avg_s_pdf = avg_s_pdf,
            # bin_list = bin_list
             )
    end


    results_dir = joinpath("results","Heis_Level_statistics_alpha")
    isdir(results_dir) || mkpath(results_dir)
    filename = joinpath(results_dir,"Heis_Level_Statistics_L$(L)_m$(m)_alpha0.2_0.2_3.0_delta$(delta)_hj$(hj)_ensemble_$(N)_all_alphas.jld2")
    @save filename Level_stat_result
    println("Saved all results to $filename")

end


L=12
#m=0.0
m=length(ARGS) >= 1 ? parse(Float64, ARGS[1])    : 0.0
alpha_list=0.2:0.2:3.0
delta_list=[0.0]#0.1,0.5,1.0]#,7.0,10.0,14.0]
#hj=length(ARGS) >= 2 ? parse(Float64, ARGS[2])    : 0.1
N=100

hj_list=[0.1,0.5]
println("Using $(nthreads()) threads ")
println(" ")
println("Computing for L=$L, m=$m, N=$N disorder realizations for alpha in $alpha_list and delta in $delta_list hj is $hj_list")
println(" ")

for hj in hj_list
    println("Computing for hj=$hj")
    for delta in delta_list

        println("Computing for delta=$delta")
        compute_and_save_level_stat_alpha_range(L, m, alpha_list, hj, delta, N)
    end
end