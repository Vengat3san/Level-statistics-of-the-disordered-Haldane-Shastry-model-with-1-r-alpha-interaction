
using Base.Threads,JLD2

include("../src/Hamiltonians.jl")
include("../src/Level_statistics.jl")
include("../src/helper_functions.jl")


function compute_and_save_all(L::Int; alpha_range=0.0:0.1:3.0, output_prefix="HS_gen_energies")
    mk_basis_dict = precompute_all_mk_basis(L)

    Threads.@threads for alpha in alpha_range

        println("Processing α = $alpha for L = $L")

        alphaf = Float64(alpha)
        dij = dij_matrix(L, alphaf)
        energies_by_mk = Dict{Tuple{Float64,Float64}, Vector{Float64}}()

        for ((m, k), (k_basis_rep, k_basis_period)) in mk_basis_dict
            
            if isempty(k_basis_rep)
                continue
            end

            m_k_states = (k_basis_rep, k_basis_period)
            H = H_S_Hamiltonian_m_k_block(L, m_k_states, dij, k)
            if size(H, 1) == 0
                continue
            end

            # Compute eigenvalues and eigenvectors
            e_vals= eigvals!(H)

            # Warn if imaginary part is not negligible
            for ev in e_vals
                if abs(imag(ev)) >= 1e-3
                    println("⚠ Warning: Imaginary part of eigenvalue ≥ 1e-3 for L=$L, α=$alphaf, sector (m=$m, k=$k)")
                    break 
                end
            end

            # Store only the real parts (sorted)
            energies_by_mk[(m, k)] = sort(real.(e_vals))
        end

        # Save results
        results_dir = joinpath( "results", "HS_gen_energies")
        isdir(results_dir) || mkpath(results_dir)
        filename = joinpath(results_dir,
            "$(output_prefix)_L$(L)_alpha$(round(alphaf, digits=2)).jld2")
            
        @save filename energies_by_mk
    end
end

 
L           = length(ARGS) >= 1 ? parse(Int, ARGS[1])    : 10
alpha_start = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 0.0
alpha_step  = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : 0.1
alpha_end   = length(ARGS) >= 4 ? parse(Float64, ARGS[4]) : 3.0

alpha_range = alpha_start:alpha_step:alpha_end
output_prefix = "HS_gen_energies"

println("Starting computation for L=$L, α in [$alpha_start, $alpha_end] with step $alpha_step")

compute_and_save_all(L; alpha_range=alpha_range, output_prefix=output_prefix)

println("Computation completed for L=$L, α in [$alpha_start, $alpha_end] with step $alpha_step")