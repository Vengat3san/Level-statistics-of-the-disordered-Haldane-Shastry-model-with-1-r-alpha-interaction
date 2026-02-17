using Base.Threads,JLD2

include("../src/Hamiltonians.jl")
include("../src/Level_statistics.jl")
include("../src/helper_functions.jl")

function compute_and_save_all(L::Int,J,Delta; output_prefix="XXZ_Heis_energies")
        mk_basis_dict = precompute_all_mk_basis(L)

        energies_by_mk = Dict()

        for ((m, k), (k_basis_rep, k_basis_period)) in mk_basis_dict
            println("Computing for sector (m=$m, k=$k)")
            
            if isempty(k_basis_rep)
                continue
            end

            m_k_states = (k_basis_rep, k_basis_period)
            H=Heis_XXZ_Hamiltonian_m_k_block(L,m_k_states,k,J,Delta)

            if size(H, 1) == 0
                continue
            end

            # Compute eigenvalues and eigenvectors
            e_vals= eigvals!(H)

            # Warn if imaginary part is not negligible
            for ev in e_vals
                if abs(imag(ev)) >= 1e-3
                    println("⚠ Warning: Imaginary part of eigenvalue ≥ 1e-3 for L=$L, sector (m=$m, k=$k)")
                    break 
                end
            end

            # Store only the real parts (sorted)
            energies_by_mk[(m, k)] = sort(real.(e_vals))
        end

        # Save results
        results_dir = joinpath( "results", "Full_spectrum_XXZ_Heis")
        isdir(results_dir) || mkpath(results_dir)
        filename = joinpath(results_dir,
            "$(output_prefix)_L$(L)_J$(J)_Delata$(Delta).jld2")
            
        @save filename energies_by_mk
end

 
L_l = [10,12,14,16,18]
J=1.0
Delta=1.0 #cos(pi/3)
output_prefix = "XXZ_Heis_energies"
for L in L_l
    println("Processing L=$L")
compute_and_save_all(L,J,Delta; output_prefix=output_prefix)
end
