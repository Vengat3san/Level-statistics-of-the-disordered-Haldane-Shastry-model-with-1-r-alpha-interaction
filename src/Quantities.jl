

function S_squared_operator_m_k_block(L::Int64, m::Float64, k::Float64)
    m_states = m_block(L, m)
    k_basis_rep = []
    k_basis_period = []

    for state in m_states
        R = check_state(state, L, k)
        if R >= 0
            push!(k_basis_rep, state)
            push!(k_basis_period, R)
        end
    end

    N_k_basis = length(k_basis_rep)
    S2 = spzeros(ComplexF64, N_k_basis, N_k_basis)

    for ind in 1:N_k_basis
        state = k_basis_rep[ind]
        for i in 0:L-2
            for j in i+1:L-1
                bit_i = (state >> i) & 1
                bit_j = (state >> j) & 1

                if bit_i == bit_j
                    S2[ind, ind] += 0.5  # 2 * (1/4)
                else
                    S2[ind, ind] += -0.5 # 2 * (-1/4)

                    # flip terms from S^+_i S^-_j + S^-_i S^+_j
                    flip_state = state ⊻ ((1 << i) | (1 << j))
                    ref_state, l = representative_state(flip_state, L)
                    flip_ind = searchsortedfirst(k_basis_rep, ref_state)

                    if flip_ind <= N_k_basis && k_basis_rep[flip_ind] == ref_state
                        coeff = 2 * 0.5 * exp(im * 2π * k * l / L) * sqrt(k_basis_period[ind] / k_basis_period[flip_ind])
                        S2[flip_ind, ind] += coeff
                    end
                end
            end
        end
    end

    # Add constant term: S² = sum_{i<j} 2 S_i·S_j + (3/4)L·𝟙
    S2 += (3/4) * L * I # I is identity operator in sparse matrix algebra
    return Hermitian(S2)
end


#=
function S_squared_operator_m_k_block(L::Int64, m::Float64, k::Float64)
    m_states = m_block(L, m)
    k_basis_rep = []
    k_basis_period = []

    for state in m_states
        R = check_state(state, L, k)
        if R >= 0
            push!(k_basis_rep, state)
            push!(k_basis_period, R)
        end
    end

    N_k_basis = length(k_basis_rep)
    S2 = spzeros(ComplexF64, N_k_basis, N_k_basis)

    for ind in 1:N_k_basis
        state = k_basis_rep[ind]
        for i in 0:L-2
            for j in i+1:L-1
                bit_i = (state >> i) & 1
                bit_j = (state >> j) & 1
                if bit_i == bit_j
                    # S^z_i S^z_j term
                    S2[ind, ind] += 0.25
                else
                    # S^z_i S^z_j term
                    S2[ind, ind] -= 0.25
                    # S^+_i S^-_j + S^-_i S^+_j term
                    flip_state = state ⊻ ((1 << i) | (1 << j))
                    ref_state, l = representative_state(flip_state, L)
                    flip_ind = searchsortedfirst(k_basis_rep, ref_state)
                    if flip_ind <= N_k_basis && k_basis_rep[flip_ind] == ref_state
                        S2[flip_ind, ind] += 0.5 * exp(im * 2π * k * l / L) * sqrt(k_basis_period[ind] / k_basis_period[flip_ind])
                    end
                end
            end
        end
    end

    return Hermitian(S2)
end
=#