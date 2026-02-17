using MKL,LinearAlgebra,Random

function H_S_Hamiltonian_m_block(L::Int64,m_states,dij)
    Len=length(m_states)

    H = zeros(Float64, Len, Len)
    
    @inbounds for (ind,state) in enumerate(m_states)
        @inbounds for i in 0:L-2
            @inbounds for j in i+1:L-1
                
                if ((state >> i) & 1) == ((state >> j) & 1)
                    H[ind, ind] += 0.25*dij[i+1,j+1]
                else
                    H[ind, ind] -= 0.25*dij[i+1,j+1]
                    flip_state= state ⊻  ((1<<i)| (1<<j)) 
                    flip_ind=searchsortedfirst(m_states, flip_state)
                    H[flip_ind, ind] += 0.5*dij[i+1,j+1]
                end
            end
        end
    end
    return Hermitian(H)
end


function H_S_Hamiltonian_m_k_block(L::Int64, m_k_states::Tuple,dij,k::Float64)
    k_basis_rep, k_basis_period = m_k_states
    N_k_basis=length(k_basis_rep)
    H = zeros(Complex{Float64}, N_k_basis, N_k_basis)

    @inbounds for ind in 1:N_k_basis
        @inbounds for i in 0:L-2
            @inbounds for j in i+1:L-1
                if ((k_basis_rep[ind] >> i) & 1) == ((k_basis_rep[ind] >> j) & 1)
                    H[ind, ind] += 0.25*dij[i+1,j+1]
                else
                    H[ind, ind] -= 0.25*dij[i+1,j+1]
                    flip_state= k_basis_rep[ind] ⊻  ((1<<i)| (1<<j)) 
                    ref_state,l=representative_state(flip_state, L)
                    flip_ind=searchsortedfirst(k_basis_rep, ref_state)
                    if  ((flip_ind <= N_k_basis) && ( k_basis_rep[flip_ind] == ref_state))
                        H[flip_ind, ind] += 0.5*dij[i+1,j+1]*exp(im*(2*pi*k*l)/L)*sqrt(k_basis_period[ind]/k_basis_period[flip_ind])
                    end                  
                end
            end
        end
    end
    return Hermitian(H)
end

function Heis_XXX_Hamiltonian_m_k_block(L::Int64,m_k_states::Tuple,k::Float64,J::Float64)
    k_basis_rep, k_basis_period = m_k_states
    N_k_basis=length(k_basis_rep)
    H = zeros(Complex{Float64}, N_k_basis, N_k_basis)

    @inbounds for ind in 1:N_k_basis
        @inbounds for i in 0:L-1
            j= (i+1)%(L)
                
                if ((k_basis_rep[ind] >> i) & 1) == ((k_basis_rep[ind] >> j) & 1)
                    H[ind, ind] += 0.25
                else

                    H[ind, ind] -= 0.25
                    flip_state= k_basis_rep[ind] ⊻  ((1<<i)| (1<<j)) 

                    ref_state,l=representative_state(flip_state, L)

                    flip_ind=searchsortedfirst(k_basis_rep, ref_state)

                    if  ((flip_ind <= N_k_basis) && ( k_basis_rep[flip_ind] == ref_state))

                        H[flip_ind, ind] += 0.5*exp(im*(2*pi*k*l)/L)*sqrt(k_basis_period[ind]/k_basis_period[flip_ind])
                        
                    end                  
                end          
        end
    end
    return J .* Hermitian(H)
end

function Heis_XXX_Hamiltonian_m_block(L::Int64,m_states,J::Float64)
    Len=length(m_states)

    H = zeros(Float64, Len, Len)
    
    @inbounds for (ind,state) in enumerate(m_states)
        @inbounds for i in 0:L-1
                j= (i+1)%(L)
                
                if ((state >> i) & 1) == ((state >> j) & 1)
                    H[ind, ind] += 0.25
                else
                    H[ind, ind] -= 0.25
                    flip_state= state ⊻  ((1<<i)| (1<<j)) 
                    flip_ind=searchsortedfirst(m_states, flip_state)
                    H[flip_ind, ind] += 0.5
                end
            end
        end
    return J .* Hermitian(H)
end

function Heis_XXZ_Hamiltonian_m_block(L::Int64,m_states,J::Float64 ,Delta::Float64)
    Len=length(m_states)

    H = zeros(Float64, Len, Len)
    
    @inbounds for (ind,state) in enumerate(m_states)
        @inbounds for i in 0:L-1
                j= (i+1)%(L)
                
                if ((state >> i) & 1) == ((state >> j) & 1)
                    H[ind, ind] += Delta*0.25
                else
                    H[ind, ind] -= Delta*0.25
                    flip_state= state ⊻  ((1<<i)| (1<<j)) 
                    flip_ind=searchsortedfirst(m_states, flip_state)
                    H[flip_ind, ind] += J*0.5
                end
            end
        end
    return Hermitian(H)
end

function Heis_XXZ_Hamiltonian_m_k_block(L::Int64,m_k_states::Tuple,k::Float64,J::Float64,Delta::Float64)
    k_basis_rep, k_basis_period = m_k_states
    N_k_basis=length(k_basis_rep)
    H = zeros(Complex{Float64}, N_k_basis, N_k_basis)

    @inbounds for ind in 1:N_k_basis
        @inbounds for i in 0:L-1
            j= (i+1)%(L)
                
                if ((k_basis_rep[ind] >> i) & 1) == ((k_basis_rep[ind] >> j) & 1)
                    H[ind, ind] += Delta*0.25
                else

                    H[ind, ind] -= Delta*0.25
                    flip_state= k_basis_rep[ind] ⊻  ((1<<i)| (1<<j)) 

                    ref_state,l=representative_state(flip_state, L)

                    flip_ind=searchsortedfirst(k_basis_rep, ref_state)

                    if  ((flip_ind <= N_k_basis) && ( k_basis_rep[flip_ind] == ref_state))

                        H[flip_ind, ind] += J*0.5*exp(im*(2*pi*k*l)/L)*sqrt(k_basis_period[ind]/k_basis_period[flip_ind])
                        
                    end                  
                end          
        end
    end
    return  Hermitian(H)
end

function z_disordered_H_S_Hamiltonian_m_block(L::Int64,m_states::Vector{Int64},alpha::Float64,delta::Float64)
    dij=disordered_dij_matrix(L,alpha,delta)
    return H_S_Hamiltonian_m_block(L,m_states,dij)
end

function h_disordered_Hamiltonian_m_block(L::Int, H, m_states, h::Float64)
    hj = h .* randn(L)  
    delta_diag = zeros(Float64,length(m_states))
    @inbounds for (ind, state) in enumerate(m_states)
        s = 0.0
        for i in 0:L-1
            s += hj[i+1] * (((state >> i) & 1) - 0.5) 
        end
        delta_diag[ind] = s
    end
    return Hermitian(H + Diagonal(delta_diag))
end


function z_h_disordered_H_S_Hamiltonian_m_block(L::Int64,m_states::Vector{Int64},alpha::Float64,delta::Float64,hj:: Float64)
    dij=disordered_dij_matrix(L,alpha,delta)

    h=hj*((pi/L)/sin(pi/L))^alpha .* randn(L)

    Len=length(m_states)

    H = zeros(Float64, Len, Len)

    @inbounds for (ind,state) in enumerate(m_states)
        H[ind, ind] += h[L] *(((state >> L-1) & 1)-0.5)
        @inbounds for i in 0:L-2
            H[ind, ind] += h[i+1] *(((state >> i) & 1)-0.5)
            @inbounds for j in i+1:L-1
                
                if ((state >> i) & 1) == ((state >> j) & 1)
                    H[ind, ind] += 0.25*dij[i+1,j+1]
                else
                    H[ind, ind] -= 0.25*dij[i+1,j+1]
                    flip_state= state ⊻  ((1<<i)| (1<<j)) 
                    flip_ind=searchsortedfirst(m_states, flip_state)
                    H[flip_ind, ind] += 0.5*dij[i+1,j+1]
                end
            end
        end
    end
    return Hermitian(H)
end

function Heis_XXX_HS_rand_Hamiltonian_m_block(L::Int64,m_states,alpha,delta::Float64)
    Len=length(m_states)
    dz= (2 .* rand(L) .- 1) .* (delta*0.5) 
    d=[exp(im*(2*pi/L)*(i+dz[i])) for i in 1:L]
    dij = zeros(Float64, L, L)
    @inbounds for i in 1:L
        j= (i+1)%(L)
        dij[i,j] = ((2 * pi / L) * (1/abs(d[i] - d[j])))^alpha
    end
    H = zeros(Float64, Len, Len)
    
    @inbounds for (ind,state) in enumerate(m_states)
        @inbounds for i in 0:L-1
            
                j= (i+1)%(L)
                
                if ((state >> i) & 1) == ((state >> j) & 1)
                    H[ind, ind] += dij[i+1,j+1]*0.25
                else
                    H[ind, ind] -= dij[i+1,j+1]*0.25
                    flip_state= state ⊻  ((1<<i)| (1<<j)) 
                    flip_ind=searchsortedfirst(m_states, flip_state)
                    H[flip_ind, ind] += dij[i+1,j+1]* 0.5
                end
            end
        end
    return Hermitian(H)
end

function Heis_XXX_rand_Hamiltonian_m_block(L::Int64,m_states,J::Float64,type::String,alpha::Float64=1.0)
    Jij = zeros(Float64, L, L)

    if type == "G"
        @inbounds for i in 1:L
            j= (i+1)%(L)
            Jij[i,j] = J*randn()
        end
    elseif type == "HS"
        dz= (2 .* rand(L) .- 1) .* (J*0.5) 
        d=[exp(im*(2*pi/L)*(i+dz[i])) for i in 1:L]
        @inbounds for i in 1:L
            j= (i+1)%(L)
            Jij[i,j] = ((2 * pi / L) * (1/abs(d[i] - d[j])))^alpha
        end

    elseif type == "U"
        @inbounds for i in 1:L
            j= (i+1)%(L)
            Jij[i,j] = 2*J*rand()
        end
    elseif type == "P"
        @inbounds for i in 1:L
            j= (i+1)%(L)
            Jij[i,j] = J
        end
    else
        error("Unknown disorder type. Use 'G' for Gaussian, 'U' for Uniform, or 'HS' for Haldane-Shastry type disorder.")        
    end

    Len=length(m_states)

    H = zeros(Float64, Len, Len)
    
    @inbounds for (ind,state) in enumerate(m_states)
        @inbounds for i in 0:L-1
                j= (i+1)%(L)
                
                if ((state >> i) & 1) == ((state >> j) & 1)
                    H[ind, ind] += Jij[i+1,j+1]*0.25
                else
                    H[ind, ind] -= Jij[i+1,j+1]*0.25
                    flip_state= state ⊻  ((1<<i)| (1<<j)) 
                    flip_ind=searchsortedfirst(m_states, flip_state)
                    H[flip_ind, ind] += Jij[i+1,j+1]*0.5
                end
            end
        end
    return Hermitian(H)
end