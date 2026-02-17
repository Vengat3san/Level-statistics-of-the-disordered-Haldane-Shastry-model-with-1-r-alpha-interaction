##===========================================##
## nup=m+(L/2)##
## m=nup-(L/2)##
##===========================================##
 
 
 function ishftc(i::Integer, shift::Integer, size::Integer=bitwidth(i))
    shift_amt = mod(shift, size)
    mask = (1 << size) - 1
    right_bits = i & mask
    rotated = ((right_bits << shift_amt) | (right_bits >> (size - shift_amt))) & mask
    result = (i & ~mask) | rotated
    return result
end

 function m_block(L::Int64, m::Float64) #Gosper Hack
    nup = Int64(m+(L ÷ 2))
    states = Int[]
    if nup > L || nup < 0
        return states
    elseif nup == 0
        return [0]
    end
    x = (1 << nup) - 1
    limit = 1 << L
    while x < limit
        push!(states, x)
        u = x & -x
        v = x + u
        x = v + (((v ⊻ x) ÷ u) >> 2)
    end
    return states
end

 function check_state(state::Int64, L::Int64,k::Float64)
    R,t=-1,state
    @inbounds for i in 1:L 
        t=ishftc(t, 1, L)
        if t < state 
            return R
        elseif t == state
            if k % (L/i) != 0
                return R
            end
            R=i
            return R
        end
    end
end

 function representative_state(state::Int64, L::Int64)
    ref_state,transl_state,l=state,state,0
    @inbounds for i in 1:L-1
        transl_state=ishftc(transl_state, 1, L)
        if transl_state < ref_state 
            ref_state=transl_state
            l=i
        end
    end
    return ref_state,l
end

 function dij_matrix(L::Int64, alpha::Float64)
    d=[exp(im*(2*pi/L)*i) for i in 1:L]
    dij = zeros(Float64, L, L)
    @inbounds for i in 1:L-1
        @inbounds for j in i+1:L
            dij[i,j] = ((2 * pi / L) * (1/abs(d[i] - d[j])))^alpha
        end
    end
    return dij
end

 function disordered_dij_matrix(L::Int64, alpha::Float64,delta::Float64)
    dz= (2 .* rand(L) .- 1) .* (delta*0.5) 
    d=[exp(im*(2*pi/L)*(i+dz[i])) for i in 1:L]
    dij = zeros(Float64, L, L)
    @inbounds for i in 1:L-1
        @inbounds for j in i+1:L
            dij[i,j] = ((2 * pi / L) * (1/abs(d[i] - d[j])))^alpha
        end
    end
    return dij
end

 function m_k_states_period(L::Int64, k::Float64,m::Float64)
    m_states = m_block(L, m)
    k_basis_rep = []
    k_basis_period = []
    @inbounds for state in m_states
        R = check_state(state, L, k)
        if R >= 0
            push!(k_basis_rep, state)
            push!(k_basis_period, R)
        end
    end
    return k_basis_rep, k_basis_period
end

 

 function precompute_all_mk_basis(L::Int)
    mk_basis_dict = Dict{Tuple{Float64,Float64}, Tuple{Vector{Int64}, Vector{Int}}}()
    for nup in 0:L
        m = Float64(nup - (L ÷ 2)) 
        for k in 0.0:(L - 1)
            k_basis_rep, k_basis_period = m_k_states_period(L, k, m)
            mk_basis_dict[(m, k)] = (k_basis_rep, k_basis_period)
        end
    end
    return mk_basis_dict
end

 function precompute_mk_basis(L::Int64, M::Array{Float64}, K::Array{Float64})
    mk_basis_dict = Dict{Tuple{Float64,Float64}, Tuple{Vector{Int64}, Vector{Int}}}()
    
    for m in M
        for k in K
            k_basis_rep, k_basis_period = m_k_states_period(L, k, m)
            mk_basis_dict[(m, k)] = (k_basis_rep, k_basis_period)
        end
    end
    return mk_basis_dict
end

function average_histograms(s_pdf_list::Vector{StatsBase.Histogram})
    N = length(s_pdf_list)
    n_bins = length(s_pdf_list[1].weights)

    edges = collect(s_pdf_list[1].edges[1]) 

    avg_counts = zeros(n_bins)
    for h in s_pdf_list
        avg_counts .+= h.weights
    end
    avg_counts ./= N

    bin_centers = 0.5 .* (edges[1:end-1] .+ edges[2:end])
    return bin_centers, avg_counts
end
