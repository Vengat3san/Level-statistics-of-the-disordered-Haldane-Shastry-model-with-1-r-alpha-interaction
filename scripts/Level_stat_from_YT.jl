using JLD2, ThreadsX

spinon_momenta(a, L) = (pi / Float64(L)) * (Float64(a) - 0.5)
energy_spinon(p, L) = 0.5 * p * (pi - p) + (pi^2) / (8.0 * Float64(L)^2)

@inline function level_stat_combined(E::AbstractVector{T}) where T<:Real
    n = length(E)
    n < 3 && return (NaN, NaN)
    tol = T(1e-10)

  
    total_r = 0.0
    count_r = 0

   
    total_r_uni = 0.0
    count_r_uni = 0

    @inbounds begin
        
        for i in 2:n-1
            s1 = E[i] - E[i-1]
            s2 = E[i+1] - E[i]
           
            if s1 > tol || s2 > tol
                total_r += (s1 < s2) ? (s1 / s2) : (s2 / s1)
                count_r += 1
            end
        end

        s_uni_prev = -1.0 
        last_val = E[1]
        
        for i in 2:n
            s_uni_curr = E[i] - last_val
            if s_uni_curr > tol
                if s_uni_prev > 0.0
                    total_r_uni += (s_uni_curr < s_uni_prev) ? (s_uni_curr / s_uni_prev) : (s_uni_prev / s_uni_curr)
                    count_r_uni += 1
                end
                s_uni_prev = s_uni_curr
                last_val = E[i]
            end
        end
    end

    res_full = count_r > 0 ? total_r / count_r : NaN
    res_uni = count_r_uni > 0 ? total_r_uni / count_r_uni : NaN
    
    return res_full, res_uni
end



@inline function is_a_dense(A::UInt64, B::UInt64)
    (count_ones(A) < count_ones(B)) && return false
    tB, tA = B, A
    while tB != 0
        lsb = tB & -tB
        tA &= -lsb 
        (count_ones(tA) < count_ones(tB)) && return false
        tB &= (tB - 1)
    end
    return true
end

@inline function calc_EP(C::UInt64, E_LUT, P_LUT)
    Es, Ps = 0.0, 0.0
    tC = C
    while tC != 0
        tz = trailing_zeros(tC)
        @inbounds Es += E_LUT[tz + 1]
        @inbounds Ps += P_LUT[tz + 1]
        tC &= tC - 1 
    end
    return Es, Ps
end

function main(L_in::Integer)
    L = UInt64(L_in)
    L_h = L >> 1
    MASK = (one(UInt64) << L) - 1
    
    LIMIT = UInt64((BigInt(4)^L_h - 1) ÷ 3) 

    E_G = -((pi^2) / 24.0) * (Float64(L) + 5.0 / Float64(L))
    P_G = -0.5 * pi * Float64(L)
    P_LUT = Float64[spinon_momenta(a, L) for a in L:-1:1]
    E_LUT = Float64[energy_spinon(p, L) for p in P_LUT]

    total_size = Int(binomial(BigInt(L), BigInt(L_h))) 
    
    E_list = Vector{Float64}(undef, total_size)
    E_list_p1 = Float64[]
    sizehint!(E_list_p1, total_size ÷ Int(L))

    idx = 1
    for nup in UInt64(L_h):-1:1
        x = (one(UInt64) << nup) - 1
        while x <= LIMIT
            x_not = ~x & MASK
            x_not_ls = (x_not << 1) | ((x_not >> (L - 1)) & MASK)
            x_rs = (x >> 1) | ((x << (L - 1)) & MASK)
            A, B = x_not & x_not_ls, x & x_rs

            if is_a_dense(A, B)
                C = A | B
                Es, Ps = calc_EP(C, E_LUT, P_LUT)
                en = E_G + Es
                @inbounds E_list[idx] = en
                if round(Int, L * mod(P_G + Ps, 2pi) / (2pi)) == 1
                    push!(E_list_p1, en)
                end
                idx += 1
            end
            
            u = x & -x
            v = x + u
            v == 0 && break
            x = v | ((v ⊻ x) >> (trailing_zeros(u) + 2))
        end
    end
    
    return E_list, E_list_p1
end

r_tilde_mean_full = Float64[]
r_tilde_mean_p1 = Float64[]
r_tilde_mean_full_unique = Float64[] 
r_tilde_mean_p1_unique = Float64[]  

Ls = 34:2:36

for L in Ls
    println("\n--- Processing L = $L ---")
    
    print("Generating energies... ")
     
    @time E, E1 = main(UInt64(L))
      
    print("Sorting full list ($(length(E)))... ")

   
    @time ThreadsX.sort!(E)

    
    @time r_full, r_uni = level_stat_combined(E)
    
    push!(r_tilde_mean_full,r_full)
    push!(r_tilde_mean_full_unique,r_uni)

    print("Sorting P=1 list ($(length(E1)))... ")
    @time ThreadsX.sort!(E1)
    @time r_p1_full, r_p1_uni = level_stat_combined(E1)
    
    push!(r_tilde_mean_p1,r_p1_full)
    push!(r_tilde_mean_p1_unique,r_p1_uni)
    
    print("Filtering unique energies (P=1)... ")
    
    E = nothing
    E1 = nothing
    GC.gc() 
    println("Memory cleared.")
end

results_dir = joinpath("results","Spinon_data_Haldane_Shastry_Model")
    isdir(results_dir) || mkpath(results_dir)
    filename = joinpath(results_dir,"Level_stat_L_34_2_36.jld2")

@save filename r_tilde_mean_full r_tilde_mean_p1 r_tilde_mean_full_unique r_tilde_mean_p1_unique



