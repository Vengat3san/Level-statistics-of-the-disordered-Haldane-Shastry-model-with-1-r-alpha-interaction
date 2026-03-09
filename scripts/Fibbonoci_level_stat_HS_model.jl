using JLD2, ThreadsX

spinon_momenta(a, L) = (pi / Float64(L)) * (Float64(a) - 0.5)

energy_spinon(a, L) = -12*a^2 + 12*a*(L+1) -6*L

@inline function level_stat_combined(E::AbstractVector{T}) where T<:Real
    n = length(E)
    n < 3 && return (NaN, NaN)

    total_r = 0.0
    count_r = 0

    total_r_uni = 0.0
    count_r_uni = 0

    @inbounds begin

        for i in 2:n-1
            s1 = E[i] - E[i-1]
            s2 = E[i+1] - E[i]

            if s1 > 0 || s2 > 0
                total_r += (s1 < s2) ? (s1 / s2) : (s2 / s1)
                count_r += 1
            end
        end

        s_uni_prev = -1.0 
        last_val = E[1]

        for i in 2:n
            s_uni_curr = E[i] - last_val
            if s_uni_curr > 0
                if s_uni_prev > 0
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

@inline function calc_EP(C::UInt64, E_LUT, P_LUT)
    Es, Ps = 0, 0.0
    tC = C
    tz_prev=-1
    diff=0
    while tC != 0
        tz_cur= trailing_zeros(tC)
        diff += (tz_cur-1) - tz_prev
        tz=tz_cur+diff
        tz_prev=tz_cur

        @inbounds Es += E_LUT[tz + 1]
        @inbounds Ps += P_LUT[tz + 1]
        tC &= tC - 1 
        
    end
    return Es, Ps
end

fibbonoci(n) = (1 <= n <= 2) ? 1 : (fibbonoci(n-1) + fibbonoci(n-2))

function main(L::Int64)
    
    E_G = -(L^3 + 5*L)
    P_G = -0.5 * pi * L
    P_LUT = Float64[spinon_momenta(a, L) for a in L:-1:1]
    E_LUT= [energy_spinon(a, L) for a in L:-1:1]

    total_size = fibbonoci(L+1)

    E_list = Vector{Int32}(undef, total_size)
    E_list_p1 = Int32[]
    sizehint!(E_list_p1, total_size ÷ L)

    idx = 1
    for nsp in Int(L/2):-1:0
       Leff=UInt64(L-nsp)
       nup =UInt64(Leff-nsp)           
    
        x = (one(UInt64) << nup) - 1
        LIMIT = one(UInt64) << Leff
        while x <= LIMIT
             
                Es, Ps = calc_EP(x, E_LUT, P_LUT)
                en = E_G + Es
                @inbounds E_list[idx] = en
                if round(Int, L * mod(P_G + Ps, 2pi) / (2pi)) == 1
                    push!(E_list_p1, en)
                end
                idx += 1
            

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

Ls = 10:2:50

for L in Ls
    println("\n--- Processing L = $L ---")

    print("Generating energies... ")

    @time E, E1 = main(Int64(L))

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
    filename = joinpath(results_dir,"Level_stat_L_10_2_50_Fibb.jld2")

@save filename r_tilde_mean_full r_tilde_mean_p1 r_tilde_mean_full_unique r_tilde_mean_p1_unique
