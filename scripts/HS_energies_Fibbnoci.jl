using JLD2, ThreadsX

spinon_momenta(a, L) = (pi / Float64(L)) * (Float64(a) - 0.5)

energy_spinon(a, L) = -12*a^2 + 12*a*(L+1) -6*L

@inline function calc_E(C::UInt64, E_LUT)
    Es = 0
    tC = C
    tz_prev=-1
    diff=0
    while tC != 0
        tz_cur= trailing_zeros(tC)
        diff += (tz_cur-1) - tz_prev
        tz=tz_cur+diff
        tz_prev=tz_cur

        @inbounds Es += E_LUT[tz + 1]
        tC &= tC - 1 

    end
    return Es
end

fibbonoci(n) = (1 <= n <= 2) ? 1 : (fibbonoci(n-1) + fibbonoci(n-2))

function main(L::Int64)

    E_G = -(L^3 + 5*L)
    P_G = -0.5 * pi * L
    P_LUT = Float64[spinon_momenta(a, L) for a in L:-1:1]
    E_LUT= [energy_spinon(a, L) for a in L:-1:1]

    total_size = fibbonoci(L+1)

    E_list = Vector{Int32}(undef, total_size)
    
    idx = 1
    for nsp in Int(L/2):-1:0
       Leff=UInt64(L-nsp)
       nup =UInt64(Leff-nsp)           

        x = (one(UInt64) << nup) - 1
        LIMIT = one(UInt64) << Leff
        while x <= LIMIT

                Es = calc_E(x, E_LUT)
                en = E_G + Es
                @inbounds E_list[idx] = en

                idx += 1

            u = x & -x
            v = x + u
            v == 0 && break
            x = v | ((v ⊻ x) >> (trailing_zeros(u) + 2))
        end
    end

    return E_list
end


Ls = 10:2:52

for L in Ls
    println("\n--- Processing L = $L ---")

    print("Generating energies... ")

    E = main(Int64(L))

    ThreadsX.sort!(E)
    unique!(E)
 
   results_dir = joinpath("results","Spinon_data_Haldane_Shastry_Model")
    isdir(results_dir) || mkpath(results_dir)
    filename = joinpath(results_dir,"Unique _energies_$(L)_Fibb.jld2")

@save filename E


    E = nothing

    GC.gc() 
    println("Memory cleared.")
end

 