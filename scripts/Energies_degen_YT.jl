using JLD2,ThreadsX

energy_spinon(a, L) = -12*a^2 + 12*a*(L+1) -6*L

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

@inline function calc_E(C::UInt64,E_LUT)
    Es = Int64(0)
    while C != zero(UInt64)
        tz = trailing_zeros(C)
        @inbounds Es += E_LUT[tz + 1]
        C &= C - 1 
    end
    return Es
end

function fibonacci(n)
  if n <= 2
    return 1
  end
  return fibonacci(n - 1) + fibonacci(n - 2)
end

struct State
    e::Int32
    s::Int32
end

function main(L::UInt64)
     L_c= Int64(L)
     L_h = UInt64(L // 2)
     MASK = (one(UInt) << L) - 1

     LIMIT = UInt64((BigInt(4)^L_h - 1) ÷ 3) 

     E_G = -(L_c^3 + 5*L_c)
     E_LUT= [energy_spinon(a, L_c) for a in L_c:-1:1]

    total_size = binomial(L, L_h) 
    
    E_S_struct = Vector{State}(undef, total_size)
    
    
    idx=1
    
    for nup in L_h:-1:1
        S = Int64(L_h - nup)
        
        x = (one(UInt64) << nup) - 1
        
        while x <= LIMIT
            
            x_not = ~x & MASK
            x_not_ls = (x_not << 1) | ((x_not >> (L - 1)) & MASK)
            x_rs = (x >> 1) | ((x << (L - 1)) & MASK)
            
            A = x_not & x_not_ls
            B = x & x_rs
            
            if is_a_dense(A, B)
            
                C = A | B
                Esp = calc_E(C,E_LUT)
                E_st= E_G + Esp
                @inbounds E_S_struct[idx] = State(E_st,S)
                idx += 1
              
            end

            u = x & -x
            v = x + u
            if v == 0 
                print("Loop broken")
                break 
            end
            x = v | ((v ⊻ x) >> (trailing_zeros(u) + 2))
        end
    end
    
    E_all = E_G + sum(E_LUT)
    @inbounds E_S_struct[idx] = State(E_all,L_h)

    ThreadsX.sort!(E_S_struct, by=x->x.e)

    n=fibonacci(L_c+1)
    eigen_vals = Vector{Float64}()
    sizehint!(eigen_vals, n)
    degeneracies = Vector{Int32}()
    sizehint!(degeneracies, n)

    @inbounds begin
        curr_e = E_S_struct[1].e
        curr_deg = 2 * E_S_struct[1].s + 1
        
        for i in 2:length(E_S_struct)
            e = E_S_struct[i].e
            s = E_S_struct[i].s

            
            if e == curr_e
                curr_deg += 2*s + 1
            else
                push!(eigen_vals, curr_e)
                push!(degeneracies, curr_deg)
                curr_e = e
                curr_deg = 2*s + 1
            end
        end

        push!(eigen_vals, curr_e)
        push!(degeneracies, curr_deg)
    end   
    return eigen_vals, degeneracies
end


for L in 34:2:36
    
    println("For L=$(L) starts")
    @time eigen_vals, degeneracies = main(UInt64(L));

    results_dir = joinpath("results","Spinon_data_Haldane_Shastry_Model")
    isdir(results_dir) || mkpath(results_dir)
    filename = joinpath(results_dir,"Unique_Energies_degen_L_$(L).jld2")
    @save filename eigen_vals degeneracies
    println("For L=$(L) end")
    
end  
