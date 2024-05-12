using Plots
includet("./tmm.jl")

θ_init = 0.0
λ = 3.0
# Silicon, Platinum,SiO2  
n = [3.4699, 4.90+11.70im ,1.5209]

function f(d)
    R,T = tmm([1.0, 3.4699, 4.90+11.70im ,1.5209 ,1.0],[0.1,0.5,d,10.0,0.1],θ_init,λ)
    return 1-T-R
end
