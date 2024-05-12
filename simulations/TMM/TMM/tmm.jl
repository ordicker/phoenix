
function θ_list(θ_i, indices)
    θ = similar(indices)
    θ[1] = θ_i
    counter = 1
    for (n_i,n_f) in zip(indices[1:end-1],indices[2:end])
        θ[counter+1] = asin(n_i*sin(θ[counter]))
        counter+=1
    end
    return θ
end

function make_t_r(n1, n2, θ1, θ2)
    r = (n2*cos(θ1)-n1*cos(θ2))/(n2*cos(θ1)+n1*cos(θ2))
    t = 2*n1*cos(θ1)/(n2*cos(θ1)+n1*cos(θ2))
    return (t,r)
end

function make_t_r_list(indices, θ)
    t_r_list = zeros(eltype(indices),2,length(indices)-1)
    counter = 1
    for (n_i,n_f,θ_i,θ_f) in zip(indices[1:end-1],indices[2:end],θ[1:end-1],θ[2:end])
        t_r_list[:,counter].= make_t_r(n_i,n_f,θ_i,θ_f)
        counter+=1
    end
    return t_r_list
end

k_list(indices,θ,λ) = 2*π.*indices.*cos.(θ)./λ

function tmm_aux(indices, d, θ, λ)
    δ_list = k_list(indices,θ,λ).*d
    M = zeros(eltype(indices),2,2)
    t_r_list = make_t_r_list(indices,θ)
    M[1,1]=1;M[2,2]=1
    for (δ, t_r) in zip(δ_list[1:end-1],eachcol(t_r_list))
        M *= [exp(-im*δ) 0; 0 exp(im*δ)]
        M *= [1 t_r[2];t_r[2] 1]./t_r[1]
    end
    δ=last(δ_list)
    M*= [exp(-im*δ) 0; 0 exp(im*δ)]
    return M
end

function tmm(indices, d, θ_init, λ)
    θ = θ_list(θ_init, indices)
    M = tmm_aux(indices, d, θ, λ)
    t = 1/M[1,1]
    r = M[2,1]/M[1,1]
    return [abs2(r), abs2(t)]
end
