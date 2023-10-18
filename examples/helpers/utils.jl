





function combinevectors(x::Vector{Vector{T}})::Vector{T} where T

    if isempty(x)
        return Vector{T}(undef, 0)
    end

    N = sum(length(x[i]) for i = 1:length(x))

    y = Vector{T}(undef,N)

    st_ind = 0
    fin_ind = 0
    for i = 1:length(x)
        st_ind = fin_ind + 1
        fin_ind = st_ind + length(x[i]) - 1

        y[st_ind:fin_ind] = x[i]
    end

    return y
end


function getPsnospininfo(As, hz2ppmfunc)

    ΩS_ppm = Vector{Vector{Float64}}(undef, length(As))

    for (n,A) in enumerate(As)

        ΩS_ppm[n] = hz2ppmfunc.( combinevectors(A.Ωs) ./ (2*π) )

        tmp = hz2ppmfunc.( A.Ωs_singlets ./ (2*π) )
        push!(ΩS_ppm[n], tmp...)
    end

    return ΩS_ppm
end

function array2matrix(X::Array{Vector{T},L})::Matrix{T} where {T,L}

    N = length(X)
    D = length(X[1])

    out = Matrix{T}(undef,D,N)
    for n = 1:N
        out[:,n] = X[n]
    end

    return out
end



##### for investigating prune combination coherences.
function allabssmaller(x, threshold)
    if maximum(abs.(x)) > threshold
        return false
    end

    return true
end

function sumabssmaller(x, threshold)
    if sum(abs.(x)) > threshold
        return false
    end

    return true
end

###### for evaluating resonance groups and singlets.

# barebones.
function evalzerophasecl1Darray(u_rad, αs::Vector{T}, Ωs::Vector{T}, λ::T)::Complex{T} where T <: AbstractFloat

    out = zero(Complex{T})
    for l = 1:length(αs)
        out += evalzerophaseclpartitionelement(u_rad, αs[l], Ωs[l], λ)
    end

    return out
end

function evalzerophaseclpartitionelement(r,
    α::T, Ω::T, λ::T)::Complex{T} where T <: AbstractFloat

    out = α/(λ+im*(r-Ω))

    return out
end

function evalzerophaseclresonancegroup(u_rad, A,
    i::Int, k::Int, λ::T)::Complex{T} where T <: AbstractFloat

    out = zero(Complex{T})
    for l in A.parts[i][k]
        out += evalzerophaseclpartitionelement(u_rad, A.αs[i][l], A.Ωs[i][l], λ)
    end

    return out
end

function evalzerophaseclsinglets(u_rad, A, λ::T)::Complex{T} where T <: AbstractFloat

    out = zero(Complex{T})
    for i = 1:length(A.Ωs_singlets)
        out += evalzerophaseclpartitionelement(u_rad, A.αs_singlets[i], A.Ωs_singlets[i], λ)
    end

    return out
end

function evalzerophaseclmixture(u_rad, As, λ::T)::Complex{T} where T <: AbstractFloat

    out = zero(Complex{T})
    for n = 1:length(As)
        for i = 1:length(As[n].parts)
            for k = 1:length(As[n].parts[i])

                out += evalzerophaseclresonancegroup(u_rad, As[n], i, k, λ)
            end
        end

        out += evalzerophaseclsinglets(u_rad, As[n], λ)
    end

    return out
end

function getqs(A, λ::T) where T <: AbstractFloat
    #
    N_sys = length(A.parts)

    qs = Vector{Vector{Function}}(undef, N_sys)
    for i = 1:N_sys

        N_groups = length(A.parts[i])
        qs[i] = collect( uu->evalzerophaseclresonancegroup(uu, A, i, k, λ) for k = 1:N_groups)
    end

    return qs
end
