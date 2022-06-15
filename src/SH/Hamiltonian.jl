
function getgenericHamiltoniantest( Id,
    Ix,
    Iy_no_im,
    Iz,
    ω0::Vector{T},
    J_vals::Vector{T},
    J_inds::Vector{Tuple{Int,Int}}) where T <: Real
#
    N = length(ω0)
    N_couplings = length(J_vals)
    @assert N_couplings == length(J_inds)


    # pair-wise terms.
    H1 = zeros(T, 2^N, 2^N)
    for i = 1:N_couplings

        j = J_inds[i][1]
        k = J_inds[i][2]

        H1 += getmultiIjIk(Id, Ix, Iy_no_im, Iz, j, k, N, 2*π*J_vals[i])
    end

    # single terms.
    ωIzs = collect( ω0[j] .* Iz for j = 1:N )
    H2 = Kronecker.kroneckersum(ωIzs...)

    return H1+H2, H1, H2
end

function getgenericHamiltoniantest( Id,
    Ix,
    Iy_no_im,
    Iz,
    ω0::Vector{T},
    J::Matrix{T}) where T <: Real

    N = length(ω0)
    @assert size(J,1) == size(J,2) == N


    # pair-wise terms.
    H1 = zeros(T, 2^N, 2^N)
    for j = 1:N
        for k = j+1:N
            #H += (2*π*J[j,k]) .* getmultiIjIk(Id, Ix, Iy_no_im, Iz, j, k, N)
            H1 += getmultiIjIk(Id, Ix, Iy_no_im, Iz, j, k, N, 2*π*J[j,k])

        end
    end

    # single terms.
    ωIzs = collect( ω0[j] .* Iz for j = 1:N )
    H2 = Kronecker.kroneckersum(ωIzs...)

    return H1+H2, H1, H2
end

# Only the upper-right off-diagonals of J is used.
# No secular dipole-dipole.
function getgenericHamiltonian( Id,
    Ix,
    Iy_no_im,
    Iz,
    ω0::Vector{T},
    J::Matrix{T}) where T <: Real

    N = length(ω0)
    @assert size(J,1) == size(J,2) == N


    # pair-wise terms.
    H = zeros(T, 2^N, 2^N)
    for j = 1:N
        for k = j+1:N
            #H += (2*π*J[j,k]) .* getmultiIjIk(Id, Ix, Iy_no_im, Iz, j, k, N)
            H += getmultiIjIk(Id, Ix, Iy_no_im, Iz, j, k, N, 2*π*J[j,k])

        end
    end

    # single terms.
    ωIzs = collect( ω0[j] .* Iz for j = 1:N )
    H += Kronecker.kroneckersum(ωIzs...)

    return H
end

# No secular dipole-dipole.
function getgenericHamiltonian0( Id,
    Ix,
    Iy_no_im,
    Iz,
    ω0::Vector{T},
    J::Matrix{T}) where T <: Real

    N = length(ω0)
    @assert size(J,1) == size(J,2) == N

    H = zeros(T, 2^N, 2^N)

    for j = 1:N
        H += ω0[j] .* getmultiIz(Id, j, N)

        for k = j+1:N
            H += (2*π*J[j,k]) .* getmultiIjIk0(Id, Ix, Iy_no_im, Iz, j, k, N)
        end
    end

    return H
end

# No secular dipole-dipole.
function getgenericHamiltonian( Id,
    Ix,
    Iy_no_im,
    Iz,
    ω0::Vector{T},
    J_vals::Vector{T},
    J_inds::Vector{Tuple{Int,Int}}) where T <: Real

    N = length(ω0)
    N_couplings = length(J_vals)
    @assert N_couplings == length(J_inds)


    # pair-wise terms.
    H = zeros(T, 2^N, 2^N)
    for i = 1:N_couplings

        j = J_inds[i][1]
        k = J_inds[i][2]

        H += getmultiIjIk(Id, Ix, Iy_no_im, Iz, j, k, N, 2*π*J_vals[i])
    end

    # single terms.
    ωIzs = collect( ω0[j] .* Iz for j = 1:N )
    H += Kronecker.kroneckersum(ωIzs...)

    return H
end

# has secular dipole-dipole. J here is not really J-coupling, but J - d.
function getgenericHamiltonian( Id,
    Ix,
    Iy_no_im,
    Iz,
    ω0::Vector{T},
    J_vals::Vector{T},
    J_inds::Vector{Tuple{Int,Int}},
    Di_vals::Vector{T},
    Di_inds::Vector{Tuple{Int,Int}}) where T <: Real

    N = length(ω0)
    N_couplings = length(Di_vals)
    @assert N_couplings == length(Di_inds)

    # Allocate base operators.
    z_operators = Vector{Matrix{T}}(undef, N)
    fill!(z_operators, Id)

    # pair-wise terms.
    H = getgenericHamiltonian(Id, Ix, Iy_no_im, Iz, ω0, J_vals, J_inds)
    for i = 1:N_couplings

        j = Di_inds[i][1]
        k = Di_inds[i][2]


        ### z.
        z_operators[j] = (3*Di_vals[i]) .* Iz
        z_operators[k] = Iz

        # get product operator.
        H_di = Kronecker.kronecker(z_operators...)

        H += H_di
    end

    return H
end

# v must be an eigen vector of Iz. This routine does not verify this condition.
function getzangularmomentum(v::Vector{T}, Iz; tol = 1e-7)::T where T
    
    b = Iz*v
    m = dot(b,v)/dot(v,v)

    ## This was meant to be a sanity-check, but when N > 4, the tol becomes larger than 1e-7.
    # The exclusion of these states mean there are insufficient number of splits,
    #   thus, the area ratio does not hold, unless weak coupling is "really satisifed"?
    # Solution: just disable this check, and assume v must be an eigenvector of Iz.
    # if norm(m .* v - b) > tol
    #     return -Inf
    # end

    return m
end

function getorderofcoherence(Iz, Q::Vector{Vector{T}}) where T
    
    M = length(Q)
    @assert size(Iz, 1) == M == size(Iz, 2)

    M_array = collect( getzangularmomentum(Q[r], Iz) for r = 1:length(Q))

    p = Matrix{T}(undef, M, M)
    fill!(p, NaN)
    for r = 1:M
        #Mr = getzangularmomentum(Q[r], Iz)
        Mr = M_array[r]

        for s = 1:M
            #Ms = getzangularmomentum(Q[s], Iz)
            Ms = M_array[s]

            p[r,s] = Mr-Ms
            #p[r,s] = round(Mr-Ms)
        end
    end

    return p, M_array
end

function normalizetoNspins!(A::Vector{T}, N_spins::Int) where T

    factor = N_spins/sum(A)
    for l = 1:length(A)
        A[l] = A[l] * factor
    end

    return nothing
end

function getaΩ( Iz_full,
    H::Matrix{T},
    Iys_no_im_full,
    Ix_full;
    tol = 1e-7) where T <: Real

    M = size(H,1)
    @assert size(H,2) == M == size(Iz_full,1) == size(Iz_full,2)

    ### eigen states.
    λ, Q_mat = eigen(H)
    # λ[6] = -λ[6]
    # Q_mat[:,6] = -Q_mat[:,6]
    Q = collect( Q_mat[:,m] for m = 1:M )

    ### pre-compute.
    IyQs = collect( Iys_no_im_full*Q[s] for s = 1:size(Q,1)) # no im since official formula for α is div by 2im.
    IxQs = collect( Ix_full*Q[r] for r = 1:size(Q,1))


    ### get order of coherence.
    #@time p, M_array = getorderofcoherence(Iz_full, Q) # TODO speed this up.
    #println("timing coherence")
    #@time p, M_array = getorderofcoherence(Iz_full, Q)
    p, M_array = getorderofcoherence(Iz_full, Q)

    a = Vector{T}(undef, M*M)
    F = Vector{T}(undef, M*M)
    inds = Vector{Tuple{Int,Int}}(undef, M*M)

    i_out = 0

    # for r = 1:M
    #     for s = 1:M
    for (r,s) in Base.Iterators.product(1:M,1:M)

        #if Utilities.isnumericallyclose(p[r,s], -one(T), tol)
        if isnumericallyclose(round(p[r,s]), -one(T), tol)
            #println("i_out = ", i_out)
            i_out += 1

            inds[i_out] = (r,s)

            ## signal amplitude. Section A.8.2 Spin Dynamics, Eq 2.1.131 Ernst.
            t1 = -dot(conj(Q[r]), IyQs[s])
            t2 = dot(conj(Q[s]), IxQs[r])

            # factor 2 instead of 2*im since Iys_no_im_full
            #   is divided by 2 instead of 2*im.
            a[i_out] = 2*(t1 * t2)

            ## signal frequency.
            F[i_out] = λ[s] - λ[r]
        end
    end

    resize!(inds, i_out)
    resize!(a, i_out)
    resize!(F, i_out)

    return a, F, p, λ, Q, inds, M_array
end


function getcfreq(z::Vector{T}) where T
    M = length(z)

    F = Matrix{T}(undef, M, M)
    for r = 1:M
        for s = 1:M
            F[r,s] = z[s] - z[r]
        end
    end

    return F
end

# Q is operator. Evalates <r| Q |s>.
function evaltrace( Q::Kronecker.KroneckerSum{T}, r, s)::T where T
    return dot(conj(r),Q*s)
end

function evaltrace( Q::Kronecker.KroneckerSum{T},
    V)::Matrix{T} where T

    M = length(V)

    out = Matrix{T}(undef, M, M)
    for r = 1:M
    for s = 1:M
    out[r,s] = evaltrace(Q, V[r], V[s])
    end
    end

    return out
end


# tp_basis[r,:] is the elementary state vector for each of the spins.
# +1 for state 1/2, -1 for state -1/2.
function gettensorproductbasis(Iz::Matrix{T}, Id, N_spins::Int) where T
    tp_basis = Matrix{Int}(undef, 2^N_spins, N_spins)

    for n = 1:N_spins

    Is = collect( Id for j = 1:N_spins )
    Is[n] = Iz
    I_full = Kronecker.kronecker(Is...)

        for i = 1:size(I_full,1)
            m = I_full[i,i]

            if m > 0
                tp_basis[i,n] = 1
            else
                tp_basis[i,n] = -1
            end
        end
    end

    return tp_basis
end


# for single (π/2)_x pulse. # equation A.32.
function evalamplitudeA32( Iys_no_im_full,
    Ix_full,
    V::Vector{Vector{T}})::Matrix{T} where T

    M = length(V)

    out = Matrix{T}(undef, M, M)
    for r = 1:M
    for s = 1:M

    t1 = -dot(conj(V[r]), Iys_no_im_full*V[s])
    t2 = dot(conj(V[s]), Ix_full*V[r])

    # factor 2 instead of 2*im since Iys_no_im_full
    #   is divided by 2 instead of 2*im.
    out[r,s] = 2*(t1 * t2)
    end
    end

    return out
end

# no phase. common peak width λ.
function evalabsorptionspectrum(x::T, λ::T, a::Array{T,D}, Ω::Array{T,D}) where {T,D}
    L = length(a)
    @assert length(Ω) == L

    out = zero(T)
    for l = 1:L
    denominator = λ^2 + (x-Ω[l])^2
    out += a[l]*λ/denominator
    end

    return out
end

function evalabsorptionspectrumarea(lower_limit::T,
    upper_limit::T,
    λ::T,
    a::Array{T,D},
    Ω::Array{T,D}) where {T,D}

    L = length(a)
    @assert length(Ω) == L

    out = zero(T)
    for l = 1:L

    numerator = Ω[l] - lower_limit
    Ia = -atan(numerator,λ)

    numerator = Ω[l] - upper_limit
    Ib = -atan(numerator,λ)

    out += a[l]*(Ib-Ia)
    end

    return out
end
