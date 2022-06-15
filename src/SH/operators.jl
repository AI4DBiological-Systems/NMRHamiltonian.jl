# matrix-form of spin-1/2 single spin operators,
#   in Zeeman basis.

### for understanding.
function getsingleI⁺()::Matrix{Float64}
    return [0 1; 0 0]
end

function getsingleI⁻()::Matrix{Float64}
    return [0 0; 1 0]
end

function getsingleIy()::Matrix{Complex{Float64}}
    return 1/(2*im) .* [0 1; -1 0]
end

function getsinglehalfId()::Matrix{Float64}
    return 0.5 .* [1 0; 0 1]
end



### actual use.

function getsingleIx()::Matrix{Float64}
    return 0.5 .* [0 1; 1 0]
end

function getsingleIynoim()::Matrix{Float64}
    return 0.5 .* [0 1; -1 0]
end

function getsingleIz()::Matrix{Float64}
    return 0.5 .* [1 0; 0 -1]
end

function getsingleId()::Matrix{Float64}
    return [1 0; 0 1]
end



#### product operators.

function getmultiIz(   Id::Matrix{T},
                        j::Int,
                        N::Int) where T <: Real

    # prepare base operators.
    # prepare base operators.
    operators = Vector{Matrix{T}}(undef, N)
    fill!(operators, Id)

    operators[j] = getsingleIz()

    # get product operator.
    #out = 2^(N-1) .* Kronecker.kronecker(operators...)
    out = Kronecker.kronecker(operators...)

    return out
end

function getmultiIjIk(   Id,
                        Ix,
                        Iy_no_im,
                        Iz,
                        j::Int,
                        k::Int,
                        N::Int,
                        c::T) where T <: Real

    ### x.
    # prepare base operators.
    x_operators = Vector{Matrix{T}}(undef, N)
    fill!(x_operators, Id)

    x_operators[j] = c .* Ix
    x_operators[k] = Ix

    # get product operator.
    x_out = Kronecker.kronecker(x_operators...)



    ### z.
    # prepare base operators.
    z_operators = Vector{Matrix{T}}(undef, N)
    fill!(z_operators, Id)

    z_operators[j] = c .* Iz
    z_operators[k] = Iz

    # get product operator.
    z_out = Kronecker.kronecker(z_operators...)



    ### y.
    # prepare base operators.
    y_operators = Vector{Matrix{T}}(undef, N)
    fill!(y_operators, Id)

    y_operators[j] = c .* Iy_no_im
    y_operators[k] = Iy_no_im

    # get product operator.
    y_out = -Kronecker.kronecker(y_operators...)

    out = x_out + y_out + z_out

    return out
end

function getmultiIjIk0(   Id::Matrix{T},
                        Ix::Matrix{T},
                        Iy_no_im::Matrix{T},
                        Iz::Matrix{T},
                        j::Int,
                        k::Int,
                        N::Int) where T <: Real

    ### x.
    # prepare base operators.
    x_operators = Vector{Matrix{T}}(undef, N)
    fill!(x_operators, Id)

    x_operators[j] = Ix
    x_operators[k] = Ix

    # get product operator.
    x_out = Kronecker.kronecker(x_operators...)



    ### z.
    # prepare base operators.
    z_operators = Vector{Matrix{T}}(undef, N)
    fill!(z_operators, Id)

    z_operators[j] = Iz
    z_operators[k] = Iz

    # get product operator.
    z_out = Kronecker.kronecker(z_operators...)



    ### y.
    # prepare base operators.
    y_operators = Vector{Matrix{T}}(undef, N)
    fill!(y_operators, Id)

    y_operators[j] = Iy_no_im
    y_operators[k] = Iy_no_im

    # get product operator.
    y_out = -Kronecker.kronecker(y_operators...)

    out = x_out + y_out + z_out

    return out
end
