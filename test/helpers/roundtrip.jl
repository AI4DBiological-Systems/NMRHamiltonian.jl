


function comparevectors(x::Vector{String}, y::Vector{String}, args...)::Bool
    return all(x .== y)
end

function comparevectors(
    x::Vector{T},
    y::Vector{T},
    zero_tol,
    )::Bool where T <: Real
    
    @assert length(x) == length(y)

    return norm(x-y) < zero_tol
end

function comparevectors(
    x::Tuple{T,T},
    y::Tuple{T,T},
    zero_tol,
    )::Bool where T <: Real

    @assert length(x) == length(y)

    discrepancy = 0
    for j in eachindex(x)
        discrepancy += (x[j][begin] - y[j][begin])^2
        discrepancy += (x[j][end] - y[j][end])^2
    end
        
    return discrepancy < zero_tol
end


function comparevectors(
    x::Vector,
    y::Vector,
    zero_tol,
    )::Bool
    
    @assert length(x) == length(y)

    if isempty(x)
        return true
    end
    
    return all( comparevectors(x[i], y[i], zero_tol) for i in eachindex(x))
end



function roundtripJSON(
    As::Vector{NMRHamiltonian.SHType{T}};
    zero_tol = 1e-12,
    ) where T <: AbstractFloat

    file_name = "test.json"

    S = NMRHamiltonian.serializemixture(As)

    NMRHamiltonian.saveasJSON(file_name, S)

    # load JSON3
    file_path = file_name
    json_string = read(file_path)
    W = JSON3.read(json_string)
    
    # round trip test.
    AS2 = NMRHamiltonian.deserializemixture(Dict(W))
    S2 = NMRHamiltonian.serializemixture(AS2)

    for (key,val) in S2
        if typeof(val) <: Vector

            @test comparevectors(val, S[key], zero_tol)
        end
    end

    return nothing
end

function roundtripBSON(
    As::Vector{NMRHamiltonian.SHType{T}};
    zero_tol = 1e-12,
    ) where T <: AbstractFloat

    file_name = "test.bson"

    S = NMRHamiltonian.serializemixture(As)

    BSON.bson(file_name, S)

    ## BSON
    W = BSON.load(file_name)

    # round trip test.
    AS2 = NMRHamiltonian.deserializemixture(Dict(W))
    S2 = NMRHamiltonian.serializemixture(AS2)

    for (key,val) in S2
        if typeof(val) <: Vector

            @test comparevectors(val, S[key], zero_tol)
        end
    end

    return nothing
end


function roundtripJSON(
    Phys::Vector{NMRHamiltonian.PhysicalParamsType{T}},
    molecule_entries::Vector{String};
    zero_tol = 1e-12,
    ) where T <: AbstractFloat

    S = NMRHamiltonian.serializephysicalparams(Phys, molecule_entries)

    file_name = "test.json"
    NMRHamiltonian.saveasJSON(file_name, S)

    # load JSON3
    file_path = file_name
    json_string = read(file_path)
    W = JSON3.read(json_string)
    
    # round trip test.
    Phys2, m_entries2 = NMRHamiltonian.deserializephysicalparams(Dict(W))
    S2 = NMRHamiltonian.serializephysicalparams(Phys2, m_entries2)

    for (key,val) in S2
        if typeof(val) <: Vector

            @test comparevectors(val, S[key], zero_tol)
        end
    end

    return nothing
end

function roundtripBSON(
    Phys::Vector{NMRHamiltonian.PhysicalParamsType{T}},
    molecule_entries::Vector{String};
    zero_tol = 1e-12,
    ) where T <: AbstractFloat    

    file_name = "test.bson"

    S = NMRHamiltonian.serializephysicalparams(Phys, molecule_entries)

    BSON.bson(file_name, S)

    ## BSON
    W = BSON.load(file_name)

    # round trip test.
    Phys2, m_entries2 = NMRHamiltonian.deserializephysicalparams(Dict(W))
    S2 = NMRHamiltonian.serializephysicalparams(Phys2, m_entries2)

    for (key,val) in S2
        if typeof(val) <: Vector

            @test comparevectors(val, S[key], zero_tol)
        end
    end

    return nothing
end