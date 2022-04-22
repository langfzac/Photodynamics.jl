#using StaticArrays, LinearAlgebra, Photodynamics
import Base: show
import Photodynamics: @copyfields

# Define traits for differentiability
abstract type Differentiability end
struct Differentiable <: Differentiability end
struct NonDifferentiable <: Differentiability end
Differentiability(::Type) = NonDifferentiable() # Default

## Abstract types for dynamical model
abstract type AbstractDynamicalModel end

## Types for Expansion based models (ie. PK20)
abstract type AbstractExpansion{T} end

struct ExpansionChain{ET<:AbstractExpansion, T<:Real}
    expansions::Vector{ET}
    body_index::Int64
    times::Vector{T}

    function ExpansionChain(expansions::AbstractVector{<:AbstractExpansion{T}}, body_index::Int) where T<:Real
        times = [ex.t0 for ex in expansions]

        ET = eltype(expansions)
        return new{ET,T}(expansions, body_index, times)
    end
end

Base.show(io::IO, ::MIME"text/plain", ec::ExpansionChain) = begin
    println(io, "ExpansionChain")
    print(io, "$(length(ec.expansions)) transits for body $(ec.body_index)")
end

Base.show(io::IO, ::MIME"text/plain", ecs::Vector{<:ExpansionChain}) = begin
    print(io, "Vector{ExpansionChain}")
    #print(io, "$(length(ec.expansions)) transits for body $(ec.body_index)")
end

get_expansion(chain::ExpansionChain, i) = chain.expansions[i]
get_times(chain::ExpansionChain) = chain.times

"""Hold the different expansions needed to compute the impact parameter"""
struct ExpansionModel{ET, DT<:Differentiability, T<:Real} <: AbstractDynamicalModel
    chains::ET
    times::Vector{T}
    bodies::Vector{Int64}
    ic::InitialConditions{T}

    function ExpansionModel(
            chains::AbstractVector{<:ExpansionChain{<:AbstractExpansion{T}}},
            times::Vector{T}, bodies::Vector{Int64},
            ic::InitialConditions{T}) where T<:Real

        ET = typeof(chains) # Nested vector types (maybe simplify this?)
        DT = typeof(Differentiability(eltype(chains))) # Every index should be the same type...
        return new{ET, DT, T}(chains, times, bodies, ic)
    end
end
Differentiability(::ExpansionModel) = NonDifferentiable()

"""Setup the dynamical model; Pre-compute the expansions"""
function ExpansionModel(ic::InitialConditions{T}, tmax::T; h=zero(T)) where T<:Real
    # Compute the nbody model and output points
    s = State(ic)
    ts = TransitSeries(tmax, ic)
    tt = TransitTiming(tmax, ic)
    if h == 0.0; h = ic.elements[2,2]/30; end  ## If needed, get a reasonable stepsize
    Integrator(h, tmax)(s, ts, tt, grad=false) ## Fix typing for gradient computation ##

    # Compute components for each transit
    expansions = compute_expansion_chains(PK20Expansion, ts)

    return ExpansionModel(expansions, ts.times, ts.bodies, ic)
end

Base.show(io::IO, ::MIME"text/plain", em::ExpansionModel) = begin
    println(io, "$(Differentiability(em))"[1:end-2] * " $(typeof(em).name.name)")
    println(io, "Number of transits: $(length(em.times))")
    print(io, "Number of bodies: $(maximum(em.bodies))")
end

## Broadcast over the expansions inside of ExpansionModel
Base.broadcastable(em::ExpansionModel) = em.chains

## Accessor functions
get_times(dyn::ExpansionModel) = dyn.times
get_bodies(dyn::ExpansionModel) = dyn.bodies
get_expansion(dyn::ExpansionModel, i, j) = dyn.chains[j-1].expansions[i]

@copyfields struct PK20Expansion{V<:AbstractVector, T<:Real} <: AbstractExpansion{T}
    xc::V # Components
    yc::V # ""
    t0::T # Expansion time

    function PK20Expansion(xc::AbstractVector{T}, yc::AbstractVector{T}, t0::T) where T<:Real
        V = typeof(xc)
        return new{V,T}(xc, yc, t0)
    end
end
Differentiability(::Type{<:PK20Expansion}) = NonDifferentiable()

struct dPK20Expansion{V<:AbstractVector, VV<:AbstractVector, T<:Real} <: AbstractExpansion{T}
    @add_PK20Expansion_fields
    #xc::V # Components
    #yc::V # ""
    dxc::VV # Derivatives (Will be vector of vectors)
    dyc::VV
    #t0::T # Expansion time

    function dPK20Expansion(
            xc::AbstractVector{T}, yc::AbstractVector{T},
            dxc::AbstractVector{<:AbstractVector{T}}, dyc::AbstractVector{<:AbstractVector{T}},
            t0::T) where T<:Real

            V = typeof(xc)
            VV = typeof(dxc)
        return new{V,VV,T}(xc, yc, dxc, dyc, t0)
    end
end
Differentiability(::Type{<:dPK20Expansion}) = Differentiable()

function compute_expansions(::Type{PK20Expansion}, ts::TransitSeries{T}) where T<:Real
    expansions = sizehint!(PK20Expansion{SVector{5, T}, T}[], length(ts.times))
    for it in eachindex(ts.times)
        t0 = ts.times[it]
        ib = ts.bodies[it]

        xc = components(@views(ts.points[ib,it,:,1]), ts.h)
        yc = components(@views(ts.points[ib,it,:,2]), ts.h)

        exp = PK20Expansion(xc, yc, t0)
        push!(expansions, exp)
    end
    return expansions
end

function compute_expansion_chains(::Type{PK20Expansion}, ts::TransitSeries{T}) where T<:Real
    expType = PK20Expansion{SVector{5, T}, T} # For type inference
    nbody = maximum(ts.bodies) - 1 # Minus the star
    expansion_chains = sizehint!(ExpansionChain{expType, T}[], nbody)
    # Loop over each body and transit time
    for ib in 2:nbody+1
        expansions = sizehint!(expType[], length(ts.times))
        for it in eachindex(ts.times)
            t0 = ts.times[it]

            xc = components(@views(ts.points[ib,it,:,1]), ts.h)
            yc = components(@views(ts.points[ib,it,:,2]), ts.h)

            # No transit of this body
            if all(xc .== 0.0) && all(yc .== 0.0); continue; end

            # Store in expansion type and push to vector
            exp = PK20Expansion(xc, yc, t0)
            push!(expansions, exp)
        end
        exp_chain = ExpansionChain(expansions, ib)
        push!(expansion_chains, exp_chain)
    end
    return expansion_chains
end


function compute_sky_position(tc::T, ex::PK20Expansion{V,T}) where {V,T<:Real}
    t = tc - ex.t0
    ts = SVector{5, T}(1.0, t, t*t, t*t*t, t*t*t*t)
    lx = dot(ex.xc, ts)
    ly = dot(ex.yc, ts)
    return lx, ly
end

@inline compute_impact_parameter(tc::T, ex::PK20Expansion{V,T}) where {V,T<:Real} = compute_impact_parameter(tc, ex.t0, ex.xc, ex.yc)
@inline compute_impact_parameter(tc::Real, ex::PK20Expansion{V, T}) where {V,T<:Real} = compute_impact_parameter(T(tc), ex)

function compute_impact_parameter(tc::T, ex::dPK20Expansion{V,VV,T}) where {V,VV,T<:Real}
    t = tc - ex.t0
    ts = SVector(1.0, t, t*t, t*t*t, t*t*t*t)
    lx = dot(ex.xc, ts)
    ly = dot(ex.yc, ts)
    b = sqrt(lx*lx + ly*ly)

    dbdlx = b == 0 ? zero(T) : lx/b; dbdly = b==0 ? zero(T) : ly/b;
    dlxdp = dot.(ex.dxc, Ref(ts))
    dlydp = dot.(ex.dyc, Ref(ts))
    grad = dbdlx*dlxdp + dbdly*dlydp
    return b, grad
end

@inline function compute_impact_parameter(em::ExpansionModel{ET, NonDifferentiable, T}, tc::T) where {ET,T<:Real}
    return compute_impact_parameter.(em, tc)
end

function compute_impact_parameter(em::ExpansionModel{ET, Differentiable, T}, tc::T) where {ET,T<:Real}
    res = compute_impact_parameter.(em, tc)
    bs = first.(res)
    grads = last.(res)
    return bs, grads
end

#=function test_expansion()
    # Get a set of Components assuming linear trajectory
    dt = 0.02
    v0 = 0.001
    xps = SVector(collect(-3*v0*dt:dt*v0:3*v0*dt)...)
    yps = @SVector zeros(7)
    xc = Photodynamics.components(xps, dt)
    yc = Photodynamics.components(yps, dt)

    # Do random gradients for now
    dxc = @SVector [@SVector randn(5) for _ in 1:6]
    dyc = @SVector [@SVector randn(5) for _ in 1:6]

    t0 = 4*dt
    ex = PK20Expansion(xc, yc, t0)
    em = ExpansionModel(SVector(ex,ex))

    dex = dPK20Expansion(xc, yc, dxc, dyc, t0)
    dem = ExpansionModel(SVector(dex,dex))

    # Non gradient version
    tc = 6*dt
    bs = compute_impact_parameter.(em, tc)

    dbs, grads = compute_impact_parameter(dem, tc)

    @assert all(bs .== dbs) "Impact parameters differ"
    return bs, grads
end

test_expansion()=#