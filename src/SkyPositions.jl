struct SkyPositions{T<:Real}
    times::Vector{T}  # Times at which positions are calculated
    points::Array{T, 3}  # [time, coordinate, body] (ntt X 3 X nbody)
    dpoints::Array{T, 5}  # [time, coordinate, body, body, parameter] (ntt X 3 X nbody X nbody X 7)

    # Internal
    s_prior::State{T}
end

function SkyPositions(times::Vector{T}, ic::InitialConditions{T}) where T<:Real
    points = zeros(T, length(times), 3, ic.nbody)
    dpoints = zeros(T, length(times), 3, ic.nbody, ic.nbody, 7)
    s_prior = State(ic)
    return SkyPositions{T}(times, points, dpoints, s_prior)
end

function (intr::Integrator)(s::State{T}, sp::SkyPositions{T}, d::Union{Derivatives,Nothing}=nothing; grad::Bool=false) where T<:Real
    @assert s.t[1] <= sp.times[1] "Integration must start before first exposure time"
    if grad && d == nothing; d = Derivatives(T, s.n); end
    if ~grad && d !== nothing; grad = true; end

    t0 = s.t[1]; h = intr.h
    nsteps = abs(round(Int64, (sp.times[end] + h - t0)/h))
    it = 1; tot = length(sp.times)
    for istep in 1:nsteps
        tc = s.t[1]
        # If we'll pass an exposure time, integrate to it instead
        while sp.times[it] < tc + h
            dt = sp.times[it] - tc
            set_state!(sp.s_prior, s)  # Save current state for after step
            intr.scheme(sp.s_prior, dt)  # Integrate to exposure time
            if grad
                save_positions!(it, sp.s_prior, sp.points, sp.dpoints)  # Save positions to structure
            else
                save_positions!(it, sp.s_prior, sp.points)
            end
            it += 1
            set_state!(sp.s_prior, s)
            # Stop if we're pass the data
            if it > tot; break; end
        end

        if grad
            intr.scheme(s,d,h)
        else
            intr.scheme(s,h)
        end
        s.t[1] = t0 + istep*intr.h
        # Stop if we're pass the data
        if it > tot; break; end
    end
end

@inline function save_positions!(it, s::State{T}, points::Array{T, 3}) where T<:Real
    points[it, :, :] .= s.x
end