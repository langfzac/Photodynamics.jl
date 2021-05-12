# NbodyGradient integration to compute points about each transit for series expansion.
import NbodyGradient: TransitOutput, InitialConditions
import NbodyGradient: check_step, set_state!

struct TransitSeries{T<:AbstractFloat} <: TransitOutput{T}
    times::Vector{T}    # Transit times, sequentially
    bodies::Vector{Int64}   # Transiting body at the ith transit time
    points::Array{T}    # [body, transit, 7, 2]

    # Internal values/arrays
    h::T                 # Stepsize for series points
    ntt::Int64           # Total number of transits
    count::Vector{Int64} # Number of transits for ith body
    s_transit::State{T}  # Duplicate state to do sub-integrations
end

function TransitSeries(times::Matrix{T}, ic::InitialConditions{T}; h=2e-2) where T<:AbstractFloat
    ntt = sum(times .> ic.t0) # Total transits, assumes no transits <= ic.t0

    # Flatten transit times and assign a body index
    count = zeros(Int64, ic.nbody)
    tmp = zeros(T,ntt,2)
    for (i,t) in enumerate(eachrow(times))
        for j in eachindex(t)
            tc = t[j]
            if tc > ic.t0
                count[i] += 1
                ind = sum(count)
                tmp[ind,1] = tc
                tmp[ind,2] = i
            end
        end
    end

    # Now sort by transit time
    tmp .= sortslices(tmp, dims=1)

    times = copy(tmp[:,1])
    bodies = copy(Int64.(tmp[:,2]))

    s_transit = State(ic)
    points = zeros(T,ic.nbody,ntt,7,2)
    TransitSeries(times, bodies, points, convert(T,h), ntt, count, s_transit)
end

"""Integrate to each transit time and record 7 points for series expansion."""
function (intr::Integrator)(s::State{T},ts::TransitSeries{T}; grad::Bool=false) where T<:AbstractFloat
    t0 = s.t[1]                             # Initial time
    tmax = maximum(ts.times)                # Final transit time
    nsteps = abs(round(Int64, tmax/intr.h)) # Total integration steps
    h = intr.h * check_step(t0, tmax+t0)    # Make sure we step in proper direction

    # Run integrator, check for transits against list of times, record sky positions
    istep = 0; itime = 1;
    for _ in 1:nsteps
        istep += 1
        tc = ts.times[itime] # Current time
        tnext = t0 + (istep * h) # Time after next step

        # Check if we will pass a transit
        if tnext >= tc
            # Save current state
            set_state!(ts.s_transit, s)

            # Integrate forward to transit time
            htrans = tc - s.t[1]
            intr.scheme(ts.s_transit, htrans)

            # Compute sky-plane positions for each body
            compute_points!(itime, ts, intr)

            itime += 1
        end

        # Now, step saved state.
        intr.scheme(s,h)
        s.t[1] = tnext

        # Break if no more transits
        if itime > ts.ntt; break; end;
    end
end

function compute_points!(itime::Int64,ts::TransitSeries{T},intr::Integrator{T}) where T<:AbstractFloat
    # Want to compute positions from transit-3h -> transit+3h
    # Integrate backward from the transit time
    for _ in 1:4; intr.scheme(ts.s_transit, -ts.h); end

    # Compute 7 sky-positions about the transit
    ibody = ts.bodies[itime]
    for n in 1:7
        intr.scheme(ts.s_transit, ts.h)
        @views ts.points[ibody,itime,n,:] .= ts.s_transit.x[1:2,ibody]
    end
end