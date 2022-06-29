# NbodyGradient integration to compute points about each transit for series expansion.

abstract type AbstractTransitTimes end
struct ComputedTimes <: AbstractTransitTimes end
struct ProvidedTimes <: AbstractTransitTimes end

struct TransitSeries{T<:Real, TT<:AbstractTransitTimes}
    times::Vector{T}      # Transit times, sequentially
    bodies::Vector{Int64} # Transiting body at the ith transit time
    points::Array{T, 4}   # [body, transit, 7, 2]
    dpoints::Array{T, 6}  # [body, transit, 7, 2, body, parameter]

    # Internal values/arrays
    h::T                  # Stepsize for series points
    ntt::Int64            # Total number of transits
    intr_times::Vector{T} # Times to integrate to anc compute points
    count::Vector{Int64}  # Number of transits for ith body
    s_prior::Vector{State{T}}
end

"""Pass a set of transit times"""
function TransitSeries(times::Matrix{T}, ic::InitialConditions{T}; h::T=T(2e-2)) where T<:Real
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

    # Compute the time of the first expansion point for each transit
    intr_times = zeros(T,length(times))
    for (i,t) in enumerate(times)
        intr_times[i] =  t - 3*h
    end

    points = zeros(T,ic.nbody,ntt,7,2)
    dpoints = zeros(T,ic.nbody,ntt,7,2,ic.nbody,7)
    s_prior = [State(ic)]
    TransitSeries{T, ProvidedTimes}(times, bodies, points, dpoints, h, ntt, intr_times, count, s_prior)
end

"""Need to compute transit times"""
function TransitSeries(tmax::T, ic::InitialConditions{T}; h::T=T(2e-2)) where T<:Real
    ntt::Int64 = 0  # Total expected transits, assumes no transits <= ic.t0
    for P in ic.elements[:,2]
        if P > 0.0 # Allow for orbital elements to be 'out-of-order'
            ntt += ceil(Int64, tmax/P) # Likely to have a buffer
        end
    end

    times = zeros(T, ntt)
    bodies = zeros(Int64, ntt)
    intr_times = zeros(T,length(times))
    points = zeros(T,ic.nbody,ntt,7,2)
    dpoints = zeros(T,ic.nbody,ntt,7,2,ic.nbody,7)
    count = zeros(Int64, ic.nbody)
    s_prior = [State(ic) for _ in 1:5]
    TransitSeries{T, ComputedTimes}(times, bodies, points, dpoints, h, ntt, intr_times, count, s_prior)
end

"""Compute 7 points about each transit time.

Integrate to right before the first transit expansion point. Save state and
integrate for 7 steps at h=ts.h. Revert to pre-transit state and continue to
next transit.
"""
function (intr::Integrator)(s::State{T},ts::TransitSeries{T, ProvidedTimes},d::Union{Derivatives,Nothing}=nothing; grad::Bool=false) where T<:Real
    if grad && d == nothing; d = Derivatives(T, s.n); end
    if ~grad && d !== nothing; grad = true; end

    # Run integrator and record sky positions for list of integration times
    nstep::Int64 = 0; t0 = s.t[1]
    s_prior = ts.s_prior[1]
    for (itime,t) in enumerate(ts.intr_times)
        ibody = ts.bodies[itime]

        # Integrate to right before the first point in the transit
        while s.t[1] < t
            # Break if we will pass a transit on this step
            if t-s.t[1] < intr.h; break; end

            if grad
                intr.scheme(s,d,intr.h)
            else
                intr.scheme(s,intr.h)
            end
            nstep += 1
            s.t[1] = t0 + nstep*intr.h
        end

        # Save the pre-transit state
        set_state!(s_prior, s)

        # Now integrate to each expansion point with step ts.h
        if grad
            compute_points!(s, ts, d, t, s_prior.t[1], ts.h, ibody, itime, intr)
        else
            compute_points!(s, ts, t, s_prior.t[1], ts.h, ibody, itime, intr)
        end

        # Return to state before transit
        set_state!(s, s_prior)
    end
end

function compute_points!(s, ts, t, t0, h, ibody, itime, intr)
    for k in 1:7
        # Step to the first point if needed
        k == 1 ? h_points = t-t0 : h_points = h
        intr.scheme(s, h_points)

        # record points for transiting body
        # Assumes transiting single star
        ts.points[ibody,itime,k,1] = s.x[1,ibody] - s.x[1,1]
        ts.points[ibody,itime,k,2] = s.x[2,ibody] - s.x[2,1]
    end
    return
end

function compute_points!(s, ts, d, t, t0, h, ibody, itime, intr)
    for k in 1:7
        # Step to the first point if needed
        k == 1 ? h_points = t-t0 : h_points = h
        intr.scheme(s, d, h_points)

        # record points for transiting body
        # Assumes transiting single star
        ts.points[ibody,itime,k,1] = s.x[1,ibody] - s.x[1,1]
        ts.points[ibody,itime,k,2] = s.x[2,ibody] - s.x[2,1]

        # Get gradients of points wrt initial orbital elements and masses
        # p body, q element
        for p in 1:s.n, q in 1:7
            ts.dpoints[ibody,itime,k,1,p,q] = s.jac_step[(ibody-1)*7+1,(p-1)*7+q] - s.jac_step[1, (p-1)*7+q]
            ts.dpoints[ibody,itime,k,2,p,q] = s.jac_step[(ibody-1)*7+2,(p-1)*7+q] - s.jac_step[2, (p-1)*7+q]
        end
    end
    return
end

"""Compute the transit times and 7 points about each.
"""
function (intr::Integrator)(s::State{T},ts::TransitSeries{T, ComputedTimes},tt::TransitOutput{T},d::Union{Derivatives,Nothing}=nothing; grad::Bool=false) where T<:Real
    if d isa Nothing; d = Derivatives(T, s.n); end
    nstates = length(ts.s_prior)

    h = intr.h; t0 = s.t[1]
    nsteps = abs(round(Int64, intr.tmax/intr.h))

    # Save the sky-acceleration at the initial conditions
    for i in tt.occs
        tt.gsave[i] = NbodyGradient.g!(i,tt.ti,s.x,s.v)
    end

    # Save initial state to each element of vector
    set_state!.(ts.s_prior, Ref(s))

    # Set a counter for which state in s_prior is nstates steps before
    state_counter = (1,1)  # (current state, furthest state)

    for i in 1:nsteps
        # Take a step with the integrator
        if grad
            intr.scheme(s,d,h)
        else
            intr.scheme(s,h)
        end
        s.t[1] = t0 + (i * h)

        # Save the current number of transits and check if any new ones occured
        prior_count = copy(tt.count)
        NbodyGradient.detect_transits!(s,d,tt,intr,grad=grad)

        # Check whether a transit did occur.
        # If so, compute 7 points around transit time.
        # This assumes the integration timestep is larger than the
        # expansion resolution (ie. ts.dt).
        count_mask = prior_count .< tt.count
        ibodies = findall(count_mask) # Indices of bodies which transited

        # Make sure the transit times are chronological
        times = @SVector T[]
        for (i,j) in zip(ibodies, @view(tt.count[ibodies]))
            times = push(times, tt.tt[i,j])
        end
        inds = sortperm(times)
        ibodies = ibodies[inds]

        # If theres more than one transit we need to revert the state for each
        s_points = ts.s_prior[last(state_counter)]
        if length(ibodies) > 1; set_state!(tt.s_prior, s_points); end

        # Loop over each transit and compute expansion points
        for (j,ibody) in enumerate(ibodies)
            # If more than one transit occured, reset state to pre-transit
            if j > 1; set_state!(s_points, tt.s_prior); end

            # Save time and body index to TransitSeries
            itime = sum(tt.count) - (length(ibodies) - j)
            ts.times[itime] = tt.tt[ibody, tt.count[ibody]]
            ts.bodies[itime] = ibody

            # Get the time of the first expansion point
            t_first = ts.times[itime] - 3*ts.h

            nstep = 0; t0_points = s_points.t[1]
            while s_points.t[1] < t_first
                t_current = s_points.t[1]
                # Break if we will pass the point on this step
                if t_first - t_current < h; break; end

                if grad
                    intr.scheme(s_points, d, h)
                else
                    intr.scheme(s_points, h)
                end
                nstep += 1
                s_points.t[1] = t0_points + (nstep * h)
            end

            # Save the pre-point state
            set_state!(tt.s_prior, s_points)

            if grad
                compute_points!(s_points, ts, d, t_first, s_points.t[1], ts.h, ibody, itime, intr)
            else
                compute_points!(s_points, ts, t_first, s_points.t[1], ts.h, ibody, itime, intr)
            end
        end

        # Shift counter by 1 and wrap at nstates
        state_counter = (
            mod1(first(state_counter) + 1, nstates),
            i < nstates ? 1 : mod1(i+2, nstates)
        )

        # reset/set to current state
        set_state!(s_points, tt.s_prior)
        set_state!(ts.s_prior[first(state_counter)], s)
    end

    # Clean up arrays/remove buffer
    while ts.bodies[end] == 0
        pop!(ts.bodies)
        pop!(ts.times)
    end
end
