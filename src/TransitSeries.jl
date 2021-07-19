# NbodyGradient integration to compute points about each transit for series expansion.

struct TransitSeries{T<:Real}
    times::Vector{T}      # Transit times, sequentially
    bodies::Vector{Int64} # Transiting body at the ith transit time
    points::Array{T, 4}      # [body, transit, 7, 2]
    dpoints::Array{T, 6}     # [body, transit, 7, 2, body, parameter]

    # Internal values/arrays
    h::T                  # Stepsize for series points
    ntt::Int64            # Total number of transits
    intr_times::Vector{T} # Times to integrate to anc compute points
    count::Vector{Int64}  # Number of transits for ith body
    s_prior::State{T}
end

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
    s_prior = State(ic)
    TransitSeries(times, bodies, points, dpoints, h, ntt, intr_times, count, s_prior)
end

"""Compute 7 points about each transit time.

Integrate to right before the first transit expansion point. Save state and
integrate for 7 steps at h=ts.h. Revert to pre-transit state and continue to
next transit.
"""
function (intr::Integrator)(s::State{T},ts::TransitSeries{T}; grad::Bool=false) where T<:Real
    if grad; d = Derivatives(T, s.n); end

    # Run integrator and record sky positions for list of integration times
    nstep::Int64 = 0; t0 = s.t[1]
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
        set_state!(ts.s_prior, s)

        # Now integrate to each expansion point with step ts.h
        for k in 1:7

            # Step to the first point if needed
            k == 1 ? h_points = t-ts.s_prior.t[1] : h_points = ts.h

            if grad
                intr.scheme(s,d,h_points)
            else
                intr.scheme(s,h_points)
            end

            # record points for transiting body
            ts.points[ibody,itime,k,1] = s.x[1,ibody]
            ts.points[ibody,itime,k,2] = s.x[2,ibody]

            if grad
                # Get gradients of points wrt initial orbital elements and masses
                # p body, q element
                for p in 1:s.n, q in 1:7
                    ts.dpoints[ibody,itime,k,1,p,q] = s.jac_step[(ibody-1)*7+1,(p-1)*7+q]
                    ts.dpoints[ibody,itime,k,2,p,q] = s.jac_step[(ibody-1)*7+2,(p-1)*7+q]
                end
            end
        end

        # Return to state before transit
        set_state!(s, ts.s_prior)
    end
end
