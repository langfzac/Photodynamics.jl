# NbodyGradient integration to compute points about each transit for series expansion.
import NbodyGradient: TransitOutput, InitialConditions, Derivatives
import NbodyGradient: check_step, set_state!

struct TransitSeries{T<:AbstractFloat} <: TransitOutput{T}
    times::Vector{T}      # Transit times, sequentially
    bodies::Vector{Int64} # Transiting body at the ith transit time
    points::Array{T}      # [body, transit, 7, 2]
    dpoints::Array{T}     # [body, transit, 7, 2, parameter]

    # Internal values/arrays
    h::T                  # Stepsize for series points
    ntt::Int64            # Total number of transits
    intr_times::Vector{T} # Times to integrate to anc compute points
    count::Vector{Int64}  # Number of transits for ith body
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

    # Now compute the integration times
    intr_times = T[]
    for t in times
        for i in -3:1:3
            push!(intr_times, t + i*h)
        end
    end

    points = zeros(T,ic.nbody,ntt,7,2)
    dpoints = zeros(T,ic.nbody,ntt,7,2,ic.nbody,7)
    TransitSeries(times, bodies, points, dpoints, convert(T,h), ntt, intr_times, count)
end

"""Integrate to each time (expanded from transit time) and record points for series expansion."""
function (intr::Integrator)(s::State{T},ts::TransitSeries{T}; grad::Bool=false) where T<:AbstractFloat

    # Run integrator and record sky positions for list of integration times
    count = 1; itime = 1; ibody = ts.bodies[itime]
    for t in ts.intr_times

        # Update itime and ibody
        if count > 7
            itime += 1
            ibody = ts.bodies[itime]
            count = 1
        end

        # Integrate to t
        intr(s,t,grad=grad)

        # record points for transiting body
        ts.points[ibody,itime,count,1] = s.x[1,ibody]
        ts.points[ibody,itime,count,2] = s.x[2,ibody]

        if grad
            # Get gradients of points wrt initial orbital elements and masses
            # p body, q element
            for p in 1:s.n, q in 1:7
                ts.dpoints[ibody,itime,count,1,p,q] = s.jac_step[(ibody-1)*7+1,(p-1)*7+q]
                ts.dpoints[ibody,itime,count,2,p,q] = s.jac_step[(ibody-1)*7+2,(p-1)*7+q]
            end
        end

        count += 1
    end
end
