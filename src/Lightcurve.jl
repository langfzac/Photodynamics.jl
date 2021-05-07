import Limbdark: transit_poly_g

struct Lightcurve
    dt     # Exposure time
    tobs   # Observed times
    fobs   # Observed flux
    eobs   # Measurement errors
    nobs   # number of flux measurements
    flux   # Computed model flux

    # Transit parameters
    u_n # Limbdark coefficients
    k   # radius ratios

    # Interal arrays/values
    dtinv

    function Lightcurve(dt, tobs, fobs, eobs, u_n, k)
        @assert (length(tobs) == length(fobs)) && (length(tobs) == length(eobs)) "Data arrays are different sizes"
        nobs = length(tobs)
        flux = zeros(nobs)
        dtinv = inv(dt)
        return new(dt,tobs,fobs,eobs,nobs,flux,u_n,k,dtinv)
    end
end

function points_of_contact_4(t0,h,points,k)
    t1 = find_zero(t -> (1.0+k-compute_impact_parameter(t,t0,h,points)), t0-h)
    t2 = find_zero(t -> (1.0-k-compute_impact_parameter(t,t0,h,points)), t0-h)
    t3 = find_zero(t -> (1.0-k-compute_impact_parameter(t,t0,h,points)), t0+h)
    t4 = find_zero(t -> (1.0+k-compute_impact_parameter(t,t0,h,points)), t0+h)
    return [t1,t2,t3,t4]
end

function points_of_contact_2(t0,h,points,k)
    t1 = find_zero(t -> (1.0+k-compute_impact_parameter(t,t0,h,points)), t0-h)
    t4 = find_zero(t -> (1.0+k-compute_impact_parameter(t,t0,h,points)), t0+h)
    return [t1,t4]
end

function compute_lightcurve!(lc, ts; tol=1e-6, maxdepth=6)

    ia = IntegralArrays(1, maxdepth, tol) # Gradients not implemented yet
    lc.flux .= 0.0 # Zero out model flux

    # Iterate over each transit time and sum Lightcurve
    for it in eachindex(ts.times[:,1])
        # check for transit
        t0 = ts.times[it]   # Get transit time
        ib = ts.bodies[it]  # Get transiting body
        b0 = compute_impact_parameter(t0,t0,ts.h, @views(ts.points[ib,it,:,:]))
        if b0 > 1.0+lc.k[ib-1]; continue; end

        # Compute points of contact
        # If grazing transit, only two points
        if (b0 + lc.k[ib-1]) > 1.0
            tc = points_of_contact_2(t0, ts.h, @views(ts.points[ib,it,:,:]), lc.k[ib-1])
        else
            tc = points_of_contact_4(t0, ts.h, @views(ts.points[ib,it,:,:]), lc.k[ib-1])
        end

        # Integrate lightcurve
        trans = transit_init(lc.k[ib-1], b0, lc.u_n, false)
        integrate_transit!(ib,it,t0,tc,trans,lc,ts,ia)
    end
    lc.flux .*= lc.dtinv # Divide by exposure time to get average flux
    return
end

function integrate_transit!(ib,it,t0,tc,trans,lc,ts,ia)
    nc = length(tc) # Number of points of contact
    dt = lc.dt

    # Integrate over each Exposure
    for i in 1:lc.nobs
        tstart = lc.tobs[i] - 0.5*dt
        tend = lc.tobs[i] + 0.5*dt

        # Check if remaining exposures are inside transit
        if tstart > tc[end]; break; end # Don't need to continue loop if passed transit
        if tend < tc[1]
            # flux == 0
            continue
        end

        # Check if points of contact are within exposure
        tlim = [tstart]
        for j in 1:nc
            if tstart < tc[j] && tc[j] < tend
                push!(tlim, tc[j])
            end
        end
        push!(tlim, tend)

        # Integrate over exposure
        for j in 1:length(tlim)-1
            integrate_timestep!(t0, tlim[j], tlim[j+1], ts.h, @views(ts.points[ib,it,:,:]), trans, ia)
            lc.flux[i] += ia.I_of_x[1]
        end
    end
    return
end

function integrate_timestep!(t0::T, a, b, h, points, trans, ia) where T<:AbstractFloat
    # Computes the flux as function of time.
    # Closure to be passed to integrator function
    transit_flux! = let trans=trans, t0=t0, points=points, h=h
        (time::T, flux::Vector{T}) -> begin
            trans.b = compute_impact_parameter(time, t0, h, points)
            flux[1] = transit_poly_g(trans)-1
        end
    end

    # Integrate transit_flux! over interval [a,b]
    integrate_simpson!(a,b,transit_flux!,ia)
end