import Photodynamics: transit_poly_g

struct Lightcurve{T<:Real}
    dt::T     # Exposure time
    tobs::Vector{T}   # Observed times
    fobs::Vector{T}   # Observed flux
    eobs::Vector{T}   # Measurement errors
    nobs::Int64       # number of flux measurements
    flux::Vector{T}   # Computed model flux

    # Transit parameters
    u_n::Vector{T} # Limbdark coefficients
    k::Vector{T}   # radius ratios

    # Interal arrays/values
    dtinv::T

    function Lightcurve(dt::T, tobs::Vector{T}, fobs::Vector{T}, eobs::Vector{T}, u_n::Vector{T}, k::Vector{T}) where T<:Real
        @assert (length(tobs) == length(fobs)) && (length(tobs) == length(eobs)) "Data arrays are different sizes"
        nobs = length(tobs)
        flux = zeros(T,nobs)
        dtinv = inv(dt)
        return new{T}(dt,tobs,fobs,eobs,nobs,flux,u_n,k,dtinv)
    end
end

function points_of_contact_4(t0::T,h::T,points::AbstractMatrix{T},k::T) where T<:Real
    t1 = find_zero(t -> (1.0+k-compute_impact_parameter(t,t0,h,points)), t0-h)
    t2 = find_zero(t -> (1.0-k-compute_impact_parameter(t,t0,h,points)), t0-h)
    t3 = find_zero(t -> (1.0-k-compute_impact_parameter(t,t0,h,points)), t0+h)
    t4 = find_zero(t -> (1.0+k-compute_impact_parameter(t,t0,h,points)), t0+h)
    return SVector(t1,t2,t3,t4)
end

function points_of_contact_2(t0::T,h::T,points::AbstractMatrix{T},k::T) where T<:Real
    t1 = find_zero(t -> (1.0+k-compute_impact_parameter(t,t0,h,points)), t0-h)
    t4 = find_zero(t -> (1.0+k-compute_impact_parameter(t,t0,h,points)), t0+h)
    return SVector(t1,t4)
end

function compute_lightcurve!(lc::Lightcurve{T}, ts::TransitSeries{T}; tol::T=1e-6, maxdepth::Int64=6) where T<:Real

    ia = IntegralArrays(1, maxdepth, tol) # Gradients not implemented yet
    lc.flux .= 0.0 # Zero out model flux

    # Make transit structure (will be updated with proper r and b later)
    trans = transit_init(lc.k[1], 0.0, lc.u_n, false)

    # Iterate over each transit time and sum Lightcurve
    for it in eachindex(ts.times)
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
        trans.r = lc.k[ib-1]
        integrate_transit!(ib,it,t0,tc,trans,lc,ts,ia)
    end
    lc.flux .*= lc.dtinv # Divide by exposure time to get average flux
    return
end

function integrate_transit!(ib::Int64,it::Int64,t0::T,tc::SVector{N,T},trans::Transit_Struct{T},lc::Lightcurve{T},ts::TransitSeries{T},ia::IntegralArrays{T}) where {N, T<:Real}
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

        # Get series expansion components
        xc = components(@views(ts.points[ib,it,:,1]), ts.h)
        yc = components(@views(ts.points[ib,it,:,2]), ts.h)
        # Integrate over exposure
        for j in 1:length(tlim)-1
            integrate_timestep!(t0, tlim[j], tlim[j+1], xc, yc, trans, ia)
            lc.flux[i] += ia.I_of_x[1]
        end
    end
    return
end

function integrate_timestep!(t0::T, a::T, b::T, xc::SVector{N,T}, yc::SVector{N,T}, trans::Transit_Struct{T}, ia::IntegralArrays{T}) where {N, T<:Real}
    # Computes the flux as function of time.
    # Closure to be passed to integrator function
    transit_flux! = let trans=trans, t0=t0, xc=xc, yc=yc
        (time::T, flux::Vector{T}) -> begin
            trans.b = compute_impact_parameter(time, t0, yc, xc)
            flux[1] = transit_poly_g(trans)-1
        end
    end

    # Integrate transit_flux! over interval [a,b]
    integrate_simpson!(a,b,transit_flux!,ia)
end