# User-level methods to compute lightcurves

struct Lightcurve{T<:Real}
    dt::T             # Exposure time
    tobs::Vector{T}   # Observed times
    fobs::Vector{T}   # Observed flux
    eobs::Vector{T}   # Measurement errors
    nobs::Int64       # number of flux measurements
    flux::Vector{T}   # Computed model flux
    dfdu::Matrix{T}   # Derivative of flux wrt limbdark coefficients
    dfdk::Matrix{T}   # Derivatives wrt the radius ratios
    dfdq0::Matrix{T}  # Derivatives wrt initial Nbody Cartesian coordinates and masses
    dfdr::Vector{T}   # Derivatives wrt stellar radius

    # Transit parameters
    u_n::Vector{T} # Limbdark coefficients
    k::Vector{T}   # radius ratios
    rstar::T       # Stellar radius

    # Interal arrays/values
    dtinv::T          # Inverse of the exposure time
    dbdq0::Vector{T}  # derivative of the impact parameter wrt the Nbody initial conditions
    n_params::Int64   # Number of nbody model parameters
    do_grad::Bool     # Compute gradients or not

    function Lightcurve(dt::T, tobs::Vector{T}, fobs::Vector{T}, eobs::Vector{T}, u_n::Vector{T}, k::Vector{T}, rstar::T, n_params::Int64=0) where T<:Real
        @assert (length(tobs) == length(fobs)) && (length(tobs) == length(eobs)) "Data arrays are different sizes"
        @assert n_params >= 0 "Number of model parameters must be a positive integer"
        n_params > 0 ? do_grad = true : do_grad = false
        nobs = length(tobs)
        flux = zeros(T,nobs)
        dtinv = inv(dt)
        dfdu = zeros(T, nobs, length(u_n))
        dfdk = zeros(T, nobs, length(k))
        dfdq0 = zeros(T, nobs, n_params)
        dbdq0 = zeros(T, n_params)
        dfdr = zeros(T, nobs)
        return new{T}(dt,tobs,fobs,eobs,nobs,flux,dfdu,dfdk,dfdq0,dfdr,u_n,k,rstar,dtinv,dbdq0,n_params,do_grad)
    end
end

"""Zero out the model arrays"""
function zero_out!(lc::Lightcurve{T}) where T<:Real
    lc.dfdu .= 0.0
    lc.dfdk .= 0.0
    lc.dfdq0 .= 0.0
    lc.dfdr .= 0.0
    lc.flux .= 0.0
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

    zero_out!(lc) # Zero out model arrays
    ia = IntegralArrays(lc.do_grad ? (lc.n_params + length(lc.u_n) + length(lc.k) + 1) : 1, maxdepth, tol)

    # Make transit structure (will be updated with proper r and b later)
    trans = transit_init(lc.k[1], 0.0, lc.u_n, lc.do_grad)

    # Iterate over each transit time and sum Lightcurve
    for it in eachindex(ts.times)
        # check for transit
        t0 = ts.times[it]   # Get transit time
        ib = ts.bodies[it]  # Get transiting body
        b0 = compute_impact_parameter(t0,t0,ts.h,@views(ts.points[ib,it,:,:]./lc.rstar))
        if b0 > 1.0+lc.k[ib-1]; continue; end

        # Compute points of contact
        # If grazing transit, only two points
        if (b0 + lc.k[ib-1]) > 1.0
            tc = points_of_contact_2(t0, ts.h, @views(ts.points[ib,it,:,:]./lc.rstar), lc.k[ib-1])
        else
            tc = points_of_contact_4(t0, ts.h, @views(ts.points[ib,it,:,:]./lc.rstar), lc.k[ib-1])
        end

        # Integrate lightcurve
        trans.r = lc.k[ib-1]
        integrate_transit!(ib,it,t0,tc,trans,lc,ts,ia)
    end
    lc.flux .*= lc.dtinv # Divide by exposure time to get average flux
    if lc.do_grad
        # Do the same for derivatives
        lc.dfdk .*= lc.dtinv
        lc.dfdu .*= lc.dtinv
        lc.dfdq0 .*= lc.dtinv
        lc.dfdr .*= lc.dtinv
    end
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
        if tend < tc[1]; continue; end

        # Check if points of contact are within exposure
        tlim = [tstart]
        for j in 1:nc
            if tstart < tc[j] && tc[j] < tend
                push!(tlim, tc[j])
            end
        end
        push!(tlim, tend)

        # Get series expansion components
        xc = components(@views(ts.points[ib,it,:,1]./lc.rstar), ts.h)
        yc = components(@views(ts.points[ib,it,:,2]./lc.rstar), ts.h)

        if lc.do_grad
            n_bodies = length(ts.count)
            dxc = [components(ts.dpoints[ib,it,:,1,k,i]./lc.rstar, ts.h) for i in 1:7, k in 1:n_bodies][:]
            dyc = [components(ts.dpoints[ib,it,:,2,k,i]./lc.rstar, ts.h) for i in 1:7, k in 1:n_bodies][:]

            # integrate over exposure
            for j in 1:length(tlim)-1
                integrate_timestep!(t0, tlim[j], tlim[j+1], xc, yc, dxc, dyc, trans, ia, lc.dbdq0, ib-1)
                lc.flux[i] += ia.I_of_f[1]
                lc.dfdq0[i,:] .+= ia.I_of_f[2:1+n_bodies*7]
                lc.dfdk[i,:] .+= ia.I_of_f[2+n_bodies*7:n_bodies*8]
                lc.dfdu[i,:] .+= trans.dgdu' * ia.I_of_f[end-trans.n-1:end-1]
                lc.dfdr[i] += ia.I_of_f[end] / lc.rstar
            end
        else
            # Integrate over exposure
            for j in 1:length(tlim)-1
                integrate_timestep!(t0, tlim[j], tlim[j+1], xc, yc, trans, ia)
                lc.flux[i] += ia.I_of_f[1]
            end
        end
    end
    return
end

function integrate_timestep!(t0::T, a::T, b::T, xc::SVector{N,T}, yc::SVector{N,T}, trans::Transit_Struct{T}, ia::IntegralArrays{T}) where {N, T<:Real}
    # Computes the flux as function of time.
    # Closure to be passed to integrator function
    transit_flux! = let trans=trans, t0=t0, xc=xc, yc=yc
        (time::T, flux::Vector{T}) -> begin
            trans.b = compute_impact_parameter(time, t0, xc, yc)
            flux[1] = transit_poly_g(trans)-1
        end
    end

    # Integrate transit_flux! over interval [a,b]
    integrate_simpson!(a,b,transit_flux!,ia)
end

function integrate_timestep!(t0::T, a::T, b::T, xc, yc, dxc, dyc, trans, ia, dbdq0, ki) where T<:Real
    n_coords = length(dbdq0)

    # Map radius ratio index to the gradient array index.
    k_ind = ki + 1 + n_coords

    # Compute the flux and derivatives
    transit_flux_grad! = let trans=trans, t0=t0, xc=xc, yc=yc, dbdq0=dbdq0, n_coords = n_coords, k_ind=k_ind
        (time::T, dflux::Vector{T}) -> begin
            trans.b = compute_impact_parameter!(time, t0, xc, yc, dxc, dyc, dbdq0)
            dflux .= 0.0
            dflux[1] = transit_poly_g!(trans)-1            # Flux at time
            dflux[2:1+n_coords] .= trans.dfdrb[2] .* dbdq0 # Gradient wrt the initial conditions
            dflux[k_ind] = trans.dfdrb[1]                  # radius ratio (derivative of others are zero)
            dflux[end-trans.n-1:end-1] .= trans.dfdg       # dfdg (greens basis)
            dflux[end] = -trans.dfdrb[2]*trans.b           # dfdrstar * rstar (divide by rstar outside)
        end
    end

    # Integrate flux and derivatives over exposure interval [a,b]
    integrate_simpson!(a,b,transit_flux_grad!,ia)
end