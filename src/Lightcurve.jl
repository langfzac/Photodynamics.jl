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
    dfdelements::Matrix{T}

    # Transit parameters
    u_n::Vector{T}   # Limbdark coefficients
    k::Vector{T}     # radius ratios
    rstar::Vector{T} # Stellar radius

    # Integral arrays/values
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
        dtinv = dt == 0.0 ? 0.0 : inv(dt)
        dfdu = zeros(T, nobs, length(u_n))
        dfdk = zeros(T, nobs, length(k))
        dfdq0 = zeros(T, nobs, n_params)
        dbdq0 = zeros(T, n_params)
        dfdr = zeros(T, nobs)
        dfdelements = zeros(T, nobs, n_params)
        return new{T}(dt,copy(tobs),copy(fobs),copy(eobs),nobs,flux,dfdu,dfdk,dfdq0,dfdr,dfdelements,copy(u_n),copy(k),[rstar],dtinv,dbdq0,n_params,do_grad)
    end
end

"""Setup Lightcurve without having data"""
function Lightcurve(dt::T, duration::T, u_n::Vector{T}, k::Vector{T}, rstar::T, n_params::Int64=0) where T<:Real
    tobs = collect(0.0:dt:duration)
    return Lightcurve(dt, tobs, zeros(T, size(tobs)), zeros(T,size(tobs)), u_n, k, rstar, n_params)
end

"""Zero out the model arrays"""
function zero_out!(lc::Lightcurve{T}) where T<:Real
    lc.dfdu .= 0.0
    lc.dfdk .= 0.0
    lc.dfdq0 .= 0.0
    lc.dfdr .= 0.0
    lc.flux .= 0.0
end

# Normalize the exposure integration by the exposure time.
function normalize!(lc::Lightcurve{T}, ia::IntegralArrays{T}) where T<:Real
    lc.flux .*= lc.dtinv # Divide by exposure time to get average flux
    if lc.do_grad
        # Do the same for derivatives
        lc.dfdk .*= lc.dtinv
        lc.dfdu .*= lc.dtinv
        lc.dfdq0 .*= lc.dtinv
        lc.dfdr .*= lc.dtinv
    end
end

# If non-integrated exposure, do nothing
@inline normalize!(lc::Lightcurve{T}, ia::T) where T<:Real = nothing

function find_transit_time(t0::T, h::T, points::AbstractMatrix{T}) where T<:Real
    tt = find_zero(t->compute_impact_parameter(t, t0, h, points), t0)
    return tt
end

function points_of_contact_4(tt::T,t0::T,h::T,points::AbstractMatrix{T},k::T) where T<:Real
    t1 = find_zero(t -> (1.0+k-compute_impact_parameter(t,t0,h,points)), tt-h)
    t2 = find_zero(t -> (1.0-k-compute_impact_parameter(t,t0,h,points)), tt-h)
    t3 = find_zero(t -> (1.0-k-compute_impact_parameter(t,t0,h,points)), tt+h)
    t4 = find_zero(t -> (1.0+k-compute_impact_parameter(t,t0,h,points)), tt+h)
    return SVector{4, T}(t1,t2,t3,t4)
end

function points_of_contact_2(t0::T,tt::T,h::T,points::AbstractMatrix{T},k::T) where T<:Real
    t1 = find_zero(t -> (1.0+k-compute_impact_parameter(t,t0,h,points)), tt-h)
    t4 = find_zero(t -> (1.0+k-compute_impact_parameter(t,t0,h,points)), tt+h)
    return SVector{2, T}(t1,t4)
end

function compute_lightcurve!(lc::Lightcurve{T}, ts::TransitSeries{T, TT}; tol::T=1e-6, maxdepth::Int64=6, body_index::Int64=0) where {T<:Real, TT<:AbstractTransitTimes}

    zero_out!(lc) # Zero out model arrays

    # Check if we're doing an integrated lightcurve
    ia = lc.dt == 0.0 ? zero(T) : IntegralArrays(lc.do_grad ? (lc.n_params + length(lc.u_n) + length(lc.k) + 2) : 1, maxdepth, tol) # Plus 2 for flux and rstar

    # Make transit structure (will be updated with proper r and b later)
    trans = transit_init(lc.k[1], zero(T), lc.u_n, lc.do_grad)

    # Iterate over each transit time and sum Lightcurve
    rstar = lc.rstar[1]
    for it in eachindex(ts.times)
        # check for transit
        t0 = ts.times[it]   # Get "data" transit time (expansion point)
        ib = ts.bodies[it]  # Get transiting body
        if ib == 0; break; end  # Break if we've reached the end of the transits

        # Check whether we're computing flux for only a single body
        # If we are, skip all other body indices.
        if body_index != 0 && body_index != ib; continue; end

        # Get the impact parameter at the transit midpoint (computed above)
        b0 = compute_impact_parameter(t0, t0, ts.h, @views(ts.points[ib,it,:,:]./rstar))
        if b0 > 1.0+lc.k[ib-1]; continue; end

        # Compute points of contact
        # If grazing transit, only two points
        if (b0 + lc.k[ib-1]) >= 1.0
            tc = points_of_contact_2(t0, t0, ts.h, @views(ts.points[ib,it,:,:]./rstar), lc.k[ib-1])
        else
            tc = points_of_contact_4(t0, t0, ts.h, @views(ts.points[ib,it,:,:]./rstar), lc.k[ib-1])
        end

        # Integrate lightcurve
        trans.r = lc.k[ib-1]
        integrate_transit!(ib,it,t0,tc,trans,lc,ts,ia)
    end

    # Normalize by exposure time, if needed.
    normalize!(lc, ia)
    return
end

function integrate_transit!(ib::Int64,it::Int64,t0::T,tc::SVector{N,T},trans::Transit_Struct{T},lc::Lightcurve{T},ts::TransitSeries{T},ia::IntegralArrays{T}) where {N, T<:Real}
    # nc = length(tc) # Number of points of contact
    dt = lc.dt
    inv_rstar = inv(lc.rstar[1])

    # Integrate over each Exposure
    # TO-DO: Only iterate over the times around the transit.
    for i in eachindex(lc.tobs)
        tstart = lc.tobs[i] - 0.5*dt
        tend = lc.tobs[i] + 0.5*dt

        # Check if remaining exposures are inside transit
        if tstart > tc[end]; break; end # Don't need to continue loop if passed transit
        if tend < tc[1]; continue; end

        # Check if points of contact are within exposure
        tlim = [tstart]
        for j in eachindex(tc)
            if tstart < tc[j] && tc[j] < tend
                push!(tlim, tc[j])
            end
        end
        push!(tlim, tend)

        # Get series expansion components
        xc = components(@views(ts.points[ib,it,:,1].*inv_rstar), ts.h)
        yc = components(@views(ts.points[ib,it,:,2].*inv_rstar), ts.h)

        if lc.do_grad
            n_bodies = length(ts.count)
            dxc = [components(ts.dpoints[ib,it,:,1,k,i].*inv_rstar, ts.h) for i in 1:7, k in 1:n_bodies][:]
            dyc = [components(ts.dpoints[ib,it,:,2,k,i].*inv_rstar, ts.h) for i in 1:7, k in 1:n_bodies][:]

            # integrate over exposure
            for j in eachindex(@view(tlim[1:end-1]))
                integrate_timestep!(t0, tlim[j], tlim[j+1], xc, yc, dxc, dyc, trans, ia, lc.dbdq0, ib-1)
                lc.flux[i] += ia.I_of_f[1]
                lc.dfdq0[i,:] .+= ia.I_of_f[2:1+n_bodies*7]
                lc.dfdk[i,:] .+= ia.I_of_f[2+n_bodies*7:n_bodies*8]
                lc.dfdu[i,:] .+= trans.dgdu' * ia.I_of_f[end-trans.n-1:end-1]
                lc.dfdr[i] += ia.I_of_f[end] * inv_rstar
            end
        else
            # Integrate over exposure
            for j in eachindex(@view(tlim[1:end-1]))
                integrate_timestep!(t0, tlim[j], tlim[j+1], xc, yc, trans, ia)
                lc.flux[i] += ia.I_of_f[1]
            end
        end
    end
    return
end

function integrate_transit!(ib::Int64,it::Int64,t0::T,tc::SVector{N,T},trans::Transit_Struct{T},lc::Lightcurve{T},ts::TransitSeries{T},ia::T) where {N, T<:Real}
    inv_rstar = inv(lc.rstar[1])

    # Compute the flux at each point
    for i in eachindex(lc.tobs)
        # Check if observation is outside of a transit
        if lc.tobs[i] > tc[end]; break; end
        if lc.tobs[i] < tc[1]; continue; end

        # Get series expansion components
        xc = components(@views(ts.points[ib,it,:,1].*inv_rstar), ts.h)
        yc = components(@views(ts.points[ib,it,:,2].*inv_rstar), ts.h)

        if lc.do_grad
            n_bodies = length(ts.count)
            dxc = [components(ts.dpoints[ib,it,:,1,k,i].*inv_rstar, ts.h) for i in 1:7, k in 1:n_bodies][:]
            dyc = [components(ts.dpoints[ib,it,:,2,k,i].*inv_rstar, ts.h) for i in 1:7, k in 1:n_bodies][:]

            compute_flux!(lc.tobs[i], t0, xc, yc, dxc, dyc, lc, trans, i, ib-1, inv_rstar)
        else
            lc.flux[i] += compute_flux(lc.tobs[i], t0, xc, yc, trans)
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

function compute_flux(tc::T, t0::T, xc::SVector{N,T}, yc::SVector{N,T}, trans::Transit_Struct{T}) where {N, T<:Real}
    # Compute flux at a particular time
    trans.b = compute_impact_parameter(tc, t0, xc, yc)
    flux = transit_poly_g(trans) - 1
    return flux
end

function integrate_timestep!(t0::T, a::T, b::T, xc::SVector{N,T}, yc::SVector{N,T}, dxc, dyc, trans::Transit_Struct{T}, ia, dbdq0, ki) where {N,T<:Real}
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

function compute_flux!(tc::T, t0::T, xc, yc, dxc, dyc, lc, trans, i, ki, inv_rstar) where T<:Real
    # Compute flux and derivatives
    trans.b = compute_impact_parameter!(tc, t0, xc, yc, dxc, dyc, lc.dbdq0)
    lc.flux[i] += transit_poly_g!(trans) - 1            # Flux
    lc.dfdq0[i,:] .+= trans.dfdrb[2] .* lc.dbdq0        # Derivative wrt initial conditions
    lc.dfdk[i,ki] += trans.dfdrb[1]                     # Radius ratio
    lc.dfdu[i,:] .+= trans.dgdu' * trans.dfdg           # Limbdarkening params
    lc.dfdr[i] += -trans.dfdrb[2]*trans.b*inv_rstar     # Stellar radius
    return
end

# There should be a way to do this using sparse matrices
"""
Transform the jacobian wrt the initial Cartesian coordinates to initial orbital
elements.
"""
transform_to_elements!(s::State{T}, lc::Lightcurve{T}) where T<:Real = mul!(lc.dfdelements, lc.dfdq0, s.jac_init);