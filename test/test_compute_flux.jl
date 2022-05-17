# Used in finite difference derivatives
function compute_timestep(params, t0, points, h, ib, it)
    rstar = params[end]
    xc = components(points[ib,it,:,1]./rstar, h)
    yc = components(points[ib,it,:,2]./rstar, h)
    trans = transit_init(params[ib - 1], zero(typeof(t0)), params[end - 2:end-1], false)
    dt = 2e-2
    return compute_flux(t0+dt, t0, xc, yc, trans)
end

# Used in finite difference derivatives
function compute_timestep_nbody(coords, s_copy, tt, trans, ic, intr, it, ib)
    s = deepcopy(s_copy)
    coords = reshape(coords, 7, s.n)
    for (i, col) in enumerate(eachcol(coords))
        s.x[:,i] .= col[1:3]
        s.v[:,i] .= col[4:6]
        s.m[i]    = col[end]
    end

    pd = compute_pd(s, ic, tt, intr, grad=false);
    rstar = big(0.00465047 * 0.1192) # Trappist-1 (Rstar/AU)
    normalize_points!(pd.points, rstar);

    t0 = pd.times[it]
    xc = components(pd.points[ib,it,:,1], pd.h)
    yc = components(pd.points[ib,it,:,2], pd.h)

    dt = 2e-2
    return compute_flux(t0+dt, t0, xc, yc, trans)
end

function test_compute_timestep_derivatives(n)
    # Setup initial conditions and compute points
    BJD = 7250.0; t0_ic = 7.0
    tmax = 2.0
    k = get_radius_ratios_trappist(n);
    u_n = get_limbdark_coeffs_trappist();
    ic = setup_ICs(n, BJD, t0_ic);
    intr = setup_integrator(ic, tmax);
    tt = compute_transit_times(ic, intr);
    pd = compute_pd(ic, tt, intr, grad=true);

    # Normalize points to stellar radius
    rstar = 0.00465047 * 0.1192 # Trappist-1 (Rstar/AU)
    normalize_points!(pd.points, rstar);
    normalize_points!(pd.dpoints, rstar);

    for it in 1:length(pd.times)
        ib = pd.bodies[it]
        t0 = pd.times[it]

        # Get a transit structure
        trans = transit_init(k[ib - 1], 0.0, u_n, false)

        # Compute x and y components for impact parameter calculation
        xc = components(pd.points[ib,it,:,1], pd.h)
        yc = components(pd.points[ib,it,:,2], pd.h)

        # Compute the flux at t0+dt
        dt = 2e-2
        flux = compute_flux(t0+dt, t0, xc, yc, trans)

        # Make bigfloat variants
        ic_big = setup_ICs(n, big(BJD), big(t0_ic));
        intr_big = setup_integrator(ic_big, big(tmax));
        tt_big = compute_transit_times(ic_big, intr_big);
        pd_big = compute_pd(ic_big, tt_big, intr_big, grad=false);

        # Make sure we get the same answer
        params = [k...,u_n...,rstar]
        flux_big = compute_timestep(big.(params), big(t0), copy(pd_big.points), pd_big.h, ib, it)
        @test flux_big ≈ flux

        # Compute finite diff derivatives with respect to the transit parameters
        grad_num = grad(central_fdm(5, 1), p -> compute_timestep(p, big(t0), copy(pd_big.points), pd_big.h, ib, it), big.(params))[1]

        # Compute analytic derivatives
        inds = [1,2,3,4,5,6,7]
        dxc = [components(pd.dpoints[ib,it,:,1,k,i], pd.h) for i in inds, k in 1:ic.nbody][:]
        dyc = [components(pd.dpoints[ib,it,:,2,k,i], pd.h) for i in inds, k in 1:ic.nbody][:]
        dbdq0 = zeros(7 * ic.nbody)
        ki = ib - 1
        trans_grad = transit_init(k[ki], 0.0, u_n, true)
        n_params = length(dbdq0) + length(k) + length(u_n) + 1 # flux
        lc = Lightcurve(0.0, [0.0], [0.0], [0.0], u_n, k, rstar, 7*(length(k)+1))
        compute_flux!(t0+dt, t0, xc, yc, dxc, dyc, lc, trans_grad, 1, ki, inv(rstar))

        # Check flux again
        @test lc.flux[1] ≈ flux

        # Check derivative wrt radius ratios
        @test isapprox(lc.dfdk[1,:], Float64.(grad_num[1:n-1]), norm=x -> maximum(abs.(x)))

        # Now the derivatives wrt limbdark coefficients
        @test isapprox(lc.dfdu[1,:], Float64.(grad_num[n:end-1]), norm=x -> maximum(abs.(x)))

        # Now the derivatives wrt the stellar radius
        @test isapprox(lc.dfdr[1], Float64.(grad_num[end]), norm=x -> maximum(abs(x)), atol=1e-8)

        # Finally, derivatives wrt the Nbody initial conditions
        s_copy = deepcopy(State(ic_big))
        params_nbody = vcat([[x...,v...,m] for (x, v, m) in zip(eachcol(s_copy.x), eachcol(s_copy.v), s_copy.m)]...)

        # Again, check flux
        trans = transit_init(k[ib - 1], 0.0, u_n, false)
        flux_nbody = compute_timestep_nbody(params_nbody, State(ic), tt, trans, deepcopy(ic), intr, it, ib)
        @test flux_nbody ≈ flux

        # Compute finite FiniteDifferences
        trans_big = transit_init(big(k[ki]), big(0.0), big.(u_n), false)
        grad_nbody_num = grad(central_fdm(12, 1), p -> compute_timestep_nbody(p, s_copy, tt_big, trans_big, ic_big, intr_big, it, ib), big.(params_nbody))[1]

        # Compare with analytic
        @test isapprox(lc.dfdq0[1,:], Float64.(grad_nbody_num), norm=x -> maximum(abs.(x)))
    end
end

@testset "Non-Integrated Flux" begin
    test_compute_timestep_derivatives(3) # This runs for a while...
end