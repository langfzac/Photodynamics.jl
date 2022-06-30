# Test the precision of the series expansion
function test_impact_parameter_accuracy(n)
    # Setup the simulation
    BJD = 7250.0; t0_ic = 7.0
    tmax = 15.0
    ic = setup_ICs(n,BJD,t0_ic)
    intr = setup_integrator(ic, tmax)
    tt = compute_transit_times(ic, intr)
    pd_provided = compute_pd(ic, tt, intr)
    pd_computed = compute_pd(ic, intr)

    # Normalize points to stellar radius
    rstar = 0.00465047 * 0.1192 # Trappist-1 (Rstar/AU)
    normalize_points!(pd_provided.points, rstar)
    normalize_points!(pd_computed.points, rstar)

    # Get radius ratios
    k = get_radius_ratios_trappist(n)

    # Test the provided times methods
    for it in 1:length(pd_provided.times)
        # Compute the transit duration
        ib = pd_provided.bodies[it]
        t0 = pd_provided.times[it];
        s = State(ic)
        els = NbodyGradient.get_orbital_elements(s, ic)
        b0 = compute_impact_parameter(t0, t0, 2e-2, pd_provided.points[ib,it,:,:])
        a = els[ib].a / rstar

        # Compute impact parameter over transit duration
        resolution = 0.0001
        bound = compute_transit_duration(b0, els[ib].P, k[ib-1], a, els[ib].I) / 2
        tc = collect(pd_provided.times[it] - bound:resolution:pd_provided.times[it] + bound)
        xc = components(@views(pd_provided.points[ib,it,:,1]), 2e-2)
        yc = components(@views(pd_provided.points[ib,it,:,2]), 2e-2)
        b_series = compute_impact_parameter.(tc, t0, Ref(xc), Ref(yc))

        # Now output the actual sky-separation from NbodyGradient
        positions = zeros(length(tc), 2)
        for (i, t) in enumerate(tc)
            s = State(ic)
            intr(s, t,;grad=false)
            positions[i,:] .= (s.x[1:2,ib] - s.x[1:2,1]) ./ rstar
        end
        b_actual = sqrt.(positions[:,1].^2 .+ positions[:,2].^2)

        @test isapprox(b_actual, b_series, norm=x -> maximum(abs.(x)))
    end
end

function test_impact_parameter_derivatives(n)
    # Test the derivatives of the impact parameter with respect
    # to the initial cartesian coordinates and masses

    # Setup the simulation
    BJD = 7250.0; t0_ic = 7.0
    tmax = 2.0
    ic = setup_ICs(n, BJD, t0_ic)
    intr = setup_integrator(ic, tmax)
    tt = compute_transit_times(ic, intr)
    pd = compute_pd(ic, tt, intr, grad=true)

    k = get_radius_ratios_trappist(n)
    rstar = get_trappist_rstar()
    normalize_points!(pd.points, rstar)
    normalize_points!(pd.dpoints, rstar)

    # Now run BigFloat precision simulations
    ic_copy = setup_ICs(n, big(BJD), big(t0_ic))
    intr_copy = setup_integrator(ic_copy, big(tmax))
    s_copy = State(ic_copy)
    tt_copy = compute_transit_times(ic_copy, intr_copy)
    pd_copy = TransitSeries(intr_copy.tmax, ic_copy)
    d_copy = NbodyGradient.Derivatives(BigFloat, ic_copy.nbody)

    # Compute the impact parameter
    # Used in finite difference derivatives
    function compute_b(coords, pd, tt, d, intr, ic, it, tc::Union{<:Real, Nothing}=nothing)
        T = eltype(coords)
        h::T = 2e-2
        s = State(ic)
        coords = reshape(coords, 7, n)
        for (i, col) in enumerate(eachcol(coords))
            s.x[:,i] .= col[1:3]
            s.v[:,i] .= col[4:6]
            s.m[i]    = col[end]
        end

        Photodynamics.zero_out!(pd)
        NbodyGradient.zero_out!(tt)
        intr(s, pd, tt, d; grad=false)
        #pd = compute_pd(s, ic, intr, grad=false)

        ib = pd.bodies[it]
        t0 = pd.times[it]
        tc isa Nothing ? tc = t0 : nothing
        points = pd.points[ib,it,:,:]
        xc = components(points[:,1], h)
        yc = components(points[:,2], h)
        return compute_impact_parameter(tc, t0, xc, yc)
    end

    # Compute the impact parameter and analytic gradient
    function compute_grad_b(coords, intr, ic, it, tc::Union{<:Real, Nothing}=nothing)
        h = 2e-2
        s = State(ic)
        coords = reshape(coords, 7, n)
        for (i, col) in enumerate(eachcol(coords))
            s.x[:,i] .= col[1:3]
            s.v[:,i] .= col[4:6]
            s.m[i]    = col[end]
        end

        pd = compute_pd(s, ic, intr, grad=true)

        ib = pd.bodies[it]
        t0 = pd.times[it]
        tc isa Nothing ? tc = t0 : nothing
        points = pd.points[ib,it,:,:]
        xc = components(points[:,1], h)
        yc = components(points[:,2], h)
        N = length(coords)
        grad = zeros(N)
        inds = [1,2,3,4,5,6,7]
        dxc = [components(pd.dpoints[ib,it,:,1,k,i], h) for i in inds, k in 1:ic.nbody][:]
        dyc = [components(pd.dpoints[ib,it,:,2,k,i], h) for i in inds, k in 1:ic.nbody][:]
        compute_impact_parameter!(tc, t0, xc, yc, dxc, dyc, grad)
        return grad
    end

    # Get orbital elements at the initial conditions
    s = State(ic)
    els = NbodyGradient.get_orbital_elements(s, ic)

    # Compare numerical and analytic derivative wrt the inital cartesian coordinates and masses
    coords = vcat([[x...,v...,m] for (x, v, m) in zip(eachcol(s_copy.x), eachcol(s_copy.v), s_copy.m)]...)
    for it in eachindex(pd.times)
        t0 = pd.times[it]
        ib = pd.bodies[it]
        b0 = compute_impact_parameter(t0, t0, 2e-2, pd.points[ib,it,:,:])
        bound = compute_transit_duration(b0, els[ib].P, k[ib-1], els[ib].a / rstar, els[ib].I) / 2
        times = range(pd.times[it]-bound, stop=pd.times[it]+bound, length=10)

        # Compute gradient at each time
        for tc in times
            grad_num = grad(central_fdm(2, 1), c -> compute_b(c, pd_copy, tt_copy, d_copy, intr_copy, ic_copy, it, big(tc)), big.(coords))[1]
            grad_analytic = compute_grad_b(Float64.(coords), intr, ic, it, tc)
            @test isapprox(asinh.(Float64.(grad_num)), asinh.(grad_analytic), norm=x -> maximum(abs.(x)))
        end
    end
end

@testset "Impact Parameter" begin
    test_impact_parameter_accuracy(8)
    test_impact_parameter_derivatives(3) # Takes a long time...
end