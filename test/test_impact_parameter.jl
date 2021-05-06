# Test the precision of the series expansion
import Photodynamics:compute_impact_parameter

function run_tests(n)
    # Setup the simulation
    BJD = 7250.0; t0_ic = 7.0
    ic = setup_ICs(n,BJD,t0_ic)
    intr = setup_integrator(ic)
    tt = compute_transit_times(ic, intr)
    pd = compute_pd(ic, tt, intr)

    # Normalize points to stellar radius
    rstar = 0.00465047 * 0.1192 # Trappist-1 (Rstar/AU)
    normalize_points!(pd.points, rstar)

    # Get radius ratios
    k = get_radius_ratios_trappist(n)

    for it in 1:length(pd.times)
        # Compute the transit duration
        ib = pd.bodies[it]
        t0 = pd.times[it];
        s = State(ic)
        els = NbodyGradient.get_orbital_elements(s, ic)
        b0 = compute_impact_parameter(t0, t0, 2e-2, pd.points[ib,it,:,:])
        a = els[ib].a / rstar

        # Compute impact parameter over transit duration
        resolution = 0.0001
        bound = els[ib].P / pi * asin(sqrt((1 + k[ib - 1])^2 + b0^2) / (a * sin(els[ib].I))) / 2
        tc = collect(pd.times[it] - bound:resolution:pd.times[it] + bound)
        b_series = compute_impact_parameter.(tc, t0, 2e-2, Ref(pd.points[ib,it,:,:]))

        # Now output the actual sky-separation from NbodyGradient
        positions = zeros(length(tc), 2)
        for (i, t) in enumerate(tc)
            s = State(ic)
            intr(s, t,;grad=false)
            positions[i,:] .= s.x[1:2,ib] ./ rstar
        end
        b_actual = sqrt.(positions[:,1].^2 .+ positions[:,2].^2)

        @test isapprox(b_actual, b_series, norm=x -> maximum(abs.(x)))
    end
end

@testset "Impact Parameter" begin
    n = 8
    run_tests(8)
end