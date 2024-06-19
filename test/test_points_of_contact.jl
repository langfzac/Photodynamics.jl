function test_points_of_contact()
    n = 8
    ic = setup_ICs(n, 0.0, 0.0)
    els = NbodyGradient.get_orbital_elements(State(ic), ic)
    k = get_radius_ratios_trappist(n)
    rstar = get_trappist_rstar()

    # Estimate the sky-plane velocity
    # Assume linear trajectory with b0=0, inclination=pi/2
    ib = 8
    dt = 2e-2
    T_dur = compute_transit_duration(0.0, els[ib].P, k[ib-1], els[ib].a/rstar, pi/2)
    V = 2/T_dur * (1 + k[ib-1]) # [rstar/day]
    
    # Generate fake trajectory
    times = collect(-3*dt:dt:3*dt)
    points = zeros(7,2)
    points[:,1] .= V.*times

    # Expect 4 points of contact
    # Compare points of contact results to analytic
    tc_analytic_4 = [-T_dur/2, -(1 - k[ib-1])/V, (1 - k[ib-1])/V, T_dur/2]
    tc_computed_4 = points_of_contact_4(0.0,0.0,dt,points,k[ib-1])
    @test all(tc_analytic_4 .≈ tc_computed_4)

    # Now a "grazing" transit -- 2 points of contact
    # Just shift trajectory in y
    a = 1 - k[ib-1]/3; c = 1 + k[ib-1]; b = sqrt(c^2 - a^2) 
    points[:,2] .= a
    tc_analytic_2 = [-b/V, b/V]
    tc_computed_2 = points_of_contact_2(0.0,0.0,dt,points,k[ib-1])
    @test all(tc_analytic_2 .≈ tc_computed_2)
end

@testset "Points of Contact" begin
    test_points_of_contact()
end