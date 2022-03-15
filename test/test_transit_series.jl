function test_integrator_methods(n)
    # Setup and run the simulation
    BJD = 7250.0; t0_ic = 7.0
    tmax = 100.0
    ic = setup_ICs(n,BJD,t0_ic)
    intr = setup_integrator(ic, tmax)
    tt = compute_transit_times(ic, intr)
    pd_provided = compute_pd(ic, tt, intr)
    pd_computed = compute_pd(ic, intr);

    @test pd_provided.times == pd_computed.times
    @test pd_provided.bodies == pd_computed.bodies

    # Computed points array is likely longer than provided
    t_range = length(pd_provided.times)
    @test pd_provided.points ≈ pd_computed.points[:,1:t_range,:,:]

    # Make sure there's no extra computed points
    @test all(pd_computed.points[:,t_range+1:end,:,:] .== 0.0)
end

@testset "Integrator Methods" begin
    test_integrator_methods(8)
end