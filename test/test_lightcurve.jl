## Test that Lightcurve constructors work
function test_lightcurve_constructors(n)
    dt = 2 / 60 / 24
    obs_duration = 100.0
    tobs = collect(0.0:dt:obs_duration)
    k = get_radius_ratios_trappist(n)
    u_n = get_limbdark_coeffs_trappist()
    rstar = get_trappist_rstar()

    lc1 = Lightcurve(dt, tobs, zeros(size(tobs)), zeros(size(tobs)), u_n, k, rstar)
    lc2 = Lightcurve(dt, obs_duration, u_n, k, rstar)

    for field in fieldnames(Lightcurve)
        @test getfield(lc1, field) == getfield(lc2, field)
    end
    return
end

@testset "Lightcurve Constructors" begin
    test_lightcurve_constructors(8)
end