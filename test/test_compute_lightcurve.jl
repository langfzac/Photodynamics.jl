## For now, just make sure the code runs ##

@testset "Lightcurve" begin
    n = 8
    BJD = 7250.0; t0_ic = 7.0
    tmax = 2.0
    ic = setup_ICs(n, BJD, t0_ic);
    intr = setup_integrator(ic, tmax);
    tt = compute_transit_times(ic, intr);
    pd = compute_pd(ic, tt, intr);

    # Normalize points to stellar radius
    rstar = 0.00465047 * 0.1192 # Trappist-1 (Rstar/AU)
    normalize_points!(pd.points, rstar);

    k = get_radius_ratios_trappist(n);
    u_n = get_limdark_coeffs_trappist();

    # Simulate lightcurve
    dt = 2 / 60 / 24; # 2 minutes in days
    obs_duration = maximum(tt.tt) + 5 * intr.h
    nobs = round(Int64, obs_duration * 30 * 24)
    tobs = collect(pd.times[1] - 10 * intr.h:obs_duration / nobs:obs_duration);
    lc = Lightcurve(dt, copy(tobs), copy(tobs), zeros(length(tobs)), u_n, k);
    compute_lightcurve!(lc, pd, tol=1e-11, maxdepth=40)

    @time compute_lightcurve!(lc, pd, tol=1e-11, maxdepth=40)
    @test 1 == 1
end