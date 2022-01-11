## For now, just make sure the code runs ##

@testset "Lightcurve" begin
    n = 8
    BJD = 7250.0; t0_ic = 7.0
    tmax = 2.0
    ic = setup_ICs(n, BJD, t0_ic);
    intr = setup_integrator(ic, tmax);
    tt = compute_transit_times(ic, intr);
    pd = compute_pd(ic, tt, intr);

    # Transit parameters
    k = get_radius_ratios_trappist(n);
    u_n = get_limdark_coeffs_trappist();
    rstar = 0.00465047 * 0.1192 # Trappist-1 (Rstar/AU)

    # Simulate lightcurve
    dt = 2 / 60 / 24; # 2 minutes in days
    obs_duration = maximum(tt.tt) + 5 * intr.h
    nobs = round(Int64, obs_duration * 30 * 24)
    tobs = collect(pd.times[1] - 10 * intr.h:obs_duration / nobs:obs_duration);
    lc = Lightcurve(dt, copy(tobs), copy(tobs), zeros(length(tobs)), u_n, k, rstar);
    lc_noint = Lightcurve(0.0, copy(tobs), copy(tobs), zeros(length(tobs)), u_n, k, rstar);
    dlc = Lightcurve(dt, copy(tobs), copy(tobs), zeros(length(tobs)), u_n, k, rstar, n*7);
    compute_lightcurve!(lc, pd, tol=1e-11, maxdepth=40)
    compute_lightcurve!(lc_noint, pd)

    @time compute_lightcurve!(lc, pd, tol=1e-11, maxdepth=40)
    @time compute_lightcurve!(lc_noint, pd)
    @time compute_lightcurve!(dlc, pd, tol=1e-11, maxdepth=40)
    @test 1 == 1
end