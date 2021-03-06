## For now, just make sure the code runs ##

@testset "Lightcurve" begin
    n = 8
    BJD = 7250.0
    t0_ic = 7.0
    tmax = 15.0
    ic = setup_ICs(n, BJD, t0_ic)
    intr = setup_integrator(ic, tmax)
    tt = compute_transit_times(ic, intr)
    pd = compute_pd(ic, tt, intr)
    pd_provided = deepcopy(pd)
    pd_computed = compute_pd(ic, intr)

    # Transit parameters
    k = get_radius_ratios_trappist(n)
    u_n = get_limbdark_coeffs_trappist()
    rstar = 0.00465047 * 0.1192 # Trappist-1 (Rstar/AU)

    # Simulate lightcurve
    dt = 2 / 60 / 24 # 2 minutes in days
    obs_duration = maximum(tt.tt) + 5 * intr.h
    nobs = round(Int64, obs_duration / dt)
    tobs = collect(pd.times[1]-10*intr.h:dt:obs_duration)
    lc = Lightcurve(dt, copy(tobs), copy(tobs), zeros(length(tobs)), u_n, k, rstar)
    lc_noint = Lightcurve(0.0, copy(tobs), copy(tobs), zeros(length(tobs)), u_n, k, rstar)
    dlc = Lightcurve(dt, copy(tobs), copy(tobs), zeros(length(tobs)), u_n, k, rstar, n * 7)
    compute_lightcurve!(lc, pd_provided, tol=1e-11, maxdepth=40)
    compute_lightcurve!(lc, pd_computed, tol=1e-11, maxdepth=40)

    if ~haskey(ENV, "CI")
        # Output computation times as naive indicator of issues
        @time compute_lightcurve!(dlc, pd_provided, tol=1e-11, maxdepth=40)
        @time compute_lightcurve!(dlc, pd_computed, tol=1e-11, maxdepth=40)
        @time compute_lightcurve!(lc_noint, pd_provided)
        @time compute_lightcurve!(lc_noint, pd_computed)
        @time compute_lightcurve!(lc, pd_provided, tol=1e-11, maxdepth=40)
        @time compute_lightcurve!(lc, pd_computed, tol=1e-11, maxdepth=40)
    end

    # Make sure different methods give same flux
    lc_p = Lightcurve(dt, copy(tobs), copy(tobs), zeros(length(tobs)), u_n, k, rstar)
    lc_c = Lightcurve(dt, copy(tobs), copy(tobs), zeros(length(tobs)), u_n, k, rstar)
    compute_lightcurve!(lc_p, pd_provided, tol=1e-11, maxdepth=40)
    compute_lightcurve!(lc_c, pd_computed, tol=1e-11, maxdepth=40)
    haskey(ENV, "CI") ? println("\nProvided - Computed: ", maximum(abs.(lc_p.flux .- lc_c.flux))) : nothing
    @test isapprox_maxabs(lc_p.flux, lc_c.flux)

    dlc_p = Lightcurve(dt, copy(tobs), copy(tobs), zeros(length(tobs)), u_n, k, rstar, n * 7)
    dlc_c = Lightcurve(dt, copy(tobs), copy(tobs), zeros(length(tobs)), u_n, k, rstar, n * 7)
    compute_lightcurve!(dlc_p, pd_provided, tol=1e-11, maxdepth=40)
    compute_lightcurve!(dlc_c, pd_computed, tol=1e-11, maxdepth=40)
    haskey(ENV, "CI") ? println("dProvided - dComputed: ", maximum(abs.(dlc_p.flux .- dlc_c.flux))) : nothing
    @test isapprox_maxabs(dlc_p.flux, dlc_c.flux)

    haskey(ENV, "CI") ? println("dProvided - Provided: ", maximum(abs.(dlc_p.flux .- lc_p.flux))) : nothing
    @test isapprox_maxabs(dlc_p.flux, lc_p.flux)
    haskey(ENV, "CI") ? println("dComputed - Computed: ", maximum(abs.(dlc_c.flux .- lc_c.flux))) : nothing
    @test isapprox_maxabs(dlc_c.flux, lc_c.flux)

    lc_noint = Lightcurve(0.0, copy(tobs), copy(tobs), zeros(length(tobs)), u_n, k, rstar)
    dlc_noint = Lightcurve(0.0, copy(tobs), copy(tobs), zeros(length(tobs)), u_n, k, rstar, n * 7)
    compute_lightcurve!(lc_noint, pd_computed, tol=1e-11, maxdepth=40)
    compute_lightcurve!(dlc_noint, pd_computed, tol=1e-11, maxdepth=40)
    haskey(ENV, "CI") ? println("dnoint - noint: ", maximum(abs.(dlc_noint.flux .- lc_noint.flux))) : nothing
    @test isapprox_maxabs(dlc_noint.flux, lc_noint.flux)
end