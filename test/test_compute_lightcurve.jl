## For now, just make sure the code runs ##

function test_lightcurve(n)
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

function compute_photodynamics(θ::AbstractVector{T}, n, integrated::Bool) where T<:Real
    t0 = zero(T)
    ind = 7 * (n - 1)

    elements = zeros(T, n, 7)
    elements[1] = 1.0
    elements[2:end, :] .= reshape(@views(θ[1:ind]), n - 1, 7)
    ic = ElementsIC(t0, n, elements)
    s = State(ic)

    u_n = θ[ind+1:ind+2]
    k = θ[ind+3:ind+4]
    rstar = θ[end]

    cadence = T(2 / 60 / 24)  # 2 minute cadence in days
    obs_duration = T(2.0)  # Duration of observations in days
    obs_times = collect(t0:cadence:obs_duration)
    dt = integrated ? cadence : zero(T)
    lc = Lightcurve(dt, obs_times, ones(T, length(obs_times)), zeros(T, length(obs_times)), u_n, k, rstar)

    intr = Integrator(T(0.05), T(0.0), obs_duration)
    ts = TransitSeries(obs_duration, ic)
    tt = TransitTiming(obs_duration, ic)
    intr(s, ts, tt; grad=false)
    compute_lightcurve!(lc, ts, tol=T(1e-8), maxdepth=20)
    return lc.flux
end

function test_jacobians(n)
    t0 = 0.0
    ic = setup_ICs(n, 0.0, t0)
    rstar = get_trappist_rstar()
    k = get_radius_ratios_trappist(n)
    u_n = get_limbdark_coeffs_trappist()

    # Setup simulated lightcurve
    cadence = 2 / 60 / 24  # 2 minute cadence in days
    obs_duration = 2.0  # Duration of observations in days
    obs_times = collect(t0:cadence:obs_duration)
    lc_noint = Lightcurve(0.0, obs_times, ones(length(obs_times)), zeros(length(obs_times)), u_n, k, rstar, 7*n)
    #lc_int = Lightcurve(cadence, obs_times, ones(length(obs_times)), zeros(length(obs_times)), u_n, k, rstar, 7*n)
 
    # Compute the dynamical model
    intr = Integrator(0.05, obs_duration)
    s = State(ic)
    ts = TransitSeries(obs_duration, ic)
    tt = TransitTiming(obs_duration, ic)
    intr(s, ts, tt; grad=true)

    # Compute the photometry
    compute_lightcurve!(lc_noint, ts)
    #compute_lightcurve!(lc_int, ts, tol=1e-8, maxdepth=20)

    # Transform the jacobian to orbital elements
    transform_to_elements!(s, lc_noint)
    #transform_to_elements!(s, lc_int)

    θ_init = [vec(ic.elements[2:end, :])..., u_n..., k..., rstar]
    jac_noint_num = jacobian(central_fdm(3,1), x->compute_photodynamics(x,n,false), big.(θ_init))[1]
    #jac_int_num = jacobian(central_fdm(3,1), x->compute_photodynamics(x,n,true), big.(θ_init))[1]

    jac_lc_inds = [7,14,1,8,2,9,3,10,4,11,5,12,6,13] .+ 7
    jac_num_inds = [1,2,3,4,5,6,7,8,9,10,11,12,13,14];

    # Check the orbital elements
    for i in eachindex(jac_lc_inds)
        @test isapprox_maxabs(jac_noint_num[:, jac_num_inds[i]], lc_noint.dfdelements[1:end-1, jac_lc_inds[i]])
        #@test isapprox_maxabs(jac_int_num[:, jac_num_inds[i]], lc_int.dfdelements[1:end, jac_lc_inds[i]])
    end

    # Check transit params
    @test isapprox_maxabs(jac_noint_num[:, 15], lc_noint.dfdu[1:end-1,1])
    @test isapprox_maxabs(jac_noint_num[:, 16], lc_noint.dfdu[1:end-1,2])
    @test isapprox_maxabs(jac_noint_num[:, 17], lc_noint.dfdk[1:end-1,1])
    @test isapprox_maxabs(jac_noint_num[:, 18], lc_noint.dfdk[1:end-1,2])
    @test isapprox_maxabs(jac_noint_num[:, 19], lc_noint.dfdr[1:end-1])

    #=@test isapprox_maxabs(jac_int_num[:, 15], lc_int.dfdu[1:end-1,1])
    @test isapprox_maxabs(jac_int_num[:, 16], lc_int.dfdu[1:end-1,2])
    @test isapprox_maxabs(jac_int_num[:, 17], lc_int.dfdk[1:end-1,1])
    @test isapprox_maxabs(jac_int_num[:, 18], lc_int.dfdk[1:end-1,2])
    @test isapprox_maxabs(jac_int_num[:, 19], lc_int.dfdr[1:end-1])=#
end

@testset "Lightcurve" begin
    test_lightcurve(8)
    test_jacobians(3)
end