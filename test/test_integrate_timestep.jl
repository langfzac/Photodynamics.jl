# Setup initial conditions and compute points
n = 3
BJD = 7250.0; t0_ic = 7.0
tmax = 2.0
ic = setup_ICs(n, BJD, t0_ic);
intr = setup_integrator(ic, tmax);
tt = compute_transit_times(ic, intr);
pd = compute_pd(ic, tt, intr, grad=true);

# Normalize points to stellar radius
rstar = 0.00465047 * 0.1192 # Trappist-1 (Rstar/AU)
normalize_points!(pd.points, rstar);

# Make a fake exposure time
dt = 2/60/24 # 2 minute cadence in days
it = 1; ib = pd.bodies[it]
t0 = pd.times[it]
texp = t0 + dt/2

# Compute x and y components for impact parameter calculation
xc = components(pd.points[ib,it,:,1], pd.h)
yc = components(pd.points[ib,it,:,2], pd.h)

# Get a transit structure
k = get_radius_ratios_trappist(n);
u_n = get_limdark_coeffs_trappist();
trans = transit_init(k[ib-1], 0.0, u_n, false)

# Integrate a timestep over a -> b
a = texp - dt/2
b = texp + dt/2
tol = 1e-15
maxdepth = 40
ia = IntegralArrays(1, maxdepth, tol)
integrate_timestep!(t0, a, b, xc, yc, trans, ia)
flux_int = ia.I_of_f[1]

# Make bigfloat variants
ic_big = setup_ICs(n, big(BJD), big(t0_ic));
intr_big = setup_integrator(ic_big, big(tmax));
tt_big = compute_transit_times(ic_big, intr_big);
pd_big = compute_pd(ic_big, tt_big, intr_big, grad=false);
normalize_points!(pd_big.points, big(rstar));
xc_big = components(pd_big.points[ib,it,:,1], big(pd.h))
yc_big = components(pd_big.points[ib,it,:,2], big(pd.h))
ia_big = IntegralArrays(1, maxdepth, big(tol))

function compute_timestep(params, t0, a, b, xc, yc, ia, ib)
    trans = transit_init(params[ib-1], zero(typeof(t0)), params[end-1:end], false)
    integrate_timestep!(t0, a, b, xc, yc, trans, ia)
    return ia.I_of_f[1]
end

# Make sure we get the same answer
params = [k...,u_n...]
flux_big = compute_timestep(big.(params), big(t0), big(a), big(b), xc_big, yc_big, ia_big, ib)
@test flux_big ≈ flux_int

# Compute derivatives with respect to the transit parameters
grad_num = grad(central_fdm(5, 1), p -> compute_timestep(p, big(t0), big(a), big(b), xc_big, yc_big, ia_big, ib), big.(params))[1]

# Compute analytic derivatives
dxc = [components(pd.dpoints[ib,it,:,1,k,i], pd.h) for i in 1:7, k in 1:ic.nbody][:]
dyc = [components(pd.dpoints[ib,it,:,2,k,i], pd.h) for i in 1:7, k in 1:ic.nbody][:]
dbdq0 = zeros(7*ic.nbody)
ki = ib-1
trans_grad = transit_init(k[ki], 0.0, u_n, true)
n_params = length(dbdq0) + length(k) + length(u_n) + 1 # flux
ia_grad = IntegralArrays(n_params, maxdepth, tol)
integrate_timestep!(t0, a, b, xc, yc, dxc, dyc, trans_grad, ia_grad, dbdq0, ki)

# Check flux again
@test ia_grad.I_of_f[1] ≈ flux_int

# Check derivative of radius ratios
@test isapprox(ia_grad.I_of_f[2+length(dbdq0):end-3], Float64.(grad_num[1:2]), norm=x->maximum(abs.(x)))

# Now the derivatives of limbdark coefficients
dfdu = trans_grad.dgdu'*ia_grad.I_of_f[end-trans_grad.n:end]
@test isapprox(dfdu, Float64.(grad_num[3:end]), norm=x->maximum(abs.(x)))

# Finally, derivatives wrt the Nbody initial conditions
function compute_timestep_nbody(params, trans, ia, ic, tmax, a, b, it, ib)
    n = ic.nbody
    ic.elements .= reshape(params, n, 7)
    ic.m .= ic.elements[:,1]
    amatrix(ic)
    intr = setup_integrator(ic, tmax);
    tt = compute_transit_times(ic, intr);
    pd = compute_pd(ic, tt, intr, grad=false);
    rstar = big(0.00465047 * 0.1192) # Trappist-1 (Rstar/AU)
    normalize_points!(pd.points, rstar);

    t0 = pd.times[it]
    xc = components(pd.points[ib,it,:,1], pd.h)
    yc = components(pd.points[ib,it,:,2], pd.h)

    integrate_timestep!(t0, a, b, xc, yc, trans, ia)
    return ia.I_of_f[1]
end

params_nbody = copy(ic.elements[:])

# Again, check flux
flux_nbody = compute_timestep_nbody(params_nbody, trans, ia, deepcopy(ic), tmax, a, b, it, ib)
@test flux_nbody ≈ flux_int

# Compute finite FiniteDifferences
# trans_big = transit_init(big(k[ki]), big(0.0), big.(u_n), false)
# This gets killed due to memory issues... Tries to allocate >10GB
# grad_nbody_num = grad(central_fdm(5,1), p->compute_timestep_nbody(p, trans_big, ia_big, ic_big, big(tmax), big(a), big(b), it, ib), big.(params_nbody))[1]

# This also causes memory issues...
#=dq = 1e-10
params_big = big.(params_nbody)
params_copy = copy(params_big)
for i in eachindex(params_big)
    params_big .= params_copy
    params_big[i] += dq
    fluxp = compute_timestep_nbody(params_big, trans_big, ia_big, ic_big, big(tmax), big(a), big(b), it, ib)
    params_big[i] -= 2dq
    fluxm = compute_timestep_nbody(params_big, trans_big, ia_big, ic_big, big(tmax), big(a), big(b), it, ib)
    grad_nbody_num = (fluxp - fluxm)/(2*dq)
    println(grad_nbody_num," ",ia_grad.I_of_f[2:1+length(dbdq0)][i])
    #@test isapprox(grad_nbody_num, ia_grad.I_of_f[2:1+length(dbdq0)][i], norm=x->maximum(abs.(x)))
end=#