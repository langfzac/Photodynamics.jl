@shorthands riverplot

# Conversions from days
const TIME_UNITS = Dict(
    "days" => 1,
    "hrs" => 24,
    "hours" => 24,
    "minutes" => 1440,
    "mins" => 1440,
    "seconds" => 86400,
    "sec" => 86400,
    "secs" => 86400
)

@recipe function f(
    ts::TransitSeries,
    ib::Int,
    texp::Real,
    u_n::Vector,
    k::Real,
    rstar::Real;
    t_units="days",
    neach = 50,
    tol = 1e-6,
    maxdepth = 6
)
    # Make a row of pixels, centered on mid point
    pixels = [-div(neach, 2):div(neach,2)-1;] .+ 0.5

    # Get linear transit times
    t_numbers, tt_linear = get_linear_transit_times(ts, ib)

    # Get exposure times around each transit and map to pixels
    t_plt = texp.*pixels
    tobs = vec(mapreduce(tt -> tt .+ t_plt, hcat, tt_linear))

    # Apply units for plotting
    t_plt = t_plt.*TIME_UNITS[t_units]

    # Compute photometry at each exposure
    ks = zeros(length(ts.count)-1)
    ks[ib-1] = k
    lc = Lightcurve(texp, tobs, similar(tobs), similar(tobs), u_n, ks, rstar)
    compute_lightcurve!(lc, ts; body_index=ib, tol=tol, maxdepth=maxdepth)

    # Make "image" from flux
    flux = permutedims(reshape(lc.flux, neach, ts.count[ib]))

    # Now setup heatmap
    planet_labels = "abcdefghijklmnopqrstuvwxyz" # Assumes star is "a" and indices are linear
    title --> "$(planet_labels[ib])"
    color --> :winter
    seriestype := :heatmap

    t_plt, t_numbers, flux
end

#=function test_riverplot()
    # Build a lightcurve model
    # First, compute the dynamics
    n = 8
    tmax = 1600.0 # Days
    ic = setup_ICs(n, 0.0, 7257.93115525, fname="extra/elements_noprior_students.txt")
    intr = setup_integrator(ic, tmax)
    ts = compute_pd(ic, intr; grad=true) # Compute dynamics

    # Now, setup and compute photometry the lightcurve
    texp = 2.0 # Exposure time in minutes
    texp_day = texp / 60 / 24 # Exposure time in days
    u_n = [0.161, 0.208]
    k = get_radius_ratios_trappist(n)

    # Star density in units of solar density
    rho_T1 = 52.5363620668614
    AU = 1.49598e13
    RSUN = 6.955e10
    fac = rho_T1^(1 // 3) * AU / RSUN

    for i in 2:n
        display(riverplot(ts, i, texp_day, u_n, k[i-1], 1 / fac, neach=100, size=(800, 1200)))
    end
end=#