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

    # Get max transit number
    ymax = length(ts.times[ts.bodies .== ib])

    # Now setup heatmap
    planet_labels = "abcdefghijklmnopqrstuvwxyz" # Assumes star is "a" and indices are linear
    title --> "$(planet_labels[ib])"
    ylabel --> "Transit Number"
    yticks --> 0:1:ymax
    xlabel --> "Time [$(t_units)]"
    color --> :winter
    seriestype := :heatmap

    t_plt, t_numbers, flux
end
