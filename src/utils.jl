# Collection of utility functions

"""
Compute a linear regression using linear least-squares.
"""
function linear_regression(x,y)
    @assert length(y) == length(x) "x and y must be the same length"
    N = length(x)
    X = [ones(N) x]

    β = (X'X) \ X'y
end

"""
Get the transit times for a particular body.
"""
function get_transit_times(ts::TransitSeries{<:Real, ComputedTimes}, ib::Integer)
    n = length(ts.count)
    # Assume only one star and the star is the first index
    @assert 1 < ib <= n "Must choose a valid body number: $(collect(2:n))"
    return ts.times[ts.bodies .== ib]
end

"""
Compute linear ephemeris for a particular transiting body.
"""
function get_linear_transit_times(ts, ib)
    n = length(ts.count)
    # Assume only one star and the star is the first index
    @assert 1 < ib <= n "Must choose a valid body number: $(collect(2:n))"

    transit_times = get_transit_times(ts, ib)
    transit_numbers = [0:1:ts.count[ib]-1;]
    t0, period_estimate = linear_regression(transit_numbers, transit_times)
    return transit_numbers, t0 .+ period_estimate .* transit_numbers
end