using DelimitedFiles, FiniteDifferences

import Photodynamics: components, transit_init
import Photodynamics: compute_impact_parameter, compute_impact_parameter!
import Photodynamics: IntegralArrays, integrate_simpson!
import Photodynamics: integrate_timestep!, compute_flux, compute_flux!
import Photodynamics: NbodyGradient.amatrix

function isapprox_maxabs(a,b; kwargs...)
    isap = isapprox(a, b; norm=x->maximum(abs.(x)), kwargs...)
    if ~isap 
        @warn maximum(@.(abs(a-b)))
    end
    return isap
end

function get_test_assertion(s)
    # Direct error message testing only in >1.8
    # Otherwise just make sure the assertion works
    if VERSION >= v"1.8"
        return s
    else
        return AssertionError
    end
end

function setup_ICs(n, BJD::T, t0::T) where T<:Real
    elements = T.(readdlm("elements.txt", ',')[1:n,:])
    elements[2:end,3] .-= BJD # Shift initial transit times
    ic = ElementsIC(t0, n, elements)
    return ic
end

function setup_integrator(ic, tmax::T) where T<:Real
    h = ic.elements[2,2] / 40
    intr = Integrator(h, T(0.0), tmax)
    return intr
end

function compute_transit_times(s, ic, intr; grad=false)
    tt = TransitTiming(intr.tmax, ic)
    intr(s, tt; grad=grad)
    return tt
end

function compute_transit_times(ic, intr; grad=false)
    s = State(ic)
    tt = TransitTiming(intr.tmax, ic)
    intr(s, tt; grad=grad)
    return tt
end

function compute_pd(ic, intr; grad=false)
    s = State(ic);
    tt = TransitTiming(intr.tmax, ic);
    pd = TransitSeries(intr.tmax, ic);
    intr(s, pd, tt;grad=grad)
    return pd
end

function compute_pd(ic, tt, intr; grad=false)
    s = State(ic);
    pd = TransitSeries(tt.tt, ic);
    intr(s, pd;grad=grad)
    return pd
end

function compute_pd(s, ic, tt, intr; grad=false)
    pd = TransitSeries(copy(tt.tt), ic);
    intr(s, pd; grad=grad)
    return pd
end

function get_radius_ratios_trappist(n)
    depth = [0.7277,0.6940,0.3566,0.4802,0.634,0.764,0.346] .* 0.01
    return sqrt.(depth)[1:n-1]
end

function get_limbdark_coeffs_trappist()
    q = [0.11235270319764341, 0.42037661035916857]#, 0.352424321959808, 0.2864053200404355]
    return [2*sqrt(q[1])*q[2],2*sqrt(q[1])*(1-2q[2])]
end

get_trappist_rstar() = 0.00465047 * 0.1192 # Trappist-1 (Rstar/AU)

normalize_points!(points, rstar) = points./=rstar
