using Photodynamics
using Test, DelimitedFiles

########################## Common test functions ###############################
function setup_ICs(n,BJD,t0)
    elements = readdlm("elements.txt", ',')[1:n,:]
    elements[2:end,3] .-= BJD # Shift initial transit times
    ic = ElementsIC(t0, n, copy(elements))
    return ic
end

function setup_integrator(ic)
    h = ic.elements[2,2] / 40
    tmax = 2.0
    intr = Integrator(h, 0.0, tmax)
    return intr
end

function compute_transit_times(ic, intr)
    s = State(ic)
    tt = TransitTiming(intr.tmax, ic)
    intr(s, tt; grad=false)
    return tt
end

function compute_pd(ic, tt, intr)
    s = State(ic);
    pd = TransitSeries(tt.tt, ic);
    intr(s, pd;grad=false)
    return pd
end

function get_radius_ratios_trappist(n)
    depth = [0.7277,0.6940,0.3566,0.4802,0.634,0.764,0.346] .* 0.01
    return sqrt.(depth)[1:n-1]
end

normalize_points!(points, rstar) = points./=rstar
################################################################################

@testset "Photodynamics.jl" begin
    print("Impact Parameter...")
    @testset "Impact Parameter" begin
        include("test_impact_parameter.jl")
    end
    println("Done.")
end
