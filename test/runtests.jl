using Photodynamics
using Test

include("common.jl")

@testset "Photodynamics.jl" begin
    print("Integration Methods... ")
    @testset "Integration Methods" begin
        include("test_transit_series.jl")
    end
    println("Done.")

    print("Impact Parameter... ")
    @testset "Impact Parameter" begin
        include("test_impact_parameter.jl")
    end
    println("Done.")

    print("Lightcurve... ")
    @testset "Lightcurve" begin
        include("test_points_of_contact.jl")
        include("test_simpson.jl")
        include("test_integrate_timestep.jl")
        include("test_compute_flux.jl")
        include("test_lightcurve.jl")
        include("test_compute_lightcurve.jl")
    end

    print("Breaking Cases...")
    @testset "Breaking Cases" begin
        include("test_warn.jl")
    end
end