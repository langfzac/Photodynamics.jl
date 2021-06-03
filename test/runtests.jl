using Photodynamics
using Test

include("common.jl")

@testset "Photodynamics.jl" begin
    print("Impact Parameter... ")
    @testset "Impact Parameter" begin
        include("test_impact_parameter.jl")
    end
    println("Done.")

    print("Lightcurve... ")
    @testset "Lightcurve" begin
        include("test_simpson.jl")
        include("test_integrate_timestep.jl")
        include("test_compute_lightcurve.jl")
    end
end