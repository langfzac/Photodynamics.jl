module Photodynamics

using Reexport
@reexport using NbodyGradient

using StaticArrays, LinearAlgebra, Roots, SpecialFunctions, Limbdark

import Limbdark: transit_poly_g, transit_poly_g!, Transit_Struct
import NbodyGradient: InitialConditions, Derivatives, ElementsIC, TransitOutput
import NbodyGradient: check_step, set_state!

export dot, compute_lightcurve!
export TransitSeries, Lightcurve
export transform_to_elements!

include("simpson.jl")
include("impact.jl")
include("TransitSeries.jl")
include("Lightcurve.jl")

"""Constant coefficients for the series expansion."""
struct Coefficients{T<:Real}
    vd::StaticVector{3,T}
    ad::StaticVector{4,T}
    jd::StaticVector{3,T}
    sd::StaticVector{4,T}

    function Coefficients(::Type{T}=Float64) where T<:Real
        vd = SVector{3,T}(1.0/60, 9/60, 45/60)
        ad = SVector{4,T}(1.0/90, -3/20, 3/2, -49/18)
        jd = SVector{3,T}(1.0/8, 1, 13/8)
        sd = SVector{4,T}(-1.0/6, 2, -13/2, 28/3)
        return new{T}(vd,ad,jd,sd)
    end
end

const COEFF = Coefficients();

end
