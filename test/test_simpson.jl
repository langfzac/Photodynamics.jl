import Photodynamics: IntegralArrays, integrate_simpson!

## Simple test of simpson integration ##
function f(x,s)
    s[1] = x^2 + 2*x
end
F(x) = x^3/3 + x^2 # Anti-derivative of f

a = -1.0; b= 1.0 # Integration bounds
ia = IntegralArrays(1,40,1e-8)
integrate_simpson!(a,b,f,ia)
@test ia.I_of_x[1] .≈ F(b) - F(a)