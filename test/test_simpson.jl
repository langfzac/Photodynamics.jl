## Simple test of simpson integration ##
function test_integration()
    function f(x,s)
        s[1] = x^2 + 2*x
    end
    F(x) = x^3/3 + x^2 # Anti-derivative of f

    a = -1.0; b= 1.0 # Integration bounds
    ia = IntegralArrays(1,40,1e-8)
    integrate_simpson!(a,b,f,ia)
    @test ia.I_of_f[1] .≈ F(b) - F(a)
end

function test_gradient_integration()
    function f(x,s)
        c1 = 1.0
        c2 = 2.0

        s[1] = c1*x^2 + c2*x
        s[2] = x^2  # dfdc1
        s[3] = x    # dfdc2
        return
    end

    function F(x)
        c1 = 1.0
        c2 = 2.0

        s = c1/3*x^3 + c2/2*x^2
        dsdc1 = x^3/3
        dsdc2 = x^2/2
        return [s, dsdc1, dsdc2]
    end

    a = -1.0; b = 1.0
    ia = IntegralArrays(2, 40, 1e-8)
    integrate_simpson!(a,b,f,ia)
    @test ia.I_of_f ≈ F(b) - F(a)
end

@testset "Simpsons Integration" begin
    test_integration()
    test_gradient_integration()
end