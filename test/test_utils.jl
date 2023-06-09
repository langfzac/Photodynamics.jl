function test_linear_regression()
    # Test we can fit a straight line
    x = collect(1:10)
    m = 3.0; b = 2.0
    y = @. m*x+b
    @test all(linear_regression(x, y) .≈ [b, m])

    # Make sure the asserts work
    x_long = collect(1:100)
    @test_throws get_test_assertion("x and y must be the same length") linear_regression(x_long, y)
end

function test_get_transit_times()
    ic = setup_ICs(8, 0.0, 0.0)
    intr = setup_integrator(ic, 100.0)
    ts = compute_pd(ic, intr)

    expr = get_test_assertion("Must choose a valid body number: $(collect(2:8))")
    @test_throws expr get_transit_times(ts, 1)
    @test_throws expr get_transit_times(ts, 10)

    @test_throws expr get_linear_transit_times(ts, 1)
    @test_throws expr get_linear_transit_times(ts, 1)
end

@testset "Utility Functions" begin
    test_linear_regression()
    test_get_transit_times()
end