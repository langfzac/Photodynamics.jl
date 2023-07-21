@testset "Plots" begin
    # Setup dynamics
    ic = get_default_ICs("trappist-1")
    intr = setup_integrator(ic, 20.0)
    ts = compute_pd(ic, intr; grad=false)

    # Test the plot recipe
    rec = apply_recipe(Dict{Symbol, Any}(), ts, 2, 0.01, [0.1,0.2], 0.01, 0.001)
    @test getfield(rec[1], 1) == Dict{Symbol, Any}(:color => :winter,
                                                   :ylabel => "Transit Number",
                                                   :title => "b",
                                                   :xlabel => "Time [days]",
                                                   :seriestype => :heatmap,
                                                   :yticks => 0:1:13)
end