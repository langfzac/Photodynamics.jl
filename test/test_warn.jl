# Test the cases where we expect the code to break, warn, or otherwise produce
# erroneous results.

"""
Test the case where we compute more transits than array elements in the
preallocated arrays.
"""
function test_to_many_transits(n)
    ic = setup_ICs(n, 0.0, 0.0)
    tmax = 10.0
    intr = Integrator(0.05, tmax)

    s = State(ic)
    tt = TransitTiming(tmax, ic)
    ts = TransitSeries(tmax - 2.0, ic)
    @test_logs (:warn, "Exceeded times array: ts.times") intr(s, ts, tt)

    s = State(ic)
    tt = TransitTiming(tmax-7.0, ic)
    ts = TransitSeries(tmax, ic)
    @test_logs (:warn, "Exceeded transit times array: tt.tt") intr(s, ts, tt)
end

@testset "Warnings" begin
    test_to_many_transits(3)
end