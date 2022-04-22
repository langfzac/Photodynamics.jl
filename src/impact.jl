# Uses method from Parviainen et al. 2020 to compute the sky-separation.

"""Compute the derivative components of the series expansion."""
function components(d::AbstractVector{T},h::T) where T<:AbstractFloat
    A = d[7]-d[1]; B = d[1]+d[7]; C = d[2]-d[6]
    D = d[2]+d[6]; E = d[5]-d[3]; F = d[3]+d[5]
    v = dot(SVector(A, C, E), COEFF.vd)/h
    a = dot(SVector(B, D, F, d[4]), COEFF.ad)/(h*h)
    j = dot(SVector(-A, -C, -E), COEFF.jd)/(h*h*h)
    s = dot(SVector(B, D, F , d[4]), COEFF.sd)/(h*h*h*h)
    return SVector(d[4],v,0.5*a,(1/6)*j,(1/24)*s)
end

"""Compute the impact parameter at t=tc using a series expansion about t0."""
function compute_impact_parameter(tc::T, t0::T, xc::AbstractVector{T}, yc::AbstractVector{T}) where T<:Real
    t = tc - t0
    ts = SVector{5, T}(1.0, t, t*t, t*t*t, t*t*t*t)
    lx = dot(xc, ts)
    ly = dot(yc, ts)
    return sqrt(lx*lx + ly*ly)
end

function compute_impact_parameter(tc::T,t0::T,h::T,points::AbstractMatrix{T}) where T<:Real
    t = tc - t0
    ts = SVector{5, T}(1.0,t,t*t,t*t*t,t*t*t*t)
    xc = components(points[:,1], h)
    yc = components(points[:,2], h)
    lx = dot(xc,ts)
    ly = dot(yc,ts)
    return sqrt(lx*lx + ly*ly)
end

function compute_impact_parameter!(tc::T, t0::T, xc, yc, dxc, dyc, grad) where T<:Real
    t = tc - t0
    ts = SVector{5, T}(1.0, t, t*t, t*t*t, t*t*t*t)
    lx = dot(xc, ts)
    ly = dot(yc, ts)
    b = sqrt(lx*lx + ly*ly)

    dbdlx = lx/b; dbdly = ly/b;
    N = length(grad)
    for p in 1:N
        dlxdp = dot(dxc[p],ts)
        dlydp = dot(dyc[p],ts)
        grad[p] = dbdlx*dlxdp + dbdly*dlydp
    end

    return b
end