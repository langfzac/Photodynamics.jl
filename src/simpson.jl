# A slight re-write of simpson_vec from github.com/rodluger/limbdark.jl

struct IntegralArrays{T<:AbstractFloat}
    g::Matrix{T}
    A::Matrix{T}
    S::Array{T,3}
    f_of_x::Vector{T}
    I_of_x::Vector{T}

    nf::Int64
    maxdepth::Int64
    tol::T

    function IntegralArrays(nf::Int64, maxdepth::Int64, tol::T) where T<:AbstractFloat
        g = zeros(T, nf, 5)
        A = zeros(T, nf, 3)
        S = zeros(T, nf, 3, maxdepth)
        f_of_x = zeros(T, nf)
        I_of_x = zeros(T, nf)
        return new{T}(g, A, S, f_of_x, I_of_x, nf, maxdepth, tol)
    end
end

function integrate_simpson!(a::T, b::T, f, ia) where T<:AbstractFloat
    third = one(T) / 3
    i = zero(T)
    m = 0; n = 0

    # 'Unpack' IntegralArrays struct
    g = ia.g
    A = ia.A
    S = ia.S
    f_of_x = ia.f_of_x
    I_of_f = ia.I_of_x
    epsilon = ia.tol
    N = ia.maxdepth
    nf = ia.nf

    # Zero out integral array
    I_of_f .= zero(T)

    f(a,f_of_x)
    @inbounds for j = 1:nf
        g[j,1] = f_of_x[j]
    end
    f(0.5 * (a + b),f_of_x)
    @inbounds for j = 1:nf
        g[j,3] = f_of_x[j]
    end
    f(b,f_of_x)
    @inbounds for j = 1:nf
        g[j,5] = f_of_x[j]
    end
    bma = b - a
    @inbounds for j = 1:nf
        A[j,1] = 0.5 * bma * (g[j,1] + 4 * g[j,3] + g[j,5])
    end
    @label AA
    d = 2^n
    h = 0.25 * bma / d

    f(a + h * (4 * m + 1),f_of_x)
    @inbounds for j = 1:nf
        g[j,2] = f_of_x[j]
    end
    f(a + h * (4 * m + 3),f_of_x)
    @inbounds for j = 1:nf
        g[j,4] = f_of_x[j]
    end
    @inbounds for j = 1:nf
        A[j,2] = h * (g[j,1] + 4 * g[j,2] + g[j,3])
        A[j,3] = h * (g[j,3] + 4 * g[j,4] + g[j,5])
    end
    maxdiff = zero(T)
    @inbounds for j = 1:nf
        maxa231 = abs((A[j,2] + A[j,3]) - A[j,1])
        if maxa231 > maxdiff
            maxdiff = maxa231
        end
    end
    if maxdiff > 3 * epsilon
        m *= 2
        n += 1
        if n > N
            @goto CC
        end
        @inbounds for j = 1:nf
            A[j,1] = A[j,2]
            S[j,1,n] = A[j,3]
            S[j,2,n] = g[j,4]
            S[j,3,n] = g[j,5]
            g[j,5] = g[j,3]
            g[j,3] = g[j,2]
        end
        @goto AA
    else
        @inbounds for j = 1:nf
            I_of_f[j] += (A[j,2] + A[j,3]) * third
        end
        m += 1
        i = a + m * bma / d
        @label BB
        if m == 2 * div(m, 2)
            m = div(m, 2)
            n -= 1
            @goto BB
        end
        if (m != 1) || (n != 0)
            @inbounds for j = 1:nf
                A[j,1] = S[j,1,n]
                g[j,1] = g[j,5]
                g[j,3] = S[j,2,n]
                g[j,5] = S[j,3,n]
            end
            @goto AA
        end
    end
    @label CC
    return
end