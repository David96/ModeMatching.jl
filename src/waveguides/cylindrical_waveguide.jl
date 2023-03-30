using SpecialFunctions, FunctionZeros, ForwardDiff, Roots, StaticArrays

#
# Helpers for bessel functions
#
function besselj_prime(n, x)
    ForwardDiff.derivative(x1 -> besselj(n, x1), x)
end

max_n::Int32 = 20
max_root::Int32 = 100
bessel_prime_zeros::Vector{Vector{Float64}} = []
function precompute_zeros()
    println("Precomputing roots of J'...")
    global max_n, max_root, bessel_prime_zeros
    bessel_prime_zeros = Vector{Vector{Float64}}(undef, max_n)
    for n1 = 0:max_n-1
        bessel_prime_zeros[n1+1] = find_zeros(x -> ForwardDiff.derivative(
                                            x1 -> besselj(n1, x1), x), (1, max_root))
    end
    println("Done.")
end
const bessel_lock = ReentrantLock()
function besselj_prime_zero(n, m)
    global max_n, max_root, bessel_prime_zeros
    lock(bessel_lock)
    if isempty(bessel_prime_zeros)
        precompute_zeros()
    end
    if n > max_n
        max_n *= 2
        precompute_zeros()
    end
    while m > length(bessel_prime_zeros[n+1])
        max_root *= 2
        precompute_zeros()
    end
    unlock(bessel_lock)
    bessel_prime_zeros[n+1][m]
end

struct CylindricalWaveguide{T<:AbstractFloat} <: Waveguide
    r::T
    f::T
    ω::T
    λ::T
    ε::T
    μ::T
    k::T
    x::T
    y::T
    z::T
    length::T
    A::T
    r_TE::T
    r_TM::T
    max_m::Int

    function CylindricalWaveguide(x, y, z, r::T, length::T,
            μ, ε, f, r_TE, r_TM, max_m) where T<:AbstractFloat
        new{T}(r, f, 2π * f, c / f, ε * ε0, μ * μ0, 2π * f * sqrt(μ * ε) / c,
               x, y, z, length, π * r^2, r_TE, r_TM, max_m)
    end
end

function mode_from_nr(g::CylindricalWaveguide, nr, n_total)
    n_TE = round(Int, n_total * g.r_TE)
    N = nr > n_TE ? (nr - n_TE - 1) : (nr - 1)
    n = div(N, g.max_m) + 1
    m = N % g.max_m + 1
    nr > n_TE ? TMMode(m, n) : TEMode(m, n)
end

function intersect(g1::CylindricalWaveguide, g2::CylindricalWaveguide)
    # TODO: off center intersections
    lb = [1e-6, 0]
    hb = [minimum([g1.r, g2.r]), 2π]
    (lb, hb, (_, _) -> 1)
end

function contains(g::CylindricalWaveguide, ρ, φ, z)
    return ρ <= g.r
end

jacobi_det(_::CylindricalWaveguide, ρ, φ, z) = ρ

k_c(g::CylindricalWaveguide, mode::TEMode) = besselj_prime_zero(mode.m, mode.n) / g.r
k_c(g::CylindricalWaveguide, mode::TMMode) = besselj_zero(mode.m, mode.n) / g.r

#
# TE Modes
#
E_spatial(g::CylindricalWaveguide, ρ, φ, _, mode::TEMode, _::Direction) = @SVector [
    cos(mode.m * φ) / ρ * besselj(mode.m, k_c(g, mode) * ρ),
    sin(mode.m * φ) * besselj_prime(mode.m, k_c(g, mode) * ρ),
    0
]
E_freq(g::CylindricalWaveguide, mode::TEMode, _::Direction) = @SVector [
    -j * g.ω * g.μ * mode.m / k_c(g, mode)^2,
    j * g.ω * g.μ / k_c(g, mode),
    0
]
H_spatial(g::CylindricalWaveguide, ρ, φ, _, mode::TEMode, _::Direction) = @SVector [
    sin(mode.m * φ) * besselj_prime(mode.m, k_c(g, mode) * ρ),
    cos(mode.m * φ) / ρ * besselj(mode.m, k_c(g, mode) * ρ),
    sin(mode.m * φ) * besselj(mode.m, k_c(g, mode) * ρ)
]
H_freq(g::CylindricalWaveguide, mode::TEMode, dir::Direction) = @SVector [
    -j * β(g, mode, dir) / k_c(g, mode),
    -j * β(g, mode, dir) * mode.m / k_c(g, mode)^2,
    1
]

#
# TM Modes
#
E_spatial(g::CylindricalWaveguide, ρ, φ, _, mode::TMMode, _::Direction) = @SVector [
    sin(mode.m * φ) * besselj_prime(mode.m, k_c(g, mode) * ρ),
    cos(mode.m * φ) / ρ * besselj(mode.m, k_c(g, mode) * ρ),
    sin(mode.m * φ) * besselj(mode.m, k_c(g, mode) * ρ)
]
E_freq(g::CylindricalWaveguide, mode::TMMode, dir::Direction) = @SVector [
    -j * β(g, mode, dir) / k_c(g, mode),
    -j * β(g, mode, dir) * mode.m / k_c(g, mode)^2,
    1
]
H_spatial(g::CylindricalWaveguide, ρ, φ, _, mode::TMMode, _::Direction) = @SVector [
    cos(mode.m * φ) / ρ * besselj(mode.m, k_c(g, mode) * ρ),
    sin(mode.m * φ) * besselj_prime(mode.m, k_c(g, mode) * ρ),
    0
]
H_freq(g::CylindricalWaveguide, mode::TMMode, _::Direction) = @SVector [
    j * g.ω * g.ε * mode.m / k_c(g, mode)^2,
    -j * g.ω * g.ε / k_c(g, mode),
    0
]
