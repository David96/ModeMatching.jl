using SpecialFunctions, FunctionZeros, ForwardDiff, Roots

#
# Helpers for bessel functions
#
function besselj_prime(n, x)
    ForwardDiff.derivative(x1 -> besselj(n, x1), x)
end

max_n = 20
max_root = 100
bessel_prime_zeros = nothing
function precompute_zeros()
    println("Precomputing roots of J'...")
    global max_n, max_root, bessel_prime_zeros
    bessel_prime_zeros = Vector{Vector{AbstractFloat}}(undef, max_n)
    for n1 = 0:max_n-1
        bessel_prime_zeros[n1+1] = find_zeros(x -> ForwardDiff.derivative(
                                            x1 -> besselj(n1, x1), x), (1, max_root))
    end
    println("Done.")
end
function besselj_prime_zero(n, m)
    global max_n, max_root, bessel_prime_zeros
    if bessel_prime_zeros === nothing
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
    bessel_prime_zeros[n+1][m]
end

struct CylindricalWaveguide <: Waveguide
    r::AbstractFloat
    f::AbstractFloat
    ω::AbstractFloat
    λ::AbstractFloat
    ε::AbstractFloat
    μ::AbstractFloat
    k::AbstractFloat
    x::AbstractFloat
    y::AbstractFloat
    z::AbstractFloat
    length::AbstractFloat

    function CylindricalWaveguide(x, y, z, r, length, μ, ε, f)
        new(r, f, 2π * f, c / f, ε * ε0, μ * μ0, 2π * f * sqrt(μ * ε) / c, x, y, z, length)
    end
end

function mode_from_nr(_::CylindricalWaveguide, nr::Integer,
        n_TE::Integer, n_TM::Integer, max_m::Integer)
    N = nr > n_TE ? (nr - n_TE - 1) : (nr - 1)
    n = div(N, max_m) + 1
    m = N % max_m + 1
    nr > n_TE ? TMMode(m, n) : TEMode(m, n)
end

function intersect(g1::CylindricalWaveguide, g2::CylindricalWaveguide)
    lb = [0, 0]
    hb = [minimum([g1.r, g2.r]), 2π]
    (lb, hb, (_, _) -> 1)
end

function contains(g::CylindricalWaveguide, ρ, φ, z)
    return ρ <= g.r
end

jacobi_det(_::CylindricalWaveguide, ρ, φ, z) = ρ

k_c(g::CylindricalWaveguide, mode::TEMode) = besselj_prime_zero(mode.n, mode.m) / g.r
k_c(g::CylindricalWaveguide, mode::TMMode) = besselj_zero(mode.n, mode.m) / g.r

#
# TE Modes
#
E_spatial(g::CylindricalWaveguide, ρ, φ, _, mode::TEMode) = [
    cos(mode.n * φ) / ρ * besselj(mode.n, k_c(g, mode) * ρ),
    sin(mode.n * φ) * besselj_prime(mode.n, k_c(g, mode) * ρ),
    0
]
E_freq(g::CylindricalWaveguide, mode::TEMode) = [
    -j * g.ω * g.μ * mode.n / k_c(g, mode)^2,
    j * g.ω * g.μ / k_c(g, mode),
    0
]
H_spatial(g::CylindricalWaveguide, ρ, φ, _, mode::TEMode) = [
    sin(mode.n * φ) * besselj_prime(mode.n, k_c(g, mode) * ρ),
    cos(mode.n * φ) / ρ * besselj(mode.n, k_c(g, mode) * ρ),
    sin(mode.n * φ) * besselj(mode.n, k_c(g, mode) * ρ)
]
H_freq(g::CylindricalWaveguide, mode::TEMode) = [
    -j * β(g, mode) / k_c(g, mode),
    -j * β(g, mode) / k_c(g, mode)^2,
    1
]

#
# TM Modes
#
E_spatial(g::CylindricalWaveguide, ρ, φ, _, mode::TMMode) = [
    sin(mode.n * φ) * besselj_prime(mode.n, k_c(g, mode) * ρ),
    cos(mode.n * φ) / ρ * besselj(mode.n, k_c(g, mode) * ρ),
    sin(mode.n * φ) * besselj(mode.n, k_c(g, mode) * ρ)
]
E_freq(g::CylindricalWaveguide, mode::TMMode) = [
    -j * β(g, mode) / k_c(g, mode),
    -j * β(g, mode) / k_c(g, mode)^2,
    1
]
H_spatial(g::CylindricalWaveguide, ρ, φ, _, mode::TMMode) = [
    cos(mode.n * φ) / ρ * besselj(mode.n, k_c(g, mode) * ρ),
    sin(mode.n * φ) * besselj_prime(mode.n, k_c(g, mode) * ρ),
    0
]
H_freq(g::CylindricalWaveguide, mode::TMMode) = [
    -j * g.ω * g.μ * mode.n / k_c(g, mode)^2,
    j * g.ω * g.μ / k_c(g, mode),
    0
]
