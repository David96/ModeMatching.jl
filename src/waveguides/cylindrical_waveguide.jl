using SpecialFunctions, FunctionZeros, StaticArrays, .BesselFunctions

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
    return ρ <= g.r && g.z <= z < g.z + g.length
end

function integral_deps(g::CylindricalWaveguide, _)
    g.r
end

jacobi_det(_::CylindricalWaveguide, ρ, φ, z) = ρ

k_c(g::CylindricalWaveguide, mode::TEMode) = besselj_prime_zero(mode.m, mode.n) / g.r
k_c(g::CylindricalWaveguide, mode::TMMode) = besselj_zero(mode.m, mode.n) / g.r

const P = 1

#
# TE Modes
#
E_spatial(g::CylindricalWaveguide, ρ, φ, _, mode::TEMode, _::Direction) = @SVector [
    (P * cos(mode.m * φ) - (1 - P) * sin(mode.m * φ)) / ρ * besselj(mode.m, k_c(g, mode) * ρ),
    (P * sin(mode.m * φ) + (1 - P) * cos(mode.m * φ)) * besselj_prime(mode.m, k_c(g, mode) * ρ),
    0
]
E_freq(g::CylindricalWaveguide, mode::TEMode, _::Direction) = @SVector [
    -j * g.ω * g.μ * mode.m / k_c(g, mode)^2,
    j * g.ω * g.μ / k_c(g, mode),
    0
]
H_spatial(g::CylindricalWaveguide, ρ, φ, _, mode::TEMode, _::Direction) = @SVector [
    (P * sin(mode.m * φ) + (1 - P) * cos(mode.m * φ)) * besselj_prime(mode.m, k_c(g, mode) * ρ),
    (P * cos(mode.m * φ) - (1 - P) * sin(mode.m * φ)) / ρ * besselj(mode.m, k_c(g, mode) * ρ),
    (P * sin(mode.m * φ) + (1 - P) * cos(mode.m * φ)) * besselj(mode.m, k_c(g, mode) * ρ)
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
    (P * sin(mode.m * φ) + (1 - P) * cos(mode.m * φ)) * besselj_prime(mode.m, k_c(g, mode) * ρ),
    (P * cos(mode.m * φ) - (1 - P) * sin(mode.m * φ)) / ρ * besselj(mode.m, k_c(g, mode) * ρ),
    (P * sin(mode.m * φ) + (1 - P) * cos(mode.m * φ)) * besselj(mode.m, k_c(g, mode) * ρ)
]
E_freq(g::CylindricalWaveguide, mode::TMMode, dir::Direction) = @SVector [
    -j * β(g, mode, dir) / k_c(g, mode),
    -j * β(g, mode, dir) * mode.m / k_c(g, mode)^2,
    1
]
H_spatial(g::CylindricalWaveguide, ρ, φ, _, mode::TMMode, _::Direction) = @SVector [
    (P * cos(mode.m * φ) - (1 - P) * sin(mode.m * φ)) / ρ * besselj(mode.m, k_c(g, mode) * ρ),
    (P * sin(mode.m * φ) + (1 - P) * cos(mode.m * φ)) * besselj_prime(mode.m, k_c(g, mode) * ρ),
    0
]
H_freq(g::CylindricalWaveguide, mode::TMMode, _::Direction) = @SVector [
    j * g.ω * g.ε * mode.m / k_c(g, mode)^2,
    -j * g.ω * g.ε / k_c(g, mode),
    0
]
