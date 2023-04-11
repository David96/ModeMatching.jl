using SpecialPolynomials

struct GaussianBeam{T} <: Waveguide where T <: AbstractFloat
    k::T
    f::T
    ω::T
    ε::T
    μ::T
    n::T
    w_0::T
    z_0::T
    z_R::T
    x::T
    y::T
    z::T
    length::T
    A::T
    max_m::Int

    function GaussianBeam(x, y, z, w_0, z_0, length, μ, ε, f, max_m)
        new{typeof(sqrt(μ*ε))}(2π * f * sqrt(μ * ε) / c, f, 2π * f, ε*ε0, μ*μ0, sqrt(ε), w_0,
            z_0, π*w_0^2*sqrt(ε)/(c/f), x, y, z, length, π, max_m)
    end
end

function mode_from_nr(g::GaussianBeam, nr, _)
    p = div(nr-1, 2*(g.max_m-1)+1)
    l = (nr-1) % (2 * (g.max_m-1) + 1) - (g.max_m-1)
    TEMMode(l, p)
end

function intersect(g1::GaussianBeam, g2::GaussianBeam)
    lb = [1e-6, 0]
    hb = [1, 2π]
    (lb, hb, (_, _) -> 1)
end
function intersect(g1::GaussianBeam, g2::CylindricalWaveguide)
    lb = [1e-6, 0]
    hb = [g2.r, 2π]
    (lb, hb, (_, _) -> 1)
end
function intersect(g1::CylindricalWaveguide, g2::GaussianBeam)
    lb = [1e-6, 0]
    hb = [g1.r, 2π]
    (lb, hb, (_, _) -> 1)
end
function contains(g::GaussianBeam, ρ, φ, z)
    g.z <= z <= g.z + g.length
end
jacobi_det(_::GaussianBeam, ρ, φ, z) = ρ

function integral_deps(g::GaussianBeam, z)
    (g, z)
end

w(g::GaussianBeam, z) = g.w_0 * sqrt(1 + (z / g.z_R)^2)
Rinv(g::GaussianBeam, z) = z / (g.z_R^2 + z^2)
ψ(g::GaussianBeam, z, mode::TEMMode) = (abs(mode.l)+2*mode.p + 1) * atan(z/g.z_R)
L(α, n, x) = Laguerre{α}(vcat(zeros(n), 1))(x)

β(g::GaussianBeam, _::Mode, dir::Direction) = Int(dir) * g.k

_A(g::GaussianBeam, ρ, φ, z, m::TEMMode, dir::Direction) = 1/w(g, z) *
        (ρ*sqrt(2)/w(g, z))^(abs(m.l)) *
        exp(-ρ^2/w(g, z)^2) *
        L(abs(m.l), m.p, 2*ρ^2/w(g, z)^2) *
        exp(-j * β(g, m, dir) * ρ^2 * Rinv(g, z) / 2) *
        exp(-j * m.l * φ) *
        exp(j*ψ(g, z, m))

E_spatial(g::GaussianBeam, ρ, φ, z, mode::TEMMode, dir::Direction) = begin
    z_rel = z - (dir == fwd ? (g.z + g.z_0) : (g.z + g.z_0 + g.length))
    @SVector [
        (g.μ / g.ε)^(1//4) * _A(g, ρ, φ, z_rel, mode, dir) * sin(φ),
        (g.μ / g.ε)^(1//4) * _A(g, ρ, φ, z_rel, mode, dir) * cos(φ),
        0,
    ]
end
E_freq(_::GaussianBeam, _::TEMMode, _::Direction) = [1, 1, 0]


H_spatial(g::GaussianBeam, ρ, φ, z, mode::TEMMode, dir::Direction) = begin
    z_rel = z - (dir == fwd ? (g.z + g.z_0) : (g.z + g.z_0 + g.length))
    @SVector [
        -(g.ε / g.μ)^(1//4) * _A(g, ρ, φ, z_rel, mode, dir) * cos(φ),
        (g.ε / g.μ)^(1//4) * _A(g, ρ, φ, z_rel, mode, dir) * sin(φ),
        0
    ]
end

H_freq(_::GaussianBeam, _::TEMMode, _::Direction) = [1, 1, 0]
