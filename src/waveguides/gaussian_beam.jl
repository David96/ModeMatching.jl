using SpecialPolynomials

struct GaussianBeam <: Waveguide
    k::AbstractFloat
    f::AbstractFloat
    ω::AbstractFloat
    ε::AbstractFloat
    μ::AbstractFloat
    n::AbstractFloat
    w_0::AbstractFloat
    z_R::AbstractFloat
    x::AbstractFloat
    y::AbstractFloat
    z::AbstractFloat
    length::AbstractFloat

    function GaussianBeam(x, y, z, w_0, length, μ, ε, f)
        new(2π * f * sqrt(μ * ε) / c, f, 2π * f, ε*ε0, μ*μ0, sqrt(ε), w_0,
            π*w_0^2*sqrt(ε)/(c/f), x, y, z, length)
    end
end

function mode_from_nr(_::GaussianBeam, nr::Integer,
        n_TE::Integer, n_TM::Integer, max_m::Integer)
    p = div(nr-1, 2*max_m+1)
    l = (nr-1) % (2 * max_m + 1) - max_m
    TEMMode(l, p)
end

function intersect(g1::GaussianBeam, g2::GaussianBeam)
    lb = [1e-6, 0]
    hb = [50e-2, 2π]
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

w(g::GaussianBeam, z) = g.w_0 * sqrt(1 + (z / g.z_R)^2)
Rinv(g::GaussianBeam, z) = z / (g.z_R^2 + z^2)
ψ(g::GaussianBeam, z, mode::TEMMode) = (abs(mode.l)+2*mode.p + 1) * atan(z/g.z_R)
L(α, n, x) = Laguerre{α}(vcat(zeros(n), 1))(x)

β(g::GaussianBeam, _::Mode, dir::Direction) = Int(dir) * g.k

E_spatial(g::GaussianBeam, ρ, φ, z, mode::TEMMode, dir::Direction) = begin
    z_rel = z - (dir == fwd ? g.z : (g.z + g.length))
    @SVector [
        0,
        1 / w(g, z_rel) *
          (ρ * sqrt(2) / w(g, z_rel))^(abs(mode.l)) * exp(-ρ^2/w(g, z_rel)^2) *
          L(abs(mode.l), mode.p, (2*ρ^2 / w(g, z_rel)^2) *
          exp(-1im * g.k * ρ^2 * 2 * Rinv(g, z_rel))) *
          exp(-1im*mode.l*φ) * exp(1im*ψ(g, z_rel, mode)),
        0,
    ]
end
E_freq(_::GaussianBeam, _::TEMMode, _::Direction) = [0, 1, 0]


H_spatial(g::GaussianBeam, ρ, φ, z, mode::TEMMode, dir::Direction) = begin
    z_rel = z - (dir == fwd ? g.z : (g.z + g.length))
    [
          -sqrt(g.ε / g.μ) / w(g, z_rel) *
          (ρ * sqrt(2) / w(g, z_rel))^(abs(mode.l)) * exp(-ρ^2/w(g, z_rel)^2) *
          L(abs(mode.l), mode.p, (2 * ρ^2 / w(g, z_rel)^2) *
          exp(-1im * g.k * ρ^2 * 2 * Rinv(g, z_rel))) *
          exp(-1im*mode.l*φ) * exp(1im*ψ(g, z_rel, mode)),
        0,
        0
    ]
end

H_freq(_::GaussianBeam, _::TEMMode, _::Direction) = [1, 0, 0]
