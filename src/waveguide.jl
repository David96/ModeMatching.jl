using HCubature

"""
    Waveguide API:
    Every waveguide implementation is expected to provide:
    - E_spatial(g::WaveguideImpl, x, y, mode::Mode)::Vector
    - E_freq(g::WaveguideImpl, mode::Mode)::Vector
    - H_spatial(g::WaveguideImpl, x, y, mode::Mode)::Vector
    - E_freq(g::WaveguideImpl, mode::Mode)::Vector
    returning three dimensional vectors of the spatial and frequency dependent parts
    of the E and H fields of the Eigenmodes.
    - z_dep(g::WaveuideImpl, mode::Mode, z)::Complex
    returning the common z dependence.
    - intersect(g1::WaveguideImpl, g2::WaveguideImpl)::Tuple
    returning (lowerbound, higherbound, mask) of the two dimensional intersection
    of the two waveguides.

    int_ExHy(g1::WaveguideImpl, g2::WaveuideImpl, mode1::Mode, mode2::Mode) and
    int_EyHx(g1::WaveguideImpl, g2::WaveuideImpl, mode1::Mode, mode2::Mode)
    can be used to provide analytical solutions to the overlap integrals.
"""
abstract type Waveguide end
abstract type Mode end

function integral(g1::Waveguide, g2::Waveguide, f)
    (lb, hb, mask) = intersect(g1, g2)
    (I, _) = hcubature(r -> f(r[1], r[2]) * mask(r[1], r[2]),
                       lb, hb;
                       rtol=1e-4,
                       #atol=1e-12,
                       maxevals=Integer(1e6))
    I
end
function int_ExHy(g1::Waveguide, g2::Waveguide, mode1::Mode, mode2::Mode)
    E1 = (x, y) -> E_spatial(g1, x, y, mode1)
    H2 = (x, y) -> H_spatial(g2, x, y, mode2)
    integral(g1, g2, (x, y) -> E1[1] * H2[2])
end
function int_EyHx(g1::Waveguide, g2::Waveguide, mode1::Mode, mode2::Mode)
    E1 = (x, y) -> E_spatial(g1, x, y, mode1)
    H2 = (x, y) -> H_spatial(g2, x, y, mode2)
    integral(g1, g2, (x, y) -> E1[2] * H2[1])
end
function scalar(g1::Waveguide, g2::Waveguide, mode1::Mode, mode2::Mode; norm=true)
    # Remove if in doubt of orthonormal system!
    if norm# && false
        if g1 == g2 && mode1 == mode2
            return Complex(1)
        elseif g1 == g2
            return Complex(0)
        end
    end
    I1 = int_ExHy(g1, g2, mode1, mode2)
    I2 = int_EyHx(g1, g2, mode1, mode2)

    E1_f = E_freq(g1, mode1)
    H2_f = H_freq(g2, mode2)
    ((norm ? C(g1, mode1) * C(g2, mode2) : 1) *
     (E1_f[1] * H2_f[2] * I1 - E1_f[2] * H2_f[1] * I2))
end
Cs = Dict()
lockCs = ReentrantLock()
function C(g::Waveguide, mode::Mode)
    hash = "$g;$mode"
    lock(lockCs)
    if !(hash in keys(Cs))
        Cm = scalar(g, g, mode, mode; norm=false)
        Cs[hash] = 1 / sqrt(Complex(Cm))
    end
    unlock(lockCs)
    Cs[hash]
end

#=
 General mode functions combining normalisation, frequency dependence, orthogonal component and
 z dependence
=#
"""
    E(g::Waveguide, x, y, z, mode::Mode[; norm=true])

    Return the vector of the E field of a certain `mode` at position `x`, `y`, `z`.
    `norm` specifies whether the modes should be normalized in regards to `scalar` (default: true)
"""
function E(g::Waveguide, x, y, z, mode::Mode; norm=true)
    (norm ? C(g, mode) : 1) * E_freq(g, mode) .* E_spatial(g, x, y, mode) * z_dep(g, mode, z)
end

"""
    H(g::Waveguide, x, y, z, mode::Mode[; norm=true])

    Return the vector of the H field of a certain `mode` at position `x`, `y`, `z`.
    `norm` specifies whether the modes should be normalized in regards to `scalar` (default: true)
"""
function H(g::Waveguide, x, y, z, mode::Mode; norm=true)
    (norm ? C(g, mode) : 1) * H_freq(g, mode) .* H_spatial(g, x, y, mode) * z_dep(g, mode, z)
end
