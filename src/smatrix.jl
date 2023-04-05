using LinearAlgebra, Base.Iterators, FLoops

struct WaveguideSetup
    waveguides::Vector{Waveguide}
    a_1::Vector
    n_modes::Vector{Int}
    n_layers::Integer

    function WaveguideSetup(waveguides, n_modes_1, a_1)
        n_modes = Vector{Int}(undef, length(waveguides))
        n_modes[1] = n_modes_1
        for (i, g)=enumerate(waveguides[2:end])
            n_modes[i+1] = round(n_modes_1 * g.A / waveguides[1].A, RoundUp)
        end
        new(waveguides, a_1, n_modes, length(waveguides))
    end
end

function power(g::Waveguide, mode::Mode, dir::Direction)
    E1(x, y, z) = E(g, x, y, z, mode, dir)
    H2(x, y, z) = H(g, x, y, z, mode, dir)
    power(g, E1, H2, 0)
end
function power(g::Waveguide, E::Function, H::Function, z)
    (lb, hb, mask) = ModeMatching.intersect(g, g)
    (I, _) = hcubature(r -> (E(r[1], r[2], z)[1] * conj(H(r[1], r[2], z)[2]) -
                             E(r[1], r[2], z)[2] * conj(H(r[1], r[2], z)[1])) *
                            mask(r[1], r[2]) * jacobi_det(g, r[1], r[2], z),
                   lb, hb;
                   rtol=1e-8)
    I / 2
end

function get_region(setup::WaveguideSetup, z)
    for (i, g) in enumerate(setup.waveguides)
        if g.z <= z < g.z + g.length
            return i
        end
    end
    NaN
end

const z_pm = -1
function Es_right(s::WaveguideSetup, region, a_i, x, y, z)
    g = s.waveguides[region]
    [a_i[region][mode] .*
     E(g, x, y, z, mode_from_nr(g, mode, s.n_modes[region]), fwd)
     for mode=1:(s.n_modes[region])]
end
function Es_left(s::WaveguideSetup, region, b_i, x, y, z)
    g = s.waveguides[region]
    [b_i[region][mode] .*
     (E(g, x, y, z, mode_from_nr(g, mode, s.n_modes[region]), bck)#= .* [1, 1, -1]=#)
     for mode=1:s.n_modes[region]]
end
function E_tot(s::WaveguideSetup, r, a_i, b_i, x, y, z)
    if r == s.n_layers
        return sum(Es_right(s, r, a_i, x, y, z))
    end
    SVector{3, ComplexF64}(sum(Es_right(s, r, a_i, x, y, z)) + sum(Es_left(s, r, b_i, x, y, z)))
end
function E(s::WaveguideSetup, a_i, b_i, x, y, z)
    r = get_region(s, z)
    if r <= length(s.waveguides)
        g = s.waveguides[r]
        if contains(g, x, y, z)
            return E_tot(s, r, a_i, b_i, x, y, z)
        end
    end
    # Outside the boundaries, the field is 0 by definition
    @SVector [0, 0, 0]
end
function Hs_right(s::WaveguideSetup, region, a_i, x, y, z)
    g = s.waveguides[region]
    [a_i[region][mode] .*
     H(g, x, y, z, mode_from_nr(g, mode, s.n_modes[region]), fwd)
     for mode=1:(s.n_modes[region])]
end
function Hs_left(s::WaveguideSetup, region, b_i, x, y, z)
    g = s.waveguides[region]
    [b_i[region][mode] .*
     (H(g, x, y, z, mode_from_nr(g, mode, s.n_modes[region]), bck)#= .* [-1, -1, 1]=#)
     for mode=1:s.n_modes[region]]
end
function H_tot(s::WaveguideSetup, r, a_i, b_i, x, y, z)
    if r == s.n_layers
        return sum(Hs_right(s, r, a_i, x, y, z))
    end
    sum(Hs_right(s, r, a_i, x, y, z)) + sum(Hs_left(s, r, b_i, x, y, z))
end
function H(s::WaveguideSetup, a_i, b_i, x, y, z)
    r = get_region(s, z)
    if r <= length(s.waveguides)
        g = s.waveguides[r]
        if contains(g, x, y, z)
            return H_tot(s, r, a_i, b_i, x, y, z)
        end
    end
    # Outside the boundaries, the field is 0 by definition
    @SVector [0, 0, 0]
end

function P_q(g, n_modes)
    Diagonal([
        begin
            m = mode_from_nr(g, mode, n_modes)
            propagation(g, m, fwd, g.z + g.length)
        end for mode=1:(n_modes)])
end

function t_r_ab(s::WaveguideSetup)
    n = length(s.waveguides)
    if n == 1
        return [Diagonal(fill(1, s.n_modes[1])), zeros(s.n_modes[1], s.n_modes[1])]
    end
    g1 = s.waveguides[1]
    g2 = s.waveguides[2]
    (G_12, G_21) = G_12_21(g1, g2, s.n_modes[1], s.n_modes[2])

    t = Array{Matrix{ComplexF64}}(undef, n, n)
    r = Array{Matrix{ComplexF64}}(undef, n, n)
    t[1, 2], r[1, 2] = t_r_12(s, 1, 2, G_12, G_21)
    t[2, 1], r[2, 1] = t_r_12(s, 2, 1, G_21, G_12)

    if n == 2
        return (t, r)
    end

    for q=2:n-1
        g1 = s.waveguides[q]
        g2 = s.waveguides[q+1]
        G_12, G_21 = G_12_21(g1, g2, s.n_modes[q], s.n_modes[q+1])
        t[q, q+1], r[q, q+1] = t_r_12(s, q, q+1, G_12, G_21)
        t[q+1, q], r[q+1, q] = t_r_12(s, q+1, q, G_21, G_12)

        P_qr = P_q(g1, s.n_modes[q])
        P_ql = P_qr
        A = r[q, 1] * P_ql * r[q, q+1] * P_qr
        tmp_inv = inv(I - A)
        r[1, q+1] = r[1, q] + t[q, 1] * P_ql * r[q, q+1] * P_qr * tmp_inv * t[1, q]
        t[1, q+1] = t[q, q+1] * P_qr * tmp_inv * t[1, q]
        r[q+1, 1] = r[q+1, q] + t[q, q+1] * P_qr * tmp_inv * r[q, 1] * P_ql * t[q+1, q]
        t[q+1, 1] = t[q, 1] * P_ql * t[q+1, q] + t[q, 1] * P_ql * r[q, q+1] * P_qr * tmp_inv *
            r[q, 1] * P_ql * t[q+1, q]
    end

    for q=n-1:-1:2
        g1 = s.waveguides[q]
        P_qr = P_q(g1, s.n_modes[q])
        P_ql = P_qr
        A = r[q, q-1] * P_ql * r[q, n] * P_qr
        tmp_inv = inv(I - A)
        r[q-1, n] = r[q-1, q] + t[q, q-1] * P_ql * r[q, n] * P_qr * tmp_inv * t[q-1, q]
        t[q-1, n] = t[q, n] * P_qr * tmp_inv * t[q-1, q]
        r[n, q-1] = r[n, q] + t[q, n] * P_qr * tmp_inv * r[q, q-1] * P_ql * t[n, q]
        t[n, q-1] = t[q, q-1] * P_ql * t[n, q] + t[q, q-1] * P_ql * r[q, n] * tmp_inv *
            r[q, q-1] * P_ql * t[n, q]
    end
    (t, r)
end

function calc_a_i(s::WaveguideSetup, t_ab, r_ab)
    n_layers = length(s.waveguides)
    if n_layers == 1
        return [s.a_1]
    elseif n_layers == 2
        aprime_1 = P_q(s.waveguides[1], s.n_modes[1]) * s.a_1
        [s.a_1, t_ab[1, 2] * aprime_1]
    else
        aprime_1 = P_q(s.waveguides[1], s.n_modes[1]) * s.a_1
        bprime_n = zeros(s.n_modes[end])
        a_qs = Vector{Vector{ComplexF64}}(undef, s.n_layers)
        a_qs[1] = s.a_1
        for q=2:s.n_layers-1
            P_qr = P_q(s.waveguides[q], s.n_modes[q])
            P_ql = P_qr
            a_qs[q] = inv(I - r_ab[q, 1] * P_ql * r_ab[q, s.n_layers] * P_qr) *
                    (t_ab[1, q] * aprime_1 + r_ab[q, 1] * P_ql * t_ab[s.n_layers, q] * bprime_n)
        end
        a_qs[end] = t_ab[1, s.n_layers] * aprime_1
        a_qs
    end
end

function calc_b_i(s::WaveguideSetup, t_ab, r_ab)
    n_layers = length(s.waveguides)
    if n_layers == 1
        return [zeros(s.n_modes[1])]
    elseif n_layers == 2
        aprime_1 = P_q(s.waveguides[1], s.n_modes[1]) * s.a_1
        [r_ab[1, 2] * aprime_1, zeros(s.n_modes[2])]
    else
        aprime_1 = P_q(s.waveguides[1], s.n_modes[1]) * s.a_1
        bprime_n = zeros(s.n_modes[end])
        b_qs = Vector{Vector{ComplexF64}}(undef, s.n_layers)
        n = s.n_layers
        b_qs[1] = r_ab[1, n] * aprime_1
        for q=2:n-1
            P_qr = P_q(s.waveguides[q], s.n_modes[q])
            P_ql = P_qr
            b_qs[q] = inv(I - r_ab[q, n] * P_qr * r_ab[q, 1] * P_ql) *
                (r_ab[q, n] * P_qr * t_ab[1, q] * aprime_1 + t_ab[n, q] * bprime_n)
        end
        b_qs[end] = zeros(s.n_modes[end])
        b_qs
    end
end

function G_12_21(g1::Waveguide, g2::Waveguide, n_modes1, n_modes2)
    @assert g1.z + g1.length ≈ g2.z "Waveguides must have an interface!"
    G_12 = Array{ComplexF64}(undef, n_modes1, n_modes2)
    G_21 = Array{ComplexF64}(undef, n_modes2, n_modes1)
    @floop for (mode1, mode2) = product(1:n_modes1, 1:n_modes2)
        m1 = mode_from_nr(g1, mode1, n_modes1)
        m2 = mode_from_nr(g2, mode2, n_modes2)
        I = scalar(g1, g2, g2.z, m1, m2)
        G_12[mode1, mode2] = I
        I = scalar(g2, g1, g2.z, m2, m1)
        G_21[mode2, mode1] = I
    end
    (G_12, G_21)
end

function t_r_12(s::WaveguideSetup, r1, r2, G_12, G_21)
    g1 = s.waveguides[r1]
    g2 = s.waveguides[r2]
    lb1, hb1, _ = intersect(g1, g1)
    lb2, hb2, _ = intersect(g2, g2)
    boundary_enlargement = any(lb2 .< lb1) || any(hb2 .> hb1)
    G = boundary_enlargement ? G_12 : G_21
    tG = transpose(G)
    S1 = Diagonal([mode_from_nr(g1, i, s.n_modes[r1]) isa TMMode ? -1. : 1.
                   for i=1:s.n_modes[r1]])
    if boundary_enlargement
        t = 2 * inv(I + tG * G) * tG
        r = S1 * (I - G * t)
        return (t, r)
    else
        t = 2 * inv(I + G * tG) * G
        r = S1 * (tG * t - I)
        return (t, r)
    end
end
