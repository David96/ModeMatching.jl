using LinearAlgebra, Base.Iterators, FLoops

struct WaveguideSetup
    waveguides::Vector{Waveguide}
    n_TE::Integer
    n_TM::Integer
    max_m::Integer
    a_1::Vector
    n_total::Integer
    n_layers::Integer

    function WaveguideSetup(waveguides, n_TE, n_TM, max_m, a_1)
        new(waveguides, n_TE, n_TM, max_m, a_1, n_TE+n_TM, length(waveguides))
    end
end

function power(g::Waveguide, mode::Mode)
    E1(x, y, z) = E(g, x, y, z, mode)
    H2(x, y, z) = H(g, x, y, z, mode)
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

function Es_right(s::WaveguideSetup, region, a_i, x, y, z)
    g = s.waveguides[region]
    [a_i[region][mode] *
     E(g, x, y, +(z - g.z), mode_from_nr(g, mode, s.n_TE, s.n_TM, s.max_m))
     for mode=1:(s.n_TM + s.n_TE)]
end
function Es_left(s::WaveguideSetup, region, b_i, x, y, z)
    g = s.waveguides[region]
    znext = g.z + g.length
    [(b_i[region][mode] *
     (E(g, x, y, -(z - znext), mode_from_nr(g, mode, s.n_TE, s.n_TM, s.max_m))
     .* ([1, 1, -1]))) # left moving -> invert z coordinate
     for mode=1:s.n_total]
end
function E_tot(s::WaveguideSetup, r, a_i, b_i, x, y, z)
    if r == s.n_layers
        return sum(Es_right(s, r, a_i, x, y, z))
    end
    sum(Es_right(s, r, a_i, x, y, z)) + sum(Es_left(s, r, b_i, x, y, z))
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
    [0, 0, 0]
end
function Hs_right(s::WaveguideSetup, region, a_i, x, y, z)
    g = s.waveguides[region]
    [a_i[region][mode] *
     H(g, x, y, +(z - g.z), mode_from_nr(g, mode, s.n_TE, s.n_TM, s.max_m))
     for mode=1:(s.n_TM + s.n_TE)]
end
function Hs_left(s::WaveguideSetup, region, b_i, x, y, z)
    g = s.waveguides[region]
    znext = g.z + g.length
    [(b_i[region][mode] *
      (H(g, x, y, -(z - znext), mode_from_nr(g, mode, s.n_TE, s.n_TM, s.max_m))
       .* ([-1, -1, 1]))) # left moving -> invert z coordinate
     for mode=1:s.n_total]
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
    [0, 0, 0]
end

function t_r_ab(s::WaveguideSetup)
    n_layers = length(s.waveguides)
    @assert 1 <= n_layers <= 3
    if n_layers == 1
        return [Diagonal(fill(1, s.n_total)), zeros(s.n_total, s.n_total)]
    end
    (G_12, G_21) = G_12_21(s.waveguides[1], s.waveguides[2], s.waveguides[2].z,
                           s.n_TE, s.n_TM, s.max_m)
    if n_layers == 2
        t_ab = Array{Matrix}(undef, 2, 2)
        r_ab = Array{Matrix}(undef, 2, 2)

        # Parameters of interface from 1st to 2nd layer
        t_ab[1, 2], r_ab[1, 2] = t_r_12(G_12, G_21)
        t_ab[2, 1], r_ab[2, 1] = t_r_12(G_21, G_12)
        return t_ab, r_ab
    end
    d_2 = s.waveguides[3].z - s.waveguides[2].z
    P_2r = P_qr(d_2, s.waveguides[2], s.n_TE, s.n_TM, s.max_m)
    P_2l = P_ql(d_2, s.waveguides[2], s.n_TE, s.n_TM, s.max_m)
    G_23, G_32 = G_12_21(s.waveguides[2], s.waveguides[3], s.waveguides[3].z,
                         s.n_TE, s.n_TM, s.max_m)
    t_r_ab(G_12, G_21, G_23, G_32, P_2r, P_2l)
end

function calc_a_i(s::WaveguideSetup, t_ab, r_ab)
    n_layers = length(s.waveguides)
    @assert 1 <= n_layers <= 3
    if n_layers == 1
        return [s.a_1]
    elseif n_layers == 2
        aprime_1 = P_qr(s.waveguides[2].z, s.waveguides[1], s.n_TE, s.n_TM, s.max_m) * s.a_1
        [s.a_1, t_ab[1, 2] * aprime_1]
    else
        aprime_1 = P_qr(s.waveguides[2].z, s.waveguides[1], s.n_TE, s.n_TM, s.max_m) * s.a_1
        bprime_3 = zeros(s.n_TE + s.n_TM)
        d_2 = s.waveguides[3].z - s.waveguides[2].z
        P_2r = P_qr(d_2, s.waveguides[2], s.n_TE, s.n_TM, s.max_m)
        P_2l = P_ql(d_2, s.waveguides[2], s.n_TE, s.n_TM, s.max_m)
        a_2 = inv(I - r_ab[2, 1] * P_2l * r_ab[2, 3] * P_2r) *
                (t_ab[1, 2] * aprime_1 + r_ab[2, 1] * P_2l * t_ab[3, 2] * bprime_3)
        a_3 = t_ab[1, 3] * aprime_1
        [s.a_1, a_2, a_3]
    end
end

function calc_b_i(s::WaveguideSetup, t_ab, r_ab)
    n_layers = length(s.waveguides)
    @assert 1 <= n_layers <= 3
    if n_layers == 1
        return [zeros(s.n_total)]
    elseif n_layers == 2
        aprime_1 = P_qr(s.waveguides[2].z, s.waveguides[1], s.n_TE, s.n_TM, s.max_m) * s.a_1
        [r_ab[1, 2] * aprime_1, zeros(s.n_total)]
    else
        aprime_1 = P_qr(s.waveguides[2].z, s.waveguides[1], s.n_TE, s.n_TM, s.max_m) * s.a_1
        bprime_3 = zeros(s.n_TE + s.n_TM)
        d_2 = s.waveguides[3].z - s.waveguides[2].z
        P_2r = P_qr(d_2, s.waveguides[2], s.n_TE, s.n_TM, s.max_m)
        P_2l = P_ql(d_2, s.waveguides[2], s.n_TE, s.n_TM, s.max_m)
        b_2 = inv(I - r_ab[2, 3] * P_2r * r_ab[2, 1] * P_2l) *
             (r_ab[2, 3] * P_2r * t_ab[1, 2] * aprime_1 + t_ab[3, 2] * bprime_3)
        [r_ab[1, 3] * aprime_1, b_2, zeros(s.n_TM + s.n_TE)]
    end
end

function G_12_21(g1::Waveguide, g2::Waveguide, z, n_TE, n_TM, max_m)
    @assert g1.z + g1.length == g2.z "Waveguides must have an interface!"
    n_total = n_TE + n_TM
    G_12 = Array{Complex}(undef, n_total, n_total)
    G_21 = Array{Complex}(undef, n_total, n_total)
    @floop for (mode1, mode2) = product(1:n_total, 1:n_total)
        m1 = mode_from_nr(g1, mode1, n_TE, n_TM, max_m)
        m2 = mode_from_nr(g2, mode2, n_TE, n_TM, max_m)
        #println("Matching $m1 with $m2...")
        I = scalar(g1, g2, z, m1, m2)
        G_12[mode1, mode2] = I
        I = scalar(g2, g1, z, m1, m2)
        G_21[mode1, mode2] = I
    end
    (G_12, G_21)
end

P_qr(d, g, n_TE, n_TM, max_m) = Diagonal([propagation(g, mode_from_nr(g, mode, n_TE, n_TM, max_m), d) for mode=1:(n_TE + n_TM)])
P_ql(d, g, n_TE, n_TM, max_m) = P_qr(d, g, n_TE, n_TM, max_m)

function t_r_12(G_12, G_21)
    t = 2 * inv(transpose(G_21) + G_12)
    r = 1/2 * (transpose(G_21) - G_12) * t
    (t, r)

    #t = 2 * inv(transpose(G_12) + G_21)
    #r = 1/2 * (transpose(G_12) - G_21) * t
    #(t, r)
end
function t_r_ab(G_12, G_21, G_23, G_32, P_2r, P_2l)
    t_ab = Array{Matrix}(undef, 3, 3)
    r_ab = Array{Matrix}(undef, 3, 3)

    # Parameters of interface from 1st to 2nd layer
    (t_ab[1, 2], r_ab[1, 2]) = t_r_12(G_12, G_21)
    (t_ab[2, 1], r_ab[2, 1]) = t_r_12(G_21, G_12)

    # Parameters of interface from 2nd to 3rd layer
    (t_ab[2, 3], r_ab[2, 3]) = t_r_12(G_23, G_32)
    (t_ab[3, 2], r_ab[3, 2]) = t_r_12(G_32, G_23)

    # Implements Numerical methods in photonics 6.109-6.112
    #display(abs.(eigvals(r_ab[2, 1] * P_2l * r_ab[2, 3] * P_2r)))
    tmp_inv = inv(I - r_ab[2, 1] * P_2l * r_ab[2, 3] * P_2r)
    r_ab[1, 3] = r_ab[1, 2] + t_ab[2, 1] * P_2l * r_ab[2, 3] * P_2r * tmp_inv * t_ab[1, 2]
    t_ab[1, 3] = t_ab[2, 3] * P_2r * tmp_inv * t_ab[1, 2]
    r_ab[3, 1] = r_ab[3, 2] + t_ab[2, 3] * P_2r * tmp_inv * r_ab[2, 1] * P_2l * t_ab[3, 2]
    t_ab[3, 1] = t_ab[2, 1] * P_2l * t_ab[3, 2] + t_ab[2, 1] * P_2l * r_ab[2, 3] * 
    P_2r * tmp_inv * r_ab[2, 1] * P_2l * t_ab[3, 2]

    (t_ab, r_ab)
end
r_1_qp1(r_1q, t_q1, P_q, r_q_qp1) = r_1q + t_q1 * P_q * r_q_qp1 * P_q * inv(I - r_q1 * P_q * r_q_qp1 * P_q) * t_1q
t_qp1_1(t_q1, P_q, t_qp1_q, r_q_qp1) = t_q1 * P_q * t_qp1_q + t_q1 * P_q * r_q_qp1 * P_q *
            inv(I - r_q1 * P_q * r_q_qp1 * P_q) * r_q1 * P_q * t_qp1_q
