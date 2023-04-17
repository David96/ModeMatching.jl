module BesselFunctions
using SpecialFunctions, FunctionZeros, ForwardDiff, Roots

export besselj_prime, besselj_prime_zero

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

end;
