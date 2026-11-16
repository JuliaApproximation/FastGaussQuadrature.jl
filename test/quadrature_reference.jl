# These references use classical polynomials and shifted-parameter derivative
# identities, independently of the default quadrature iterations.
module QuadratureReference
using SpecialFunctions
using Test

function polynomials(family, n, a, b, x)
    p = [one(x)]
    n == 0 && return p
    if family == :jacobi
        ap, bp = a + 1, b + 1
        sumparameters = ap + bp
        push!(p, ((a - b) + sumparameters * x) / 2)
        for k in 2:n
            s, sm2 = 2(k - 1) + sumparameters, 2(k - 2) + sumparameters
            push!(
                p, (
                    (s - 1) * (s * sm2 * x + (a - b) * (a + b)) * p[end] -
                        2 * (k - 2 + ap) * (k - 2 + bp) * s * p[end - 1]
                ) /
                    (2k * (k - 2 + sumparameters) * sm2)
            )
        end
    elseif family == :laguerre
        push!(p, 1 + a - x)
        for k in 1:(n - 1)
            push!(p, ((2k + 1 + a - x) * p[end] - (k + a) * p[end - 1]) / (k + 1))
        end
    else
        push!(p, 2x)
        for k in 1:(n - 1)
            push!(p, 2x * p[end] - 2k * p[end - 1])
        end
    end
    return p
end

function derivative(family, n, a, b, x)
    family == :jacobi && return (n + a + b + 1) / 2 * polynomials(family, n - 1, a + 1, b + 1, x)[end]
    family == :laguerre && return -polynomials(family, n - 1, a + 1, b, x)[end]
    return 2n * polynomials(family, n - 1, a, b, x)[end]
end

function nodeweight(family, n, a, b, seed; bits)
    return setprecision(bits + 128) do
        a, b, x = BigFloat(a), BigFloat(b), BigFloat(seed)
        converged = false
        for _ in 1:30
            step = polynomials(family, n, a, b, x)[end] / derivative(family, n, a, b, x)
            x -= step
            # The degree-n recurrence has a rounding floor proportional to n.
            # The 128 guard bits keep it below the requested output precision.
            if abs(step) <= 8n * eps(BigFloat) * max(one(x), abs(x))
                converged = true
                break
            end
        end
        converged || error("Independent $family reference did not converge")
        d = derivative(family, n, a, b, x)
        w = if family == :jacobi
            c = exp(
                (a + b + 1) * log(BigFloat(2)) + loggamma(n + a + 1) +
                    loggamma(n + b + 1) - loggamma(n + a + b + 1) - loggamma(BigFloat(n + 1))
            )
            c / ((1 - x) * (1 + x) * d^2)
        elseif family == :laguerre
            gamma(n + a + 1) / factorial(big(n)) / (x * d^2)
        else
            BigFloat(2)^(n + 1) * factorial(big(n)) * sqrt(big(π)) / d^2
        end
        return x, w
    end
end

# Sturm counts also detect skipped or duplicated roots, which Newton refinement
# from the computed nodes alone would not detect.
function roots_below(family, n, a, b, x)
    p = polynomials(family, n, a, b, x)
    previous, changes = 1, 0
    for k in 1:n
        current = sign(family == :laguerre && isodd(k) ? -p[k + 1] : p[k + 1])
        if !iszero(current)
            changes += current != previous
            previous = current
        end
    end
    return n - changes
end

# Conservative weight-only envelope of the frozen native GST budget. At a
# root U''=0, so the coordinate sensitivity is |-a/(1-x)+b/(1+x)|. Reserve
# native recurrence/Taylor amplitude error proportional to n, logarithmic
# evaluation error, and normalization error separately (GST 2020, §2.1.1).
# P'(1)/P(1)=sum(1/(1-x_i)) gives the endpoint bounds below. Interlacing with
# P' bounds the ordinary roots after an extreme root is refined by the series.
# The entropy bound sum(-p_i*log(p_i))≤log(n) bounds normalization. Coefficient
# 8 was frozen before independent reserved-case validation; it is not an ulp
# target. Detailed pointwise budgets remain in the developer evidence.
function native_jacobi_weight_budget(n, a, b, reference)
    bits = precision(a)
    lower_series, upper_series = b < typeof(b)(-0.9), a < typeof(a)(-0.9)
    return setprecision(max(256, 2bits + 128)) do
        a, b, wr = BigFloat(a), BigFloat(b), BigFloat(reference)
        e = BigFloat(2)^(1 - bits)
        degree = BigFloat(n)
        D = degree * (degree + a + b + 1)
        lower, upper = -1 + 2 * (b + 1) / D, 1 - 2 * (a + 1) / D
        z, nu = min(a, b), degree + (a + b + 1) / 2
        G = abs(z) * max(abs(log1p(-lower^2)), abs(log1p(-upper^2))) +
            abs(a - z) * max(abs(log1p(-lower)), abs(log1p(-upper))) +
            abs(b - z) * max(abs(log1p(lower)), abs(log1p(upper)))
        KX = (lower_series || upper_series) ? 32 * (abs(a) + abs(b)) : zero(e)
        if n > lower_series + upper_series
            Dderivative = (degree - 1) * (degree + a + b + 2)
            lo = lower_series ? -1 + 2 * (b + 2) / Dderivative : lower
            hi = upper_series ? 1 - 2 * (a + 2) / Dderivative : upper
            KX = max(
                KX, 8 * (abs(a) + abs(b)) + 8 * (degree / nu) *
                    (abs(a) * sqrt((1 + hi) / (1 - hi)) + abs(b) * sqrt((1 - lo) / (1 + lo)))
            )
        end
        gamma_a, gamma_b, gamma_sum = loggamma(a + 1), loggamma(b + 1), loggamma(a + b + 2)
        logpower = (a + b + 1) * log(BigFloat(2))
        logmass = logpower + gamma_a + gamma_b - gamma_sum
        mass = 8 * (1 + abs(logpower) + abs(gamma_a) + abs(gamma_b) + abs(gamma_sum))
        local_error = 16 * (degree + 1) + KX + 8 * (1 + G + abs(log(wr) - logmass))
        mean_error = 16 * (degree + 1) + KX + 8 * (1 + G + log(degree))
        denominator = 1 - (degree - 1) * e / 2
        allowance = e * (mass + local_error + mean_error + 8) +
            (degree - 1) * e / (2 * denominator)
        denominator > 0 && allowance < BigFloat(1) / 8 || error("Native GST budget is inconclusive")
        return allowance
    end
end

function check(family, n, a, b, x, w; indices = eachindex(x))
    bits = precision(BigFloat)
    tolerance = BigFloat(2)^(-bits + 16)
    @test length(x) == length(w) == n
    @test eltype(x) == eltype(w) == BigFloat
    @test all(z -> precision(z) == bits, x)
    @test all(z -> precision(z) == bits, w)
    @test issorted(x) && all(>(0), diff(x))
    @test all(isfinite, x) && all(isfinite, w) && all(>(0), w)
    @testset "$family, n: $n, i: $i" for i in indices
        xr, wr = nodeweight(family, n, a, b, x[i]; bits)
        @test isapprox(x[i], xr; atol = tolerance, rtol = tolerance)
        # Relative error is essential for the tiny outer Laguerre weights.
        @test isapprox(w[i], wr; atol = 0, rtol = tolerance)
        setprecision(bits + 128) do
            aa, bb, z = BigFloat(a), BigFloat(b), BigFloat(x[i])
            delta = BigFloat(2)^(-bits ÷ 2) * max(one(z), abs(z))
            @test roots_below(family, n, aa, bb, z - delta) == i - 1
            @test roots_below(family, n, aa, bb, z + delta) == i
        end
    end
    return nothing
end
end
