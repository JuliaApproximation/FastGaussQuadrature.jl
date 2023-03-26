var documenterSearchIndex = {"docs":
[{"location":"misc/#Misc.","page":"Misc.","title":"Misc.","text":"","category":"section"},{"location":"misc/#About-the-logo","page":"Misc.","title":"About the logo","text":"","category":"section"},{"location":"misc/","page":"Misc.","title":"Misc.","text":"(Image: )","category":"page"},{"location":"misc/","page":"Misc.","title":"Misc.","text":"The center points of circles are generated from nodes(x-axis) and weights(y-axis) of Gaussian quadrature.\nNumber of points : n = 15\ntextcolorCB3C33texttextbfRed : gausslegendre(n)\ntextcolor389826texttextbfGreen : gausschebyshev(n, 3)\ntextcolor9558B2texttextbfPurple : gaussjacobi(n, 5/2, 1/2)\nThese colors are from julia-logo-graphics.\nThe logo is generated with Luxor.jl, and its source code is here.","category":"page"},{"location":"misc/#Other-docstrings-from-private-methods","page":"Misc.","title":"Other docstrings from private methods","text":"","category":"section"},{"location":"misc/","page":"Misc.","title":"Misc.","text":"Modules = [FastGaussQuadrature]\nPrivate = true\nPublic = false","category":"page"},{"location":"misc/#FastGaussQuadrature.AIRY_ROOTS","page":"Misc.","title":"FastGaussQuadrature.AIRY_ROOTS","text":"The first 11 roots of the Airy function in Float64 precision https://mathworld.wolfram.com/AiryFunctionZeros.html\n\nExamples\n\njulia> zeros = airy.(FastGaussQuadrature.AIRY_ROOTS);\n\njulia> all(zeros .< 1e-14)\ntrue\n\n\n\n\n\n","category":"constant"},{"location":"misc/#FastGaussQuadrature.BESSELJ0_ROOTS","page":"Misc.","title":"FastGaussQuadrature.BESSELJ0_ROOTS","text":"First twenty roots of Bessel function J_0 in Float64. https://mathworld.wolfram.com/BesselFunctionZeros.html\n\nExamples\n\njulia> zeros = besselj0.(FastGaussQuadrature.BESSELJ0_ROOTS);\n\njulia> all(zeros .< 1e-14)\ntrue\n\n\n\n\n\n","category":"constant"},{"location":"misc/#FastGaussQuadrature.BESSELJ1_ON_BESSELJ0_ROOTS","page":"Misc.","title":"FastGaussQuadrature.BESSELJ1_ON_BESSELJ0_ROOTS","text":"Values of Bessel function J_1 on first ten roots of Bessel function J_0.\n\nExamples\n\njulia> roots = approx_besselroots(0,10);\n\njulia> (besselj1.(roots)).^2 ≈ FastGaussQuadrature.BESSELJ1_ON_BESSELJ0_ROOTS\ntrue\n\n\n\n\n\n","category":"constant"},{"location":"misc/#FastGaussQuadrature.PIESSENS_C","page":"Misc.","title":"FastGaussQuadrature.PIESSENS_C","text":"Coefficients of Chebyshev series approximations for the zeros of the Bessel functions\n\nj_nu s approx sum_k^n C_ks T_k(fracnu-23)\n\nwhere j_nu s is a s-th zero of Bessel function J_nu, T_k are Chebyshev polynomials and C_k s is the coefficients.\n\n\n\n\n\n","category":"constant"},{"location":"misc/#FastGaussQuadrature.evalLaguerreRec-Tuple{Any, Any, Any}","page":"Misc.","title":"FastGaussQuadrature.evalLaguerreRec","text":"Evaluate the orthonormal associated Laguerre polynomial with positive leading coefficient, as well as its derivative, in the point x using the recurrence relation.\n\n\n\n\n\n","category":"method"},{"location":"misc/#FastGaussQuadrature.feval_asy1-Tuple{Integer, Float64, Float64, AbstractVector, Any}","page":"Misc.","title":"FastGaussQuadrature.feval_asy1","text":"Evaluate the interior asymptotic formula at x = cos(t). Assumption:\n\nlength(t) == n ÷ 2\n\n\n\n\n\n","category":"method"},{"location":"misc/#FastGaussQuadrature.gausslaguerre_GW-Tuple{Any, Any}","page":"Misc.","title":"FastGaussQuadrature.gausslaguerre_GW","text":"Calculate Gauss-Laguerre nodes and weights from the eigenvalue decomposition of the Jacobi matrix.\n\n\n\n\n\n","category":"method"},{"location":"misc/#FastGaussQuadrature.gausslaguerre_asy-Tuple{Integer, Any}","page":"Misc.","title":"FastGaussQuadrature.gausslaguerre_asy","text":"Compute the Gauss-Laguerre rule using explicit asymptotic expansions for the nodes and weights. Optional parameters are:\n\nreduced: compute a reduced quadrature rule, discarding all points and weights as soon as the weights underflow\nT: the order of the expansion. Set T=-1 to determine the order adaptively depending on the size of the terms in the expansion\nrecompute: if a crude measure of the error is larger than a tolerance, the point and weight are recomputed using the (slower) recursion+newton approach, yielding more reliable accurate results.\n\n\n\n\n\n","category":"method"},{"location":"misc/#FastGaussQuadrature.gausslaguerre_rec-Tuple{Any, Any}","page":"Misc.","title":"FastGaussQuadrature.gausslaguerre_rec","text":"Compute Gauss-Laguerre rule based on the recurrence relation, using Newton iterations on an initial guess.\n\n\n\n\n\n","category":"method"},{"location":"besselroots/#Roots-of-Bessel-function","page":"Roots of Bessel function","title":"Roots of Bessel function","text":"","category":"section"},{"location":"besselroots/","page":"Roots of Bessel function","title":"Roots of Bessel function","text":"Since SpecialFunctions.jl doesn't have a method to calculate roots of Bessel function, we implemented approx_besselroots.","category":"page"},{"location":"besselroots/","page":"Roots of Bessel function","title":"Roots of Bessel function","text":"approx_besselroots(ν::Real, n::Integer)","category":"page"},{"location":"besselroots/#FastGaussQuadrature.approx_besselroots-Tuple{Real, Integer}","page":"Roots of Bessel function","title":"FastGaussQuadrature.approx_besselroots","text":"approx_besselroots(ν::Real, n::Integer) -> Vector{Float64}\n\nReturn the first n roots of Bessel function. Note that this function is only 12-digits of precision.\n\nJ_nu(x) = sum_m=0^inftyfrac(-1)^jGamma(nu+j+1)j left(fracx2right)^2j+nu\n\nExamples\n\njulia> ν = 0.3;\n\njulia> roots = approx_besselroots(ν, 10);\n\njulia> zeros = (x -> besselj(ν, x)).(roots);\n\njulia> all(zeros .< 1e-12)\ntrue\n\n\n\n\n\n","category":"method"},{"location":"besselroots/","page":"Roots of Bessel function","title":"Roots of Bessel function","text":"This method approx_besselroots is used to calculate gaussjacobi and gausslaguerre.","category":"page"},{"location":"reference/#References","page":"References","title":"References","text":"","category":"section"},{"location":"reference/","page":"References","title":"References","text":"[1] I. Bogaert, \"Iteration-free computation of Gauss-Legendre quadrature nodes and weights\", SIAM J. Sci. Comput., 36(3), A1008-A1026, 2014.","category":"page"},{"location":"reference/","page":"References","title":"References","text":"[2] A. Glaser, X. Liu, and V. Rokhlin. \"A fast algorithm for the calculation of the roots of special functions.\" SIAM J. Sci. Comput., 29 (2007), 1420-1438.","category":"page"},{"location":"reference/","page":"References","title":"References","text":"[3] N. Hale and A. Townsend, \"Fast and accurate computation of Gauss-Legendre and Gauss-Jacobi quadrature nodes and weights\", SIAM J. Sci. Comput., 2012.","category":"page"},{"location":"reference/","page":"References","title":"References","text":"[4] J. C. Mason and D. C. Handscomb, \"Chebyshev Polynomials\", CRC Press, 2002.","category":"page"},{"location":"reference/","page":"References","title":"References","text":"[5] P. Opsomer, (in preparation).","category":"page"},{"location":"reference/","page":"References","title":"References","text":"[6] A. Townsend, The race for high order Gauss-Legendre quadrature, in SIAM News, March 2015.","category":"page"},{"location":"reference/","page":"References","title":"References","text":"[7] A. Townsend, T. Trogdon, and S. Olver, \"Fast computation of Gauss quadrature nodes and weights on the whole real line\", to appear in IMA Numer. Anal., 2014.","category":"page"},{"location":"reference/","page":"References","title":"References","text":"[8] M. Vanlessen, \"Strong asymptotics of Laguerre-Type orthogonal polynomials and applications in Random Matrix Theory\", Constr. Approx., 25:125-175, 2007.","category":"page"},{"location":"benchmark/#Benchmark","page":"Benchmark","title":"Benchmark","text":"","category":"section"},{"location":"benchmark/","page":"Benchmark","title":"Benchmark","text":"Here we compute 100000 nodes and weights of the Gauss rules. Try a million or ten million.","category":"page"},{"location":"benchmark/","page":"Benchmark","title":"Benchmark","text":"using LinearAlgebra, BenchmarkTools, FastGaussQuadrature\n@btime gausschebyshev(100000);\n@btime gausslegendre(100000);\n@btime gaussjacobi(100000, 0.9, -0.1);\n@btime gaussradau(100000);\n@btime gausslobatto(100000);\n@btime gausslaguerre(100000);\n@btime gausshermite(100000);","category":"page"},{"location":"benchmark/","page":"Benchmark","title":"Benchmark","text":"The paper [1] computed a billion Gauss-Legendre nodes. So here we will do a billion + 1.","category":"page"},{"location":"benchmark/","page":"Benchmark","title":"Benchmark","text":"julia> @btime gausslegendre(1000000001);\n  23.363 s (10 allocations: 22.35 GiB)","category":"page"},{"location":"benchmark/","page":"Benchmark","title":"Benchmark","text":"(The nodes near the endpoints coalesce in 16-digits of precision.)","category":"page"},{"location":"gaussquadrature/#Gaussian-Quadrature","page":"Gaussian Quadrature","title":"Gaussian Quadrature","text":"","category":"section"},{"location":"gaussquadrature/#Gauss-Legendre-quadrature","page":"Gaussian Quadrature","title":"Gauss-Legendre quadrature","text":"","category":"section"},{"location":"gaussquadrature/","page":"Gaussian Quadrature","title":"Gaussian Quadrature","text":"Gauss quadrature for the weight function w(x) = 1.","category":"page"},{"location":"gaussquadrature/","page":"Gaussian Quadrature","title":"Gaussian Quadrature","text":"For n  5: Use an analytic expression.\nFor n  60: Use Newton's method to solve P_n(x)=0. Evaluate Legendre polynomials P_n and their derivatives P_n by 3-term recurrence. Weights are related to P_n.\nFor n  60: Use asymptotic expansions for the Legendre nodes and weights [1].","category":"page"},{"location":"gaussquadrature/","page":"Gaussian Quadrature","title":"Gaussian Quadrature","text":"gausslegendre(n::Integer)","category":"page"},{"location":"gaussquadrature/#FastGaussQuadrature.gausslegendre-Tuple{Integer}","page":"Gaussian Quadrature","title":"FastGaussQuadrature.gausslegendre","text":"gausslegendre(n::Integer) -> x, w  # nodes, weights\n\nReturn nodes x and weights w of Gauss-Legendre quadrature.\n\nint_-1^1 f(x) dx approx sum_i=1^n w_i f(x_i)\n\nExamples\n\njulia> x, w = gausslegendre(3);\n\njulia> f(x) = x^4;\n\njulia> I = dot(w, f.(x));\n\njulia> I ≈ 2/5\ntrue\n\n\n\n\n\n","category":"method"},{"location":"gaussquadrature/#Gauss-Hermite-quadrature","page":"Gaussian Quadrature","title":"Gauss-Hermite quadrature","text":"","category":"section"},{"location":"gaussquadrature/","page":"Gaussian Quadrature","title":"Gaussian Quadrature","text":"Gauss quadrature for the weight function w(x) = exp(-x^2) on the real line.","category":"page"},{"location":"gaussquadrature/","page":"Gaussian Quadrature","title":"Gaussian Quadrature","text":"For n  200: Use Newton's method to solve H_n(x)=0. Evaluate Hermite polynomials H_n and their derivatives H_n by three-term recurrence.\nFor n  200: Use Newton's method to solve H_n(x)=0. Evaluate H_n and H_n by a uniform asymptotic expansion, see [7].","category":"page"},{"location":"gaussquadrature/","page":"Gaussian Quadrature","title":"Gaussian Quadrature","text":"The paper [7] also derives an O(n) algorithm for generalized Gauss-Hermite nodes and weights associated to weight functions of the form exp(-V(x)), where V(x) is a real polynomial.","category":"page"},{"location":"gaussquadrature/","page":"Gaussian Quadrature","title":"Gaussian Quadrature","text":"gausshermite(n::Integer)","category":"page"},{"location":"gaussquadrature/#FastGaussQuadrature.gausshermite-Tuple{Integer}","page":"Gaussian Quadrature","title":"FastGaussQuadrature.gausshermite","text":"gausshermite(n::Integer) -> x, w  # nodes, weights\n\nReturn nodes x and weights w of Gauss-Hermite quadrature.\n\nint_-infty^+infty f(x) exp(-x^2) dx approx sum_i=1^n w_i f(x_i)\n\nExamples\n\njulia> x, w = gausshermite(3);\n\njulia> f(x) = x^4;\n\njulia> I = dot(w, f.(x));\n\njulia> I ≈ 3(√π)/4\ntrue\n\n\n\n\n\n","category":"method"},{"location":"gaussquadrature/#Gauss-Laguerre-quadrature","page":"Gaussian Quadrature","title":"Gauss-Laguerre quadrature","text":"","category":"section"},{"location":"gaussquadrature/","page":"Gaussian Quadrature","title":"Gaussian Quadrature","text":"Gauss quadrature for the weight function w(x) = exp(-x) on 0+infty)","category":"page"},{"location":"gaussquadrature/","page":"Gaussian Quadrature","title":"Gaussian Quadrature","text":"For n  128: Use the Golub-Welsch algorithm.\nFor method=GLR: Use the Glaser-Lui-Rohklin algorithm. Evaluate Laguerre polynomials L_n and their derivatives L_n by using Taylor series expansions near roots generated by solving the second-order differential equation that L_n satisfies, see [2].\nFor n  128: Use a Newton procedure on Riemann-Hilbert asymptotics of Laguerre polynomials, see [5], based on [8]. There are some heuristics to decide which expression to use, it allows a general weight w(x) = x^alpha exp(-q_m x^m) and this is O(sqrtn) when allowed to stop when the weights are below the smallest positive floating point number.","category":"page"},{"location":"gaussquadrature/","page":"Gaussian Quadrature","title":"Gaussian Quadrature","text":"gausslaguerre(n::Integer)","category":"page"},{"location":"gaussquadrature/#FastGaussQuadrature.gausslaguerre-Tuple{Integer}","page":"Gaussian Quadrature","title":"FastGaussQuadrature.gausslaguerre","text":"gausslaguerre(n::Integer) -> x, w  # nodes, weights\n\nReturn nodes x and weights w of Gauss-Laguerre quadrature.\n\nint_0^+infty f(x) e^-x dx approx sum_i=1^n w_i f(x_i)\n\nExamples\n\njulia> x, w = gausslaguerre(3);\n\njulia> f(x) = x^4;\n\njulia> I = dot(w, f.(x));\n\njulia> I ≈ 24\ntrue\n\n\n\n\n\n","category":"method"},{"location":"gaussquadrature/","page":"Gaussian Quadrature","title":"Gaussian Quadrature","text":"gausslaguerre(n::Integer, α::Real)","category":"page"},{"location":"gaussquadrature/#FastGaussQuadrature.gausslaguerre-Tuple{Integer, Real}","page":"Gaussian Quadrature","title":"FastGaussQuadrature.gausslaguerre","text":"gausslaguerre(n::Integer, α::Real) -> x, w  # nodes, weights\n\nReturn nodes x and weights w of generalized Gauss-Laguerre quadrature.\n\nint_0^+infty f(x) x^alpha e^-x dx approx sum_i=1^n w_i f(x_i)\n\nExamples\n\njulia> x, w = gausslaguerre(3, 1.0);\n\njulia> f(x) = x^4;\n\njulia> I = dot(w, f.(x));\n\njulia> I ≈ 120\ntrue\n\nOptionally, a reduced quadrature rule can be computed. In that case, only those points and weights are computed for which the weight does not underflow in the floating point precision type. Supply the optional argument reduced = true.\n\nThough the code is generic, heuristical choices on the choice of the algorithm are based on achieving machine precision accuracy only for Float64 type. In case the default choice of algorithm does not achieve the desired accuracy, the user can manually invoke the following routines:\n\ngausslaguerre_GW: computation based on Golub-Welsch\ngausslaguerre_rec: computation based on Newton iterations applied to evaluation  using the recurrence relation\ngausslaguerre_asy: the asymptotic expansions\n\n\n\n\n\n","category":"method"},{"location":"gaussquadrature/#Gauss-Chebyshev-quadrature","page":"Gaussian Quadrature","title":"Gauss-Chebyshev quadrature","text":"","category":"section"},{"location":"gaussquadrature/","page":"Gaussian Quadrature","title":"Gaussian Quadrature","text":"There are four kinds of Gauss-Chebyshev quadrature rules, corresponding to four weight functions:","category":"page"},{"location":"gaussquadrature/","page":"Gaussian Quadrature","title":"Gaussian Quadrature","text":"1st kind, weight function w(x) = 1sqrt1-x^2\n2nd kind, weight function w(x) = sqrt1-x^2\n3rd kind, weight function w(x) = sqrt(1+x)(1-x)\n4th kind, weight function w(x) = sqrt(1-x)(1+x)","category":"page"},{"location":"gaussquadrature/","page":"Gaussian Quadrature","title":"Gaussian Quadrature","text":"They are all have explicit simple formulas for the nodes and weights [4].","category":"page"},{"location":"gaussquadrature/","page":"Gaussian Quadrature","title":"Gaussian Quadrature","text":"gausschebyshev(n::Integer, kind::Integer)","category":"page"},{"location":"gaussquadrature/#FastGaussQuadrature.gausschebyshev-Tuple{Integer, Integer}","page":"Gaussian Quadrature","title":"FastGaussQuadrature.gausschebyshev","text":"gausschebyshev(n::Integer) -> x, w  # nodes, weights\ngausschebyshev(n::Integer, 1) -> x, w  # nodes, weights\n\nReturn nodes x and weights w of Gauss-Chebyshev quadrature of the 1st kind.\n\nint_-1^1 fracf(x)sqrt1-x^2 dx approx sum_i=1^n w_i f(x_i)\n\nExamples\n\njulia> x, w = gausschebyshev(3);\n\njulia> f(x) = x^4;\n\njulia> I = dot(w, f.(x));\n\njulia> I ≈ 3π/8\ntrue\n\n\n\ngausschebyshev(n::Integer, 2) -> x, w  # nodes, weights\n\nReturn nodes x and weights w of Gauss-Chebyshev quadrature of the 2nd kind.\n\nint_-1^1 f(x)sqrt1-x^2 dx approx sum_i=1^n w_i f(x_i)\n\nExamples\n\njulia> x, w = gausschebyshev(3, 2);\n\njulia> f(x) = x^4;\n\njulia> I = dot(w, f.(x));\n\njulia> I ≈ π/16\ntrue\n\n\n\ngausschebyshev(n::Integer, 3) -> x, w  # nodes, weights\n\nReturn nodes x and weights w of Gauss-Chebyshev quadrature of the 3rd kind.\n\nint_-1^1 f(x)sqrtfrac1+x1-x dx approx sum_i=1^n w_i f(x_i)\n\nExamples\n\njulia> x, w = gausschebyshev(3, 3);\n\njulia> f(x) = x^4;\n\njulia> I = dot(w, f.(x));\n\njulia> I ≈ 3π/8\ntrue\n\n\n\ngausschebyshev(n::Integer, 4) -> x, w  # nodes, weights\n\nReturn nodes x and weights w of Gauss-Chebyshev quadrature of the 4th kind.\n\nint_-1^1 f(x)sqrtfrac1-x1+x dx approx sum_i=1^n w_i f(x_i)\n\nExamples\n\njulia> x, w = gausschebyshev(3, 4);\n\njulia> f(x) = x^4;\n\njulia> I = dot(w, f.(x));\n\njulia> I ≈ 3π/8\ntrue\n\n\n\n\n\n","category":"method"},{"location":"gaussquadrature/#Gauss-Jacobi-quadrature","page":"Gaussian Quadrature","title":"Gauss-Jacobi quadrature","text":"","category":"section"},{"location":"gaussquadrature/","page":"Gaussian Quadrature","title":"Gaussian Quadrature","text":"Gauss quadrature for the weight functions w(x) = (1-x)^alpha(1+x)^beta, alphabeta  -1.","category":"page"},{"location":"gaussquadrature/","page":"Gaussian Quadrature","title":"Gaussian Quadrature","text":"For n  100: Use Newton's method to solve P_n(x)=0. Evaluate P_n and P_n by three-term recurrence.\nFor n  100: Use Newton's method to solve P_n(x)=0. Evaluate P_n and P_n by an asymptotic expansion (in the interior of -11) and the three-term recurrence O(n^-2) close to the endpoints. (This is a small modification to the algorithm described in [3].)\nFor max(alphabeta)  5: Use the Golub-Welsch algorithm requiring O(n^2) operations.","category":"page"},{"location":"gaussquadrature/","page":"Gaussian Quadrature","title":"Gaussian Quadrature","text":"gaussjacobi(n::Integer, α::Real, β::Real)","category":"page"},{"location":"gaussquadrature/#FastGaussQuadrature.gaussjacobi-Tuple{Integer, Real, Real}","page":"Gaussian Quadrature","title":"FastGaussQuadrature.gaussjacobi","text":"gaussjacobi(n::Integer, α::Real, β::Real) -> x, w  # nodes, weights\n\nReturn nodes x and weights w of Gauss-Jacobi quadrature for exponents α and β.\n\nint_-1^1 f(x) (1-x)^alpha(1+x)^beta dx approx sum_i=1^n w_i f(x_i)\n\nExamples\n\njulia> x, w = gaussjacobi(3, 1/3, -1/3);\n\njulia> f(x) = x^4;\n\njulia> I = dot(w, f.(x));\n\njulia> I ≈ 268π/729(√3)\ntrue\n\n\n\n\n\n","category":"method"},{"location":"gaussquadrature/#Gauss-Radau-quadrature","page":"Gaussian Quadrature","title":"Gauss-Radau quadrature","text":"","category":"section"},{"location":"gaussquadrature/","page":"Gaussian Quadrature","title":"Gaussian Quadrature","text":"Gauss quadrature for the weight function w(x)=1, except the endpoint -1 is included as a quadrature node.","category":"page"},{"location":"gaussquadrature/","page":"Gaussian Quadrature","title":"Gaussian Quadrature","text":"The Gauss-Radau nodes and weights can be computed via the (01) Gauss-Jacobi nodes and weights[3].","category":"page"},{"location":"gaussquadrature/","page":"Gaussian Quadrature","title":"Gaussian Quadrature","text":"gaussradau(n::Integer)","category":"page"},{"location":"gaussquadrature/#FastGaussQuadrature.gaussradau-Tuple{Integer}","page":"Gaussian Quadrature","title":"FastGaussQuadrature.gaussradau","text":"gaussradau(n::Integer) -> x, w  # nodes, weights\n\nReturn nodes x and weights w of Gauss-Radau quadrature.\n\nint_-1^1 f(x) dx approx sum_i=1^n w_i f(x_i)\n\nExamples\n\njulia> x, w = gaussradau(3);\n\njulia> f(x) = x^4;\n\njulia> I = dot(w, f.(x));\n\njulia> I ≈ 2/5\ntrue\n\nNote that the first node is fixed at -1.\n\njulia> x, w = gaussradau(3);\n\njulia> x[1]\n-1.0\n\n\n\n\n\n","category":"method"},{"location":"gaussquadrature/#Gauss-Lobatto-quadrature","page":"Gaussian Quadrature","title":"Gauss-Lobatto quadrature","text":"","category":"section"},{"location":"gaussquadrature/","page":"Gaussian Quadrature","title":"Gaussian Quadrature","text":"Gauss quadrature for the weight function w(x)=1, except the endpoints -1 and 1 are included as nodes.","category":"page"},{"location":"gaussquadrature/","page":"Gaussian Quadrature","title":"Gaussian Quadrature","text":"The Gauss-Lobatto nodes and weights can be computed via the (11) Gauss-Jacobi nodes and weights[3].","category":"page"},{"location":"gaussquadrature/","page":"Gaussian Quadrature","title":"Gaussian Quadrature","text":"gausslobatto(n::Integer)","category":"page"},{"location":"gaussquadrature/#FastGaussQuadrature.gausslobatto-Tuple{Integer}","page":"Gaussian Quadrature","title":"FastGaussQuadrature.gausslobatto","text":"gausslobatto(n::Integer) -> x, w  # nodes, weights\n\nReturn nodes x and weights w of Gauss-Lobatto quadrature.\n\nint_-1^1 f(x) dx approx sum_i=1^n w_i f(x_i)\n\nExamples\n\njulia> x, w = gausslobatto(4);\n\njulia> f(x) = x^4;\n\njulia> I = dot(w, f.(x));\n\njulia> I ≈ 2/5\ntrue\n\nNote that the both ends of nodes are fixed at -1 and 1.\n\njulia> x, w = gausslobatto(4);\n\njulia> x[1], x[end]\n(-1.0, 1.0)\n\n\n\n\n\n","category":"method"},{"location":"#FastGaussQuadrature.jl","page":"Home","title":"FastGaussQuadrature.jl","text":"","category":"section"},{"location":"#Abstract","page":"Home","title":"Abstract","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"FastGaussQuadrature.jl is a Julia package to compute n-point Gauss quadrature nodes and weights to 16-digit accuracy and in O(n) time. So far the package includes gausschebyshev(), gausslegendre(), gaussjacobi(), gaussradau(), gausslobatto(), gausslaguerre(), and gausshermite(). This package is heavily influenced by Chebfun.","category":"page"},{"location":"","page":"Home","title":"Home","text":"An introduction to Gauss quadrature can be found here. For a quirky account on the history of computing Gauss-Legendre quadrature, see [6].","category":"page"},{"location":"#Our-Aims","page":"Home","title":"Our Aims","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The fastest Julia code for Gauss quadrature nodes and weights (without tabulation).\nChange the perception that Gauss quadrature rules are expensive to compute.","category":"page"},{"location":"#First-example","page":"Home","title":"First example","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"To check an integral","category":"page"},{"location":"","page":"Home","title":"Home","text":"int_-1^1 x^4 dx = frac25","category":"page"},{"location":"","page":"Home","title":"Home","text":"by numerically, try following code.","category":"page"},{"location":"","page":"Home","title":"Home","text":"using FastGaussQuadrature, LinearAlgebra\nx, w = gausslegendre(3)\nf(x) = x^4\nI = dot(w, f.(x))\nI ≈ 2/5","category":"page"}]
}
