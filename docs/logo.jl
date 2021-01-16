using FastGaussQuadrature, LinearAlgebra, Luxor
logocolors = Luxor.Colors.JULIA_LOGO_COLORS

@svg begin
    n = 15
    Drawing(500,500, "docs/src/assets/logo.svg")
    origin()
    unitlength = 210
    translate(0, unitlength)

    Y = (n+1/2)/π

    x_lgd, w_lgd = gausslegendre(n)
    x_cbs, w_cbs = gausschebyshev(n,3)
    x_jcb, w_jcb = gaussjacobi(n,5/2,1/2)

    c_lgd = logocolors[1]  # Red
    c_cbs = logocolors[2]  # Green
    c_jcb = logocolors[4]  # Purple

    p_jcb = [Point(unitlength*x_jcb[i], -Y*unitlength*w_jcb[i]) for i in 1:n]
    p_cbs = [Point(unitlength*x_cbs[i], -Y*unitlength*w_cbs[i]) for i in 1:n]
    p_lgd = [Point(unitlength*x_lgd[i], -Y*unitlength*w_lgd[i]) for i in 1:n]

    l_jcb = [norm(p_jcb[i+1]-p_jcb[i]) for i in 1:n-1]
    l_cbs = [norm(p_cbs[i+1]-p_cbs[i]) for i in 1:n-1]
    l_lgd = [norm(p_lgd[i+1]-p_lgd[i]) for i in 1:n-1]

    r_jcb = zeros(n)
    r_cbs = zeros(n)
    r_lgd = zeros(n)

    r_jcb[end] = norm(l_jcb[end])/2 * 0.8
    r_cbs[end] = norm(l_cbs[end])/2 * 0.7
    r_lgd[end] = norm(l_lgd[end])/2 * 0.75

    for i in reverse(1:n-1) r_jcb[i] = l_jcb[i] - r_jcb[i+1] end
    for i in reverse(1:n-1) r_cbs[i] = l_cbs[i] - r_cbs[i+1] end
    for i in reverse(1:n-1) r_lgd[i] = l_lgd[i] - r_lgd[i+1] end

    sethue(c_cbs)
    for i in 5:n circle(p_cbs[i], r_cbs[i], :fill) end
    sethue(c_jcb)
    for i in 1:n circle(p_jcb[i], r_jcb[i], :fill) end
    sethue(c_lgd)
    for i in 1:n circle(p_lgd[i], r_lgd[i], :fill) end
    sethue(c_cbs)
    for i in 1:4 circle(p_cbs[i], r_cbs[i], :fill) end

    finish()
    preview()
end

run(`convert -density 256x256 -background transparent docs/src/assets/logo.svg -define icon:auto-resize -colors 256 docs/src/assets/favicon.ico`)
