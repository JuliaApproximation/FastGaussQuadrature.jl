# Finite native Laguerre GST engineering model

Frozen v1, 2026-11-16, before the broad candidate/reference grid. This
model covers the measured native GST routes after removal of the transient
two-part square correction. It is not a theorem for all orders/parameters and
does not certify unchanged REC/ASY routes. Coefficients are operational reserves
derived from the formulas below, not fitted to observed candidate errors.

Write `e=eps(T)`, `a=alpha`, `delta=a+1`, `nu=4n+2a+2`, `x=z^2` at an independent
root, and `mu=Gamma(delta)`. Inputs are exactly the original represented `a`.
All model quantities use reference precision. Define

```text
M(s) = 8*(2 + abs(loggamma(s)) + s*abs(log(s)))
C(x) = 16x + 16n*sqrt(x/nu)
Bx(x) = C(x) + x
K(x) = abs((a+1/2)/x - 1)
G(x) = abs(a+1/2)*abs(log(x/pivot)) + abs(x-pivot)
L(x,w) = 16(n+1) + K(x)*C(x) + 8*(1+G(x)+abs(log(w/mu))).
```

Here `pivot` is the exact mathematical startup x-coordinate used by GST: if
`abs(a)>1/2`, `sqrt(a^2-1/4)`; otherwise the published lower-bound expression in
`gausslaguerre_gst`, evaluated independently in reference arithmetic. It only
sets a common logarithmic scale. The formula is documented in the source from
GST (2019), equation (30); no production intermediate is read by the model.

`C` comes from a z allowance `8e*(z+n/sqrt(nu))`: coordinate-scale rounding plus
accumulation of ordinary phase errors over at most n roots. Multiplication by
`2z` gives the internal x allowance; final `z*z` rounding adds `e*x` for nodes.
This deliberately conservative propagation assumption is challenged by the
grid; it is not a rigorous stability proof near every turning point.

The normal-form equation has `y''=0` at a root, so the exact local derivative
of log-weight with respect to x is `(a+1/2)/x-1`. The derivative amplitude's
linear accumulation reserve is `8(n+1)e`, doubled by the inverse square.
`G` is the magnitude of the direct native log/exponential expression. The
log-probability term reserves scale-dependent exp/log evaluation. `M(s)`
includes a native gamma allowance and the effect of rounding the shifted
argument: `s*abs(digamma(s)) <= 1+s*abs(log(s))` is sufficient here.

## Ordinary normalization: delta > 1/8

Rigorous analytic root bounds avoid using unobserved production roots:

```text
l = delta/n
u = 2n+a-1 + 2sqrt((n-1)*(n-1+a))
Kglobal = 16*(abs(a+1/2)+u)
          + 16n/sqrt(nu)*(abs(a+1/2)/sqrt(l)+sqrt(u))
Gglobal = abs(a+1/2)*max(abs(log(l/pivot)),abs(log(u/pivot)))
          + max(abs(l-pivot),abs(u-pivot))
Mean = 16(n+1)+Kglobal+8*(1+Gglobal+log(n))
gamma_sum = (n-1)*e/(2*(1-(n-1)*e/2))
tau = reduced ? min(1, n*10floatmin(T)/mu) : 0
```

`l` follows from `sum(1/x_i)=n/delta`; `u` follows from a row-sum bound on the
classical Jacobi matrix. `Kglobal` bounds `K*C` and `Gglobal` bounds the direct
log-expression scale on that interval. The entropy inequality
`sum(p_i*abs(log(p_i)))<=log(n)` bounds the weighted normalization mean.
`tau` conservatively accounts for the omitted mass in reduced calculation;
the existing stopping proof bounds each omitted weight below the cutoff.

Require absolute node error at most `e*Bx` and normal-weight relative error at
most

```text
e*(M(delta)+L(x,w)+Mean+16) + gamma_sum + tau/(1-tau).
```

The fixed 16e includes ordinary summation/exp, normalization and final output
rounding reserves. For public Float64 half-integer n>=128 rules, the first
2/12 Bessel-overwritten points are separately preserved controls; this GST
formula is applied only to the remaining points.

## Endpoint route: delta <= 1/8

This route anchors the remaining weights at the second root, rather than
normalizing their sum. It needs a separate model. Let `t=n*x1/delta` and

```text
F(t) = sum(c_k, k=0:n),
c_0=1, c_1=-t,
c_k=c_(k-1)*((k-1-n)/k)*(delta*t/n)/(delta+k-1).
AF=sum(abs(c_k)); AD=sum(k*abs(c_k)/t).
EF=8(n+1)*AF; ED=8(n+1)*AD;
Ct=EF/abs(F'(t))+8abs(t).
logp=sum(log1p(delta/k), k=1:n-1).
```

The sums and derivative identities are reference quantities. A decreasing
alternating tail may be bounded rather than evaluated to degree n. The envelope
still allows n rounded term operations; near the boundary the actual series
uses many fewer. For the production Float64 beta difference, use
`LP=8*(1+abs(logbeta(n,delta))+abs(log(delta))+abs(logp))`.
For other native formats, use `LP=8*(n+abs(logp))`, allowing rounding through
the existing sum of log1p terms.

The first root has absolute allowance
`e*(delta*Ct/n+4x1)`. Its relative-weight allowance is

```text
e*(M(delta)+LP+2ED/abs(F'(t))
   +Ct*abs(1/t+2F''(t)/F'(t))+8).
```

For the second-root anchor define independently

```text
P(x2) = 1F1(1-n;delta+1;x2),
q_0=1, q_k=q_(k-1)*((k-n)/k)*x2/(delta+k),
chi2=sum(abs(q_k))/abs(P(x2)),
Anchor=M(delta+1)+LP+16(n+1)*chi2
       +abs(1/x2+2P'(x2)/P(x2))*(C(x2)+x2)+8.
```

The condition `chi2` explicitly measures cancellation in the finite anchor
series. Its node sensitivity includes the rounded x2 used by production.
For each ordinary root `i>=2`, use node allowance `e*Bx(x_i)` and relative
weight allowance `e*(Anchor+L(x_i,w_i)+L(x2,w2)+8)`. This follows the actual
`anchor*exp(lograw_i-lograw_2)` route; no fictitious normalization mean is used.

## Acceptance and limitations

All relative allowances must be finite and below 1/8. A nonpositive summation
denominator, tau>=1, unavailable finite mass, unresolved series condition, or
failed independent reference makes the model inconclusive, not passed. A lost
required normal weight, wrong indexed root, or invalid full/reduced semantics
is a failure independently of any budget. For subnormal/true-underflow weights,
use absolute allowance `abs(wref)*relative_allowance + one local target spacing`,
and report local-spacing errors and their range limitations separately.

Candidate iteration/truncation must fit this model without changing existing
tighter regression tests. Existing REC, alpha=0 ASY, explicit n<=2, and hybrid
Bessel points remain exact-before controls. Before GW rows are historical
comparison context, not native-GST certification. Original measured errors and
all representation failures remain saved. No below-24-bit calculation is judged
directly with this model: the Float32 result and final conversion are separate.
