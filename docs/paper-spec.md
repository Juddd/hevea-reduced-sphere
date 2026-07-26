# Reduced-sphere paper specification

Source: Bartzos, Borrelli, Denis, Lazarus, Rohmer, Thibert, *An Explicit
Isometric Reduction of the Unit Sphere into an Arbitrarily Small Ball*, FoCM
(2017), DOI `10.1007/s10208-017-9360-1`.  Local immutable source hash is listed
in `tests/fixtures/SHA256SUMS`.

## Domain and metrics

- `D = S^1 x (-yInfinity,yInfinity)`, with `x` periodic modulo `2 pi`.
- Round parametrization `h(x,y)=(cos(y)cos(x),cos(y)sin(x),sin(y))`.
- Equation (2.1): `I_h = [[cos^2(y),0],[0,1]]`.
- Equation (2.2), for `B=[[E,F],[F,G]]`:
  `rho1(B)=E+F`, `rho2(B)=E-F`, `rho3(B)=G-E`, using
  `ell1=(dx+dy)/sqrt(2)`, `ell2=(-dx+dy)/sqrt(2)`, `ell3=dy`.

## Initial short embedding

- The rotational ribbon is `f0(x,y)=(X(y)cos(x),X(y)sin(x),Z(y))`.
- Conditions (C1)-(C4): immersion/embedding, image inside the requested ball,
  C1 attachment to the translated polar caps, and strict positivity of every
  primitive coordinate needed by the iteration.
- Equations (3.5)-(3.7):
  `cos^2(y)>X(y)^2`, `sin^2(y)+X(y)^2>X'(y)^2+Z'(y)^2`, and the endpoint
  position/derivative data
  `X(+-yInf)=cos(yInf)`, `X'(-yInf)=sin(yInf)`,
  `X'(+yInf)=-sin(yInf)`, `Z(-yInf)=eta-sin(yInf)`,
  `Z(+yInf)=sin(yInf)-eta`, `Z'(+-yInf)=cos(yInf)`.
- The implementation described by the authors uses a Hermite cubic spline
  reparametrized by a second Hermite cubic spline; its coordinates have degree
  nine.  The coefficients are not printed in the paper.

## One-dimensional convex integration

- Equation (4.8): integrate
  `r (cos(theta) t + sin(theta) n)` along each image flow curve.
- Equation (4.9): `alpha=J0Inverse(norm(df(w))/r)`, on the first monotone branch
  of `J0` ending at its first positive root.
- Equation (4.10): use the shifted flow parameter `s`, not an arbitrary grid
  coordinate; `theta=alpha cos(2 pi N s)`.
- Equations (4.11)-(4.12): the averaged tangent equals the old derivative and
  the positional perturbation is `O(1/N)`.
- Lemmas 2-3: affine lower phase and the arithmetic closure condition are
  required both for quotient periodicity and transverse derivative control.

## Primitive directions and characteristic fields

- Equation (5.13): `r^2=norm(df(w_i))^2+rho_i`.
- Equation (5.14): `mu_i=I_f+rho_i(I_target-I_f) ell_i tensor ell_i`.
- Equation (5.15): kernel directions
  `v1=(-1,1)/sqrt(2)`, `v2=-(1,1)/sqrt(2)`, `v3=(-1,0)`.
- `w_i` is uniquely fixed by `ell_i(w_i)=1` and
  `mu_i(w_i,v_i)=0`.  Since `ell_i(v_i)=0`, the latter equals orthogonality in
  the current induced metric along `v_i`.
- Equation (5.16): the actual flow parameter is
  `s=ell_i(gamma_p(s)-(0,0))`.  Ridge lines are parallel to `v_i`.

## Nested ribbons, target metrics, and gluing

- Equation (6.18): `I_f0 <= I_k,i <= I_l,j <= I_h` in lexicographic order.
- Figure 6: `lambda_k,i` equals one through the previous ribbon boundary,
  transitions to zero by the midpoint, and is zero on the later gluing half.
  Conversely, `chi_k,i` is zero through the midpoint and transitions to one at
  the current ribbon boundary.  Their transition interiors do not overlap.
- Equation (6.19): build `F_k,i` by high-accuracy integration on complete flow
  lines from the south boundary of `D_k,i`.
- Lemma 7: `2 pi N_k,i ell_i(e_x)` is integral.
- Equation (6.20): `f_k,i=(1-chi)F_k,i+chi f_k,i-1`; outside `D_k,i`, keep `f0`.
- Equation (7.21): each step corrects exactly one primitive coordinate relative
  to its altered target metric; a negative coordinate is a hard precondition
  failure, never something to clamp.

## Paper numerical fixture

Grid: `4000 x 20000`.  Corrugation numbers `4.72,31.96,334.92`; visible ridge
counts `21,142,997`.  The machine-readable table is
`tests/fixtures/paper_metrics.json`.  All averages and suprema must be computed
on the paper's stated domain/target for that row, not on an invented global
metric.

## Code responsibility map

| Paper responsibility | Required module/test |
|---|---|
| (2.1)-(2.2) metrics/primitives | independent `reference_metrics`, production metric module |
| (3.5)-(3.7), C1-C4 | profile generator and dense certificate |
| (4.8)-(4.12) Bessel corrugation | shared Hevea numerical kernel |
| Lemmas 2-5, (5.13)-(5.16) | flow/phase closure/inverse-flow tests |
| (6.18)-(6.20), Figure 6 | target schedule, disjoint lambda/chi tests, cap seam tests |
| (7.21), Stage theorem prerequisites | per-step positive primitive/Jacobian hard gates |
| Figure 9/table | paper-grid verifier, mesh audit, fixed-camera WL render |
