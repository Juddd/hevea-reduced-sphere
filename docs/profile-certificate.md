# Certified initial reduced-sphere profile

The production profile is the symmetric degree-thirteen Hermite ansatz in
`src/reduced_sphere_profile.hpp`.  Endpoint values and derivatives fix the
constant and linear-in-`t=1-(y/yInfinity)^2` terms.  The ten remaining
coefficients are

```text
-0.52164358694220492, -0.024827807054374828,
 1.0452244736858278,  -1.3656657598369506,
 0.54300716429149642, -0.38685452062281694,
 0.31457715720002261, -0.61976108254554241,
 1.1664678885754221,  -0.73680735226233729
```

They were selected by the deterministic C++ envelope search in
`tools/stage1_envelope_model.cpp`.  The search jointly constrains the author's
public WRL mode-21 latitude envelope, its axisymmetric mean profile, the
paper's initial metric row, C1--C4, all three static primitive/characteristic
margins and the containing-ball radius.  The associated first-stage
transition parameters are

```text
previous = 0.63199298223094058
outer    = 1.2792974481304009
power    = 0.60979123180328654
split    = 0.80193065649632189
fraction = 0.237
```

For the machine-readable author fixture
`tests/fixtures/reference_stage1_mode21_envelope.csv`, the independent model
has envelope RMS/max errors `0.0012763432/0.0023301835`.  Its initial profile
has metric mean/max `0.8500041/1.1930015`, radius `0.5035085`, WRL
axisymmetric-profile C0 RMS/max `0.0176735/0.0209893`, static primitive minima
`0.1098067,0.00749035,0.0330070`, and zero search penalty.  As a concrete
high-latitude check at `|y|=1.10`, the author WRL mode-21 amplitude is
`0.00212205`; the actual `4000x4000` full-flow mesh gives `0.00205424`.  This
latitude-envelope constraint is what prevents corrugation energy from being
concentrated at the middle while disappearing too early near the caps.

The production ribbon fractions are
`0.41148336662952684,0.83293586428607069,0.89850213833742598,0.96908999242130689`.
The first lambda/chi transition uses the power and split above; the other two
use the default smooth step.  Environment-variable overrides remain available
for diagnostics, but these values are the no-environment production defaults.

Two actual-flow probes distinguish structural stability from finite-frequency
metric convergence:

- At `4000x4000`, all three stages complete with final primitive/Jacobian
  minima `0.0051292/0.690848`, radius `0.504109`, and eight fresh VTK outputs.
  Its Stage-3 target mean `0.0635959` is deliberately not used as a paper-grid
  certificate: four latitude samples per 997-cycle corrugation under-resolve
  the final derivative.
- At `4000x8000`, Stage-3 round mean/max are `0.661384/1.02816`, target
  mean/max are `0.0382096/0.207837`, radius is `0.504062`, and final
  primitive/Jacobian minima are `0.00394808/0.690908`.  This convergence probe
  passes the Phase-5 numerical bands and demonstrates that the elevated
  `4000x4000` target mean is a sampling artifact rather than a geometry
  failure.

The fresh `4000x20000` production certificate uses executable SHA256
`f53b78395bb4159abc020336d351ab46a52ae205b5cc1ffc7dd042c234d31731`
and config hash
`0b76d983209473c0b8fc92d53506bb8d928c72d7fb055afc58b752638e086360`.
Its independent stage rows `(round mean, round max; target mean, target max)`
are

```text
Stage 1  (0.799190, 1.081687; 0.106374, 0.213958)
Stage 2  (0.737414, 1.019307; 0.089828, 0.184377)
Stage 3  (0.658453, 1.022083; 0.035475, 0.218496)
```

The final primitive/Jacobian minima are `0.00315251/0.69072`, the containing
radius is `0.504041`, and the independent/native metric deltas are at most
`3.29e-6`.  The mesh audit reports a closed consistently oriented sphere with
Euler characteristic 2, no non-finite or degenerate cells; all binary/log/
manifest SHA256 checks pass.

The third characteristic family uses the equator as its affine integration
transversal and integrates toward both ribbon boundaries.  This changes only
the characteristic-wise integration constant; the prescribed derivative,
phase, target metric and 997 visible ridges are unchanged.

The stable identifier is
`profile-20260725-author-envelope-paper-flow-certified-v5`.
`tools/certify_profile --dense 1000000` is the authoritative dense initial
certificate; `tests/run_fast_e2e.sh` checks actual-flow structural stability,
and the paper-grid verifier certifies the finite corrugations.  Defects tend to
zero at `+-yInfinity`, necessarily, because the ribbon is C1-attached to the
isometric translated caps; no uniform positive lower bound is claimed on that
open-domain limit.

The paper does not publish its spline coefficients or transition parameters.
The public WRL constrains the final surface but does not uniquely identify
those hidden inputs.  This is therefore a constrained, numerically certified
reproduction within the paper's printed tolerances, not a claim that the
authors' unpublished coefficient vector has been recovered exactly.
