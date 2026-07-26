# Hévéa numerical-kernel map

| Upstream flat-torus location | Responsibility | Sphere-neutral replacement |
|---|---|---|
| `utils/include/UTI_utils.h::J0_inv` | inverse first monotone branch of J0 | `src/hevea_numeric.hpp::j0Inverse` |
| `iso/.../ISO_integrate.h::flow` | integrate characteristic scalar coordinate | `integrateDp54` with explicit cylinder field |
| `ISO_h_s` + `integrate_h_hai` | integrate corrugation vector along each flow line | generic `State<N>` adaptive integrator |
| `CYL_embedding::gluing` | periodic cylinder seam correction | sphere layer must instead use paper `chi` cap/ribbon gluing |
| `CYL_embedding::cyl_to_torus` | inverse flow and doubly-periodic reparametrization | Phase 4 cylinder inverse-flow, x-periodic/y-bounded |

The original Hairer wrapper requests `rtol=atol=1e-10`.  The extracted layer
defaults to Dormand-Prince 5(4), `rtol=1e-11`, `atol=1e-12`, reports accepted
and rejected steps, and is verified against analytic solutions.  It deliberately
does not include any `TOR_embedding` or y-periodic indexing.
