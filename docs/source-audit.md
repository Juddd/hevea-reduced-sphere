# Source and supplementary-material audit

Audit date: 2026-07-23.  The goal was to find the actual reduced-sphere program
or the omitted degree-nine profile and ribbon parameters before attempting a
reimplementation.

| Resource/query | Result |
|---|---|
| `https://hevea-project.fr/indexSphere.html` | HTTP 200; project landing page, no source link |
| `pageSpherePublication.html` | HTTP 200; paper/publication material only |
| `pageSphereImages.html` | HTTP 200; rendered images only |
| `codeSphere/` | HTTP 404 |
| `codeSphere/hevea-sphere.tgz` | HTTP 404 |
| `hevea-sphere.tgz` | HTTP 404 |
| HAL `hal-01374730` | HTTP 200; publication record/PDF, no code attachment found |
| GitHub repository search `hevea sphere convex integration` | zero repositories |
| GitHub code search exact paper title | four bibliographic/index matches, no implementation |
| GitHub code search author names and sphere | four bibliographic/index matches, no implementation |
| GitHub code search paper numbers | no relevant implementation match |
| Local `/home/yode` archive/source search | no sphere tgz/zip/cpp; only flat-torus Hevea source |
| Tavily search skill | unavailable without OAuth token; fallback sources above used |

The public Francis Lazarus GitHub pages link only to DOI and
`hevea-project.fr/pdfSphere/focm-revised-2.pdf`.  The repository search did not
expose a separate code repository.  This is evidence of non-discovery, not a
claim that no private or unindexed source exists.

## Licensing consequence

The vendored flat-torus Hevea source is GPL-3.0-or-later and may be reused under
that license with attribution.  No sphere source was found, so no unknown code
is copied.  Phase 2 may extract the generic numerical mechanisms from the GPL
Hevea tree; Phase 3 must either discover authoritative sphere coefficients or
produce and label a reproducible constrained reconstruction.
