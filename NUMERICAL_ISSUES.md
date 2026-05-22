# Numerical Issues: CES Effective Distance

Last updated: 2026-05-21

## Negative-theta log-sum-exp

`ces_effective_distance()` computes the CES power mean in log space for
`theta < 0`:

```text
log(sum_xy w_xy d_xy^theta) = logsumexp(log(w_xy) + theta log(d_xy))
```

For the empirically relevant range `theta in [-7, -3]` and floored pairwise
distances `d in [0.5, 20000]` km, the log-distance contribution spans about
`[-69.3, 4.85]` at `theta = -7`. That range is well inside double-precision
limits. The implementation now includes `log(w_xy)` inside the log-sum-exp,
which also avoids avoidable underflow/overweighting when the closest pair has
very small probability mass.

## Distance floor and harmonic-dominated means

The function applies a default raster-resolution floor:

```text
d_floored = max(d_raw, 0.5 km)
```

This is not numerically neutral when raw distances fall below 0.5 km. For
non-overlapping supports whose pairwise distances are all at least 0.5 km,
the floor introduces no bias. For overlapping support, shared cells, or
sub-grid separations, the floor deliberately regularizes the singularity in
`d^theta` for `theta < 0`.

The consequence is that the implemented limit is:

```text
lim_{theta -> -infinity} d_eff(theta) = min_xy max(d_raw_xy, 0.5 km)
```

not the literal raw minimum when any pairwise distance is zero. Thus an
overlapping-support pair converges to 0.5 km under the default floor, not 0 km.
That is a modeling choice tied to raster resolution and should be reported in
any sensitivity table for the harmonic-dominated specifications.

## Planar fallback

The default path uses haversine great-circle distances and is the correct
choice for cross-country work and high-latitude comparisons. The planar
fallback is only for small regions. It now uses the shortest wrapped longitude
difference and scales longitude by `cos(latitude_mid)`, so small antimeridian
crossings at high latitude do not blow up into near-global distances.
