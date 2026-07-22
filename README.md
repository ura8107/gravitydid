# jwdid

An R port of Stata's `jwdid` extended two-way fixed-effects estimator.

Use `jwdid_aggte()` for post-estimation aggregation. The shorter name
`aggte()` is intentionally not exported because it collides with the
unrelated function exported by the `did` package. Stata users can use
`estat()`, and typed wrappers such as `jwdid_event()` are also available.
