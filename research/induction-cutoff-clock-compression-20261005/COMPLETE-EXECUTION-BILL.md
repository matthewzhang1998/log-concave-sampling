# Exact callback and complete replay bill

This note distinguishes a finite scalar quadrature implementation, a finite original-VALUE source implementation, and a native sampler instantiation. The first two exist in the pinned imported packages; the last requires the same admitted native wrapper, explicit source floors and current ledger as before. No native sampler was run in this research package.

## Clock callback

Call `prepare_certified_rule(N,K,D,eta,epsilon,W,caller_weight=...)` in `cutoff_free_sector_gauss.py`. It returns exact parameter choices, certified rational one-interval Gauss nodes/weights, and the dyadic interval atlas. Iterate `sector_nodes(base,panels)` twice, swapping the two original tree-edge roles on the second pass. Returned `(t_large,t_small,w)` already includes the Jacobian; do not multiply by r again.

No F, tensor entry, expectation, derivative of g, or native law is an argument to the builder. Memory need not store all Q nodes: the scalar node iterator is lazy. The single-axis rule and the panel list are finite. Every emitted node still requires its original source execution.

## Original VALUE source

The pinned `active_probe_source.py` source is the literal finite callback used by the input construction. At order k, center z, private width t and normalization A, it evaluates

    (A t)^(-1) 2^(-k) sum_(eps in {-1,1}^k) product eps_j
        [g(z+a t(x+sum eps_j P_j))-g(z+a t sum eps_j P_j)],
    a=(k+1)^(-1/2).

It uses exactly 2^(k+1) original gradient VALUES before duplicate caching. Every pair is same-center and on the same complete bank. The optional radial preparation applies its explicit Jacobian to this field, preserves its exact private curl, and changes no VALUE count. The existing floating implementation is diagnostic; certified scalar geometry and the separately required native/source Sobolev floors remain real obligations.

For each node, use the original complete retained B and the sealed c,R,delta disintegration and simultaneous private-heat allocation. Generate all required complete V/E banks and query residual roots. Original injection rows and old six-packet residual definitions are retained. Build the marked force tree with its actual selected source spines, physical readouts, positive native fills, and untouched keep. The complex OU representation is never executed.

## Complete census

For a fixed ordered triple of original cubic histories there are 243 new hit assignments across three replica trees. The maximal raw main bill at one node/history is 44 original VALUES, with total primitive order 7 and maximum local order 3.

For a fixed complete old mixed history there are 648 new hit assignments across three replica trees. The maximal raw main bill is 72 VALUES, with total primitive order 10 and maximum local order 4. If the six-packet's nine OLD covariance hit pairs have not already been included in the old history definition, multiply its census by nine; the executable checks 9*648=5,832 expanded mixed hits. It does not multiply an already complete old census by nine again.

Let H3 be the number of complete ordered cubic old-history triples, including all literal finite clock/physical/source versions. If all wrappers on the actual descriptors have repetition bound Mmax, the new raw-source part is at most

    44 * 243 * H3 * Q_tree * Mmax.

For the mixed target with H336 complete old-history choices (including old six-packet hit labels), it is at most

    72 * 648 * H336 * Q_tree * Mmax.

Use the exact sum over local wrapper repetitions when they differ; these maxima are merely upper bounds. Repeated original clock labels remain repeated descriptors with literal weight products, not independent draws.

## Replay cannot be dropped

For every source invocation list its complete version key:

- original g/service version and actual requested VALUE floor
- source order, same-center sign/anchor identity and radial preparation
- all captured caller values and strict ancestor versions
- original and new private widths, roots and frozen clock labels
- selected/native pair, filter, readout, positive gap and numerical floor versions

A cached value is reusable only for an identical complete key. For a changed key, re-evaluate all dependent original VALUES and their same-center anchors, all affected caller/ancestor captures, and every native wrapper invocation at the requested floor. Scalar root preparation, physical rows, residual banks, fills, final keep, and retained internal replay storage have their own actual costs. No tensor oracle replaces this dependency graph.

Consequently the literal work expression is

    sum_(actual complete native source invocations v)
        2^(k_v+1) M_native(v)
    + all old complete-version/capture calls
    + all changed-key ancestor and anchor replays
    + scalar/vector geometry, roots, readout and keep work.

The scalar quadrature theorem changes the number of new bridge descriptors to Q_tree. If an admitted per-descriptor complete old/native/replay service costs C_complete at the chosen actual floors, its contribution is multiplied by the actual finite hit/history/node census; C_complete is not replaced by one. The same requested floor used to certify the final error must appear in that service's bill. A local fixed-rank LOW30 logarithmic floor dependence is usable only under its stated guards.

For eta=alpha^q_eta and epsilon=alpha^p, this new Q_tree contributes no inverse-alpha exponent at fixed rank. An old C_complete=alpha^-c, an inherited history count, a native inverse-width/gap/readout factor, and the full mixed caller eta^-chi still retain their own exponents. If these depend on p, they must remain in the final max-plus ledger. Thus this is a genuine removal of the quadrature cutoff exponent, not a claim that the complete original-VALUE exponent is zero or globally sublinear in p.

## Guard dependencies

The original target envelope W is fixed before choosing Q_tree. New native inverse readout shares may depend on total packet count. When their normalization cancels in the selected target coefficient, charge them in the native source/caller/error/replay guards rather than feeding them circularly into W. If a different implementation places an uncancelled Q-dependent normalization in its target, it must solve and certify its own parameter-selection fixed point.

Innovation gaps may reach the accuracy-dependent zeta. The displayed Gaussian program only square-roots them. A native wrapper that introduces an inverse gap must include that cost, which can be accuracy-dependent, in its complete bill. The node theorem alone does not certify that no such inverse is used.

No common-carrier native theorem is supplied here. Independently owned native groups retain their actual readout allocations and inverse products. Common-carrier reuse remains conditional on the separately required carrier-fixed replacement result.
