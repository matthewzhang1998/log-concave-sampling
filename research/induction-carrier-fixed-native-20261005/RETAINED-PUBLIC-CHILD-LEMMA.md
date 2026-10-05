# Integrated retained-public child-tree lemma

This spells out the only strengthened-looking interface in the root-seed argument. It is a property of the original LOW30 spine coupling; no internal Gaussian innovation is newly retained.

Let Q be the tuple of original marked publics in an oriented descendant subtree. Include the primitive distinguished leaf passive in Q. Let I^a and I^f be its actual and ideal forward outputs, using the same finite source fields and coefficient record E. There exists a coupling kernel C_E(q,d i_a,d i_f) with the correct two conditional output laws and

    integral |i_a-i_f|^2 C_E(q,d i_a,d i_f) dgamma(q) <= d(E)^2,
    d(E) <= Gamma sqrt(D) sum_{v in child tree}
                 rho_v^b product_{a strict ancestor within child} rho_a.

The physical public tuple is fixed in the coupling. The integrated norm, rather than a uniform bound at each q, is the assertion.

## Construction and estimate

1. Follow the distinguished spine from its Gaussian leaf passive R to the child-tree output. Freeze the spine's own nonleaf publics and all side-subtree outputs. They are independent of R. Conditional on these frozen variables, the ideal channels have fixed matrices; their inputs have standard Gaussian marginals after integration over R and prior channel innovations. At each chronological step apply the native conditional prior kernel at the actual input, couple to the ideal kernel, and use a common fresh ideal innovation to propagate the existing incoming error. Its error envelope is controlled exactly by LOW30's composition recurrence. R is never changed. None of the spine's own publics changes. The side subtrees and their publics are left intact during this stage.

2. Collapse the ideal spine to its terminal pair (R,Y). Conditional on the frozen side outputs and own publics, its cross matrix A is the ordered product of its channel matrices. The terminal representation can be taken as

       Y=A R+(I-AA^T)^(1/2) Z,

   with Z independent of the entire side-subtree banks and all public primitives. This is a new coupling representation of the completed ideal spine, not retention of its intermediate outputs. All intermediate spine outputs and innovations are consumed.

3. For each side subtree use the same construction recursively, keeping its own public tuple fixed. Its joint coupling depends only on that subtree's publics and private randomness (and E). In particular it does not depend on R, Z, or the enclosing spine's other independent publics. The changed terminal cross matrix A is uniformly input-to-matrix-HS Lipschitz, with the actual product of strict-ancestor source factors. The independent R converts its HS change to output L2. The fresh Z and the gapped covariance-root inequality do the same for the completion. These are exactly the two dimension-sharp estimates in LOW30's tree-replacement proof. No conditioning argument treats R as Gaussian after fixing its value: independence is used after the coupling has been constructed, in the integrated L2 estimate.

4. Sum the spine and side costs in the original spine-before-side order. All original marked publics survived identically. Source-private innovations, internal spine outputs, and the terminal collapse innovation remain private. The same ancestor-attenuated prior floor follows, with the fixed tree/path constants in Gamma.

## Appending a correlated physical carrier

Let T be the root compiler's primitive source-zero Gaussian bank, independent of the child tree, and let the physical packet carrier be H=L(T,Q,other affine publics). Thus H can correlate with Q but reads no child-private innovation. Under the common-carrier construction the marginal law of (T,Q,other affine publics) remains the original standard Gaussian bank. For fixed H and those public primitives, draw (I^a,I^f) using C_E(Q). This keeps H identical and gives the correct actual and ideal child kernels, because the child's private conditional law depends only on Q and E, not on T or H.

The integrated cost is exactly the same:

    E_{H,V} integral |i_a-i_f|^2 C_E(Q(H,V),d i_a,d i_f)
      = integral |i_a-i_f|^2 C_E(q,d i_a,d i_f) dgamma(q)
      <= d(E)^2.

Across packets, choose these couplings independently conditional on E,H and their independent packet null/public banks. No cross-packet public independence conditional on H is claimed. Every individual packet marginal remains the one used by the bound.

This lemma fails if the purported carrier reads an internal child innovation that is consumed during the spine comparison. The marked-leaf root's literal zero-amplitude map rules out that read in the current application: its own carrier is independent of its incoming passive, and all nonroot outputs enter only its vanishing residual. This check is structural, not an equality-in-law argument.
