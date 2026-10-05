# Exact whole-bank connected tree census with unchanged private heat

2026-10-05. A coefficient-level source map for higher shared-bank cumulants. It supplies a positive Gaussian replica geometry and native bridge graph; it does not supply the missing uniform higher-clock/heat ledger.

## 1. The complete old bank is one variable

Let B be the full standard Gaussian coefficient bank at the boundary, with zero source-zero physical readout. Let F_i(B) be finitely privately heated coefficient tensors, with their original scalar clocks, affine query ancestry and physical slots retained. If several old histories share B, expand their SUM before applying the formula below. Do not replace it by a sum of separate covariances/cumulants.

For r>=2 the joint Gaussian cumulant has the connected forest representation

    kappa_B(F_1,...,F_r)
      = sum_(labelled trees T on {1,...,r})
          integral_[0,1]^(r-1)
            E_(B_1,...,B_r; K_T(t) tensor I)
            [ product_((i,j) in T) (D_(B_i) dot D_(B_j))
                product_i F_i(B_i) ] dt.

Here K_ii=1 and K_ij is the MINIMUM edge parameter on the unique i-j tree path. All tensor physical slots are retained. Each derivative contraction is in the COMPLETE original bank, using its actual injection row at whichever original gradient occurrence is hit.

The coefficient 1/r! required by a logarithmic generating function is inserted separately once. Derivative product rules enumerate every hit history, including repeated hits at one occurrence. They are never interpreted as tensor oracle queries.

## 2. Why the replica covariance and weights are positive

For 0<=u<=1 delete every tree edge with t_e<u. Let Pi_u be the resulting partition of replicas. Its connectivity matrix has entry one within a block and zero otherwise; it is positive semidefinite. The exact identity

    K_T(t) = integral_0^1 connectivity(Pi_u) du

holds entrywise, including the diagonal. Hence K_T is positive semidefinite. It is a known r-by-r scalar covariance matrix; use its known scalar root, or the nested coalescent block representation after sorting the edge parameters. No force covariance or unknown tensor is inverted.

Let c=min_e t_e. Since the tree remains connected for u<=c,

    K_T(t)-c 11* >= 0.

Thus a symmetric disintegration preserving the ACTUAL retained old B is

    B_i = sqrt(c) B + (R_T V)_i,
    R_T R_T* = K_T(t)-c 11*.

V is a fresh r-copy complete standard bank independent of B. Each B_i has its exact original marginal law, and cross-correlation is exactly K_ij. This extends the two-copy symmetric bridge without independently replacing the old service's B. The original caller means and every original private force shield sigma_v are unchanged.

For r=2 this gives the prior formula B_+=sqrt(c)B+sqrt(1-c)U_+, B_-=sqrt(c)B+sqrt(1-c)U_-. The derivative row is the original old injection A, not sqrt(c)A: the formula represents a smoothed derivative, not the derivative of a smoothed function.

## 3. Finite derivation of the connected forest identity

The finite forest interpolation identity for a smooth function of pair parameters is

    f(1) = sum_forests F integral_[0,1]^|F|
              (product_(e in F) partial_e) f(x^F(t)) dt,

where x^F_ij=min-path t for vertices in one component and zero across components. One elementary proof on monomials is probabilistic. For f(x)=product_e x_e^(n_e), give each support edge an independent random weight with density n_e t^(n_e-1). A term is nonzero precisely when F is a spanning forest of the support graph. Conditional on its tree weights, the event that F is the maximum-weight spanning forest is exactly

    X_e <= min_(path_F(e)) t  for every nonforest edge e.

Its conditional probability is the remaining monomial factor in the forest integral. These mutually exclusive maximum-spanning-forest events partition probability one. Linearity gives the identity for polynomials; smooth approximation with the finite required derivative norms gives its finite-order form.

Apply it to Gaussian expectations with off-diagonal replica covariances scaled first by a small lambda>0, so that the entire pair-parameter cube has positive definite covariance. Gaussian covariance differentiation contributes D_i dot D_j, with no extra 1/2 because the two symmetric off-diagonal entries are differentiated together. The identity extends to 0<=lambda<1 by analyticity along the positive forest covariances; privately heated coefficient functions have the required finite Gaussian moments and derivatives. Dominated convergence takes lambda to one. Alternatively, the same finite formula follows after a positive diagonal regularization and its limit.

Insert the forest expansion into the partition/Mobius formula for the joint cumulant. A disconnected forest factors over its components, and its partition coefficient is zero. Exactly the connected forests, hence labelled trees, remain. This proves the displayed representation.

At a requested finite current cutoff, only finitely many cumulant orders and derivative histories are kept. No infinite polynomial moment-generating series is executed or assumed to converge.

## 4. Native source graph and marked spanning tree

For each tensor coefficient F_i retain its original force graph and marked-spanning-tree certificate. The tree interpolation differentiates at its original occurrences. A hit raises C_k to C_(k+1) and exposes one extra bank slot. Contract each pair of hit slots along the corresponding replica tree edge.

The union of the original certified spanning trees and these r-1 bridge edges is itself a spanning tree on all original force occurrences. Its remaining leaves were leaves of an original tree and still have their physical marks. A hit at an original marked C0 leaf therefore becomes a marked C1 vertex with a bank edge; it is admitted by the new invariant. Multiple hits simply increase its adapter order and valence together.

The global all-proper-cut and one-Hilbert theorem applies. The finite original-gradient producer uses those C_k VALUE sources at the correlated replica centers and their UNCHANGED original private shields. Cut/opening Gaussians are a different dynamic bank from the replica coefficient roots. Known signs, permutations, clock weights and inverse readout normalizations are placed at the selected root.

This is a concrete source map for the full shared-bank cumulant, including cross-history terms. It does not average a random tensor by a Monte Carlo bank, and it does not erase any old observer.

## 5. Exact remaining analytical guard

The r-copy positive geometry alone does not establish a useful uniform endpoint estimate. Derivative hits can accumulate at one original force. Its normalized factor has the actual inverse shield power t_v^(-k_v), and its complete caller first has one additional t_v^-1. Repeated original clocks must remain repeated; independent-node bounds cannot be substituted.

The two-copy symmetric redistribution proved in the pinned source pays the tight squared-clock one-hit sum. A higher tree requires an explicit multi-copy redistribution/scale decomposition that pays all corresponding derivative powers while preserving K_T(t) and the original heat. Using only the smallest gap 1-max_e t_e generally throws away the distinct tree-clock resources. Merely stating that K_T is positive is not that estimate.

A deliberate common heat floor gives a finite conservative certificate, with all inverse powers entered into root amplitudes, radius guards, first paths, native priors and precision. But the mismatch between that heated coefficient and the original target remains a separate current and must be restored at the claimed grade. This note does not claim that restoration or an arbitrary-order sampler.

The pure conditional queue for a certified native packet doubles its root surplus. Whole-bank cumulants add root excesses of their arguments, but a complete mixed queue must still prove that every argument contributes strictly positive excess AFTER its actual heat and observer losses. A source-zero nonzero old-bank physical row invalidates that claim and can generate infinitely many same-grade ranks, as the retained-public obstruction demonstrates.

## 6. Equality scope and honest finite bill

The symmetric retained-B disintegration realizes the displayed cumulant only AFTER averaging the retained B together with the new residual replicas. It is not a pointwise-in-old-B cancellation. Its conditional random coefficient is a new old-bank-dependent field, and its mixed cumulants with the original service must remain in the complete joined ledger. This distinction is exactly where the four-term decorated correction is stronger: that correction matches its original conditional coefficient pointwise in B.

Fix a finite derivative/history list and an absolute quadrature tolerance. A finite positive scalar cubature exists by ordinary quadrature of the smooth privately shielded finite integrands. Split into ordering sectors of the tree parameters if using the min-path representation. Let N_T be the ACTUAL certified positive node count with its error bound; do not label N_T a public-log quantity before proving the needed endpoint/analyticity estimate. A simple finite grid is a legitimate existence fallback but its possibly large cost must remain visible.

If history h differentiates r original coefficient graphs into a final graph with N_h force vertices and M_h internal edges, its per-node producer bill is

    sum_(v in h) Q_C(k_v),b(v)
       + Q_complete_old_captures + Q_known_replica_root
       + Q_known_physical_rows + Q_readout + Q_scalar_encoding,

plus M_h-N_h+1 cut-edge D-roots, every physical row, every complete original source/pair/filter/private tape and the reserved keep. Multiply by the full ordered derivative-hit list, original shared-clock/history combinations, labelled tree count and certified N_T. For r labelled distinct arguments the tree count is r^(r-2); repeated histories are reduced only under a fixed verified multiplicity convention. Add all required old source-version rebuilds and ancestor replay. The replica root is a known r-by-r scalar matrix acting on the full original bank and is charged in both scalar work and Gaussian storage.

The substantive source-current estimates spend one Hilbert factor through the marked-spanning-tree theorem. Numerical cubature may be assigned a stricter conservative absolute tolerance using its literal finite tensor bound; any extra dimension dependence in that numerical bill remains explicit. No hidden inverse-heat Monte Carlo sampling is asserted, and no unverified eventual-sublinear high-accuracy exponent follows.
