# Complete linear-cost outer milestone

2026-10-04. The direct-mean outer theorem and its elementary source specialization have passed independent review.

## Admitted result

Let V be C2 with alpha I <= Hess V <= beta I, kappa=beta/alpha, 0<epsilon<=1, and supplied initialization |grad V(x0)|<=sqrt(alpha d). Set L=10+log((d+1)kappa/epsilon) and D=d+L. For each fixed j>=2, let P=2j-1.

There is a finite positive original-query sampler, using same/comparable posterior heats, with

    sqrt(beta) W2(Law X_out,pi)<=epsilon,
    Q<=Lambda_P {kappa^(2P)
           +kappa^(3-1/P)(D/epsilon^2)^(1-1/(2P))}.

Here Lambda_P is a fixed power of public/precision logarithms times constants which may depend rapidly on the FIXED requested order P. This is not uniform as P grows with the input size. Initialization if not supplied is a separate cost. The original finite numerical-query model is retained. Arithmetic, storage and Gaussian generation are not constant-time queries.

For arbitrary fixed requested target T>=3, choose j=ceil((T+1)/2). The actual accumulated order P=2j-1 is at least T, and the complete per-phase inverse-heat exponent P-1<T+1 is linear in T.

## Exactly what “linear” measures

- Source state-law order: R=j+1/2.
- Complete hidden source and force-mean cost: A^-2(j-1)=A^-(P-1), including all independent bank copies, paired coarse sides, actual ancestor histories and finite numerical work.
- Force approximation: ||ghat-m_A||_2 <=Lambda_P sqrt(D) A^(j-1/2), plus separately budgeted caller-dependent numerical error.
- With A comparable to H^2 up to public logarithms, one physical phase costs error Lambda_P sqrt(D) h H^P.
- Its direct grade triple is (I,P_*,F_*)=(P,P,P+2); no future-rank improvement is used.
- The arbitrary-order smoothed kinetic reference has stationary local defect h^(2j+1) plus its separately bounded collocation interpolation term. Its accumulated capacity 2j exceeds P.
- Full chronological mixing/reset accumulation is Lambda_P kappa sqrt(D) H^P.
- There are Lambda_P kappa/H phases, hence total work is Lambda_P kappa H^-(2P-1). The stated rate is the result of the actual H choice, including the kappa^(2P) conditioning branch.

Thus this proves a complete any-fixed-order outer construction with a linear inverse-heat query exponent. It establishes no algorithmic improvement over known/specialized samplers. Its dimension/accuracy exponent tends to 1, and its conditioning dependence also worsens with order. The specialized Exact32 result is much better in the displayed dimension rate. This is an architectural completion milestone and a baseline for future cost reductions.

## Why this route closes

The source's complete multilevel mean is already an executable vector with actual strong error A^j and a stable cost fixed point. Use that mean directly at each hidden exact-score Picard target; pay its whole fluctuation in ordinary L2 coupling. No outer covariance/even correction, marked cap, proxy, retained current, root reserve or future-rank argument is needed.

The selected hidden history has the original physical h^2 gradient rows and captured once-per-phase momentum refresh. No sampled ancestor becomes a deterministic caller. Every source bank expands its full hierarchy. The elementary specialization starts from a finite-mode Gaussian and uses quarter maps only, so no inherited marked-domain restriction D<=A^-b is present. All small-heat guards depend only on fixed order and public logarithms, apart from the explicitly factored sqrt(D) errors.

A real finite-program omission was found and repaired: the unanchored quarter predictor carries A^(M+1)|grad V(b)| even at an exact mode seed. Its fixed iteration count is chosen from the original numerical/caller budget; this residual is retained and deterministic zero trajectories are not falsely declared exact. An erroneous extra sqrt(D) on the recursive physical state-law error was also removed. The admitted versions below contain both repairs.

## Admitted proof and audit pins

Under direct-mean/:

- ANY-ORDER-DIRECT-MEAN-OUTER.md
  f7514fbaf777a573d262416759f2776dcce7bbc909741ceeb08a6276671785f4
- ELEMENTARY-SEED-AND-COMPLETE-SOURCE.md
  47e48da4e92885a48ce681048594d0dc8efbb1f71338e9b6df759f32053819cd
- independent-audit/INDEPENDENT-DIRECT-MEAN-OUTER-AUDIT.md
  1696d5059b9ca4afdf31f4b07aee2fd1a8105605ccc118c328e6c0265215baf6

Author check_direct_mean_outer.py: 6629 passing assertions.
Independent check: 213 grouped checks, including 1200 scalar finite-reference contraction cases and an exact symbolic regression example for the repaired predictor residual.

The audit directly inspected the original LOW30 all-fixed-order collocation theorem, the conditional decoder identity, physical history chronology, contraction/reset and original numerical-caller assumptions. The diagnostics supplement the analytical proof; they are not a full nonlinear implementation of the sampler.

## Publication scope

The unaudited auxiliary parameter and grammar discussion in the original milestone is omitted from this publication copy. Neither branch is used by the admitted direct-mean theorem. Its full proof, elementary specialization, independent audit, and numerical/initialization qualifications are retained above. Original and public hashes are separately recorded in the inventory.
