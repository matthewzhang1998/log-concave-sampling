# Independent covariance-test tensor audit

2026-10-05. Verdict: **PASS for the finite-dimensional real tensor lemma.** The source-specific variance theorem and positive grouped producer in Sections 2–3 remain **OPEN**.

Audited source: `../COVARIANCE-TEST-PORT-REDUCTION.md`, SHA-256 `e8424b630508dab9ed042f9275a1f414d662a011cfa2e66f240b3f573b4bcd9f`. This supersedes the initially supplied `5a30772efcb0874cae742ef21809b7e1343520475bb8c2074ce93c5377c29575`; the revisions correctly restrict the automatic bound `h ≤ K√D` to rank `r ≥ 2` and explicitly specify real tensor spaces. This audit did not edit the source.

## 1. Operator Cauchy–Schwarz and every cut

Use real tensor spaces, `r ≥ 1`, and nonnegative `K,v`. For a fixed external cut, local reshapes have types `A: X_A → Y_A` and `B: X_B → Y_B`. A local empty row-index set means a **row** vector; a full row-index set means a **column** vector. Singleton scalar spaces must not be silently identified with the full old tensor space.

The rectangular operator inequality is valid:

    ||E[A ⊗ B]|| ≤ ||E[AA*]||^(1/2) ||E[B*B]||^(1/2).

Indeed, factor `A ⊗ B = (A ⊗ I_YB)(I_XA ⊗ B)` and apply operator Cauchy–Schwarz to the averaged product. Its second orientation is

    ||E[A ⊗ B]|| ≤ ||E[A*A]||^(1/2) ||E[BB*]||^(1/2).

Permuting row or column slots does not change operator norms. Apply these inequalities to the centered local reshapes:

- **Both old blocks split:** both factors are proper local cuts. The pointwise `2K` bounds give `4K²`.
- **Exactly one whole old block:** if it is a column, choose its covariance Gram `E[AA*]`; if it is a row, choose `E[A*A]`. Equivalently interchange the roles of A and B when necessary. The covariance Gram is bounded by `v²I`, and the split factor contributes at most `4K²I`, giving `2Kv`.
- **Both old blocks whole:** a proper global cut places them on opposite sides. The flattening is the covariance matrix, or its transpose up to slot permutations, so its norm is at most `v²`. Placing both whole blocks on the same side is an improper global cut and is excluded.

These cases exhaust all proper cuts. Using the other Gram of an uncut block gives the scalar `h²`, not the covariance norm; that orientation is still a valid inequality but loses the claimed gain. The source's “or reverse orientation” is therefore essential and sufficient.

## 2. Hilbert bound, rank-one edge, and necessity

With `z = vec(F−EF)`, the balanced matrix is `Σ = E[zzᵀ] ≥ 0`. Thus

    ||T||HS² = Tr(Σ²) ≤ ||Σ||op Tr(Σ) ≤ v² h².

For `r ≥ 2`, isolate one physical slot in **uncentered** F. Its flattening has at most D singular values, each at most K, so

    h² = E||F||HS² − ||EF||HS² ≤ E||F||HS² ≤ DK².

For `r = 1`, there is no proper local cut, and only the two orientations of the balanced global cut remain. The main covariance and Hilbert bounds still hold with the separately declared h. The revised source states this correctly.

The scalar-test inequality with constant b is exactly equivalent to `||Σ||op ≤ b`. Therefore it is necessary and sufficient for a balanced-cut bound with that same constant. A bound on all cuts by `M = max(4K²,2Kv,v²)` only entails the scalar-test bound with M, not necessarily with the smaller `v²`; the source's qualifier “with its corresponding constant” handles this distinction.

Proper local cuts alone cannot supply dimension-free balanced control: for a Rademacher sign ε and `F = ε Σ_i e_i⊗e_i⊗e_i`, every proper local cut has norm 1, while the balanced covariance norm is D. The source's Gaussian `sin(G)` version gives `D Var(sin G) = D(1−e^(−2))/2`. Neither example establishes anything about the source-qualified higher-heat-jet genealogy.

## 3. Optional strengthening and one scope clarification

The stated constants are conservative. For any proper uncentered local reshape A,

    E[Abar Abar*] = E[AA*] − (EA)(EA)* ≤ K²I,
    E[Abar* Abar] = E[A*A] − (EA)*(EA) ≤ K²I.

Using these bounds in the same argument improves the all-cut bound to `max(K²,Kv,v²)`. This is optional and does not change the open source-specific obligations.

The final source now explicitly says **real**, resolving the only additional scope clarification raised in this audit: `E[Fbar ⊗ Fbar]` need not be positive semidefinite over a complex tensor space. For a complex version, use the conjugated second block and the Hermitian variance convention. No correction is needed when real tensors are understood, as in the Gaussian-gradient setting.

## 4. Reproducible check and limits

Run `python check_covariance_tensor.py` from this directory. The checker verifies the pinned source hash; enumerates 394 proper cuts across random real rank-1 through rank-4 families and the rank-3 diagonal-sign family; checks both operator Cauchy–Schwarz orientations, both sets of case constants, and the Hilbert bounds; and checks empty/full local Gram orientation on an isotropic old tensor space. `checks.json` records **PASS**. Numerical checks supplement the proof, rather than establish the arbitrary-rank statement.

Sections 2–3 of the source explicitly preserve both remaining obligations: uniform source-qualified scalar-test covariance control for the complete higher-heat-jet coefficient, and a positive sign-programmable grouped VALUE/current producer respecting shields, descendants, costs, and observer boundaries. No such theorem, producer, recurrence, or eventual-sublinear conclusion is proved or claimed by this audit.
