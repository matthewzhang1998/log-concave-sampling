# Independent repaired-source audit

Audited: `NONLINEAR-ANCESTRY-BRIDGE-REPAIR.md`, sections 3–8, against the pinned parent source's definitions of the canonical history and half-variance adapter. This audit does not recertify the native imported compilers.

## Verdict

No substantive mathematical defect found in the shifted bias proof, the noncommuting first/curl ports, the stated normalized bounds, the conditional energy/origin bounds, or the VALUE-error floor. The order-three bias qualification is genuinely stronger than unrestricted C2: bounded block dimension and Lip(Dg_block)<=BA are used exactly where required.

Resolved interface clarifications (independently rechecked in the revised theorem):

1. The revised theorem explicitly displays the linear Taylor discrepancy as the conditional expectation E[Dg(S)(F1-F2-g(g(w))) | x,v,w]. The unconditioned bracket is correctly used as an L2 coupling bound.
2. The revised theorem explicitly distinguishes the raw anchored caller bounds ACx and A from the half-variance-normalized bounds sqrt(2)ACx and sqrt(2)A used by the native caller guards.
3. The revised theorem explicitly defines Lambda as the collected guarded completion constant, absorbing the fixed factors from the imported bracket and source-energy bound.

All three points are resolved. No outstanding defect was found in the audited scope. The added main-theorem Section 4.1 agrees with the independently proved fixed-(D,h) corollary in Section 6 of this audit.

## 1. Shifted cancellation and dimension dependence

I checked the expansion at the executed S=x-g(v)+g(g(w)), rather than replacing Dg(S) by Dg(x). The exact bridge's linear term is

Dg(S)[g(v)-E(F1 | x,v,w)],

while the canonical term is

Dg(S)[g(v)-g(g(w))-E(F2 | x,v,w)].

Their true-minus-repaired difference is therefore the conditional expectation stated above. On the actual history, ||F1-F2||2<=A^2 sqrt(d), and ||g(g(w))||2<=A^2 sqrt(3/8) sqrt(d). Multiplication by Dg(S) yields precisely A^3(1+k)sqrt(d). No Hessian commutation or unshifted response approximation enters this step.

For a d-dimensional block, the Gaussian L4 factor is [d(d+2)]^(1/4). The actual displacement bounds are Au_A times that factor for B_shift-F2 and Ad0 times that factor for each bridge difference. Lip(Dg)<=BA then produces the two remainder coefficients BA^3u_A^2/2 and BA^3d0^2/2, each multiplied by sqrt(d(d+2)). The centered difference has the stated factor 1/2: its two Taylor remainder bounds are averaged, not added without that average.

Positive weights and conditional Jensen suffice; shared L across nodes does not invalidate either estimate. Summing squared block norms uses sqrt(d(d+2))<=sqrt(b+2)sqrt(d), exactly yielding the displayed C_(b,B,A). The bridge's operator error is applied only to g and then multiplied by the bounded random matrix Dg(S), giving delta_br A^2 sqrt(D). The outer norm bound C_F also checks directly.

## 2. Exact noncommuting derivatives and curl

Differentiating the actual centered difference gives C_j S_y+K_j(D_j)_y with the two stated matrix averages. All four S_y and (D_j)_y rows are correct, including the order JW. I found no omitted S-dependence, no missing D_j path, and no illicit commutation.

The key norm bounds use the PSD interval: ||H-H0||<=A, ||(H_j^+-H_j^-)/2||<=A/2. Thus both ||Hbar|| and ||Hbar-H0|| are <=3A/2. Applying these to the complete rows gives exactly the source's Cx,Cn,Cm bounds.

For the lifted curl, the leading G block is a sum of symmetric matrices Hbar-H0. The remaining G skew is bounded by beta(3A eta+3A^2). The off-diagonal block row contributes A sqrt(Cn^2+Cm^2+A^2). This is the operator norm of a skew block matrix, so no extra factor two is needed for that off-diagonal contribution.

Half-variance normalization multiplies both first and curl by sqrt(2). After dividing first by A and curl by A^2, the displayed upper bounds are nondecreasing in A and beta. At A=1/2 and beta=sqrt(3)/2,

(normalized first/A)^2 <= 32231/2048 < 16,

normalized curl/A^2 <= (sqrt(2)/16)[39sqrt(3)+sqrt(2381)]
                       =10.283584122805527... <12.

Hence ell_E=4A and a_seed=3A safely satisfy the declared normalized radius and curl product bound. The first estimate is relatively close to 4A, but it is valid with a strict exact-rational margin.

The independent executable check also compared the complete forward derivative to centered finite differences for a nonlinear two-dimensional family with noncommuting Hessians, with two shared-root inner and outer nodes. Across 100 random points, the largest operator discrepancy was 2.964e-11. This is a supporting diagnostic, not a substitute for the preceding calculation.

## 3. Energy and actual caller origins

A direct deterministic inequality is

|E| <= A^2 sum_i omega_i [2|v_i|+A|w_i|+sum_j p_j|U_ij|].

Given Z=z, the conditional Gaussian means are t_i z/2, t_i z/4, and r_j t_i z. Their noise standard-deviation factors are bounded by 1/sqrt(2), k, and 1. The two exact first moments therefore give the coefficient

2(1/2)(1/2) + A(1/4)(1/2) + (1/2)(1/2)
=3/4+A/8

on |z|. The noise coefficient is sqrt(2)+Ak+1. This independently recovers the quoted conditional Lp energy bound without assuming private independence across nodes.

Setting all private roots to zero gives the actual captured program, including the nonzero U_ij=r_j t_i z. Removing the Gaussian terms in the preceding inequality proves the stated origin bound. Differentiating the captured execution and subtracting it doubles the raw caller estimate but does not alter private first/curl. The baseline assertions check by the same argument. Literal zeros require the stated coherent zero-site/numerical conventions; a nonzero caller origin is never erased.

## 4. Absolute VALUE floors

Let Delta_b<=nu_b+A nu_w, Delta_S<=nu_v+Delta_b, and Delta_Dj<=nu_v+nu_Uj. At the actual perturbed sites, the base g(S) contributes nu_S+A Delta_S. Each centered difference contributes

(nu_+j+nu_-j)/2 + A Delta_S + A Delta_Dj.

Mass one then yields exactly

nu_S + sum_j p_j(nu_+j+nu_-j)/2
+3A nu_v+2A nu_b+2A^2 nu_w+A sum_j p_j nu_Uj.

Thus the uniform floor (2+6A+2A^2)nu is valid. The extra baseline nu for E and the doubling for a separately evaluated captured origin are also valid. None of these estimates divides by energy or an A-dependent small quantity.

## 5. Completion interface and exclusions

The anchored residual energy remains O(A^2)(|z|+sqrt(D)). With the declared parameters, the imported near-gradient bracket is

16A^2+64A^3+64A^(5/2)=O(A^2).

Consequently the claimed order-four completion allowance follows **conditional on the stated native service contract and all its guards**, with universal constants absorbed appropriately. The analysis does not prove those services, remove their guards, or validate their finite-mode/precision/replay implementation. The source VALUE counts and the 4D residual private dimension match the executed graph.

Companion diagnostic: `independent_repaired_source_checks.py`. It certifies the exact endpoint first constant and checks noncommuting derivative/curl rows; it explicitly does not certify native imported compilers.

## 6. Optional fixed-source little-o corollary

The following additional analytical conclusion is valid without the block-curvature qualification. Fix D and an anchored C1 gradient h with 0<=Dh<=I, and set g_A=A h. Choose any of the positive mass-one bridge rules with delta_br<=A. Then

||psi_br,A-psi2,A||_(L2 gamma_D)=o(A^2) as A decreases to zero.

There is no uniform rate in h or D. This is a bias corollary, not an extension of native compiler eligibility.

To check the nontrivial uniformity over A-dependent nodes, put R=(|x|^2+|N|^2+|M|^2+|L|^2)^(1/2). Every real bridge row has coefficient norm one, so |U_s|<=R for every s. Also |v|<=R/sqrt(2), |w|<=kR. For A<=1/2, all points S and S±tD_s, 0<=t<=1, lie in the ball of radius 3AR about x. Define the local derivative oscillation

omega_x(r)=sup_{y,z in B(x,r)} ||Dh(y)-Dh(z)||.

It tends to zero as r tends to zero and is bounded by 2. The centered Taylor remainder therefore satisfies, simultaneously for every s,

|R_A,s|/A^2 <= d0 R omega_x(3AR).

Dominated convergence in L2 proves that the right side has norm tending to zero. Hence any A-dependent positive mass-one combination of node remainders is o(A^2).

For the true history, define G0=int e^-t|X_t|dt and G1=int t e^-t|X_t|dt. Both are in L2 (indeed every fixed finite Lp). The actual F2 obeys |F2|<=A(G0+A G1), while |B_shift|<=A(|v|+A|w|). Thus |B_shift-F2|/A is bounded by a fixed L2 random variable, and every point on its shifted Taylor segment tends to x. The same derivative-oscillation argument and dominated convergence make the true remainder o(A^2).

The exact linear mismatch and bridge quadrature error are O(A^3 sqrt(D)). This proves the corollary. Applying a finite outer rule preserves it provided delta_out=o(A); the main theorem's delta_out<=A^3 is more than sufficient. Absolute numerical floors must still be retained separately, or be chosen o(A^2) to include them in the little-o conclusion.
