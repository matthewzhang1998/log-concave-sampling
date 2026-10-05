# Independent audit: source-qualified uniform covariance countertest

2026-10-05.

## Verdict

PASS for the proof in `../MODULATED-C2-SIX-TREE-COUNTEREXAMPLE.md`, SHA-256
`8d0086e37b9cb1273ff483dc3dd9ddbeee297c9e93b7675b4aad11f2d6bd4dc3`.

The printed bound

    Var <H_n,J_n> >= n / 2^40,  n>=128,

holds for the source, exact illustrative clocks, normalized tensor convention, and selected G1/G1 center-hit history specified there, uniformly over its complete Price interval. Every proper cut is bounded by a numerical constant, while H_n is an HS-unit symmetric six-tensor. This is a genuine bounded-PSD-Hessian original-gradient counterexample to the proposed unweighted shield-uniform covariance scalar-test theorem.

The result is not an impossibility theorem for the general positive-law program. In particular no lower bound for an arbitrary original-clock weighted aggregate, arbitrary signed correction-history sum, or redesigned common-heat/coarse-observer grouping follows merely from this fixed-node result.

## Exact genealogy and normalizations

The relative source path `../../../EXACT-FIVE-CLOCK-NATIVE-SOURCE-ADAPTER.md` in the main note resolves correctly from its directory. The original half-shield splitting is retained exactly. With retained Z=0 and old q=s=1/2, the G0,G1 center row is `(r2 sqrt(3)/4, r2 sqrt(3)/2)`, the W2 coefficient is sigma, and the other two coarse-block coefficients are zero. Thus

    Var(q2_i)=(15/16)r2²+sigma²=15/16-7/(8n).

The entire correlated 5D-dimensional Price pair produces independent coordinate pairs, each with precisely this variance and correlation t. The independent low-coordinate pair supplies the common modulation. No private shield is enlarged, and no original ancestor is resampled.

The selected G1/G1 injection is

    gamma_n=(3/4)r2²=(3/4)(1-2/n),

which is nonvanishing and at least 1/2 in the claimed range. It is not the W2-only injection, whose squared coefficient sigma² would invalidate the stated growing-variance conclusion. Selecting a center product-rule hit is legitimate even though the same G1 block also changes the right leaf query: those leaf hits are different histories. Summing center-center injections gives v_n instead, also nonvanishing.

The center C2 normalization `sigma² P_sigma D³g/A` is correct: the old center C1 coefficient is `sigma P_sigma D²g/A`, so one center query derivative costs precisely `1/sigma` after conversion to C2. Combining both hits with the old cubic scalar restores `A^7 beta`, beta=`w²/sigma^4`, in the convention of the sealed tight bridge. Four unchanged leaf C0 coefficients each contribute kappa in the comparison tree, giving kappa^4, and the two center jets contribute eta² d_n². The all-i six-slot component has only internal index i or 0. The mixed term's sigma² factor and sign in (12) are exact.

Symmetrization leaves the H_n pairing unchanged. Fixed nonzero calibration/readout conventions alter the numerical prefactor only; the printed `2^-40` is for the explicit tensor convention (10), not every arbitrarily rescaled convention.

## Constants and analytic bounds

The Hessian diagonal perturbation has norm at most eta and its off-diagonal arrow norm is at most eta sigma sqrt(n)=eta. With eta=1/8, kappa=1/2, this gives the printed [1/4,3/4] sandwich globally.

The heat multiplier is exactly `exp[-h²(1+sigma^-2)/2]` for every perturbation term, including mixed and 00 entries. The actual leaf perturbation norm is therefore at most `(1/4)exp[-(n+1)/8]`.

For the printed variance proof:

- `v_n>=7/8`, `mu_n>=1/3`, and `|nu_n|<=1/2` hold for n>=16 over the whole t interval.
- `sd(cos² U)=(1-exp(-4v_n))/sqrt(8)>=1/3` is exact and sufficient.
- `||cos U cos V-cos² U||_2<=||V-U||_2<=2 n^-3/2` is valid, and implies `sd(M)>=1/4`.
- `sd(mu M+sigma² nu N)>=1/12-1/(2n)>=1/32` holds in the claimed range.
- `c_n>=2^-14` follows from gamma>=1/2, kappa^4 eta²=1/1024, and `d_n²>=exp(-2)>1/8`.
- Hence `sd<H,J0>>=2^-19 sqrt(n)`.

The 16-placement center-tensor cut estimate is sound. For each placement with k zero indices, the coefficient is bounded by eta d_n sigma^k; the only possible sqrt(n) factor arises when the repeated positive indices all lie on one side of the cut. For k>=1 it is absorbed by sigma^k sqrt(n)<=1. For k=0 both sides necessarily carry a repeated positive index. The all-zero term is separately at most eta n sigma^4<=eta. This proves the uniform cut bound 16eta<=2.

The one-edge contraction lemma therefore applies after all leaf linear maps. In each of the four telescoping differences, both center bounds contribute 2, one leaf contributes its error, and every other leaf has norm at most one. This gives the printed all-proper-cut error `4exp[-(n+1)/8]`. Isolating a physical index gives HS error at most that quantity times sqrt(n+1), and the factor 6 in (18) is conservative.

At n=128, `6exp[-129/8]` is approximately 5.96e-7, below `2^-20` (approximately 9.54e-7), and decreases thereafter. Standard deviation is 1-Lipschitz in L2 because subtraction of the expectation is an orthogonal projection. Thus there is no extra factor 2 when converting the pointwise error to the final centered-scalar estimate. The claimed final lower bound follows.

The asymptotic coefficient in (17) is also correct: the averaged high-coordinate product tends to 1/2, and `Var(cos²(sqrt(v)Z))=(1-exp(-4v))²/8`, giving the printed denominator 32.

## Clock support and weighted scope

The illustrative numbers q=s=1/2, r1=r3=1/sqrt(2), sigma²=1/n are legitimate parameters of the analytical five-clock target. The proof does not establish that these exact values occur in any particular frozen finite quadrature. The main note correctly states the robust extension: take any actual endpoint center-node sequence sigma->0, set n=floor(sigma^-2), keep the other clock nodes in fixed interior compact sets, and require A/sigma²->0. Then n sigma²<=1 and approaches one; all source guards survive and the same positive limiting variance remains. The explicit `n/2^40` need not survive unchanged for arbitrary different interior clock constants, but a positive constant times n does.

For the illustrative A=n^-3=sigma^6, a cutoff sigma_min=A^b with b>1/6 includes this center scale. A separately restricted node family excluding the required scale is outside this counterexample. More generally any endpoint scale sigma comparable to A^a with 0<a<min(b,1/2) suffices for the robust version.

The uniform t proof survives positive normalized Price averaging on its declared common G,H family. It also survives a randomly selected independent Price clock. It should not be confused with summing independently resampled quadrature histories, where squared weights and the actual bank architecture must be retained.

At a fixed chosen original node, in the exact normalization C=A^7 beta J, scalar restoration is rigorous:

    Var <H_n,C_n> >= A_n^14 beta_n² n / 2^40.

If fixed readout constants remain outside that equation, their square also remains in this lower bound. Deriving an original-clock integrated lower bound, a no-go theorem for all grouped clocks, or a statement about a different positive producer requires the actual node weights, cross-node correlations, and all signs. Those conclusions are not supplied here.

A redesigned grouped common-heat/coarse-observer calculation remains a possible route, as do original-clock-weighted energy bounds with genuinely additional information. The invalid step would be to obtain such a route merely by invoking the now-falsified uniform covariance scalar-test inequality.

## Restricted separable test and conditional grouping (Sections 7–8)

The added Sections 7–8 also pass. For a coordinate-separable potential with deterministic retained caller and coordinatewise isotropic Gaussian ancestry, all primitive tensors are diagonal in the separating orthonormal basis. The scalar normalized center coefficient obeys the dimension-free estimate

    |sigma² partial_q² P_sigma(u_i″/A)(q)|
      = |E[(Z²-1)(u_i″/A)(q+sigma Z)]|
      <= E|Z²-1| <= sqrt(2).

This follows from the bounded original Hessian and the exact second Hermite integration-by-parts formula. Leaf factors and the declared injection are bounded, so |j_i|<=K. Different coordinate banks remain independent under the full Price pair. Thus the variance identity and bound (20) are valid, including after an orthogonal change to the separating basis. The deterministic retained-caller condition matters: an additional randomly shared caller would require conditioning or its own covariance analysis.

For the nonseparable example, conditioning on U,V indeed leaves independent bounded high-coordinate summands in the ideal test (12). Their conditional variance after n^(-1/2) normalization is O(1). The order-n contribution is precisely the variance of the displayed conditional mean. Retaining that mean as an explicit coherent current is an identified grouping obligation, not a completed positive realization. The main note makes that distinction correctly.

The selected-injection variance limit from (17) evaluates to approximately 2.163278239334616e-9, consistently with the parent calculation.

## Diagnostics and files

`check_modulated_jet_audit.py` ran successfully with 514 deterministic checks. These include Hessian eigenvalue/norm tests; a separate Fourier-frequency construction of the fourth derivative tensor; every proper center cut for several small dimensions; all six-tree proper cuts and the actual-leaf telescoping bounds; direct all-i contraction identities; and exact Gaussian conditional-variance calculations throughout representative complete Price panels. No Monte Carlo is used for the variance diagnostics.

`modulated_jet_audit_checks.json` records every result and source SHA-256 hashes, including the exact main-note revision audited. The executable and this audit supplement rather than replace the proof.

`INDEPENDENT-MODULATED-JET-AUDIT.md` gives an independent derivation and a separately bounded observation about pooling derivative hits at one fixed node. That observation is not needed for the main result and does not establish any original-clock aggregate or arbitrary signed-history claim.
