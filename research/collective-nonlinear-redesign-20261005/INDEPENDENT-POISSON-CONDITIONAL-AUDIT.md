# Independent audit: bounded Poisson positive conditional law

2026-10-05. Scope: `POISSON-POSITIVE-CONDITIONAL-KERNEL.md` only. No sampler-composition, Gibbs/mixing, cold-start, or pathwise derivative theorem is asserted here.

## Verdict

The exact-real conditional sampler and its finite positive cap are correct as stated. The acceptance identity, query-count bound, deterministic caps, W2 delta^(1/4) constant, caller-law estimate, and quadratic whole-block obstruction all pass independent verification. The numerical bound also survives, subject to the coupling-language correction and explicit implementation qualifications below.

## Checked identities and constants

1. For v=y-c, the fundamental theorem of calculus gives
   R(y)=integral_0^1 [h(c+sv)-h(c)] dot v ds.
   The Hessian interval makes 0<=p(s,v)<=1. With S density 2s and Lambda=A|v|^2/2, Lambda E p=R exactly. Thus Poisson-mark survival is exp(-R). Multiplying this by the N(-h_c,I) proposal density cancels the tangent term and gives the genuine conditional target.

2. The Gaussian integral is exactly
   E exp[-A|Z-ell|^2/2]=(1+A)^(-d/2) exp[-A e^2/(2(1+A))].
   The expected all-marks count over the unlimited sampler is E K/alpha, even though K and success are dependent within one trial: reaching a trial depends only on preceding independent trials. E K=A(d+e^2)/2. Early rejection only decreases executed original-gradient queries.

3. At t=log 2, E exp(tK)=(1-A)^(-d/2) exp[A e^2/(2(1-A))]. For 0<A<=1/2, -log(1-A)/2<=A and 1/[2(1-A)]<=1. The displayed L therefore gives an overflow probability <=delta/2. The displayed N gives no acceptance probability <=delta/2. Pre-drawing all N proposals legitimizes the union bound even though the actual sampler stops early.

4. For the target potential W, strong convexity gives |b-c|<=|grad W(c)|=e. Integration by parts gives E|Y-b|^2<=d and E|Y-b|^4<=(d+2)E|Y-b|^2<=d(d+2). Since [d(d+2)]^(1/4)<=sqrt(d+1), the fourth norm about c is <=sqrt(d+1)+e. On a cap, the capped output is c; otherwise it is the unlimited exact output. Cauchy-Schwarz therefore gives precisely W2<=[sqrt(d+1)+e]delta^(1/4). No independence between the cap event and exact output is needed.

5. At fixed source and projection, sup_z |h_a(z)-h_a'(z)|<=A|a-a'|. Synchronous Langevin contraction is a valid proof of W2(pi_a,pi_a')<=A|a-a'|, without being an implementation or mixing claim. Adding two independently certified cap errors gives the stated 2epsilon term.

6. For g(x)=Ax, a=c=0 and d=D, p=1 for every nondegenerate mark. Exactly the K=0 proposals succeed; every failure executes exactly one original gradient mark query. The expected mark-query count is therefore (1+A)^(D/2)-1. This is exp(Theta(AD)) for A<=1/2, hence exp(Theta(1/A)) at D=Theta(A^(-2)). This is an obstruction for this unmodified proposal mechanism, not a general lower bound for all gradient samplers.

## Numerical correction and qualifications

- Correct the sentence claiming that the common-Y branch of maximal coupling has the exact displacement/count laws. Its unnormalized proposal measure is min(q,q_hat), which is dominated by q. Conditional on equality, it is generally not q. Domination suffices for the displayed bounds: exact rejection probability at |v|<rho is <=A rho^2/2, and the unconditioned expected number of marks with S<s0 is <=A(d+e^2)s0^2/2 per proposal. Thus equation (8) does not need changed constants.

- Gaussian mean total variation is <=nu0/2 per proposal, and the mark coin perturbation away from the cutoffs is <= (nu+nu0)/(A s0 rho); clipping is nonexpansive relative to a true p in [0,1]. Uniform source-error guarantees over all executed adaptive queries are required. Failure probabilities of any probabilistic source-error guarantees must be added explicitly.

- If caps are computed from an inexact anchor query, use a certified e_upper>=e, for example |c+h_hat_c|+nu0. Freeze the same N and L in the comparison. The exact e in an analysis formula must not become an unavailable numerical input.

- Keep the charged-count cap convention unchanged when skipping small marks. In particular, skipped marks still count against the all-K budget. Otherwise the comparison needs an additional explanation of changed cap decisions.

- The numerical fourth-moment estimate is valid when every proposal is an exact unit Gaussian plus a mean perturbation bounded by nu0. The exact capped output has fourth moment at most that of the unlimited target, and the numerical output is dominated in fourth moment by the sum over N proposals. This proves the displayed (1+N^(1/4))[sqrt(d+1)+e+nu0] p_num^(1/4) floor. Small TV error of an arbitrary approximate Gaussian sampler alone does not bound its fourth moment or W2 error; scalar-law approximations need a separate moment or direct transport certificate. The text's separate scalar/probability-law error bill must preserve this condition.

## Caller and readset qualifications

- A genuine conditional-law port conditions on the incoming caller/anchor transcript, then uses fresh independent proposal/count/mark/coin randomness. Common random numbers may couple evaluations at different callers. Selecting a caller using that same kernel randomness does not automatically preserve the asserted conditional law; an adaptive composition would need its own argument.

- Every executed h value is obtained from an original full gradient value g(a+Pz); no Hessian or potential query is used. For an auditable source readset, explicitly retain the returned original vector g(x), its original location and source version, not only the projected h value, together with the already listed stochastic and decision records. Replay requirements in the packet correctly include changed proposal means and all affected descendants.

- Rejection decisions are discontinuous in caller/source data. The distinction between a conditional-law port and a pathwise first/curl source port is essential and correctly stated. Neither the W2 bounds nor the finite readset licenses appending the raw transcript to the observed output law.
