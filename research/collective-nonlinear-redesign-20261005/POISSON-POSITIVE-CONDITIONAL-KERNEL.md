# A bounded gradient-only positive conditional kernel, and its dimension obstruction

2026-10-05. This packet proves ONE conditional draw mechanism. It does not start a Gibbs sampler, a mixing argument, an outer cold-start program, or a path/Picard implementation. A coordinate-by-coordinate use makes at least D nontrivial conditional updates and cannot meet a sublinear-in-D update-count goal.

## 1. Conditional target and exact finite VALUE mechanism

Let g=grad U on R^D, U in C2 and 0<=Dg<=A I, 0<A<=1/2. Let P:R^d->R^D be a fixed public isometry and let a be a live caller in the orthogonal complement of its range. The target is the genuine block-conditional law

    pi_a(dz) proportional to exp[-|z|^2/2-U(a+Pz)] dz.

Define h_a(z)=P^T g(a+Pz). Each h VALUE is ONE original full-g VALUE followed by a known projection; no d-by-d Hessian or matrix reconstruction is queried. The block Hessian interval is 0<=Dh_a<=A I_d.

Choose a known anchor c in R^d, without claiming it is a mode. Capture h_c=h_a(c), and compute

    ell=c+h_c, e=|ell|,
    alpha0=(1+A)^(-d/2) exp[-A e^2/(2(1+A))].         (1)

The residual ell is not erased. There is no mode-finding iteration in this lemma.

Expose the incoming caller/source/anchor transcript first. The following proposal, count, mark and coin randomness is fresh conditional on that entire transcript. Common-random-number couplings used in proofs do not authorize reusing an already exposed bank as fresh randomness.

For each independent proposal:

1. Draw Z~N(0,I_d), set Y=-h_c+Z and v=Y-c=Z-ell.
2. Put Lambda=A|v|^2/2. Draw K~Poisson(Lambda).
3. For each mark draw independent S with density 2s on (0,1) (S=sqrt(V), V uniform), and an independent uniform B. Query h_a(c+Sv) and set

       p(S,v)=[h_a(c+Sv)-h_c] dot v / [A S |v|^2].   (2)

   If B<=p, reject this proposal immediately. If every mark survives, accept Y. At Lambda=0 accept directly.

The original Hessian interval gives 0<=p<=1. Every site is a real original-gradient VALUE. There is no potential-value, integral, expectation, derivative, or source-dependent covariance oracle.

## 2. Exact law and acceptance bound

The convex tangent gap is

    R(y)=U(a+Py)-U(a+Pc)-h_c dot(y-c)
        =integral_0^1 [h_a(c+sv)-h_c] dot v ds,
    0<=R(y)<=A|v|^2/2.

Its appearance is only in the proof. Poisson thinning gives

    P(accept |Y=y)=exp[-Lambda E_S p(S,v)]=exp[-R(y)].

The proposal density q(y) is proportional to exp[-|y+h_c|^2/2]. Multiplication by exp[-R(y)] cancels its linear tilt exactly, leaving exp[-|y|^2/2-U(a+Py)] times a constant. Thus the first accepted proposal is exactly pi_a.

Moreover

    alpha=E_q exp[-R(Y)]
       >= E exp[-A|Z-ell|^2/2]=alpha0.                (3)

This explicitly handles imperfect centering. For d=1 and e<=1, acceptance is bounded away from zero uniformly over A<=1/2. For arbitrary d, the bound is useful only when A(d+e^2) is controlled.

The anchor costs one original VALUE, retained and reused. The expected sum of all Poisson counts per proposal is A(d+e^2)/2, so the expected original leaf count of the unlimited exact algorithm is at most

    1 + A(d+e^2)/(2 alpha0).                          (4)

Early rejection can reduce the actual count. Proposal/scalar-randomness work costs O(d) per trial, plus O(D) for forming each original full-space query unless the embedding has cheaper known structure. These arithmetic costs are separate from original-query counts.

## 3. Deterministic finite caps and a positive error certificate

For a chosen delta in (0,1), set

    N=ceil[log(2/delta)/alpha0],
    L=ceil[{N A(d+e^2)+log(2/delta)}/log 2].           (5)

Try at most N proposals. Also keep a budget of L Poisson marks, charging all K marks of a proposed trial before evaluating any, even if early rejection would use fewer. If the next K would exceed the budget, or no proposal has been accepted by N trials, output the anchor c. Thus the output is an ordinary positive probability law with at most 1+L original VALUES. A capped Poisson categorical draw with values 0,...,L and an overflow atom is sufficient; no unbounded loop is required for count generation.

Proof of cap probability: failed acceptance trials cost at most

    (1-alpha)^N <=exp(-alpha0 N)<=delta/2.

For independent pre-drawn proposals Y_j, let K_j be their Poisson counts. At t=log 2,

    E exp(t K_j)
      =(1-A)^(-d/2)exp[A e^2/(2(1-A))]
      <=exp[A(d+e^2)].

The last inequality uses A<=1/2. Chernoff therefore gives

    P(sum_(j=1)^N K_j>L)
      <=exp[N A(d+e^2)-L log 2]<=delta/2.

Couple the capped algorithm to the unlimited exact one with the same entire proposal/mark stream. They agree unless a cap fires, an event of probability at most delta. The conditional target has Hessian at least I. If b is its mode, strong monotonicity gives |b-c|<=e. Integration by parts around b yields

    E_pi |Y-b|^2<=d,
    E_pi |Y-b|^4<=d(d+2).

The second inequality follows by integrating the divergence of |y-b|^2(y-b); the convex remainder of the potential is nonnegative in the radial direction. Consequently

    ||Y-c||_4 <=sqrt(d+1)+e.

Cauchy-Schwarz on the cap event then proves

    W2(pi_a, pi_(a,capped))
      <=[sqrt(d+1)+e] delta^(1/4).                    (6)

For target epsilon, choose delta=min(1/2,[epsilon/(sqrt(d+1)+e)]^4). Equations (1),(5) are explicit finite order/cost constants. No undefined accuracy-dependent expectation remains.

## 4. Caller law port, readset, and important derivative boundary

For two complement callers a,a', the same original Hessian bound implies

    sup_z |h_a(z)-h_(a')(z)|<=A|a-a'|.

The targets are 1-strongly log-concave. Synchronous Langevin coupling, used only as a proof, yields

    W2(pi_a,pi_(a'))<=A|a-a'|.                        (7)

Indeed the drift difference is bounded by -|Delta|+A|a-a'|, and taking its stationary limit gives (7). If both cap errors are bounded by epsilon, the capped laws have

    W2(pi_(a,capped),pi_(a',capped))
      <=A|a-a'|+2epsilon.

For a fixed caller region, choose a public upper bound e_max and use it in (1),(5),(6), so the number of slots and cap version are frozen. An anchor c(a) may be used, but its ACTUAL h_a(c(a)) and ell(a) remain live. Exact target law does not depend on the anchor; cap constants do.

The complete readset includes the full original g(a+Pc) vector and every full g(a+P(c+Sv)) vector, not only their projected h values. Retain its original full-space location, caller, source version, captured anchor, proposal Gaussian, Poisson count, mark, acceptance coin and decision. Changing a caller, source, anchor or coefficient changes proposal means and later sites; all affected descendants must be replayed. One may preallocate the finite capped random stream and use common random numbers, but acceptance boundaries remain discontinuous.

THIS IS A CONDITIONAL-LAW PORT, NOT A BOUNDED PATHWISE FIRST/CURL SOURCE PORT. Acceptance/rejection can switch the output proposal under arbitrarily small caller changes. A finite HVP sweep through one fixed decision branch does not certify a global derivative. The existing smooth near-gradient completion compiler cannot consume this kernel without a separate law-level theorem. No raw clock, decision, source record or proposal transcript is appended as an observer under (6) or (7).

## 5. Numerical floors are not free

The exact-real VALUE model above has a literal finite call cap. A numerical implementation must separately price Gaussian/probability arithmetic, Poisson tails, original VALUE errors and rare denominators in (2). Positivity-preserving clipping alone is not a proof of exact acceptance.

Here is one explicit error-accounting option, demonstrating that no inverse-denominator error is silently discarded. Work with a proposal anchor estimate of error at most nu0 and all mark source estimates of error at most nu. Choose the caps using a certified e_upper>=|c+h_c|, for example |c+h_hat_c|+nu0; using the uncorrected noisy residual is not certified. There are at most N proposals and L marks. Choose rho,s0>0. Automatically accept proposals with |v|<rho and skip marks with S<s0, declaring these approximation changes. Skipped marks still count against the original K-based cap; the numerical shortcut does not enlarge the finite readset budget. On all remaining marks, calculate (2) and clip its approximation to [0,1]. Coupling to exact proposal/coins gives an extra failure-probability allowance bounded by

    p_num <= N nu0/2
       + N A rho^2/2
       + N A(d+e^2) s0^2/2
       + L(nu+nu0)/(A s0 rho),                       (8)

in addition to separately recorded scalar/probability-law errors. The first term is the total-variation allowance for N unit-Gaussian mean perturbations. The next terms bound exact rejections skipped at small displacement or small marks. The last bounds the sum of coin-probability perturbations away from the cutoffs; clipping cannot increase distance from the true p in [0,1]. Under the proposal maximal coupling used for the first term, the common-Y branch has subprobability density min(q,q_hat), dominated by the exact q. That domination, rather than an incorrect claim of conditional equality of laws, justifies the following cutoff bounds with the original e (or its certified upper bound).

Assuming exact Gaussian proposal draws (or a separately certified fourth-moment bound for their numerical substitute), the numerical output can be one of at most N Gaussian proposals or c, so a conservative fourth-moment bound about c is N^(1/4)[sqrt(d+1)+e+nu0]. Hence an additional W2 floor of at most

    (1+N^(1/4))[sqrt(d+1)+e+nu0] p_num^(1/4)

covers these disagreements by Cauchy-Schwarz. This is intentionally conservative. Source precision and known-scalar sampling work needed to achieve the displayed parameters must be included in a complete numerical bill; they are not counted as zero work. Alternatively retain only the exact-real theorem and declare a separate implementation floor.

## 6. Source-valid obstruction for the whole-vector tangent proposal

The mechanism is a general-class whole-block positive law, but its dimension bill is genuinely bad. Take the original anchored quadratic source

    g(x)=A x, U(x)=A|x|^2/2,
    a=0, c=0, d=D.

Then every p in (2) equals one, so a proposal is accepted exactly when K=0, and

    alpha=(1+A)^(-D/2).

Every rejected proposal makes one original mark query before rejection. Thus its expected original mark count is exactly

    alpha^(-1)-1=(1+A)^(D/2)-1.                      (9)

No analytical looseness caused this cost. For the two-shell dimension scale D~A^(-2), even this elementary quadratic costs exp(Theta(1/A)) in this unmodified tangent-proposal mechanism. It cannot meet the desired dimension regime.

A d=1 conditional has controlled acceptance, but independently updating all D coordinates requires D updates and is outside the sublinear dimension goal. Blocks with d of order 1/A control acceptance at the cost of order AD block updates for a full sweep; no mixing or independence claim follows. This packet does not pursue that sampler.

The useful reusable output is therefore exact nonlinear positive acceptance from original gradients, together with a finite-cap conditional-law/caller port. To become a competitive whole-vector mechanism it needs a genuinely better collective proposal/envelope or a law-level composition that avoids both (9) and coordinate-count cost. Merely relabeling the tangent proposal as a nonlinear resummation does not achieve that.
