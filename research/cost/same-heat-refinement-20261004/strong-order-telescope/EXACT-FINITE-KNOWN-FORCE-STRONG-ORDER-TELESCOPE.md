# Exact finite strong-order telescope for known-center forces

2026-10-04. New construction relative to the reviewed quarter/short known-center source. This is an actual source coupling, not a Wasserstein coupling inferred from posterior-law accuracy.

## Result

A finite strong telescope is available for the **known-center original-force mean**. At the canonical force radius r=A, its level-j actual paired difference has uncentered Lp energy

    ||Delta_j||_p <= Lambda_(J,p) sqrt(Dim) A^j,  j>=1.

Its complete fresh first is only asserted to be O(A), and its caller first O(sqrt(A)). No small derivative is inferred from that small VALUE difference.

There are two lawful grid schedules, with different costs:

1. **Minimal matched schedule:** level j executes j quarter/short maps, each with local finite allowance A^(j+1/2). Its complete original-query cost is C_j<=Lambda_J A^(-(j-1)_+), and its posterior-law order is j+1/2. The complete CW7 known seed is overqualified at the first few levels; this is harmless. All coarse grids remain their own actual grids.
2. **Preserve the seven-order seed offset:** use allowance A^(j+7). Then level j has posterior-law order j+7 and C_j<=Lambda_J A^(-(j+11/2)). The same paired force energy A^j holds, but the extra numerical accuracy does not improve this suffix-coupling energy. It makes the telescope's mean bill unnecessarily larger.

For the first schedule, an exact-unbiased, positive buffered empirical kernel for E[H_J|theta], with normalized W2 allowance sqrt(Dim) delta, costs

    Q_mean <= Lambda_J [1 + A^(-(J-1)_+) + A^2/delta].       (R1)

This includes one or more complete pairs at every level, all deterministic zeros/anchors, and all original-query replays. Fixed target constants and public logarithms are suppressed only in Lambda_J. The full private tape can have an inverse-A number of complete seed copies and must not be described as polylogarithmic.

If H_J must itself be a known posterior source of order R=J+1/2, and the mean is inserted physically with sqrt(A), choosing delta=A^(R-1/2)=A^J gives exponent R-3/2. If only that physically inserted force mean must approximate the exact posterior force mean to A^R, the r=A readout gives an additional A in the sampler-bias estimate: J>=R-3/2 suffices, and the half-integer-order schedule gives exponent R-5/2. For arbitrary non-grid R, round J upward and retain that rounding in the exponent. These are query, not dense-arithmetic, exponents.

**This does not yet construct a linear COMPLETE growing-depth hidden family.** In particular, a marginal buffered-mean law of order delta is not an actual cross-order statistic coupling of order delta. Section 10 gives a concrete square-root gap and the missing hidden-source contract.

## 1. Frozen inputs and the actual finite levels

Fix A, the original potential, the original exposed caller theta (including the known physical center y(theta)), maximum level J, required moment list, all grids, iteration counts, original finite numerical versions, and row encodings before sampling or differentiating theta. Use one fixed complete finite known CW7 seed X_*(A,y;W_*). Its law allowance is Lambda sqrt(Dim) A^7 and its complete query cost has heat exponent zero. Its actual firsts, original zero, moment profiles, and retained-gradient/finite-restoration ports are the reviewed input contract.

The same sufficiently restored seed version is used in every level. Tighter finite arithmetic/restoration of this fixed seed is a paid polynomial-logarithmic choice under that input port, not a replacement of its law floor or a free finest deterministic grid.

Write Phi_j(x;Z,L) for ONE level-j quarter map followed by ONE protected short map, using exactly the formulas (Q) and (S) of the reviewed source. The short duration is always h=A^(19/30), independently of j. Its emitted value includes the full L dependence. No auxiliary protected projection is used in this telescope.

For a general offset b in [1/2,7], set

    R_j=j+b,
    N_q,j=max(1,ceil(Lambda_J A^(-(R_j-3/2)))),
    N_s,j=max(1,ceil(Lambda_J A^(-(R_j-17/5)))).

Choose finite M_q,j>=1 so M_q,j+3/2>=R_j, and M_s,j so

    1/2 +(34/15)(M_s,j+1)>=R_j.

Allocate finite arithmetic, original-node and original-zero restorations to the same or smaller local allowance. All occurrences of a named level use its same fixed versions. Constants/strict allocations are chosen after enumerating the complete finite construction. The meaningful real finite circuit is differentiated using its literal original HVPs; a discontinuously rounded machine implementation does not acquire those first bounds by a VALUE claim.

The actual level is

    X_0=X_*,
    X_j=Phi_j(...Phi_j(Phi_j(X_*;Z_1,L_1);Z_2,L_2)...;Z_j,L_j).

Thus all j maps in X_j use level-j grids. This intentionally differs from appending one increasingly accurate map to each previous finite level. That latter definition would preserve a cheap earlier local numerical error after too few contractions and does not have the same law declaration automatically.

Each level's actual finite grid is fixed globally. In particular, when X_(j-1) occurs as the coarse side of a pair it still uses Phi_(j-1), never Phi_j.

For definiteness the terminal force is

    H_j(theta;W_j)=(r/sqrt(A)) grad V(X_j(theta;W_j)).

A fixed common physical anchor may be subtracted from every H_j without changing Delta_j. Any caller-dependent anchor remains an actual original computation with its full caller path. It is not a full-joint constant merely because it is independent of the private record.

## 2. Same-grid contraction is a finite VALUE identity

Let K_h=1-cos(h). Positivity and exact row sums in the finite quarter recurrence imply, for M_q>=1 and the same Z,

    |Q_j(x;Z)-Q_j(x';Z)| <= A/(1-A) |x-x'|.

Indeed the direct endpoint cos(pi/2) input coefficient is exactly zero, while all intermediate differences are bounded by |x-x'|/(1-A). This argument concerns the actual finite recurrences and is independent of their quadrature accuracy.

Similarly the actual finite short recurrence, at the same L, obeys

    |S_j(u;L)-S_j(u';L)| <= 1/(1-A K_h) |u-u'|.

Consequently every Phi_j has the same safe input contraction

    lambda_A=A/[(1-A)(1-A K_h)] <= C A <1.             (C)

The exact quarter/short reference has the analogous bound and preserves the same Q_(A,y). Use exact endpoint affine coefficients or an explicitly budgeted normalized encoding; an accidentally nonzero rounded cos(pi/2) is not silently discarded in (C).

Analytical centering at the exact posterior mode m is used for estimates only. The executable recurrence stays the original unanchored one, or uses a finite mode with its exact residual retained. Seed and finite-trajectory moment profiles give

    ||Phi_j(U;Z,L)-Phi_exact(U;Z,L)||_p
       <= Lambda_(J,p) sqrt(Dim) A^(R_j)              (L)

at every actual input profile used below. The same statement includes its explicitly allocated finite numerical allowance. No modulus of the Hessian beyond the original bounded-Hessian hypotheses is used in the VALUE bound.

Applying (C), (L), and exact stationarity at each of the j steps proves

    W2(Law(X_j|theta),Q_(A,y))
      <= lambda_A^j epsilon_* + C epsilon_j
      <= Lambda_J sqrt(Dim)[A^(7+j)+A^(j+b)].         (Law)

For b in [1/2,7] this is the stated order j+b. The same moment estimates hold along the intermediate finite trajectories. Caller-uniformity is not inferred when the imported profiles are only integrated against the caller law.

## 3. The actual suffix pair and the coarse-grid bill

One paired bank at level j samples one COMPLETE W_* and j independent momentum pairs eta_i=(Z_i,L_i). It computes

    X_j^+ = Phi_j( ... Phi_j(X_*;eta_1)...;eta_j),
    X_(j-1)^- = Phi_(j-1)(...Phi_(j-1)(X_*;eta_2)...;eta_j).

For j=1 the second value is X_* itself. Define the executed difference

    Delta_j=H_j^+ - H_(j-1)^-.

The fine side has precisely the named finite level-j marginal; the coarse side precisely the named level-(j-1) marginal. Sharing the last j-1 momenta changes neither marginal. Sharing the seed is legal within this complete pair. It does not promote the seed to the external caller.

After the first fine map let d_1 be the fine-minus-seed state difference. The actual moment profiles give

    ||d_1||_p <= Lambda_(J,p) sqrt(Dim) sqrt(A).

At each remaining common momentum pair compare first at a common fine grid, then pay the grid discrepancy at the coarse state:

    ||d_i||_p <= lambda_A ||d_(i-1)||_p
                 + Lambda_(J,p) sqrt(Dim)
                     [A^(R_j)+A^(R_(j-1))].

The second term follows by comparing BOTH different finite maps to the same exact map at that same state and same momentum. Thus, for j>=2,

    ||X_j^+-X_(j-1)^-||_p
      <= Lambda_(J,p) sqrt(Dim)
          [sqrt(A) A^(j-1) + A^(j-1+b)].             (P)

This is an uncentered source bound. It is proved on the actual shared Gaussian tape. It is not obtained from W2(Law(X_j),Law(X_(j-1))).

The last coarse discretization discrepancy is not multiplied by A^(j-1): it enters at the last common map. Formula (P) retains that fact. At b=1/2 it has exactly the same order as the contracted initial discrepancy; at b=7 it is much smaller. A hypothetical common finest grid on both sides would eliminate this term by changing the coarse marginal and its price, which is not this construction.

Bounded original Hessian and the literal readout now give

    ||Delta_j||_p <= Lambda_(J,p) sqrt(Dim)
                       r[A^(j-1)+A^(j+b-3/2)].      (F)

For every b>=1/2 this is at most Lambda sqrt(Dim) r A^(j-1). At r=A it is Lambda sqrt(Dim) A^j.

## 4. Exact expectation telescope and complete independence

Because the pair's marginal versions are the named ones,

    E[Delta_j|theta]=E[H_j|theta]-E[H_(j-1)|theta].

Therefore the finite identity is exactly

    E[H_J|theta]=E[H_0|theta]+sum_(j=1)^J E[Delta_j|theta]. (T)

No limiting source, limiting grid, unknown mean query, or independent fine/coarse subtraction is used in (T).

Different empirical pairs, including pairs at different levels, are independent COMPLETE banks conditional only on the original theta. Each draws a new complete W_* as well as its own complete eta list. The old seed's internal histories, observations, hidden statistics, decoder calls, finite zeros, and aliases remain inside that copied seed. Sharing an expensive sampled seed across these banks and then conditioning on it would define a different conditional mean problem.

The base bank for H_0 is likewise independent of all pairs. The terminal numerical Gaussian keep is drawn independently after all signal banks and is unread by them. All later law comparisons discussed here concern this complete output statistic plus the original theta, not its independently exposed internal sample records.

## 5. Real zeros, anchors, and actual firsts

Let z_j(theta)=H_j(theta;0) be the actual all-private-zero computation of the NAMED finite level j. For a paired zero, both sides are their own actual zero versions, hence

    Delta_j(theta;0)=z_j-z_(j-1)

exactly. Set S_0=H_0-z_0 and S_j=Delta_j-(z_j-z_(j-1)). Then all S_j vanish on their actual zero and

    E[H_J|theta]=z_J+sum_(j=0)^J E[S_j|theta].         (TZ)

All z_j and their caller derivatives are executed and charged. They may be cached only at the identical original caller, version, and semantic key. Computing only z_J in the final displayed readout does not remove the intermediate z_j needed by a literally centered implementation; those evaluations are in the zero bill. Alternatively the uncentered implementation using H_0 and Delta_j has the same algebraic zero identity and can retain its explicit offsets.

The centered variances used in the mean proof satisfy

    ||S_j-E S_j||_2=||Delta_j-E Delta_j||_2<=||Delta_j||_2

for j>=1, independent of the size of the zero. Thus the buffered-law proof does not need an unsupported zero estimate. If an uncentered retained-deletion bound for S_j is requested, its zero allowance must also be proved. Under the reviewed seed's actual centered-state and posterior-mode profiles this can be done directly: the inequality

    |X_*(theta;0)-m| <= E[|X_*(theta;0)-X_*||theta]
                        +E[|X_*-m||theta]

supplies the seed-zero envelope; the deterministic zero recurrences and their separately restored numerical errors then give the zero version of (P). This is a direct deterministic-zero argument, not evaluation of a Gaussian Lp inequality at zero. Any caller-integrated rather than uniform input qualification remains in force.

The source's original finite first recurrences prove, on the full paired tape,

    Lip_private(H_j)<=Lambda_J r,
    Lip_private(Delta_j)<=Lambda_J r,
    caller(H_j), caller(Delta_j)<=Lambda_J r/sqrt(A).

Subtracting actual zeros does not change private firsts. Caller firsts of centered sources also include the zero derivatives and have the same scale up to fixed-J constants. The common seed columns are differentiated as shared columns within a pair. The coarse and fine derivative sweeps are both executed, with shared seed adjoints accumulated before its saved sweep.

No bound Lip(Delta_j)=O(r A^(j-1)) is claimed. Even a tiny displacement between two terminal states does not force their Hessians to be close at that power under only bounded continuous Hessian. For example, V''(x)=lambda+kappa cos(kx), with lambda>|kappa| and lambda+|kappa|<=1, has a uniformly bounded Hessian. Points separated by pi/k have a vanishing VALUE difference but an order-one Hessian difference. The resulting common-terminal-momentum derivative of a force difference can stay of order r.

## 6. One positive buffered mean and its finite law proof

Choose integers N_0,...,N_J>=1 and v>0. Execute independent complete samples S_(j,i) as above and an untouched independent Z~N(0,I_d). Return

    Y=z_J+sum_(j=0)^J (1/N_j)sum_(i=1)^N_j S_(j,i)+sqrt(v) Z. (M)

The exact mean is E[Y|theta]=E[H_J|theta] by (TZ). The reference law is N(E[H_J|theta],v I_d). There is no unknown covariance calculation, covariance subtraction, clipping, or executed unknown expectation in (M).

Let e_j=||S_j-E[S_j|theta]||_2 and ell_j=Lip_private S_j. The complete independent-bank Gaussian Riesz argument in W29 applies level by level to the single interpolation

    E[H_J|theta]+t sum_j T_j+sqrt(v)Z,
    T_j=N_j^(-1)sum_i(S_(j,i)-E[S_j|theta]).

The Riesz and derivative rows of the j-th independent bank are zero in every other bank. Thus its tensor contribution is the same as the one-bank contribution, and summing the L2 velocity bounds gives

    W2(Law(Y|theta),N(E[H_J|theta],v I_d))
       <= (1/(2sqrt(v))) sum_(j=0)^J ell_j e_j/N_j.   (B)

One common numerical keep is sufficient; separate variance shares are unnecessary. This proof uses the complete paired source as its C1 VALUE input and does not require Delta_j to be a gradient or have small curl. It uses the original W29 one-energy estimate, not an asserted variance-matching theorem.

For the canonical normalization, factor out sqrt(Dim) and public/caller envelopes. The admitted bounds are

    e_0<=Lambda sqrt(Dim) A, ell_0<=Lambda A,
    e_j<=Lambda sqrt(Dim) A^j, ell_j<=Lambda A  (j>=1).

Writing a_j=ell_j e_j/(2sqrt(v) sqrt(Dim)), safe canonical envelopes are

    a_0<=Lambda A^2,  a_j<=Lambda A^(j+1).

The target sqrt(Dim) delta is met by sum a_j/N_j<=delta. At caller-dependent envelopes use their full admitted integrated profiles when applying (B); do not replace a product of correlated path quantities by a product of means.

The emitted signal's private first is bounded by

    [sum_j ell_j^2/N_j]^(1/2),

plus the independent keep's sqrt(v) first. Its caller derivative has only the coherent sum bound O_J(r/sqrt(A)), with actual anchor terms. Neither its caller derivative nor its common numerical bias receives a 1/sqrt(N_j) gain. Finite numerical errors common to every copy are budgeted separately, or absorbed by targeting the same exact named finite H_j as in (T).

## 7. Optimal empirical allocation, including mandatory samples

Let c_j be the COMPLETE query price of one base sample or paired sample, including its prescribed original numerical version. Ignore deterministic setup for this optimization and minimize

    sum_j c_j N_j  subject to  N_j>=1, sum_j a_j/N_j<=delta.

The exact continuous constrained optimizer is

    N_j=max(1,sqrt(lambda a_j/c_j)),                  (A1)

where lambda>=0 is chosen to saturate the constraint unless all N_j=1 already suffice. On an active set I, its explicit value is

    sqrt(lambda)=sum_(j in I) sqrt(a_j c_j)
                   / [delta-sum_(j not in I) a_j].

The active set must agree with lambda a_j>c_j. Rounding each continuous N_j upward is legal and costs at most an extra sum c_j.

A simpler allocation with the same sharp heat powers is

    B=sum_j sqrt(a_j c_j),
    N_j=max(1,ceil((B/delta)sqrt(a_j/c_j))).          (A2)

It gives the completely explicit bounds

    sum a_j/N_j<=delta,
    sum c_j N_j<=B^2/delta+sum c_j.                 (A3)

Cauchy-Schwarz gives the universal lower bound B^2/delta for this empirical error certificate; N_j>=1 gives the additional lower bound sum c_j. Thus (A2) is within a factor two of the larger continuous lower bound, and (A1) is the exact continuous optimum. No exact integer-optimality claim is made for ceiling.

For b=1/2, c_0<=Lambda and c_j<=Lambda A^(-(j-1)_+). Every j>=1 then has

    a_j c_j<=Lambda A^2,

and the base has the same bound. Hence B<=Lambda_J A and

    N_0=O_J(max(1,A^2/delta)),
    N_j=O_J(max(1,A^(j+1)/delta)), j>=1,
    Q_mean<=Lambda_J[sum c_j+A^2/delta].            (A4)

This proves (R1), after adding the actual deterministic zero/setup bill, which is O_J(sum c_j).

For b=7 instead, c_j<=Lambda A^(-(j+11/2)) and a_j c_j<=Lambda A^(-9/2). Therefore

    Q_mean<=Lambda_J[A^(-(J+11/2))+A^(-9/2)/delta].  (A5)

This is still linear in the target order but has a worse fixed offset. Its extra cost cannot be erased by merely declaring that the fine and coarse final distributions are both high order.

## 8. Complete query, FIRST, arithmetic, and tape ledger

One level-j sample performs:

- one COMPLETE seed call at heat A;
- j M_q,j N_q,j original quarter-gradient calls;
- j M_s,j N_s,j original short-gradient calls;
- its terminal original gradient and every required anchor/zero/numerical call;
- its actual known coefficient operations.

One level-j pair performs BOTH its level-j and its level-(j-1) transcripts. The identical seed VALUE may be shared once within that pair. No later original gradient is reused unless its physical point, tape, version, and semantic key actually coincide. The safe bound is c_j<=C_j+C_(j-1) plus explicitly excluded setup, and this has the displayed level-j heat exponent.

All N_j complete pairs are paid. The cost is additive over levels and banks. There is no N_j^2 original-query term, but there is also no replacement of an N_j-fold complete seed bill by one shared sampled seed. Deterministic zeros at a fixed caller may be evaluated once and reused legally; a zero at a new private child center is a new complete occurrence.

The pair's private Gaussian record has

    n_pair,j=n_seed+2 j d,

with the coarse side selecting the suffix coordinates. It retains both primal finite transcripts. Deterministic grid nodes add query/storage cost, not Gaussian variables. The complete mean record has

    n_mean=N_0 n_seed +sum_(j=1)^J N_j(n_seed+2 j d)+d_keep,

before any separately specified auxiliary records. These are actual copies. The potential inverse-A growth cannot be removed from a later consumer's input dimension without a literal retained deletion and its source energy/caller/origin certificate.

Each forward or adjoint sweep uses original HVPs at every executed gradient node. Shared-node adjoints accumulate before a single saved shared sweep. Discarded tapes require paid replay. Dense known product-integration work can cost O(j M_q,j N_q,j^2+j M_s,j N_s,j^2) per sample; this report does not give that arithmetic the smaller original-query exponent. Gaussian generation, storage, coefficient normalization, and precision bit lengths also have their real counts.

The reviewed weighted normalization may be applied separately to every deterministic group and inherited seed block. It preserves original-query identities and avoids a spurious sqrt(N_q) graph-radius loss. It does not eliminate the N_j independent empirical Gaussian records or make their later retained/proxy source return automatic.

## 9. Exact target, posterior-force bias, and half-order conventions

The empirical identity (T) targets the finite E H_J exactly in expectation. If the actual objective is the exact known-center posterior force mean, bounded original Hessian gives the additional deterministic bias

    |E H_J - (r/sqrt(A)) E_(X~Q_(A,y)) grad V(X)|
       <= Lambda sqrt(Dim) r A^(R_J-1/2).

After the physical sqrt(A) insertion this is Lambda sqrt(Dim) r A^(R_J). At r=A and the minimal schedule R_J=J+1/2, it is A^(J+3/2), while buffered physical error is sqrt(A) delta.

Thus the two different choices in the Result are intentional:

- requesting an order-R posterior source as well as its force mean fixes J>=R-1/2;
- requesting only a physical force-mean allowance A^R requires J>=R-3/2 and delta<=A^(R-1/2).

Neither changes the exact finite telescope. For half-integer R the second choice makes the mandatory finest-pair and empirical terms both have exponent R-5/2. For arbitrary R use the actual chosen integer J in (R1), rather than suppressing a heat-power rounding loss. Low targets use the actual max(1,...) counts and the completed seed; negative exponents do not mean fractional oracle calls.

## 10. What is still needed for depth-R hidden centers

The construction above handles the means of the specified known-center force hierarchy. It does not show that two COMPLETE hidden-level sources at adjacent order have the same strong difference. Their private hidden-statistic programs, their empirical banks, and their decoder children are part of each source occurrence; none becomes a cheap known caller by declaration.

The elementary source

    Y_N=mu+(A/N)sum_(i=1)^N W_i+Z

already separates the relevant statements. Its marginal law is exactly N(mu,1+A^2/N), so

    W2(Law(Y_N),N(mu,1))=sqrt(1+A^2/N)-1<=A^2/(2N).

But on the literal shared-Z source coupling its deviation from mu+Z has L2 norm A/sqrt(N). This is a square-root, not a same-order, gain. Even for this linear fixture, the marginal mean theorem cannot justify deleting or strongly coupling the empirical source at its smaller 1/N law error. An alternative executable coupling exploiting its known covariance would be a different construction. No unknown covariance correction is available for a general force source here.

The same problem appears when adjacent orders use nested sample counts. Shared initial samples give the exact variance

    E|bar S_N-bar S_M|^2=e^2 |1/N-1/M|

for iid centered S and M>=N, not an e*ell/N squared bound. Full-bank independence for each marginal must still hold after any cross-order coupling is chosen.

A sufficient next theorem would have to supply an ACTUAL pair of COMPLETE hidden statistics (T_j,T_(j-1)) at the original caller with:

1. Exactly the named finite marginal versions, including each level's own sample counts, zeros, variance shares, and complete old-source bills.
2. An uncentered strong statistic gap at the level power required by the proposed hidden force telescope, with separately proved deterministic-zero gap and ordinary actual first/caller bounds.
3. Shared actual decoder observation noise and a legal actual coupling of each pair of differently accurate known-center children, retaining their different grids. If decoder lengths differ, the extra prefix, mode, and observation-covariance changes are explicitly coupled and priced.
4. Retained-state, full original-gradient/proxy/finite-restoration, active-tape, origin, caller and protected-row returns for the complete hidden output. A small endpoint law or a late protected L alone does not establish these returns.
5. A level AND hidden-depth recurrence on complete work. It must charge every lower-depth source in every paired bank. A depth proportional to R cannot be put into a fixed-order prefactor when it changes an A exponent.

Given such a statistic pair, the actual decoder center derivative is useful: with a common observation pool and matched steps, the center multiplier is L_child/2=1/2+O(A), and a geometric sum transfers the actual statistic gap plus actual child gaps. This transfer uses the finite child first, not posterior W2 contraction. One can choose a common sufficiently long decoder length if its logarithmic repetitions and changed observation covariance are actually fixed and paid, but one cannot choose a universal finest original-gradient grid for free.

The known-source pair proved here supplies the child-pair VALUE component of that possible argument. It does not supply the strong hidden-statistic pair or the complete hidden graph return. In particular, (R1) is not a recurrence for a growing-depth hidden family and is not an o(R) original-query result.

## 11. Sources and checks

Read and used:

- `protected-short-refresh/DETERMINISTIC-KNOWN-CENTER-SOURCE-WITH-PROTECTED-REFRESH.md`, SHA256 44cfd4b2d5ef98b86df8b824c1b98217281547eb7cf19c3d3ff05f095d816036.
- Its independent short-refresh and weighted-extension audits.
- `p-weak-mean/WEAK-29-ACTUAL-REMAINDER-AND-RETAINED-MEAN-RECONSTRUCTION.md`, especially the complete-bank buffered empirical proof and first/cost ledgers in Sections 5, 6, 9, 10.
- `no-copy-rank-20261004/host-observer-audit/ACTUAL-CW7-WEAK-E-HOST-SEPARATOR-AND-COMPARISON-ORDER.md`, for the original-caller/complete-statistic law boundary and actual decoder comparison.
- `protected-short-refresh/CONSUMER-DEPTH-AND-ACTUAL-REPLAY-AUDIT.md`, for the complete lower-depth replay issue.

`check_strong_order_telescope.py` writes `strong_order_telescope_checks.json`. It checks 36 finite quadratic two-grid suffix pairs and their contractions/law scalings; the actual all-zero/mean telescope; an explicit non-attenuated coarse-grid floor fixture; exact continuous constrained allocations and legal rounding; rational heat powers; and the shared-keep square-root gap. All checks pass. These are finite diagnostics in addition to the analytic proof, not an independent audit of every imported old seed or a proof of the missing hidden-family theorem.

No source file outside this new directory was changed. No external publication or push was performed.
