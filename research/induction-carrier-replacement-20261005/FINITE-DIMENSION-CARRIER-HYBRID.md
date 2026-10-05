# A proved carrier-compatible native hybrid by proof-only polynomial taming

2026-10-05. This closes the ACTUAL shared-carrier replacement at a deliberately explicit, dimension-dependent smallness guard. It does not claim the desired public-log-only dimension dependence. No executed source is changed, no conditional mean oracle is added, and no Gaussian is subtracted from an abstract LAW. The extra comparator below is used only in a coupling proof.

## 1. Exact imported input and conclusion

Keep the actual shared Gaussian geometry of COMMON-CARRIER-AMALGAMATION.md:

    P_h=L_h^T H/v+Pi_h V_h,   Pi_h=I-L_h^T L_h/v,
    H~N(0,v I_D),   L_h L_h^T=v I_D,
    X_h^a=H+u_h^a(H,Z_h),   X_h^p=H+u_h^p(H,Z_h).

Here Z_h contains the null and private tapes. Conditional on the retained external labels E, the base Z_h banks are independent across h AND independent of H; the stated Gaussian construction has exactly this product structure. The retained labels include B,Y as required; every assertion below is conditional on them with the stated uniform envelope, or integrated with its justified product-integrability envelope. The old component reads E but no new H,Z tape. No outside observer reads any individual packet's tape.

Import for fixed finite native orders:

(i) u_h^a is globally H-Lipschitz with constant s_h/sqrt(v), and ||u_h^a||_2 <= s_h sqrt(D).
(ii) u_h^p is polynomial of degree at most m in x=H/sqrt(v), with arbitrary measurable Z_h coefficients. Its combined value/operator-first envelope obeys [(||u_h^p||_12/sqrt(D))^12 + ||D_x u_h^p||_(L12;op)^12]^(1/12) <= s_h. Enlarging the original two native constants by 2^(1/12) suffices.
(iii) The completed, ordinary marginal native comparison BEFORE the group’s final untouched keep, with its external labels retained, gives W2(Law(X_h^a),Law(X_h^p)) <= e_h. This includes the actual-prior floor and strong ideal-forest polynomial remainder. It need not retain H.
(iv) s_h>=0 is a deterministic envelope, proportional to the actual root amplitude; at root zero both residuals vanish pointwise.

LOW30's raw-path theorem and the bounded radial-pullback C_k source give (i) through the complete graph paths. Its polynomial lemma gives (ii), after the Gaussian coisometry above. Its ordinary tree replacement plus strong polynomial comparison gives (iii); their separate floors are summed, not silently upgraded to a joint comparison. For the pre-keep qualification: LOW30 tree replacement couples the distinguished root spine conditional on its publics and side inputs before the final readout/buffer; its internal completion Gaussians remain part of the spine. The strong forest-to-polynomial comparison uses an identical final keep, which cancels from that strong difference. Taking the triangle of these pre-keep comparisons gives (iii); no post-keep W2 estimate is deconvolved. The finite degree m and the full p=12 constants must be taken from the actual compiled version. Enlarge s_h to contain all these constants. Set S=sum s_h and E0=sum e_h.

For any K>=1 and R=sqrt(D)+ell with ell>=0, define

    delta = K^(-5) + exp(-5 ell^2/24)
              + 2K(1+D+R^2)^(m/2) exp(-ell^2/4),          (1)
    L = max{1/sqrt(v),
            2K(1+D+R^2)^((m-1)/2)/sqrt(v)}. (2)

For m=0 omit the second term in the maximum. Then the ACTUAL original shared-carrier group and the original ideal polynomial group obey

    W2(group_actual,group_polynomial)
      <= (1+LS) E0 + (2+2LS) delta S sqrt(D)
                      + 2 L S^2 sqrt(D).                 (3)

Both groups may be translated by the unchanged old component and have the same untouched Gaussian keep appended. Equation (3) is a genuine environment-stable native hybrid. In particular, the missing global-Lipschitz condition is proved for every mixed comparison stage; it is not inferred from Gaussian Lp firsts.

## 2. Explicit polynomial taming lemma

Write x=H/sqrt(v). Expand the finite polynomial in the orthonormal multivariate Hermite basis phi_alpha(x)=He_alpha(x)/sqrt(alpha!):

    u_h^p(sqrt(v)x,Z)=sum_(|alpha|<=m) c_alpha(Z) phi_alpha(x),
    A(Z)^2=sum_alpha |c_alpha(Z)|^2=E_x |u_h^p|^2.

The coefficients and conditional energies below are analytic objects, not executed source/tensor or derivative oracles. Put

    A0(Z)=A(Z),
    A1(Z)^2=E_x ||D_x u_h^p||op^2,
    Astar(Z)=max(A0(Z)/sqrt(D),A1(Z)).

Jensen and input (ii) give ||Astar||_12 <= s_h. Let

    M_h=K s_h sqrt(D),
    q_h(Z)=min(1,K s_h/Astar(Z)),
    pi_R(x)=x min(1,R/|x|),
    u_h^t(H,Z)=q_h(Z)u_h^p(sqrt(v)pi_R(H/sqrt(v)),Z).

Use q=1 when Astar=0; if s_h=0 all residuals are zero. The Euclidean ball projection is 1-Lipschitz. This comparator is a measurable positive pushforward with the ORIGINAL Gaussian inputs. No differentiability in Z is needed for its use in the marginal hybrid. It is never fed into the ideal current theorem.

Here are dimension-explicit kernel bounds, valid for |x|<=R:

    sum_(|alpha|<=j) phi_alpha(x)^2 <=4(1+D+R^2)^j,
    sum_(|alpha|<=m) |grad phi_alpha(x)|^2
       <=4(D+m)(1+D+R^2)^(m-1).                          (4)

For the first bound use Mehler's positive generating series with rho=(1+D+R^2)^(-1):

    sum_alpha rho^|alpha| phi_alpha(x)^2
       =(1-rho^2)^(-D/2) exp(rho|x|^2/(1+rho)) <4.

The truncated sum is at most rho^(-j) times this series. For the second bound, grad_i phi_alpha=sqrt(alpha_i)phi_(alpha-e_i); summing gives

    sum_(|beta|<=m-1)(D+|beta|)phi_beta(x)^2.

The sharper operator-valued bound uses the Jacobian polynomial J(x)=D_x u_h^p, whose degree is at most m-1. If J=sum_beta C_beta phi_beta, matrix-valued Parseval gives

    sum_beta C_beta^T C_beta=E_x J(x)^T J(x).

For any unit input v, Cauchy-Schwarz gives

    |J(x)v|^2 <= [sum_beta phi_beta(x)^2]
                          v^T[sum_beta C_beta^T C_beta]v.

Consequently ||J(x)||op <=2(1+D+R^2)^((m-1)/2) A1(Z) on the ball. Since q_h A1<=K s_h, composition with the 1-Lipschitz projection proves the GLOBAL bound

    Lip_H(u_h^t) <= L s_h.                                (5)

This is deliberately conservative: the Hermite evaluation kernel pays a polynomial dimension factor. The matrix Parseval step avoids a further Hilbert/dimension loss in the derivative. It makes no claim that such a factor is necessary.

## 3. Paid strong bias, including the cap coefficient tail

The coefficient cap costs

    ||u_h^p-q_h u_h^p||_2^2
       <= E[A0^2 1_(Astar>K s_h)]
       <= D E[Astar^2 1_(Astar>K s_h)]
       <= D(K s_h)^(-10) E[Astar^12]
       <= s_h^2 D K^(-10).

For a standard D-Gaussian x, P(|x|>sqrt(D)+ell)<=exp(-ell^2/2). On that event, Holder with p=12 bounds the original polynomial part by

    ||u_h^p 1_tail||_2 <= s_h sqrt(D) exp(-5 ell^2/24).

The projected polynomial has pointwise norm at most

    2 M_h(1+D+R^2)^(m/2).

Its tail L2 contribution is therefore at most the last term of (1) times s_h sqrt(D). Combining these estimates gives

    d_h:=||u_h^p-u_h^t||_2 <= delta s_h sqrt(D),
    ||u_h^t||_2 <= s_h sqrt(D)+d_h.                       (6)

The same-root coupling keeps H identical pointwise on both sides; the displayed L2 bias is integrated over H and is not a uniform conditional-in-H estimate. The taming is NOT claimed to preserve root-linear currents, caller firsts or curl. We pay its whole strong bias twice around the native replacement and restore the original polynomial group before using any current identity.

## 4. One-packet coupling against every mixed environment

At an arbitrary hybrid stage the unchanged remainder contains actual or tamed packets only. It therefore has the deterministic GLOBAL carrier-Lipschitz bound

    Lip_H R_-h <= L S_-h,   S_-h=sum_(g!=h) s_g.

Triangle inequality supplies an ordinary marginal coupling of X_h^a=H_a+u_h^a and X_h^t=H_t+u_h^t with error epsilon_h<=e_h+d_h. Disintegrate the carrier conditional on the packet output on each side. This retains each pair's correct joint law; it does NOT assert H_a=H_t. Since both carriers have covariance vI,

    ||H_a-H_t||_2
      <= epsilon_h+||u_h^a||_2+||u_h^t||_2
      <= epsilon_h+2s_h sqrt(D)+d_h.

Draw all unchanged packet tapes Z_-h independently of this coupling, and use the SAME Z_-h on both sides. Their independence conditional on the carrier is precisely the executed read-set rule. Each side then has the correct full hybrid law. The carrier displacement changes the rest by at most L S_-h times that displacement. Thus

    W2(hybrid_before,hybrid_after)
      <=(1+L S_-h)(e_h+d_h)
           +L S_-h[2s_h sqrt(D)+d_h].                    (7)

Sum (7), using sum_h s_h S_-h <= S^2. This compares all actual packets to all tamed packets at cost

    (1+LS)E0 +(1+2LS)sum d_h +2 L S^2 sqrt(D).

Finally restore all ORIGINAL polynomial residuals on the same common H and Z bank, at cost sum d_h. Equation (3) follows. In particular no intermediate polynomial residual with an unproved global first appears in the marginal telescope.

## 5. Actual triple-cubic admission and the dimension-dependent guard

Use the original full-v common-carrier execution, with eight nonroots alpha^(1/16) and roots c_h b_h alpha^(17/2). All fixed normalization and p=12/path constants are included in s_h, so

    S <= B alpha^(17/2),    B=sum_h Gamma_h b_h.

The actual source list, original-bank dependencies and force-hit labels are unchanged. Choose K=alpha^(-1/2), and choose the smallest integer ell>=1 for which the sum of the last two terms of (1) is <=alpha^(5/2). Such an ell exists and has

    ell^2=O_m(log(2+D)+log(1/alpha)),
    delta<=2alpha^(5/2).

The finite radial threshold is explicit: it is found by evaluating inequality (1), not by appealing to an unspecified tail event. Define

    C_dim=2(1+D+R^2)^((m-1)/2)/sqrt(v),

so L=max(v^(-1/2),alpha^(-1/2)C_dim). Then the extra hybrid row is at most

    2 max(v^(-1/2),alpha^(-1/2)C_dim)
                  B^2 alpha^17 sqrt(D).                 (8)

The original ideal common-carrier current row remains C_* B^2 alpha^17 sqrt(D). The new strong taming bias is <=4(1+LS)B alpha^11 sqrt(D), and choose E0<=c alpha^11 sqrt(D) before this comparison. All these are beyond alpha^10 PROVIDED their literal finite guards pass; an explicit sufficient new guard is

    LS<=1,
    2 L B^2 alpha^7 <= c_hybrid.                         (9)

For the alpha-dependent part of L, the latter is

    2 C_dim B^2 alpha^(13/2) <= c_hybrid.                (10)

For each FIXED finite D, fixed graph/current order, fixed m and fixed v>0, with the stated public-log envelope B, these hold for sufficiently small alpha. This gives an ACTUAL positive finite original-VALUE triple-cubic group carrying the exact finite alpha^9 target with error Lambda alpha^10 sqrt(D), together with the already admitted old-bank alpha^10 join and paid cubature/old-target debts.

The crucial qualification: C_dim grows roughly as D^((m-1)/2), up to logarithmic radius factors. Equation (10) is not a public-log-D condition. This theorem therefore does not establish the desired dimension-uniform/high-dimensional admission window. It supplies an actual finite-D closure of the precise native-hybrid gap, not the stronger claim previously left conditional.

## 6. Execution, costs, curl, read sets and limitations

No taming operation or Hermite coefficient is executed. The sampler is exactly the original common-carrier finite VALUE program with its actual known source-zero Gaussian arithmetic. Consequently:

- Executed VALUE count, native repetitions, root count, memory, scalar arithmetic, covariance gaps and complete replay are unchanged from that program. Tighter e_h floors have their real fixed-order native/logarithmic cost, plus any original-oracle precision cost.
- The whole shared H and all null/private tapes are owned by the group. No separate packet value, carrier or internal source tape can be exported. Every old/external read must still be in the original coefficient record and full observer ledger.
- Actual first/caller bounds are the original complete raw-path sums, with every 1/t caller injection and numerical Sobolev floor. The proof-only comparator's large L does not become an executed terminal first.
- Private curl remains exact for the mathematical radial-pullback VALUE sources. Terminal curl still follows only from the actual literal-gradient baseline and the full executed residual first. It is not inferred from (3), and the nonsmooth proof-only projection is never part of that terminal source.
- No finite mean, old target heat debt, cubature discrepancy, source bias, calibration/Sobolev error or replay is erased by the coupling.

The missing PUBLIC-LOG dimension port is now sharply separated: replace the dimension-expensive polynomial taming by a globally carrier-Lipschitz approximation with dimension-sharp bias and only public-log enlargement, or prove a different environment-stable replacement. Ordinary Gaussian Lp firsts alone were not used to assert that stronger statement.
