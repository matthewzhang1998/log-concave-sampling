# A finite three-C0 positive reserve for the reduced mixed covariance

2026-10-05. Constructive bounded-order return. Independent audit requested.

## Result and exact scope

Let g=grad U, g(0)=0 and 0<=Dg<=A I, with only the admitted C2 regularity. Write F=Dg analytically, R2 f=integral_0^1 r P_r f dr, B=R2 F, and C_H=2R2(B^2). The target here is exactly

    T(Z)=-2 Sym R2[C_H F](Z),

where Sym M=(M+M*)/2. There is a finite original-g-VALUE positive Gaussian-root program R_T(Z), using only the existing fixed-order native pair compiler, with

    ||W2(Law(R_T(Z)|Z),N(0,v0 I+T(Z)))||_(L2(standard Z))
        <= C delta A^3 sqrt(D)+Lambda A^5 sqrt(D)
                         +C A^6 sqrt(D)+e_absolute.

Here v0>0 is a fixed allocated variance; delta is the positive-clock Hermite multiplier error; and all actual normalized-radius/native guards below are imposed. Choosing delta<=A and fixed pair order b=5 gives the requested Lambda A^4 sqrt(D) grade. The pair/clock/numerical floors are absolute, never divided by a realized small source energy. The resulting total original VALUE count is a fixed public-log polynomial and the actual carrier-subtracted private/captured-Z first is at most Lambda A.

This realizes a signed covariance correction through an everywhere positive pushforward. No Hessian or coefficient matrix is executed as a producer, no expectation oracle is used, and no finite path grid is used. The exact reduction from the old mixed K current to this T remains an independent analytical input. This note does not claim the separate same-carrier mean or final full endpoint join.

## 1. Finite target and one-energy clock certificate

Use the same finite positive dyadic Gauss rule Q as in the admitted square-clock covariance service, with nodes in (0,1), weights beta_i>0,

    sum beta_i=1, sum beta_i r_i=1/2,
    sup_m |sum beta_i r_i^m-1/(m+1)|<=delta.

It has N_clock=O(log^2(1/delta)) nodes, with a recorded positive terminal gap. Define only analytically

    R2,Q f=sum beta_i r_i P_(r_i)f,
    B_Q=R2,Q F, C_Q=2R2,Q(B_Q^2),
    T_Q=-2 Sym R2,Q[C_Q F].

The multiplier bound implies ||R2,Q-R2||_(L2->L2)<=delta and both resolvents have norm 1/2. Their positivity gives ||B||op,||B_Q||op<=A/2 and ||C_H||op,||C_Q||op<=A^2/4. With one HS energy in each telescoping term,

    ||B_Q-B||2,HS<=delta A sqrt(D),
    ||B_Q^2-B^2||2,HS<=delta A^2 sqrt(D),
    ||C_Q-C_H||2,HS<=(3/2)delta A^2 sqrt(D),
    ||T_Q-T||2,HS<=2delta A^3 sqrt(D).                 (1)

The final line follows by multiplying the covariance error by F in operator norm, then replacing the outer resolvent. It is an integrated standard-Gaussian-Z certificate, not a deterministic-Z modulus for the unsmoothed original F.

Expanding the four clocks gives

    T_Q(Z)=-4 sum_j w_j Sym E[F(U_j) F(V_j) F(Y_j)|Z], (2)

where a node j=(q,s,r,t) has w_j=beta_q beta_s beta_r beta_t q s r t>0, and

    Y=qZ+c_q G0, X=sY+c_s G1,
    U=rX+c_r G2, V=tX+c_t G3,
    c_a=sqrt(1-a^2).

The original G0,G1,G2,G3 in this display are analytical query genealogy, not extra executed path inputs. In particular Y is genuinely shared with both leaves and X is shared by the leaves. There are N=N_q N_s N_r N_t nodes and sum_j w_j=1/16. Nothing interchanges noncommuting Hessian factors.

## 2. A known coherent shield split, including the F(Y) occurrence

Conditional on Z, the three original sites (Y,U,V) have known means

    mZ=(q,rsq,tsq) tensor Z

and a scalar 3-by-3 covariance K=L L*, where

    L=[[c_q,0,0,0],
       [rs c_q,r c_s,c_r,0],
       [ts c_q,t c_s,0,c_t]].

Let h=min(c_q^2,c_r^2,c_t^2) and sigma^2=h/10. Dropping the independent shared-X innovation gives K>=L0 L0*, with

    L0=[[c_q,0,0],[rs c_q,c_r,0],[ts c_q,0,c_t]].

The inverse has squared Frobenius norm

    ||L0^-1||HS^2
       =c_q^-2+(1+r^2s^2)c_r^-2+(1+t^2s^2)c_t^-2
       <=5/h.

Consequently K>=h I/5=2 sigma^2 I, so K-sigma^2 I is a known positive definite matrix. Execute one owned standard W_j in R^(3D), form the three coarse centers

    (a_Y,a_U,a_V)=mZ+[(K-sigma^2 I)^(1/2) tensor I_D]W_j.

The square root is of a known 3-by-3 scalar matrix, not a source-dependent covariance. Its versions/precision are frozen. The full coarse injection norm is bounded by a numerical constant; each site has variance at most one. Adding independent sigma-Gaussian shields to these three centers restores exactly the joint law in (2). This split does not enlarge any original heat or decorrelate original sites.

The dyadic panel weight estimates used in the square-clock service give

    sum_j w_j/sigma_j <= C,
    sum_j w_j/sigma_j^2 <= C(1+log(1/h_min)).           (3)

For the first inequality use 1/sigma<=sqrt(10)(c_q^-1+c_r^-1+c_t^-1), then factor the positive four-clock sum. Each one-clock sum of beta_a a/c_a is bounded: a dyadic endpoint panel of gap h has weight h and 1/c=O(h^-1/2). The terminal midpoint contributes O(sqrt(h_min)). The second inequality follows identically from 1/sigma^2<=10(c_q^-2+c_r^-2+c_t^-2), with O(1) per panel. Only the first bound is needed for the substantive reserve error and final first below.

## 3. Literal VALUE sources and the exact selected coefficients

For i in {Y,U,V}, freeze a_i and the known sigma and define

    f_i(u)=[g(a_i+sigma u)-g(a_i)]/(A sigma).          (4)

Every f_i is a genuine square gradient with actual private first<=1, f_i(0)=0 by same-query reuse, ||f_i||_Lp<=C_p sqrt(D), and captured-center first<=2/sigma. The origin g(a_i) is one real original VALUE call, captured separately for each vertex/node. It is reused only at that exact complete key. A changes only as a declared public bound; it is frozen in every derivative.

Use the unfiltered native C0 source

    S_i(x,z)=(sqrt(2)/s0) f_i((x+z)/sqrt(2)),          (5)

where s0 is the fixed selected-pair convention. Both partial Jacobians are symmetric; its complete private first is at most kappa0=sqrt(2)/s0 and its captured-center first at most C/sigma. There are no external probes, no separate first-chaos calibration, and no derivative source leaf.

The selected x-Jacobian mean satisfies exactly

    s0 E D_x S_i=H_i(a_i),
    H_i(a_i)=E_N F(a_i+sigma N)/A.                    (6)

Thus H_i is symmetric PSD with op norm<=1 and HS norm<=sqrt(D). Formula (6) is analytical only. Its actual implementation is (4)-(5) substituted throughout the finite pair VALUE graph.

At any fixed center, gradient symmetry, one Gaussian integration by parts and first-chaos Bessel give

    ||D_a H_i[a']||HS<=|a'|/sigma.                   (7)

Indeed the derivative matrix equals (A sigma)^-1 E[(F(a+sigma N)a')N*]. This identity uses only a smoothed derivative of F, justified distributionally and by C2 approximation, and never calls an original third derivative. The known complete coarse/retained-Z injections must be included when using (7).

## 4. Three calls, correct orientation, and positive readout

Allocate v_path=v0/2 and one untouched final variance v_free=v0/2. Give each of the N nodes path variance nu=v_path/N and put

    c=a=sqrt(nu/2), rho_Y=rho_V=A,
    rho_U=-4 w_j A/nu.                              (8)

With a fresh standard passive p_j, execute three independent COMPLETE order-b pair banks chronologically, all retaining exactly Z and their node's coarse W_j:

    y1=P_b(rho_Y S_Y; p_j),
    y2=P_b(rho_V S_V; y1),
    y3=P_b(rho_U S_U; y2),
    R_j=c y3+a p_j.                                 (9)

The y1,y2 values have no observer other than their immediate distinguished-input consumer. Keep p_j until the readout in (9). The full path, not individual unqualified marginals, is replaced using the native chronological-envelope theorem. Conditional on the coarse centers, the ideal path has three standard Gaussian channels and its root-to-passive cross matrix is

    rho_U rho_V rho_Y H_U H_V H_Y
       =-4 w_j A^3 H_U H_V H_Y/nu.                  (10)

In particular the middle factor is H_V and the unsmoothed original site's factor is at the leaf; reversing an isolated nonsymmetric product would be wrong. The scalar root sign is an actual source amplitude, never a negative probability.

The ideal readout (9), conditional on the coarse roots, is exactly a centered Gaussian with covariance

    nu I-4w_j A^3 Sym(H_U H_V H_Y).                  (11)

This is an exact Gaussian-channel identity, not a covariance truncation. There is no omitted degree-six covariance root feedback: the channel innovations already preserve the standard marginals exactly. All native gaps are checked before execution at the ACTUAL radii kappa0 |rho_i|, with the imported public-log factors, not just at A.

Return the positive finite program

    R_T(Z)=sum_j R_j+sqrt(v_free) z_free,              (12)

with independent complete node banks conditional on Z, and with z_free independent of all of them. Its source-zero carrier is the literal

    G_T=sum_j[c G_(P_b,U,j)+a p_j]+sqrt(v_free)z_free.

Here G_(P_b,U,j) is the admitted pair's recorded affine Gaussian source-zero carrier. It is independent of its incoming passive and of the source labels; its row has unit covariance by the native variance allocation. Therefore G_T is exactly N(0,v0 I). The child carriers are still present in the full tape, even though they have no source-zero path to the output. No Gaussian introduced by a law coupling is substituted for an executed carrier.

## 5. Prior errors, owned-coarse mixture, and one energy

The actual complete normalized source radii are r_i=kappa0 |rho_i|, with the fixed source/block constants enlarged as prescribed by LOW30. Its chronological pointwise-envelope composition gives a node error bounded by

    c Lambda sqrt(D) [r_U^b+r_U r_V^b+r_U r_V r_Y^b]. (13)

All inverse nu and sum-over-nodes factors are recorded public-log factors. At fixed b=5 and the actual small-radius guards, the sum of (13) is at most Lambda A^5 sqrt(D). A larger fixed b can place this separate prior floor at any higher fixed grade. There is no force-dependent inverse-energy normalization.

After these COMPLETE comparisons, condition on all W_j. The sum in (12) is exactly the Gaussian mixture with covariance

    S(W;Z)=v0 I+Delta(W;Z),
    Delta=-4A^3 sum_j w_j Sym M_j,
    M_j=H_U H_V H_Y.                                (14)

By the coherent shield split and independence of the three coefficient means within a fixed coarse record,

    E_W Delta=T_Q(Z).

This means integration of the whole original conditional product in (2); it is not strong sampling of H_i or M_j. In particular it does not split a shared original Y or X into independent centers.

We now pay the actual coarse first. From (7), bounded center injections, the product rule, and op(H_i)<=1,

    Lip_(W_j->HS) M_j<=C/sigma_j,
    ||M_j||HS<=sqrt(D).

The full independent W bank therefore obeys

    Lip_(W->HS) Delta
      <=C A^3 (sum_j w_j^2/sigma_j^2)^(1/2)
      <=C A^3 sum_j w_j/sigma_j<=C A^3,
    ||Delta-E Delta||_(L2(W);HS)<=C A^3 sqrt(D),
    ||Delta||op<=4 A^3 sum_j w_j=A^3/4.              (15)

Thus the potentially large inverse shield has been paid by its actual positive clock weight, not suppressed. Impose A^3/4<=v0/4. To apply the previously proved PSD covariance-mixture lemma to this SIGNED correction, write

    S=(v0/2)I+C'(W), C'=v0/2 I+Delta>=v0/4 I.

The centered energy and Lipschitz constant of C' are exactly those of Delta. The fixed-gap one-energy mixture theorem gives, uniformly in captured Z,

    W2(Law(N(0,S(W;Z))|Z),N(0,v0 I+T_Q(Z)))
         <=C_v0 A^6 sqrt(D).                        (16)

This uses one matrix HS energy and one Gaussian-to-HS first, never two dimension-sized vector energies. Alternatively node-by-node comparison gives a weighted w_j^2/sigma_j bill with public-log inverse-share factors; (14)-(16) avoid those unnecessary denominators by retaining the fixed global buffer.

Both v0 I+T_Q and v0 I+T have a fixed positive gap. The Gaussian-root Lipschitz bound and (1) cost at most C_v0 delta A^3 sqrt(D). Combining this with (13) and (16) proves the stated result. Only the last clock-restoration step is integrated over standard Z; the finite-target/prior/mixture service is uniform in captured Z.

## 6. Actual private/caller first and absolute VALUE precision

The source-path theorem, rather than the output-law bound, supplies each pair's actual residual first Lambda r_i, its captured-center first Lambda |rho_i|/sigma, and its absolute path norm to original g VALUES Lambda |rho_i|/(A sigma). The pair's carrier is independent of its incoming passive. Hence a nonroot carrier reaches y3 only through the residual first of every strict ancestor.

After the readout c, a safe per-node residual/private/coarse/captured-Z bound is

    Lambda c |rho_U|(1+sigma^-1)(1+A+A^2),            (17)

including the root's own tape, both descendant carriers/tapes, all three center paths, and the leaf passive residual path. The direct a p_j belongs to the known carrier and is not charged as a residual. Summing (17), using (3) and nu=v0/(2N), gives

    Lip_(complete private tape,Z)(R_T-G_T)
         <=Lambda A sum_j [w_j/sqrt(nu)](1+sigma_j^-1)
         <=Lambda A.                               (18)

The factor nu^-1/2=O(sqrt(N)) is explicitly part of the public-log majorant. Shared Z paths are added; owned independent banks are not falsely identified. An earlier physical caller receives its literal recorded chain-rule row through Z and any original-g caller dependence. No derivative of a W2 statement is taken.

A safe total absolute path norm from uniform original-g VALUE errors to (12) is

    V_g<=Lambda sum_j [w_j/(sqrt(nu) sigma_j)]
                              (1+A+A^2)<=Lambda.     (19)

An original VALUE floor eps_g<=eps_abs/[4(1+V_g)] therefore prices the total VALUE substitution error. Record and charge every individual 1/(A sigma), coefficient, native response width, origin subtraction and scalar row; (19) is their executed downstream sum, not permission to erase individual precision requirements. Allocate arithmetic/known-covariance-root/weight errors through the complete frozen absolute graph majorant and their Gaussian moment profiles. The finite endpoint gap is a known fixed power of the requested accuracy: inverse sigma changes bit precision logarithmically and does not create a replication count. Origin and zero-input paired evaluations reuse the identical primitive record to preserve literal zero. Original/conditional finite-mode work and any external caller-only anchors are additive, as in the imported service; they are not hidden inside (19).

## 7. Complete original-VALUE and Gaussian-tape count

Let Q_P,b(L;r) be LOW30's LITERAL fully expanded count of raw-source VALUE invocations in one order-b selected-pair call, including every nested mean/pair call, finite filter, signed response, stage, origin, numerical-restoration replay and copied ancestor. The public log L at each call includes its actual radius r, nu, sigma, source normalization 1/(A sigma), all coefficient/response widths, active dimensions, desired absolute precision, and declared caller scales. This is the enumerated count in LOW30's serial-query recurrence, not a unit-cost abstract action.

For one packet, cache the three original anchors g(a_i) once after W_j,Z are fixed. Every invocation of (5) then uses ONE shifted original-g VALUE plus that recorded anchor. An explicit safe primal bound is

    Q_packet <=3+Q_P,b(L_U;r_U)+Q_P,b(L_V;r_V)
                                      +Q_P,b(L_Y;r_Y).            (20)

Without the allowed identical-anchor cache, replace the three Q terms by twice themselves. Executions at zero can reuse their anchor; (20) does not rely on that further saving. Captured origins inside the pair implementation are already included in Q_P,b and cannot be silently removed. Original mode/caller-only/numerical external source work remains additive. All changes of a source argument replay the whole underlying pair source graph; the three pair calls are serial, and (20) counts the children before their consumers.

For clarity the expanded three-depth recurrence is

    Q_y1=Q_P,b(L_Y;r_Y),
    Q_y2=Q_y1+Q_P,b(L_V;r_V),
    Q_y3=Q_y2+Q_P,b(L_U;r_U),
    Q_g,total<=3N+sum_j Q_y3,j+Q_external.             (21)

The child y1 or y2 is a recorded vector input, not a function provider replayed at each parent source evaluation. In contrast EVERY internal source invocation of a pair is expanded and counted by Q_P,b. This distinction is the literal program in (9).

The imported finite enumeration gives, for fixed b,

    Q_P,b(L;r)<=C_b(1+L)^(K_b),
    Q_g,total<=C_b N(1+L)^(K_b)+Q_external
              =O_b((1+L)^(8+K_b))+Q_external.        (22)

The polynomial exponent K_b is the fixed enumerated LOW30 compiler exponent; no growing-order or inverse-A count is claimed. This explicitly displays the four-clock factor and all three complete pair coefficients.

Let d_P,b(L;r) be the actual private Gaussian dimension of each fully expanded pair call, excluding its incoming passive and external retained labels. The exact declared root-dimension upper bill is

    d_total =D(4N+1)+sum_j [d_P,b(L_U;r_U)
                          +d_P,b(L_V;r_V)+d_P,b(L_Y;r_Y)].         (23)

Here 3DN is the owned coarse W bank, DN the leaf passives, and D the untouched global buffer. Source shields and native filter/mean/pair fills are all inside the corresponding d_P terms; no additional analytical query-tree Gaussians are sampled. Thus d_total/D is a fixed public-log polynomial. Discarded primal records incur the full replay bill. A requested first or adjoint in one direction uses at most a fixed multiple of (21) original HVPs at the recorded original-g VALUE sites, through all original source ancestors. No HVP is a primal/producer leaf and no saved HVP is differentiated.

## 8. Exact quadratic check and source pins

For g(x)=H x with fixed symmetric PSD H, (4) is H u/A for every center and every sigma. Thus H_i=H/A exactly, all coarse randomness disappears from (14), and

    B=H/2, C_H=H^2/4, T=T_Q=-H^3/4.

The ideal complete reserve therefore has EXACT covariance v0 I-H^3/4, without any degree-six covariance correction. Only the separately declared finite pair/numerical approximation error remains. This is a full noncommuting path calibration at general nodes and a full matrix, not scalar-only, quadratic check.

Pinned source: LOW30 `research-source/High Acc Ideas/ai-bucket/30_low_acc.tex`, SHA256 `7767ac389d6b677eb9517b2585f7a24cb9f85ede47d099a653c635b7be750fb8`, especially `b27:compiler:pair`, `b27:nc:eq:envelope`, `b27:nc:eq:composition`, `b27:nc:eq:vertex-source`, `b27:nc:tree-replacement`, `b27:compiler:paths`, and `b27:compiler:serial`. These actual interfaces were read for this construction.

Other inputs:

- `../../clock-compression/POSITIVE-SQUARE-CLOCK-COVARIANCE-SERVICE.md`: positive clock rule and original-VALUE/endpoint conventions.
- `../GAUSSIAN-COVARIANCE-MIXTURE-ONE-ENERGY-LAW.md`: fixed-gap PSD mixture theorem; the signed shift is proved in Section 5 here.
- `../../order-reentry/RECTANGULAR-FIRST-COEFFICIENT-MEAN-AND-GRAM-RETURN.md`: complete-bank/owned-coarse and first-versus-law rules.
- `../../rank3-continuation/EXACT-FIVE-CLOCK-NATIVE-SOURCE-ADAPTER.md` and its independent native audit: native C0 normalization, positive shares and strict-ancestor path conventions.

The only new consumer step is the coherent three-site shield split followed by a three-C0 chronological chain and a signed covariance-mixture shift. There is no remaining constructive gate inside this reduced target once the stated imported fixed-order compiler/clock guards and ordinary absolute restoration floors are met. Independent review of this specialization is still required before calling it admitted.
