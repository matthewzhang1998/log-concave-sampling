# Exact K3 native curl split and a paid feedback remainder

This is a new analytical supplement to the actual finite source c8b20bf4 and source audit83914986. It isolates a smaller coefficient target without claiming an executor for that target.

## 1. Notation and literal complete derivative

All maps below act on the full master record (X,Y,Z). Let C and P denote the original ancestor and recorded terminal rows, and T the Z selector. Write

    S=PW, q0=CW,
    q1=CW+aM*[gt(S)-gt(S+epsilon Z)],
    t0=S+aM g0(q0), t1=S+aM g0(q1),
    A0=Ht(t0), A1=Ht(t1), B0=H0(q0), B1=H0(q1),
    L0=Ht(S), L1=Ht(S+epsilon Z).

Each displayed Hessian is analytical notation for an original primitive derivative at its LITERAL shared query. It is not an executed derivative oracle. All have operator norm at most one.

Direct differentiation of the finite source E=r[gt(t0)-gt(t1)] gives

    DE = r(A0-A1)P
         + r a(A0 M B0-A1 M B1)C
         - r a² A1 M B1 M*(L0-L1)P
         + r a² epsilon A1 M B1 M*L1 T.                 (1)

The final sign agrees with the separately checked positive D_ZE product. No query or source version is changed.

## 2. Full recorded curl

Put F=P*E, a square source on the master Gaussian space. For any matrix N define SkewDiff(N)=N-N*, with no factor one-half. The first term in P*DE is symmetric. Therefore

    Curl F = C_lead + C_feedback,

where

    C_lead = r a SkewDiff(P*(A0 M B0-A1 M B1)C),
    C_feedback = -r a² SkewDiff(P*A1 M B1 M*(L0-L1)P)
                  +r a² epsilon SkewDiff(P*A1 M B1 M*L1 T).  (2)

The contraction bounds give the UNIFORM complete-master estimate

    ||C_feedback||op <=(4+2epsilon) r a² <=6 r a².       (3)

This does not use small Hessian differences, a derivative of the small E VALUE, or an independence approximation. In the affine hostile-twin specialization C_lead vanishes while the feedback terms can produce the true auxiliary covariance; they are paid here at their actual higher-order allowance rather than asserted to be zero.

## 3. Whole-source one-energy orientation consequence

Let R_F=D(-L)^-1(F-EF), with the Gaussian number-operator convention. Its full matrix-HS L2 norm is at most e=||F-EF||₂, equal to the original physical centered energy because P is a coisometry. The exact forward-orientation identity is

    O_F=Cov(F)-Sym B_F=-Sym E[R_F Curl F].

Define the leading constant target O_lead=-Sym E[R_F C_lead]. Equation(3), Cauchy-Schwarz and the intact whole-source energy give

    ||O_F-O_lead||HS <=6 r a² e.                       (4)

No surrogate source energy replaces e. The ordinary matrix operator bound is also O(r² a²), using the actual O(r) row frame of R_F. With r=A, a=A^g and kappa0=ra, (4) is O(a kappa0 e), smaller than b kappa0 e whenever b=A^zeta and zeta<g.

At an independently gapped constant Gaussian covariance endpoint, (4) is directly payable by covariance transport. It does not require differentiating the error. A producer using random retained roots must still establish its own complete joint law, coefficient frame and consumer chronology BEFORE this constant-target replacement can be used.

## 4. What remains to execute

The leading word is the difference of A0 M B0 and A1 M B1, at the SAME original and feedback-shifted ancestor queries. The marked R_F remains that of the complete actual E, including both terminal/twin records. Independent primitive heating of these matrices is not licensed by (4).

In particular C_lead need not be a closed two-form by itself: closedness is guaranteed only for C_lead+C_feedback. Thus closed-curl matrix-test or root-derivative frames cannot be imported for the leading word merely from its operator bound. The supplement reduces the coefficient target and quantifies a payable remainder; it does not solve the native common-query extraction or its positive retained return.
