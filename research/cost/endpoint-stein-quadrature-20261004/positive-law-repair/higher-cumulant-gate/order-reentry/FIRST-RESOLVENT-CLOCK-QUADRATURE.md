# Positive clock quadrature also controls one resolvent derivative

2026-10-04. Analytical Hilbert-operator lemma. No tensor field is made into an oracle and no new law compiler is claimed.

The existing positive dyadic rule Q_delta0 approximates R_k=(k-L_OU)^(-1), k>=1, by a sum of P_r with weights w r^(k-1). Its ordinary Hermite-multiplier error is at most delta0:

    sup_(n>=0)|Q_delta0(r^(n+k-1))-1/(n+k)|<=delta0.       (1)

For the supplied dyadic panel Gauss rules plus final midpoint there is also an absolute constant C such that

    Q_delta0(r^m)<=C/(m+1),  m>=0.                      (2)

Indeed a panel at distance [a,2a] from one has positive weights of total mass a and nodes r<=1-a. Its entire contribution is at most a exp(-ma). Summing over dyadic a is at most C/(m+1). The terminal interval [1-h,1] has one node 1-h/2 and weight h; its contribution h exp(-mh/2) has the same bound. This argument is independent of the number of nodes inside each panel. One computable choice is C=6: for m>=1 the dyadic-panel sum is at most 2/m because x exp(-x)<=2 integral_(x/2)^x exp(-t)dt and those dyadic intervals are disjoint; the terminal term is at most 2/(e m)<1/m. Thus the full sum is at most 3/m<=6/(m+1); at m=0 its mass is one. Consequently C' below may be set to 7.

Let e_n denote the error multiplier in (1). Positivity, (2), and the exact resolvent coefficient imply

    |e_n|<=min(delta0,C'/(n+1)).                        (3)

Gaussian chaos isometry for a Hilbert-valued F gives

    ||D(R_k-Q_delta0)F||_L2(HS)^2
          =sum_(n>=1) n |e_n|^2 ||F_n||^2.

Using n|e_n|^2<=C' delta0 in (3),

    ||D(R_k-Q_delta0)F||_L2(HS)
        <=sqrt(C' delta0) ||F||_L2.                    (4)

Thus requesting delta0<=delta^2/C' supplies derivative error at most delta with only O_k(log^2(1/delta)) positive nodes. The same finite nodes approximate both the value and first resolvent derivative. This applies to a whole Hilbert/tensor-valued source field in Gaussian L2 and introduces no input/output dimension factor.

More generally a chaos Sobolev multiplier (n+1)^a, 0<=a<1, has error at most C_a delta0^(1-a), by (3). The restriction a<1 matters: a finite interior-node rule cannot approximate the twice-differentiated resolvent in uniform L2 operator norm. Its multiplier has n|e_n| tending one at arbitrarily high chaos because the quadrature side decays exponentially and 1/(n+k) does not. A C2 source may nevertheless have extra structural first information; that is a separate norm and cannot be inferred from this quadrature lemma.

Application to the connected-cumulant hierarchy: if its complete coefficient source F_k is already controlled in Gaussian Hilbert L2, then both kappa_k=R_k F_k and Dkappa_k have polylogarithmic positive-clock approximations. For B=Dv=R_2(Dg), DB uses ONE derivative of this resolvent applied to the already bounded original Hessian field Dg, so (4) applies analytically without requiring a Hessian continuity modulus. Neither Dg nor DB is therefore an allowed producer instruction. A native original-VALUE action realizing their required contractions, its ownership and product bounds, remains necessary.
