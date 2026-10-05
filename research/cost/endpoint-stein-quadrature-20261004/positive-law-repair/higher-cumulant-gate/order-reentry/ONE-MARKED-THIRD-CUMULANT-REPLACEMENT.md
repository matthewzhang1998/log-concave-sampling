# One-marked third-cumulant replacement on the same Gaussian tape

2026-10-04. Dimension-safe tensor estimate for joining the additive-H cubic packet to higher stationary-resolvent references. It does not identify their means or covariances.

Let F,G be D-output C1 maps of the SAME standard Gaussian tape, with private first at most L. Put E=F-G and e=||E-EE||2. Center F,G,E separately, retaining their common tape. Then

    ||kappa_3(F)-kappa_3(G)||HS <=6 L^2 e.             (1)

To prove this, first let U,V be any centered Gaussian-image vectors with private first at most L. For a matrix M with ||M||HS=1, Gaussian Poincare gives

    Var(U^T M V)
       <=E|DU^T M V+DV^T M^T U|^2
       <=2L^2(E|MV|^2+E|M^TU|^2)
       <=4L^4.

The last inequality uses Cov(U),Cov(V)<=L^2 I. Consequently the operator from scalar centered Gaussian L2 to the matrix Hilbert space,

    a -> E[a(U tensor V-E(U tensor V))],

has norm at most 2L^2. Apply it to each coordinate of the centered marked E and sum squared coordinates. This gives

    ||E[E tensor U tensor V]||HS<=2L^2 e.

The exact telescoping tensor identity is

    F_c^(tensor3)-G_c^(tensor3)
      =E_c tensor F_c tensor F_c
          +G_c tensor E_c tensor F_c
          +G_c tensor G_c tensor E_c.

Each term has the same bound, proving (1). Arbitrary slot symmetrization is a contraction. The proof differentiates the OTHER two sources only once and does not differentiate the energy-small E as a small-first source.

Conditional application is valid at any captured caller, with the actual conditional L,e profiles. Integrated L2 caller error follows by coupling at the same caller and Minkowski. Source independence is neither required nor allowed to replace the specified shared record.

For the analytical stationary-resolvent references, H=int g(X_r)dr and any fixed-depth force F_j have first O(A) on their complete conditional Gaussian path, while F_j-H has L2 profile O(A^2(|Z|+sqrt(D))). Therefore

    ||kappa_3(F_j|Z)-kappa_3(H|Z)||HS
       <=C_j A^4(|Z|+sqrt(D)).                       (2)

At a fresh standard Z the integrated bound is C_j A^4 sqrt(D). Hence a cubic packet for the SAME H conditional cumulant suffices through fourth law order for F_j, once all mean/covariance branches match F_j and the analytical full-law comparison is supplied.

There is no analogous assertion that Cov(F_j-H) is automatically A^4 sqrt(D). From its only-known first O(A) and energy O(A^2 sqrt(D)), the dimension-safe covariance bound is O(A^3 sqrt(D)). That self-covariance must be retained or corrected separately. A small VALUE mark cannot be differentiated to invent a smaller first.
