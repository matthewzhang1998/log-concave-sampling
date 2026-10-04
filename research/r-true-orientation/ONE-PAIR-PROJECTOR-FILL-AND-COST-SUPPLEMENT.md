# Explicit cheap Gaussian fills for the projector one-pair fork

New implementation supplement to held addendum0b3e0dcd. Assume D=D*=D^2 is the known projector used there. For V=[D;I]/sqrt(2) and U=[I,-D]/sqrt(2), no dense matrix square root is needed.

For independent standard Z1,Z2, a square root of I-VV* acts as

    Z_in,1=(I-D)Z1 + D(Z1-Z2)/2,
    Z_in,2=(I-D)Z2/sqrt(2) + D(Z2-Z1)/2.

On the range of D this is the orthogonal projection onto the antisymmetric pair; on its orthogonal complement it is diag(I,I/sqrt(2)). The two spaces are orthogonal, so its covariance is exactly I-VV*.

The output fill is simply

    (I-D) Z_out/sqrt(2),

because I-UU*=(I-D)/2. Both formulas retain the literal shared Gaussian labels inside the input fill. Independent coordinate draws are not substituted for its correlated range-D components.

For B=I+aD, B^-1=I-aD/(1+a). Applying B, B^-1 or the two fills therefore costs only a fixed number of applications of the actual known D, plus vector arithmetic and Gaussian generation. For the survived rank-one counterexample D=e1e1*, this is explicitly linear in the output dimension, with no unknown eigenbasis or coefficient matrix. A general supplied projector still carries its actual D-application/storage cost.

The pulled-back original VALUE is h_B(x)=B*g(Bx). Its anchor h_B(cX) can be cached at the same clock/root and reused across identical source calls; its first/adjoint contributions must be accumulated. The marked f-clock calls g(B(cX+vz)) and g(cX+vz) on their literal common input, with the matching source-zero values. Shared identical anchors may be cached, but distinct active source points and finite iteration versions cannot.

This supplement does not change the source or target of0b3e0dcd, its covariance/current comparison, or its retained unmarked-body price. It supplies an explicit known-fill implementation for the projection-D specialization rather than leaving Q_fill abstract.
