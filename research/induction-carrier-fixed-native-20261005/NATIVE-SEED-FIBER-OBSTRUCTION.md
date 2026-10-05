# The carrier-fixed high-order prior is not the theorem being proved

The order-two seed in the pinned LOW30 source already distinguishes its ordinary conditional prior from a source-zero-carrier-fixed prior. This is an exact native example, not merely an arbitrary rotated Gaussian construction.

Take dimension D, one selected square-gradient source f(u)=r u, a fixed 0<s0<1/4 and c0=sqrt(1-s0^2). The literal pair seed uses

    D_x(U)=r[s0 x+(c0-1)U].

Its higher-chaos tail is zero; the finite tail filter can be chosen to vanish exactly on this linear source. Collapsing the seed's independent baseline Gaussian summands gives

    Y_seed = H + s0 r x + (c0-1)r U,
    H,U independent N(0,I_D).

The selected-pair reference is

    Y_ref = s0 r x + sqrt(1-s0^2 r^2) H_ref.

With x retained, the ordinary conditional W2 error is exactly

    [sqrt(1+(1-c0)^2 r^2)-sqrt(1-s0^2 r^2)] sqrt(D)
      = (1-c0)r^2 sqrt(D)+O(r^4 sqrt(D)).

If the physical source-zero carrier must be identical, H_ref=H, the two conditional residual kernels instead have integrated W2 distance

    sqrt((1-c0)^2 r^2
          +(1-sqrt(1-s0^2 r^2))^2) sqrt(D)
      >= (1-c0)|r| sqrt(D).

The actual seed residual has a global complete first proportional to |r|, and the source is a genuine linear gradient. Thus the native order-two conditional floor cannot be upgraded from O(r^2 sqrt(D)) to a carrier-fixed residual floor of that order.

The new root-seed bridge uses precisely the feature visible here: the O(r) conditional residual fluctuation is centered in a private bank, and its effect on the full smoothed output is O(r^2 sqrt(D)). It keeps that fluctuation until a valid grouped Riesz calculation with the original one keep. It never calls the O(r) conditional fluctuation an O(r^2) strong or fiber error.

This counterexample only addresses the asserted metric upgrade for the native seed. It is not a no-go theorem for a differently designed carrier-preserving compiler, and it does not refute the environment-stable root-seed group bound.
