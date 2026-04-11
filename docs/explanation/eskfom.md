# ESKFOM Explanation

## Perturbation Convention

I just use right(local)-perturbation as same as [4] and [7], but left(global)-perturbation maybe more efficiency according to (Li and Mourikis, 2012).
In fact, i've tried combining two kinds of perturbation before via a template named `template<typename Perturbation = eststack::Perturbation::Right>`, but it's too annoying for appropriate coding and unnecessary:<, but you have to focus on left/right invariant and equivariant properties of lie group in Invariant Kalman Filter/Equivariant Filter.

## Reset Jacobian Derivation

we need some tricks about lie bracket $[\cdot,\cdot]$ to derive how small adjoint related to wedge mapping, here provides derivation in $\boldsymbol{SO}(3)$:

![reset jacobian derivation](./assets/reset_jacobian_derivation.jpg)

**NOTE:** if manifold is $\boldsymbol{SE}(3)$ or $\boldsymbol{SE}_2(3)$, the skew-symmetric matrix would not equal to small adjoint matrix since $(\cdot)^{\wedge}: \mathbb{R}^{n} \to \mathbb{R}^{m\times m}$, but $ \boldsymbol{ad}_{(\cdot)}: \mathbb{R}^n \to \mathbb{R}^{n\times n}$, here we define $m$ is dimension of lie group (`manif::Dim`) and $n$ is dimension of lie algebra (`manif::DoF`).