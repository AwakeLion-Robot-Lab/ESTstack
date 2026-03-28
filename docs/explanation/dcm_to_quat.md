# DCM to Quaternion Explanation

## Quaternion Convention
paper[14] use JPL quaternion, but `Eigen` and `manif` prefer Hamilton quaternion, also matrix $\mathbf{K}$ have different definition between [14] and [13]. $\mathbf{K} = \begin{bmatrix}
\mathbf{\sigma} & \mathbf{z}^{\top} \\
\mathbf{z} & \mathbf{S}-\mathbf{\sigma I}_3 \\
\end{bmatrix} $ in [13], but [14] is denoted by $\mathbf{K} = \begin{bmatrix}
\mathbf{S}-\mathbf{\sigma I}_3 & \mathbf{z} \\
\mathbf{z}^{\top} & \mathbf{\sigma} \\
\end{bmatrix} $. 
**For convenient, we prefer Hamilton quaternion here.**

## Sarabandi Performance

![sarabandi eta performance table](./assets/sarabandi_eta_performance_table.png)

As paper[12] said: 

> The highest errors are obtained for rotation axes corresponding to points contained in the spherical triangles centered in each octant. Now, if we increase the threshold, these spherical triangles, where the highest errors are obtained, are reduced. For η = 0, they have almost completely disappeared. This gives a visual justification of why the optimal threshold is located at η = 0.

So I set `eta` default as `0.0`, you can change eta via this table.

## Davenport’s q-Method Derivation

![q method derivation page 1](./assets/q_method_derivation_p1.jpg)
![q method derivation page 2](./assets/q_method_derivation_p2.jpg)

