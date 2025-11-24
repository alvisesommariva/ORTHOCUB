# ORTHOCUB
ORTHOCUB: integral and  differential cubature rules by orthogonal moments: Integration and differentiation by orthogonal moments


In what follows, we mention the routines that are stored in this project.

• cheap_startup: given the algebraic degree of precision ade, it computes a low cardinality 
cubature rule rule_ref on the reference square [−1, 1]^2 or cube [−1, 1]^3 w.r.t.
tensorial-Chebyshev weight (17) or (18).
Next, calling the subroutine dCHEBVAND_orthn, it evaluates the tensorial-orthonormal
Chebyshev basis {ψ_j}_j of total degree ade at the nodes, storing the results in the Van-
dermonde matrix V_ref. The sequence of indices (h, k) or (h, k, l),
defining the ordering of the basis elements, is listed in the matrix basis_indices.
We stress that this routine in independent of the functional L to be used and is only related
to the ADE of the reference cubature rule.

• cheap_rule: given the matrix dbox of dimension 2 × d in which the i column con-
tains the extrema of the bounding box w.r.t. the i-th variable, the moments m storing
L(ψ j({Λ−1}(· − C))), as well as rule_ref, V_ref computed by cheap_startup, it determines 
the nodes X and the weights w of the ORTHOCUB cubature rule.

• dCHEBVAND_orthn evaluates the Vandermonde matrix on the set X ⊂ Ω ⊆ B, relatively to
the scaled tensorial-Chebyshev {φ_j}_j. As input, it is requires the variable
dbox describing B and the sequence of degree indices basis_indices.

• vandermonde_jacobi allows the evaluation of the tensorial-Jacobi polynomials. 
This general purpose function takes into account a specific normalization of
the basis, e.g. monic, orthonormal as well classical ones depending on the weight.

• cub_square_mpx, cub_cube_mpx are low cardinality cubature rules, respectively on
the unit-square [−1, 1]^2 and on the cube [−1, 1]^3 relatively to the tensorial-Chebyshev
measure.

• mom_derivative_2D, mom_derivative_3D produce the derivatives of order with re-
spect to a variable specified by the user of the polynomial basis {φ_k} on a specific 2D or
3D mesh. Depending on the mesh, it is possible to define a suitable bounding box B via
the variable dbox.

As for the examples regarding numerical cubature, the algorithms require some additional functions.

• compute_spline_boundary is useful in the case of the spline-curvilinear domains. It
provides the description of its boundary, given as inputs the abscissae of the vertices XV
and the relative ordinates YV , a variable spline_parms storing the order of the piecewise
splines and the vertices involved. In case the spline is cubic, spline_type fixes the additional 
conditions (e.g. natural, periodic). As output, the routine supplies the vector of
structures Sx and Sy that determine the piecewise splines x = x(t), y = y(t) parametrizing
the boundary.

• spline_chebmom, given the polynomial degree n, the vectors of structures Sx, Sy, the
basis ordering basis_indices, computes the moments of the basis {ψ_j}_j and, if required,
of {φ_j}_j. This purpose is obtained by applying Green theorem, requiring the numerical
integration of a certain piecewise polynomial function on the boundary.

• QMC_union_balls, given the vectors of centers and radii defining the balls Bi, i =
1, 2, . . ., computes card Halton points XB in the bounding box B of the integration domain 
Ω = ∪_i B_i. Next, after the application of an in-domain function, provides as output
the nodes and the weights of a QMC rule in Ω.


5.2 Description of the demos

In order to replicate the examples, we have implemented the following demos, used for the
numerical experiments above.

• demo_cubature_splines illustrates the basic usage of the ORTHOCUB technique to
approximate definite integrals on a spline-curvilinear domain Ω1 (see Figure 1).
• demo_cubature_ade_splines, fixed a degree of precision n equal to 2, 4, . . . , 16, approximates
via the rules proposed in this work 100 integrals of random polynomials of
the form
p_n(x, y) = (c_0 + c_1 x + c_2 y)^n.
Each exact integral is computed via Green theorem. Finally it plots the
single relative errors and the value of their geometric mean.

• demo_cubature_sumweights_spline first computes ORTHOCUB rules with ADE equal
to 2, 4, . . . , 16 for integration on Ω1 and then their stability ratios as in Table 2.

• demo_cubature_weights_spline, first computes ORTHOCUB rules with ADE equal to
2, 4, . . . , 16 for integration on Ω1, and then sorts the weights, illustrating their distribution
as in Figure 1.

• demo_cubature_QMC shows the basic usage of the ORTHOCUB technique to approximate
definite integrals on the domain Ω2 that is union of 5 balls (see Figure 2).

• demo_cubature_ade_QMC determines via QMC_union_balls two QMC rules SL, SH
on Ω2, based respectively on 105 and 106 in the bounding box B. Next, fixed a degree
of precision n equal to 2, 4, . . . , 16, approximates via the rules proposed in this work 100
integrals of random polynomials of the form
p_n(x, y, z) = (c_0 + c_1 x + c_2 y * c_3 z)^n.
In particular, takes into account SH (p_n) as reference value, while as approximation of
the moments m uses the vector whose i-th component is SL(φ_i).
Finally it plots the single relative errors and the value of their geometric mean.

• demo_cubature_sumweights_QMC computes ORTHOCUB rules with ADE equal to 2, 4, . . . , 16
for QMC integration on Ω1 and their stability ratios as in Table 2.

• demo_cubature_weights_QMC computes ORTHOCUB rules with ADE equal to 2, 4, . . . , 16
for QMC integration on Ω2, sorts the weights and illustrates their distribution as in Figure
2.

• demo_derivative_ade_2D computes by the ORTHOCUB technique described above the
first and second order derivatives of random polynomials of the form
p_n(x, y) = (c_0 + c_1 x + c_2 y)^n
on a mesh of the first 10000 Halton points {P_k} in [−1, 1]^2, for n = 2, 4, . . . , 16. For
any exponent n and point of mesh Pk, we determine by an ORTHOCUB rule on [−1, 1]^2
defined by mom_derivative_2D the numerical approximation of these derivatives on P_k.
Since the exact value is known analytically, we can easily establish the relative error in
norm 2.
Integration and differentiation by orthogonal moments

• demo_derivative_ade_3D computes by the ORTHOCUB technique described above the
first and second order derivatives of random polynomials of the form
p_n(x, y, z) = (c_0 + c_1 x + c_2 y + c_3 z)^n
on a mesh of the first 10000 Halton points {P_k} in [−1, 1]^3, for n = 2, 4, . . . , 16. For
any exponent n and point of mesh P_k, we determine by an ORTHOCUB rule on [−1, 1]^3
defined by mom_derivative_3D the numerical approximation of these derivatives on P_k.
Since the exact value is known analytically, we can easily establish the relative error in
norm 2.

For more details see the paper:

"ORTHOCUB: integral and  differential cubature rules by orthogonal moments", 2025,
by L. Rinaldi, A. Sommariva, M. Vianello.
