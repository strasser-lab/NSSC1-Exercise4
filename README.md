The discretization of the Laplace equation

$$
\Delta u = \frac{\partial u}{\\partial x} + \frac{\partial u}{\partial y} = 0
$$

with a five point stencil is given by:

$$
u_{i-1,j}^{(n+1)} -2u_{i,j}^{(n+1)} +u_{i+1,j}^{(n+1)} + u_{i,j-1}^{(n+1)} -2u_{i,j}^{(n+1)} + u_{i,j+1}^{(n+1)} = 0
$$

This equation is also used to check the solution and calculate the residual.
When we rearrange the equation above, we obtain the update rule for a point:

$$
u_{i,j}^{(n+1)} =\frac{u_{i-1,j}^{(n)} + u_{i,j-1}^{(n)} +u_{i+1,j}^{(n)} + u_{i,j+1}^{(n)}}{4}
$$

This equation can be interpreted as calculating the mean value of the neighbours at that point .

The solutions of the two functions `Jacobi2D_A` and `Jacobi2D_C` match up very closely.

![](plots/Jacobi2D_A.png)
![](plots/Jacobi2D_C.png)
