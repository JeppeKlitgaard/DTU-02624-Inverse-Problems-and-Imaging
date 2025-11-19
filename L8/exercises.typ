== 1

// Consider a finite difference operator $D$ defined as :
// $
//   D f(x) = (f(x+h) - f(x-h))/h
// $

// Where $h$ is the grid spacing arising from the discretisation.

// $
// f(t_j) = (f(t_j + h/2) - f(t_j - h/2))/h = (f(x_j) - f(x_(j-1)))/h
// $

We consider the double distance version, corresponding to difference operator $D$:
$
  D f(x) = (f(x+h) - f(x-h))/(2h)
$

Which when discretised on a grid $x_j = j h$ becomes:
$
  f(x_j) = (f(x_(j+1)) - f(x_(j-1)))/(2h)
$

Corresponding to a system of equations:
$
  A = 1/(2h) mat(delim: "[",
  0, 1, ,,,;
  -1, 0, 1,,;
  ,⋱, ⋱, ⋱,;
  ,, -1, 0, 1;
  ,,, -1, 0;
  )
$

This is different from given answer, which is suppose to be:
$
  hat(L)_1 =
  1/(2h)
  mat(delim: "[",
  -1, 0, 1, ;
  , ⋱, ⋱, ⋱,;
  ,, -1, 0, 1;
  )
$

=== Null Space

Solving for the null space involves:
$
  hat(L)_1 f = 0
$

Given the repeating nature of the matrix, we can restrict ourselves to
consider only a small subset of the problem.

By inspection, we find that any vector of the form:
$
  f = mat(f_"odd,1", f_"even,1", f_"odd,2", f_"even,2" ...)
$
Will satisfy the equation for:
$
  f_"odd,i" &&= c_1\
  f_"even,i" &&= c_2
$

From which we recover the expected solution by considering the case where the signs of $c_1$ and $c_2$ are the same and the case where they are different:
$
  cal(N)(hat(L)_1) = "span"{mat(1, 1, 1, 1, ...), mat(1, -1, 1, -1, ...)}
$

=== Why delicate?

If we have noise, then the first difference operator will be very sensitive to the noise. Suppression of this noise was exactly what we were trying to reduce.

This will mean we accidentally may include some components just because the noise happened to be high for them.

Answer from Kim:

From null space we see that our regularisation scheme ignores levels and the highest frequency modes (where the sign changes every step). This is delicate because noise often has high frequency components, so these will not be penalised by the regularisation and may end up dominating the solution.

== 2

