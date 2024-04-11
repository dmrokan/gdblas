# GDBlas
This native [Godot](https://github.com/godotengine/godot) extension provides real and complex matrix algebra. It uses data structures and matrix iterators of [Eigen](https://gitlab.com/libeigen/eigen) library. Also, there is a demo project with numerous tests and displacement simulation of a flexible structure. Its mathematical model can be found in my [PhD thesis](https://www.proquest.com/openview/28b57f84e375831a4f1ae27be456ba2d/1?pq-origsite=gscholar&cbl=2026366&diss=y) (Chapter 6).

**An example**:
```gdscript
var gbl = GDBlas.new()
var A = gbl.new_mat(3, 2)
var b = gbl.new_mat(2, 1)
A.from_array([ [1, 2], [3, 4], [5, 6] ])
b.from_array([ [1], [-1] ])
var c = A.prod(b)
print(c.to_array())
c.abs()
print(c.to_array())
c.log()
print(c.to_array())
c.add(3)
print(c.to_array())
```
will print out
```
>>> [ [-1], [-1], [-1] ]
>>> [ [1], [1], [1] ]
>>> [ [0], [0], [0] ]
>>> [ [3], [3], [3] ]
```

# Objects
## GDBlas
Reference counted base class which is used to create new matrices

### Functions
- `new_mat(p_rows, p_cols = -1)`: Creates new real matrix, Usage:
```gdscript
var gbl = GDBlas.new()
var A = gbl.new_mat(3, 12) # Creates a 3 by 12 real matrix
var B = gbl.new_mat(3) # Creates a 3 by 3 real matrix
var C = gbl.new_mat(Vector2i(4, 2)) # Creates a 4 by 2 real matrix
```
- `new_complex_mat(p_rows, p_cols = -1)`: Creates new complex matrix
- `linspace(p_start, p_end, p_count)`: Creates a column vector of linearly spaced values
```gdscript
var gbl = GDBlas.new()
var A = gbl.linspace(0, 1, 3) # Creates a 3 by 1 matrix with entries [ [0], [0.5], [1] ]
```

## GDBlasMat
Reference counted real or complex matrix object. A real matrix returns enries as a `float` and complex matrix as `Vector2`.

### Functions
- `resize(m, n)`: Resizes matrix to m by n
- `copy()`: Creates a copy of matrix.
```gdscript
var gbl = GDBlas.new()
var A = gbl.new_mat(3)
var B = A.copy()
```
- `dimension()`: Returns the size of matrix as a `Vector2i` object
- `get(i, j, m = -1, n = -1)`: Get a matrix entry or a submatrix of size m by n starting ith row and jth column. Returns `0` on success or error value.
```gdscript
var gbl = GDBlas.new()
var A = gbl.new_mat(3)
var a00: float = A.get(0, 0)
var Asub: GDBlasMat = A.get(1, 0, 2, 2) # Returns 2 by 2 sub matrix
var Ac = gbl.new_complex_mat(3)
var ac00: Vector2 = A.get(0, 0)
```
- `set(val, i = -1, j = -1)`: Set a matrix entry or a submatrix of size m by n starting ith row and jth column. Returns `0` on success or error value.
```gdscript
var gbl = GDBlas.new()
var A = gbl.new_mat(3)
var a00: float = A.set(1, 0, 0)
var Asub = gbl.new_mat(2)
A.set(Asub, 1, 0)
```
- `add(val)`: Adds a number or a matrix of same dimension. Returns `0` on success or error value.
```gdscript
var gbl = GDBlas.new()
var A = gbl.new_mat(3)
A.add(1)
var B = A.copy()
B.add(A)
```
- `mul(val)`: Multiplies by a number or a matrix of same dimension. Returns `0` on success or error value.
- `div(val)`: Divides by a number or a matrix of same dimension. Returns `0` on success or error value.
**NOTE**: Can not add complex matrices to a real matrix
- `sub(val)`: Subtracts a number or a matrix of same dimension. Returns `0` on success or error value.
- `transpose()`: Transposes the matrix.  Returns `0` on success or error value.
- `hermitian()`: Hermitian transpose of matrix. Returns `0` on success or error value.
- `is_eq(other, p_eps, p_norm_type)`: Checks if matrices are equal, meaning that norm of their differences are less than `p_eps` . Returns `true` or `false`
```gdscript
var gbl = GDBlas.new()
var A = gbl.new_mat(3)
var B = A.copy()
if A.is_eq(B):
    print("They are equal")
```
- `fill(val)`: Sets all matrix entries to `val`. `val` can be a number or `Vector2` if it is a complex matrix. Returns `0` on success.
- `eye(val)`: Sets all diagonal matrix entries to `val`. `val` can be a number or `Vector2` if it is a complex matrix. Returns `0` on success.
- `reset()`: Sets all entries to `0`.
- `conj()`: Conjugates all matrix entries. Returns `0` on success.
- `real(matrix)`: Returns or sets the real part of complex matrix. It is equivalent to `set(matrix, 0, 0)` on real case.
- `imag(matrix)`: Returns or sets the imaginary part of complex matrix. It returns a matrix of zeros for real matrix case.
- `prod(matrix)`: Returns product of matrices.
```gdscript
var gbl = GDBlas.new()
var A = gbl.new_mat(3)
var B = gbl.new_mat(3, 5)
var C = A.prod(B)
```
It is equivalent to `C=AB`, column count of `A` and row count of `B` must be equal.
- `inv()`: Computes the inverse of matrix. It can only be applied to square matrices. It will return the inverse matrix or `null` if matrix is singular.
- `from_array(arr)`: Sets matrix entries according to the values on 2 dimensional array `arr`. Return `0` on success.
```gdscript
var gbl = GDBlas.new()
var A = gbl.new_mat(3, 2)
A.from_array([ [1, 2], [3, 4], [5, 6] ])
```
- `to_array()`: Returns a 2 dimensional array filled with matrix entries.
```gdscript
var gbl = GDBlas.new()
var A = gbl.new_complex_mat(2, 2)
A.eye(Vector2(1, -1))
print(A.to_array())
```
will print
```
>>> [ [(1,-1), (0,0)], [(0,0), (1,-1)] ]
```
- `integrate(axis = -1)`: Calculates sum of matrix entries on given axis. `axis=0`: sums over rows, `axis=1`: sums of cols, `axis=-1`: (default) sum of all entries. Returns a `GDBlasMat` object containing the sum of values.
- `mean(axis = -1)`: Calculates mean of matrix entries on given axis. `axis=0`: mean over rows, `axis=1`: mean of cols, `axis=-1`: (default) mean of all entries. Returns a `GDBlasMat` object containing the means.
- `min(axis = -1)`: Finds the minimum of matrix entries on given axis. `axis=0`: min over rows, `axis=1`: min of cols, `axis=-1`: (default) min of all entries. Returns a `GDBlasMat` object containing the minimums.
- `max(axis = -1)`: Finds the maximum of matrix entries on given axis. `axis=0`: max over rows, `axis=1`: max of cols, `axis=-1`: (default) max of all entries. Returns a `GDBlasMat` object containing the maximums.
- `argmin(axis = -1)`: Finds the index of minimum of matrix entries on given axis. `axis=0`: min over rows, `axis=1`: min of cols, `axis=-1`: (default) min of all entries. Returns an `Array` of `Vector2i` containing the index o minimums.
- `argmax(axis = -1)`: Finds the index of maximum of matrix entries on given axis. `axis=0`: max over rows, `axis=1`: max of cols, `axis=-1`: (default) max of all entries. Returns an `Array` of `Vector2i` containing the index o minimums.
- `norm(norm_type)`: Computes $L1$, $L_{\infty}$ or Frobenius norm of matrix. Accepted arguments are `GDBlas.NORM_1`, `GDBlas.NORM_INF` or `GDBlas.NORM_FRO`. Returns `float`.
- `eval_ode(p_f: Callable, p_dt: float, p_max_step: float = 1e-2)`: Evaluates the ordinary differential equation (ODE) defined in `p_f` for an amount of time given by `p_dt` starting from the current value of matrix. It can be called on only real n by 1 matrices (equivalent to a column vector). Returns the step count (how many times the ODE function is evaluated) or a negative value on error.
```gdscript
var A: GDBlasMat = null
func ode_fx(x: GDBlasMat, t: float):
	return a + args[0]

func some_func():
	var gbl = GDBlas.new()
	A = gbl.new_mat(2)
	A.eye(-1)
	var x = gbl.new_mat(2, 1)
	x.fill(1)

	x.eval_ode(ode_fx, 0.5, 1e-3) # Writes final value at t = 0.5 into x itself
```

### Elementwise functions
A list of implemented math functions are given below. They operate elementwise on the matrix and modifies matrix itself instead of creating a copy. You can visit C++ stdlib documentation for mathematical meaning of these functions.

- `f(p_func: Callable, p_args: Array)`: Applies `p_func` on each matrix entry and writes the result in place.
```gdscript
func add_1(a):
	return a + 1

func add_const(a, args):
	return a + args[0]

func some_func():
	var gbl = GDBlas.new()
	var A = gbl.new_mat(2, 2)
	A.f(add_1)
	A.f(add_const, [ 3 ])
```
**NOTE**: If the matrix is complex, first argument of callable is a `Vector2` and its return type must also be `Vector2`.

- `sin()`
- `cos()`
- `abs()`
- `exp()`
- `log()`
- `log10()`
- `log2()`
- `sqrt()`
- `cbrt()`
- `tan()`
- `asin()`
- `acos()`
- `atan()`
- `sinh()`
- `cosh()`
- `tanh()`
- `atanh()`
- `erf()`
- `erfc()`
- `tgamma()`
- `lgamma()`
- `ceil()`
- `floor()`
- `trunc()`
- `round()`

## Build from source
```sh
git clone https://github.com/dmrokan/gdblas.git
cd gdblas
git submodule update --init --recursive
```
**Build boost**
```sh
cd boost
./bootstrap.sh
./b2 headers
```
You can visit [Boost wiki](https://github.com/boostorg/wiki/wiki) for more information.

**Build extension**
```sh
scons platform=<platform> target_path=<target_path> target_name=libgdblas
```
You can visit Godot's [build system](https://docs.godotengine.org/en/stable/contributing/development/compiling/introduction_to_the_buildsystem.html) documentation for more information.
