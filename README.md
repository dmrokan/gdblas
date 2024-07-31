# GDBlas
This native [Godot](https://github.com/godotengine/godot) extension provides real and complex matrix algebra. It uses data structures and matrix iterators of [Eigen](https://gitlab.com/libeigen/eigen) library, also includes ODE solver based on [ODEINT](https://github.com/headmyshoulder/odeint-v2).

In version `1.4.0`, [BoostC++ Geometry](https://www.boost.org/doc/libs/1_85_0/libs/geometry/doc/html/index.html) algorithms are added.

## Demos
1. There is a demo project in `demo` directory which includes numerous tests and displacement simulation of a flexible structure. Its mathematical model can be found in my [PhD thesis](https://www.proquest.com/openview/28b57f84e375831a4f1ae27be456ba2d/1?pq-origsite=gscholar&cbl=2026366&diss=y) (Chapter 6).
2. 3D demo project based on a [Godot example](https://github.com/godotengine/godot-demo-projects/releases/download/4.2-31d1c0c/3d_occlusion_culling_mesh_lod.zip). It displays filtered versions of the texture in `Viewport` as shown below.

![Demo3D screenshot](docs/demo3d_screenshot.jpg?raw=true)

https://gist.github.com/user-attachments/assets/aa9c1056-6da8-4c55-9c86-5937d3cd315c

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

# Classes
## GDBlas
Reference counted base class which is used to create new matrices

### Methods
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
- `mat_to_image_data(p_mat_array: Array, p_channel_width: int = 1)`: Places the entries of `GDBlasMat` objects in `p_mat_array` into a `PackedByteArray` which matches the data structure returned from `Image::get_data()`.
```gdscript
var gbl = GDBlas.new()
var dim = Vector2i(480, 640)
var R = gbl.new_mat(dim)
var G = gbl.new_mat(dim)
var B = gbl.new_mat(dim)

# fill and process R, G, B matrices

var pack: PackedByteArray = gbl.mat_to_image_data([ R, G, B ])
# An RGB8 formatted Image object can be created by using the data in 'pack'
```
- `area(p_polygon: PackedVector2Array)`: Calculates the area of polygon. Points must be ordered in clockwise (CW) direction.

### Boost Geometry

#### Data structures

List of Boost Geometry [data structures](https://www.boost.org/doc/libs/1_85_0/libs/geometry/doc/html/geometry/reference/models.html) and their representations in GDScript to be used by `GDBlas`'s bindings of Boost Geometry [algorithms](https://www.boost.org/doc/libs/1_85_0/libs/geometry/doc/html/geometry/reference/algorithms.html).

- `model::point` &equiv; `Vector2`
- `model::linestring` &equiv; `PackedVector2Array`
- `model::ring` &equiv; `PackedVector2Array` if the first and last values of the array are equal.
- `model::polygon` &equiv; `Array` of `PackedVector2Array` where the first entry represents the outer ring and the other entries represent inner rings.
  - `ring_type model::polygon::outer` &equiv; `PackedVector2Array` where the first and last values of the array are equal.
  - `ring_type model::polygon::inners[i]` &equiv; `PackedVector2Array` where the first and last values of the array are equal.
- `model::box` &equiv; `Rect2`

#### Algorithms

- [`area`](https://www.boost.org/doc/libs/1_85_0/libs/geometry/doc/html/geometry/reference/algorithms/area.html)
- [`buffer`](https://www.boost.org/doc/libs/1_85_0/libs/geometry/doc/html/geometry/reference/algorithms/buffer.html)
- [`centroid`](https://www.boost.org/doc/libs/1_85_0/libs/geometry/doc/html/geometry/reference/algorithms/centroid.html)
- [`clear`](https://www.boost.org/doc/libs/1_85_0/libs/geometry/doc/html/geometry/reference/algorithms/clear.html)
- [`closest_points`](https://www.boost.org/doc/libs/1_85_0/libs/geometry/doc/html/geometry/reference/algorithms/closest_points.html)
- [`convert`](https://www.boost.org/doc/libs/1_85_0/libs/geometry/doc/html/geometry/reference/algorithms/convert.html)
- [`convex_hull`](https://www.boost.org/doc/libs/1_85_0/libs/geometry/doc/html/geometry/reference/algorithms/convex_hull.html)
- [`correct`](https://www.boost.org/doc/libs/1_85_0/libs/geometry/doc/html/geometry/reference/algorithms/correct.html)
- [`covered_by`](https://www.boost.org/doc/libs/1_85_0/libs/geometry/doc/html/geometry/reference/algorithms/covered_by.html)
- [`crosses`](https://www.boost.org/doc/libs/1_85_0/libs/geometry/doc/html/geometry/reference/algorithms/crosses.html)
- [`densify`](https://www.boost.org/doc/libs/1_85_0/libs/geometry/doc/html/geometry/reference/algorithms/densify.html)
- [`difference`](https://www.boost.org/doc/libs/1_85_0/libs/geometry/doc/html/geometry/reference/algorithms/difference.html)
- [`discrete_frechet_distance`](https://www.boost.org/doc/libs/1_85_0/libs/geometry/doc/html/geometry/reference/algorithms/discrete_frechet_distance.html)
- [`discrete_hausdorff_distance`](https://www.boost.org/doc/libs/1_85_0/libs/geometry/doc/html/geometry/reference/algorithms/discrete_hausdorff_distance.html)
- [`disjoint`](https://www.boost.org/doc/libs/1_85_0/libs/geometry/doc/html/geometry/reference/algorithms/disjoint.html)
- [`distance`](https://www.boost.org/doc/libs/1_85_0/libs/geometry/doc/html/geometry/reference/algorithms/distance.html)
- [`envelope`](https://www.boost.org/doc/libs/1_85_0/libs/geometry/doc/html/geometry/reference/algorithms/envelope.html)
- [`equals`](https://www.boost.org/doc/libs/1_85_0/libs/geometry/doc/html/geometry/reference/algorithms/equals.html)
- [`intersection`](https://www.boost.org/doc/libs/1_85_0/libs/geometry/doc/html/geometry/reference/algorithms/intersection.html)
- [`intersects`](https://www.boost.org/doc/libs/1_85_0/libs/geometry/doc/html/geometry/reference/algorithms/intersects.html)
- [`is_empty`](https://www.boost.org/doc/libs/1_85_0/libs/geometry/doc/html/geometry/reference/algorithms/is_empty.html)
- [`is_simple`](https://www.boost.org/doc/libs/1_85_0/libs/geometry/doc/html/geometry/reference/algorithms/is_simple.html)
- [`is_valid`](https://www.boost.org/doc/libs/1_85_0/libs/geometry/doc/html/geometry/reference/algorithms/is_valid.html)
- [`length`](https://www.boost.org/doc/libs/1_85_0/libs/geometry/doc/html/geometry/reference/algorithms/length.html)
- [`overlaps`](https://www.boost.org/doc/libs/1_85_0/libs/geometry/doc/html/geometry/reference/algorithms/overlaps.html)
- [`perimeter`](https://www.boost.org/doc/libs/1_85_0/libs/geometry/doc/html/geometry/reference/algorithms/perimeter.html)
- [`relation`](https://www.boost.org/doc/libs/1_85_0/libs/geometry/doc/html/geometry/reference/algorithms/relation.html)
- [`reverse`](https://www.boost.org/doc/libs/1_85_0/libs/geometry/doc/html/geometry/reference/algorithms/reverse.html)
- [`simplify`](https://www.boost.org/doc/libs/1_85_0/libs/geometry/doc/html/geometry/reference/algorithms/simplify.html)
- [`sym_difference`](https://www.boost.org/doc/libs/1_85_0/libs/geometry/doc/html/geometry/reference/algorithms/sym_difference.html)
- [`touches`](https://www.boost.org/doc/libs/1_85_0/libs/geometry/doc/html/geometry/reference/algorithms/touches.html)
- [`transform`](https://www.boost.org/doc/libs/1_85_0/libs/geometry/doc/html/geometry/reference/algorithms/transform.html)
- [`union_`](https://www.boost.org/doc/libs/1_85_0/libs/geometry/doc/html/geometry/reference/algorithms/union_.html)
- [`unique`](https://www.boost.org/doc/libs/1_85_0/libs/geometry/doc/html/geometry/reference/algorithms/unique.html)
- [`within`](https://www.boost.org/doc/libs/1_85_0/libs/geometry/doc/html/geometry/reference/algorithms/within.html)

#### Return types and errors

The algorithms above may not be used all combinations of `model`s. In such cases and error message is emitted and function returns an invalid value. The map of return types and invalid values are given below.

- `bool` return types: Returns `int`. `1` = `true`, `0` = `false`, negative value on error.
- `int` return types: Returns `int`, negative value on error.
- `double` return types: Returns `float`, `NaN` on error.
- `model::polygon` return types: Returns `Array` of `PackedVector2Array`, empty `Array` on error.
- `model::ring` return types: Returns `PackedVector2Array`, empty `PackedVector2Array` on error.
- `model::line` return types: Returns `PackedVector2Array`, empty `PackedVector2Array` on error.
- `model::point` return types: Returns `Vector2`, `Vector2(NaN, NaN)` on error.
- `model::box` return types: Returns `Rect2`, `Rect2(NaN, NaN, NaN, NaN)` on error.

#### Examples

You can check the tests in (demo/GDBlasTest.gd) for their usage.

## GDBlasMat
Reference counted real or complex matrix object. A real matrix returns enries as a `float` and complex matrix as `Vector2`.

### Methods
- `resize(m: Variant, n: int = -1)`: Resizes matrix to m by n if both are integer. `n` is not required if `m` is `Vector2i`.
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
- `argmin(axis = -1)`: Finds the index of minimum of matrix entries on given axis. `axis=0`: min over rows, `axis=1`: min of cols, `axis=-1`: (default) min of all entries. Returns an `Array` of `Vector2i` containing the indices of minimums.
- `argmax(axis = -1)`: Finds the index of maximum of matrix entries on given axis. `axis=0`: max over rows, `axis=1`: max of cols, `axis=-1`: (default) max of all entries. Returns an `Array` of `Vector2i` containing the indices of minimums.
- `norm(norm_type)`: Computes $L1$, $L_{\infty}$ or Frobenius norm of matrix. Accepted arguments are `GDBlas.NORM_1`, `GDBlas.NORM_INF` or `GDBlas.NORM_FRO`. Returns `float`.
- `eval_ode(p_f: Callable, p_dt: float, p_max_step: float = 1e-2)`: Evaluates the ordinary differential equation (ODE) defined in `p_f` for an amount of time given by `p_dt` starting from the current value of matrix. It can be called on only real n by 1 matrices (equivalent to a column vector). Returns the step count (how many times the ODE function is evaluated) or a negative value on error.
```gdscript
var A: GDBlasMat = null
func ode_fx(x: GDBlasMat, t: float):
	return A.prod(x)

func some_func():
	var gbl = GDBlas.new()
	A = gbl.new_mat(2)
	A.eye(-1)
	var x = gbl.new_mat(2, 1)
	x.fill(1)

	x.eval_ode(ode_fx, 0.5, 1e-3) # Writes final value at t = 0.5 into x itself
```
- `conv(p_other: GDBlasMat, p_same: bool = false)`: Computes convolution of matrices. If `p_same` is `true`, returns the central part of the result.
```gdscript
var gbl = GDBlas.new()
A = gbl.new_mat(m1, n1)
A.from_array( ... ) # Fill with values.
B = gbl.new_mat(m2, n2)
B.from_array( ... ) # Fill with values.
var C = A.conv(B)
assert(C.dimension() == Vector2i(m1 + m2 - 1, n1 + n2 -1))
var D = B.conv(A, 'same')
assert(D.dimension() == Vector2i(m2, n2))
```
- `pack(p_component: int = GDBlas.BOTH_COMPONENTS)`: Packs matrix entries into a `PackedFloat64Array` in row major format. Argument `p_component` can take values `GDBlas.REAL_COMPONENT`, `GDBlas.IMAG_COMPONENT` or `GDBlas.BOTH_COMPONENTS`. If both components of a complex matrix is packed imaginary part of each entry is placed right after the real component in the `PackedFloat64Array`. Returns `PackedFloat64Array`.
- `unpack(p_packed_data: PackedFloat64Array, p_component: int = GDBlas.BOTH_COMPONENTS, p_step: int = 1, p_offset: int = 0)`: Unpacks the data in `p_packed_data` into the matrix. If `p_step = n`, each nth entry in the `p_packed_data` placed into the matrix starting from the entry indexed by `p_offset`. Number of elements in the array divided by `p_step` must match the matrix dimension.
- `downsample(p_factor_m: int, p_factor_n: int, p_filter: GDBlasMat = null)`: Returns a new matrix constructed by picking rows and columns whose index satisfy `row_index % p_factor_m == 0` and `col_index % p_factor_n == 0`. If `p_filter` provided, the matrix is filtered (by convolving) by the the coefficients of `p_filter` before down sampling.

### Elementwise functions
A list of implemented math functions are given below. They operate elementwise on the matrix and modifies matrix itself instead of creating a copy. You can visit C++ stdlib documentation for mathematical meaning of these functions.

- `f(p_func: Callable, p_args: Array = null, p_indexed: bool = false)`: Applies `p_func` on each matrix entry and writes the result in place. If `p_indexed` is `true`, `p_func` can have additional argumens which gives the row and column number of the matrix entry.
```gdscript
func add_1(a):
	return a + 1

func add_const(a, args: Array):
	return a + args[0]

func add_const_2(a, args: Array, i: int, j: int):
	if i < 1 and j < 1:
		return a + args[0]
	else:
		return a

func some_func():
	var gbl = GDBlas.new()
	var A = gbl.new_mat(2, 2)
	A.f(add_1)
	A.f(add_const, [ 3 ])
	A.f(add_const_2, [ 3 ], true)
```
**NOTE**: If the matrix is complex, first argument of callable is a `Vector2` and its return type must also be `Vector2`.

- `sin()`

```gdscript
var gbl = GDBlas.new()
var vec = gbl.linspace(-0.5, 0.5, 256)
vec.mul(2 * PI)
vec.sin() # Calculate sine of vec and write in place.
```

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
