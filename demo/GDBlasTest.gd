extends Node2D

var _has_geometry_funcs: bool = false

func C(real: float, imag: float = 0.0) -> Vector2:
	return Vector2(real, imag)

func _is_real_number(v) -> bool:
	if v is int or v is float:
		return true

	return false

func _is_complex_number(v) -> bool:
	if v is Vector2:
		return true

	return false

func _is_number(v) -> bool:
	return _is_real_number(v) or _is_complex_number(v)

func cmp_scalar(s1, s2, eps: float = 1e-7) -> bool:
	if not _is_number(s1) or not _is_number(s2):
		return false

	if _is_real_number(s1) and _is_complex_number(s2):
		s1 = Vector2(s1, 0)
	if _is_real_number(s2) and _is_complex_number(s1):
		s2 = Vector2(s2, 0)

	if _is_complex_number(s1) and _is_complex_number(s2) and (s1 - s2).length() < eps:
		return true
	elif _is_real_number(s1) and _is_real_number(s2) and abs(s1 - s2) < eps:
		return true

	return false

func cmp_mat(mat1: Array, mat2: Array) -> bool:
	if mat1.size() != mat2.size():
		return false

	for i in range(mat1.size()):
		var row1 = mat1[i]
		var row2 = mat2[i]
		if row1.size() != row2.size():
			return false
		for j in range(row1.size()):
			var val1 = row1[j]
			var val2 = row2[j]

			if not cmp_scalar(val1, val2):
				return false

	return true

func cmp_vec2_array(a1: PackedVector2Array, a2: PackedVector2Array, eps: float = 1e-7):
	var s = a1.size()
	if s != a1.size():
		return false

	for i in range(s):
		if not cmp_scalar(a1[i], a2[i], eps):
			return false

	return true

func random_mat(m: int, n: int, real: bool = true) -> Array:
	var mat = Array()
	mat.resize(m)

	for i in range(m):
		var row = Array()
		row.resize(n)

		for j in range(n):
			if real:
				row[j] = randf_range(-100, 100)
			else:
				row[j] = Vector2(randf_range(-100, 100), randf_range(-100, 100))

		mat[i] = row

	return mat

func get_real_or_imag(mat: Array, real: bool = true) -> Array:
	var out = Array()
	out.resize(mat.size())

	for i in range(out.size()):
		var row = Array()
		row.resize(mat[i].size())

		for j in range(row.size()):
			if real:
				row[j] = mat[i][j].x
			else:
				row[j] = mat[i][j].y

		out[i] = row

	return out

func new_mat(m: int, n: int, val = 0, real: bool = true) -> Array:
	var mat = Array()
	mat.resize(m)

	for i in range(m):
		var row = Array()
		row.resize(n)

		for j in range(n):
			if real:
				if _is_complex_number(val):
					print("ERR: _is_complex_number(val)")
					return Array()

				row[j] = val
			else:
				if _is_real_number(val):
					val = Vector2(val, 0)
				row[j] = val

		mat[i] = row

	return mat

# Called when the node enters the scene tree for the first time.
func _ready():
	var gbl = GDBlas.new()
	if gbl.has_method("area"):
		_has_geometry_funcs = true

	test01()
	test02()
	test02_1()
	test03()
	test04()
	test05()
	test06()
	test07()
	test08()
	test09()
	test10()
	test11()
	test12()
	test13()
	test14()
	test15()
	test16()
	test17()
	test18()
	test19()
	test20()
	test21()
	test22()
	test23()

	if _has_geometry_funcs:
		test24()
		test25()
		test26()
		test27()
		test28()
		test29()
		test30()
		test31()
		test32()
		test33()
		test34()
		test35()
		test36()
		test37()
		test38()
		test39()
		test40()
		test41()
		test42()
		test43()
		test44()
		test45()
		test46()
		test47()
		test48()
		test49()
		test50()
		test51()
		test52()
		test53()
		test54()
		test55()
	else:
		print("Geometry functions are not enabled. Skipping...")

func test01():
	print("test01")

	var gbl = GDBlas.new()
	var mat1 = gbl.new_mat(2, 2)

	var mat1_arr = [ [0.0, 0.0], [0.0, 0.0] ]
	var tmp = mat1.to_array()
	assert(cmp_mat(tmp, mat1_arr))
	assert(mat1.dimension() == Vector2i(2, 2))

	mat1.resize(3, 4)
	mat1_arr = [ [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0] ]
	tmp = mat1.to_array()
	assert(cmp_mat(tmp, mat1_arr))
	assert(mat1.dimension() == Vector2i(3, 4))

	var mat2 = gbl.new_complex_mat(2, 3)
	var mat2_arr = [ [C(0.0), C(0.0), C(0.0)], [C(0.0), C(0.0), C(0.0)] ]
	tmp = mat2.to_array()
	assert(cmp_mat(tmp, mat2_arr))
	assert(mat2.dimension() == Vector2i(2, 3))

	mat2.resize(3, 4)
	mat2_arr = [
		[C(0.0), C(0.0), C(0.0), C(0.0)],
		[C(0.0), C(0.0), C(0.0), C(0.0)],
		[C(0.0), C(0.0), C(0.0), C(0.0)]
	]
	tmp = mat2.to_array()
	assert(tmp == mat2_arr)
	assert(mat2.dimension() == Vector2i(3, 4))

	assert(mat2.set(C(1, -1), 4, 4) != 0)

	mat1.set(1, 1, 0)
	mat1.set(-0.3, 0, 1)
	mat1.set(1e3, 2, 3)
	mat1_arr[1][0] = 1
	mat1_arr[0][1] = -0.3
	mat1_arr[2][3] = 1e3
	tmp = mat1.to_array()
	assert(cmp_scalar(mat1.get(1, 0), 1))
	assert(cmp_scalar(mat1.get(0, 1), -0.3))
	assert(cmp_scalar(mat1.get(2, 3), 1e3))
	assert(cmp_mat(tmp, mat1_arr))

	mat2.set(C(1, -1), 1, 0)
	mat2.set(C(2, -0.3), 0, 1)
	mat2.set(C(1e3, 1e-2), 2, 3)
	mat2_arr[1][0] = C(1, -1)
	mat2_arr[0][1] = C(2, -0.3)
	mat2_arr[2][3] = C(1e3, 1e-2)
	tmp = mat2.to_array()
	assert(cmp_scalar(mat2.get(1, 0), C(1, -1)))
	assert(cmp_scalar(mat2.get(0, 1), C(2, -0.3)))
	assert(cmp_scalar(mat2.get(2, 3), C(1e3, 1e-2)))
	assert(cmp_mat(tmp, mat2_arr))

	var mat3 = mat1.copy()
	assert(cmp_mat(mat3.to_array(), mat1.to_array()))
	var mat4 = mat2.copy()
	assert(cmp_mat(mat4.to_array(), mat2.to_array()))

	mat3.sub(mat1)
	assert(
		cmp_mat(mat3.to_array(), gbl.new_mat(mat3.dimension()).to_array())
	)

func test02():
	print("test02")

	var gbl = GDBlas.new()
	var mat1 = gbl.new_mat(Vector2i(2, 4))
	var mat1_arr = [ [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0] ]
	assert(cmp_mat(mat1.to_array(), mat1_arr))

	mat1_arr = [ [0.2, 3.0, -2.0, 0.0], [0.0, 0.0, 1e-3, 0.0], [0.1, 4.0, 0.3, 10.0] ]
	mat1.from_array(mat1_arr)
	assert(mat1.dimension() == Vector2i(3, 4))
	assert(cmp_mat(mat1.to_array(), mat1_arr))

	mat1 = gbl.new_complex_mat(Vector2i(2, 4))
	mat1_arr = [ [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0] ]
	assert(cmp_mat(mat1.to_array(), mat1_arr))

	mat1_arr = [ [0.2, 3.0, -2.0, 0.0], [0.0, 0.0, 1e-3, 0.0], [0.1, 4.0, 0.3, 10.0] ]
	mat1.from_array(mat1_arr)
	assert(mat1.dimension() == Vector2i(3, 4))
	assert(cmp_mat(mat1.to_array(), mat1_arr))

	mat1.T()
	assert(mat1.dimension() == Vector2i(4, 3))
	mat1.T()
	assert(mat1.dimension() == Vector2i(3, 4))
	assert(cmp_mat(mat1.to_array(), mat1_arr))

	mat1_arr = random_mat(5, 5, false)
	mat1.from_array(mat1_arr)
	var mat2 = mat1.copy()
	assert(not mat2.sub(mat1))
	assert(cmp_mat(mat2.to_array(), new_mat(mat2.dimension().x, mat2.dimension().y)))
	mat2 = mat1.copy()
	mat2.T()
	mat2.H()
	assert(not mat2.sub(mat1))
	assert(not cmp_mat(mat2.to_array(), new_mat(mat2.dimension().x, mat2.dimension().y)))
	mat2 = mat1.copy()
	mat2.H()
	mat2.H()
	assert(not mat2.sub(mat1))
	assert(cmp_mat(mat2.to_array(), new_mat(mat2.dimension().x, mat2.dimension().y)))

	mat2 = mat1.copy()
	var mat3 = mat1.copy()
	mat2.H()
	mat2.H()
	assert(not mat2.add(mat1))
	assert(not mat3.mul(4))
	assert(not mat3.div(Vector2(2, 0)))
	assert(cmp_mat(mat2.to_array(), mat3.to_array()))

func _test02_1_set_mat(a, i: int, j: int):
	return i + j + 1

func test02_1():
	print("test02_1")

	var gbl = GDBlas.new()
	var A = gbl.new_mat(6)
	var A11 = gbl.new_mat(4, 2)
	A.fill(-1)
	A11.f(_test02_1_set_mat, null, true)
	A.set(A11, 1, 1)
	var A11_cmp = A.get(1, 1, 4, 2)
	assert(A11.is_eq(A11_cmp))

	var B = gbl.new_complex_mat(6)
	B.fill(-1)
	B.set(A11, 1, 1)
	var B11_cmp = B.get(1, 1, 4, 2)
	assert(A11.is_eq(B11_cmp))
	var B11 = gbl.new_complex_mat(3, 3)
	B11.fill(Vector2(10, 20))
	assert(B.set(B11, 4, 4) != 0)
	B.set(B11, 3, 3)
	B11_cmp = B.get(3, 3, 3, 3)
	assert(B11.is_eq(B11_cmp))

func test03():
	print("test03")

	var gbl = GDBlas.new()
	var mat1 = gbl.new_mat()
	assert(mat1.dimension() == Vector2i(1, 1))
	var mat1_arr = random_mat(4, 7)
	mat1.from_array(mat1_arr)
	var mat2 = gbl.new_mat()
	var mat2_arr = random_mat(4, 7)
	mat2.from_array(mat2_arr)
	assert(mat2.add(mat1) == 0)
	mat1.T()
	assert(mat2.add(mat1) != 0)
	mat2.T()
	assert(mat2.add(mat1) == 0)
	assert(mat2.sub(mat1) == 0)
	assert(mat2.sub(mat1) == 0)
	mat2.T()
	assert(cmp_mat(mat2.to_array(), mat2_arr))

	var mat1c = gbl.new_complex_mat()
	assert(mat1c.dimension() == Vector2i(1, 1))
	var mat1c_arr = random_mat(4, 7, false)
	mat1c.from_array(mat1c_arr)
	var mat2c = gbl.new_complex_mat()
	var mat2c_arr = random_mat(4, 7, false)
	mat2c.from_array(mat2c_arr)
	assert(mat2c.add(mat1c) == 0)
	mat1c.T()
	assert(mat2c.add(mat1c) != 0)
	mat2c.T()
	assert(mat2c.add(mat1c) == 0)
	assert(mat2c.sub(mat1c) == 0)
	assert(mat2c.sub(mat1c) == 0)
	mat2c.T()
	assert(cmp_mat(mat2c.to_array(), mat2c_arr))

func test04():
	print("test04")

	var gbl = GDBlas.new()
	var mat1 = gbl.new_mat(12)
	assert(mat1.dimension() == Vector2i(12, 12))
	assert(mat1.eye(-0.224) == 0)
	var mat2 = mat1.copy()
	assert(mat2.dimension() == Vector2i(12, 12))
	assert(cmp_mat(mat1.to_array(), mat2.to_array()))
	mat1.T()
	assert(cmp_mat(mat1.to_array(), mat2.to_array()))
	mat1.H()
	assert(cmp_mat(mat1.to_array(), mat2.to_array()))
	var mat1_arr = random_mat(mat1.dimension().x, mat1.dimension().y)
	mat1.from_array(mat1_arr)
	mat2 = mat1.copy()
	assert(cmp_mat(mat1.to_array(), mat2.to_array()))
	mat1.T()
	assert(not cmp_mat(mat1.to_array(), mat2.to_array()))
	mat1.H()
	assert(cmp_mat(mat1.to_array(), mat2.to_array()))

func test05():
	print("test05")

	var gbl = GDBlas.new()
	var mat1 = gbl.new_complex_mat(12)
	assert(mat1.dimension() == Vector2i(12, 12))
	assert(mat1.eye(Vector2(-0.224, 10)) == 0)
	var mat2 = mat1.copy()
	assert(mat2.dimension() == Vector2i(12, 12))
	assert(cmp_mat(mat1.to_array(), mat2.to_array()))
	mat1.T()
	assert(cmp_mat(mat1.to_array(), mat2.to_array()))
	mat1.H()
	assert(not cmp_mat(mat1.to_array(), mat2.to_array()))
	mat1.conj()
	assert(cmp_mat(mat1.to_array(), mat2.to_array()))

	var mat1_arr = random_mat(mat1.dimension().x, mat1.dimension().y, false)
	mat1.from_array(mat1_arr)
	mat2 = mat1.copy()
	assert(cmp_mat(mat1.to_array(), mat2.to_array()))
	mat1.T()
	assert(not cmp_mat(mat1.to_array(), mat2.to_array()))
	mat1.H()
	assert(not cmp_mat(mat1.to_array(), mat2.to_array()))
	mat1.conj()
	assert(cmp_mat(mat1.to_array(), mat2.to_array()))
	assert(mat1.is_eq(mat2, 1e-16, GDBlas.NORM_FRO))
	assert(mat2.is_eq(mat1))

func test06():
	print("test06")

	var gbl = GDBlas.new()
	var mat1 = gbl.new_mat()
	assert(mat1.dimension() == Vector2i(1, 1))
	var mat1_arr = random_mat(4, 7)
	mat1.from_array(mat1_arr)
	var mat2 = gbl.new_mat()
	var mat2_arr = random_mat(4, 7)
	mat2.from_array(mat2_arr)
	assert(mat2.mul(mat1) == 0)
	mat1.T()
	assert(mat2.mul(mat1) != 0)
	mat2.T()
	assert(mat2.mul(mat1) == 0)
	assert(mat2.div(mat1) == 0)
	assert(mat2.div(mat1) == 0)
	mat2.T()
	assert(cmp_mat(mat2.to_array(), mat2_arr))

	mat2.mul(2)
	mat2.div(2)
	assert(cmp_mat(mat2.to_array(), mat2_arr))
	assert(mat2.mul(C(1e3, 0.002)) != 0)
	assert(mat2.div(C(1e3, 0.002)) != 0)
	mat2.mul(0.002)
	assert(not cmp_mat(mat2.to_array(), mat2_arr))
	mat2.div(0.002)
	assert(cmp_mat(mat2.to_array(), mat2_arr))
	var sc = gbl.new_mat()
	sc.set(-7)
	mat2.mul(sc)
	assert(not cmp_mat(mat2.to_array(), mat2_arr))
	mat2.div(sc)
	assert(cmp_mat(mat2.to_array(), mat2_arr))

	var mat1c = gbl.new_complex_mat()
	assert(mat1c.dimension() == Vector2i(1, 1))
	var mat1c_arr = random_mat(4, 7, false)
	mat1c.from_array(mat1c_arr)
	var mat2c = gbl.new_complex_mat()
	var mat2c_arr = random_mat(4, 7, false)
	mat2c.from_array(mat2c_arr)
	assert(mat2c.mul(mat1c) == 0)
	mat1c.T()
	assert(mat2c.mul(mat1c) != 0)
	mat2c.T()
	assert(mat2c.mul(mat1c) == 0)
	assert(mat2c.div(mat1c) == 0)
	assert(mat2c.div(mat1c) == 0)
	mat2c.T()
	assert(cmp_mat(mat2c.to_array(), mat2c_arr))

	mat2c.mul(2)
	mat2c.div(2)
	assert(cmp_mat(mat2c.to_array(), mat2c_arr))
	mat2c.mul(C(1e3, 0.002))
	assert(not cmp_mat(mat2c.to_array(), mat2c_arr))
	mat2c.div(C(1e3, 0.002))
	assert(cmp_mat(mat2c.to_array(), mat2c_arr))
	mat2c.mul(C(1e3, 0.002))
	assert(not cmp_mat(mat2c.to_array(), mat2c_arr))
	mat2c.div(C(1e3, 0.002))
	assert(cmp_mat(mat2c.to_array(), mat2c_arr))
	var scc = gbl.new_complex_mat()
	scc.set(C(2.33, -7))
	mat2c.mul(scc)
	assert(not cmp_mat(mat2c.to_array(), mat2c_arr))
	mat2c.div(scc)
	assert(cmp_mat(mat2c.to_array(), mat2c_arr))

func test07():
	print("test07")

	var gbl = GDBlas.new()
	var mat1 = gbl.new_complex_mat()
	var mat1_arr = random_mat(20, 27, false)
	mat1.from_array(mat1_arr)
	assert(mat1.dimension() == Vector2i(20, 27))
	var mat1re = mat1.real()
	var mat1im = mat1.imag()
	assert(mat1re.dimension() == Vector2i(20, 27))
	assert(mat1im.dimension() == Vector2i(20, 27))
	var mat1re_arr = get_real_or_imag(mat1_arr, true)
	var mat1im_arr = get_real_or_imag(mat1_arr, false)
	assert(cmp_mat(mat1re.to_array(), mat1re_arr))
	assert(cmp_mat(mat1im.to_array(), mat1im_arr))

	mat1.real(mat1im)
	mat1.imag(mat1re)
	var mat1re2 = mat1.real()
	var mat1im2 = mat1.imag()
	assert(cmp_mat(mat1re2.to_array(), mat1im_arr))
	assert(cmp_mat(mat1im2.to_array(), mat1re_arr))
	assert(mat1re.is_eq(mat1im2))
	assert(mat1im.is_eq(mat1re2))

func test08():
	print("test08")

	var gbl = GDBlas.new()
	var a = gbl.new_mat()
	var A = [
		[ 1, 2 ],
		[ 3, 4 ],
		[ 5, 6 ],
		[ 7, 8 ],
	]
	a.from_array(A)
	var b = gbl.new_mat()
	b.from_array(A)
	b.T()
	var c = a.prod(b)
	var C_ = [
		[  5,  11,  17,  23],
 		[ 11,  25,  39,  53],
 		[ 17,  39,  61,  83],
 		[ 23,  53,  83, 113]
	]
	assert(cmp_mat(c.to_array(), C_))

	var d = b.prod(a)
	var D = [
		[84, 100],
		[100, 120],
	]
	assert(cmp_mat(d.to_array(), D))

	c.set(d, 0, 0)
	c.set(d, 2, 2)
	assert(c.set(d, 3, 3) != 0)
	var cinv = c.inv()
	var cinv_arr = [
		[-0.04584144,  0.05921693,  0.3099465 , -0.27565661],
		[ 0.05921693, -0.05881161, -0.26333495,  0.23407101],
 		[ 0.3099465 , -0.26333495,  0.46648022, -0.33183366],
 		[-0.27565661,  0.23407101, -0.33183366,  0.2343142 ]
	]
	assert(cmp_mat(cinv.to_array(), cinv_arr))

	var eye = gbl.new_mat(cinv.dimension())
	eye.eye(1)
	var eye_hat = c.prod(cinv)
	assert(eye.is_eq(eye_hat, 1e-12))
	eye_hat = cinv.prod(c)
	assert(eye.is_eq(eye_hat, 1e-12))

	var E = [
		[  5,  11,  17,  -17],
 		[ 11,  25,  39,  -39],
 		[ 17,  39,  61,  -61],
 		[ 23,  53,  83, -83]
	]
	var e = gbl.new_mat()
	e.from_array(E)
	assert(e.inv() == null)

	var F = [
		[  5,  11],
 		[ 11,  25],
 		[ 17,  39],
 		[ 23,  53]
	]
	var f = gbl.new_mat()
	f.from_array(F)
	assert(f.inv() == null)

func test09():
	print("test09")

	var gbl = GDBlas.new()
	var a = gbl.new_complex_mat(4, 2)

	var areim = gbl.new_mat(a.dimension())
	var A = [
		[ 1, 2 ],
		[ 3, 4 ],
		[ 5, 6 ],
		[ 7, 8 ],
	]
	areim.from_array(A)
	a.real(areim)
	areim.mul(-1)
	a.imag(areim)
	var b = a.copy()
	b.H()

	var c = a.prod(b)
	var C_ = [
		[ 10,  22,  34,  46],
 		[ 22,  50,  78, 106],
 		[ 34,  78, 122, 166],
 		[ 46, 106, 166, 226]
	]
	assert(cmp_mat(c.to_array(), C_))

	var d = b.prod(a)
	var D = [
		[168, 200],
		[200, 240],
	]
	assert(cmp_mat(d.to_array(), D))

	d.imag(d.real())

	var dinv = d.inv()
	var dinv_arr = [
		[ C(0.375, -0.375), C(-0.3125, 0.3125) ],
		[ C(-0.3125, 0.3125), C(0.2625, -0.2625) ]
	]
	assert(cmp_mat(dinv.to_array(), dinv_arr))

	var eye = gbl.new_mat(dinv.dimension())
	eye.eye(1)
	var eye_hat = d.prod(dinv)
	assert(eye.is_eq(eye_hat, 1e-12))
	eye_hat = dinv.prod(d)
	assert(eye.is_eq(eye_hat, 1e-12))

func test10():
	print("test10")

	var gbl = GDBlas.new()
	var mat1 = gbl.new_mat()
	var mat1_arr = [
		[ 10,  22,  34,  -46],
 		[ 22,  50,  78, 106],
 		[ 340,  78, 122, 166],
 		[ 46, 10, 16, 226],
 		[ -460, -100, -160, -26]
	]
	mat1.from_array(mat1_arr)
	assert(mat1.integrate().get(0, 0) == 534)
	assert(cmp_mat(mat1.integrate(0).to_array(), [ [ -42, 60, 90, 426 ] ]))
	assert(cmp_mat(mat1.integrate(1).to_array(), [ [20], [256], [706], [298], [-746] ]))
	assert(cmp_mat(mat1.mean().to_array(), [ [26.7] ]))
	assert(cmp_mat(mat1.mean(0).to_array(), [ [ -8.4, 12, 18, 85.2 ] ]))
	assert(cmp_mat(mat1.mean(1).to_array(), [ [5], [64], [176.5], [74.5], [-186.5] ]))
	assert(cmp_mat(mat1.min().to_array(), [ [-460] ]))
	assert(cmp_mat(mat1.min(0).to_array(), [ [ -460, -100, -160, -46 ] ]))
	assert(cmp_mat(mat1.min(1).to_array(), [ [-46], [22], [78], [10], [-460] ]))
	assert(cmp_mat(mat1.max().to_array(), [ [340] ]))
	assert(cmp_mat(mat1.max(0).to_array(), [ [ 340, 78, 122, 226 ] ]))
	assert(cmp_mat(mat1.max(1).to_array(), [ [34], [106], [340], [226], [-26] ]))
	assert(mat1.argmax() == [ Vector2i(2, 0) ])
	assert(mat1.argmax(0) == [ Vector2i(2, 0), Vector2i(2, 1), Vector2i(2, 2), Vector2i(3, 3) ])
	assert(mat1.argmax(1) == [ Vector2i(0, 2), Vector2i(1, 3), Vector2i(2, 0), Vector2i(3, 3), Vector2i(4, 3) ])
	assert(mat1.argmin() == [ Vector2i(4, 0) ])
	assert(mat1.argmin(0) == [ Vector2i(4, 0), Vector2i(4, 1), Vector2i(4, 2), Vector2i(0, 3) ])
	assert(mat1.argmin(1) == [ Vector2i(0, 3), Vector2i(1, 0), Vector2i(2, 1), Vector2i(3, 1), Vector2i(4, 0) ])

func test11():
	print("test11")

	var gbl = GDBlas.new()
	var mat1 = gbl.new_complex_mat()
	var mat1_arr = [
		[ 10,  22,  34,  -46],
 		[ 22,  50,  78, 106],
 		[ 340,  78, 122, 166],
 		[ 46, 10, 16, 226],
 		[ -460, -100, -160, -26]
	]
	mat1.from_array(mat1_arr)
	mat1.imag(mat1.real())
	mat1.H()
	assert(mat1.integrate().get(0, 0) == C(534, -534))
	assert(cmp_mat(mat1.integrate(0).to_array(), [ [ C(20, -20), C(256, -256), C(706, -706), C(298, -298), C(-746, 746) ] ]))
	assert(cmp_mat(mat1.integrate(1).to_array(), [ [C(-42, 42)], [C(60, -60)], [C(90, -90)], [C(426, -426)] ]))
	assert(cmp_mat(mat1.mean().to_array(), [ [C(26.7, -26.7)] ]))
	assert(cmp_mat(mat1.mean(0).to_array(), [ [ C(5, -5), C(64, -64), C(176.5, -176.5), C(74.5, -74.5), C(-186.5, 186.5) ] ]))
	assert(cmp_mat(mat1.mean(1).to_array(), [ [C(-8.4, 8.4)], [C(12, -12)], [C(18, -18)], [C(85.2, -85.2)] ]))
	assert(mat1.min() == null)
	assert(mat1.max() == null)

func test12():
	print("test12")

	var gbl = GDBlas.new()
	var mat1 = gbl.new_mat()
	var mat1_arr = [
		[ 10,  -22,  34,  46],
 		[ -22,  50,  -78, 106],
 		[ 340,  -78, 122, 166],
 		[ -46, 17, 16, 2.26],
 		[ 46, 110, -16, 26],
 		[ -46, 10, 16, 6]
	]
	mat1.from_array(mat1_arr)
	assert(cmp_scalar(mat1.norm(GDBlas.NORM_1), 510))
	assert(cmp_scalar(mat1.norm(GDBlas.NORM_INF), 706))
	assert(cmp_scalar(mat1.norm(GDBlas.NORM_FRO), 456.734176))

	var mat1c = gbl.new_complex_mat()
	mat1c.from_array(mat1_arr)
	mat1c.imag(mat1c.copy().real())
	mat1c.H()
	assert(cmp_scalar(mat1c.norm(GDBlas.NORM_1), 998.434775))
	assert(cmp_scalar(mat1c.norm(GDBlas.NORM_INF), 721.2489168))
	assert(cmp_scalar(mat1c.norm(GDBlas.NORM_FRO), 645.9196662))

func gdblas_cos(a):
	if a < 0.74:
		return cos(a)

	return -1

func gdblas_sin(a, args):
	if a < 0.74:
		return sin(a + args[0])

	return -1

func test13():
	print("test13")

	var gbl = GDBlas.new()
	var x = gbl.linspace(0, 1, 6)
	x.T()
	assert(cmp_mat(x.to_array(), [ [0, 0.2, 0.4, 0.6, 0.8, 1] ]))

	var y = x.copy(); y.pow(0.3)
	assert(cmp_mat(y.to_array(), [ [0., 0.61703386, 0.75965779, 0.8579172 , 0.93524845, 1.] ]))

	y = x.copy(); y.sin()
	assert(cmp_mat(y.to_array(), [ [0.        , 0.19866933, 0.38941834, 0.56464247, 0.71735609,
	   0.84147098] ]))

	y = x.copy(); y.cos()
	assert(cmp_mat(y.to_array(), [ [1.        , 0.98006658, 0.92106099, 0.82533561, 0.69670671,
	   0.54030231] ]))

	y = x.copy(); y.f(gdblas_cos)
	assert(cmp_mat(y.to_array(), [ [1.        , 0.98006658, 0.92106099, 0.82533561, -1, -1] ]))

	y = x.copy(); y.f(gdblas_sin, [ PI / 2])
	assert(cmp_mat(y.to_array(), [ [1.        , 0.98006658, 0.92106099, 0.82533561, -1, -1] ]))

	y = x.copy(); y.mul(-1); y.abs()
	assert(cmp_mat(y.to_array(), [ [0, 0.2, 0.4, 0.6, 0.8, 1] ]))

	y = x.copy(); y.exp()
	assert(cmp_mat(y.to_array(), [ [1.        , 1.22140276, 1.4918247 , 1.8221188 , 2.22554093,
	   2.71828183] ]))

	y = x.copy(); y.add(1); y.log()
	assert(cmp_mat(y.to_array(), [ [0.        , 	0.18232156, 0.33647224, 0.47000363, 0.58778666,
	   0.69314718] ]))

	y = x.copy(); y.add(1); y.log10()
	assert(cmp_mat(y.to_array(), [ [0.        , 0.07918125, 0.14612804, 0.20411998, 0.25527251,
	   0.30103   ] ]))

	y = x.copy(); y.add(1); y.log2()
	assert(cmp_mat(y.to_array(), [ [0.        , 0.26303441, 0.48542683, 0.67807191, 0.84799691,
	   1.        ] ]))

	y = x.copy(); y.sqrt()
	assert(cmp_mat(y.to_array(), [ [0.        , 0.4472136 , 0.63245553, 0.77459667, 0.89442719,
	   1.        ] ]))

	y = x.copy(); y.cbrt()
	assert(cmp_mat(y.to_array(), [ [0.        , 0.58480355, 0.7368063 , 0.84343267, 0.92831777,
	   1.        ] ]))

func gdblas_exp1(a: Vector2) -> Vector2:
	return Vector2(cos(a.y), sin(a.y))

func gdblas_exp2(a: Vector2, args: Array) -> Vector2:
	return Vector2(cos(a.y), sin(a.y)) * args[0] + args[1]

func test14():
	print("test14")

	var gbl = GDBlas.new()
	var x = gbl.linspace(0, 1, 6)
	var z = gbl.new_complex_mat(x.dimension())
	z.imag(x)
	z.H()

	assert(cmp_mat(z.to_array(), [ [C(0, 0), C(0, -0.2), C(0, -0.4), C(0, -0.6), C(0, -0.8), C(0, -1)] ]))

	var y1 = z.copy(); y1.f(gdblas_exp1)
	var y2 = z.copy(); y2.exp()
	assert(y1.is_eq(y2, 1e-7))

	y1 = z.copy(); y1.f(gdblas_exp2, [ 2, C(1, -1) ])
	y2 = z.copy(); y2.exp(); y2.mul(2); y2.add(C(1, -1))
	assert(y1.is_eq(y2, 1e-6))

var A: GDBlasMat = null
func ode_fx(x: GDBlasMat, t: float):
	return A.prod(x)

func test15():
	print("test15")

	var gbl = GDBlas.new()
	var x = gbl.new_mat(5, 1)
	if not x.has_method("eval_ode"):
		print("ODE not enabled, skipping...")
		return

	x.fill(1.0)
	A = gbl.new_mat(5)
	A.eye(-0.01)

	x.eval_ode(ode_fx, 1)

	var xend = x.copy()
	xend.fill(0.99)
	assert(x.is_eq(xend, 1e-7))

func _test16_set_mat(a, i: int, j: int):
	return i + j;

func test16():
	print("test16")

	var gbl = GDBlas.new()
	var A = gbl.new_mat(4, 3)
	A.f(_test16_set_mat, null, true)
	var B = A.copy()
	B.mul(-1)

	var C = A.conv(B)
	var C1 = gbl.new_mat(C.dimension())
	assert(C1.from_array(
		[
			[0, 0, -1, -4, -4],
			[0, -2, -8, -14, -12],
			[-1, -8, -24, -32, -25],
			[-4, -20, -52, -60, -44],
			[-10, -32, -69, -68, -46],
			[-12, -34, -68, -62, -40],
			[-9, -24, -46, -40, -25]]
	) == 0)
	assert(C.is_eq(C1))

	var D = A.conv(B, true)
	var D1 = gbl.new_mat(A.dimension())
	assert(D1.from_array(
		[ [-8, -24, -32], [-20, -52, -60], [-32, -69, -68], [-34, -68, -62] ]
	) == 0)
	assert(D.is_eq(D1))

	assert(cmp_mat(D.conv(C1, true).to_array(), [
		[2408, 6076, 9424],
		[5756, 12696, 17824],
		[10126, 20168, 26178],
		[13552, 24937, 30392]
	]))

func test17():
	print("test17")

	var gbl = GDBlas.new()
	var A = gbl.new_mat(3, 2)
	A.from_array([
		[ 1, 2 ],
		[ 3, 4 ],
		[ 5, 6 ]
	])
	var packA: PackedFloat64Array = A.pack()
	var A2 = gbl.new_mat(A.dimension())
	A2.unpack(packA)
	assert(A.is_eq(A2))

	var B = gbl.new_complex_mat(3, 2)
	B.real(A)
	B.imag(A)
	B.H()
	var packB: PackedFloat64Array = B.pack()
	assert(packB == PackedFloat64Array([1, -1, 3, -3, 5, -5, 2, -2, 4, -4, 6, -6]))
	var packBre: PackedFloat64Array = B.pack(GDBlas.REAL_COMPONENT)
	assert(packBre == PackedFloat64Array([1, 3, 5, 2, 4, 6]))
	var packBim: PackedFloat64Array = B.pack(GDBlas.IMAG_COMPONENT)
	assert(packBim == PackedFloat64Array([-1, -3, -5, -2, -4, -6]))

	var B2 = gbl.new_complex_mat(B.dimension())
	B2.unpack(packB)
	assert(B.is_eq(B2))
	var B2re = gbl.new_mat(B.dimension())
	B2re.unpack(packBre)
	assert(B.real().is_eq(B2re))
	var B2im = gbl.new_mat(B.dimension())
	B2im.unpack(packBim)
	assert(B.imag().is_eq(B2im))
	B2.fill(0)
	B2.unpack(packBre, GDBlas.REAL_COMPONENT)
	B2.unpack(packBre, GDBlas.IMAG_COMPONENT)
	B2.T()
	B.H()
	assert(B.is_eq(B2))

func test18():
	print("test18")

	var gbl = GDBlas.new()
	var A = gbl.new_mat(3, 2)
	A.from_array([
		[ 1, 2 ],
		[ 3, 4 ],
		[ 5, 6 ]
	])
	var packA: PackedFloat32Array = PackedFloat32Array(Array(range(1, 7)))
	var A2 = gbl.new_mat(A.dimension())
	A2.unpack(packA)
	assert(A2.is_eq(A))

func test19():
	print("test19")

	var gbl = GDBlas.new()
	var A = gbl.new_mat(3, 2)
	A.from_array([
		[ 1, 2 ],
		[ 3, 4 ],
		[ 5, 6 ]
	])
	var packA: PackedByteArray = PackedByteArray(Array(range(1, 7)))
	var A2 = gbl.new_mat(A.dimension())
	A2.unpack(packA)
	assert(A2.is_eq(A))

func test20():
	print("test20")

	var gbl = GDBlas.new()
	var A = gbl.new_mat(3, 2)
	A.from_array([
		[ 1, 2 ],
		[ 3, 4 ],
		[ 5, 6 ]
	])
	var packA: PackedByteArray = [ 1, -1, 2, -2, 3, -3, 4, -4, 5, -5, 6, -6 ]
	var A2 = gbl.new_mat(A.dimension())
	A2.unpack(packA, GDBlas.REAL_COMPONENT, 2)
	assert(A2.is_eq(A))

func _test21_set_mat(a: float, i: int, j: int):
	return i + j + a + 1;

func test21():
	print("test21")
	const dim = Vector2i(3, 5)
	var gbl = GDBlas.new()
	var R = gbl.new_mat(dim)
	R.f(_test21_set_mat, null, true)
	var G = gbl.new_mat(dim)
	G.f(_test21_set_mat, null, true)
	G.f(_test21_set_mat, null, true)
	var B = gbl.new_mat(dim)
	B.f(_test21_set_mat, null, true)
	B.f(_test21_set_mat, null, true)
	B.f(_test21_set_mat, null, true)

	var pack: PackedByteArray = gbl.mat_to_image_data([ R, G, B ])

	var Ra = R.to_array()
	var Ga = G.to_array()
	var Ba = B.to_array()
	var pack2: PackedByteArray = PackedByteArray()
	for i in range(len(Ra)):
		for j in range(len(Ra[0])):
			var r = Ra[i][j]
			var g = Ga[i][j]
			var b = Ba[i][j]
			pack2.append(r)
			pack2.append(g)
			pack2.append(b)

	assert(pack == pack2)

func _test22_set_mat(a, args: Array, i: int, j: int):
	if args[0] or (i % 3 == 0 and j % 2 == 0):
		if args[1]:
			return Vector2(2, -2)
		else:
			return 2

	if args[1]:
		return Vector2(3, -3)
	else:
		return 2

func test22():
	print("test22")

	var gbl = GDBlas.new()
	var size = Vector2i(72, 50)
	var A = gbl.new_mat(size)
	A.f(_test22_set_mat, [ false, false ], true)
	var Ad = A.downsample(3, 2)
	var Ad2 = gbl.new_mat(Ad.dimension())
	Ad2.f(_test22_set_mat, [ true, false ], true)
	assert(Ad.is_eq(Ad2))

	var B = gbl.new_complex_mat(size)
	B.f(_test22_set_mat, [ false, true ], true)
	var Bd = B.downsample(3, 2)
	var Bd2 = gbl.new_complex_mat(Bd.dimension())
	Bd2.f(_test22_set_mat, [ true, true ], true)
	assert(Bd.is_eq(Bd2))

	var C = gbl.new_complex_mat(11, 13)
	var Cd = C.downsample(3, 2)
	assert(Cd.dimension() == Vector2i(3, 6))

func test23():
	print("test23")

	var gbl = GDBlas.new()
	var A = gbl.new_mat(6)
	A.fill(1)
	var filter = gbl.new_mat()
	filter.from_array([
		[ 0, 1, 0 ],
		[ 1, 1, 1],
		[ 0, 1, 0],
	])
	filter.div(filter.integrate())

	var Ad = A.downsample(2, 2, filter)
	assert(cmp_mat(Ad.to_array(), [[0.6, 0.8, 0.8], [0.8, 1, 1], [0.8, 1, 1]]))

	var B = gbl.new_complex_mat(6)
	B.real(A)
	B.imag(A)
	B.H()
	var Bd = B.downsample(2, 2, filter)
	var Bdr = Bd.real()
	var Bdi = Bd.imag()
	Bdi.mul(-1)
	assert(cmp_mat(Bdr.to_array(), [[0.6, 0.8, 0.8], [0.8, 1, 1], [0.8, 1, 1]]))
	assert(cmp_mat(Bdi.to_array(), [[0.6, 0.8, 0.8], [0.8, 1, 1], [0.8, 1, 1]]))

func test24():
	print("test24")

	var gbl = GDBlas.new()
	var polygon: PackedVector2Array = PackedVector2Array()
	polygon.append(Vector2(0, 0))
	polygon.append(Vector2(0, 7))
	polygon.append(Vector2(4, 2))
	polygon.append(Vector2(2, 0))
	polygon.append(Vector2(0, 0))
	var area: float = gbl.area(polygon)
	assert(cmp_scalar(area, 16))

	polygon.clear()
	polygon.append(Vector2(0, 0))
	polygon.append(Vector2(0, 5))
	polygon.append(Vector2(5, 5))
	polygon.append(Vector2(5, 0))
	polygon.append(Vector2(0, 0))

	area = gbl.area(polygon)
	assert(cmp_scalar(area, 25))

var _line_buffer: Array;
var _poly_colors: PackedColorArray;
func test25():
	print("test25")

	var gbl = GDBlas.new()
	var polygon: PackedVector2Array = PackedVector2Array()
	polygon.append(Vector2(30, 30))
	polygon.append(Vector2(30, 600))
	polygon.append(Vector2(1080, 600))
	polygon.append(Vector2(1080, 30))
	_line_buffer = gbl.buffer(polygon, 15, 16, 16, 16)

	var buffer: PackedVector2Array = PackedVector2Array([
		C(1065, 585), C(45, 585), C(45, 30), C(43.85819, 24.25975),
		C(40.6066, 19.3934), C(35.74025, 16.14181), C(30, 15), C(24.25975, 16.14181),
		C(19.3934, 19.3934), C(16.14181, 24.25975), C(15, 30), C(15, 600),
		C(16.14181, 605.7402), C(19.3934, 610.6066), C(24.25975, 613.8582),
		C(30, 615), C(1080, 615), C(1085.74, 613.8582), C(1090.607, 610.6066),
		C(1093.858, 605.7402), C(1095, 600), C(1095, 30), C(1093.858, 24.25975),
		C(1090.607, 19.3934), C(1085.74, 16.14181), C(1080, 15), C(1074.26, 16.14181),
		C(1069.393, 19.3934), C(1066.142, 24.25975), C(1065, 30), C(1065, 585)
	])

	assert(cmp_vec2_array(buffer, _line_buffer[0], 1e-3))

	_poly_colors.resize(_line_buffer[0].size())
	for i in range(_poly_colors.size()):
		_poly_colors[i] = Color.from_hsv((randi() % 20) / 30.0, 0.75, 0.5)

func test26():
	print("test26")

	var gbl = GDBlas.new()
	var polygon: Array = [
		PackedVector2Array([C(2, 1.3), C(2.4, 1.7), C(2.8, 1.8), C(3.4, 1.2), C(3.7, 1.6),
							C(3.4, 2), C(4.1, 3), C(5.3, 2.6), C(5.4, 1.2), C(4.9, 0.8),
							C(2.9, 0.7), C(2, 1.3)]),
		PackedVector2Array([C(4.0, 2.0), C(4.2, 1.4), C(4.8, 1.9), C(4.4, 2.2), C(4.0, 2.0)])
	]

	var C: Vector2 = gbl.centroid(polygon)
	assert(cmp_scalar(C, Vector2(4.04663, 1.6349), 1e-4))

func test27():
	print("test27")

	var gbl = GDBlas.new()
	var polygon: Array = [
		PackedVector2Array([
			C(2, 1.3), C(2.4, 1.7), C(2.8, 1.8), C(3.4, 1.2), C(3.7, 1.6),
			C(3.4, 2), C(4.1, 3), C(5.3, 2.6), C(5.4, 1.2), C(4.9, 0.8), C(2.9, 0.7), C(2, 1.3)
		]),
		PackedVector2Array([ C(4.0, 2.0), C(4.2, 1.4), C(4.8, 1.9), C(4.4, 2.2), C(4.0, 2.0) ])
	]
	var line: PackedVector2Array = PackedVector2Array([ C(5, 0), C(2, 0) ])
	var point: Vector2 = Vector2(4.3, 1.9)

	var polygon2: Array = [
		PackedVector2Array([ C(3, 3), C(4, 3), C(4, 4), C(3, 4) ])
	]

	var closest = gbl.closest_points(polygon, point)
	assert(cmp_scalar(closest[0], C(4.2, 2.1), 1e-4))
	assert(cmp_scalar(closest[1], C(4.3, 1.9), 1e-4))

	closest = gbl.closest_points(point, polygon)
	assert(cmp_scalar(closest[0], C(4.3, 1.9), 1e-4))
	assert(cmp_scalar(closest[1], C(4.2, 2.1), 1e-4))

	closest = gbl.closest_points(polygon, line)
	assert(cmp_scalar(closest[0], C(2.9, 0.7), 1e-4))
	assert(cmp_scalar(closest[1], C(2.9, 0), 1e-4))

	closest = gbl.closest_points(line, polygon)
	assert(cmp_scalar(closest[0], C(2.9, 0), 1e-4))
	assert(cmp_scalar(closest[1], C(2.9, 0.7), 1e-4))

	polygon = [
		PackedVector2Array([ C(1, 1), C(2, 1), C(2, 2), C(1, 2) ])
	]

	closest = gbl.closest_points(polygon, polygon2)
	assert(cmp_scalar(closest[0], C(2, 2), 1e-4))
	assert(cmp_scalar(closest[1], C(3, 3), 1e-4))

	closest = gbl.closest_points(polygon2, polygon)
	assert(cmp_scalar(closest[0], C(3, 3), 1e-4))
	assert(cmp_scalar(closest[1], C(2, 2), 1e-4))

func test28():
	print("test28")

	var gbl = GDBlas.new()
	var polgon: Array = [
		[ C(2, 1.3), C(2.4, 1.7), C(2.8, 1.8), C(3.4, 1.2), C(3.7, 1.6),
		  C(3.4, 2), C(4.1, 3), C(5.3, 2.6), C(5.4, 1.2), C(4.9, 0.8), C(2.9, 0.7), C(2, 1.3) ]
	]

	var hull: PackedVector2Array = gbl.convex_hull(polgon)
	var hull2: PackedVector2Array = PackedVector2Array([
		C(2, 1.3), C(2.4, 1.7), C(4.1, 3), C(5.3, 2.6), C(5.4, 1.2), C(4.9, 0.8), C(2.9, 0.7), C(2, 1.3)
	])

	assert(cmp_vec2_array(hull, hull2))

func test29():
	print("test29")

	var gbl = GDBlas.new()
	var polygon: Array = [
		[ C(0, 0), C(10, 10), C(0, 9) ],
		[ C(1, 2), C(4, 6), C(2, 8), C(1, 2) ],
	]

	var area: float = gbl.area(polygon)
	assert(cmp_scalar(area, -7))

	var corrected_polygon: Array = gbl.correct(polygon)

	area = gbl.area(corrected_polygon)
	assert(cmp_scalar(area, 38))

func test30():
	print("test30")

	var gbl = GDBlas.new()
	var polygon1: Array = [
		[ C(0, 2), C(0, 3), C(2, 4), C(1, 2), C(0, 2) ]
	]
	var polygon2: Array = [
		[ C(0, 4), C(3, 4), C(2, 2), C(0, 1), C(0, 4) ]
	]
	var polygon3: Array = [
		[ C(-1, -1), C(-3, -4), C(-7, -7), C(-4, -3), C(-1, -1) ]
	]

	assert(gbl.covered_by(polygon1, polygon2) == 1)
	assert(gbl.covered_by(polygon1, polygon3) == 0)
	assert(gbl.covered_by(polygon1, polygon1) == 1)

func test31():
	print("test31")

	var gbl = GDBlas.new()
	var polygon: Array = [
		[ C(0, 0), C(0, 3), C(3, 3), C(3, 0), C(0, 0) ]
	]
	var line1: PackedVector2Array = PackedVector2Array([ C(1, 1), C(2, 2), C(4, 4) ])
	var line2: PackedVector2Array = PackedVector2Array([ C(1, 1), C(1, 2), C(1, 3) ])

	assert(gbl.crosses(polygon, line1) == 1)
	assert(gbl.crosses(line1, polygon) == 1)
	assert(gbl.crosses(polygon, line2) == 0)

func test32():
	print("test32")

	var gbl = GDBlas.new()
	var polygon: Array = [
		PackedVector2Array([ C(0, 0), C(0, 10), C(10, 10), C(10, 0), C(0, 0) ]),
		PackedVector2Array([ C(1, 1), C(4, 1), C(4, 4), C(1, 4), C(1, 1) ]),
	]

	var densified_polygon: Array = gbl.densify(polygon, 6)

	var densified_polygon2: Array = [
		PackedVector2Array([
			C(0, 0), C(0, 5), C(0, 10), C(5, 10), C(10, 10), C(10, 5), C(10, 0), C(5, 0), C(0, 0)
		]),
		PackedVector2Array([ C(1, 1), C(4, 1), C(4, 4), C(1, 4), C(1, 1) ]),
	]

	assert(cmp_vec2_array(densified_polygon[0], densified_polygon2[0]))
	assert(cmp_vec2_array(densified_polygon[1], densified_polygon2[1]))

	var ring: PackedVector2Array = PackedVector2Array([ C(0, 0), C(0, 10), C(10, 10), C(10, 0), C(0, 0) ])
	var densified_ring: PackedVector2Array = gbl.densify(ring, 6)

	var densified_ring2: PackedVector2Array = PackedVector2Array([
		C(0, 0), C(0, 5), C(0, 10), C(5, 10), C(10, 10), C(10, 5), C(10, 0), C(5, 0), C(0, 0)
	])

	assert(cmp_vec2_array(densified_ring, densified_ring2))

	var line: PackedVector2Array = PackedVector2Array([ C(0, 0), C(0, 10), C(10, 10), C(10, 0) ])
	var densified_line: PackedVector2Array = gbl.densify(line, 6)

	var densified_line2: PackedVector2Array = PackedVector2Array([
		C(0, 0), C(0, 5), C(0, 10), C(5, 10), C(10, 10), C(10, 5), C(10, 0), C(5, 0)
	])

	assert(cmp_vec2_array(densified_line, densified_line2))

func test33():
	print("test33")

	var gbl = GDBlas.new()
	var polygon1: Array = [
		PackedVector2Array([
			C(2, 1.3), C(2.4, 1.7), C(2.8, 1.8), C(3.4, 1.2), C(3.7, 1.6), C(3.4, 2),
			C(4.1, 3), C(5.3, 2.6), C(5.4, 1.2), C(4.9, 0.8), C(2.9, 0.7), C(2, 1.3)
		]),
		PackedVector2Array([ C(4.0, 2.0), C(4.2, 1.4), C(4.8, 1.9), C(4.4, 2.2), C(4.0, 2.0) ])
	]
	var polygon2: Array = [ PackedVector2Array([
		C(4.0, -0.5), C(3.5, 1.0), C(2.0, 1.5), C(3.5, 2.0), C(4.0, 3.5),
		C(4.5, 2.0), C(6.0, 1.5), C(4.5, 1.0), C(4.0, -0.5)
	]) ]

	var difference: Array = gbl.difference(polygon1, polygon2)

	var areas: Dictionary = {
		0: 0.02375,
		1: 0.542951,
		2: 0.0149697,
		3: 0.226855,
		4: 0.839424,
	}

	var i: int = 0;
	for diff in difference:
		var area: float = gbl.area(diff)
		assert(cmp_scalar(areas[i], area, 1e-4))
		i += 1

	difference = gbl.difference(polygon2, polygon1)

	areas = {
		0: 0.525154,
		1: 0.015,
		2: 0.181136,
		3: 0.340083,
		4: 0.128798,
		5: 0.307778,
	}

	i = 0;
	for diff in difference:
		var area: float = gbl.area(diff)
		assert(cmp_scalar(areas[i], area, 1e-4))
		i += 1

func test34():
	print("test34")

	var gbl = GDBlas.new()
	var line1: PackedVector2Array = PackedVector2Array([ C(0, 0), C(1, 1), C(1, 2), C(2, 1), C(2, 2) ])
	var line2: PackedVector2Array = PackedVector2Array([ C(1, 0), C(0, 1), C(1, 1), C(2, 1), C(3, 1) ])

	var dist: float = gbl.discrete_frechet_distance(line1, line2)
	assert(cmp_scalar(dist, 1.41421, 1e-4))

func test35():
	print("test35")

	var gbl = GDBlas.new()
	var line1: PackedVector2Array = PackedVector2Array([ C(0, 0), C(1, 1), C(1, 2), C(2, 1), C(2, 2) ])
	var line2: PackedVector2Array = PackedVector2Array([ C(1, 0), C(0, 1), C(1, 1), C(2, 1), C(3, 1) ])

	var dist: float = gbl.discrete_hausdorff_distance(line1, line2)
	assert(cmp_scalar(dist, 1.0, 1e-4))

func test36():
	print("test36")

	var gbl = GDBlas.new()

	var polygon1: Array = [
		PackedVector2Array([ C(0,  2),  C(-2,  0),  C(-4,  2),  C(-2,  4),  C(0,  2) ])
	]
	var polygon2: Array = [
		PackedVector2Array([ C(2,  2),  C(4,  4),  C(6,  2),  C(4,  0),  C(2,  2) ])
	]
	var polygon3: Array = [
		PackedVector2Array([ C(0,  2),  C(2,  4),  C(4,  2),  C(2,  0),  C(0,  2) ])
	]

	var ring: PackedVector2Array = PackedVector2Array([ C(2,  2),  C(4,  4),  C(6,  2),  C(4,  0),  C(2,  2) ])
	var point: Vector2 = Vector2(-2, 2)

	assert(gbl.disjoint(polygon1, polygon2) == 1)
	assert(gbl.disjoint(polygon1, polygon3) == 0)
	assert(gbl.disjoint(polygon1, ring) == 1)
	assert(gbl.disjoint(ring, polygon1) == 1)
	assert(gbl.disjoint(polygon1, point) == 0)
	assert(gbl.disjoint(point, polygon1) == 0)

func test37():
	print("test37")

	var gbl = GDBlas.new()
	var line: PackedVector2Array = PackedVector2Array()
	for x in range(0.0, 5.0, 1.0):
		for y in range(0.0, 5.0, 1.0):
			line.append(Vector2(x, y))

	var point: Vector2 = Vector2(1.4, 2.6)
	var dist: float = gbl.comparable_distance(line, point)

	assert(cmp_scalar(dist, 0.00235, 1e-4))

func test38():
	print("test38")

	var gbl = GDBlas.new()
	var polygon: Array = [
		PackedVector2Array([
			C(2, 1.3), C(2.4, 1.7), C(2.8, 1.8), C(3.4, 1.2), C(3.7, 1.6),
			C(3.4, 2), C(4.1, 3), C(5.3, 2.6), C(5.4, 1.2), C(4.9, 0.8), C(2.9, 0.7), C(2, 1.3)
		]),
		PackedVector2Array([ C(4.0, 2.0), C(4.2, 1.4), C(4.8, 1.9), C(4.4, 2.2), C(4.0, 2.0) ])
	]
	var line: PackedVector2Array = PackedVector2Array([ C(0, 0), C(0, 3) ])
	var point: Vector2 = Vector2(1, 2)

	var dist: float = gbl.distance(point, polygon)
	assert(cmp_scalar(dist, 1.22066, 1e-5))

	dist = gbl.distance(point, line)
	assert(cmp_scalar(dist, 1, 1e-5))

func test39():
	print("test39")

	var gbl = GDBlas.new()
	var polygon: Array = [
		PackedVector2Array([
			C(2, 1.3), C(2.4, 1.7), C(2.8, 1.8), C(3.4, 1.2), C(3.7, 1.6), C(3.4, 2),
			C(4.1, 3), C(5.3, 2.6), C(5.4, 1.2), C(4.9, 0.8), C(2.9, 0.7), C(2, 1.3)
		]),
		PackedVector2Array([ C(4.0, 2.0), C(4.2, 1.4), C(4.8, 1.9), C(4.4, 2.2), C(4.0, 2.0) ])
	]

	var envelope: Rect2 = gbl.envelope(polygon)
	assert(cmp_scalar(envelope.position, C(2, 0.7)))
	assert(cmp_scalar(envelope.size, C(3.4, 2.3)))

func test40():
	print("test40")

	var gbl = GDBlas.new()
	var ring1: PackedVector2Array = PackedVector2Array([ C(0, 0), C(0, 5), C(5, 5), C(5, 0), C(0, 0) ])
	var ring2: PackedVector2Array = PackedVector2Array([ C(5, 0), C(0, 0), C(0, 5), C(5, 5), C(5, 0) ])

	assert(gbl.equals(ring1, ring2) == 1)

func test41():
	print("test41")

	var gbl = GDBlas.new()
	var polygon1: Array = [
		PackedVector2Array([
			C(2, 1.3), C(2.4, 1.7), C(2.8, 1.8), C(3.4, 1.2), C(3.7, 1.6), C(3.4, 2),
			C(4.1, 3), C(5.3, 2.6), C(5.4, 1.2), C(4.9, 0.8), C(2.9, 0.7), C(2, 1.3)
		]),
		PackedVector2Array([ C(4.0, 2.0), C(4.2, 1.4), C(4.8, 1.9), C(4.4, 2.2), C(4.0, 2.0) ])
	]
	var polygon2: Array = [ PackedVector2Array([
		C(4.0, -0.5), C(3.5, 1.0), C(2.0, 1.5), C(3.5, 2.0), C(4.0, 3.5),
		C(4.5, 2.0), C(6.0, 1.5), C(4.5, 1.0), C(4.0, -0.5)
	]) ]

	var intersection: Array = gbl.intersection(polygon1, polygon2)

	var areas: Dictionary = {
		0: 2.50205,
	}

	var i: int = 0;
	for intr in intersection:
		var area: float = gbl.area(intr)
		assert(cmp_scalar(areas[i], area, 1e-4))
		i += 1

	intersection = gbl.intersection(polygon2, polygon1)

	areas = {
		0: 2.50205,
	}

	i = 0;
	for intr in intersection:
		var area: float = gbl.area(intr)
		assert(cmp_scalar(areas[i], area, 1e-4))
		i += 1

func test42():
	print("test42")

	var gbl = GDBlas.new()
	var line1: PackedVector2Array = PackedVector2Array([ C(1, 1), C(2, 2), C(3, 3) ])
	var line2: PackedVector2Array = PackedVector2Array([ C(2, 1), C(1, 2), C(4, 0) ])
	var line3: PackedVector2Array = PackedVector2Array([ C(0.5, 0.5), C(0.5, 2.5), C(4, 4) ])

	assert(gbl.intersects(line1, line2) == 1)
	assert(gbl.intersects(line1, line3) == 0)

func test43():
	print("test43")

	var gbl = GDBlas.new()
	var polygon: Array = [
		PackedVector2Array([
			C(2, 1.3), C(2.4, 1.7), C(2.8, 1.8), C(3.4, 1.2), C(3.7, 1.6), C(3.4, 2),
			C(4.1, 3), C(5.3, 2.6), C(5.4, 1.2), C(4.9, 0.8), C(2.9, 0.7), C(2, 1.3)
		]),
		PackedVector2Array([ C(1.9, 1.3), C(4.2, 1.4), C(4.8, 1.9), C(4.4, 2.2), C(4.0, 2.0) ])
	]
	var ring: PackedVector2Array = PackedVector2Array([ C(0, 0), C(10, 0), C(10, 10), C(0, 0), C(-10, 0) ])

	assert(gbl.is_simple(polygon) == 1)
	assert(gbl.is_simple(ring) == 0)

func test44():
	print("test44")

	var gbl = GDBlas.new()
	var polygon: Array = [
		PackedVector2Array([ C(0, 0), C(0, 10), C(10, 10), C(10, 0) ]),
		PackedVector2Array([ C(0, 0), C(9, 2), C(9, 1), C(0, 0) ]),
		PackedVector2Array([ C(0, 0), C(2, 9), C(1, 9), C(0, 0) ]),
	]

	var failure: int = gbl.is_valid(polygon)
	assert(failure != 0)
	var could_be_fixed: bool = failure == 20 or failure == 22
	assert(could_be_fixed)
	if could_be_fixed:
		polygon = gbl.correct(polygon)
		failure = gbl.is_valid(polygon)
		assert(failure == 0)

func test45():
	print("test45")

	var gbl = GDBlas.new()
	var line: PackedVector2Array = PackedVector2Array([ C(0, 0), C(1, 1), C(4, 8), C(3, 2) ])

	var len: float = gbl.length(line)
	assert(cmp_scalar(len, 15.1127, 1e-4))

func test46():
	print("test46")

	var gbl = GDBlas.new()
	var polygon1: Array = [ PackedVector2Array([ C(0, 0), C(0, 4), C(4, 4), C(4, 0), C(0, 0) ]) ]
	var polygon2: Array = [ PackedVector2Array([ C(2, 2), C(2, 6), C(6, 7), C(6, 1), C(2, 2) ]) ]
	var polygon3: Array = [ PackedVector2Array([ C(-1, -1), C(-3, -4), C(-7, -7), C(-4, -3), C(-1, -1) ]) ]

	assert(gbl.overlaps(polygon1, polygon2) == 1)
	assert(gbl.overlaps(polygon1, polygon3) == 0)

func test47():
	print("test47")

	var gbl = GDBlas.new()
	var polygon: Array = [ PackedVector2Array([ C(0, 0), C(3, 4), C(5, -5), C(-2, -4), C(0, 0) ]) ]

	var perim: float = gbl.perimeter(polygon)
	assert(cmp_scalar(perim, 25.7627, 1e-4))

func test48():
	print("test48")

	var gbl = GDBlas.new()
	var polygon: Array = [
		PackedVector2Array([C(2, 1.3), C(2.4, 1.7), C(2.8, 1.8), C(3.4, 1.2), C(3.7, 1.6),
							C(3.4, 2), C(4.1, 3), C(5.3, 2.6), C(5.4, 1.2), C(4.9, 0.8),
							C(2.9, 0.7), C(2, 1.3)]),
		PackedVector2Array([C(4.0, 2.0), C(4.2, 1.4), C(4.8, 1.9), C(4.4, 2.2), C(4.0, 2.0)])
	]

	var rel: String = gbl.relation(Vector2(4, 1), polygon)
	assert(rel == "0FFFFF212")

func test49():
	print("test49")

	var gbl = GDBlas.new()
	var polygon: Array = [
		PackedVector2Array([ C(0, 0), C(0, 9), C(10, 10), C(0, 0) ]),
		PackedVector2Array([ C(1, 2), C(4, 6), C(2, 8), C(1, 2) ]),
	]
	var ring: PackedVector2Array = PackedVector2Array([ C(0, 0), C(8, 8), C(0, 9), C(0, 0) ])

	var area: float = gbl.area(polygon)
	assert(cmp_scalar(area, 38))
	area = gbl.area(ring)
	assert(cmp_scalar(area, -36))

	polygon = gbl.reverse(polygon)
	area = gbl.area(polygon)
	assert(cmp_scalar(area, -38))
	ring = gbl.reverse(ring)
	area = gbl.area(ring)
	assert(cmp_scalar(area, 36))

func test50():
	print("test50")

	var gbl = GDBlas.new()
	var line: PackedVector2Array = PackedVector2Array([
		C(1.1, 1.1), C(2.5, 2.1), C(3.1, 3.1), C(4.9, 1.1), C(3.1, 1.9)
	])
	var simple_line: PackedVector2Array = PackedVector2Array([
		C(1.1, 1.1), C(3.1, 3.1), C(4.9, 1.1), C(3.1, 1.9)
	])

	line = gbl.simplify(line, 0.5)
	assert(cmp_vec2_array(line, simple_line))

func test51():
	print("test51")

	var gbl = GDBlas.new()
	var polygon1: Array = [
		PackedVector2Array([
			C(2, 1.3), C(2.4, 1.7), C(2.8, 1.8), C(3.4, 1.2), C(3.7, 1.6), C(3.4, 2),
			C(4.1, 3), C(5.3, 2.6), C(5.4, 1.2), C(4.9, 0.8), C(2.9, 0.7), C(2, 1.3)
		]),
		PackedVector2Array([ C(4.0, 2.0), C(4.2, 1.4), C(4.8, 1.9), C(4.4, 2.2), C(4.0, 2.0) ])
	]
	var polygon2: Array = [ PackedVector2Array([
		C(4.0, -0.5), C(3.5, 1.0), C(2.0, 1.5), C(3.5, 2.0), C(4.0, 3.5),
		C(4.5, 2.0), C(6.0, 1.5), C(4.5, 1.0), C(4.0, -0.5)
	]) ]

	var difference: Array = gbl.sym_difference(polygon1, polygon2)

	var areas: Dictionary = {
		0: 0.542951,
		1: 0.525154,
		2: 0.0149697,
		3: 0.181136,
		4: 0.226855,
		5: 0.839424,
		6: 0.02375,
		7: 0.340083,
		8: 0.015,
		9: 0.128798,
		10: 0.307778,
	}

	var i: int = 0;
	for diff in difference:
		var area: float = gbl.area(diff)
		assert(cmp_scalar(areas[i], area, 1e-4))
		i += 1

func test52():
	print("test52")

	var gbl = GDBlas.new()
	var polygon: Array = [
		PackedVector2Array([
			C(0, 0), C(0, 3), C(2, 3), C(2, 2), C(1, 2), C(1, 1),
			C(2, 1), C(2, 2), C(3, 2), C(3, 0), C(0, 0)
		])
	]

	assert(gbl.touches(polygon) == 1)

	var polygon1: Array = [
		PackedVector2Array([ C(0, 0), C(0, 4), C(4, 4), C(4, 0), C(0, 0) ])
	]
	var polygon2: Array = [
		PackedVector2Array([ C(0, 0), C(0, -4), C(-4, -4), C(-4, 0), C(0, 0) ])
	]
	var polygon3: Array = [
		PackedVector2Array([ C(1, 1), C(0, -4), C(-4, -4), C(-4, 0), C(1, 1) ])
	]

	assert(gbl.touches(polygon1, polygon2) == 1)
	assert(gbl.touches(polygon1, polygon3) == 0)
	assert(gbl.touches(polygon2, polygon3) == 0)


func test53():
	print("test53")

	var gbl = GDBlas.new()
	var polygon1: Array = [
		PackedVector2Array([
			C(2, 1.3), C(2.4, 1.7), C(2.8, 1.8), C(3.4, 1.2), C(3.7, 1.6), C(3.4, 2),
			C(4.1, 3), C(5.3, 2.6), C(5.4, 1.2), C(4.9, 0.8), C(2.9, 0.7), C(2, 1.3)
		]),
		PackedVector2Array([ C(4.0, 2.0), C(4.2, 1.4), C(4.8, 1.9), C(4.4, 2.2), C(4.0, 2.0) ])
	]
	var polygon2: Array = [ PackedVector2Array([
		C(4.0, -0.5), C(3.5, 1.0), C(2.0, 1.5), C(3.5, 2.0), C(4.0, 3.5),
		C(4.5, 2.0), C(6.0, 1.5), C(4.5, 1.0), C(4.0, -0.5)
	]) ]

	var union: Array = gbl.union_(polygon1, polygon2)

	var areas: Dictionary = {
		0: 5.64795,
	}

	var i: int = 0;
	for u in union:
		var area: float = gbl.area(u)
		assert(cmp_scalar(areas[i], area, 1e-4))
		i += 1

func test54():
	print("test54")

	var gbl = GDBlas.new()
	var polygon: Array = [ PackedVector2Array([
		C(0, 0), C(0, 0), C(0, 5), C(5, 5), C(5, 5), C(5, 5),
		C(5, 0), C(5, 0), C(0, 0), C(0, 0), C(0, 0), C(0, 0)
	]) ]

	var unique: Array = gbl.unique(polygon)
	var unique2: Array = [ PackedVector2Array([ C(0, 0), C(0, 5), C(5, 5), C(5, 0), C(0, 0) ]) ]
	assert(cmp_vec2_array(unique[0], unique2[0]))

func test55():
	print("test55")

	var gbl = GDBlas.new()
	var polygon: Array = [
		PackedVector2Array([C(2, 1.3), C(2.4, 1.7), C(2.8, 1.8), C(3.4, 1.2), C(3.7, 1.6),
							C(3.4, 2), C(4.1, 3), C(5.3, 2.6), C(5.4, 1.2), C(4.9, 0.8),
							C(2.9, 0.7), C(2, 1.3)]),
		PackedVector2Array([C(4.0, 2.0), C(4.2, 1.4), C(4.8, 1.9), C(4.4, 2.2), C(4.0, 2.0)])
	]

	assert(gbl.within(Vector2(4, 1), polygon) == 1)
	assert(gbl.within(Vector2(40, 1), polygon) == 0)

func _process(delta):
	pass

func _draw():
	if _has_geometry_funcs:
		draw_polygon(_line_buffer[0], _poly_colors)
		draw_polyline(_line_buffer[0], Color.SADDLE_BROWN, 6)
