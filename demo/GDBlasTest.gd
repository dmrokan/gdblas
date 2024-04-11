extends Node

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
	test01()
	test02()
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

# Called every frame. 'delta' is the elapsed time since the previous frame.
func _process(delta):
	pass

