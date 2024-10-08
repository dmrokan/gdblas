#include "gdblas_mat.h"

#if defined(GDBLAS_WITH_ODE)
#define BOOST_NO_EXCEPTIONS
#include <boost/numeric/odeint.hpp>
#endif

using namespace godot;

typedef GDBlasMat::scalar_t scalar_t;

void GDBlasMat::_bind_methods() {
	ClassDB::bind_method(D_METHOD("resize", "p_variant", "p_cols"), &GDBlasMat::resize, DEFVAL(-1));
	ClassDB::bind_method(D_METHOD("copy"), &GDBlasMat::copy);
	ClassDB::bind_method(D_METHOD("dimension"), &GDBlasMat::size);
	ClassDB::bind_method(D_METHOD("get", "p_i", "p_j", "p_m", "p_n"),
						 &GDBlasMat::get, DEFVAL(-1), DEFVAL(-1));
	ClassDB::bind_method(D_METHOD("set", "p_value", "p_i", "p_j"), &GDBlasMat::set, DEFVAL(-1), DEFVAL(-1));
	ClassDB::bind_method(D_METHOD("add", "p_other"), &GDBlasMat::add);
	ClassDB::bind_method(D_METHOD("sub", "p_other"), &GDBlasMat::sub);
	ClassDB::bind_method(D_METHOD("T"), &GDBlasMat::transpose);
	ClassDB::bind_method(D_METHOD("H"), &GDBlasMat::hermitian);
	ClassDB::bind_method(D_METHOD("is_eq", "p_other", "p_eps", "p_norm_type"),
						 &GDBlasMat::is_eq, DEFVAL(GDBlasMat::EPS), DEFVAL(GDBlasMat::NORM_1));
	ClassDB::bind_method(D_METHOD("mul", "p_other"), &GDBlasMat::mul);
	ClassDB::bind_method(D_METHOD("div", "p_other"), &GDBlasMat::div);
	ClassDB::bind_method(D_METHOD("fill", "p_value"), &GDBlasMat::fill);
	ClassDB::bind_method(D_METHOD("eye", "p_value"), &GDBlasMat::eye, DEFVAL(1.0));
	ClassDB::bind_method(D_METHOD("reset"), &GDBlasMat::reset);
	ClassDB::bind_method(D_METHOD("conj"), &GDBlasMat::conj);
	ClassDB::bind_method(D_METHOD("real", "p_matrix"), &GDBlasMat::real, DEFVAL(Variant()));
	ClassDB::bind_method(D_METHOD("imag", "p_matrix"), &GDBlasMat::imag, DEFVAL(Variant()));
	ClassDB::bind_method(D_METHOD("prod", "p_other"), &GDBlasMat::prod);
	ClassDB::bind_method(D_METHOD("inv"), &GDBlasMat::inv);
	ClassDB::bind_method(D_METHOD("to_array"), &GDBlasMat::to_array);
	ClassDB::bind_method(D_METHOD("from_array", "p_array"), &GDBlasMat::from_array);
	ClassDB::bind_method(D_METHOD("sin"), &GDBlasMat::sin);
	ClassDB::bind_method(D_METHOD("cos"), &GDBlasMat::cos);
	ClassDB::bind_method(D_METHOD("abs"), &GDBlasMat::abs);
	ClassDB::bind_method(D_METHOD("exp"), &GDBlasMat::exp);
	ClassDB::bind_method(D_METHOD("log"), &GDBlasMat::log);
	ClassDB::bind_method(D_METHOD("log10"), &GDBlasMat::log10);
	ClassDB::bind_method(D_METHOD("log2"), &GDBlasMat::log2);
	ClassDB::bind_method(D_METHOD("sqrt"), &GDBlasMat::sqrt);
	ClassDB::bind_method(D_METHOD("cbrt"), &GDBlasMat::cbrt);
	ClassDB::bind_method(D_METHOD("tan"), &GDBlasMat::tan);
	ClassDB::bind_method(D_METHOD("asin"), &GDBlasMat::asin);
	ClassDB::bind_method(D_METHOD("acos"), &GDBlasMat::acos);
	ClassDB::bind_method(D_METHOD("atan"), &GDBlasMat::atan);
	ClassDB::bind_method(D_METHOD("sinh"), &GDBlasMat::sinh);
	ClassDB::bind_method(D_METHOD("cosh"), &GDBlasMat::cosh);
	ClassDB::bind_method(D_METHOD("atanh"), &GDBlasMat::atanh);
	ClassDB::bind_method(D_METHOD("erf"), &GDBlasMat::erf);
	ClassDB::bind_method(D_METHOD("erfc"), &GDBlasMat::erfc);
	ClassDB::bind_method(D_METHOD("tgamma"), &GDBlasMat::tgamma);
	ClassDB::bind_method(D_METHOD("lgamma"), &GDBlasMat::lgamma);
	ClassDB::bind_method(D_METHOD("ceil"), &GDBlasMat::ceil);
	ClassDB::bind_method(D_METHOD("floor"), &GDBlasMat::floor);
	ClassDB::bind_method(D_METHOD("trunc"), &GDBlasMat::trunc);
	ClassDB::bind_method(D_METHOD("round"), &GDBlasMat::round);
	ClassDB::bind_method(D_METHOD("pow", "p_exponent"), &GDBlasMat::pow);
	ClassDB::bind_method(D_METHOD("integrate", "p_axis"), &GDBlasMat::integrate, DEFVAL(-1));
	ClassDB::bind_method(D_METHOD("norm", "p_norm_type"), &GDBlasMat::norm);
	ClassDB::bind_method(D_METHOD("mean", "p_axis"), &GDBlasMat::mean, DEFVAL(-1));
	ClassDB::bind_method(D_METHOD("min", "p_axis"), &GDBlasMat::min, DEFVAL(-1));
	ClassDB::bind_method(D_METHOD("max", "p_axis"), &GDBlasMat::max, DEFVAL(-1));
	ClassDB::bind_method(D_METHOD("argmin", "p_axis"), &GDBlasMat::argmin, DEFVAL(-1));
	ClassDB::bind_method(D_METHOD("argmax", "p_axis"), &GDBlasMat::argmax, DEFVAL(-1));
	ClassDB::bind_method(D_METHOD("f", "p_func", "p_args", "p_indexed"),
						 &GDBlasMat::unary_func, DEFVAL(Variant()), DEFVAL(false));
	ClassDB::bind_method(D_METHOD("conv", "p_other", "p_same"), &GDBlasMat::conv, DEFVAL(false));
	ClassDB::bind_method(D_METHOD("pack", "p_component"), &GDBlasMat::pack, DEFVAL(BOTH_COMPONENTS));
	ClassDB::bind_method(D_METHOD("unpack", "p_packed_data", "p_component", "p_step", "p_offset"),
						 &GDBlasMat::unpack, DEFVAL(BOTH_COMPONENTS), DEFVAL(1), DEFVAL(0));
	ClassDB::bind_method(D_METHOD("downsample", "p_factor_m", "p_factor_n", "p_filter"),
						 &GDBlasMat::downsample, DEFVAL(Variant()));

#if defined(GDBLAS_WITH_ODE)
	ClassDB::bind_method(D_METHOD("eval_ode", "p_f", "p_dt", "p_max_step"), &GDBlasMat::eval_ode,
						 DEFVAL(1e-2));
#endif
}

Variant GDBlasMat::resize(Variant p_m, int n) {
	int m = 0;
	int type = p_m.get_type();
	if (type == Variant::INT) {
		m = p_m;
	} else if (type == Variant::VECTOR2I) {
		Vector2i tmp = p_m;
		m = tmp.x;
		n = tmp.y;
	} else {
		return ERR_INVALID_INPUT;
	}

	if (m <= 0 || n <= 0) {
		GDBLAS_ERROR("(%dx%d) is not a valid matrix dimension", m, n);

		return ERR_INVALID_DIM;
	}

	_resize_implementation(m, n);

	return 0;
}

Variant GDBlasMat::copy() {
	if (_is_zero_dim()) {
		GDBLAS_ERROR("Matrix is dimensionless");

		return ERR_INVALID_DIM;
	}

	int error = 0;

	Ref< GDBlasMat > mat = _copy_implementation(&error);

	if (error)
		return Variant();

	return Variant(mat);
}

Variant GDBlasMat::size() {
	Dimension d = _size();

	Vector2i s(d.m, d.n);

	return s;
}

Variant GDBlasMat::get(int i, int j, int m, int n) {
	if (_is_zero_dim()) {
		GDBLAS_ERROR("Matrix is dimensionless");

		return ERR_INVALID_DIM;
	}

	Dimension d = _size();

	if (i >= (int) d.m || j >= (int) d.n || i < 0 || j < 0) {
		GDBLAS_ERROR("i >= d.m || j >= d.n || i < 0 || j < 0");

		return ERR_INVALID_INDEX;
	}

	if (m <= 0 || n <= 0) {
		if (get_type() == BLAS_COMPLEX_MATRIX) {
			complex_t c = ctx()->getc(i, j);

			return _complex_to_variant(c);
		}

		return ctx()->get(i, j);
	} else if (i + m > (int) d.m || j + n > (int) d.n) {
		GDBLAS_ERROR("Can not get submatrix of size (%dx%d) from "
					  "matrix of size (%dx%d) at index (%dx%d)", m, n,
					  d.m, d.n, i, j);

		return Variant();
	}

	int error = 0;
	Ref< GDBlasMat > mat = new_mat(m, n, get_type(), &error);
	if (error)
		return Variant();

	error = mat->ctx()->get(ctx(), i, j);
	if (error)
		return Variant();

	return Variant(mat);
}

int GDBlasMat::_set_scalar(Variant &val, int i, int j) {
	Variant::Type t = val.get_type();

	if (t != Variant::INT && t != Variant::FLOAT && t != Variant::VECTOR2)
		return ERR_INVALID_INPUT;

	if (get_type() == BLAS_COMPLEX_MATRIX) {
		complex_t c = _variant_to_complex(val);

		return ctx()->set(c, i, j);
	} else if (!_is_real_number(val)) {
		GDBLAS_ERROR("Can not set: not a number");

		return ERR_INVALID_INPUT;
	}

	return ctx()->set(val, i, j);
}

int GDBlasMat::_set_submatrix(const Variant &val, int i, int j) {
	GDBlasMat *o = _cast(val);

	if (o == nullptr || this == o) {
		GDBLAS_ERROR("o == nullptr || this == o");

		return ERR_INVALID_INPUT;
	}

	Dimension d1 = _size();
	Dimension d2 = o->_size();

	if (i + d2.m > d1.m || j + d2.n > d1.n) {
		GDBLAS_ERROR("Can not set submatrix of size (%dx%d) into "
					  "matrix of size (%dx%d) at index (%dx%d)", d2.m, d2.n,
					  d1.m, d1.n, i, j);

		return ERR_INVALID_INDEX;
	}

	return ctx()->set(o->ctx(), i, j);
}

Variant GDBlasMat::set(Variant val, int i, int j) {
	if (_is_zero_dim()) {
		GDBLAS_ERROR("Matrix is dimensionless");

		return ERR_INVALID_DIM;
	}

	i = GDBLAS_MAX(0, i);
	j = GDBLAS_MAX(0, j);

	Dimension d = _size();

	if (i >= (int) d.m || j >= (int) d.n) {
		GDBLAS_ERROR("i >= d.m || j >= d.n");

		return ERR_INVALID_INDEX;
	}

	int error = _set_scalar(val, i, j);
	if (error) {
		GDBlasMat *o = _cast(val);
		if (o == nullptr || this == o) {
			GDBLAS_ERROR("o == nullptr || this == o");

			return ERR_INVALID_INPUT;
		}

		if (_cmp_size(o))
			error = _eq_implementation(o);
		else
			error = _set_submatrix(val, i, j);
	}

	return error;
}

Variant GDBlasMat::add(Variant other) {
	if (_is_zero_dim()) {
		GDBLAS_ERROR("Matrix is dimensionless");

		return ERR_INVALID_DIM;
	}

	if (_is_real_number(other)) {
		scalar_t val = other;

		auto add_func_r = [&](scalar_t a) { return a + val; };
		auto add_func_rc = [&](const complex_t &a) { return a + val; };

		int type = get_type();
		if (type == BLAS_MATRIX) {
			ctx()->get_real_mdata()->elementwise_func(add_func_r);
		} else if (type == BLAS_COMPLEX_MATRIX) {
			ctx()->get_complex_mdata()->elementwise_func(add_func_rc);
		} else {
			return ERR_INVALID_TYPE;
		}

		return 0;
	} else if (_is_complex_number(other) && get_type() == BLAS_COMPLEX_MATRIX) {
		complex_t val = _variant_to_complex(other);

		auto add_func_cc = [&](const complex_t &a) { return a + val; };

		ctx()->get_complex_mdata()->elementwise_func(add_func_cc);

		return 0;
	}

	GDBlasMat *o = _cast(other);

	if (o == nullptr || this == o) {
		GDBLAS_ERROR("o == nullptr || this == o");

		return ERR_INVALID_INPUT;
	}

	if (!_can_add(o)) {
		return ERR_INCOMPATIBLE_DIM;
	}

	return _add_implementation(o);
}

Variant GDBlasMat::sub(Variant other) {
	if (_is_zero_dim()) {
		GDBLAS_ERROR("Matrix is dimensionless");

		return ERR_INVALID_DIM;
	}

	if (_is_real_number(other)) {
		scalar_t val = other;

		auto add_func_r = [&](scalar_t a) { return a - val; };
		auto add_func_rc = [&](const complex_t &a) { return a - val; };

		int type = get_type();
		if (type == BLAS_MATRIX) {
			ctx()->get_real_mdata()->elementwise_func(add_func_r);
		} else if (type == BLAS_COMPLEX_MATRIX) {
			ctx()->get_complex_mdata()->elementwise_func(add_func_rc);
		} else {
			return ERR_INVALID_TYPE;
		}

		return 0;
	} else if (_is_complex_number(other) && get_type() == BLAS_COMPLEX_MATRIX) {
		complex_t val = _variant_to_complex(other);

		auto add_func_cc = [&](const complex_t &a) { return a - val; };

		ctx()->get_complex_mdata()->elementwise_func(add_func_cc);

		return 0;
	}

	GDBlasMat *o = _cast(other);

	if (o == nullptr || this == o) {
		GDBLAS_ERROR("o == nullptr || this == o");

		return ERR_INVALID_INPUT;
	}

	if (!_can_add(o)) {
		return ERR_INCOMPATIBLE_DIM;
	}

	return _sub_implementation(o);
}

void GDBlasMat::transpose() {
	if (_is_zero_dim()) {
		GDBLAS_ERROR("Matrix is dimensionless");
	}

	_transpose_implementation();
}

void GDBlasMat::hermitian() {
	if (_is_zero_dim()) {
		GDBLAS_ERROR("Matrix is dimensionless");
	}

	_hermitian_implementation();
}

Variant GDBlasMat::is_eq(Variant other, scalar_t p_eps, int p_norm_type) {
	if (_is_zero_dim()) {
		GDBLAS_ERROR("Matrix is dimensionless");

		return false;
	}

	GDBlasMat *o = _cast(other);

	if (o == nullptr || this == o) {
		GDBLAS_ERROR("o == nullptr || this == o");

		return false;
	}

	if (!_cmp_size(o)) {
		return false;
	}

	return _is_eq_implementation(o, p_eps, p_norm_type);
}

Variant GDBlasMat::mul(Variant other) {
	if (_is_zero_dim()) {
		GDBLAS_ERROR("Matrix is dimensionless");

		return ERR_INVALID_DIM;
	}

	Variant::Type t = other.get_type();

	if (t == Variant::INT || t == Variant::FLOAT) {
		scalar_t s = other;

		return _muls_implementation(s);
	}
	else if (t == Variant::VECTOR2) {
		complex_t c = _variant_to_complex(other);

		return _muls_implementation(c);
	}
	else {
		GDBlasMat *o = _cast(other);

		if (o == nullptr || this == o) {
			GDBLAS_ERROR("o == nullptr || this == o");

			return ERR_INVALID_INPUT;
		}

		if (o->_is_scalar()) {
			if (o->get_type() == BLAS_MATRIX)
				return _muls_implementation(o->ctx()->get(0, 0));
			else if (o->get_type() == BLAS_COMPLEX_MATRIX)
				return _muls_implementation(o->ctx()->getc(0, 0));
		}
		else if (!_cmp_size(o)) {
			GDBLAS_ERROR("Matrix dimensions do not match");

			return ERR_INCOMPATIBLE_DIM;
		}

		auto func_rrr = [&](const scalar_t &a, const scalar_t &b) { return a * b; };
		auto func_ccr = [&](const complex_t &a, const scalar_t &b) { return a * b; };
		auto func_ccc = [&](const complex_t &a, const complex_t &b) { return a * b; };

		int type1 = get_type();
		int type2 = o->get_type();
		if (type1 == BLAS_MATRIX && type2 == BLAS_MATRIX) {
			ctx()->get_real_mdata()->elementwise_func(func_rrr, o->ctx()->get_real_mdata()->matrix());
		}
		else if (type1 == BLAS_COMPLEX_MATRIX && type2 == BLAS_MATRIX) {
			ctx()->get_complex_mdata()->elementwise_func(func_ccr, o->ctx()->get_real_mdata()->matrix());
		}
		else if (type1 == BLAS_COMPLEX_MATRIX && type2 == BLAS_COMPLEX_MATRIX) {
			ctx()->get_complex_mdata()->elementwise_func(func_ccc, o->ctx()->get_complex_mdata()->matrix());
		}
		else {
			return ERR_INVALID_TYPE;
		}

		return 0;
	}

	return ERR_INVALID_INPUT;
}

Variant GDBlasMat::div(Variant other) {
	if (_is_zero_dim()) {
		GDBLAS_ERROR("Matrix is dimensionless");

		return ERR_INVALID_DIM;
	}

	Variant::Type t = other.get_type();

	if (t == Variant::INT || t == Variant::FLOAT) {
		scalar_t s = other;

		return _muls_implementation(1.0 / s);
	}
	else if (t == Variant::VECTOR2) {
		complex_t c = _variant_to_complex(other);

		c = std::conj(c) / std::norm(c);

		return _muls_implementation(c);
	}
	else {
		GDBlasMat *o = _cast(other);

		if (o == nullptr || this == o) {
			GDBLAS_ERROR("o == nullptr || this == o");

			return ERR_INVALID_INPUT;
		}

		if (o->_is_scalar()) {
			if (o->get_type() == BLAS_MATRIX) {
				return _muls_implementation(1.0 / o->ctx()->get(0, 0));
			}
			else if (o->get_type() == BLAS_COMPLEX_MATRIX) {
				complex_t c(1.0 / o->ctx()->getc(0, 0));

				return _muls_implementation(c);
			}
		}
		else if (!_cmp_size(o)) {
			GDBLAS_ERROR("Matrix dimensions do not match");

			return ERR_INCOMPATIBLE_DIM;
		}

		auto func_rrr = [&](const scalar_t &a, const scalar_t &b) { return a / b; };
		auto func_ccr = [&](const complex_t &a, const scalar_t &b) { return a / b; };
		auto func_ccc = [&](const complex_t &a, const complex_t &b) { return a / b; };

		int type1 = get_type();
		int type2 = o->get_type();
		if (type1 == BLAS_MATRIX && type2 == BLAS_MATRIX) {
			ctx()->get_real_mdata()->elementwise_func(func_rrr, o->ctx()->get_real_mdata()->matrix());
		}
		else if (type1 == BLAS_COMPLEX_MATRIX && type2 == BLAS_MATRIX) {
			ctx()->get_complex_mdata()->elementwise_func(func_ccr, o->ctx()->get_real_mdata()->matrix());
		}
		else if (type1 == BLAS_COMPLEX_MATRIX && type2 == BLAS_COMPLEX_MATRIX) {
			ctx()->get_complex_mdata()->elementwise_func(func_ccc, o->ctx()->get_complex_mdata()->matrix());
		}
		else {
			return ERR_INVALID_TYPE;
		}

		return 0;
	}

	return ERR_INVALID_INPUT;
}

Variant GDBlasMat::fill(Variant other) {
	if (_is_zero_dim()) {
		GDBLAS_ERROR("Matrix is dimensionless");

		return ERR_INVALID_DIM;
	}

	Variant::Type t = other.get_type();

	if (t == Variant::INT || t == Variant::FLOAT) {
		scalar_t s = other;

		return _fill_implementation(s);
	}
	else if (t == Variant::VECTOR2) {
		complex_t c = _variant_to_complex(other);

		return _fill_implementation(c);
	}

	return ERR_INVALID_INPUT;
}

Variant GDBlasMat::eye(Variant p_val) {
	if (_is_zero_dim()) {
		GDBLAS_ERROR("Matrix is dimensionless");

		return ERR_INVALID_DIM;
	}

	Variant::Type t = p_val.get_type();

	if (t == Variant::INT || t == Variant::FLOAT) {
		scalar_t s = p_val;

		int result = _fill_implementation(0.0);
		if (result)
			return result;

		return _fill_implementation(s, true);
	}
	else if (t == Variant::VECTOR2) {
		complex_t c = _variant_to_complex(p_val);

		int result = _fill_implementation(0.0);
		if (result)
			return result;

		return _fill_implementation(c, true);
	}

	return ERR_INVALID_INPUT;
}

Variant GDBlasMat::reset() {
	return _fill_implementation(0.0);
}

void GDBlasMat::conj() {
	if (_is_zero_dim()) {
		GDBLAS_ERROR("Matrix is dimensionless");
	}

	_conj_implementation();
}

Variant GDBlasMat::real(Variant p_matrix) {
	int type = p_matrix.get_type();

	if (type == Variant::NIL) {
		int error = 0;

		if (_is_zero_dim()) {
			GDBLAS_ERROR("Matrix is dimensionless");

			return ERR_INVALID_DIM;
		}

		Ref< GDBlasMat > mat = _real_implementation(&error);
		if (error)
			return Variant();

		return Variant(mat);
	}

	GDBlasMat *o = _cast(p_matrix);

	if (o == nullptr || this == o) {
		GDBLAS_ERROR("o == nullptr || this == o");

		return ERR_INVALID_INPUT;
	}

	if (!_cmp_size(o)) {
		return ERR_INCOMPATIBLE_DIM;
	}

	return _set_real_implementation(o);
}

Variant GDBlasMat::imag(Variant p_matrix) {
	if (get_type() != BLAS_COMPLEX_MATRIX)
		return Variant();

	int type = p_matrix.get_type();

	if (type == Variant::NIL) {
		int error = 0;

		if (_is_zero_dim()) {
			GDBLAS_ERROR("Matrix is dimensionless");

			return ERR_INVALID_DIM;
		}

		Ref< GDBlasMat > mat = _imag_implementation(&error);
		if (error)
			return Variant();

		return Variant(mat);
	}

	GDBlasMat *o = _cast(p_matrix);

	if (o == nullptr || this == o) {
		GDBLAS_ERROR("o == nullptr || this == o");

		return ERR_INVALID_INPUT;
	}

	if (!_cmp_size(o)) {
		return ERR_INCOMPATIBLE_DIM;
	}

	return _set_imag_implementation(o);
}

Variant GDBlasMat::prod(Variant other) {
	if (_is_zero_dim()) {
		GDBLAS_ERROR("Matrix is dimensionless");

		return ERR_INVALID_DIM;
	}

	GDBlasMat *o = _cast(other);

	if (o == nullptr || this == o) {
		GDBLAS_ERROR("o == nullptr || this == o");

		return ERR_INVALID_INPUT;
	}

	Dimension d1 = _size();
	Dimension d2 = o->_size();

	if (d1.n != d2.m) {
		GDBLAS_ERROR("Can not apply dot product on matrices with sizes "
					  "(%dx%d) and (%dx%d)", d1.m, d1.n, d2.m, d2.n);

		return ERR_INCOMPATIBLE_DIM;
	}

	int error = 0;
	Ref< GDBlasMat > mat = _prod_implementation(o, &error);
	if (error)
		return Variant();

	return Variant(mat);
}

Variant GDBlasMat::inv() {
	if (_is_zero_dim()) {
		GDBLAS_ERROR("Matrix is dimensionless");

		return Variant();
	}

	if (!_is_square()) {
		GDBLAS_ERROR("Matrix must be square");

		return Variant();
	}

	int error = 0;
	Ref< GDBlasMat > mat = _inv_implementation(&error);
	if (error)
		return Variant();

	return Variant(mat);
}

Variant GDBlasMat::to_array() {
	Array arr;

	if (_is_zero_dim()) {
		GDBLAS_ERROR("Matrix is dimensionless");

		return Variant(arr);
	}

	Dimension d = _size();
	int type = get_type();

	arr.resize(d.m);

	for (s_t i = 0; i < d.m; ++i) {
		Array row;
		row.resize(d.n);

		for (s_t j = 0; j < d.n; ++j) {
			if (type == BLAS_MATRIX) {
				row[j] = ctx()->get(i, j);
			} else if (type == BLAS_COMPLEX_MATRIX) {
				row[j] = _complex_to_variant(ctx()->getc(i, j));
			}
		}

		arr[i] = row;
	}

	return Variant(arr);
}

Variant GDBlasMat::from_array(Array p_array) {
	s_t m = GDBLAS_SIZE_CAST(p_array.size());
	s_t n = 0;
	bool resized = false;
	int type = get_type();

	if (!m) {
		return ERR_INVALID_DIM;
	}

	for (s_t i = 0; i < m; ++i) {
		Variant &v = p_array[i];

		if (v.get_type() != Variant::ARRAY) {
			return ERR_INVALID_TYPE;
		}

		Array row = v;

		if (!n) {
			n = GDBLAS_SIZE_CAST(row.size());
		} else if (n != GDBLAS_SIZE_CAST(row.size())) {
			return ERR_INVALID_DIM;
		}

		if (n && !resized) {
			resize((int)m, (int)n);
			resized = true;
		}

		for (s_t j = 0; j < n; ++j) {
			Variant &entry = row[j];
			int error = 0;

			if (type == BLAS_MATRIX && _is_real_number(entry)) {
				scalar_t val = entry;
				error = ctx()->set(val, i, j);
			} else if (type == BLAS_COMPLEX_MATRIX) {
				error = ctx()->set(_variant_to_complex(entry), i, j);
			}

			if (error)
				return error;
		}
	}

	return 0;
}

#define MATH_FUNC_GEN(func) do { \
	int type = get_type(); \
	if (type == BLAS_MATRIX) { \
		ctx()->get_real_mdata()->elementwise_func(static_cast<scalar_t(*)(scalar_t)>(func)); \
	} else if (type == BLAS_COMPLEX_MATRIX) { \
		ctx()->get_complex_mdata()->elementwise_func(static_cast<complex_t(*)(const complex_t&)>(func)); \
	} else { \
		return ERR_INVALID_TYPE; \
	} \
} while (0)

#define MATH_FUNC_GEN_REAL(func) do { \
	int type = get_type(); \
	if (type == BLAS_MATRIX) { \
		ctx()->get_real_mdata()->elementwise_func(static_cast<scalar_t(*)(scalar_t)>(func)); \
	} else { \
		return ERR_INVALID_TYPE; \
	} \
} while (0)

#include "gdblas_math_extra.inc"

Variant GDBlasMat::sin() {
	MATH_FUNC_GEN(std::sin);

	return 0;
}

Variant GDBlasMat::cos() {
	MATH_FUNC_GEN(std::cos);

	return 0;
}

Variant GDBlasMat::abs() {
	MATH_FUNC_GEN(std::abs2);

	return 0;
}

Variant GDBlasMat::exp() {
	MATH_FUNC_GEN(std::exp);

	return 0;
}

Variant GDBlasMat::log() {
	MATH_FUNC_GEN(std::log);

	return 0;
}

Variant GDBlasMat::log10() {
	MATH_FUNC_GEN(std::log10);

	return 0;
}

Variant GDBlasMat::log2() {
	MATH_FUNC_GEN(std::log2);

	return 0;
}

Variant GDBlasMat::sqrt() {
	MATH_FUNC_GEN(std::sqrt);

	return 0;
}

Variant GDBlasMat::cbrt() {
	MATH_FUNC_GEN(std::cbrt);

	return 0;
}

Variant GDBlasMat::tan() {
	MATH_FUNC_GEN(std::tan);

	return 0;
}

Variant GDBlasMat::asin() {
	MATH_FUNC_GEN(std::asin);

	return 0;
}

Variant GDBlasMat::acos() {
	MATH_FUNC_GEN(std::acos);

	return 0;
}

Variant GDBlasMat::atan() {
	MATH_FUNC_GEN(std::atan);

	return 0;
}

Variant GDBlasMat::sinh() {
	MATH_FUNC_GEN(std::sinh);

	return 0;
}

Variant GDBlasMat::cosh() {
	MATH_FUNC_GEN(std::cosh);

	return 0;
}

Variant GDBlasMat::tanh() {
	MATH_FUNC_GEN(std::tanh);

	return 0;
}

Variant GDBlasMat::atanh() {
	MATH_FUNC_GEN(std::atanh);

	return 0;
}

Variant GDBlasMat::erf() {
	MATH_FUNC_GEN_REAL(std::erf);

	return 0;
}

Variant GDBlasMat::erfc() {
	MATH_FUNC_GEN_REAL(std::erfc);

	return 0;
}

Variant GDBlasMat::tgamma() {
	MATH_FUNC_GEN_REAL(std::tgamma);

	return 0;
}

Variant GDBlasMat::lgamma() {
	MATH_FUNC_GEN_REAL(std::lgamma);

	return 0;
}

Variant GDBlasMat::ceil() {
	MATH_FUNC_GEN_REAL(std::ceil);

	return 0;
}

Variant GDBlasMat::floor() {
	MATH_FUNC_GEN_REAL(std::floor);

	return 0;
}

Variant GDBlasMat::trunc() {
	MATH_FUNC_GEN_REAL(std::trunc);

	return 0;
}

Variant GDBlasMat::round() {
	MATH_FUNC_GEN_REAL(std::round);

	return 0;
}

Variant GDBlasMat::pow(Variant p_exponent) {
	scalar_t exp_r;
	complex_t exp_c;
	bool is_complex_exp = false;
	if (_is_real_number(p_exponent)) {
		exp_r = p_exponent;
	} else if (_is_complex_number(p_exponent)) {
		exp_c = _variant_to_complex(p_exponent);
		is_complex_exp = true;
	} else {
		return ERR_INVALID_INPUT;
	}

	auto pow_func_r = [&](scalar_t a) { return std::pow(a, exp_r); };
	auto pow_func_rc = [&](const complex_t &a) { return std::pow(a, exp_r); };
	auto pow_func_cc = [&](const complex_t &a) { return std::pow(a, exp_c); };

	int type = get_type();
	if (type == BLAS_MATRIX) {
		ctx()->get_real_mdata()->elementwise_func(pow_func_r);
	} else if (type == BLAS_COMPLEX_MATRIX) {
		if (is_complex_exp)
			ctx()->get_complex_mdata()->elementwise_func(pow_func_cc);
		else
			ctx()->get_complex_mdata()->elementwise_func(pow_func_rc);
	} else {
		return ERR_INVALID_TYPE;
	}

	return 0;
}

Variant GDBlasMat::integrate(int axis) {
	if (_is_zero_dim()) {
		GDBLAS_ERROR("Matrix is dimensionless");

		return Variant();
	}

	int error = 0;
	Ref< GDBlasMat > mat = _integrate_implementation(axis, &error);
	if (error)
		return Variant();

	return Variant(mat);
}

Variant GDBlasMat::norm(int p_norm_type) {
	switch (p_norm_type) {
	case NORM_1:
	case NORM_INF:
	case NORM_FRO:
		break;
	default:
		GDBLAS_ERROR("Undefined norm type");

		return Variant();
	}

	if (_is_zero_dim()) {
		GDBLAS_ERROR("Matrix is dimensionless");

		return Variant();
	}

	int error = 0;
	scalar_t l = _norm_implementation(p_norm_type, &error);
	if (error)
		return Variant();

	return Variant(l);
}

Variant GDBlasMat::mean(int axis) {
	if (_is_zero_dim()) {
		GDBLAS_ERROR("Matrix is dimensionless");

		return Variant();
	}

	int error = 0;
	Ref< GDBlasMat > mat = _mean_implementation(axis, &error);
	if (error)
		return Variant();

	return Variant(mat);
}

Variant GDBlasMat::min(int axis) {
	if (_is_zero_dim()) {
		GDBLAS_ERROR("Matrix is dimensionless");

		return Variant();
	}

	int error = 0;
	Ref< GDBlasMat > mat = _min_max_implementation(axis, &error, true);
	if (error)
		return Variant();

	return Variant(mat);
}

Variant GDBlasMat::max(int axis) {
	if (_is_zero_dim()) {
		GDBLAS_ERROR("Matrix is dimensionless");

		return Variant();
	}

	int error = 0;
	Ref< GDBlasMat > mat = _min_max_implementation(axis, &error, false);
	if (error)
		return Variant();

	return Variant(mat);
}

Variant GDBlasMat::argmin(int axis) {
	if (_is_zero_dim()) {
		GDBLAS_ERROR("Matrix is dimensionless");

		return Variant();
	}

	int error = 0;
	Array arg = _arg_min_max_implementation(axis, &error, true);

	return arg;
}

Variant GDBlasMat::argmax(int axis) {
	if (_is_zero_dim()) {
		GDBLAS_ERROR("Matrix is dimensionless");

		return Variant();
	}

	int error = 0;
	Array arg = _arg_min_max_implementation(axis, &error, false);

	return arg;
}

Variant GDBlasMat::unary_func(Callable p_func, Variant p_args, bool p_indexed) {
	int type = get_type();
	int args_type = p_args.get_type();

	if (args_type != Variant::ARRAY && args_type != Variant::NIL) {
		return ERR_INVALID_INPUT;
	}

	Array args = p_args;

	auto func_rr = [&](scalar_t &a) -> scalar_t {
		if (args_type == Variant::NIL)
			return p_func.call(a);
		else if (args_type == Variant::ARRAY)
			return p_func.call(a, args);

		return 0;
	};

	auto func_rr_ij = [&](scalar_t &a, int i, int j) -> scalar_t {
		if (args_type == Variant::NIL)
			return p_func.call(a, i, j);
		else if (args_type == Variant::ARRAY)
			return p_func.call(a, args, i, j);

		return 0;
	};

	auto func_cc = [&](complex_t &a) -> complex_t {
		if (args_type == Variant::NIL)
			return _variant_to_complex(p_func.call(_complex_to_variant(a)));
		else if (args_type == Variant::ARRAY)
			return _variant_to_complex(p_func.call(_complex_to_variant(a), args));

		return 0.0;
	};

	auto func_cc_ij = [&](complex_t &a, int i, int j) -> complex_t {
		if (args_type == Variant::NIL)
			return _variant_to_complex(p_func.call(_complex_to_variant(a), i, j));
		else if (args_type == Variant::ARRAY)
			return _variant_to_complex(p_func.call(_complex_to_variant(a), args, i, j));

		return 0.0;
	};

	Variant tmp;
	if (args_type == Variant::ARRAY) {
		if (p_indexed)
			tmp = p_func.call(1.0, args, 0, 0);
		else
			tmp = p_func.call(1.0, args);
	} else if (args_type == Variant::NIL) {
		if (p_indexed)
			tmp = p_func.call(1.0, 0, 0);
		else
			tmp = p_func.call(1.0);
	}

	if (!((_is_real_number(tmp) && type == BLAS_MATRIX)
		  || (_is_complex_number(tmp) && type == BLAS_COMPLEX_MATRIX))) {
		GDBLAS_ERROR("Callable must return number for real matrix and Vector2 for complex matrix");
		return ERR_INVALID_TYPE;
	}

	if (type == BLAS_MATRIX) {
		if (p_indexed)
			ctx()->get_real_mdata()->elementwise_func_indexed(func_rr_ij);
		else
			ctx()->get_real_mdata()->elementwise_func(func_rr);
	} else if (type == BLAS_COMPLEX_MATRIX) {
		if (p_indexed)
			ctx()->get_complex_mdata()->elementwise_func_indexed(func_cc_ij);
		else
			ctx()->get_complex_mdata()->elementwise_func(func_cc);
	} else {
		return ERR_INVALID_TYPE;
	}

	return 0;
}

Variant GDBlasMat::conv(Variant p_other, bool p_same) {
	if (_is_zero_dim()) {
		GDBLAS_ERROR("Matrix is dimensionless");

		return Variant();
	}

	GDBlasMat *o = _cast(p_other);

	if (o == nullptr || this == o) {
		GDBLAS_ERROR("o == nullptr || this == o");

		return Variant();
	}

	if (o->_is_zero_dim()) {
		GDBLAS_ERROR("Matrix is dimensionless");

		return Variant();
	}

	int error = 0;
	Ref<GDBlasMat> mat = _conv_implementation(o, p_same, &error);
	if (error)
		return Variant();

	return Variant(mat);
}

Variant GDBlasMat::pack(int p_component) {
	if (_is_zero_dim()) {
		GDBLAS_ERROR("Matrix is dimensionless");

		return Variant();
	}

	PackedFloat64Array packed_data;

	int error = _pack_implementation(packed_data, p_component, 1);
	if (error) {
		return PackedFloat64Array();
	}

	return packed_data;
}

Variant GDBlasMat::unpack(Variant p_packed_data, int p_component, int p_step, int p_offset) {
	if (_is_zero_dim()) {
		GDBLAS_ERROR("Matrix is dimensionless");

		return ERR_INVALID_DIM;
	}

	int64_t packed_size = 0;

	Variant::Type packed_type = p_packed_data.get_type();
	if (packed_type == Variant::PACKED_BYTE_ARRAY) {
		packed_size = ((PackedByteArray)p_packed_data).size();
	} else if (packed_type == Variant::PACKED_FLOAT32_ARRAY) {
		packed_size = ((PackedFloat32Array)p_packed_data).size();
	} else if (packed_type == Variant::PACKED_FLOAT64_ARRAY) {
		packed_size = ((PackedByteArray)p_packed_data).size();
	} else {
		GDBLAS_ERROR("Packed data type must be BYTE, FLOAT32 or FLOAT64 array");

		return ERR_INVALID_INPUT;
	}

	if (p_step < 1 || packed_size % p_step != 0) {
		GDBLAS_ERROR("Step size must be a positive integer divisor of input array size");

		return ERR_INVALID_INPUT;
	}

	packed_size /= p_step;

	int type = get_type();
	Dimension d = _size();
	if (type == BLAS_MATRIX && d.m * d.n != packed_size) {
		GDBLAS_ERROR("Can not unpack %ld values into a matrix of size (%ux%u)", packed_size, d.m, d.n);

		return ERR_INVALID_DIM;
	} else if (type == BLAS_COMPLEX_MATRIX && p_component == BOTH_COMPONENTS
			   && 2 * d.m * d.n != packed_size) {
		GDBLAS_ERROR("Can not unpack %ld values into a complex matrix of size (%ux%u)",
					 packed_size, d.m, d.n);

		return ERR_INVALID_DIM;
	}

	if (packed_type == Variant::PACKED_FLOAT64_ARRAY) {
		PackedFloat64Array packed_data = p_packed_data;
		return _unpack_implementation(packed_data, p_component, p_step);
	} else if (packed_type == Variant::PACKED_FLOAT32_ARRAY) {
		if (type == BLAS_COMPLEX_MATRIX && p_component == BOTH_COMPONENTS) {
			GDBLAS_ERROR("Can not unpack FLOAT32 array into complex matrix");

			return ERR_INVALID_INPUT;
		}
		PackedFloat32Array packed_data = p_packed_data;
		return _unpack_implementation<float>(packed_data, p_component, p_step);
	} else if (packed_type == Variant::PACKED_BYTE_ARRAY) {
		if (type == BLAS_COMPLEX_MATRIX && p_component == BOTH_COMPONENTS) {
			GDBLAS_ERROR("Can not unpack BYTE array into complex matrix");

			return ERR_INVALID_INPUT;
		}
		PackedByteArray packed_data = p_packed_data;
		return _unpack_implementation<uint8_t>(packed_data, p_component, p_step, p_offset);
	}

	return ERR_INVALID_INPUT;
}

Variant GDBlasMat::downsample(int p_factor_m, int p_factor_n, Variant p_filter) {
	if (p_factor_m < 1 || p_factor_n < 1) {
		GDBLAS_ERROR("Downsampling factor must be a positive integer");

		return Variant();
	}

	GDBlasMat *filter = _cast(p_filter);

	Ref<GDBlasMat> output = _downsample_implementation(p_factor_m, p_factor_n, filter);

	return Variant(output);
}

#if defined(GDBLAS_WITH_ODE)
Variant GDBlasMat::eval_ode(Callable p_f, double p_dt, double p_max_step) {
	typedef std::vector<scalar_t> state_type;

	Dimension d = _size();

	if (!d.m || d.n != 1) {
		GDBLAS_ERROR("ODE can be evaluated on only column vectors");

		return ERR_INVALID_DIM;
	}

	if (get_type() != BLAS_MATRIX) {
		GDBLAS_ERROR("ODE can be evaluated on only real vectors");

		return ERR_INVALID_TYPE;
	}

	int64_t result = 0;

#define VEC_TMP(m) m->ctx()->get_real_mdata()->matrix().col(0)

	auto fx_wrapper = [&](const state_type &x, state_type &dx, double t) -> void {
		Variant tmp = p_f.call(Variant(this), t);

		GDBlasMat *dx_mat = _cast(tmp);
		if (dx_mat == nullptr) {
			GDBLAS_ERROR("ODE Callable must return a matrix of same dimension as its first argument");

			result = ERR_INVALID_INPUT;

			return;
		}

		if (dx_mat == this) {
			GDBLAS_ERROR("dx matrix can not be x itself");

			result = ERR_INVALID_INPUT;

			return;
		}

		dx.assign(VEC_TMP(dx_mat).begin(), VEC_TMP(dx_mat).end());
	};

	state_type x_tmp(VEC_TMP(this).begin(), VEC_TMP(this).end());
	double end_time = m_prev_end_time + p_dt;
	size_t step_count = boost::numeric::odeint::integrate(fx_wrapper, x_tmp,
														  m_prev_end_time, end_time, p_max_step);

	if (step_count) {
		auto it1 = VEC_TMP(this).begin();
		auto it2 = x_tmp.begin();
		for (; it1 != VEC_TMP(this).end() && it2 != x_tmp.end(); ++it1, ++it2) {
			*it1 = *it2;
		}
	}

#undef VEC_TMP

	if (!result)
		result = step_count;

	return result;
}
#endif

/********************* Implementations ********************/

void GDBlasMat::_resize_implementation(s_t m, s_t n) {
	if (get_type() == BLAS_MATRIX) {
		ctx()->resize<Context::RealMatrixData>(m, n);
		ctx()->clear<Context::RealMatrixData>();
	} else if (get_type() == BLAS_COMPLEX_MATRIX) {
		ctx()->resize<Context::ComplexMatrixData>(m, n);
		ctx()->clear<Context::ComplexMatrixData>();
	}
}

int GDBlasMat::_add_implementation(GDBlasMat *other) {
	if (get_type() == BLAS_MATRIX && other->get_type() == BLAS_MATRIX) {
		ctx()->add<Context::RealMatrixData, Context::RealMatrixData>(other->ctx());
	} else if (get_type() == BLAS_COMPLEX_MATRIX && other->get_type() == BLAS_MATRIX) {
		ctx()->add<Context::ComplexMatrixData, Context::RealMatrixData>(other->ctx());
	} else if (get_type() == BLAS_COMPLEX_MATRIX && other->get_type() == BLAS_COMPLEX_MATRIX) {
		ctx()->add<Context::ComplexMatrixData, Context::ComplexMatrixData>(other->ctx());
	} else {
		return ERR_INVALID_TYPE;
	}

	return 0;
}

int GDBlasMat::_sub_implementation(GDBlasMat *other) {
	if (get_type() == BLAS_MATRIX && other->get_type() == BLAS_MATRIX) {
		ctx()->sub<Context::RealMatrixData, Context::RealMatrixData>(other->ctx());
	} else if (get_type() == BLAS_COMPLEX_MATRIX && other->get_type() == BLAS_MATRIX) {
		ctx()->sub<Context::ComplexMatrixData, Context::RealMatrixData>(other->ctx());
	} else if (get_type() == BLAS_COMPLEX_MATRIX && other->get_type() == BLAS_COMPLEX_MATRIX) {
		ctx()->sub<Context::ComplexMatrixData, Context::ComplexMatrixData>(other->ctx());
	} else {
		return ERR_INVALID_TYPE;
	}

	return 0;
}

void GDBlasMat::_transpose_implementation() {
	if (get_type() == BLAS_MATRIX) {
		ctx()->transpose<Context::RealMatrixData>();
	} else if (get_type() == BLAS_COMPLEX_MATRIX) {
		ctx()->transpose<Context::ComplexMatrixData>();
	}
}

void GDBlasMat::_hermitian_implementation() {
	if (get_type() == BLAS_MATRIX) {
		ctx()->transpose<Context::RealMatrixData>();
	} else if (get_type() == BLAS_COMPLEX_MATRIX) {
		_conj_implementation();
		ctx()->transpose<Context::ComplexMatrixData>();
	}
}

int GDBlasMat::_eq_implementation(GDBlasMat *other) {
	if (get_type() == BLAS_MATRIX && other->get_type() == BLAS_MATRIX) {
		ctx()->eq<Context::RealMatrixData, Context::RealMatrixData>(other->ctx());
	} else if (get_type() == BLAS_COMPLEX_MATRIX && other->get_type() == BLAS_COMPLEX_MATRIX) {
		ctx()->eq<Context::ComplexMatrixData, Context::ComplexMatrixData>(other->ctx());
	} else if (get_type() == BLAS_COMPLEX_MATRIX && other->get_type() == BLAS_MATRIX) {
		ctx()->eq<Context::ComplexMatrixData, Context::RealMatrixData>(other->ctx());
	} else {
		return ERR_INVALID_TYPE;
	}

	return 0;
}

bool GDBlasMat::_is_eq_implementation(GDBlasMat *other, scalar_t eps, int norm_type) {
	bool result = false;

	if (get_type() == BLAS_MATRIX && other->get_type() == BLAS_MATRIX) {
		result = ctx()->is_eq<Context::RealMatrixData, Context::RealMatrixData>(other->ctx(),
				eps, norm_type);
	} else if (get_type() == BLAS_COMPLEX_MATRIX && other->get_type() == BLAS_COMPLEX_MATRIX) {
		result = ctx()->is_eq<Context::ComplexMatrixData, Context::ComplexMatrixData>(other->ctx(),
				eps, norm_type);
	} else if (get_type() == BLAS_COMPLEX_MATRIX && other->get_type() == BLAS_MATRIX) {
		result = ctx()->is_eq<Context::ComplexMatrixData, Context::RealMatrixData>(other->ctx(),
				eps, norm_type);
	} else if (get_type() == BLAS_MATRIX && other->get_type() == BLAS_COMPLEX_MATRIX) {
		result = other->ctx()->is_eq<Context::ComplexMatrixData, Context::RealMatrixData>(ctx(),
				eps, norm_type);
	}

	return result;
}

int GDBlasMat::_muls_implementation(scalar_t s) {
	if (get_type() == BLAS_MATRIX) {
		ctx()->muls<Context::RealMatrixData, scalar_t>(s);
	} else if (get_type() == BLAS_COMPLEX_MATRIX) {
		ctx()->muls<Context::ComplexMatrixData, scalar_t>(s);
	} else {
		return ERR_INVALID_TYPE;
	}

	return 0;
}

int GDBlasMat::_muls_implementation(complex_t s) {
	if (get_type() == BLAS_COMPLEX_MATRIX) {
		ctx()->muls<Context::ComplexMatrixData, complex_t>(s);
	} else {
		return ERR_INVALID_TYPE;
	}

	return 0;
}

Ref<GDBlasMat> GDBlasMat::_copy_implementation(int *error) {
	Dimension d = _size();

	Ref<GDBlasMat> mat = GDBlasMat::new_mat(d.m, d.n, get_type(), error);

	mat->_eq_implementation(this);

	return mat;
}

int GDBlasMat::_fill_implementation(scalar_t s, bool diag) {
	if (get_type() == BLAS_MATRIX) {
		if (s == 0.0 && !diag)
			ctx()->clear<Context::RealMatrixData>();
		else
			ctx()->fill<Context::RealMatrixData, scalar_t>(s, diag);
	} else if (get_type() == BLAS_COMPLEX_MATRIX) {
		if (s == 0.0 && !diag) {
			ctx()->clear<Context::ComplexMatrixData>();
		} else {
			complex_t c;
			c.real(s);

			ctx()->fill<Context::ComplexMatrixData, complex_t>(c, diag);
		}
	} else {
		return ERR_INVALID_TYPE;
	}

	return 0;
}

int GDBlasMat::_fill_implementation(complex_t &s, bool diag) {
	if (get_type() == BLAS_COMPLEX_MATRIX) {
		ctx()->fill<Context::ComplexMatrixData, complex_t>(s, diag);
	} else {
		return ERR_INVALID_TYPE;
	}

	return 0;
}

void GDBlasMat::_conj_implementation() {
	if (get_type() == BLAS_COMPLEX_MATRIX) {
		auto lconj = [](complex_t &a) { return std::conj(a); };
		ctx()->get_complex_mdata()->elementwise_func(lconj);
	}
}

Ref<GDBlasMat> GDBlasMat::_real_implementation(int *error) {
	Dimension d = _size();

	Ref<GDBlasMat> mat = new_mat(d.m, d.n, BLAS_MATRIX, error);
	if (*error)
		return mat;

	int type = get_type();
	if (type == BLAS_COMPLEX_MATRIX)
		ctx()->set_real_or_imag(mat->ctx(), true, true);
	else if (type == BLAS_MATRIX)
		mat->ctx()->eq<Context::RealMatrixData, Context::RealMatrixData>(ctx());

	return mat;
}

int GDBlasMat::_set_real_implementation(GDBlasMat *other) {
	int type1 = get_type();
	int type2 = other->get_type();

	if (type1 == BLAS_COMPLEX_MATRIX && type2 == BLAS_MATRIX) {
		ctx()->set_real_or_imag(other->ctx(), true, false);
	} else if (type1 == BLAS_MATRIX && type2 == BLAS_MATRIX) {
		ctx()->eq<Context::RealMatrixData, Context::RealMatrixData>(other->ctx());
	} else {
		return ERR_INVALID_TYPE;
	}

	return 0;
}

Ref<GDBlasMat> GDBlasMat::_imag_implementation(int *error) {
	Dimension d = _size();

	Ref<GDBlasMat> mat = new_mat(d.m, d.n, BLAS_MATRIX, error);
	if (*error)
		return mat;

	*error = 0;

	if (get_type() == BLAS_COMPLEX_MATRIX)
		ctx()->set_real_or_imag(mat->ctx(), false, true);
	else
		*error = ERR_INVALID_TYPE;

	return mat;
}

int GDBlasMat::_set_imag_implementation(GDBlasMat *other) {
	if (get_type() == BLAS_COMPLEX_MATRIX && other->get_type() == BLAS_MATRIX) {
		ctx()->set_real_or_imag(other->ctx(), false, false);
	} else {
		return ERR_INVALID_TYPE;
	}

	return 0;
}

Ref<GDBlasMat> GDBlasMat::_prod_implementation(GDBlasMat *other, int *error) {
	Dimension d1 = _size();
	Dimension d2 = other->_size();

	int type = get_type();
	if (other->get_type() == BLAS_COMPLEX_MATRIX) {
		type = BLAS_COMPLEX_MATRIX;
	}

	Ref<GDBlasMat> mat = new_mat(d1.m, d2.n, type, error);
	if (*error)
		return mat;

	mat->ctx()->prod(this->ctx(), other->ctx());

	return mat;
}

Ref<GDBlasMat> GDBlasMat::_inv_implementation(int *error) {
	Dimension d = _size();

	Ref<GDBlasMat> mat = new_mat(d.m, d.n, get_type(), error);
	if (*error)
		return mat;

	*error = ctx()->inv(mat->ctx());

	return mat;
}

Ref<GDBlasMat> GDBlasMat::_integrate_implementation(int axis, int *error) {
	Ref<GDBlasMat> mat = new_mat(1, 1, get_type(), error);
	if (*error)
		return nullptr;

	*error = ctx()->integrate(mat->ctx(), axis);

	return mat;
}

Ref<GDBlasMat> GDBlasMat::_mean_implementation(int axis, int *error) {
	Ref<GDBlasMat> mat = new_mat(1, 1, get_type(), error);
	if (*error)
		return nullptr;

	*error = ctx()->mean(mat->ctx(), axis);

	return mat;
}

Ref<GDBlasMat> GDBlasMat::_min_max_implementation(int axis, int *error, bool is_min) {
	if (get_type() == BLAS_COMPLEX_MATRIX) {
		*error = ERR_INVALID_TYPE;

		return nullptr;
	}

	Ref<GDBlasMat> mat = new_mat(1, 1, get_type(), error);
	if (*error)
		return nullptr;

	*error = ctx()->min_max(mat->ctx(), is_min, axis);

	return mat;
}

Array GDBlasMat::_arg_min_max_implementation(int axis, int *error, bool is_min) {
	std::vector<index_t> arg{ INVALID_INDEX };

	if (get_type() == BLAS_MATRIX) {
		*error = ctx()->arg_min_max(arg, is_min, axis);
	} else {
		GDBLAS_ERROR("Invalid type");

		*error = ERR_INVALID_TYPE;
	}

	Array result;
	result.resize(arg.size());
	s_t i = 0;
	for (auto iter = arg.begin(); iter != arg.end(); ++iter) {
		result[i++] = Vector2i((*iter)[0], (*iter)[1]);
	}

	return result;
}

scalar_t GDBlasMat::_norm_implementation(int norm_type, int *error) {
	scalar_t l = ctx()->norm(norm_type, error);

	return l;
}

Ref<GDBlasMat> GDBlasMat::_conv_implementation(GDBlasMat *other, bool same, int *error) {
	Dimension d1 = _size();
	Dimension d2 = other->_size();
	Dimension d3;
	*error = 0;

	if (same) {
		d3.m = d1.m;
		d3.n = d1.n;
	} else {
		d3.m = d1.m + d2.m - 1;
		d3.n = d1.n + d2.n - 1;
	}

	int type1 = get_type();
	int type2 = get_type();
	int type3 = type1;
	if (type3 == BLAS_MATRIX)
		type3 = type2;

	Ref<GDBlasMat> mat = new_mat(d3.m, d3.n, type3, error);
	if (*error)
		return nullptr;

	*error = mat->ctx()->conv(ctx(), other->ctx(), same);

	return mat;
}

Ref<GDBlasMat> GDBlasMat::_downsample_implementation(int factor_m, int factor_n, GDBlasMat *filter) {
	int error = 0;
	Dimension d = _size();
	s_t m = d.m / (s_t)factor_m;
	s_t n = d.n / (s_t)factor_n;
	Ref<GDBlasMat> input = this;

	if (filter == nullptr && factor_m == 1 && factor_n == 1) {
		Ref<GDBlasMat> output = _copy_implementation(&error);
		if (error)
			return nullptr;
		else
			return output;
	}

	Ref<GDBlasMat> mat = new_mat(m, n, get_type(), &error);
	if (error)
		return nullptr;

	if (filter != nullptr) {
		int type1 = get_type();
		int type2 = filter->get_type();
		if (type2 == BLAS_COMPLEX_MATRIX && type1 != BLAS_COMPLEX_MATRIX) {
			GDBLAS_ERROR("Can not filter real matrix with complex filter coefficients");

			return nullptr;
		}

		input = _conv_implementation(filter, true, &error);
		if (error)
			return nullptr;
	}

	if (factor_m == 1 && factor_n == 1) {
		Ref<GDBlasMat> output = _copy_implementation(&error);
		if (error)
			return nullptr;
		else
			return output;
	}

	error = input->ctx()->downsample(mat->ctx(), factor_m, factor_n);
	if (error)
		return nullptr;

	return mat;
}
