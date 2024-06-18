#ifndef GDBLAS_MAT_H_
#define GDBLAS_MAT_H_

#define EIGEN_MAX_ALIGN_BYTES 16

#include <cmath>
#include <complex>
#include <vector>

#include <godot_cpp/classes/ref.hpp>
#include <godot_cpp/core/class_db.hpp>
#include <godot_cpp/variant/array.hpp>
#include <godot_cpp/variant/packed_float64_array.hpp>
#include <godot_cpp/variant/variant.hpp>

#include <Eigen/Dense>

#include "defs.h"
#include "utils.h"

#define GDBLAS_NaN std::nan("nan")
#define GDBLAS_MAX(a, b) ((a) > (b) ? (a) : (b))
#define GDBLAS_MIN(a, b) ((a) < (b) ? (a) : (b))
#define GDBLAS_SIZE_CAST(a) ((s_t)GDBLAS_MIN((s_t)GDBLAS_MAX(a, 0), std::numeric_limits<s_t>::max()))

namespace godot {

class GDBlasMat : public RefCounted {
	friend class GDBlas;

	GDCLASS(GDBlasMat, RefCounted)

public:
	typedef size_t s_t;
	typedef double scalar_t;
	typedef std::complex<scalar_t> complex_t;
	typedef Eigen::MatrixX<scalar_t> Matrix;
	typedef Eigen::MatrixX<complex_t> ComplexMatrix;
	typedef std::array<int, 2> index_t;

	template <typename PT = scalar_t>
	struct entry_packed_t {
		static_assert(sizeof(PT) <= sizeof(scalar_t),
					  "sizeof packed type must be less than or equal to sizeof(scalar_t)");

		PT packed_entry;

		operator scalar_t() const {
			return packed_entry;
		}
	};

	template <int W, typename PT = scalar_t>
	struct _entry_real_t {
		scalar_t components[W];

		typedef entry_packed_t<PT> _packed_t;
		void operator=(const _packed_t &other) {
			*components = other;
		}
		operator entry_packed_t<PT>() const {
			entry_packed_t<PT> entry;
			entry.packed_entry = *components;
			return entry;
		}
	};

	typedef _entry_real_t<1> entry_real_t;
	typedef _entry_real_t<2> entry_complex_t;

	enum ErrorDefinitions {
		ERR_GENERAL = std::numeric_limits<short>::lowest(),
		ERR_INVALID_DIM,
		ERR_SINGULAR_MAT,
		ERR_INVALID_INDEX,
		ERR_INVALID_TYPE,
		ERR_INVALID_INPUT,
		ERR_INCOMPATIBLE_DIM,
	};

	enum NormType {
		NORM_1 = 1,
		NORM_INF,
		NORM_FRO,
	};

	enum NumberComponent {
		REAL_COMPONENT = 1,
		IMAG_COMPONENT = 2,
		BOTH_COMPONENTS = REAL_COMPONENT | IMAG_COMPONENT,
	};

	struct Dimension {
		s_t m;
		s_t n;

		Dimension(s_t p_m, s_t p_n) :
				m(p_m), n(p_n) {}
		Dimension() :
				Dimension(0, 0) {}
	};

	static constexpr scalar_t EPS = 1e-16;
	static constexpr index_t INVALID_INDEX{ -1, -1 };

	static GDBlasMat *_cast(Variant &v) {
		Object *tmp = v;
		if (tmp == nullptr)
			return nullptr;

		return Object::cast_to<GDBlasMat>(tmp);
	}

	template <class T>
	static void ignore_unused(T &) {}

protected:
	static void _bind_methods();

	_ALWAYS_INLINE_ bool _is_real_number(Variant &v) {
		Variant::Type t = v.get_type();

		if (t == Variant::FLOAT || t == Variant::INT)
			return true;

		return false;
	}

	_ALWAYS_INLINE_ bool _is_complex_number(Variant &v) {
		Variant::Type t = v.get_type();

		if (t == Variant::VECTOR2)
			return true;

		return false;
	}

	_ALWAYS_INLINE_ Variant _complex_to_variant(complex_t &c) {
		Vector2 vec(c.real(), c.imag());

		return Variant(vec);
	}

	_ALWAYS_INLINE_ Variant _complex_to_variant(complex_t &&c) {
		return _complex_to_variant(c);
	}

	_ALWAYS_INLINE_ complex_t _variant_to_complex(Variant &v) {
		Variant::Type t = v.get_type();
		complex_t c;

		if (t == Variant::FLOAT || t == Variant::INT) {
			c.real(v);
		} else if (t == Variant::VECTOR2) {
			Vector2 vec = v;
			c.real(vec.x);
			c.imag(vec.y);
		} else {
			c.real(GDBLAS_NaN);
		}

		return c;
	}

	_ALWAYS_INLINE_ complex_t _variant_to_complex(Variant &&c) {
		return _variant_to_complex(c);
	}

	_ALWAYS_INLINE_ Dimension _size() {
		Dimension d = ctx()->dimension();

		return d;
	}

	_ALWAYS_INLINE_ bool _is_zero_dim() {
		Dimension d = _size();

		if (d.m == 0 || d.n == 0)
			return true;

		return false;
	}

	_ALWAYS_INLINE_ bool _cmp_size(GDBlasMat *m) {
		if (m == nullptr) {
			return false;
		}

		Dimension d1 = _size();
		Dimension d2 = m->_size();

		return d1.m == d2.m && d1.n == d2.n;
	}

	_ALWAYS_INLINE_ bool _is_scalar() {
		Dimension d = _size();

		if (d.m == 1 && d.n == 1)
			return true;

		return false;
	}

	_ALWAYS_INLINE_ bool _is_square() {
		Dimension d = _size();

		return d.m == d.n;
	}

	_ALWAYS_INLINE_ bool _is_vector() {
		Dimension d = _size();

		return d.m == 1 || d.n == 1;
	}

	_ALWAYS_INLINE_ bool _can_add(GDBlasMat *m) {
		if (m == nullptr) {
			return false;
		}

		if (get_type() == BLAS_MATRIX && m->get_type() == BLAS_COMPLEX_MATRIX) {
			GDBLAS_ERROR("Can not add complex to real matrix");

			return false;
		}

		if (_is_zero_dim() || m->_is_zero_dim())
			return false;

		if (_cmp_size(m))
			return true;

		Dimension d1 = _size();
		Dimension d2 = m->_size();

		GDBLAS_ERROR("Invalid operation on matrices with size (%ux%u) and (%ux%u)",
				d1.m, d1.n, d2.m, d2.n);

		return false;
	}

	_ALWAYS_INLINE_ bool _can_prod(GDBlasMat *m) {
		if (m == nullptr) {
			return false;
		}

		if (_is_zero_dim() || m->_is_zero_dim())
			return false;

		Dimension d1 = _size();
		Dimension d2 = m->_size();

		if (d1.n == d2.m)
			return true;

		GDBLAS_ERROR("Invalid operation on matrices with size (%ux%u) and (%ux%u)",
				d1.m, d1.n, d2.m, d2.n);

		return false;
	}

public:
	enum {
		BLAS_MATRIX = 1,
		BLAS_COMPLEX_MATRIX
	};

	struct Context {
		template <typename T>
		struct BaseMatrixDataContext {
			T m;
		};

		template <typename T, typename U>
		struct BaseMatrixData {
			BaseMatrixDataContext<T> _ctx;

			BaseMatrixData() {
				GDBLAS_V_DEBUG("Created real ctx: %p", this);
			}

			virtual ~BaseMatrixData() {
				GDBLAS_V_DEBUG("Deleted real ctx: %p", this);
			}

			virtual BaseMatrixDataContext<T> &ctx() {
				return _ctx;
			}

			_ALWAYS_INLINE_ T &matrix() {
				return ctx().m;
			}

			_ALWAYS_INLINE_ Dimension dimension() {
				Dimension d(matrix().rows(), matrix().cols());

				return d;
			}

			_ALWAYS_INLINE_ void resize(s_t m, s_t n) {
				matrix().resize(m, n);
			}

			_ALWAYS_INLINE_ U get(s_t i, s_t j) {
				return matrix().coeff(i, j);
			}

			_ALWAYS_INLINE_ T get(s_t i, s_t j, s_t m, s_t n) {
				return matrix().block(i, j, m, n);
			}

			_ALWAYS_INLINE_ void set(U &v, s_t i, s_t j) {
				matrix().coeffRef(i, j) = v;
			}

			_ALWAYS_INLINE_ void set(const T &o, s_t i, s_t j) {
				matrix().block(i, j, o.rows(), o.cols()) = o;
			}

			_ALWAYS_INLINE_ void clear() {
				matrix().setZero();
			}

			_ALWAYS_INLINE_ void transpose() {
				matrix().transposeInPlace();
			}

			void fill(U &s, bool diag = false) {
				for (Eigen::Index i = 0; i < matrix().rows(); ++i) {
					if (diag) {
						if (i >= matrix().cols())
							break;

						matrix().coeffRef(i, i) = s;
					} else {
						matrix().row(i).setConstant(s);
					}
				}
			}

			int inv(Eigen::MatrixX<U> &inverse) {
				Eigen::FullPivLU<T> lu(matrix());

				if (!lu.isInvertible())
					return ERR_SINGULAR_MAT;

				inverse = lu.inverse();

				return 0;
			}

			template <typename Func>
			void elementwise_func(Func &&func) {
				T &m = matrix();

				for (Eigen::Index i = 0; i < m.rows(); ++i) {
					for (auto iter = m.row(i).begin(); iter != m.row(i).end(); ++iter) {
						*iter = func(*iter);
					}
				}
			}

			template <typename Func>
			void elementwise_func_indexed(Func &&func) {
				T &m = matrix();

				for (Eigen::Index i = 0; i < m.rows(); ++i) {
					Eigen::Index j = 0;
					for (auto iter = m.row(i).begin(); iter != m.row(i).end(); ++iter) {
						*iter = func(*iter, i, j++);
					}
				}
			}

			template <typename Func, typename Arg>
			void elementwise_func(Func &&func, Arg &mat) {
				for (Eigen::Index i1 = 0, i2 = 0; i1 < matrix().rows() && i2 < mat.rows(); ++i1, i2++) {
					auto iter1 = matrix().row(i1).begin();
					auto iter2 = mat.row(i2).begin();

					for (; iter1 != matrix().row(i1).end() && iter2 != mat.row(i2).end(); ++iter1, ++iter2) {
						*iter1 = func(*iter1, *iter2);
					}
				}
			}

			template <typename Func, typename Iter>
			U _vector_accumulator(Func &func, Iter &iter, s_t idx1) {
				U l(0);

				s_t idx2 = 0;
				for (auto entry = iter.begin(); entry != iter.end(); ++entry) {
					l = func(l, *entry, idx1, idx2++);
				}

				return l;
			}

			template <typename Func, typename Iter>
			U _vector_accumulator(Func &func, Iter &&iter, s_t idx1) {
				return _vector_accumulator(func, iter, idx1);
			}

			template <typename Func, typename ToScalarFunc>
			int _matrix_accumulator(Func &func, T &out, int axis, ToScalarFunc &to_scalar) {
				T &m = matrix();

				if (axis <= 0) {
					out.resize(1, m.cols());
					for (Eigen::Index idx1 = 0; idx1 < out.cols(); ++idx1) {
						out.coeffRef(0, idx1) = _vector_accumulator(func, m.col(idx1), idx1);
					}
				} else if (axis == 1) {
					out.resize(m.rows(), 1);
					for (Eigen::Index idx1 = 0; idx1 < out.rows(); ++idx1) {
						out.coeffRef(idx1, 0) = _vector_accumulator(func, m.row(idx1), idx1);
					}
				}

				if (axis < 0) {
					U l = _vector_accumulator(to_scalar, out.row(0), 0);

					out.resize(1, 1);
					out(0, 0) = l;
				}

				return 0;
			}

			template <typename Func>
			int _matrix_accumulator(Func &func, T &out, int axis) {
				auto lsum = [](U &a, U &b, s_t i, s_t j) {
					ignore_unused(i);
					ignore_unused(j);

					return a + b;
				};

				return _matrix_accumulator(func, out, axis, lsum);
			}

			int integrate(T &out, int axis) {
				auto lsum = [](U &a, U &b, s_t i, s_t j) {
					ignore_unused(i);
					ignore_unused(j);

					return a + b;
				};

				return _matrix_accumulator(lsum, out, axis);
			}

			int mean(T &out, int axis) {
				Dimension d = dimension();

				scalar_t s = 1.0;

				if (axis == 0)
					s /= d.m;
				else if (axis == 1)
					s /= d.n;
				else
					s /= d.m * d.n;

				int error = integrate(out, axis);
				if (error)
					return error;

				out *= s;

				return 0;
			}

			template <typename Func>
			int _minmax(T &out, int axis, Func &&func, scalar_t fill_with) {
				Dimension d = dimension();

				std::vector<U> max_val(1);
				if (axis == 0)
					max_val.resize(d.n);
				else if (axis == 1)
					max_val.resize(d.m);

				std::fill(max_val.begin(), max_val.end(), fill_with);

				auto lmax = [&](U &a, U &b, s_t i, s_t j) {
					ignore_unused(a);
					ignore_unused(j);
					auto &m = max_val[GDBLAS_MIN(max_val.size() - 1, i)];

					if (func(std::real(b), std::real(m)))
						m = b;

					return m;
				};

				return _matrix_accumulator(lmax, out, axis, lmax);
			}

			_ALWAYS_INLINE_ int max(T &out, int axis) {
				return _minmax(
						out, axis, [](scalar_t &&a, scalar_t &&b) -> bool { return a > b; },
						std::numeric_limits<scalar_t>::lowest());
			}

			_ALWAYS_INLINE_ int min(T &out, int axis) {
				return _minmax(
						out, axis, [](scalar_t &&a, scalar_t &&b) -> bool { return a < b; },
						std::numeric_limits<scalar_t>::max());
			}

			template <typename Func>
			int _argminmax(std::vector<index_t> &arg_minmax, int axis, Func &&func, scalar_t fill_with) {
				Dimension d = dimension();

				T out;
				std::vector<U> minmax_val(1);

				if (axis == 0) {
					minmax_val.resize(d.n);
					arg_minmax.resize(d.n);
				} else if (axis == 1) {
					minmax_val.resize(d.m);
					arg_minmax.resize(d.m);
				}

				std::fill(minmax_val.begin(), minmax_val.end(), fill_with);
				std::fill(arg_minmax.begin(), arg_minmax.end(), INVALID_INDEX);

				auto lminmax = [&](U &a, U &b, s_t i, s_t j) {
					ignore_unused(a);
					auto &m = minmax_val[GDBLAS_MIN(minmax_val.size() - 1, i)];
					index_t &idx = arg_minmax[GDBLAS_MIN(arg_minmax.size() - 1, i)];

					if (func(std::real(b), std::real(m))) {
						m = b;
						if (axis <= 0) {
							idx[0] = j;
							idx[1] = i;
						} else if (axis == 1) {
							idx[0] = i;
							idx[1] = j;
						}
					}

					return m;
				};

				return _matrix_accumulator(lminmax, out, axis, lminmax);
			}

			_ALWAYS_INLINE_ int argmax(std::vector<index_t> &arg_max, int axis) {
				return _argminmax(
						arg_max, axis, [](scalar_t &&a, scalar_t &&b) -> bool { return a > b; },
						std::numeric_limits<scalar_t>::lowest());
			}

			_ALWAYS_INLINE_ int argmin(std::vector<index_t> &arg_min, int axis) {
				return _argminmax(
						arg_min, axis, [](scalar_t &&a, scalar_t &&b) -> bool { return a < b; },
						std::numeric_limits<scalar_t>::max());
			}

			U l1_norm(int *error) {
				U l = 0.0;
				T tmp;

				auto labssum = [](U &a, U &b, s_t i, s_t j) {
					ignore_unused(i);
					ignore_unused(j);

					return a + std::abs(b);
				};

				*error = _matrix_accumulator(labssum, tmp, 0);
				if (*error)
					return l;

				auto cmp = [](U &a, U &b) {
					return std::real(a) < std::real(b);
				};

				auto ltmp = std::max_element(tmp.row(0).begin(), tmp.row(0).end(), cmp);

				return *ltmp;
			}

			U linf_norm(int *error) {
				U l = 0.0;
				T tmp;

				auto labssum = [](U &a, U &b, s_t i, s_t j) {
					ignore_unused(i);
					ignore_unused(j);

					return a + std::abs(b);
				};

				*error = _matrix_accumulator(labssum, tmp, 1);
				if (*error)
					return l;

				auto cmp = [](U &a, U &b) {
					return std::real(a) < std::real(b);
				};

				auto ltmp = std::max_element(tmp.col(0).begin(), tmp.col(0).end(), cmp);

				return *ltmp;
			}

			U fro_norm(int *error) {
				U l = 0.0;
				T tmp;

				auto lsqsum = [](U &a, U &b, s_t i, s_t j) {
					ignore_unused(i);
					ignore_unused(j);

					complex_t c = b;
					return a + std::norm(c);
				};

				*error = _matrix_accumulator(lsqsum, tmp, -1);
				if (*error)
					return l;

				return std::sqrt(tmp(0, 0));
			}

			template <typename T1, typename T2>
			int conv(T1 *mat1, T2 *mat2, bool same) {
				Dimension d1 = mat1->dimension();
				Dimension d2 = mat2->dimension();
				Dimension d3 = dimension();

				int m1 = d1.m;
				int n1 = d1.n;
				int m2 = d2.m;
				int n2 = d2.n;
				int i1_0 = same ? m2 / 2 : 0;
				int i2_0 = same ? n2 / 2 : 0;
				int i1_1 = same ? m1 + m2 / 2 : d3.m;
				int i2_1 = same ? n1 + n2 / 2 : d3.n;
				int yi1_offset = same ? m2 / 2 : 0;
				int yi2_offset = same ? n2 / 2 : 0;

				for (int i1 = i1_0; i1 < i1_1; ++i1) {
					for (int i2 = i2_0; i2 < i2_1; ++i2) {
						int j1_0 = GDBLAS_MAX(i1 - m1 + 1, 0);
						int j2_0 = GDBLAS_MAX(i2 - n1 + 1, 0);
						int j1_1 = GDBLAS_MIN(m2, i1 + 1);
						int j2_1 = GDBLAS_MIN(n2, i2 + 1);

						U yi = 0;

						for (int j1 = j1_0; j1 < j1_1; ++j1) {
							for (int j2 = j2_0; j2 < j2_1; ++j2) {
								yi += mat1->get(i1 - j1, i2 - j2) * mat2->get(j1, j2);
							}
						}

						set(yi, i1 - yi1_offset, i2 - yi2_offset);
					}
				}

				return 0;
			}

			// Can cast std::complex to a double pointer to access real and imag
			// parts as a C array. Please, check references below:
			// https://www.fftw.org/fftw3_doc/Complex-numbers.html#Complex-numbers
			// https://www.open-std.org/jtc1/sc22/wg21/docs/papers/2002/n1388.pdf
			template <typename T1, typename T2 = T1, typename PT, typename Func>
			int pack_unpack(Func &&func, PT &packed_data, bool imag, int step, int offset = 0) {
				const int imag_offset = imag;

				T1 *pdata = reinterpret_cast<T1 *>(packed_data.ptrw());
				pdata += offset;

				for (Eigen::Index i = 0; i < matrix().rows(); ++i) {
					for (auto iter = matrix().row(i).begin(); iter != matrix().row(i).end(); ++iter) {
						U &val = *iter;
						func(pdata, reinterpret_cast<T2 *>(&val) + imag_offset);
						pdata += step;
					}
				}

				return 0;
			}

			template <typename T1, typename T2 = T1, typename PT>
			int pack(PT &packed_data, bool imag = false, int step = 1, int offset = 0) {
				return pack_unpack<T1, T2>([](T1 *a, T2 *b) { *a = *b; }, packed_data, imag, step, offset);
			}

			template <typename T1, typename T2 = T1, typename PT>
			int unpack(PT &packed_data, bool imag = false, int step = 1, int offset = 0) {
				return pack_unpack<T1, T2>([](T1 *a, T2 *b) { *b = *a; }, packed_data, imag, step, offset);
			}

			void downsample(T &output_mat, int factor_m, int factor_n) {
				s_t m = output_mat.rows();
				s_t n = output_mat.cols();
				output_mat = matrix()(Eigen::seqN(0, m, factor_m), Eigen::seqN(0, n, factor_n));
			}
		};

		struct ComplexMatrixData : public BaseMatrixData<ComplexMatrix, complex_t> {
			ComplexMatrixData() {
				GDBLAS_V_DEBUG("Created complex ctx: %p", this);
			}

			~ComplexMatrixData() {
				GDBLAS_V_DEBUG("Deleted complex ctx: %p", this);
			}

			BaseMatrixDataContext<ComplexMatrix> &ctx() {
				return _ctx;
			}
		};

		typedef BaseMatrixData<Matrix, scalar_t> RealMatrixData;
		int type;

		union {
			RealMatrixData *r;
			ComplexMatrixData *c;
		} mdata;

		Context(int p_type) :
				type(p_type) {
			switch (type) {
				case BLAS_MATRIX: {
					mdata.r = new RealMatrixData;
				} break;
				case BLAS_COMPLEX_MATRIX: {
					mdata.c = new ComplexMatrixData;
				} break;
				default:
					mdata.r = nullptr;
					mdata.c = nullptr;
					break;
			}

			GDBLAS_V_DEBUG("Created ctx: %p", this);
		}

		Context() :
				Context(BLAS_MATRIX) {}

		~Context() {
			if (type == BLAS_MATRIX)
				delete get_real_mdata();
			else if (type == BLAS_COMPLEX_MATRIX)
				delete get_complex_mdata();

			GDBLAS_V_DEBUG("Deleted ctx: %p", this);
		}

		_ALWAYS_INLINE_ void get_mdata_ptr(RealMatrixData **ptr) {
			*ptr = mdata.r;
		}

		_ALWAYS_INLINE_ void get_mdata_ptr(ComplexMatrixData **ptr) {
			*ptr = mdata.c;
		}

		template <typename MT>
		_ALWAYS_INLINE_ MT *get_mdata() {
			static_assert(std::is_base_of<RealMatrixData, MT>::value || std::is_base_of<ComplexMatrixData, MT>::value,
					"T is not derived from MatrixData");

			MT *ptr = nullptr;
			get_mdata_ptr(&ptr);

			return ptr;
		}

		_ALWAYS_INLINE_ RealMatrixData *get_default_mdata() {
			return get_mdata<RealMatrixData>();
		}

		_ALWAYS_INLINE_ RealMatrixData *get_real_mdata() {
			return get_default_mdata();
		}

		_ALWAYS_INLINE_ ComplexMatrixData *get_complex_mdata() {
			return get_mdata<ComplexMatrixData>();
		}

		Dimension dimension() {
			if (type == BLAS_MATRIX)
				return get_real_mdata()->dimension();
			else if (type == BLAS_COMPLEX_MATRIX)
				return get_complex_mdata()->dimension();

			return Dimension();
		}

		template <typename T1>
		_ALWAYS_INLINE_ void resize(s_t m, s_t n) {
			return get_mdata<T1>()->resize(m, n);
		}

		template <typename T1>
		_ALWAYS_INLINE_ void clear() {
			get_mdata<T1>()->clear();
		}

		_ALWAYS_INLINE_ complex_t getc(s_t i, s_t j) {
			if (type == BLAS_COMPLEX_MATRIX) {
				return get_complex_mdata()->get(i, j);
			}

			complex_t c;
			c.real(get_real_mdata()->get(i, j));

			return c;
		}

		_ALWAYS_INLINE_ scalar_t get(s_t i, s_t j) {
			if (type == BLAS_COMPLEX_MATRIX)
				return GDBLAS_NaN;

			return get_real_mdata()->get(i, j);
		}

		_ALWAYS_INLINE_ int set(complex_t &v, s_t i, s_t j) {
			if (type == BLAS_COMPLEX_MATRIX) {
				get_complex_mdata()->set(v, i, j);

				return 0;
			}

			return ERR_INVALID_INDEX;
		}

		_ALWAYS_INLINE_ int set(complex_t &&v, s_t i, s_t j) {
			return set(v, i, j);
		}

		_ALWAYS_INLINE_ int set(scalar_t v, s_t i, s_t j) {
			if (type == BLAS_COMPLEX_MATRIX) {
				complex_t c;
				c.real(v);

				return set(c, i, j);
			}

			get_real_mdata()->set(v, i, j);

			return 0;
		}

		int get(Context *other, s_t i, s_t j) {
			Dimension d = dimension();

			if (type == BLAS_COMPLEX_MATRIX && other->type == BLAS_COMPLEX_MATRIX) {
				get_complex_mdata()->matrix() = other->get_complex_mdata()->get(i, j, d.m, d.n);
			} else if (type == BLAS_COMPLEX_MATRIX && other->type == BLAS_MATRIX) {
				get_complex_mdata()->matrix() = other->get_real_mdata()->get(i, j, d.m, d.n);
			} else if (type == BLAS_MATRIX && other->type == BLAS_MATRIX) {
				get_real_mdata()->matrix() = other->get_real_mdata()->get(i, j, d.m, d.n);
			} else {
				return ERR_INVALID_TYPE;
			}

			return 0;
		}

		int set(Context *other, s_t i, s_t j) {
			if (type == BLAS_COMPLEX_MATRIX && other->type == BLAS_COMPLEX_MATRIX) {
				get_complex_mdata()->set(other->get_complex_mdata()->matrix(), i, j);
			} else if (type == BLAS_COMPLEX_MATRIX && other->type == BLAS_MATRIX) {
				get_complex_mdata()->set(other->get_real_mdata()->matrix(), i, j);
			} else if (type == BLAS_MATRIX && other->type == BLAS_MATRIX) {
				get_real_mdata()->set(other->get_real_mdata()->matrix(), i, j);
			} else {
				GDBLAS_ERROR("Can not place a complex entry to a real matrix");

				return ERR_INVALID_TYPE;
			}

			return 0;
		}

		template <typename T1, typename T2>
		_ALWAYS_INLINE_ void add(Context *other) {
			get_mdata<T1>()->matrix() += other->get_mdata<T2>()->matrix();
		}

		template <typename T1, typename T2>
		_ALWAYS_INLINE_ void sub(Context *other) {
			get_mdata<T1>()->matrix() -= other->get_mdata<T2>()->matrix();
		}

		template <typename T1, typename T2>
		void eq(Context *other) {
			get_mdata<T1>()->matrix() = other->get_mdata<T2>()->matrix();
		}

		template <typename T1, typename T2>
		bool is_eq(Context *other, scalar_t eps, int norm_type) {
			T1 tmp_mat_data;

			Dimension d = get_mdata<T1>()->dimension();
			tmp_mat_data.resize(d.m, d.n);
			tmp_mat_data.matrix() = get_mdata<T1>()->matrix() - other->get_mdata<T2>()->matrix();

			int error = 0;
			scalar_t l = std::numeric_limits<scalar_t>::max();
			if (norm_type == NORM_1)
				l = std::real(tmp_mat_data.l1_norm(&error));
			else if (norm_type == NORM_INF)
				l = std::real(tmp_mat_data.linf_norm(&error));
			else if (norm_type == NORM_FRO)
				l = std::real(tmp_mat_data.fro_norm(&error));

			if (!error && l < eps)
				return true;

			return false;
		}

		template <typename T1, typename U>
		void muls(U s) {
			get_mdata<T1>()->matrix() *= s;
		}

		template <typename T1>
		void transpose() {
			get_mdata<T1>()->transpose();
		}

		template <typename T1, typename U>
		void fill(U &s, bool diag = false) {
			get_mdata<T1>()->fill(s, diag);
		}

		void real(Context *other, bool to) {
			if (to)
				other->get_real_mdata()->matrix() = get_complex_mdata()->matrix().real();
			else
				get_complex_mdata()->matrix().real() = other->get_real_mdata()->matrix();
		}

		void imag(Context *other, bool to) {
			if (to)
				other->get_real_mdata()->matrix() = get_complex_mdata()->matrix().imag();
			else
				get_complex_mdata()->matrix().imag() = other->get_real_mdata()->matrix();
		}

		void set_real_or_imag(Context *other, bool p_real, bool to) {
			if (p_real)
				real(other, to);
			else if (!p_real)
				imag(other, to);
		}

		void linspace(scalar_t start, scalar_t end, s_t count) {
			scalar_t dx = (end - start) / (scalar_t)(count - 1);

			Matrix &m = get_real_mdata()->matrix();

			for (Eigen::Index j = 0; j < m.cols(); ++j) {
				scalar_t x = start;

				auto iter = m.col(j).begin();
				for (; iter != m.col(j).end() - 1; ++iter) {
					*iter = x;
					x += dx;
				}
				*iter = end;
			}
		}

		void prod(Context *a, Context *b) {
			if (type == BLAS_MATRIX) {
				get_real_mdata()->matrix() = a->get_real_mdata()->matrix() * b->get_real_mdata()->matrix();
			} else if (type == BLAS_COMPLEX_MATRIX) {
				ComplexMatrixData *ml = get_complex_mdata();

				if (a->type == BLAS_COMPLEX_MATRIX && b->type == BLAS_COMPLEX_MATRIX) {
					ml->matrix() = a->get_complex_mdata()->matrix() * b->get_complex_mdata()->matrix();
				} else if (a->type == BLAS_MATRIX && b->type == BLAS_COMPLEX_MATRIX) {
					ml->matrix() = a->get_real_mdata()->matrix() * b->get_complex_mdata()->matrix();
				}
			}
		}

		int inv(Context *a) {
			int result = 0;

			if (type == BLAS_MATRIX) {
				result = get_real_mdata()->inv(a->get_real_mdata()->matrix());
			} else if (type == BLAS_COMPLEX_MATRIX) {
				result = get_complex_mdata()->inv(a->get_complex_mdata()->matrix());
			}

			return result;
		}

		int integrate(Context *a, int axis = -1) {
			if (type == BLAS_MATRIX && a->type == BLAS_MATRIX) {
				return get_real_mdata()->integrate(a->get_real_mdata()->matrix(), axis);
			} else if (type == BLAS_COMPLEX_MATRIX && a->type == BLAS_COMPLEX_MATRIX) {
				return get_complex_mdata()->integrate(a->get_complex_mdata()->matrix(), axis);
			}

			return ERR_INVALID_TYPE;
		}

		int mean(Context *a, int axis = -1) {
			if (type == BLAS_MATRIX && a->type == BLAS_MATRIX) {
				return get_real_mdata()->mean(a->get_real_mdata()->matrix(), axis);
			} else if (type == BLAS_COMPLEX_MATRIX && a->type == BLAS_COMPLEX_MATRIX) {
				return get_complex_mdata()->mean(a->get_complex_mdata()->matrix(), axis);
			}

			return ERR_INVALID_TYPE;
		}

		int min_max(Context *a, bool is_min, int axis = -1) {
			if (type == BLAS_MATRIX && a->type == BLAS_MATRIX) {
				if (is_min)
					return get_real_mdata()->min(a->get_real_mdata()->matrix(), axis);
				else
					return get_real_mdata()->max(a->get_real_mdata()->matrix(), axis);
			} else if (type == BLAS_COMPLEX_MATRIX && a->type == BLAS_COMPLEX_MATRIX) {
				if (is_min)
					return get_complex_mdata()->min(a->get_complex_mdata()->matrix(), axis);
				else
					return get_complex_mdata()->max(a->get_complex_mdata()->matrix(), axis);
			}

			return ERR_INVALID_TYPE;
		}

		int arg_min_max(std::vector<index_t> &arg, bool is_min, int axis = -1) {
			if (type == BLAS_MATRIX) {
				if (is_min)
					return get_real_mdata()->argmin(arg, axis);
				else
					return get_real_mdata()->argmax(arg, axis);
			} else if (type == BLAS_COMPLEX_MATRIX) {
				if (is_min)
					return get_complex_mdata()->argmin(arg, axis);
				else
					return get_complex_mdata()->argmax(arg, axis);
			}

			return ERR_INVALID_TYPE;
		}

		scalar_t norm(int norm_type, int *error) {
			scalar_t l = 0.0;

			if (type == BLAS_MATRIX) {
				if (norm_type == NORM_1)
					l = get_real_mdata()->l1_norm(error);
				else if (norm_type == NORM_INF)
					l = get_real_mdata()->linf_norm(error);
				else if (norm_type == NORM_FRO)
					l = get_real_mdata()->fro_norm(error);
				else
					*error = ERR_INVALID_INPUT;
			} else if (type == BLAS_COMPLEX_MATRIX) {
				if (norm_type == NORM_1)
					l = get_complex_mdata()->l1_norm(error).real();
				else if (norm_type == NORM_INF)
					l = get_complex_mdata()->linf_norm(error).real();
				else if (norm_type == NORM_FRO)
					l = get_complex_mdata()->fro_norm(error).real();
				else
					*error = ERR_INVALID_INPUT;
			}

			return l;
		}

		int conv(Context *a, Context *b, bool same) {
			int error = 0;

			if (type == BLAS_MATRIX) {
				if (a->type == BLAS_MATRIX && b->type == BLAS_MATRIX)
					error = get_real_mdata()->conv(a->get_real_mdata(), b->get_real_mdata(), same);
			} else if (type == BLAS_COMPLEX_MATRIX) {
				if (a->type == BLAS_COMPLEX_MATRIX && b->type == BLAS_MATRIX)
					error = get_complex_mdata()->conv(a->get_complex_mdata(), b->get_real_mdata(), same);
				else if (a->type == BLAS_COMPLEX_MATRIX && b->type == BLAS_COMPLEX_MATRIX)
					error = get_complex_mdata()->conv(a->get_complex_mdata(), b->get_complex_mdata(), same);
				else if (a->type == BLAS_MATRIX && b->type == BLAS_COMPLEX_MATRIX)
					error = get_complex_mdata()->conv(a->get_real_mdata(), b->get_complex_mdata(), same);
				else
					error = ERR_INVALID_TYPE;
			} else {
				error = ERR_INVALID_TYPE;
			}

			return error;
		}

		int downsample(Context *output, int factor_m, int factor_n) {
			if (type == BLAS_MATRIX && output->type == BLAS_MATRIX)
				get_real_mdata()->downsample(output->get_real_mdata()->matrix(), factor_m, factor_n);
			else if (type == BLAS_COMPLEX_MATRIX && output->type == BLAS_COMPLEX_MATRIX)
				get_complex_mdata()->downsample(output->get_complex_mdata()->matrix(), factor_m, factor_n);
			else
				return ERR_INVALID_TYPE;

			return 0;
		}
	};

	Context *m_ctx;

#ifdef GDBLAS_WITH_ODE
	double m_prev_end_time;
#endif

	GDBlasMat() :
			m_ctx(nullptr) {
#ifdef GDBLAS_WITH_ODE
		m_prev_end_time = 0.0;
#endif

		GDBLAS_V_DEBUG("Created GDBlasMat: %lu", get_instance_id());
	}

	~GDBlasMat() {
		delete m_ctx;

		GDBLAS_V_DEBUG("Deleted GDBlasMat: %lu", get_instance_id());
	}

	_ALWAYS_INLINE_ virtual Context *ctx() {
		if (m_ctx == nullptr) {
			m_ctx = new Context(get_type());
		}

		return m_ctx;
	}

	_ALWAYS_INLINE_ bool init(int m, int n, int type) {
		if (m_ctx != nullptr) {
			return false;
		}

		m_ctx = new Context(type);

		if (m_ctx != nullptr) {
			_resize_implementation(m, n);

			return true;
		}

		return false;
	}

	_ALWAYS_INLINE_ int get_type() {
		return ctx()->type;
	}

	static Ref<GDBlasMat> new_mat(s_t m, s_t n, int type, int *error) {
		Ref<GDBlasMat> mat;
		mat.instantiate();

		if (!mat->init(m, n, type)) {
			mat.unref();

			*error = ERR_GENERAL;
		}

		*error = 0;

		return mat;
	}

	static Ref<GDBlasMat> _linspace_implementation(scalar_t start, scalar_t end,
			int count, int *error) {
		if (count < 2) {
			*error = ERR_INVALID_DIM;

			return nullptr;
		}

		Ref<GDBlasMat> mat = new_mat(count, 1, BLAS_MATRIX, error);
		if (*error)
			return mat;

		mat->ctx()->linspace(start, end, count);

		return mat;
	}

	template <typename ET = scalar_t, typename T>
	int _pack_implementation(T &packed_data, int component, int step, int offset = 0, bool resize = true) {
		Dimension d = _size();
		s_t packed_size = d.m * d.n;
		int type = get_type();
		int result = 0;

		typedef entry_packed_t<ET> _packed_t;
		typedef _entry_real_t<1, ET> _entry_t;

		if (type == BLAS_MATRIX) {
			if (resize)
				packed_data.resize(packed_size);

			if (component == IMAG_COMPONENT) {
				return 0;
			}

			result = ctx()->get_real_mdata()->pack<_packed_t, _entry_t>(packed_data, false,
					step, offset);
		} else if (type == BLAS_COMPLEX_MATRIX) {
			if (component == BOTH_COMPONENTS)
				packed_size *= 2;

			if (resize)
				packed_data.resize(packed_size);

			if (component == REAL_COMPONENT)
				result = ctx()->get_complex_mdata()->pack<_packed_t, _entry_t>(packed_data, false,
						step, offset);
			else if (component == IMAG_COMPONENT)
				result = ctx()->get_complex_mdata()->pack<_packed_t, _entry_t>(packed_data, true,
						step, offset);
			else if (component == BOTH_COMPONENTS)
				result = ctx()->get_complex_mdata()->pack<entry_complex_t>(packed_data, false,
						step, offset);
		} else {
			return ERR_INVALID_TYPE;
		}

		return result;
	}

	template <typename ET = scalar_t, typename T>
	int _unpack_implementation(T &packed_data, int component, int step, int offset = 0) {
		int type = get_type();
		int result = 0;

		typedef entry_packed_t<ET> _packed_t;
		typedef _entry_real_t<1, ET> _entry_t;

		if (type == BLAS_MATRIX) {
			if (component == IMAG_COMPONENT) {
				return 0;
			}

			result = ctx()->get_real_mdata()->unpack<_packed_t, _entry_t>(packed_data, false,
					step, offset);
		} else if (type == BLAS_COMPLEX_MATRIX) {
			if (component == REAL_COMPONENT)
				result = ctx()->get_complex_mdata()->unpack<_packed_t, _entry_t>(packed_data, false,
						step, offset);
			else if (component == IMAG_COMPONENT)
				result = ctx()->get_complex_mdata()->unpack<_packed_t, _entry_t>(packed_data, true,
						step, offset);
			else if (component == BOTH_COMPONENTS)
				result = ctx()->get_complex_mdata()->unpack<entry_complex_t>(packed_data, false,
						step, offset);
		} else {
			return ERR_INVALID_TYPE;
		}

		return result;
	}

	void _resize_implementation(s_t m, s_t n);
	int _add_implementation(GDBlasMat *other);
	int _sub_implementation(GDBlasMat *other);
	void _transpose_implementation();
	void _hermitian_implementation();
	int _eq_implementation(GDBlasMat *other);
	bool _is_eq_implementation(GDBlasMat *other, scalar_t eps = EPS, int norm_type = NORM_FRO);
	int _muls_implementation(scalar_t s);
	int _muls_implementation(complex_t s);
	Ref<GDBlasMat> _copy_implementation(int *error);
	int _fill_implementation(scalar_t s, bool diag = false);
	int _fill_implementation(complex_t &s, bool diag = false);
	void _conj_implementation();
	Ref<GDBlasMat> _real_implementation(int *error);
	int _set_real_implementation(GDBlasMat *other);
	Ref<GDBlasMat> _imag_implementation(int *error);
	int _set_imag_implementation(GDBlasMat *other);
	Ref<GDBlasMat> _prod_implementation(GDBlasMat *other, int *error);
	Ref<GDBlasMat> _inv_implementation(int *error);
	Ref<GDBlasMat> _integrate_implementation(int axis, int *error);
	Ref<GDBlasMat> _mean_implementation(int axis, int *error);
	Ref<GDBlasMat> _min_max_implementation(int axis, int *error, bool is_min = true);
	Array _arg_min_max_implementation(int axis, int *error, bool is_min = true);
	scalar_t _norm_implementation(int norm_type, int *error);
	Ref<GDBlasMat> _conv_implementation(GDBlasMat *other, bool same, int *error);
	int _set_scalar(Variant &val, int i, int j);
	int _set_submatrix(Variant &val, int i, int j);
	Ref<GDBlasMat> _downsample_implementation(int factor_m, int factor_n,
			GDBlasMat *filter = nullptr);

	Variant resize(Variant p_m, int n);
	Variant copy();
	Variant size();
	Variant get(int i, int j, int m = -1, int n = -1);
	Variant set(Variant val, int i = -1, int j = -1);
	Variant add(Variant other);
	Variant sub(Variant other);
	void transpose();
	void hermitian();
	Variant is_eq(Variant other, scalar_t p_eps = EPS, int p_norm_type = NORM_FRO);
	Variant mul(Variant other);
	Variant div(Variant other);
	Variant fill(Variant other);
	Variant eye(Variant p_val);
	Variant reset();
	void conj();
	Variant real(Variant p_matrix);
	Variant imag(Variant p_matrix);
	Variant prod(Variant other);
	Variant inv();
	Variant to_array();
	Variant from_array(Array p_array);
	Variant sin();
	Variant cos();
	Variant abs();
	Variant exp();
	Variant log();
	Variant log10();
	Variant log2();
	Variant sqrt();
	Variant cbrt();
	Variant tan();
	Variant asin();
	Variant acos();
	Variant atan();
	Variant sinh();
	Variant cosh();
	Variant tanh();
	Variant atanh();
	Variant erf();
	Variant erfc();
	Variant tgamma();
	Variant lgamma();
	Variant ceil();
	Variant floor();
	Variant trunc();
	Variant round();
	Variant pow(Variant p_exponent);
	Variant integrate(int axis = -1);
	Variant norm(int p_norm_type);
	Variant mean(int axis = -1);
	Variant min(int axis = -1);
	Variant max(int axis = -1);
	Variant argmin(int axis = -1);
	Variant argmax(int axis = -1);
	Variant unary_func(Callable p_func, Variant p_args, bool p_indexed = false);
	Variant conv(Variant p_other, bool p_same = false);
	Variant pack(int p_component);
	Variant unpack(Variant p_packed_data, int p_component = BOTH_COMPONENTS,
			int p_step = 1, int p_offset = 0);
	Variant downsample(int p_factor_m, int p_factor_n, Variant p_filter);

#ifdef GDBLAS_WITH_ODE
	Variant eval_ode(Callable p_f, double p_dt, double p_max_step = 1e-2);
#endif
};
} //namespace godot

#endif // GDBLAS_MAT_H_
