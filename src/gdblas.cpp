#include "gdblas.h"

#include "gdblas_geometry.h"

#define GDBLAS_BIND_CONSTANT(m_const, m_name) \
	godot::ClassDB::bind_integer_constant(get_class_static(), "", #m_name, m_const);

using namespace godot;

void GDBlas::_bind_methods() {
	ClassDB::bind_method(D_METHOD("new_mat", "p_rows", "p_cols"), &GDBlas::new_mat, DEFVAL(-1), DEFVAL(-1));
	ClassDB::bind_method(D_METHOD("new_complex_mat", "p_rows", "p_cols"),
						 &GDBlas::new_complex_mat, DEFVAL(-1), DEFVAL(-1));
	ClassDB::bind_method(D_METHOD("linspace", "p_start", "p_end", "p_count"), &GDBlas::linspace);
	ClassDB::bind_method(D_METHOD("get_version"), &GDBlas::get_version);
	ClassDB::bind_method(D_METHOD("mat_to_image_data", "p_mat_array", "p_channel_width"),
						 &GDBlas::mat_to_image_data, DEFVAL(1));

#ifdef GDBLAS_WITH_GEOMETRY

	ClassDB::bind_method(D_METHOD("area", "p_polygon"), &GDBlas::geom_area);
	ClassDB::bind_method(D_METHOD("buffer", "p_geom", "p_buffer_distance", "p_points_per_join",
								  "p_points_per_end", "p_points_per_circle"),
						 &GDBlas::geom_buffer, DEFVAL(2.0), DEFVAL(4), DEFVAL(4), DEFVAL(4));
	ClassDB::bind_method(D_METHOD("centroid", "p_geom"), &GDBlas::geom_centroid);
	ClassDB::bind_method(D_METHOD("closest_points", "p_geom1", "p_geom2"), &GDBlas::geom_closest_points);
	ClassDB::bind_method(D_METHOD("convex_hull", "p_geom"), &GDBlas::geom_convex_hull);
	ClassDB::bind_method(D_METHOD("correct", "p_geom"), &GDBlas::geom_correct);
	ClassDB::bind_method(D_METHOD("covered_by", "p_geom1", "p_geom2"), &GDBlas::geom_covered_by);
	ClassDB::bind_method(D_METHOD("crosses", "p_geom1", "p_geom2"), &GDBlas::geom_crosses);
	ClassDB::bind_method(D_METHOD("densify", "p_geom", "p_max_distance"), &GDBlas::geom_densify);
	ClassDB::bind_method(D_METHOD("difference", "p_geom1", "p_geom2"), &GDBlas::geom_difference);
	ClassDB::bind_method(D_METHOD("discrete_frechet_distance", "p_geom1", "p_geom2"),
						 &GDBlas::geom_discrete_frechet_distance);
	ClassDB::bind_method(D_METHOD("discrete_hausdorff_distance", "p_geom1", "p_geom2"),
						 &GDBlas::geom_discrete_hausdorff_distance);
	ClassDB::bind_method(D_METHOD("disjoint", "p_geom1", "p_geom2"), &GDBlas::geom_disjoint);
	ClassDB::bind_method(D_METHOD("distance", "p_geom1", "p_geom2"), &GDBlas::geom_distance);
	ClassDB::bind_method(D_METHOD("comparable_distance", "p_geom1", "p_geom2"), &GDBlas::geom_comparable_distance);
	ClassDB::bind_method(D_METHOD("envelope", "p_geom"), &GDBlas::geom_envelope);
	ClassDB::bind_method(D_METHOD("equals", "p_geom1", "p_geom2"), &GDBlas::geom_equals);
	ClassDB::bind_method(D_METHOD("intersection", "p_geom1", "p_geom2"), &GDBlas::geom_intersection);
	ClassDB::bind_method(D_METHOD("intersects", "p_geom1", "p_geom2"), &GDBlas::geom_intersects);
	ClassDB::bind_method(D_METHOD("is_simple", "p_geom"), &GDBlas::geom_is_simple);
	ClassDB::bind_method(D_METHOD("is_valid", "p_geom"), &GDBlas::geom_is_valid);
	ClassDB::bind_method(D_METHOD("length", "p_geom"), &GDBlas::geom_length);
	ClassDB::bind_method(D_METHOD("overlaps", "p_geom1", "p_geom2"), &GDBlas::geom_overlaps);
	ClassDB::bind_method(D_METHOD("perimeter", "p_geom"), &GDBlas::geom_perimeter);
	ClassDB::bind_method(D_METHOD("relation", "p_geom1", "p_geom2"), &GDBlas::geom_relation);
	ClassDB::bind_method(D_METHOD("reverse", "p_geom"), &GDBlas::geom_reverse);
	ClassDB::bind_method(D_METHOD("simplify", "p_geom", "p_max_distance"), &GDBlas::geom_simplify);
	ClassDB::bind_method(D_METHOD("sym_difference", "p_geom1", "p_geom2"), &GDBlas::geom_sym_difference);
	ClassDB::bind_method(D_METHOD("touches", "p_geom1", "p_geom2"), &GDBlas::geom_touches, DEFVAL(Variant()));
	ClassDB::bind_method(D_METHOD("union_", "p_geom1", "p_geom2"), &GDBlas::geom_union_);
	ClassDB::bind_method(D_METHOD("unique", "p_geom1"), &GDBlas::geom_unique);
	ClassDB::bind_method(D_METHOD("within", "p_geom1", "p_geom2"), &GDBlas::geom_within);

#endif // GDBLAS_WITH_GEOMETRY

	GDBLAS_BIND_CONSTANT(GDBlasMat::ERR_GENERAL, ERR_GENERAL);
	GDBLAS_BIND_CONSTANT(GDBlasMat::ERR_INVALID_DIM, ERR_INVALID_DIM);
	GDBLAS_BIND_CONSTANT(GDBlasMat::ERR_SINGULAR_MAT, ERR_SINGULAR_MAT);
	GDBLAS_BIND_CONSTANT(GDBlasMat::ERR_INVALID_INDEX, ERR_INVALID_INDEX);
	GDBLAS_BIND_CONSTANT(GDBlasMat::ERR_INVALID_TYPE, ERR_INVALID_TYPE);
	GDBLAS_BIND_CONSTANT(GDBlasMat::ERR_INVALID_INPUT, ERR_INVALID_INPUT);
	GDBLAS_BIND_CONSTANT(GDBlasMat::ERR_INCOMPATIBLE_DIM, ERR_INCOMPATIBLE_DIM);

	GDBLAS_BIND_CONSTANT(GDBlasMat::NORM_1, NORM_1);
	GDBLAS_BIND_CONSTANT(GDBlasMat::NORM_INF, NORM_INF);
	GDBLAS_BIND_CONSTANT(GDBlasMat::NORM_FRO, NORM_FRO);
	GDBLAS_BIND_CONSTANT(GDBlasMat::REAL_COMPONENT, REAL_COMPONENT);
	GDBLAS_BIND_CONSTANT(GDBlasMat::IMAG_COMPONENT, IMAG_COMPONENT);
	GDBLAS_BIND_CONSTANT(GDBlasMat::BOTH_COMPONENTS, BOTH_COMPONENTS);
}

GDBlas::GDBlas() {
	GDBLAS_V_DEBUG("Created GDBlas: %lu", get_instance_id());
}

GDBlas::~GDBlas() {
	GDBLAS_V_DEBUG("Deleted GDBlas: %lu", get_instance_id());
}

GDBlasMat::Dimension GDBlas::_new_mat_dimension(Variant p_rows, int p_cols, int *error) {
	GDBlasMat::Dimension d;

	Variant::Type t1 = p_rows.get_type();
	int rows = -1;
	int cols = -1;
	if (t1 == Variant::INT) {
		rows = p_rows;
		rows = GDBLAS_MAX(rows, 1);
		cols = p_cols;
	}
	else if (t1 == Variant::VECTOR2I) {
		Vector2i dim = p_rows;
		rows = dim.x;
		cols = p_cols < 0 ? dim.y : p_cols;
	}
	else {
		GDBLAS_ERROR("Invalid parameter");

		*error = GDBlasMat::ERR_INVALID_INPUT;
		return d;
	}

	if (cols < 0)
		cols = rows;

	if (rows <= 0 || cols <= 0) {
		GDBLAS_ERROR("Invalid dimension: (%dx%d)", rows, cols);

		*error = GDBlasMat::ERR_INVALID_DIM;

		return d;
	}

	d.m = rows;
	d.n = cols;

	return d;
}

Variant GDBlas::new_mat(Variant p_rows, int p_cols = -1) {
	int error = 0;
	GDBlasMat::Dimension d = _new_mat_dimension(p_rows, p_cols, &error);
	if (error)
		return Variant();

	Ref< GDBlasMat > mat = GDBlasMat::new_mat(d.m, d.n, GDBlasMat::BLAS_MATRIX, &error);

	if (error)
		return Variant();

	return Variant(mat);
}

Variant GDBlas::new_complex_mat(Variant p_rows, int p_cols = -1) {
	int error = 0;
	GDBlasMat::Dimension d = _new_mat_dimension(p_rows, p_cols, &error);
	if (error)
		return Variant();

	Ref<GDBlasMat> mat = GDBlasMat::new_mat(d.m, d.n, GDBlasMat::BLAS_COMPLEX_MATRIX, &error);

	if (error)
		return Variant();

	return Variant(mat);
}

Variant GDBlas::linspace(GDBlasMat::scalar_t p_start, GDBlasMat::scalar_t p_end, int p_count) {
	int error = 0;
	Ref< GDBlasMat > mat = GDBlasMat::_linspace_implementation(p_start, p_end, p_count, &error);
	if (error) {
		return Variant();
	}

	return Variant(mat);
}

Variant GDBlas::mat_to_image_data(Array p_mat_array, int p_channel_width) {
	if (p_channel_width < 1) {
		GDBLAS_ERROR("Invalid channel width");

		return Variant();
	}

	if (p_channel_width > 2) {
		GDBLAS_ERROR("Channeld width must be less than or equalt to 2");

		return Variant();
	}

	int64_t channel_count = p_mat_array.size();

	PackedByteArray output;

	GDBlasMat::Dimension d;
	for (int64_t i = 0; i < channel_count; ++i) {
		GDBlasMat *mat = GDBlasMat::_cast(p_mat_array[i]);
		if (mat == nullptr) {
			GDBLAS_ERROR("Array entries must be GDBlasMat type");

			return Variant();
		}

		if (mat->get_type() != GDBlasMat::BLAS_MATRIX) {
			GDBLAS_ERROR("All input matrices must be real valued");

			return Variant();
		}

		GDBlasMat::Dimension di = mat->_size();

		if ((d.m && d.m != di.m) || (d.n && d.n != di.n)) {
			GDBLAS_ERROR("All input matrices must have the same dimension");

			return Variant();
		}

		d.m = di.m;
		d.n = di.n;
		if (output.size() == 0) {
			GDBlasMat::s_t total_size = d.m * d.n * channel_count * p_channel_width;
			output.resize(total_size);
		}

		int error = 0;
		if (p_channel_width == 1) {
			error = mat->_pack_implementation<uint8_t>(output, GDBlasMat::REAL_COMPONENT,
													   channel_count, i, false);
		} else if (p_channel_width == 2) {
			error = mat->_pack_implementation<uint16_t>(output, GDBlasMat::REAL_COMPONENT,
														channel_count, i, false);
		}

		if (error) {
			GDBLAS_ERROR("'mat_to_image_data' exits with code %d", error);

			return Variant();
		}
	}

	return output;
}

#ifdef GDBLAS_WITH_GEOMETRY

static bool swap_geom(Variant **geom1, Variant **geom2) {
	Variant *tmp = *geom1;
	*geom1 = *geom2;
	*geom2 = tmp;

	return true;
}

template <typename C>
static bool swap_geom(Variant **geom1, Variant **geom2, C &&condition) {
	Variant *tmp1 = *geom1;
	Variant *tmp2 = *geom2;

	bool cond = condition(*tmp1, *tmp2);

	if (cond) {
		*geom1 = *geom2;
		*geom2 = tmp1;

		return true;
	}

	return false;
}

template <typename T, typename F>
static T PRLV_geom_unary(const Variant &geom, const F &&func, const T &&error_val) {
	int type = geom.get_type();

	switch (type) {
	case Variant::ARRAY: {
		const Array array = geom;
		if (array.is_empty()) {
			GDBLAS_ERROR("Empty polygon geometry");
			break;
		}
		return func(array);
	}
	case Variant::PACKED_VECTOR2_ARRAY: {
		const PackedVector2Array array = geom;
		if (array.is_empty()) {
			GDBLAS_ERROR("Empty line/ring geometry");
			break;
		}
		return func(array);
	}
	case Variant::VECTOR2: {
		const Vector2 vec = geom;
		return func(vec);
	}
	default:
		break;
	}

	return error_val;
}

template <typename T, typename F>
static T PRL_geom_unary(const Variant &geom, const F &&func, const T &&error_val) {
	int type = geom.get_type();

	switch (type) {
	case Variant::ARRAY: {
		const Array array = geom;
		if (array.is_empty()) {
			GDBLAS_ERROR("Empty polygon geometry");
			break;
		}
		return func(array);
	}
	case Variant::PACKED_VECTOR2_ARRAY: {
		const PackedVector2Array array = geom;
		if (array.is_empty()) {
			GDBLAS_ERROR("Empty line/ring geometry");
			break;
		}
		return func(array);
	}
	default:
		break;
	}

	return error_val;
}

inline static int64_t gbl_type_(const Variant &v1, const Variant &v2) {
	return v1.get_type() | ((int64_t) v2.get_type() << 32);
}

#define GBL_TYPE(t1, t2) (((int64_t) Variant:: t1) | (((int64_t) Variant:: t2) << 32))

template <typename T, typename F>
static T PRL_geom_binary_commutative(Variant &geom1, Variant &geom2,
									 const F &&func, const T &&error_val) {
	Variant *g1 = &geom1;
	Variant *g2 = &geom2;
	int64_t type = gbl_type_(*g1, *g2);

	if (type == GBL_TYPE(PACKED_VECTOR2_ARRAY, ARRAY))
		swap_geom(&g1, &g2);

	switch (type) {
	case GBL_TYPE(ARRAY, ARRAY): {
		const Array array1 = *g1;
		const Array array2 = *g2;
		if (array1.is_empty() || array2.is_empty()) {
			GDBLAS_ERROR("Empty polygon geometry");
			break;
		}
		return func(array1, array2);
	}
	case GBL_TYPE(PACKED_VECTOR2_ARRAY, PACKED_VECTOR2_ARRAY): {
		const PackedVector2Array array1 = *g1;
		const PackedVector2Array array2 = *g2;
		if (array1.is_empty() || array2.is_empty()) {
			GDBLAS_ERROR("Empty line/ring geometry");
			break;
		}
		return func(array1, array2);
	}
	case GBL_TYPE(PACKED_VECTOR2_ARRAY, ARRAY):
	case GBL_TYPE(ARRAY, PACKED_VECTOR2_ARRAY): {
		const Array array1 = *g1;
		const PackedVector2Array array2 = *g2;
		if (array1.is_empty()) {
			GDBLAS_ERROR("Empty polygon geometry");
			break;
		}
		if (array2.is_empty()) {
			GDBLAS_ERROR("Empty line/ring geometry");
			break;
		}
		return func(array1, array2);
	}
	default:
		break;
	}

	return error_val;
}

template <typename T, typename F>
static T PRLV_geom_binary_commutative(Variant &geom1, Variant &geom2,
									  const F &&func, const T &&error_val) {
	Variant *g1 = &geom1;
	Variant *g2 = &geom2;
	int64_t type = gbl_type_(*g1, *g2);

	bool cond = type == GBL_TYPE(PACKED_VECTOR2_ARRAY, ARRAY) ||
		type == GBL_TYPE(VECTOR2, ARRAY) ||
		type == GBL_TYPE(VECTOR2, PACKED_VECTOR2_ARRAY);

	if (cond)
		swap_geom(&g1, &g2);

	switch (type) {
	case GBL_TYPE(ARRAY, ARRAY): {
		const Array array1 = *g1;
		const Array array2 = *g2;
		if (array1.is_empty() || array2.is_empty()) {
			GDBLAS_ERROR("Empty polygon geometry");
			break;
		}
		return func(array1, array2);
	}
	case GBL_TYPE(PACKED_VECTOR2_ARRAY, PACKED_VECTOR2_ARRAY): {
		const PackedVector2Array array1 = *g1;
		const PackedVector2Array array2 = *g2;
		if (array1.is_empty() || array2.is_empty()) {
			GDBLAS_ERROR("Empty line/ring geometry");
			break;
		}
		return func(array1, array2);
	}
	case GBL_TYPE(PACKED_VECTOR2_ARRAY, ARRAY):
	case GBL_TYPE(ARRAY, PACKED_VECTOR2_ARRAY): {
		const Array array1 = *g1;
		const PackedVector2Array array2 = *g2;
		if (array1.is_empty()) {
			GDBLAS_ERROR("Empty polygon geometry");
			break;
		}
		if (array2.is_empty()) {
			GDBLAS_ERROR("Empty line/ring geometry");
			break;
		}
		return func(array1, array2);
	}
	case GBL_TYPE(VECTOR2, ARRAY):
	case GBL_TYPE(ARRAY, VECTOR2): {
		const Array array = *g1;
		const Vector2 vec = *g2;

		return func(array, vec);
	}
	case GBL_TYPE(VECTOR2, PACKED_VECTOR2_ARRAY):
	case GBL_TYPE(PACKED_VECTOR2_ARRAY, VECTOR2): {
		const PackedVector2Array array = *g1;
		const Vector2 vec = *g2;

		return func(array, vec);
	}
	default:
		break;
	}

	return error_val;
}

template <typename T, typename F>
static T PRL_geom_binary_noncommutative(Variant &geom1, Variant &geom2,
										 const F &&func, const T &&error_val) {
	Variant *g1 = &geom1;
	Variant *g2 = &geom2;
	int64_t type = gbl_type_(*g1, *g2);

	bool cond = type == GBL_TYPE(PACKED_VECTOR2_ARRAY, ARRAY);

	if (cond)
		swap_geom(&g1, &g2);

	switch (type) {
	case GBL_TYPE(ARRAY, ARRAY): {
		const Array array1 = *g1;
		const Array array2 = *g2;
		if (array1.is_empty() || array2.is_empty()) {
			GDBLAS_ERROR("Empty polygon geometry");
			break;
		}
		return func(array1, array2, cond);
	}
	case GBL_TYPE(PACKED_VECTOR2_ARRAY, PACKED_VECTOR2_ARRAY): {
		const PackedVector2Array array1 = *g1;
		const PackedVector2Array array2 = *g2;
		if (array1.is_empty() || array2.is_empty()) {
			GDBLAS_ERROR("Empty line/ring geometry");
			break;
		}
		return func(array1, array2, cond);
	}
	case GBL_TYPE(PACKED_VECTOR2_ARRAY, ARRAY):
	case GBL_TYPE(ARRAY, PACKED_VECTOR2_ARRAY): {
		const Array array1 = *g1;
		const PackedVector2Array array2 = *g2;
		if (array1.is_empty()) {
			GDBLAS_ERROR("Empty polygon geometry");
			break;
		}
		if (array2.is_empty()) {
			GDBLAS_ERROR("Empty line/ring geometry");
			break;
		}
		return func(array1, array2, cond);
	}
	default:
		break;
	}

	return error_val;
}

template <typename T, typename F>
static T PRLV_geom_binary_noncommutative(Variant &geom1, Variant &geom2,
										 const F &&func, const T &&error_val) {
	Variant *g1 = &geom1;
	Variant *g2 = &geom2;
	int64_t type = gbl_type_(*g1, *g2);

	bool cond = type == GBL_TYPE(PACKED_VECTOR2_ARRAY, ARRAY) ||
		type == GBL_TYPE(VECTOR2, ARRAY) ||
		type == GBL_TYPE(VECTOR2, PACKED_VECTOR2_ARRAY);

	if (cond)
		swap_geom(&g1, &g2);

	switch (type) {
	case GBL_TYPE(ARRAY, ARRAY): {
		const Array array1 = *g1;
		const Array array2 = *g2;
		if (array1.is_empty() || array2.is_empty()) {
			GDBLAS_ERROR("Empty polygon geometry");
			break;
		}
		return func(array1, array2, cond);
	}
	case GBL_TYPE(PACKED_VECTOR2_ARRAY, PACKED_VECTOR2_ARRAY): {
		const PackedVector2Array array1 = *g1;
		const PackedVector2Array array2 = *g2;
		if (array1.is_empty() || array2.is_empty()) {
			GDBLAS_ERROR("Empty line/ring geometry");
			break;
		}
		return func(array1, array2, cond);
	}
	case GBL_TYPE(PACKED_VECTOR2_ARRAY, ARRAY):
	case GBL_TYPE(ARRAY, PACKED_VECTOR2_ARRAY): {
		const Array array1 = *g1;
		const PackedVector2Array array2 = *g2;
		if (array1.is_empty()) {
			GDBLAS_ERROR("Empty polygon geometry");
			break;
		}
		if (array2.is_empty()) {
			GDBLAS_ERROR("Empty line/ring geometry");
			break;
		}
		return func(array1, array2, cond);
	}
	case GBL_TYPE(VECTOR2, ARRAY):
	case GBL_TYPE(ARRAY, VECTOR2): {
		const Array array = *g1;
		const Vector2 vec = *g2;

		return func(array, vec, cond);
	}
	case GBL_TYPE(VECTOR2, PACKED_VECTOR2_ARRAY):
	case GBL_TYPE(PACKED_VECTOR2_ARRAY, VECTOR2): {
		const PackedVector2Array array = *g1;
		const Vector2 vec = *g2;

		return func(array, vec, cond);
	}
	default:
		break;
	}

	return error_val;
}

#undef GBL_TYPE

Variant GDBlas::geom_area(Variant p_geom) {
	return PRL_geom_unary(p_geom, [](const auto &geom) -> GDBlasGeometry::scalar_t {
		return GDBlasGeometry::area(geom);
	}, GDBLAS_NaN);
}

Variant GDBlas::geom_buffer(Variant p_geom, double p_buffer_distance, int p_points_per_join,
							int p_points_per_end, int p_points_per_circle) {
	return PRL_geom_unary(p_geom, [&](const auto &geom) -> Array {
		return GDBlasGeometry::buffer(geom, p_buffer_distance, p_points_per_join, p_points_per_end,
									  p_points_per_circle);
	}, Array());
}

Variant GDBlas::geom_centroid(Variant p_geom) {
	return PRL_geom_unary(p_geom, [](const auto &geom) -> Vector2 {
		return GDBlasGeometry::centroid(geom);
	}, Vector2(GDBLAS_NaN, GDBLAS_NaN));
}

Variant GDBlas::geom_closest_points(Variant p_geom1, Variant p_geom2) {
	return PRLV_geom_binary_noncommutative(p_geom1, p_geom2, [](const auto &geom1, const auto &geom2,
																bool is_swapped) -> PackedVector2Array {
		return GDBlasGeometry::closest_points(geom1, geom2, !is_swapped);
	}, PackedVector2Array());
}

Variant GDBlas::geom_convex_hull(Variant p_geom) {
	return PRL_geom_unary(p_geom, [](const auto &geom) -> PackedVector2Array {
		return GDBlasGeometry::convex_hull(geom);
	}, PackedVector2Array());
}

Variant GDBlas::geom_correct(Variant p_geom) {
	return PRL_geom_unary(p_geom, [](const auto &geom) -> Variant {
		return GDBlasGeometry::correct(geom);
	}, Variant());
}

Variant GDBlas::geom_covered_by(Variant p_geom1, Variant p_geom2) {
	return PRLV_geom_binary_noncommutative(p_geom1, p_geom2, [](const auto &geom1, const auto &geom2,
																bool is_swapped) -> int {
		return GDBlasGeometry::covered_by(geom1, geom2, !is_swapped);
	}, -1);
}

Variant GDBlas::geom_crosses(Variant p_geom1, Variant p_geom2) {
	return PRL_geom_binary_commutative(p_geom1, p_geom2, [](const auto &geom1, const auto &geom2) -> int {
		return GDBlasGeometry::crosses(geom1, geom2);
	}, -1);
}

Variant GDBlas::geom_densify(Variant p_geom, GDBlasMat::scalar_t p_max_distance) {
	return PRL_geom_unary(p_geom, [&](const auto &geom) {
		return GDBlasGeometry::densify(geom, p_max_distance);
	}, Variant());
}

Variant GDBlas::geom_difference(Variant p_geom1, Variant p_geom2) {
	return PRL_geom_binary_noncommutative(p_geom1, p_geom2, [](const auto &geom1, const auto &geom2,
															   bool is_swapped) -> Array {
		return GDBlasGeometry::difference(geom1, geom2, !is_swapped);
	}, Array());
}

Variant GDBlas::geom_discrete_frechet_distance(Variant p_geom1, Variant p_geom2) {
	int type1 = p_geom1.get_type();
	int type2 = p_geom2.get_type();

	if (type1 == Variant::PACKED_VECTOR2_ARRAY && type2 == Variant::PACKED_VECTOR2_ARRAY) {
		const PackedVector2Array array1 = p_geom1;
		const PackedVector2Array array2 = p_geom2;

		if (array1.size() < 1 || array2.size() < 1) {
			GDBLAS_ERROR("Empty line geometry");

			return GDBLAS_NaN;
		}

		return GDBlasGeometry::discrete_frechet_distance(array1, array2);
	} else {
		GDBLAS_ERROR("Undefined pair of geometries, only supports line");
	}

	return GDBLAS_NaN;
}

Variant GDBlas::geom_discrete_hausdorff_distance(Variant p_geom1, Variant p_geom2) {
	int type1 = p_geom1.get_type();
	int type2 = p_geom2.get_type();

	if (type1 == Variant::PACKED_VECTOR2_ARRAY && type2 == Variant::PACKED_VECTOR2_ARRAY) {
		const PackedVector2Array array1 = p_geom1;
		const PackedVector2Array array2 = p_geom2;

		if (array1.size() < 1 || array2.size() < 1) {
			GDBLAS_ERROR("Empty line geometry");

			return GDBLAS_NaN;
		}

		return GDBlasGeometry::discrete_hausdorff_distance(array1, array2);
	} else {
		GDBLAS_ERROR("Undefined pair of geometries, only supports line");
	}

	return GDBLAS_NaN;
}

Variant GDBlas::geom_disjoint(Variant p_geom1, Variant p_geom2) {
	return PRLV_geom_binary_commutative(p_geom1, p_geom2, [](const auto &geom1, const auto &geom2) -> int {
		return GDBlasGeometry::disjoint(geom1, geom2);
	}, -1);
}

static Variant _geom_distance(Variant &p_geom1, Variant &p_geom2, bool p_comparable) {
	return PRLV_geom_binary_commutative(p_geom1, p_geom2, [&](const auto &geom1, const auto &geom2) -> GDBlasGeometry::scalar_t {
		return GDBlasGeometry::distance(geom1, geom2, p_comparable);
	}, GDBLAS_NaN);
}

Variant GDBlas::geom_distance(Variant p_geom1, Variant p_geom2) {
	return _geom_distance(p_geom1, p_geom2, false);
}

Variant GDBlas::geom_comparable_distance(Variant p_geom1, Variant p_geom2) {
	return _geom_distance(p_geom1, p_geom2, true);
}

Variant GDBlas::geom_envelope(Variant p_geom) {
	return PRL_geom_unary(p_geom, [](const auto &geom) -> Rect2 {
		return GDBlasGeometry::envelope(geom);
	}, GDBLAS_NaN_RECT);
}

Variant GDBlas::geom_equals(Variant p_geom1, Variant p_geom2) {
	int type1 = p_geom1.get_type();
	int type2 = p_geom2.get_type();

	if (type1 == Variant::ARRAY && type2 == Variant::ARRAY) {
		Array array1 = p_geom1;
		Array array2 = p_geom2;

		return GDBlasGeometry::equals(array1, array2);
	} else if (type1 == Variant::PACKED_VECTOR2_ARRAY && type2 == Variant::PACKED_VECTOR2_ARRAY) {
		PackedVector2Array array1 = p_geom1;
		PackedVector2Array array2 = p_geom2;

		return GDBlasGeometry::equals(array1, array2);
	} else {
		GDBLAS_ERROR("Undefined geometry");
	}

	return -1;
}

Variant GDBlas::geom_intersection(Variant p_geom1, Variant p_geom2) {
	return PRL_geom_binary_commutative(p_geom1, p_geom2, [](const auto &geom1, const auto &geom2) -> Array {
		return GDBlasGeometry::intersection(geom1, geom2);
	}, Array());
}

Variant GDBlas::geom_intersects(Variant p_geom1, Variant p_geom2) {
	return PRL_geom_binary_commutative(p_geom1, p_geom2, [](const auto &geom1, const auto &geom2) -> int {
		return GDBlasGeometry::intersects(geom1, geom2);
	}, -1);
}

Variant GDBlas::geom_is_simple(Variant p_geom) {
	return PRL_geom_unary(p_geom, [](const auto &geom) -> int {
		return GDBlasGeometry::is_simple(geom);
	}, -1);
}

Variant GDBlas::geom_is_valid(Variant p_geom) {
	return PRL_geom_unary(p_geom, [](const auto &geom) -> int {
		return GDBlasGeometry::is_valid(geom);
	}, -1);
}

Variant GDBlas::geom_length(Variant p_geom) {
	return PRL_geom_unary(p_geom, [](const auto &geom) -> GDBlasGeometry::scalar_t {
		return GDBlasGeometry::length(geom);
	}, GDBLAS_NaN);
}

Variant GDBlas::geom_overlaps(Variant p_geom1, Variant p_geom2) {
	return PRL_geom_binary_commutative(p_geom1, p_geom2, [](const auto &geom1, const auto &geom2) -> int {
		return GDBlasGeometry::overlaps(geom1, geom2);
	}, -1);
}

Variant GDBlas::geom_perimeter(Variant p_geom) {
	return PRL_geom_unary(p_geom, [](const auto &geom) -> GDBlasGeometry::scalar_t {
		return GDBlasGeometry::perimeter(geom);
	}, GDBLAS_NaN);
}

Variant GDBlas::geom_relation(Variant p_geom1, Variant p_geom2) {
	return PRLV_geom_binary_noncommutative(p_geom1, p_geom2, [](const auto &geom1, const auto &geom2,
																bool is_swapped) -> String {
		return GDBlasGeometry::relation(geom1, geom2, !is_swapped);
	}, String());
}

Variant GDBlas::geom_reverse(Variant p_geom) {
	return PRL_geom_unary(p_geom, [](const auto &geom) -> Variant {
		return GDBlasGeometry::reverse(geom);
	}, Variant());
}

Variant GDBlas::geom_simplify(Variant p_geom, GDBlasMat::scalar_t p_max_distance) {
	return PRL_geom_unary(p_geom, [&](const auto &geom) -> Variant {
		return GDBlasGeometry::simplify(geom, p_max_distance);
	}, Variant());
}

Variant GDBlas::geom_sym_difference(Variant p_geom1, Variant p_geom2) {
	return PRL_geom_binary_noncommutative(p_geom1, p_geom2, [](const auto &geom1, const auto &geom2,
															   bool is_swapped) -> Array {
		return GDBlasGeometry::sym_difference(geom1, geom2, !is_swapped);
	}, Array());
}

Variant GDBlas::geom_touches(Variant p_geom1, Variant p_geom2) {
	int type2 = p_geom2.get_type();

	if (type2 == Variant::NIL) {
		return PRL_geom_unary(p_geom1, [&](const auto &geom) -> int {
			return GDBlasGeometry::touches(geom);
		}, -1);
	}

	return PRL_geom_binary_commutative(p_geom1, p_geom2, [](const auto &geom1, const auto &geom2) -> int {
		return GDBlasGeometry::touches(geom1, geom2);
	}, -1);
}

Variant GDBlas::geom_union_(Variant p_geom1, Variant p_geom2) {
	return PRL_geom_binary_commutative(p_geom1, p_geom2, [](const auto &geom1, const auto &geom2) -> Array {
		return GDBlasGeometry::union_(geom1, geom2);
	}, Array());
}

Variant GDBlas::geom_unique(Variant p_geom) {
	return PRL_geom_unary(p_geom, [](const auto &geom) -> Variant {
		return GDBlasGeometry::unique(geom);
	}, Variant());
}

Variant GDBlas::geom_within(Variant p_geom1, Variant p_geom2) {
	return PRLV_geom_binary_noncommutative(p_geom1, p_geom2, [](const auto &geom1, const auto &geom2,
																bool is_swapped) -> int {
		return GDBlasGeometry::within(geom1, geom2, !is_swapped);
	}, -1);
}

#endif // GDBLAS_WITH_GEOMETRY

Variant GDBlas::get_version() {
	return GDBLAS_VERSION;
}
