#include "gdblas.h"


using namespace godot;


#define GDBLAS_BIND_CONSTANT(m_const, m_name) \
	godot::ClassDB::bind_integer_constant(get_class_static(), "", #m_name, m_const);



void GDBlas::_bind_methods() {
	ClassDB::bind_method(D_METHOD("new_mat", "p_rows", "p_cols"), &GDBlas::new_mat, DEFVAL(-1), DEFVAL(-1));
	ClassDB::bind_method(D_METHOD("new_complex_mat", "p_rows", "p_cols"),
						 &GDBlas::new_complex_mat, DEFVAL(-1), DEFVAL(-1));
	ClassDB::bind_method(D_METHOD("linspace", "p_start", "p_end", "p_count"), &GDBlas::linspace);
	ClassDB::bind_method(D_METHOD("get_version"), &GDBlas::get_version);

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

	Ref< GDBlasMat > mat = GDBlasMat::new_mat(d.m, d.n, GDBlasMat::BLAS_COMPLEX_MATRIX, &error);

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

Variant GDBlas::get_version() {
	return GDBLAS_VERSION;
}
