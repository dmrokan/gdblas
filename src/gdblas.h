#ifndef GDBLAS_H_
#define GDBLAS_H_

#include <godot_cpp/core/class_db.hpp>
#include <godot_cpp/variant/array.hpp>

#include "gdblas_mat.h"

namespace godot {

class GDBlas : public RefCounted {
	GDCLASS(GDBlas, RefCounted)

protected:
	static void _bind_methods();

	GDBlasMat::Dimension _new_mat_dimension(Variant p_rows, int p_cols, int *error);

public:
	GDBlas();
	~GDBlas();

	Variant new_mat(Variant p_rows, int p_cols);
	Variant new_complex_mat(Variant p_rows, int p_cols);
	Variant linspace(GDBlasMat::scalar_t p_start, GDBlasMat::scalar_t p_end, int p_count);
	Variant mat_to_image_data(Array p_mat_array, int p_channel_width = 1);

#ifdef GDBLAS_WITH_GEOMETRY

	Variant geom_area(Variant p_geom);
	Variant geom_buffer(Variant p_geom, double p_buffer_distance = 2.0, int p_points_per_join = 4,
			int p_points_per_end = 4, int p_points_per_circle = 4);
	Variant geom_centroid(Variant p_geom);
	Variant geom_closest_points(Variant p_geom1, Variant p_geom2);
	Variant geom_convex_hull(Variant p_geom);
	Variant geom_correct(Variant p_geom);
	Variant geom_covered_by(Variant p_geom1, Variant p_geom2);
	Variant geom_crosses(Variant p_geom1, Variant p_geom2);
	Variant geom_densify(Variant p_geom, GDBlasMat::scalar_t p_max_distance);
	Variant geom_difference(Variant p_geom1, Variant p_geom2);
	Variant geom_discrete_frechet_distance(Variant p_geom1, Variant p_geom2);
	Variant geom_discrete_hausdorff_distance(Variant p_geom1, Variant p_geom2);
	Variant geom_disjoint(Variant p_geom1, Variant p_geom2);
	Variant geom_distance(Variant p_geom1, Variant p_geom2);
	Variant geom_comparable_distance(Variant p_geom1, Variant p_geom2);
	Variant geom_envelope(Variant p_geom);
	Variant geom_equals(Variant p_geom1, Variant p_geom2);
	Variant geom_intersection(Variant p_geom1, Variant p_geom2);
	Variant geom_intersects(Variant p_geom1, Variant p_geom2);
	Variant geom_is_simple(Variant p_geom);
	Variant geom_is_valid(Variant p_geom);
	Variant geom_length(Variant p_geom);
	Variant geom_overlaps(Variant p_geom1, Variant p_geom2);
	Variant geom_perimeter(Variant p_geom);
	Variant geom_relation(Variant p_geom1, Variant p_geom2);
	Variant geom_reverse(Variant p_geom);
	Variant geom_simplify(Variant p_geom, GDBlasMat::scalar_t p_max_distance);
	Variant geom_sym_difference(Variant p_geom1, Variant p_geom2);
	Variant geom_touches(Variant p_geom1, Variant p_geom2 = Variant());
	Variant geom_union_(Variant p_geom1, Variant p_geom2);
	Variant geom_unique(Variant p_geom);
	Variant geom_within(Variant p_geom1, Variant p_geom2);

#endif // GDBLAS_WITH_GEOMETRY

	Variant get_version();
};

} //namespace godot

#endif // GDBLAS_H_
