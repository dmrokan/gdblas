#ifndef GDBLAS_GEOMETRY_H_
#define GDBLAS_GEOMETRY_H_

#include <godot_cpp/core/class_db.hpp>
#include <godot_cpp/variant/packed_vector2_array.hpp>
#include <godot_cpp/variant/rect2.hpp>
#include <godot_cpp/variant/string.hpp>

#define GDBLAS_NaN_RECT Rect2(std::nan("nan"), std::nan("nan"), std::nan("nan"), std::nan("nan"))

namespace godot {

namespace GDBlasGeometry {
typedef double scalar_t;

scalar_t area(const Array &array);
scalar_t area(const PackedVector2Array &array);

Array buffer(const PackedVector2Array &array, double buffer_distance,
		size_t points_per_join, size_t points_per_end, size_t points_per_circle);

Vector2 centroid(const Array &array);

PackedVector2Array closest_points(const Array &geom1, const Array &geom2, bool left_to_right = true);
PackedVector2Array closest_points(const Array &geom1, const Vector2 &geom2, bool left_to_right = true);
PackedVector2Array closest_points(const Array &geom1, const PackedVector2Array &geom2, bool left_to_right = true);
PackedVector2Array closest_points(const PackedVector2Array &geom1, const Vector2 &geom2, bool left_to_right = true);

PackedVector2Array convex_hull(const Array &array);
PackedVector2Array convex_hull(const PackedVector2Array &array);

Array correct(const Array &array);
PackedVector2Array correct(const PackedVector2Array &array);

int covered_by(const Array &geom1, const Array &geom2, bool left_to_right = true);
int covered_by(const Array &geom1, const PackedVector2Array &geom2, bool left_to_right = true);
int covered_by(const PackedVector2Array &geom1, const PackedVector2Array &geom2, bool left_to_right = true);
int covered_by(const Array &geom1, const Vector2 &geom2, bool left_to_right = true);
int covered_by(const PackedVector2Array &geom1, const Vector2 &geom2, bool left_to_right = true);

int crosses(const Array &geom1, const Array &geom2);
int crosses(const Array &geom1, const PackedVector2Array &geom2);
int crosses(const PackedVector2Array &geom1, const PackedVector2Array &geom2);

Array densify(const Array &geom, scalar_t max_distance);
PackedVector2Array densify(const PackedVector2Array &geom, scalar_t max_distance);

Array difference(const Array &geom1, const Array &geom2, bool left_to_right = true);
Array difference(const Array &geom1, const PackedVector2Array &geom2, bool left_to_right = true);
Array difference(const PackedVector2Array &geom1, const PackedVector2Array &geom2, bool left_to_right = true);

scalar_t discrete_frechet_distance(const PackedVector2Array &geom1, const PackedVector2Array &geom2);

scalar_t discrete_hausdorff_distance(const PackedVector2Array &geom1, const PackedVector2Array &geom2);

int disjoint(const Array &geom1, const Array &geom2);
int disjoint(const Array &geom1, const PackedVector2Array &geom2);
int disjoint(const Array &geom1, const Vector2 &geom2);
int disjoint(const PackedVector2Array &geom1, const PackedVector2Array &geom2);
int disjoint(const PackedVector2Array &geom1, const Vector2 &geom2);

scalar_t distance(const Array &geom1, const Array &geom2, bool comparable = false);
scalar_t distance(const Array &geom1, const PackedVector2Array &geom2, bool comparable = false);
scalar_t distance(const Array &geom1, const Vector2 &geom2, bool comparable = false);
scalar_t distance(const PackedVector2Array &geom1, const PackedVector2Array &geom2, bool comparable = false);
scalar_t distance(const PackedVector2Array &geom1, const Vector2 &geom2, bool comparable = false);

Rect2 envelope(const Array &geom);
Rect2 envelope(const PackedVector2Array &geom);

int equals(const Array &geom1, const Array &geom2);
int equals(const PackedVector2Array &geom1, const PackedVector2Array &geom2);

Array intersection(const Array &geom1, const Array &geom2);
Array intersection(const Array &geom1, const PackedVector2Array &geom2);
Array intersection(const PackedVector2Array &geom1, const PackedVector2Array &geom2);

int intersects(const Array &geom1, const Array &geom2);
int intersects(const Array &geom1, const PackedVector2Array &geom2);
int intersects(const PackedVector2Array &geom1, const PackedVector2Array &geom2);

int is_simple(const Array &geom);
int is_simple(const PackedVector2Array &geom);

int is_valid(const Array &geom);
int is_valid(const PackedVector2Array &geom);

scalar_t length(const Array &geom);
scalar_t length(const PackedVector2Array &geom);

int overlaps(const Array &geom1, const Array &geom2);
int overlaps(const Array &geom1, const PackedVector2Array &geom2);
int overlaps(const PackedVector2Array &geom1, const PackedVector2Array &geom2);

scalar_t perimeter(const Array &geom);
scalar_t perimeter(const PackedVector2Array &geom);

String relation(const Array &geom1, const Array &geom2, bool left_to_right = true);
String relation(const Array &geom1, const PackedVector2Array &geom2, bool left_to_right = true);
String relation(const PackedVector2Array &geom1, const PackedVector2Array &geom2, bool left_to_right = true);
String relation(const Array &geom1, const Vector2 &geom2, bool left_to_right = true);
String relation(const PackedVector2Array &geom1, const Vector2 &geom2, bool left_to_right = true);

Array reverse(const Array &geom);
PackedVector2Array reverse(const PackedVector2Array &geom);

Array simplify(const Array &geom, scalar_t max_distance);
PackedVector2Array simplify(const PackedVector2Array &geom, scalar_t max_distance);

Array sym_difference(const Array &geom1, const Array &geom2, bool left_to_right = true);
Array sym_difference(const Array &geom1, const PackedVector2Array &geom2, bool left_to_right = true);
Array sym_difference(const PackedVector2Array &geom1, const PackedVector2Array &geom2, bool left_to_right = true);

int touches(const Array &geom);
int touches(const PackedVector2Array &geom);
int touches(const Array &geom1, const Array &geom2);
int touches(const Array &geom1, const PackedVector2Array &geom2);
int touches(const PackedVector2Array &geom1, const PackedVector2Array &geom2);

Array union_(const Array &geom1, const Array &geom2);
Array union_(const Array &geom1, const PackedVector2Array &geom2);
Array union_(const PackedVector2Array &geom1, const PackedVector2Array &geom2);

Array unique(const Array &geom);
PackedVector2Array unique(const PackedVector2Array &geom);

int within(const Array &geom1, const Array &geom2, bool left_to_right = true);
int within(const Array &geom1, const PackedVector2Array &geom2, bool left_to_right = true);
int within(const Array &geom1, const Vector2 &geom2, bool left_to_right = false);
int within(const PackedVector2Array &geom1, const PackedVector2Array &geom2, bool left_to_right = true);
int within(const PackedVector2Array &geom1, const Vector2 &geom2, bool left_to_right = false);
} // namespace GDBlasGeometry

} // namespace godot

#endif // GDBLAS_GEOMETRY_H_
