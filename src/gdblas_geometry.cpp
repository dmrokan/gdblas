#include "gdblas_geometry.h"

#include <vector>
#include <variant>
#include <cmath>

#define BOOST_NO_EXCEPTIONS
#include <boost/geometry.hpp>
#include <boost/polygon/polygon.hpp>
#include <boost/foreach.hpp>

namespace bg = boost::geometry;

using namespace godot;

using point_t = bg::model::d2::point_xy<GDBlasGeometry::scalar_t>;
using polygon_t = bg::model::polygon<point_t>;
using polygon_set_t = std::vector<point_t>;
using ring_t = bg::model::ring<point_t>;
using mpolygon_t = bg::model::multi_polygon<polygon_t>;
using line_t = boost::geometry::model::linestring<point_t>;
using segment_t = boost::geometry::model::segment<point_t>;
using box_t = boost::geometry::model::box<point_t>;

namespace godot {

namespace GDBlasGeometry {

static polygon_set_t vector2_array_to_polygon_set(const PackedVector2Array &array) {
	polygon_set_t set;
	set.resize(array.size());

	size_t i = 0;
	for (auto it = set.begin(); it != set.end(); ++it, ++i) {
		const Vector2 &entry = array[i];
		(*it).set<0>(entry.x);
		(*it).set<1>(entry.y);
	}

	return set;
}

inline static bool is_vector2_array_ring(const PackedVector2Array &array) {
	size_t size = array.size();

	if (size < 2)
		return false;

	return array[0] == array[size-1];
}

template <typename ARG>
static int vector2_array_to_geom(const PackedVector2Array &array, ARG &geom) {
	if (array.is_empty())
		return -1;

	geom.resize(array.size());

	size_t i = 0;
	for (auto it = boost::begin(geom); it != boost::end(geom); ++it, ++i) {
		const Vector2 &entry = array[i];
		bg::set<0>(*it, entry.x);
		bg::set<1>(*it, entry.y);
	}

	return i;
}

template <typename ARG>
static int geom_to_vector2_array(const ARG &geom, PackedVector2Array &array) {
	if (geom.size() < 1)
		return -1;

	array.resize(geom.size());

	int i = 0;
	for (auto it = boost::begin(geom); it != boost::end(geom); ++it, ++i) {
		const Vector2 entry(bg::get<0>(*it), bg::get<1>(*it));

		array[i] = entry;
	}

	return array.size();
}

static int polygon_to_array(const polygon_t &polygon, Array &array) {
	auto &outer = polygon.outer();
	if (outer.size() < 1)
		return -1;

	PackedVector2Array outer_array;
	int size = geom_to_vector2_array(outer, outer_array);

	array.append(outer_array);

	auto &inners = polygon.inners();
	for (auto it = boost::begin(inners); it != boost::end(inners); ++it) {
		if ((*it).size() < 1)
			continue;

		PackedVector2Array inners_i;
		size = geom_to_vector2_array(*it, inners_i);
		array.append(inners_i);
	}

	return array.size();
}

static int array_to_polygon(const Array &array, polygon_t &polygon) {
	const PackedVector2Array ai = array[0];
	if (ai.size() < 1)
		return -1;

	int result = vector2_array_to_geom(ai, polygon.outer());
	if (result < 0)
		return result;

	int inner_size = array.size() - 1;
	if (inner_size < 1) {
		return 1;
	}

	polygon.inners().resize(inner_size);

	int i;
	for (i = 1; i < array.size(); ++i) {
		const PackedVector2Array ai2 = array[i];
		if (ai2.size() < 1)
			continue;

		result = vector2_array_to_geom(ai2, polygon.inners()[i-1]);
		if (result < 0)
			return result;
	}

	return i;
}

inline static int convert_from_godot_type_to_boost_type(const Array &array, polygon_t &polygon) {
	return array_to_polygon(array, polygon);
}

inline static int convert_from_godot_type_to_boost_type(const PackedVector2Array &array, line_t &line) {
	return vector2_array_to_geom(array, line);
}

inline static int convert_from_godot_type_to_boost_type(const PackedVector2Array &array, ring_t &ring) {
	return vector2_array_to_geom(array, ring);
}

inline static int convert_from_godot_type_to_boost_type(const Vector2 &vec, point_t &point) {
	point.set<0>(vec.x);
	point.set<1>(vec.y);

	return 0;
}

inline static int convert_from_boost_type_to_godot_type(const ring_t &ring, PackedVector2Array &vec) {
	return geom_to_vector2_array(ring, vec);
}

inline static int convert_from_boost_type_to_godot_type(const line_t &line, PackedVector2Array &vec) {
	return geom_to_vector2_array(line, vec);
}

inline static int convert_from_boost_type_to_godot_type(const polygon_t &polygon, Array &array) {
	return polygon_to_array(polygon, array);
}

inline static int convert_from_boost_type_to_godot_type(const point_t &point, Vector2 &vec) {
	vec.x = point.get<0>();
	vec.y = point.get<1>();

	return 0;
}

inline static int convert_from_boost_type_to_godot_type(const box_t &box, Rect2 &rect) {
	const point_t &min_corner = box.min_corner();
	const point_t &max_corner = box.max_corner();

	const Vector2 pos(bg::get<0>(min_corner), bg::get<1>(min_corner));
	const Vector2 end(bg::get<0>(max_corner), bg::get<1>(max_corner));

	rect.set_position(pos);
	rect.set_size(end - pos);

	return 0;
}

template <typename T, typename ARG>
static scalar_t area(const ARG &arg) {
	T geom{};

	int result = convert_from_godot_type_to_boost_type(arg, geom);
	if (result < 0)
		return std::nan("nan");

	scalar_t area = bg::area(geom);

	return area;
}

template <typename T1, typename T2, typename ARG1, typename ARG2>
static PackedVector2Array closest_points(const ARG1 &arg1, const ARG2 &arg2, bool left_to_right) {
	T1 geom1{};

	int result = convert_from_godot_type_to_boost_type(arg1, geom1);
	if (result < 0)
		return Array();

	T2 geom2{};
	result = convert_from_godot_type_to_boost_type(arg2, geom2);
	if (result < 0)
		return Array();

	segment_t seg{};

	if (left_to_right)
		bg::closest_points(geom1, geom2, seg);
	else
		bg::closest_points(geom2, geom1, seg);

	Array output;
	output.append(Vector2(bg::get<0, 0>(seg), bg::get<0, 1>(seg)));
	output.append(Vector2(bg::get<1, 0>(seg), bg::get<1, 1>(seg)));

	return output;
}

template <typename T, typename ARG>
static PackedVector2Array convex_hull(const ARG &arg) {
	T geom{};

	int result = convert_from_godot_type_to_boost_type(arg, geom);
	if (result < 0)
		return PackedVector2Array();

	ring_t hull{};
	bg::convex_hull(geom, hull);

	PackedVector2Array output;
	result = convert_from_boost_type_to_godot_type(hull, output);
	if (result < 0)
		return PackedVector2Array();

	return output;
}

template <typename T, typename ARG>
static ARG correct(const ARG &arg) {
	T geom{};

	int result = convert_from_godot_type_to_boost_type(arg, geom);
	if (result < 0)
		return ARG();

	bg::correct(geom);

	ARG output;
	result = convert_from_boost_type_to_godot_type(geom, output);
	if (result < 0)
		return ARG();

	return output;
}

template <typename T1, typename T2, bool B = true, typename ARG1, typename ARG2>
static int covered_by(const ARG1 &geom1, const ARG2 &geom2, bool left_to_right) {
	T1 g1{};
	int result = convert_from_godot_type_to_boost_type(geom1, g1);
	if (result < 0)
		return result;

	T2 g2{};
	result = convert_from_godot_type_to_boost_type(geom2, g2);
	if (result < 0)
		return result;

	if constexpr (B) {
		result = static_cast<int>(left_to_right ? bg::covered_by(g1, g2) : bg::covered_by(g2, g1));
	} else {
		result = static_cast<int>(left_to_right ? bg::covered_by(g1, g2) : -1);
	}

	return result;
}

template <typename T1, typename T2, typename ARG1, typename ARG2>
static int crosses(const ARG1 &arg1, const ARG2 &arg2) {
	T1 geom1{};
	int result = convert_from_godot_type_to_boost_type(arg1, geom1);
	if (result < 0)
		return result;

	T2 geom2{};
	result = convert_from_godot_type_to_boost_type(arg2, geom2);
	if (result < 0)
		return result;

	result = static_cast<int>(bg::crosses(geom1, geom2));

	return result;
}

template <typename T, typename ARG>
static ARG densify(const ARG &arg, scalar_t max_distance) {
	T geom{};
	int result = convert_from_godot_type_to_boost_type(arg, geom);
	if (result < 0)
		return ARG();

	T output_geom{};
	bg::densify(geom, output_geom, max_distance);

	ARG output;
	result = convert_from_boost_type_to_godot_type(output_geom, output);
	if (result < 0)
		return ARG();

	return output;
}

template <typename T1, typename T2, bool B = true, typename ARG1, typename ARG2>
static Array difference(const ARG1 &arg1, const ARG2 &arg2, bool left_to_right) {
	T1 geom1{};
	int result = convert_from_godot_type_to_boost_type(arg1, geom1);
	if (result < 0)
		return Array();

	T2 geom2{};
	result = convert_from_godot_type_to_boost_type(arg2, geom2);
	if (result < 0)
		return Array();

	std::list<T1> collection;

	if constexpr (B) {
		left_to_right ? bg::difference(geom1, geom2, collection) : bg::difference(geom2, geom1, collection);
	} else {
		left_to_right ? bg::difference(geom1, geom2, collection) : (void) 0;
	}

	Array output;

	BOOST_FOREACH(T1 const &g, collection) {
		ARG1 out_i;
		result = convert_from_boost_type_to_godot_type(g, out_i);
		if (result < 0)
			continue;

		output.append(out_i);
    }

	return output;
}
template <typename T1, typename T2, typename ARG1, typename ARG2>
static scalar_t discrete_frechet_distance(const ARG1 &arg1, const ARG2 &arg2) {
	T1 geom1{};
	int result = convert_from_godot_type_to_boost_type(arg1, geom1);
	if (result < 0)
		return std::nan("nan");

	T2 geom2{};
	result = convert_from_godot_type_to_boost_type(arg2, geom2);
	if (result < 0)
		return std::nan("nan");

	return bg::discrete_frechet_distance(geom1, geom2);
}

template <typename T1, typename T2, typename ARG1, typename ARG2>
static scalar_t discrete_hausdorff_distance(const ARG1 &arg1, const ARG2 &arg2) {
	T1 geom1{};
	int result = convert_from_godot_type_to_boost_type(arg1, geom1);
	if (result < 0)
		return std::nan("nan");

	T2 geom2{};
	result = convert_from_godot_type_to_boost_type(arg2, geom2);
	if (result < 0)
		return std::nan("nan");

	return bg::discrete_hausdorff_distance(geom1, geom2);
}

template <typename T1, typename T2, typename ARG1, typename ARG2>
static int disjoint(const ARG1 &arg1, const ARG2 &arg2) {
	T1 geom1{};
	int result = convert_from_godot_type_to_boost_type(arg1, geom1);
	if (result < 0)
		return result;

	T2 geom2{};
	result = convert_from_godot_type_to_boost_type(arg2, geom2);
	if (result < 0)
		return result;

	result = static_cast<int>(bg::disjoint(geom1, geom2));

	return result;
}

template <typename T1, typename T2, typename ARG1, typename ARG2>
static scalar_t distance(const ARG1 &arg1, const ARG2 &arg2, bool comparable) {
	T1 geom1{};
	int result = convert_from_godot_type_to_boost_type(arg1, geom1);
	if (result < 0)
		return result;

	T2 geom2{};
	result = convert_from_godot_type_to_boost_type(arg2, geom2);
	if (result < 0)
		return result;

	if (comparable)
		return bg::comparable_distance(geom1, geom2);

	return bg::distance(geom1, geom2);
}

template <typename T, typename ARG>
static Rect2 envelope(const ARG &arg) {
	T geom{};
	int result = convert_from_godot_type_to_boost_type(arg, geom);
	if (result < 0)
		return GDBLAS_NaN_RECT;

	box_t output{};
	bg::envelope(geom, output);

	Rect2 rect;
	result = convert_from_boost_type_to_godot_type(output, rect);
	if (result < 0)
		return GDBLAS_NaN_RECT;

	return rect;
}

template <typename T, typename ARG>
static int equals(const ARG &arg1, const ARG &arg2) {
	T geom1{};
	int result = convert_from_godot_type_to_boost_type(arg1, geom1);
	if (result < 0)
		return result;

	T geom2{};
	result = convert_from_godot_type_to_boost_type(arg2, geom2);
	if (result < 0)
		return result;

	result = static_cast<int>(bg::equals(geom1, geom2));

	return result;
}

template <typename T1, typename T2, typename ARG1, typename ARG2>
static Array intersection(const ARG1 &arg1, const ARG2 &arg2) {
	T1 geom1{};
	int result = convert_from_godot_type_to_boost_type(arg1, geom1);
	if (result < 0)
		return Array();

	T2 geom2{};
	result = convert_from_godot_type_to_boost_type(arg2, geom2);
	if (result < 0)
		return Array();

	std::deque<T2> collection;
	bg::intersection(geom1, geom2, collection);

	Array output;

	BOOST_FOREACH(T2 const &g, collection) {
		ARG2 out_i;
		result = convert_from_boost_type_to_godot_type(g, out_i);
		if (result < 0)
			continue;

		output.append(out_i);
    }

	return output;
}

template <typename T, typename ARG>
static int intersects(const ARG &arg) {
	T geom{};
	int result = convert_from_godot_type_to_boost_type(arg, geom);
	if (result < 0)
		return result;

	result = static_cast<int>(bg::intersects(geom));

	return result;
}

template <typename T1, typename T2, typename ARG1, typename ARG2>
static int intersects(const ARG1 &arg1, const ARG2 &arg2) {
	T1 geom1{};
	int result = convert_from_godot_type_to_boost_type(arg1, geom1);
	if (result < 0)
		return result;

	T2 geom2{};
	result = convert_from_godot_type_to_boost_type(arg2, geom2);
	if (result < 0)
		return result;

	result = static_cast<int>(bg::intersects(geom1, geom2));

	return result;
}

template <typename T, typename ARG>
static int is_simple(const ARG &arg) {
	T geom{};
	int result = convert_from_godot_type_to_boost_type(arg, geom);
	if (result < 0)
		return result;

	result = static_cast<int>(bg::is_simple(geom));

	return result;
}
template <typename T, typename ARG>
static int is_valid(const ARG &arg) {
	T geom{};
	int result = convert_from_godot_type_to_boost_type(arg, geom);
	if (result < 0)
		return result;

	bg::validity_failure_type failure;
	bg::is_valid(geom, failure);

	return static_cast<int>(failure);
}

template <typename T, typename ARG>
static scalar_t length(const ARG &arg) {
	T geom{};
	scalar_t result = convert_from_godot_type_to_boost_type(arg, geom);
	if (result < 0)
		return std::nan("nan");

	return bg::length(geom);
}

template <typename T1, typename T2, typename ARG1, typename ARG2>
static int overlaps(const ARG1 &arg1, const ARG2 &arg2) {
	T1 geom1{};
	int result = convert_from_godot_type_to_boost_type(arg1, geom1);
	if (result < 0)
		return result;

	T2 geom2{};
	result = convert_from_godot_type_to_boost_type(arg2, geom2);
	if (result < 0)
		return result;

	result = static_cast<int>(bg::overlaps(geom1, geom2));

	return result;
}

template <typename T, typename ARG>
static scalar_t perimeter(const ARG &arg) {
	T geom{};
	scalar_t result = convert_from_godot_type_to_boost_type(arg, geom);
	if (result < 0)
		return std::nan("nan");

	return bg::perimeter(geom);
}

template <typename T1, typename T2, typename ARG1, typename ARG2>
static String relation(const ARG1 &arg1, const ARG2 &arg2, bool left_to_right) {
	T1 geom1{};
	int result = convert_from_godot_type_to_boost_type(arg1, geom1);
	if (result < 0)
		return String();

	T2 geom2{};
	result = convert_from_godot_type_to_boost_type(arg2, geom2);
	if (result < 0)
		return String();

	bg::de9im::matrix matrix = left_to_right ? bg::relation(geom1, geom2) : bg::relation(geom2, geom1);

	String output(matrix.str().c_str());

	return output;
}

template <typename T, typename ARG>
static ARG reverse(const ARG &arg) {
	T geom{};
	scalar_t result = convert_from_godot_type_to_boost_type(arg, geom);
	if (result < 0)
		return ARG();

	bg::reverse(geom);

	ARG output;
	result = convert_from_boost_type_to_godot_type(geom, output);
	if (result < 0)
		return ARG();

	return output;
}

template <typename T, typename ARG>
static ARG simplify(const ARG &arg, scalar_t max_distance) {
	T geom{};
	int result = convert_from_godot_type_to_boost_type(arg, geom);
	if (result < 0)
		return ARG();

	T output_geom{};
	bg::simplify(geom, output_geom, max_distance);

	ARG output;
	result = convert_from_boost_type_to_godot_type(output_geom, output);
	if (result < 0)
		return ARG();

	return output;
}

template <typename T1, typename T2, typename ARG1, typename ARG2>
static Array sym_difference(const ARG1 &arg1, const ARG2 &arg2, bool left_to_right) {
	T1 geom1{};
	int result = convert_from_godot_type_to_boost_type(arg1, geom1);
	if (result < 0)
		return Array();

	T2 geom2{};
	result = convert_from_godot_type_to_boost_type(arg2, geom2);
	if (result < 0)
		return Array();

	std::list<T1> collection;
	left_to_right ? bg::sym_difference(geom1, geom2, collection) :
		bg::sym_difference(geom2, geom1, collection);

	Array output;

	BOOST_FOREACH(T1 const &g, collection) {
		ARG1 out_i;
		result = convert_from_boost_type_to_godot_type(g, out_i);
		if (result < 0)
			continue;

		output.append(out_i);
    }

	return output;
}

template <typename T, typename ARG>
static int touches(const ARG &arg) {
	T geom{};
	scalar_t result = convert_from_godot_type_to_boost_type(arg, geom);
	if (result < 0)
		return result;

	result = static_cast<int>(bg::touches(geom));

	return result;
}

template <typename T1, typename T2, typename ARG1, typename ARG2>
static int touches(const ARG1 &arg1, const ARG2 &arg2) {
	T1 geom1{};
	int result = convert_from_godot_type_to_boost_type(arg1, geom1);
	if (result < 0)
		return result;

	T2 geom2{};
	result = convert_from_godot_type_to_boost_type(arg2, geom2);
	if (result < 0)
		return result;

	result = static_cast<int>(bg::touches(geom1, geom2));

	return result;
}

template <typename T1, typename T2, typename ARG1, typename ARG2>
static Array union_(const ARG1 &arg1, const ARG2 &arg2) {
	T1 geom1{};
	int result = convert_from_godot_type_to_boost_type(arg1, geom1);
	if (result < 0)
		return Array();

	T2 geom2{};
	result = convert_from_godot_type_to_boost_type(arg2, geom2);
	if (result < 0)
		return Array();

	std::list<T1> collection;
	bg::union_(geom1, geom2, collection);

	Array output;

	BOOST_FOREACH(T1 const &g, collection) {
		ARG1 out_i;
		result = convert_from_boost_type_to_godot_type(g, out_i);
		if (result < 0)
			continue;

		output.append(out_i);
    }

	return output;
}

template <typename T, typename ARG>
static ARG unique(const ARG &arg) {
	T geom{};
	scalar_t result = convert_from_godot_type_to_boost_type(arg, geom);
	if (result < 0)
		return ARG();

	bg::unique(geom);

	ARG output;
	result = convert_from_boost_type_to_godot_type(geom, output);
	if (result < 0)
		return ARG();

	return output;
}

template <typename T1, typename T2, bool B = true, typename ARG1, typename ARG2>
static int within(const ARG1 &arg1, const ARG2 &arg2, bool left_to_right) {
	T1 geom1{};
	int result = convert_from_godot_type_to_boost_type(arg1, geom1);
	if (result < 0)
		return result;

	T2 geom2{};
	result = convert_from_godot_type_to_boost_type(arg2, geom2);
	if (result < 0)
		return result;

	if constexpr (B) {
		result = left_to_right ? static_cast<int>(bg::within(geom1, geom2)) :
			static_cast<int>(bg::within(geom2, geom1));
	} else {
		result = left_to_right ? static_cast<int>(bg::within(geom1, geom2)) : -1;
	}

	return result;
}

} // namespace GDBlasGeometry

} // godot

GDBlasGeometry::scalar_t GDBlasGeometry::area(const Array &array) {
	return area<polygon_t>(array);
}

GDBlasGeometry::scalar_t GDBlasGeometry::area(const PackedVector2Array &array) {
	if (is_vector2_array_ring(array))
		return area<ring_t>(array);

	return std::nan("nan");
}

Array GDBlasGeometry::buffer(const PackedVector2Array &array, double buffer_distance,
							 size_t points_per_join, size_t points_per_end, size_t points_per_circle) {
	using end_flat_t = bg::strategy::buffer::end_flat;
	using end_round_t = bg::strategy::buffer::end_round;
	using join_miter_t = bg::strategy::buffer::join_miter;
	using join_round_t = bg::strategy::buffer::join_round;

	int strategy_flags = 0;

    bg::strategy::buffer::distance_symmetric<scalar_t> distance_strategy{ buffer_distance };
    bg::strategy::buffer::point_circle circle_strategy{ points_per_circle };
    bg::strategy::buffer::side_straight side_strategy{};

	std::variant<end_round_t, end_flat_t> end_strategy;
	if (points_per_end > 0) {
		end_strategy.emplace<0>(points_per_end);
		strategy_flags |= 1;
	} else {
		end_strategy.emplace<1>();
	}

	std::variant<join_round_t, join_miter_t> join_strategy;
	if (points_per_join > 0) {
		join_strategy.emplace<0>(points_per_join);
		strategy_flags |= 2;
	} else {
		join_strategy.emplace<1>();
	}

    mpolygon_t mpoly;

	const polygon_set_t vec = vector2_array_to_polygon_set(array);
    line_t ls(vec.begin(), vec.end());

	switch (strategy_flags) {
	case 1:
		bg::buffer(ls, mpoly, distance_strategy, side_strategy, std::get<join_miter_t>(join_strategy),
				   std::get<end_round_t>(end_strategy), circle_strategy);
		break;
	case 2:
		bg::buffer(ls, mpoly, distance_strategy, side_strategy, std::get<join_round_t>(join_strategy),
				   std::get<end_flat_t>(end_strategy), circle_strategy);
		break;
	case 3:
		bg::buffer(ls, mpoly, distance_strategy, side_strategy, std::get<join_round_t>(join_strategy),
				   std::get<end_round_t>(end_strategy), circle_strategy);
		break;
	default:
		bg::buffer(ls, mpoly, distance_strategy, side_strategy, std::get<join_miter_t>(join_strategy),
				   std::get<end_flat_t>(end_strategy), circle_strategy);
		break;
	}

	Array result;

	for (auto it = boost::begin(mpoly); it != boost::end(mpoly); ++it) {
		polygon_to_array(*it, result);
	}

	return result;
}

Vector2 GDBlasGeometry::centroid(const Array &array) {
	polygon_t polygon;
	array_to_polygon(array, polygon);

	point_t c{};
	bg::centroid(polygon, c);

	return Vector2(bg::get<0>(c), bg::get<1>(c));
}

PackedVector2Array GDBlasGeometry::closest_points(const Array &geom1, const Array &geom2, bool left_to_right) {
	return closest_points<polygon_t, polygon_t>(geom1, geom2, left_to_right);
}

PackedVector2Array GDBlasGeometry::closest_points(const Array &geom1, const Vector2 &geom2, bool left_to_right) {
	return closest_points<polygon_t, point_t>(geom1, geom2, left_to_right);
}

PackedVector2Array GDBlasGeometry::closest_points(const Array &geom1, const PackedVector2Array &geom2, bool left_to_right) {
	if (is_vector2_array_ring(geom2))
		return closest_points<polygon_t, ring_t>(geom1, geom2, left_to_right);

	return closest_points<polygon_t, line_t>(geom1, geom2, left_to_right);
}

PackedVector2Array GDBlasGeometry::closest_points(const PackedVector2Array &geom1, const Vector2 &geom2, bool left_to_right) {
	if (is_vector2_array_ring(geom1))
		return closest_points<ring_t, point_t>(geom1, geom2, left_to_right);

	return closest_points<line_t, point_t>(geom1, geom2, left_to_right);
}

PackedVector2Array GDBlasGeometry::convex_hull(const Array &array) {
	return convex_hull<polygon_t>(array);
}

PackedVector2Array GDBlasGeometry::convex_hull(const PackedVector2Array &array) {
	if (is_vector2_array_ring(array))
		return convex_hull<ring_t>(array);

	return convex_hull<line_t>(array);
}

Array GDBlasGeometry::correct(const Array &array) {
	return correct<polygon_t>(array);
}

PackedVector2Array GDBlasGeometry::correct(const PackedVector2Array &array) {
	if (is_vector2_array_ring(array))
		return correct<ring_t>(array);

	return correct<line_t>(array);
}

int GDBlasGeometry::covered_by(const Array &geom1, const Array &geom2, bool left_to_right) {
	return covered_by<polygon_t, polygon_t>(geom1, geom2, left_to_right);
}

int GDBlasGeometry::covered_by(const Array &geom1, const PackedVector2Array &geom2, bool left_to_right) {
	if (is_vector2_array_ring(geom2))
		return covered_by<polygon_t, ring_t, false>(geom1, geom2, left_to_right);

	return -1;
}

int GDBlasGeometry::covered_by(const PackedVector2Array &geom1, const PackedVector2Array &geom2, bool left_to_right) {
	int cond = static_cast<int>(is_vector2_array_ring(geom1));
	cond |= is_vector2_array_ring(geom2) ? 2 : 0;

	switch (cond) {
	case 0:
		return covered_by<line_t, line_t, true>(geom1, geom2, left_to_right);
	case 2:
		return covered_by<line_t, ring_t, false>(geom1, geom2, left_to_right);
	case 3:
		return covered_by<ring_t, ring_t, true>(geom1, geom2, left_to_right);
	default:
		break;
	}

	return -1;
}

int GDBlasGeometry::covered_by(const Array &geom1, const Vector2 &geom2, bool left_to_right) {
	return covered_by<point_t, polygon_t, false>(geom2, geom1, !left_to_right);
}

int GDBlasGeometry::covered_by(const PackedVector2Array &geom1, const Vector2 &geom2, bool left_to_right) {
	if (is_vector2_array_ring(geom1))
		return covered_by<point_t, ring_t, false>(geom2, geom1, !left_to_right);

	return covered_by<point_t, line_t, false>(geom2, geom1, !left_to_right);
}

int GDBlasGeometry::crosses(const Array &geom1, const Array &geom2) {
	return crosses<polygon_t, polygon_t>(geom1, geom2);
}

int GDBlasGeometry::crosses(const Array &geom1, const PackedVector2Array &geom2) {
	if (is_vector2_array_ring(geom2))
		return crosses<polygon_t, ring_t>(geom1, geom2);

	return crosses<polygon_t, line_t>(geom1, geom2);
}

int GDBlasGeometry::crosses(const PackedVector2Array &geom1, const PackedVector2Array &geom2) {
	int cond = static_cast<int>(is_vector2_array_ring(geom1));
	cond |= is_vector2_array_ring(geom2) ? 2 : 0;

	switch (cond) {
	case 0:
		return crosses<line_t, line_t>(geom1, geom2);
	case 1:
		return crosses<ring_t, line_t>(geom1, geom2);
	case 2:
		return crosses<line_t, ring_t>(geom1, geom2);
	case 3:
		return crosses<ring_t, ring_t>(geom1, geom2);
	default:
		break;
	}

	return -1;
}

Array GDBlasGeometry::densify(const Array &geom, scalar_t max_distance) {
	return densify<polygon_t>(geom, max_distance);
}

PackedVector2Array GDBlasGeometry::densify(const PackedVector2Array &geom, scalar_t max_distance) {
	if (is_vector2_array_ring(geom))
		return densify<ring_t>(geom, max_distance);

	return densify<line_t>(geom, max_distance);
}

Array GDBlasGeometry::difference(const Array &geom1, const Array &geom2, bool left_to_right) {
	return difference<polygon_t, polygon_t>(geom1, geom2, left_to_right);
}

Array GDBlasGeometry::difference(const Array &geom1, const PackedVector2Array &geom2, bool left_to_right) {
	if (is_vector2_array_ring(geom2))
		return difference<polygon_t, ring_t>(geom1, geom2, left_to_right);

	return difference<polygon_t, line_t, false>(geom1, geom2, left_to_right);
}

Array GDBlasGeometry::difference(const PackedVector2Array &geom1, const PackedVector2Array &geom2, bool left_to_right) {
	int cond = static_cast<int>(is_vector2_array_ring(geom1));
	cond |= is_vector2_array_ring(geom2) ? 2 : 0;

	switch (cond) {
	case 0:
		return difference<line_t, line_t, true>(geom1, geom2, left_to_right);
	case 2:
		return difference<line_t, ring_t, false>(geom1, geom2, left_to_right);
	case 3:
		return difference<ring_t, ring_t, true>(geom1, geom2, left_to_right);
	default:
		break;
	}

	return Array();
}

GDBlasGeometry::scalar_t GDBlasGeometry::discrete_frechet_distance(const PackedVector2Array &geom1, const PackedVector2Array &geom2) {
	if (is_vector2_array_ring(geom1) || is_vector2_array_ring(geom2))
		return std::nan("nan");

	return discrete_frechet_distance<line_t, line_t>(geom1, geom2);
}

GDBlasGeometry::scalar_t GDBlasGeometry::discrete_hausdorff_distance(const PackedVector2Array &geom1, const PackedVector2Array &geom2) {
	return discrete_hausdorff_distance<line_t, line_t>(geom1, geom2);
}

int GDBlasGeometry::disjoint(const Array &geom1, const Array &geom2) {
	return disjoint<polygon_t, polygon_t>(geom1, geom2);
}

int GDBlasGeometry::disjoint(const Array &geom1, const PackedVector2Array &geom2) {
	if (is_vector2_array_ring(geom2))
		return disjoint<polygon_t, ring_t>(geom1, geom2);

	return disjoint<polygon_t, line_t>(geom1, geom2);
}

int GDBlasGeometry::disjoint(const Array &geom1, const Vector2 &geom2) {
	return disjoint<polygon_t, point_t>(geom1, geom2);
}

int GDBlasGeometry::disjoint(const PackedVector2Array &geom1, const PackedVector2Array &geom2) {
	int cond = static_cast<int>(is_vector2_array_ring(geom1));
	cond |= is_vector2_array_ring(geom2) ? 2 : 0;

	switch (cond) {
	case 0:
		return disjoint<line_t, line_t>(geom1, geom2);
	case 1:
		return disjoint<ring_t, line_t>(geom1, geom2);
	case 2:
		return disjoint<line_t, ring_t>(geom1, geom2);
	case 3:
		return disjoint<ring_t, ring_t>(geom1, geom2);
	default:
		break;
	}

	return std::nan("nan");
}

int GDBlasGeometry::disjoint(const PackedVector2Array &geom1, const Vector2 &geom2) {
	if (is_vector2_array_ring(geom1))
		return disjoint<ring_t, point_t>(geom1, geom2);

	return disjoint<line_t, point_t>(geom1, geom2);
}

GDBlasGeometry::scalar_t GDBlasGeometry::distance(const Array &geom1, const Array &geom2, bool comparable) {
	return distance<polygon_t, polygon_t>(geom1, geom2, comparable);
}

GDBlasGeometry::scalar_t GDBlasGeometry::distance(const Array &geom1, const PackedVector2Array &geom2, bool comparable) {
	if (is_vector2_array_ring(geom2))
		return distance<polygon_t, ring_t>(geom1, geom2, comparable);

	return distance<polygon_t, line_t>(geom1, geom2, comparable);
}

GDBlasGeometry::scalar_t GDBlasGeometry::distance(const Array &geom1, const Vector2 &geom2, bool comparable) {
	return distance<polygon_t, point_t>(geom1, geom2, comparable);
}

GDBlasGeometry::scalar_t GDBlasGeometry::distance(const PackedVector2Array &geom1, const PackedVector2Array &geom2, bool comparable) {
	int cond = static_cast<int>(is_vector2_array_ring(geom1));
	cond |= is_vector2_array_ring(geom2) ? 2 : 0;

	switch (cond) {
	case 0:
		return distance<line_t, line_t>(geom1, geom2, comparable);
	case 1:
		return distance<ring_t, line_t>(geom1, geom2, comparable);
	case 2:
		return distance<line_t, ring_t>(geom1, geom2, comparable);
	case 3:
		return distance<ring_t, ring_t>(geom1, geom2, comparable);
	default:
		break;
	}

	return std::nan("nan");
}

GDBlasGeometry::scalar_t GDBlasGeometry::distance(const PackedVector2Array &geom1, const Vector2 &geom2, bool comparable) {
	if (is_vector2_array_ring(geom1))
		return distance<ring_t, point_t>(geom1, geom2, comparable);

	return distance<line_t, point_t>(geom1, geom2, comparable);
}

Rect2 GDBlasGeometry::envelope(const Array &geom) {
	return envelope<polygon_t>(geom);
}

Rect2 GDBlasGeometry::envelope(const PackedVector2Array &geom) {
	if (is_vector2_array_ring(geom))
		return envelope<ring_t>(geom);

	return envelope<line_t>(geom);
}

int GDBlasGeometry::equals(const Array &geom1, const Array &geom2) {
	return equals<polygon_t>(geom1, geom2);
}

int GDBlasGeometry::equals(const PackedVector2Array &geom1, const PackedVector2Array &geom2) {
	int cond = is_vector2_array_ring(geom1);
	cond |= is_vector2_array_ring(geom2) ? 2 : 0;

	switch (cond) {
	case 0:
		return equals<line_t>(geom1, geom2);
	case 3:
		return equals<ring_t>(geom1, geom2);
	default:
		break;
	}

	return -1;
}

Array GDBlasGeometry::intersection(const Array &geom1, const Array &geom2) {
	return intersection<polygon_t, polygon_t>(geom1, geom2);
}

Array GDBlasGeometry::intersection(const Array &geom1, const PackedVector2Array &geom2) {
	if (is_vector2_array_ring(geom2))
		return intersection<polygon_t, ring_t>(geom1, geom2);

	return intersection<polygon_t, line_t>(geom1, geom2);
}

Array GDBlasGeometry::intersection(const PackedVector2Array &geom1, const PackedVector2Array &geom2) {
	int cond = static_cast<int>(is_vector2_array_ring(geom1));
	cond |= is_vector2_array_ring(geom2) ? 2 : 0;

	switch (cond) {
	case 0:
		return intersection<line_t, line_t>(geom1, geom2);
	case 3:
		return intersection<ring_t, ring_t>(geom1, geom2);
	default:
		break;
	}

	return Array();
}

int GDBlasGeometry::intersects(const Array &geom1, const Array &geom2) {
	if (&geom1 == &geom2)
		return intersects<polygon_t>(geom1);

	return intersects<polygon_t, polygon_t>(geom1, geom2);
}

int GDBlasGeometry::intersects(const Array &geom1, const PackedVector2Array &geom2) {
	if (is_vector2_array_ring(geom2))
		return intersects<polygon_t, ring_t>(geom1, geom2);

	return intersects<polygon_t, line_t>(geom1, geom2);
}

int GDBlasGeometry::intersects(const PackedVector2Array &geom1, const PackedVector2Array &geom2) {
	int cond = static_cast<int>(is_vector2_array_ring(geom1));

	if (&geom1 == &geom2) {
		if (cond)
			return intersects<ring_t>(geom1);
		else
			return intersects<line_t>(geom1);
	}

	cond |= is_vector2_array_ring(geom2) ? 2 : 0;

	switch (cond) {
	case 0:
		return intersects<line_t, line_t>(geom1, geom2);
	case 1:
		return intersects<ring_t, line_t>(geom1, geom2);
	case 2:
		return intersects<line_t, ring_t>(geom1, geom2);
	case 3:
		return intersects<ring_t, ring_t>(geom1, geom2);
	default:
		break;
	}

	return -1;
}

int GDBlasGeometry::is_simple(const Array &geom) {
	return is_simple<polygon_t>(geom);
}

int GDBlasGeometry::is_simple(const PackedVector2Array &geom) {
	if (is_vector2_array_ring(geom))
		return is_simple<ring_t>(geom);

	return is_simple<line_t>(geom);
}

int GDBlasGeometry::is_valid(const Array &geom) {
	return is_valid<polygon_t>(geom);
}

int GDBlasGeometry::is_valid(const PackedVector2Array &geom) {
	if (is_vector2_array_ring(geom))
		return is_valid<ring_t>(geom);

	return is_valid<line_t>(geom);
}

GDBlasGeometry::scalar_t GDBlasGeometry::length(const Array &geom) {
	return length<polygon_t>(geom);
}

GDBlasGeometry::scalar_t GDBlasGeometry::length(const PackedVector2Array &geom) {
	if (is_vector2_array_ring(geom))
		return length<ring_t>(geom);

	return length<line_t>(geom);
}

int GDBlasGeometry::overlaps(const Array &geom1, const Array &geom2) {
	return overlaps<polygon_t, polygon_t>(geom1, geom2);
}

int GDBlasGeometry::overlaps(const Array &geom1, const PackedVector2Array &geom2) {
	return overlaps<polygon_t, ring_t>(geom1, geom2);
}

int GDBlasGeometry::overlaps(const PackedVector2Array &geom1, const PackedVector2Array &geom2) {
	if (is_vector2_array_ring(geom1) || is_vector2_array_ring(geom2))
		return -1;

	return overlaps<line_t, line_t>(geom1, geom2);
}

GDBlasGeometry::scalar_t GDBlasGeometry::perimeter(const Array &geom) {
	return perimeter<polygon_t>(geom);
}

GDBlasGeometry::scalar_t GDBlasGeometry::perimeter(const PackedVector2Array &geom) {
	if (is_vector2_array_ring(geom))
		return perimeter<ring_t>(geom);

	return perimeter<line_t>(geom);
}

String GDBlasGeometry::relation(const Array &geom1, const Array &geom2, bool left_to_right) {
	return relation<polygon_t, polygon_t>(geom1, geom2, left_to_right);
}

String GDBlasGeometry::relation(const Array &geom1, const PackedVector2Array &geom2, bool left_to_right) {
	if (is_vector2_array_ring(geom2))
		return relation<polygon_t, ring_t>(geom1, geom2, left_to_right);

	return relation<polygon_t, ring_t>(geom1, geom2, left_to_right);
}

String GDBlasGeometry::relation(const PackedVector2Array &geom1, const PackedVector2Array &geom2, bool left_to_right) {
	int cond = static_cast<int>(is_vector2_array_ring(geom1));
	cond |= is_vector2_array_ring(geom2) ? 2 : 0;

	switch (cond) {
	case 0:
		return relation<line_t, line_t>(geom1, geom2, left_to_right);
	case 1:
		return relation<ring_t, line_t>(geom1, geom2, left_to_right);
	case 2:
		return relation<line_t, ring_t>(geom1, geom2, left_to_right);
	case 3:
		return relation<ring_t, ring_t>(geom1, geom2, left_to_right);
	default:
		break;
	}

	return String();
}

String GDBlasGeometry::relation(const Array &geom1, const Vector2 &geom2, bool left_to_right) {
	return relation<polygon_t, point_t>(geom1, geom2, left_to_right);
}

String GDBlasGeometry::relation(const PackedVector2Array &geom1, const Vector2 &geom2, bool left_to_right) {
	if (is_vector2_array_ring(geom1))
		return relation<ring_t, point_t>(geom1, geom2, left_to_right);

	return relation<line_t, point_t>(geom1, geom2, left_to_right);
}

Array GDBlasGeometry::reverse(const Array &geom) {
	return reverse<polygon_t>(geom);
}

PackedVector2Array GDBlasGeometry::reverse(const PackedVector2Array &geom) {
	if (is_vector2_array_ring(geom))
		return reverse<ring_t>(geom);

	return reverse<line_t>(geom);
}

Array GDBlasGeometry::simplify(const Array &geom, scalar_t max_distance) {
	return simplify<polygon_t>(geom, max_distance);
}

PackedVector2Array GDBlasGeometry::simplify(const PackedVector2Array &geom, scalar_t max_distance) {
	if (is_vector2_array_ring(geom))
		return simplify<ring_t>(geom, max_distance);

	return simplify<line_t>(geom, max_distance);
}

Array GDBlasGeometry::sym_difference(const Array &geom1, const Array &geom2, bool left_to_right) {
	return sym_difference<polygon_t, polygon_t>(geom1, geom2, left_to_right);
}

Array GDBlasGeometry::sym_difference(const Array &geom1, const PackedVector2Array &geom2, bool left_to_right) {
	if (is_vector2_array_ring(geom2))
		return sym_difference<polygon_t, ring_t>(geom1, geom2, left_to_right);

	return Array();
}

Array GDBlasGeometry::sym_difference(const PackedVector2Array &geom1, const PackedVector2Array &geom2, bool left_to_right) {
	if (is_vector2_array_ring(geom1) || is_vector2_array_ring(geom2))
		return Array();

	return sym_difference<line_t, line_t>(geom1, geom2, left_to_right);
}

int GDBlasGeometry::touches(const Array &geom) {
	return touches<polygon_t>(geom);
}

int GDBlasGeometry::touches(const PackedVector2Array &geom) {
	if (is_vector2_array_ring(geom))
		return touches<ring_t>(geom);

	return touches<line_t>(geom);
}

int GDBlasGeometry::touches(const Array &geom1, const Array &geom2) {
	return touches<polygon_t, polygon_t>(geom1, geom2);
}

int GDBlasGeometry::touches(const Array &geom1, const PackedVector2Array &geom2) {
	if (is_vector2_array_ring(geom2))
		return touches<polygon_t, ring_t>(geom1, geom2);

	return touches<polygon_t, line_t>(geom1, geom2);
}

int GDBlasGeometry::touches(const PackedVector2Array &geom1, const PackedVector2Array &geom2) {
	int cond = static_cast<int>(is_vector2_array_ring(geom1));
	cond |= is_vector2_array_ring(geom2) ? 2 : 0;

	switch (cond) {
	case 0:
		return touches<line_t, line_t>(geom1, geom2);
	case 1:
		return touches<ring_t, line_t>(geom1, geom2);
	case 2:
		return touches<line_t, ring_t>(geom1, geom2);
	case 3:
		return touches<ring_t, ring_t>(geom1, geom2);
	default:
		break;
	}

	return -1;
}

Array GDBlasGeometry::union_(const Array &geom1, const Array &geom2) {
	return union_<polygon_t, polygon_t>(geom1, geom2);
}

Array GDBlasGeometry::union_(const Array &geom1, const PackedVector2Array &geom2) {
	if (is_vector2_array_ring(geom2))
		return union_<polygon_t, ring_t>(geom1, geom2);

	return Array();
}

Array GDBlasGeometry::union_(const PackedVector2Array &geom1, const PackedVector2Array &geom2) {
	int cond = static_cast<int>(is_vector2_array_ring(geom1));
	cond |= is_vector2_array_ring(geom2) ? 2 : 0;

	switch (cond) {
	case 0:
		return union_<line_t, line_t>(geom1, geom2);
	case 3:
		return union_<ring_t, ring_t>(geom1, geom2);
	default:
		break;
	}

	return Array();
}

Array GDBlasGeometry::unique(const Array &geom) {
	return unique<polygon_t>(geom);
}

PackedVector2Array GDBlasGeometry::unique(const PackedVector2Array &geom) {
	if (is_vector2_array_ring(geom))
		return unique<ring_t>(geom);

	return unique<line_t>(geom);
}

int GDBlasGeometry::within(const Array &geom1, const Array &geom2, bool left_to_right) {
	return within<polygon_t, polygon_t>(geom1, geom2, left_to_right);
}

int GDBlasGeometry::within(const Array &geom1, const PackedVector2Array &geom2, bool left_to_right) {
	if (is_vector2_array_ring(geom2))
		return within<polygon_t, ring_t>(geom1, geom2, left_to_right);

	return within<line_t, polygon_t, false>(geom2, geom1, !left_to_right);
}

int GDBlasGeometry::within(const Array &geom1, const Vector2 &geom2, bool left_to_right) {
	return within<point_t, polygon_t, false>(geom2, geom1, !left_to_right);
}

int GDBlasGeometry::within(const PackedVector2Array &geom1, const PackedVector2Array &geom2, bool left_to_right) {
	int cond = static_cast<int>(is_vector2_array_ring(geom1));
	cond |= is_vector2_array_ring(geom2) ? 2 : 0;

	switch (cond) {
	case 0:
		return within<line_t, line_t>(geom1, geom2, left_to_right);
	case 2:
		return within<line_t, ring_t, false>(geom1, geom2, left_to_right);
	case 3:
		return within<ring_t, ring_t>(geom1, geom2, left_to_right);
	default:
		break;
	}

	return -1;
}

int GDBlasGeometry::within(const PackedVector2Array &geom1, const Vector2 &geom2, bool left_to_right) {
	if (is_vector2_array_ring(geom1))
		return within<point_t, ring_t, false>(geom2, geom1, !left_to_right);

	return within<point_t, line_t, false>(geom2, geom1, !left_to_right);
}
