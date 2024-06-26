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

	Variant get_version();
};

} //namespace godot

#endif // GDBLAS_H_
