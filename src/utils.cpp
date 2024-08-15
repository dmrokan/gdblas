#include "utils.h"

#include <cstdio>
#include <cstdarg>

#include <godot_cpp/variant/utility_functions.hpp>

static constexpr size_t MAX_MESSAGE_SIZE = 512;

void *gdblas_alloc(size_t s) {
	if (!s)
		return NULL;

	return calloc(1, s);
}

void gdblas_free(void *ptr) {
	if (ptr) {
		free(ptr);
	}
}

void gdblas_error(const char *fmt, ...) {
	va_list args;
	using namespace godot;

	va_start(args, fmt);
	PackedByteArray bytes;
	bytes.resize(MAX_MESSAGE_SIZE);
	vsnprintf(reinterpret_cast<char *>(bytes.ptrw()), MAX_MESSAGE_SIZE, fmt, args);
	va_end(args);

	UtilityFunctions::printerr(bytes.get_string_from_ascii());

#if defined(GDBLAS_DEBUG_PRINT)
	fprintf(stderr, "ERR: ");
	va_start(args, fmt);
	vfprintf(stderr, fmt, args);
	va_end(args);
#endif
}

void gdblas_warn(const char *fmt, ...) {
#if defined(GDBLAS_DEBUG_PRINT)
	fprintf(stderr, "WRN: ");
	va_list args;
	va_start(args, fmt);
	vfprintf(stderr, fmt, args);
	va_end(args);
#endif
}

void gdblas_debug(const char *fmt, ...) {
#if defined(GDBLAS_DEBUG_PRINT)
	fprintf(stdout, "DBG: ");
	va_list args;
	va_start(args, fmt);
	vfprintf(stdout, fmt, args);
	va_end(args);
#endif
}
