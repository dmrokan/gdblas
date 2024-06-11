#ifndef DEFS_H_
#define DEFS_H_

#include <cstdio>

#define GDBLAS_SET_VERSION(a, b, c) (a << 20 | b << 10 | c)

static constexpr uint32_t GDBLAS_VERSION = GDBLAS_SET_VERSION(1, 3, 2);

#ifdef GDBLAS_DEBUG_PRINT

#define GDBLAS_ERROR(...)                                     \
	do {                                                      \
		gdblas_error(__VA_ARGS__, nullptr);                   \
		fprintf(stderr, ": ( %s:%d )\n", __FILE__, __LINE__); \
	} while (0);

#define GDBLAS_WARN(...)                                      \
	do {                                                      \
		gdblas_warn(__VA_ARGS__, nullptr);                    \
		fprintf(stderr, ": ( %s:%d )\n", __FILE__, __LINE__); \
	} while (0);

#define GDBLAS_DEBUG(...)                                     \
	do {                                                      \
		gdblas_debug(__VA_ARGS__, nullptr);                   \
		fprintf(stdout, ": ( %s:%d )\n", __FILE__, __LINE__); \
	} while (0);

#define GDBLAS_V_DEBUG(...)                                   \
	do {                                                      \
		gdblas_debug(__VA_ARGS__, nullptr);                   \
		fprintf(stdout, ": ( %s:%d )\n", __FILE__, __LINE__); \
	} while (0);

#else // GDBLAS_DEBUG_PRINT

#define GDBLAS_ERROR(...)                   \
	do {                                    \
		gdblas_error(__VA_ARGS__, nullptr); \
	} while (0);

#define GDBLAS_WARN(...) \
	do {                 \
	} while (0);

#define GDBLAS_DEBUG(...) \
	do {                  \
	} while (0);

#define GDBLAS_V_DEBUG(...) \
	do {                    \
	} while (0);

#endif

#endif // DEF_H_
