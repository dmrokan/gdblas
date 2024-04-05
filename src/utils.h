#ifndef UTILS_H_
#define UTILS_H_

#include <stdint.h>
#include <stdlib.h>

void *gdblas_alloc(size_t s);
void gdblas_free(void *ptr);
void gdblas_error(const char *fmt, ...);
void gdblas_warn(const char *fmt, ...);
void gdblas_debug(const char *fmt, ...);

#endif // UTILS_H_
