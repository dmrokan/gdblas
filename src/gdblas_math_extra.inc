namespace std {
typedef GDBlasMat::complex_t complex_t;
typedef GDBlasMat::scalar_t scalar_t;
static constexpr scalar_t GDBlas_log_2 = 0.693147180559945309417232121458;

_ALWAYS_INLINE_ complex_t abs2(const complex_t &a) {
	complex_t b;
	b.real(std::abs(a));
	b.imag(0);

	return b;
}

_ALWAYS_INLINE_ scalar_t abs2(scalar_t a) {
	return std::abs(a);
}

_ALWAYS_INLINE_ complex_t cbrt(const complex_t &a) {
	scalar_t mag = std::abs(a);
	scalar_t arg = std::arg(a);

	return std::polar(std::cbrt(mag), arg / 3.0);
}

_ALWAYS_INLINE_ complex_t log2(const complex_t &a) {
	complex_t b = std::log(a) / GDBlas_log_2;

	return b;
}
} //namespace std
