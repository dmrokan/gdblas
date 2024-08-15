#if defined(GDBLAS_WITH_ODE) || defined(GDBLAS_WITH_GEOMETRY)

#define BOOST_NO_EXCEPTIONS
#include <boost/exception/diagnostic_information.hpp>
#include <boost/throw_exception.hpp>
#include <godot_cpp/variant/string.hpp>
#include <godot_cpp/classes/ref.hpp>

#include "gdblas.h"

void boost::throw_exception(std::exception const & e, boost::source_location const & loc) {
	const std::string info = boost::diagnostic_information(e);
	const godot::String error_info(info.c_str());
	godot::Ref<godot::GDBlas> tmp;
	tmp.instantiate();
	tmp->emit_signal("gdblas_boost_exception", error_info);
}

void boost::throw_exception(std::exception const & e) {
	const source_location loc{};
	throw_exception(e, loc);
}

#endif
