
if (NOT EMSCRIPTEN)
	library_with_tests("minimum_linear_colgen"
		"minimum_core"
		"minimum_linear")
else()
	# Hack: Leave out minimum_linear for Emscripten because
	# there is a conflict with statically linked protobuf.
	#
	# Emscripten does not do dynamic linking anyway; it
	# saves everything until the end.
	library_with_tests("minimum_linear_colgen"
		"minimum_core")
	set_property(TARGET minimum_linear_colgen PROPERTY INTERFACE_LINK_LIBRARIES "")
endif()

add_subdirectory(examples)
