
add_library(lbfgsb
	lbfgsb.c
	linesearch.c
	subalgorithms.c
	print.c
	linpack.c
	miniCBLAS.c
	timer.c)
set_property(TARGET lbfgsb PROPERTY FOLDER "Third-party/lbfgsb")

add_executable(lbfgsb_driver1
	driver1.c)
target_link_libraries(lbfgsb_driver1 lbfgsb ${C_MATH_LIBRARY})
set_property(TARGET lbfgsb_driver1 PROPERTY FOLDER "Third-party/lbfgsb")
add_test(NAME lbfgsb_driver1 COMMAND lbfgsb_driver1)
