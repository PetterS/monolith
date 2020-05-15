
one_target_per_file(
    absl_strings
	absl_str_format
	minimum_algorithms
	minimum_core
	minimum_linear
	minimum_linear_colgen)

add_test(NAME minimum_linear_colgen_cutting_stock COMMAND cutting_stock default)
add_test(
	NAME minimum_linear_colgen_shift_scheduling_3
	COMMAND shift_scheduling_colgen --nosave_solution --nouse_first_order_solver --num_solutions=2 --verify_solution_value=1300 Instance3.txt
	WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/minimum/linear/data/shift_scheduling
)
add_test(
	NAME minimum_linear_colgen_shift_scheduling_gent
	COMMAND shift_scheduling_colgen --nosave_solution --verify_solution_value=305 not_used.txt --read_gent_instance
	WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/minimum/linear/data/shift_scheduling
)
add_test(
	NAME minimum_linear_colgen_shift_scheduling_primal_dual
	COMMAND shift_scheduling_colgen --nosave_solution --use_first_order_solver --num_solutions=2 --verify_solution_value=610 Instance1.txt
	WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/minimum/linear/data/shift_scheduling
)
add_test(
	NAME minimum_linear_colgen_task_scheduling_1
	COMMAND task_scheduling_colgen --verify_solution_value=46 data_6_48_80_66.dat
	WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/minimum/linear/data/ptask
)
add_test(
	NAME minimum_linear_colgen_task_scheduling_sat
	COMMAND task_scheduling_colgen --sat data_2_24_40_33.dat
	WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/minimum/linear/data/ptask
)
add_test(
	NAME minimum_linear_colgen_shift_generation
	COMMAND shift_generation_colgen --min_rmp_iterations=15 --verify_solution=0.5)
