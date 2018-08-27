
include_directories(${OpenCV_INCLUDE_DIRS})

one_target_per_file(
	minimum_core
	${OpenCV_LIBS})
