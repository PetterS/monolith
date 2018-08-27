# Petter Strandmark

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_library(glucose
	core/Solver.cc
	simp/SimpSolver.cc)

set_property(TARGET glucose PROPERTY FOLDER "Third-party/Glucose")
