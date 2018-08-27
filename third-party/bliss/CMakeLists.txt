
add_library(
	bliss
	bliss/defs.cc
	bliss/graph.cc
	bliss/partition.cc
	bliss/orbit.cc
	bliss/uintseqhash.cc
	bliss/heap.cc
	bliss/timer.cc
	bliss/utils.cc
	bliss/bliss_C.cc)

set_property(TARGET bliss PROPERTY FOLDER "Third-party")
