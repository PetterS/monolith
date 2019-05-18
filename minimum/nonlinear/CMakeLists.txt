
library_with_tests("minimum_nonlinear"
	"minimum_core"
	"absl::hash"
	"absl::flat_hash_set"
	${C_MATH_LIBRARY}
	"lbfgsb"
	"meschach")

add_subdirectory(data)
add_subdirectory(examples)
