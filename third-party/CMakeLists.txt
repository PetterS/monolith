
if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_COMPILER_IS_GNUCXX)
	# No warnings for C. The Meschach library contains really old code.
	set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -w")
	# No warnings for C++. Too noisy with Coin code and too many to fix atm.
	set(CMAKE_CXX_FLAGS   "${CMAKE_CXX_FLAGS} -w")
endif ()

# Absl flags.
set(ABSL_STD_CXX_FLAG "" CACHE STRING "c++ std flag" FORCE)

add_subdirectory(abseil-cpp)
add_subdirectory(bliss)
add_subdirectory(cadical)
add_subdirectory(Catch)
add_subdirectory(cleaneling)
add_subdirectory(coin)
add_subdirectory(duktape)
add_subdirectory(gflags)
add_subdirectory(glpk)
add_subdirectory(glucose-syrup)
add_subdirectory(L-BFGS-B-C)
add_subdirectory(meschach)
add_subdirectory(minisat)
add_subdirectory(minisatp)
add_subdirectory(protobuf/cmake)
add_subdirectory(scs)
add_subdirectory(sqlite3)
