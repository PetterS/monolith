cmake_minimum_required(VERSION 3.12)
project(Minimum C CXX)

include(CheckCCompilerFlag)

enable_testing()
set(CTEST_OUTPUT_ON_FAILURE TRUE)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)  # For Clang-based tools.
set_property(GLOBAL PROPERTY USE_FOLDERS ON)
if (NOT CMAKE_BUILD_TYPE)
  message("-- No build type specified; defaulting to RelWithDebInfo.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING
    "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
    FORCE)
else()
	message("-- Requested build type is ${CMAKE_BUILD_TYPE}.")
endif ()
set(CMAKE_POSITION_INDEPENDENT_CODE TRUE)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if (MSVC)
  # Use the latest C++ standard.
  add_definitions("/std:c++latest")	
endif()

include(cmake/LibraryWithTests.cmake)
include(cmake/OneTargetPerFile.cmake)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/archive)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR}/archive)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR}/bin)

#
# OpenMP
#
option(OPENMP
       "Enable multi-threading (requires OpenMP)"
       ON)
if (OPENMP)
  find_package(OpenMP)
  if (OPENMP_FOUND)
    message("-- Found OpenMP.")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    add_definitions(-DUSE_OPENMP)
    add_definitions(-DEIGEN_HAS_OPENMP=1)
  else()
    message("-- Can't find OpenMP. Continuing without it.")
  endif()
else ()
  message("-- OpenMP is disabled. Continuing without it.")
endif ()

set(C_MATH_LIBRARY "")

#
# OpenCV
#
find_package(OpenCV QUIET)
if (OpenCV_FOUND)
	message("-- Found OpenCV.")
else ()
	message("-- Did not find OpenCV (optional). Set OpenCV_DIR to the location of OpenCVConfig.cmake if desired.")
endif ()

#
# Clang and GCC settings
#
if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_COMPILER_IS_GNUCXX)
	# Run cmake with
	#   cmake -DCMAKE_CXX_FLAGS="-march=native -mtune=native -mavx"
	# to enable processor-specific optimizations with GCC.

	# TODO: Enable sign-compare warning.
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wno-sign-compare -Wno-unused-function")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-sign-compare -Wno-unused-function")	

	if (CMAKE_BUILD_TYPE STREQUAL "Debug" AND NOT EMSCRIPTEN)
		if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
			message("-- Debug mode enabled for Clang; Adding sanitize flags.")
			set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=undefined,address -fno-sanitize-recover=undefined -fno-sanitize=float-divide-by-zero,unsigned-integer-overflow -fno-omit-frame-pointer -fsanitize-coverage=edge,trace-pc-guard,indirect-calls,trace-cmp,trace-div,trace-gep")
			set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=undefined,address -fno-sanitize-recover=undefined -fno-sanitize=float-divide-by-zero,unsigned-integer-overflow -fno-omit-frame-pointer -fsanitize-coverage=edge,trace-pc-guard,indirect-calls,trace-cmp,trace-div,trace-gep")
		elseif(CMAKE_COMPILER_IS_GNUCXX)
			# Gcov does not seem to play nicely with sanitize failures. Add Gcov only
			# for GCC.
			message("-- Debug mode enabled for Gcc; adding support for Gcov.")
			set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fprofile-arcs -ftest-coverage")
			set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
		endif()
	endif ()

	set(C_MATH_LIBRARY "m")

	# Clang and GCC doesn't print colored diagnostics when invoked from Ninja
	if (CMAKE_GENERATOR STREQUAL "Ninja")
		if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
			set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fcolor-diagnostics")
			set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fcolor-diagnostics")
		endif()
		if (CMAKE_COMPILER_IS_GNUCXX)
			set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fdiagnostics-color=always")
			set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fdiagnostics-color=always")
		endif()
	endif()
endif ()

# Some third-party C libraries use pthread. For example, Clang 3.6 requires
# this flag (but not GCC 5.2). Just add it if it is available.
check_c_compiler_flag(-pthread HAS_PTHREAD)
if (HAS_PTHREAD)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pthread")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
endif()

#
# MSVC settings
#
if (MSVC)
  # Disable deprecation warning for standard functions.
  add_definitions("/wd4996")
  # To get rid of annoying min and max macros if windows.h
  # is included.
  add_definitions("-DNOMINMAX=1")
  # Macro redefinitions.
  add_definitions("/wd4005")
  # Type conversions.
  add_definitions("/wd4244")
  add_definitions("/wd4267")
  # Signed/unsigned.
  add_definitions("/wd4018")
  add_definitions("/wd4146")
  # https://github.com/google/protobuf/tree/master/cmake
  add_definitions("/wd4251")
  # Non dll-interface class 'google::protobuf::Message' used as base for dll-interface class
  add_definitions("/wd4275")

  # Workaround for ICE in Eigen.
  add_definitions("-DEIGEN_DONT_PARALLELIZE")

  # Or-tools can fail in debug mode without this.
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj")
endif ()


#
# Emscripten settings
#
if (EMSCRIPTEN)
	# Make Emscripten output html pages.
	# set(CMAKE_EXECUTABLE_SUFFIX ".html")

	# To be able to detect Emscripten in code.
	add_definitions("-DEMSCRIPTEN=1")

	set(EMSCRIPTEN_FEATURES "")
	# We want exceptions.
	set(EMSCRIPTEN_FEATURES "${EMSCRIPTEN_FEATURES} -s DISABLE_EXCEPTION_CATCHING=0")
	# We want Webassembly.
	set(EMSCRIPTEN_FEATURES "${EMSCRIPTEN_FEATURES} -s WASM=1 -s ALLOW_MEMORY_GROWTH=1")
	# Filesystem support.
	set(EMSCRIPTEN_FEATURES "${EMSCRIPTEN_FEATURES} -s FORCE_FILESYSTEM=1 -lidbfs.js")
	# Expose callMain since we use it in our client.
	set(EMSCRIPTEN_FEATURES "${EMSCRIPTEN_FEATURES} -s EXTRA_EXPORTED_RUNTIME_METHODS='[\"callMain\"]'")
	if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
		# Do not remove whitespace etc. in the JS parts (easier debugging).
		set(EMSCRIPTEN_FEATURES "${EMSCRIPTEN_FEATURES} -g1")
	endif()
	# We don't want the page to hang while executing.
	# There will be some effort to make this work; e.g. creating files.
	set(EMSCRIPTEN_FEATURES "${EMSCRIPTEN_FEATURES} -s PROXY_TO_WORKER=1 --proxy-to-worker")
	# Add preamble.
	set(EMSCRIPTEN_FEATURES "${EMSCRIPTEN_FEATURES} --pre-js ${CMAKE_CURRENT_SOURCE_DIR}/emscripten/pre.js")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${EMSCRIPTEN_FEATURES}")
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${EMSCRIPTEN_FEATURES}")
	set(CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS} ${EMSCRIPTEN_FEATURES}")

	# Need an external protoc.
	if (NOT PROTOC)
		message(FATAL_ERROR "PROTOC needs to be set for Emscripten." )
	endif()

	# Extra things to make various detections in third-party/ work.
	set(SIZEOF_INT_P_run_result 4)
	set(GFLAGS_INTTYPES_FORMAT C99)
endif ()

find_package(BZip2)
find_package(ZLIB)
macro(link_to_coin TARGET)
  target_link_libraries(${TARGET} Cbc CbcSolver Cgl Clp CoinUtils Osi OsiClp)
  
  if (BZIP2_FOUND)
    target_link_libraries(${TARGET} ${BZIP2_LIBRARIES})
  endif ()
  
  if (ZLIB_FOUND)
    target_link_libraries(${TARGET} ${ZLIB_LIBRARIES})
  endif ()
endmacro()

#
# Source code
#
include_directories(
	${CMAKE_CURRENT_SOURCE_DIR}

	SYSTEM ${PROJECT_BINARY_DIR}  # Generated files (e.g. protocol buffers).

	SYSTEM ${PROJECT_BINARY_DIR}/third-party/gflags/include
	# SYSTEM ${PROJECT_BINARY_DIR}/third-party/or-tools
	SYSTEM third-party/abseil-cpp
	SYSTEM third-party/bliss
	SYSTEM third-party/Catch
	SYSTEM third-party/Eigen
	SYSTEM third-party/FADBAD++
	SYSTEM third-party/glpk/src
	SYSTEM third-party/glucose-syrup
	SYSTEM third-party/L-BFGS-B-C/src
	SYSTEM third-party/meschach
	SYSTEM third-party/minisat
	# SYSTEM third-party/or-tools/src
	SYSTEM third-party/protobuf/src
	SYSTEM third-party/scs
	SYSTEM third-party/sqlite3
	SYSTEM third-party/sym-ildl
)


# For minisat.
add_definitions(-DNO_GMP -D_FILE_OFFSET_BITS=64)
# For or-tools.
add_definitions("-D_SILENCE_STDEXT_HASH_DEPRECATION_WARNINGS=1")
# All of our compilers support shared_ptr. This definition
# is for protobuf.
add_definitions("-DUTIL_GTL_USE_STD_SHARED_PTR=1")

message("--")
message("-- Main configuration done; adding subdirectories.")
message("--")

add_subdirectory(third-party)
add_subdirectory(minimum)
add_subdirectory(misc)
add_subdirectory(version)
if (OpenCV_FOUND)
	add_subdirectory(vision)
endif()
