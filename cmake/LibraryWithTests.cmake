

function (library_with_tests NAME)
	file(GLOB HEADER_FILES *.h)
	file(GLOB CPP_FILES *.cpp)
	file(GLOB TEST_FILES *_test.cpp)
	file(GLOB FUZZ_FILES *_fuzz.cpp)
	file(GLOB COMPILE_FAIL_FILES *_compilefail.cpp)
	file(GLOB PROTO_FILES *.proto)

	foreach (TEST_FILE ${TEST_FILES})
		list(REMOVE_ITEM CPP_FILES ${TEST_FILE})
	endforeach()
	foreach (FUZZ_FILE ${FUZZ_FILES})
		list(REMOVE_ITEM CPP_FILES ${FUZZ_FILE})
	endforeach()
	foreach (COMPILE_FAIL_FILE ${COMPILE_FAIL_FILES})
		list(REMOVE_ITEM CPP_FILES ${COMPILE_FAIL_FILE})
	endforeach()

	string(REGEX REPLACE ".*/master/" "" ide_path ${CMAKE_CURRENT_SOURCE_DIR})

	foreach(PROTO_NAME ${PROTO_FILES})
		get_filename_component(only_name ${PROTO_NAME} NAME)
		set(proto_file_in_binary_dir "${CMAKE_CURRENT_BINARY_DIR}/${only_name}")
		string(REPLACE ".proto" ".pb.cc" output_cc_file ${proto_file_in_binary_dir})
		string(REPLACE ".proto" ".pb.h" output_h_file ${proto_file_in_binary_dir})
		string(TOUPPER ${NAME} library_name_uppercase)
		if (EMSCRIPTEN OR DEFINED PROTOC)
			add_custom_command(
				OUTPUT ${output_cc_file} ${output_h_file}
				COMMAND ${PROTOC} ${PROTO_NAME} --cpp_out=dllexport_decl=${library_name_uppercase}_API:${CMAKE_BINARY_DIR} --python_out=${CMAKE_BINARY_DIR} --proto_path=${CMAKE_SOURCE_DIR} --proto_path=${CMAKE_SOURCE_DIR}/third-party/protobuf/src
				DEPENDS ${PROTO_NAME}
				WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
			)
		else()
			add_custom_command(
				OUTPUT ${output_cc_file} ${output_h_file}
				COMMAND ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/protoc ${PROTO_NAME} --cpp_out=dllexport_decl=${library_name_uppercase}_API:${CMAKE_BINARY_DIR} --python_out=${CMAKE_BINARY_DIR} --proto_path=${CMAKE_SOURCE_DIR} --proto_path=${CMAKE_SOURCE_DIR}/third-party/protobuf/src
				DEPENDS ${PROTO_NAME} protoc
				WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
			)
		endif()
		list(APPEND library_pb_files ${output_cc_file} ${output_h_file})
		if (MSVC)
			set_source_files_properties(
				${output_cc_file}
				PROPERTIES
				COMPILE_DEFINITIONS "${library_name_uppercase}_API=__declspec(dllexport)"
			)
		else()
			set_source_files_properties(
				${output_cc_file}
				PROPERTIES
				COMPILE_DEFINITIONS ${library_name_uppercase}_API=
			)
		endif()
	endforeach()


	if (CPP_FILES)
		add_library(
			${NAME}
			SHARED
			${CPP_FILES}
			${HEADER_FILES}
			${library_pb_files})
		set_property(TARGET ${NAME} PROPERTY FOLDER ${ide_path})
		foreach(DEP ${ARGN})
			target_link_libraries(${NAME} ${DEP})
		endforeach()
	endif()

	foreach(TEST_FILE ${TEST_FILES})
		get_filename_component(TEST_NAME ${TEST_FILE} NAME_WE)
		add_executable(
			${TEST_NAME}
			${TEST_NAME}.cpp
			${HEADER_FILES})
		target_link_libraries(${TEST_NAME} catch)
		if (CPP_FILES)
		   		target_link_libraries(${TEST_NAME} ${NAME})
		endif()
		foreach(DEP ${ARGN})
			target_link_libraries(${TEST_NAME} ${DEP})
		endforeach()
		set_property(TARGET ${TEST_NAME} PROPERTY FOLDER ${ide_path}/Tests)
		add_test(
			NAME ${NAME}_${TEST_NAME}
			COMMAND ${TEST_NAME}
			WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
	endforeach()

	# https://stackoverflow.com/a/30191576/116186
	#
	# But these tests are currently not used since they take several
	# seconds to run each.
	foreach(COMPILE_FAIL_FILE ${COMPILE_FAIL_FILES})
		get_filename_component(TEST_NAME ${COMPILE_FAIL_FILE} NAME_WE)
		add_executable(
			${TEST_NAME}
			${TEST_NAME}.cpp
			${HEADER_FILES})
		if (CPP_FILES)
		   		target_link_libraries(${TEST_NAME} ${NAME})
		endif()
		foreach(DEP ${ARGN})
			target_link_libraries(${TEST_NAME} ${DEP})
		endforeach()
		set_property(TARGET ${TEST_NAME} PROPERTY FOLDER ${ide_path}/CompileFail)


		# Avoid building these targets normally
		set_target_properties(${TEST_NAME} PROPERTIES
			EXCLUDE_FROM_ALL TRUE
			EXCLUDE_FROM_DEFAULT_BUILD TRUE)
		add_test(
			NAME ${NAME}_${TEST_NAME}
			COMMAND ${CMAKE_COMMAND} --build . --target ${TEST_NAME} --config $<CONFIGURATION>
			WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
		set_tests_properties(${NAME}_${TEST_NAME} PROPERTIES WILL_FAIL TRUE)
	endforeach()

	if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang" AND CMAKE_BUILD_TYPE STREQUAL "Debug")
		foreach(FUZZ_FILE ${FUZZ_FILES})
			get_filename_component(FUZZ_NAME ${FUZZ_FILE} NAME_WE)
			add_executable(
				${FUZZ_NAME}
				${FUZZ_NAME}.cpp
				${HEADER_FILES})
			target_link_libraries(${FUZZ_NAME} "-fsanitize=fuzzer")
			if (CPP_FILES)
			   		target_link_libraries(${FUZZ_NAME} ${NAME})
			endif()
			foreach(DEP ${ARGN})
				target_link_libraries(${FUZZ_NAME} ${DEP})
			endforeach()
			set_property(TARGET ${FUZZ_NAME} PROPERTY FOLDER ${ide_path}/Fuzz)

			message("-- Added fuzz target ${FUZZ_NAME}.")
		endforeach()
	endif()

	message("-- Added ${NAME}.")
endfunction()
