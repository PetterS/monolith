set(output_file "${CMAKE_CURRENT_BINARY_DIR}/git-revision-internal.h")
string(TIMESTAMP VERSION "%Y-%m-%d %H:%M")

find_program(git_executable NAMES git git.exe git.cmd)
if (git_executable)
	find_file(git_dir NAMES .git PATHS ${PROJECT_SOURCE_DIR} NO_DEFAULT_PATH)
	if (git_dir)
		set(git_state "${git_dir}/logs/HEAD")		
		message("-- Will write version information from ${git_dir} to ${output_file}")
		add_custom_command(
			OUTPUT ${output_file}
	        COMMAND ${CMAKE_COMMAND} -E echo_append "#define REVISION " > ${output_file}
	        COMMAND ${git_executable} log -1 "--pretty=format:%h %ai" >> ${output_file}
	        DEPENDS ${git_state}
	        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
	        VERBATIM
	        COMMENT "Generating git revision file."
	    )
	else()
		message("-- No .git directory found.")
		file(WRITE ${output_file} "#define REVISION \"${VERSION} (CMake)\"")
	endif ()
else ()
	message("-- Git is not available.")
	file(WRITE ${output_file} "#define REVISION \"${VERSION} (CMake)\"")
endif ()

add_library(version version.cpp version.h ${output_file})
set_property(TARGET version PROPERTY FOLDER Version)

add_executable(version_test version_test.cpp)
target_link_libraries(version_test version catch)
set_property(TARGET version_test PROPERTY FOLDER Version)
add_test(NAME version_test COMMAND version_test)
