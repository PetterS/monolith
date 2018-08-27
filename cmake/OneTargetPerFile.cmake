

function (one_target_per_file)
	file(GLOB HEADER_FILES *.h)
	file(GLOB CPP_FILES *.cpp)

	foreach (example_file ${CPP_FILES})
		get_filename_component(name ${example_file} NAME_WE)
		add_executable(
			${name}
			${name}.cpp
			${HEADER_FILES})
		foreach(dep ${ARGN})
        	target_link_libraries(${name} ${dep})
        endforeach()
		string(REGEX REPLACE ".*/master/" "" ide_path ${CMAKE_CURRENT_SOURCE_DIR})
		set_property(TARGET ${name} PROPERTY FOLDER ${ide_path})
		message("-- Added ${ide_path}/${name}.")
	endforeach()
endfunction()
