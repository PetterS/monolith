
add_library(sqlite3 SHARED sqlite3.c)
if (UNIX)
	target_link_libraries(sqlite3 dl)
endif (UNIX)
if (MSVC)
	target_compile_definitions(sqlite3 PRIVATE "SQLITE_API=__declspec(dllexport)")
endif (MSVC)

add_executable(sqlite3_shell shell.c)
target_link_libraries(sqlite3_shell sqlite3)

set_property(TARGET sqlite3 PROPERTY FOLDER "Third-party")
set_property(TARGET sqlite3_shell PROPERTY FOLDER "Third-party")
