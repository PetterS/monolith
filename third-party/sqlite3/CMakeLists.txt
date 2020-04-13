
add_library(sqlite3 SHARED sqlite3.c)
if (UNIX)
	target_link_libraries(sqlite3 dl)
endif (UNIX)
if (MSVC)
	target_compile_definitions(sqlite3 PRIVATE "SQLITE_API=__declspec(dllexport)")
endif (MSVC)

add_executable(sqlite3_shell shell.c)
target_link_libraries(sqlite3_shell sqlite3)

find_path(Readline_INCLUDE_DIR
    NAMES readline/readline.h
)
find_library(Readline_LIBRARY
    NAMES readline
)
find_library(Ncurses_LIBRARY
    NAMES ncurses
)
if (Readline_INCLUDE_DIR AND Readline_LIBRARY AND Ncurses_LIBRARY)
	message("-- Building Sqlite shell with readline support.")
	target_include_directories(sqlite3_shell PRIVATE ${Readline_INCLUDE_DIR})
	target_link_libraries(sqlite3_shell ${Readline_LIBRARY} ${Ncurses_LIBRARY})
	target_compile_definitions(sqlite3_shell PRIVATE "-DHAVE_READLINE")
endif()

set_property(TARGET sqlite3 PROPERTY FOLDER "Third-party")
set_property(TARGET sqlite3_shell PROPERTY FOLDER "Third-party")
