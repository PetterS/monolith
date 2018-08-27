
add_library(catch main_test_program.cpp catch.hpp)
set_property(TARGET catch PROPERTY FOLDER "Third-party")
