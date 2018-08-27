# Petter Strandmark

add_definitions("-DNLGL")
add_library(cleaneling
	libcleaneling.cc
	options.h
	cleaneling.h)

set_property(TARGET cleaneling PROPERTY FOLDER "Third-party")
