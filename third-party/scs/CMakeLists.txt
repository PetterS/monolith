
add_library(
	scs
	scs/common.c
	scs/cones.c
	scs/cs.c
	scs/ctrlc.c
	scs/linAlg.c
	scs/private.c
	scs/scs.c
	scs/scs_version.c
	scs/util.c
	scs/amd_1.c
	scs/amd_2.c
	scs/amd_aat.c
	scs/amd_control.c
	scs/amd_defaults.c
	scs/amd_dump.c
	scs/amd_global.c
	scs/amd_info.c
	scs/amd_order.c
	scs/amd_postorder.c
	scs/amd_post_tree.c
	scs/amd_preprocess.c
	scs/amd_valid.c
	scs/ldl.c
	scs/SuiteSparse_config.c
)

set_property(TARGET scs PROPERTY FOLDER "Third-party")
