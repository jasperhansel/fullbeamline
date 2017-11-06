!+
! This module was originally in obj_src/xlinbcg.f90.
! This was a problem since this forced compiliation of the xlinbcg program and
! then the cmake system would put the program in the library.
! The solution was to separate the module from the program.
!-


	MODULE xlinbcg_data
	USE nrtype
	TYPE(sprs2_dp) :: sa
	END MODULE xlinbcg_data

