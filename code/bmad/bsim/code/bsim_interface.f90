module bsim_interface
  
  interface
     subroutine writefile(in_file, parameters)
       use precision_def
       implicit none
       real(rp), dimension(1:,1:), intent(in) ::  parameters
       character(*), intent(in) ::  in_file
     end subroutine writefile
  end interface
  
  interface
     subroutine gaussian_dist (ele, mode, coupling, min_sig, coord)
       use bmad_struct, only: normal_modes_struct, coord_struct, ele_struct, rp
       implicit none
       type (normal_modes_struct) mode
       type (coord_struct), allocatable :: coord(:)
       type (ele_struct) ele
       real(rp) min_sig
       real(rp) coupling
     end subroutine gaussian_dist
  end interface
  
  interface
     subroutine histogram (ele, coord, in_file, sig,a_out)
       use bmad_struct, only: coord_struct, ele_struct, rp
       implicit none
       type (coord_struct) coord(:)
       type (ele_struct) ele
       real(rp) sig(3), a_out(3)
       character(*) in_file
     end subroutine histogram
  end interface
  
integer, private :: private_dummy ! This is to suppress the ranlib "has no symbols" message

end module bsim_interface


