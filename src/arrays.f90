module arrays

  use iso_fortran_env, only : dp => real64, i4 => int32

  implicit none

  complex(dp), allocatable, dimension(:,:,:) :: u
  real(dp), allocatable, dimension(:) :: beta 
  real(dp), allocatable, dimension(:) :: plq_action,top_char, avr_top, err_top, avr_pion, err_pion
  real(dp), allocatable, dimension(:,:) :: slb_top_char
  real(dp), allocatable, dimension(:,:) :: pion_correlator
  real(dp) :: avr_action,err_action
  integer(i4) :: bins
  
end module arrays
