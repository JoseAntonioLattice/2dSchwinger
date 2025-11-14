module arrays

  use iso_fortran_env, only : dp => real64, i4 => int32

  implicit none

  complex(dp), allocatable, dimension(:,:,:) :: u
  real(dp), allocatable, dimension(:) :: beta 
  real(dp), allocatable, dimension(:) :: plq_action, avr_top, err_top
  real(dp), allocatable, dimension(:,:):: top_char

  real(dp) :: avr_action,err_action
  integer(i4) :: bins
  
end module arrays
