#if defined(PARALLEL)
#define CODIM ,codimension[:]
#else
#define CODIM 
#endif
module arrays

  use iso_fortran_env, only : dp => real64, i4 => int32

  implicit none

  complex(dp), allocatable, dimension(:,:,:) CODIM :: u
  real(dp), allocatable, dimension(:) :: beta 
  real(dp), allocatable, dimension(:) CODIM :: plq_action,top_char
  real(dp), allocatable, dimension(:) :: avr_top, err_top, avr_pion, err_pion
  real(dp), allocatable, dimension(:,:) CODIM :: slb_top_char
  real(dp), allocatable, dimension(:,:) CODIM :: pion_correlator
  real(dp) :: avr_action,err_action
  integer(i4) :: bins
  
end module arrays
