module parameters

  use iso_fortran_env, only : dp => real64, i4 => int32

  implicit none

  integer(i4) :: Lx, Ly
  integer(i4) :: N_thermalization
  integer(i4) :: N_measurements
  integer(i4) :: N_skip
  real(dp) :: beta_i, beta_f
  integer(i4) :: n_beta
  real(dp) :: m0
  real(dp) :: tol
  integer(i4) :: max_iter
  namelist /parametersfile/ Lx,Ly,N_thermalization,N_measurements,N_skip, &
       beta_i, beta_f, n_beta,m0,tol, max_iter
contains

  subroutine read_input()

    integer(i4) :: inunit
    character(99) :: inputfilename
    
    read(*,"(a)") inputfilename
    print*, 'User typed: ', trim(inputfilename)
    open(newunit = inunit,file = trim(inputfilename), status = 'old')
    read(inunit, nml = parametersfile)
    write(*,nml = parametersfile)
  end subroutine read_input

  
end module parameters
