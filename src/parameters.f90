#if defined(PARALLEL)
#define IFPARALLEL if(this_image()==1)then
#define ENDIFPARALLEL endif
#elif !defined(PARALLEL)
#define IFPARALLEL 
#define ENDIFPARALLEL 
#endif

module parameters

  use iso_fortran_env, only : dp => real64, i4 => int32

  implicit none

  integer(i4) :: L(2), Lx, Ly
  integer(i4) :: N_thermalization
  integer(i4) :: N_measurements
  integer(i4) :: N_skip
  real(dp) :: beta_i, beta_f
  integer(i4) :: n_beta
  real(dp) :: m0
  real(dp) :: tol
  integer(i4) :: max_iter
  integer(i4) :: MD_steps
  real(dp) :: trajectory_length
  namelist /parametersfile/ L,N_thermalization,N_measurements,N_skip, &
       beta_i, beta_f, n_beta,m0,tol, max_iter, MD_steps, trajectory_length
contains

  subroutine read_input()

    integer(i4) :: inunit
    character(99) :: inputfilename

#ifdef PARALLEL
    if(this_image()==1) then
#endif
       read(*,"(a)") inputfilename
       print*, 'User typed: ', trim(inputfilename)
       open(newunit = inunit,file = trim(inputfilename), status = 'old')
       read(inunit, nml = parametersfile)
       
       write(*,nml = parametersfile)
#ifdef PARALLEL
    endif
    call co_broadcast(L,source_image=1)
    call co_broadcast(N_thermalization,source_image=1)
    call co_broadcast(N_measurements,source_image=1)
    call co_broadcast(N_skip,source_image=1)
    call co_broadcast(beta_i,source_image=1)
    call co_broadcast(beta_f,source_image=1)
    call co_broadcast(n_beta,source_image=1)
    call co_broadcast(m0,source_image=1)
    call co_broadcast(tol,source_image=1)
    call co_broadcast(max_iter,source_image=1)
    call co_broadcast(MD_steps,source_image=1)
    call co_broadcast(trajectory_length,source_image=1)
#endif
  end subroutine read_input

  
end module parameters
