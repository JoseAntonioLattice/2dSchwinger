#if defined(PARALLEL)
#define IFPARALLEL if(this_image()==1)then
#define ENDIFPARALLEL endif
#elif !defined(PARALLEL)
#define IFPARALLEL 
#define ENDIFPARALLEL
#endif

program U1_2d
  !use statistics
  use pbc
  use parameters
  use arrays
  use dynamics
  use number2string
  use starts
  implicit none

  integer :: i_b, j, k
  call random_init(.false.,.true.)
  
  call read_input
  call set_memory(u,L,beta,beta_i,beta_f,n_beta,plq_action,top_char,slb_top_char,pion_correlator,n_measurements)
  allocate(avr_top(Lx),err_top(Lx))
  call hot_start(u)

  IFPARALLEL
  open( unit = 10, file = 'data/data.dat', status = 'unknown')
  open( unit = 20, file = 'data/pion_correlator.dat', status = 'unknown')
  open( unit = 30, file = 'data/topological_charge.dat', status = 'unknown')
  open( unit = 40, file = 'data/polyakov_correlator.dat', status = 'unknown')
  ENDIFPARALLEL

  do i_b = 1, n_beta
     if(read_config) then
        call read_configs(U,beta(i_b))
     else   
        call initialization(u,plq_action,top_char,slb_top_char,pion_correlator,beta(i_b),N_thermalization,N_measurements, N_skip)
     end if
     
     IFPARALLEL
     print*, beta(i_b), avr(plq_action)/product(L), std_Err(plq_action)/product(L), avr(top_char),std_err(top_char)
     do j = 0, L(1)/2-1
        write(40,*) j,avr(poly_corr(j,:)%re),std_err(poly_corr(j,:)%re)
     end do
     write(40,"(2/)")
     flush(40)
     ENDIFPARALLEL
  end do
  
contains
  
  function avr(x)
    real(dp), intent(in) :: x(:)
    real(dp) :: avr

    avr = sum(x)/size(x)
  end function avr
  
  function var(x)
    real(dp), intent(in) :: x(:)
    real(dp) :: var, avg

    avg = avr(x)
    var = sum((x-avg)**2)/(size(x)-1)
  end function var
  
  function std_err(x)
    real(dp), intent(in) :: x(:)
    real(dp) :: std_err

    std_err = sqrt(var(x)/size(x))
  end function std_err
  
end program U1_2d
