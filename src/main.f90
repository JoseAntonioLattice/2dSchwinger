program U1_2d

  !use statistics
  use pbc
  use parameters
  use arrays
  use dynamics
  
  implicit none

  integer :: i_b, j

  call read_input
  call set_memory(u,[Lx,Ly],beta,beta_i,beta_f,n_beta,plq_action,top_char,slb_top_char,pion_correlator,n_measurements)
  allocate(avr_top(Lx),err_top(Lx))
  !print*, beta
  !call check_CG()
  !print*,sqrt(2/beta)
  call hot_start(u)
  open( unit = 10, file = 'data/data.dat', status = 'unknown')
  open( unit = 20, file = 'data/pion_correlator.dat', status = 'unknown')
  open( unit = 30, file = 'data/topological_charge.dat', status = 'unknown')
  !go to 100
  do i_b = 1, n_beta
     call initialization(u,plq_action,top_char,slb_top_char,pion_correlator,beta(i_b),N_thermalization,N_measurements, N_skip)
     print*, beta(i_b), avr(plq_action), std_Err(plq_action), avr(top_char),std_err(top_char)
     
     do j = 1, Lx
        write(20,*) j-1,avr(pion_correlator(j,:)),std_err(pion_correlator(j,:))
        write(30,*) j, avr(slb_top_char(j,:)), std_err(slb_top_char(j,:))
     end do
     write(20,*) Lx,avr(pion_correlator(1,:)),std_err(pion_correlator(1,:))
     write(20,"(2/)")
     write(30,"(2/)")

     !print*,beta(i_b), avr_action, err_action, avr_top, err_top
     write(10,*) beta(i_b), avr(plq_action), std_Err(plq_action), avr(top_char),std_err(top_char),&
          avr(top_char**2)/(Lx*Ly),std_err(top_char**2)/(Lx*Ly) 
     flush(10)
     flush(20)
     flush(30)
  end do
  100 print*, "Goodbye World!"
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
