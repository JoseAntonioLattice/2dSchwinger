program U1_2d

  !use statistics
  use pbc
  use parameters
  use arrays
  use dynamics
  
  implicit none

  integer :: i_b, j

  call read_input
  call set_memory(u,[Lx,Ly],beta,beta_i,beta_f,n_beta,plq_action,top_char,n_measurements)
  allocate(avr_top(Lx),err_top(Lx))
  !print*, beta
  !call check_CG()
  !print*,sqrt(2/beta)
  call hot_start(u,Lx)
  open( unit = 10, file = 'data/data.dat', status = 'unknown')
  open( unit = 20, file = 'data/top.dat', status = 'unknown')
  do i_b = 1, n_beta
     call initialization(u,plq_action,top_char,beta(i_b),N_thermalization,N_measurements, N_skip)
     print*, beta(i_b), avr(plq_action), std_Err(plq_action)
     !write(20,*) beta(i_b), avr(plq_action), std_Err(plq_action)
     !do j = 1, L
     !   call max_jackknife_error_2(top_char(j,:),avr_top(j),err_top(j),bins)
     !   write(20,*) beta(i_b),avr(top_char(Lx,:))/Lx**2,std_err(top_char(Lx,:))/Lx**2
     !end do
     !write(20,"(2/)")


     !print*,beta(i_b), avr_action, err_action, avr_top, err_top
     write(10,*) beta(i_b), avr(plq_action), std_Err(plq_action)!, avr_top, err_top
     flush(10)
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
