#include "inc"
module dynamics
  use iso_fortran_env, only : dp => real64, i4 => int32
  use pbc
  use starts
  use hmc
  use constants, only : i, pi
  use observables
  use configurations
  use number2string
  implicit none
  complex(dp), dimension(2,2,2) :: sigma
#ifdef PARALLEL
  integer(i4),  dimension(2) :: cores 
#endif
contains
 
  subroutine set_memory(u,L,beta,betai,betaf,nbeta, plqaction,top_char,slb_top_char,pion_correlator,n_measurements)
    use parameters, only : Lx, Ly
    use arrays, only : poly_corr
    complex(dp), dimension(:,:,:), allocatable CODIM2 :: u
    integer(i4), intent(in) :: L(2)
    real(dp), allocatable, dimension(:) :: beta
    real(dp), allocatable, dimension(:,:) CODIM2 :: slb_top_char, pion_correlator
    real(dp), intent(in) :: betai, betaf
    real(dp), allocatable, dimension(:) CODIM2 :: plqaction, top_char
    integer(i4), intent(in) :: nbeta, n_measurements
    real(dp), allocatable, dimension(:) :: beta_copy
    integer(i4) :: i_beta

    Lx = L(1)
    Ly = L(2)
#ifdef PARALLEL

    if(this_image() == 1) then
       print*, "Enter cores arrays"
       read(*,*) cores
       print*, "user typed # cores", cores
       
    end if
    call co_broadcast(cores,source_image=1)
    sync all
    
    Lx = Lx/cores(1)
    Ly = Ly/cores(2)
    call set_pbc([Lx,Ly],cores)
    
#else
    call set_pbc([Lx,Ly])
#endif
    allocate(Up(DIM)CDIM)
    allocate(psi(DIM)CDIM)
    allocate(phi(DIM)CDIM)
    allocate(cg_p(DIM)CDIM)
    allocate(cg_Ddagp(DIM)CDIM)
    allocate(chi(DIM)CDIM)
    
#ifdef PARALLEL   
    !print*, this_image(), left, right, up, down, left_down, left_up, right_down, right_up
    print*, this_image(), ",", im([1,1],1), im([1,1],2)
    print*, this_image(), ",", ip([Lx,Ly],1), ip([Lx,Ly],2)
#endif
    allocate(u(DIM)CDIM)
    allocate(beta(nbeta))

    if(nbeta/=1)then
       do i_beta = 1, nbeta 
          beta(i_beta) = betai + (betaf - betai)/(nbeta-1)*(i_beta-1)
       end do
    else
       beta(1) = betai
    end if
    sigma(1,:,:) = reshape([0.0_dp,1.0_dp,1.0_dp,0.0_dp],[2,2])
    sigma(2,:,:) = reshape([(0.0_dp,0.0_dp),i,-i,(0.0_dp,0.0_dp)],[2,2])


    allocate(plqaction(n_measurements)CDIM)

    allocate(top_char(N_measurements)CDIM)
    allocate(slb_top_char(L(1),n_measurements)CDIM)
    allocate(pion_correlator(L(1),n_measurements)CDIM)
    allocate(poly_corr(0:L(1)/2-1,n_measurements)CDIM)
  end subroutine set_memory
  
  subroutine initialization(u,plqaction,top_char,slb_top_char,pion_correlator,beta,N_thermalization,N_measurements, N_skip)
    use parameters, only : L, Lx, Ly, save_config
    use arrays, only : poly_corr
    complex(dp), dimension(DIM) CODIM, intent(inout) :: u
    real(dp), dimension(:) CODIM, intent(out) :: plqaction, top_char
    integer(i4), intent(in) :: N_thermalization, N_measurements, N_skip
    real(dp), dimension(:,:) CODIM, intent(out) :: slb_top_char, pion_correlator
    real(dp), intent(in) :: beta
    real(dp) :: top(Lx,Ly)
    integer(i4) :: i_skip, i_sweeps, j

  
    call thermalization(u, N_thermalization, beta)
    
    do i_sweeps = 1, N_measurements
       do i_skip = 1, n_skip
          call sweeps(u,beta)
       end do
       !plqaction(i_sweeps) = action(u)
       !poly_corr(0:,i_sweeps) = correlation_polyakov(u)
       if(save_config) call save_configuration(U,beta)       
#ifdef PARALLEL
       sync all
#endif
    end do
    
  end subroutine initialization


  subroutine read_configs(U,beta)
    use parameters
    use arrays, only : plq_action, poly_corr
    complex(dp), dimension(DIM) CODIM, intent(inout) :: U
    real(dp), intent(in) :: beta
    integer :: k
    logical :: condition
    character(:), allocatable :: filename
    
    k = 0
    do 
       k = k + 1
       filename = "data/configurations/m0="//real2str(m0,nint(log(abs(m0))+1),4)//"/Lx="//int2str(L(1))// &
            "/Lt="//int2str(L(2))//"/beta="//real2str(beta,nint(log(abs(beta))+1),4)//"/U_"//int2str(k)//".bin"
       inquire(file = filename, exist = condition)
       if(condition)then
          call read_configuration(U,filename)
       else
          exit
       end if
       plq_action(k) = action(u)
       poly_corr(0:,k) = correlation_polyakov(u)
    end do
    
  end subroutine read_configs
  
  subroutine thermalization(u,N_thermalization,beta)
    use parameters, only : Lx, Ly
    complex(dp), dimension(DIM) CODIM, intent(inout) :: u
    integer(i4) :: N_thermalization
    real(dp) :: beta
    integer(i4) :: i_sweeps

    do i_sweeps = 1, N_thermalization
       call sweeps(u,beta)
    end do
    
  end subroutine thermalization

  subroutine sweeps(u,beta)
    use parameters, only : Lx, Ly
    complex(dp), dimension(DIM) CODIM, intent(inout) :: u
    real(dp), intent(in) :: beta
    integer(i4) :: x,y,mu

    call hmc_1(U,beta)
       
  end subroutine sweeps


end module dynamics
