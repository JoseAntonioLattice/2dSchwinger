program main
  use format
  use iso_fortran_env, only : dp => real64
  implicit none

  integer, parameter :: Lx = 32
  integer, parameter :: Lt = 8

  integer, parameter :: nconfig = 1000

  real(dp) :: beta = 6.0_dp
  real(dp) :: m0 = 100.0_dp
  
  complex(dp), dimension(2,Lx,Lt) :: U
  character(:), allocatable :: filename

  integer :: i

  real(dp) :: action(nconfig)
  
  call set_pbc([Lx,Lt])

  do i = 1, nconfig
     filename = "data/configurations/m0="//real2str(m0,4)// &
          "/Lx="//int2str(Lx)//"/Lt="//int2str(Lt)//&
          "/beta="//real2str(beta)// &
          "/U_"//int2str(i)//".bin"
     call read_configuration(U,filename)
     a_action(i) = action(U)
  end do
  
contains

  subroutine read_configuration(U,filename)
    use number2string
    complex(dp), intent(out) :: U(:,:,:)
    character(*), intent(in) :: filename 
    integer(i4) :: un
    
    open(newunit = un, file = filename, access = "sequential", form = "unformatted")  
    read(un) U
    close(un)
       
  end subroutine read_configuration

  function plaquette(u,x)
    use pbc, only : ip, im
    complex(dp) :: plaquette
    complex(dp), dimension(:,:,:), intent(in) :: u
    integer(i4), intent(in) :: x(2)
    integer(i4), dimension(2) :: x2, x3


    x2 = ip(x,1)
    x3 = ip(x,2)
    
    plaquette = U(1,x(1),x(2)) * U(2,x2(1),x2(2)) * &
          conjg(U(1,x3(1),x3(2))) * conjg(U(2,x(1),x(2)))
    
  end function plaquette

  function action(u)
    real(dp) :: action
    complex(dp), dimension(:,:,:), intent(in) :: u
    integer(i4) :: x,y

    action = 0.0_dp
    
    do x = 1, size(u(1,:,1))
       do y = 1,  size(u(1,1,:))
          action = action + real(plaquette(u,[x,y]))
       end do
    end do

  end function action
  
  function correlation_polyakov(U) result(corr_poly)
    complex(dp), intent(in) :: U(:,:,:)
    integer(i4) :: x, r, t
    complex(dp), dimension(0:size(U(1,:,1))/2-1) :: corr_poly
    complex(dp), dimension(size(U(1,:,1))) :: polyakov_loop
    integer :: Lx

    Lx = size(U(1,:,1))
    do x = 1, Lx
       polyakov_loop(x) = product(U(2,x,:))
    end do
    
    corr_poly = 0.0_dp
    do t = 0, Lx/2 - 1
       do x = 1, Lx
          r = mod(x-1+t,Lx)+1
          corr_poly(t) = corr_poly(t)+polyakov_loop(x)*conjg(polyakov_loop(r))
       end do
    end do
    corr_poly = corr_poly/Lx
  end function correlation_polyakov

  subroutine topological_charge_density(u,top)
    complex(dp), dimension(DIM), intent(in) :: u
    real(dp), intent(out) :: top(size(u(1,:,1)),size(u(1,1,:))) 
    integer(i4) :: i, j, Lx, Ly
    complex(dp) :: plq

    Lx = size(u(1,:,1))
    Ly = size(u(1,1,:))
    
    do i = 1, Lx
       do j = 1, Ly
          plq = plaquette(u,[i,j])
          top(i,j) = atan2(plq%im,plq%re)!/(2*pi)
       end do
    end do
       
  end subroutine topological_charge_density

  function slab_top_char(top,ix)
    real(dp), intent(in) :: top(:,:)
    integer(i4), intent(in) :: ix
    real(dp) :: slab_top_char
    integer(i4) :: i, j
    complex(dp) :: plq

    slab_top_char = (sum(top(1:ix,:))/(2*pi))**2!*sum(top(ix+1:,:))
    
  end function slab_top_char


  
!!$  subroutine wilson_flow(U,beta)
!!$    use parameters, only : Lx, Ly
!!$    complex(dp), intent(inout) :: U(DIM)
!!$    real(dp), intent(in) :: beta
!!$    real(dp), parameter :: epsilon = 1.0E-3_dp
!!$    integer, parameter :: Nt = 100
!!$    integer :: x,y,t
!!$    
!!$    do t = 1, Nt
!!$       do x = 1, Lx
!!$          do y = 1, Ly
!!$             Up(1,x,y) = exp(-epsilon*beta*real(U(1,x,y)*conjg(staples(U,[x,y],1))))*U(1,x,y)
!!$             Up(2,x,y) = exp(-epsilon*beta*real(U(2,x,y)*conjg(staples(U,[x,y],2))))*U(2,x,y)
!!$          end do
!!$       end do
!!$       SYNC(Up)
!!$       U = Up
!!$    end do

  end subroutine wilson_flow


end program main
