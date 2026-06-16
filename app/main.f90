program main
  !use format
  use number2string
  use pbc
  use iso_fortran_env, only : dp => real64, i4 => int32
  implicit none

  integer, parameter :: Lx = 16
  integer, parameter :: Lt = 16
  integer, parameter :: Nt = 100
  integer, parameter :: top_sector = 0

  real(dp), dimension(Lx,Lt) :: Q

  integer, parameter :: nconfig = 1000
  real(dp), parameter :: pi = acos(-1.0_dp)

  real(dp) :: beta = 3.0_dp
  real(dp) :: m0 = 100.0_dp
  real(dp) :: epsilon = 1.0E-2_dp
  
  complex(dp), dimension(2,Lx,Lt) :: U, V
  character(:), allocatable :: filename

  integer :: i, un, count, ix, it, t

  real(dp) :: a_action(nconfig), slb(nconfig,Lx,0:Nt), E(nconfig,0:Nt)

  INTEGER :: Q0, Qf
  

  
  call set_pbc([Lx,Lt])

  
  
  open(newunit = un, file = "energy_Q0.dat")
  open(unit = 69, file = "q0.dat")
  open(unit = 42, file = "slab_Q0.dat")
  count = 0
  slb = 0.0_dp
  do i = 1, nconfig
     filename = "data/configurations/m0="//real2str(m0,3,4)// &
          "/Lx="//int2str(Lx)//"/Lt="//int2str(Lt)//&
          "/beta="//real2str(beta,1,4)// &
          "/U_"//int2str(i)//".bin"
     call read_configuration(U,filename)

     V = U
     Q0 = nint(abs(top_char(U)))
     if(Q0 == top_Sector)then
        !count = count + 1
        call wilson_flow_rk4(V,beta,Nt,epsilon,.false.)
        Qf = nint(abs(top_char(V)))
        
        if( Q0 == Qf ) then
           count = count + 1
           Print*, count, i, Q0, Qf, Q0==Qf
           call wilson_flow_rk4(U,beta,Nt,epsilon,.true.)
        end if
     END if
        

        
  end do


  
  do t = 0, Nt
     do ix = 1, Lx
        write(42,*) 1.0_dp*ix/(Lx), avr(slb(1:count,ix,t)), std_err(slb(1:count,ix,t))
     end do
     write(42,"(2/)")
     print*, t*epsilon, avr(E(1:count,t)),std_err(E(1:count,t))
  end do
  print*, count
  
contains

  function avr(x)
    real(dp), intent(in) :: x(:)
    real(dp) :: avr

    avr = sum(x)/size(x)
    
  end function avr

  function var(x)
    real(dp), intent(in) :: x(:)
    real(dp) :: var

    var = sum((x - avr(x))**2)/(size(x)-1) 
    
  end function var

  function std_err(x)
    real(dp), intent(in) :: x(:)
    real(dp) :: std_err

    std_err = sqrt(var(x)/size(x))
    
  end function std_err
  
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
    !use pbc, only : ip, im
    complex(dp) :: plaquette
    complex(dp), dimension(:,:,:), intent(in) :: u
    integer(i4), intent(in) :: x(2)
    integer(i4), dimension(2) :: x2, x3


    x2 = ip(x,1)
    x3 = ip(x,2)
    
    plaquette = U(1,x(1),x(2)) * U(2,x2(1),x2(2)) * &
          conjg(U(1,x3(1),x3(2))) * conjg(U(2,x(1),x(2)))
    
  end function plaquette

  function staples(u,x,mu)
    complex(dp) :: staples
    complex(dp), dimension(:,:,:), intent(in) :: u
    integer(i4), intent(in) :: x(2), mu
    integer(i4), dimension(2) :: x2, x3, x4, x5, x6

    integer(i4) :: nu
    
    if ( mu == 1 ) then
       nu = 2
    elseif( mu == 2)then
       nu = 1
    end if

    x2 = ip(x,nu)
    x3 = ip(x,mu)
    x4 = im(x,nu)
    x5 = x4
    x6 = im(x3,nu)
    
    staples = u(nu,x(1),x(2)) * u(mu,x2(1), x2(2)) * conjg( u(nu,x3(1), x3(2)) ) + &
         conjg( u(nu,x4(1),x4(2)) ) * u(mu,x5(1), x5(2)) * u(nu,x6(1), x6(2))
    
  end function staples
  
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

    function energy(u)
    real(dp) :: energy
    complex(dp), dimension(:,:,:), intent(in) :: u
    integer(i4) :: x,y, Lx, Ly

    Lx = size(u(1,:,1))
    Ly = size(u(1,1,:))
    
    energy = 0.0_dp
    
    do x = 1, size(u(1,:,1))
       do y = 1,  size(u(1,1,:))
          energy = energy + real(plaquette(u,[x,y]))
       end do
    end do

    energy = Lx*Ly-energy
    
  end function energy
  
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
    complex(dp), dimension(:,:,:), intent(in) :: u
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

  function top_char(U)
    complex(dp), dimension(:,:,:), intent(in) :: u
    integer(i4) :: i, j, Lx, Ly
    complex(dp) :: plq
    real(dp) :: top_char
    
    Lx = size(u(1,:,1))
    Ly = size(u(1,1,:))
    top_char = 0.0_dp
    do i = 1, Lx
       do j = 1, Ly
          plq = plaquette(u,[i,j])
          top_char = top_char + atan2(plq%im,plq%re)
       end do
    end do
    top_char = top_char/(2*pi)
  end function top_char

  function top_char_arr(U)
    complex(dp), dimension(:,:,:), intent(in) :: u
    integer(i4) :: i, j, Lx, Ly
    complex(dp) :: plq
    real(dp) :: top_char_arr(size(U(1,:,1)),size(U(1,1,:)))
    
    Lx = size(u(1,:,1))
    Ly = size(u(1,1,:))
  
    do i = 1, Lx
       do j = 1, Ly
          plq = plaquette(u,[i,j])
          top_char_arr(i,j) =  atan2(plq%im,plq%re)!/(2*pi)
       end do
    end do
    !top_char = top_char!/(2*pi)
  end function top_char_arr


  function slab_top_char(top,ix)
    real(dp), intent(in) :: top(:,:)
    integer(i4), intent(in) :: ix
    real(dp) :: slab_top_char
    integer(i4) :: i, j
    complex(dp) :: plq
    real(dp), parameter :: pi = acos(-1.0_dp)

    slab_top_char = (sum(top(1:ix,:))/(2*pi))**2!*sum(top(ix+1:,:))
    
  end function slab_top_char

 
  subroutine wilson_flow(U,beta,Nt,epsilon) !result(Up)
    complex(dp), intent(inout) :: U(:,:,:)
    complex(dp) :: Up(2,size(U(1,:,1)),size(U(1,1,:)))
    real(dp), intent(in) :: beta
    real(dp), intent(in) :: epsilon
    integer, intent(in) :: Nt 
    integer :: x,y,t, Lx, Ly
    !real(dp) :: slb(size(u(1,:,1)))

    Lx = size(U(1,:,1))
    Ly = size(U(1,1,:))

    Q = top_char_arr(U)
    E(count,0) = energy(U)
    
    write(un,*) 0.0_dp, energy(U), sum(top_char_arr(U))/(2*pi)
    do ix = 1, Lx
       slb(count,ix,0) = slab_top_char(Q,ix)
    end do
    do t = 1, Nt
       do x = 1, Lx
          do y = 1, Ly
             Up(1,x,y) = exp( cmplx(0.0_dp, -epsilon*beta*aimag(U(1,x,y)*conjg(staples(U,[x,y],1))), dp) ) * U(1,x,y)
             Up(2,x,y) = exp( cmplx(0.0_dp, -epsilon*beta*aimag(U(2,x,y)*conjg(staples(U,[x,y],2))), dp) ) * U(2,x,y)
          end do
       end do
       U = Up
       Q = top_char_arr(U)
       
       E(count,t) = energy(U)
       write(un,*) t*epsilon, energy(U), sum(top_char_arr(U))/(2*pi)
       do ix = 1, Lx
          slb(count,ix,t) = slab_top_char(Q,ix)
       end do
    end do
    write(un,"(2/)") 

  end subroutine wilson_flow



  subroutine wilson_flow_rk4(U,beta,Nt,epsilon,condition) !result(Up)
    complex(dp), intent(inout) :: U(:,:,:)
    complex(dp), dimension(2,size(U(1,:,1)),size(U(1,1,:))) :: W0, W1, W2, W3, Z0, Z1, Z2, B
    real(dp), intent(in) :: beta
    real(dp), intent(in) :: epsilon
    integer, intent(in) :: Nt
    logical,intent(in) :: condition
    integer :: x,y,t, Lx, Ly
    
    !real(dp) :: slb(size(u(1,:,1)))

    Lx = size(U(1,:,1))
    Ly = size(U(1,1,:))

    
    if(condition) then
       Q = top_char_arr(U)
       E(count,0) = energy(U)
       write(un,*) 0.0_dp, energy(U), sum(Q)/(2*pi)
       do ix = 1, Lx
          slb(count,ix,0) = slab_top_char(Q,ix)
          write(69,*) 0.0_dp, i, slb(count,ix,0)
       end do
    end if
    do t = 1, Nt

       W0 = U
       do x = 1, Lx
          do y = 1, Ly
             Z0(1,x,y) = epsilon*Zeta(W0,[x,y],1,beta)
             Z0(2,x,y) = epsilon*Zeta(W0,[x,y],2,beta)
             
          end do
       end do

       W1 = exp(0.25_dp*Z0)*W0

       do x = 1, Lx
          do y = 1, Ly
             Z1(1,x,y) = epsilon*Zeta(W1,[x,y],1,beta)
             Z1(2,x,y) = epsilon*Zeta(W1,[x,y],2,beta)
          end do
       end do

       B = 8.0_dp/9.0_dp*Z1-17.0_dp/36.0_dp*Z0
       W2 = exp(B)*W1

       do x = 1, Lx
          do y = 1, Ly
             Z2(1,x,y) = epsilon*Zeta(W2,[x,y],1,beta)
             Z2(2,x,y) = epsilon*Zeta(W2,[x,y],2,beta)
          end do
       end do

       W3 = exp(0.75_dp*Z2)*W1
       
       U = W3
       
       if(condition) then
          Q = top_char_arr(U)
          E(count,t) = energy(U)
          write(un,*) t*epsilon, E(count,t), sum(Q)/(2*pi)
          !write(*,*) t*epsilon, E(count,t), sum(Q)/(2*pi)
          do ix = 1, Lx
             slb(count,ix,t) = slab_top_char(Q,ix)
             write(69,*) i, slb(count,ix,t)
          end do
       end if
       
       
    end do
    if(condition) write(un,"(2/)") 

  end subroutine wilson_flow_rk4

  function Zeta(U,x,mu,beta)
    complex(dp) :: Zeta
    complex(dp), intent(in) :: U(:,:,:)
    integer, intent(in) :: x(2), mu
    real(dp) :: beta
    Zeta = cmplx( 0.0_dp, -beta*aimag(U(mu,x(1),x(2))*conjg(staples(U,x,mu))), dp)
  end function Zeta

end program main
