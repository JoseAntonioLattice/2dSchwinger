#ifdef PARALLEL
#define DIM 2,0:Lx+1,0:Ly+1
#define CODIM [*]
#define CODIM2 ,codimension[:]
#define CODIM3 ,codimension[*]
#define ALLOC ,allocatable
#define SYNC(A) sync all; \
  A(:,  Lx+1,1:Ly  ) = A(:,1    ,1:Ly)[right]; \
  A(:,0     ,1:Ly  ) = A(:,   Lx,1:Ly)[left]; \
  A(:,1:Lx  ,  Ly+1) = A(:, 1:Lx,1   )[up]; \
  A(:,1:Lx  ,0     ) = A(:, 1:Lx,  Ly)[down]; \
  sync all; \
  A(:,Lx+1, Ly+1) = A(:, 1 , 1 )[right_up]; \
  A(:,Lx+1, 0   ) = A(:, 1 , Ly)[right_down]; \
  A(:,0   , Ly+1) = A(:, Lx, 1 )[left_up]; \
  A(:,0   , 0   ) = A(:, Lx, Ly)[left_down]; \
  sync all
#else
#define DIM 2,Lx,Ly
#define CODIM
#define CODIM2
#define CODIM3
#define ALLOC
#define SYNC(A)
#endif

program test
  use iso_fortran_env, only: dp => real64, i4 => int32
  complex(dp), dimension(:,:,:), allocatable CODIM2 :: U, chi, phi, pp
  integer :: Lx, Ly
  integer, dimension(2), parameter :: L = [16,16]
  real(dp) :: m0 = 10.0_dp
  integer, allocatable, dimension(:) :: ip1, im1, ip2, im2, sgnp, sgnm
  complex(dp), dimension(:,:,:), allocatable CODIM2 :: U_global, chi_global
  complex(dp), parameter :: i = (0.0_dp,1.0_dp)
  
#ifdef PARALLEL
  integer :: cores(2)
  integer, allocatable, dimension(:) :: ip1_c, im1_c, ip2_c, im2_c
  integer :: left, right, up, down, left_down, left_up, right_down, right_up
#endif

  
  Lx = L(1)
  Ly = L(2)
  
#ifdef PARALLEL
  if(this_image() == 1) then
     print*, "Enter cores array: "
     read(*,*) cores
     print*, "user typed # cores: ", cores
     
  end if
  call co_broadcast(cores,source_image=1)
    
  Lx = Lx/cores(1)
  Ly = Ly/cores(2)
  call set_pbc([Lx,Ly],cores)
  sync all
#else
  call set_pbc([Lx,Ly])
#endif

  print*, Lx, Ly
  allocate(U(DIM)CODIM)
  allocate(chi(DIM)CODIM)
  allocate(phi(DIM)CODIM)
  
  ALLOCATE(U_global(2,L(1),L(2))CODIM)
  ALLOCATE(chi_global(2,L(1),L(2))CODIM)
  allocate(pp(2,Lx,Ly)CODIM)
  call read_fields(U,chi)
  sync all
  call save_configuration(U,1.0_dp,1)
  SYNC(U)
  SYNC(chi)
  !call check_Dirac(U,chi)

  phi(:,1:Lx,1:Ly) = D(chi,U); SYNC(phi)
 
  pp = conjugate_gradient(phi,U)
  print*, "CG done!"
#ifdef PARALLEL
  sync all
  if(this_image() == 1) then
#endif
  print*, "phi", phi(1,1,1)
  print*, "(DDdag)^-1", pp(1,1,1)
#ifdef PARALLEL
  endif
#endif  
contains

  subroutine check_Dirac(U,chi)
    complex(dp), dimension(DIM), intent(in) CODIM3 :: U, chi
    integer :: x, y
    
    phi(:,1:Lx,1:Ly) = D(chi,U); SYNC(phi)

#ifdef PARALLEL
    if(this_image()==1)then
#endif
    do x = 1, Lx
       do y = 1, Ly
          print"(a,i3,',',x,a,i3,',',x,a,'(',f8.4,',',f8.4,')')", "x=",x,"y=",y,"Dchi=", phi(1,x,y)
       end do
    end do
#ifdef PARALLEL
    endif    
#endif
  end subroutine check_Dirac

  
  subroutine read_fields(U,chi)
    complex(dp), dimension(DIM), intent(out) CODIM3 :: U, chi
    integer :: ix, ex, iy, ey, a(2)
    
#ifdef PARALLEL
    sync all
    
    if(this_image() == 1) then
#endif
       open(unit = 69, file = "u.dat"  ,status="old")
       open(unit = 71, file = "chi.dat",status="old")
       read(69,*) U_global
       read(71,*) chi_global
       print*, U_global(1,1,1)      
#ifdef PARALLEL
    end if
    call co_broadcast(U_global,source_image=1)
    call co_broadcast(chi_global,source_image=1)
    sync all
    a = get_index_array(this_image(),cores)
    ix = L(1)/cores(1)*(a(1)-1)+1
    ex = L(1)/cores(1)*a(1)
    iy = L(2)/cores(2)*(a(2)-1)+1
    ey = L(2)/cores(2)*a(2)

    print"('image: ',i0,',',x,'x in','[',i3,',',i3,']',3x,'y in','[',i3,',',i3,']')", this_image(), ix,ex,iy,ey
    
      U(:,1:Lx,1:Ly) =   U_global(:,ix:ex,iy:ey)!; SYNC(U)
    chi(:,1:Lx,1:Ly) = chi_global(:,ix:ex,iy:ey)!; SYNC(chi)
    
#else
    U = U_global
    chi = chi_global
#endif
  end subroutine read_fields

  
#ifndef PARALLEL
  subroutine set_pbc(L)
#elif defined(PARALLEL)
  subroutine set_pbc(L,cores)
    integer(i4), dimension(2), intent(in) :: cores
#endif
    integer(i4), dimension(2), intent(in) :: L
    integer :: a(2)

    allocate(sgnp(L(1)),sgnm(L(1)))
#ifndef PARALLEL    
    allocate(ip1(L(1)), im1(L(1)))
    allocate(ip2(L(2)), im2(L(2)))
    call set_periodic_bounds(ip1,im1,L(1))
    call set_periodic_bounds(ip2,im2,L(2))
#elif defined(PARALLEL)
    print*, this_image(), "cores", cores
    print*, this_image(), "inside set_pbc()"
    allocate(ip1_c(cores(1)), im1_c(cores(1)))
    allocate(ip2_c(cores(2)), im2_c(cores(2)))
    call set_periodic_bounds(ip1_c,im1_c,cores(1))
    call set_periodic_bounds(ip2_c,im2_c,cores(2))

    a = get_index_array(this_image(),cores)

    left  = get_index(im_core(a,1),cores)
    right = get_index(ip_core(a,1),cores)
    up    = get_index(ip_core(a,2),cores)
    down  = get_index(im_core(a,2),cores)

    left_down  = get_index(im_core(im_core(a,1),2),cores)
    left_up    = get_index(ip_core(im_core(a,1),2),cores)
    right_down = get_index(im_core(ip_core(a,1),2),cores)
    right_up   = get_index(ip_core(ip_core(a,1),2),cores)

#endif
    sgnp = 1
    sgnm = 1
#if defined(PARALLEL)
    if( a(1) == cores(1) ) then
#endif
       sgnp(L(1)) = -1
#if defined(PARALLEL)
    endif
    if( a(1) == 1 ) then
#endif       
       sgnm(1) = -1
#if defined(PARALLEL)       
    endif
    print"(a,i2,',',x,*(i2,x))", "image: ",this_image(),sgnp
    print"(a,i2,',',x,*(i2,x))", "image: ",this_image(),sgnm
#endif
  end subroutine set_pbc

  subroutine set_periodic_bounds(ip_array,im_array,L)
    integer(i4), intent(in) :: L
    integer, dimension(L), intent(out) :: ip_array, im_array
    integer(i4) :: i

    do i = 1, L
       ip_array(i) = i + 1
       im_array(i) = i - 1
    end do
    ip_array(L) = 1
    im_array(1) = L
    
  end subroutine set_periodic_bounds
  
  function ip(x, mu)
    integer(i4), dimension(2), intent(in) :: x
    integer(i4) :: mu
    integer(i4), dimension(2) :: ip

    ip = x

#ifndef PARALLEL
    select case(mu)
    case(1)
       ip(mu) = ip1(x(mu))
    case(2)
       ip(mu) = ip2(x(mu))
    end select
#elif defined(PARALLEL)
    ip(mu) = x(mu) + 1
#endif
  end function ip

  function im(x, mu)
    integer(i4), dimension(2), intent(in) :: x
    integer(i4) :: mu
    integer(i4), dimension(2) :: im
    
    im = x

#ifndef PARALLEL
    select case(mu)
    case(1)
       im(mu) = im1(x(mu))
    case(2)
       im(mu) = im2(x(mu))
    end select
#elif defined(PARALLEL)
    im(mu) = x(mu) - 1
#endif
  end function im


  function im_core(x,mu)
    integer(i4) :: im_core(2)
    integer(i4), intent(in) :: x(2), mu
    
    im_core = x
    
    select case(mu)
    case(1)
       im_core(mu) = im1_c(x(mu))
    case(2)
       im_core(mu) = im2_c(x(mu))
    end select
    
  end function im_core
  
  function ip_core(x,mu)
    integer(i4) :: ip_core(2)
    integer(i4), intent(in) :: x(2), mu
    
    ip_core = x
    
    select case(mu)
    case(1)
       ip_core(mu) = ip1_c(x(mu))
    case(2)
       ip_core(mu) = ip2_c(x(mu))
    end select
    
  end function ip_core

  function get_index(vector,length_lattice)

    integer :: get_index

    integer, dimension(:), intent(in) :: vector
    integer :: dimension_space
    integer, intent(in) :: length_lattice(:)

    integer :: suma, prod
    integer :: i

    dimension_space = size(length_lattice)

    suma = vector(1)
    prod = 1
    if( dimension_space > 1)then
       do i = 2, dimension_space
          suma = suma + (vector(i) - 1) * product(length_lattice(1:i-1))
       end do
    end if

    get_index = suma


  end function get_index


  function get_index_array(idx,L) result(vector)
    integer, intent(in) :: idx
    integer :: d
    integer, intent(in) :: L(:)

    integer, dimension(size(L)) :: vector

    integer :: i, n, modx, suma, prod1, prod2

    d = size(L)
    
    modx = mod(idx,L(1))
    vector(1) = modx
    if(modx == 0) vector(1) = L(1)

    suma = vector(1)
    do i = 2, d
      modx = mod(idx,product(L(1:i)))
      if (i > 2) suma = suma + product(L(1:i-2))*(vector(i-1)-1)
      vector(i) = (modx - suma)/product(L(1:i-1)) + 1
      if(modx == 0) vector(i) = L(i)
    end do

  end function get_index_array

  function D(phi,U)
    complex(dp), dimension(DIM), intent(in) :: U, phi
    complex(dp), dimension(2,Lx,Ly) :: D
    integer(i4) :: x1, x2, mu
    integer(i4), dimension(2) :: x, xm1, xm2,xp1, xp2

    do x1 = 1, Lx
       do x2 = 1, Ly
                        ! Including periodic boundary conditions
          x = [x1,x2]   ! x = (x_1,x_2) 
          xm1 = im(x,1) ! (x_1 - 1,    x_2)
          xm2 = im(x,2) ! (x_1,    x_2 - 1)
          xp1 = ip(x,1) ! (x_1 + 1,    x_2)
          xp2 = ip(x,2) ! (x1,     x_2 + 1)
          D(1,x(1),x(2)) = (m0+2.0_dp) * phi(1,x(1),x(2)) - 0.5*( &
               U(1,x(1),x(2))*sgnp(x1) *(phi(1,xp1(1),xp1(2)) -  phi(2,xp1(1),xp1(2))) + &
               U(2,x(1),x(2))          *(phi(1,xp2(1),xp2(2)) +i*phi(2,xp2(1),xp2(2))) + &
               conjg(U(1,xm1(1),xm1(2)))*sgnm(x1)*(phi(1,xm1(1),xm1(2)) +  phi(2,xm1(1),xm1(2))) + &
               conjg(U(2,xm2(1),xm2(2)))         *(phi(1,xm2(1),xm2(2)) -i*phi(2,xm2(1),xm2(2))) &
               )
          D(2,x(1),x(2)) = (m0+2.0_dp) * phi(2,x(1),x(2)) - 0.5*( &
               U(1,x(1),x(2))*sgnp(x1) *(  -phi(1,xp1(1),xp1(2)) + phi(2,xp1(1),xp1(2))) + &
               U(2,x(1),x(2))          *(-i*phi(1,xp2(1),xp2(2)) + phi(2,xp2(1),xp2(2))) + &
               conjg(U(1,xm1(1),xm1(2)))*sgnm(x1)*(  phi(1,xm1(1),xm1(2)) + phi(2,xm1(1),xm1(2))) + &
               conjg(U(2,xm2(1),xm2(2)))         *(i*phi(1,xm2(1),xm2(2)) + phi(2,xm2(1),xm2(2))) &
               )
       end do
    end do
    
  end function D

    function Ddagger(phi,U)
    complex(dp), dimension(DIM), intent(in) :: U, phi
    complex(dp), dimension(2,Lx,Ly) :: Ddagger
    integer(i4) :: x(2), xm1(2), xm2(2),xp1(2), xp2(2)
    integer(i4) :: x1, x2
    
    do x1 = 1, Lx
       do x2 = 1, Ly
          x = [x1,x2]
          xm1 = im(x,1)
          xm2 = im(x,2)
          xp1 = ip(x,1)
          xp2 = ip(x,2)
          Ddagger(1,x(1),x(2)) = (m0+2.0_dp) * phi(1,x(1),x(2)) - 0.5*( &
               conjg(U(1,xm1(1),xm1(2)))*sgnm(x1)*(phi(1,xm1(1),xm1(2)) -  phi(2,xm1(1),xm1(2))) + &
               conjg(U(2,xm2(1),xm2(2)))         *(phi(1,xm2(1),xm2(2)) +i*phi(2,xm2(1),xm2(2))) + &
               U(1,x(1),x(2))*sgnp(x1)*(phi(1,xp1(1),xp1(2)) +  phi(2,xp1(1),xp1(2))) + &
               U(2,x(1),x(2))         *(phi(1,xp2(1),xp2(2)) -i*phi(2,xp2(1),xp2(2))) &
               )
          Ddagger(2,x(1),x(2)) = (m0+2.0_dp) * phi(2,x(1),x(2)) - 0.5*( &
               conjg(U(1,xm1(1),xm1(2)))*sgnm(x1)*(  -phi(1,xm1(1),xm1(2)) + phi(2,xm1(1),xm1(2))) + &
               conjg(U(2,xm2(1),xm2(2)))         *(-i*phi(1,xm2(1),xm2(2)) + phi(2,xm2(1),xm2(2))) + &
               U(1,x(1),x(2))*sgnp(x1)*(  phi(1,xp1(1),xp1(2)) + phi(2,xp1(1),xp1(2))) + &
               U(2,x(1),x(2))         *(i*phi(1,xp2(1),xp2(2)) + phi(2,xp2(1),xp2(2))) &
               )
       end do
    end do
        
  end function Ddagger

  function conjugate_gradient(phi,U) result(x) 
    integer, parameter :: max_iter = 15
    real(dp), parameter :: tol = 1.0E-6_dp
    complex(dp), dimension(DIM), intent(in) :: U, phi
    complex(dp), dimension(2,Lx,Ly) :: x
    complex(dp), dimension(2,Lx,Ly) :: r,Ap
    complex(dp), allocatable, dimension(:,:,:) CODIM2 :: p, Ddagp
    complex(dp) :: alpha, pAp
    real(dp)    :: phi_norm2, beta, rr, rr_new, pAp_r, pAp_i
    integer(i4) :: k

    allocate(p(DIM)CODIM)
    allocate(Ddagp(DIM)CODIM)
        
    x = phi(:,1:Lx,1:Ly)
    Ddagp(:,1:Lx,1:Ly) = Ddagger(phi,U); SYNC(Ddagp) 
    Ap = D(Ddagp,U)
    r = x - Ap
#ifdef PARALLEL
  if(this_image() == 1) then
#endif
     print*, "Inside CG phi",x(1,1,1)
#ifdef PARALLEL
  endif
#endif
    p(:,1:Lx,1:Ly) = r; SYNC(p)

    rr = real(sum(r*conjg(r)),dp)
    phi_norm2 = real(sum(phi(:,1:Lx,1:Ly)*conjg(phi(:,1:Lx,1:Ly))),dp)

#ifdef PARALLEL
    call co_sum(rr)
    call co_sum(phi_norm2)
    if(this_image() == 1) then
#endif
       !print*,"r", r
       print*, 0, "rr", rr
#ifdef PARALLEL
    end if
#endif
    do k = 1, max_iter
       Ddagp(:,1:Lx,1:Ly) = Ddagger(p,U); SYNC(Ddagp) 
       Ap = D(Ddagp,U)
       pAp_r = sum(real (p(:,1:Lx,1:Ly)*conjg(Ap)))
       pAp_i = sum(aimag(p(:,1:Lx,1:Ly)*conjg(Ap)))
#ifdef PARALLEL
       call co_sum(pAp_r)
       call co_sum(pAp_i)
#endif
       pAp = cmplx(pAp_r,pAp_i)
       alpha = rr/pAp
       
       x = x + alpha*p(:,1:Lx,1:Ly)
       r = r - alpha*Ap

       rr_new = real(sum(r*conjg(r)))
#ifdef PARALLEL
       call co_sum(rr_new)
       if(this_image()==1)then
#endif
          print*, k,"rr_new",rr_new,"phi_norm2",phi_norm2,"x", x(1,1,1)
#ifdef PARALLEL
       end if
       sync all
#endif
          
       if( rr_new < tol*phi_norm2 ) exit
       beta = rr_new/rr
       p(:,1:Lx,1:Ly) = r + beta*p(:,1:Lx,1:Ly); SYNC(p)
       rr = rr_new
   
    end do
  end function conjugate_gradient

  subroutine save_configuration(U,beta,n)
    use number2string
    !use parameters, only : Lx, Ly, L
    complex(dp), intent(in) :: U(DIM)
    real(dp), intent(in) :: beta
    integer(i4), intent(in) :: n
    integer(i4) :: un, ix,ex, iy, ey, a(2)
    complex(dp), allocatable :: U_global(:,:,:)[:]
    
    
#ifndef PARALLEL
    allocate(U_global(2,L(1),L(2)))
    U_global = U
#else
   
    allocate(U_global(2,L(1),L(2))[*])
  
    a = get_index_array(this_image(),cores)
    ix = L(1)/cores(1)*(a(1)-1)+1
    ex = L(1)/cores(1)*a(1)
    iy = L(2)/cores(2)*(a(2)-1)+1
    ey = L(2)/cores(2)*a(2)
    sync all
    U_global(:,ix:ex,iy:ey)[1] = U(:,1:Lx,1:Ly)
    sync all
    if( this_image()==1 ) then
#endif
       !open(newunit = un, file = "U"//int2str(n)//".bin", access = "sequential", form = "unformatted")
       open(newunit = un, file = "U"//int2str(n)//".dat")
       write(un,*) U_global
       close(un)
#ifdef PARALLEL
    endif
#endif
    
  end subroutine save_configuration


  
end program test
