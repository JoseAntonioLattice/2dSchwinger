#define DIM :
#ifdef PARALLEL
#define  DIM 0:
#endif

module pbc

  use iso_fortran_env, only : i4 => int32
  use indices
  implicit none

  integer, allocatable, dimension(:) :: ip1, im1, ip2, im2, sgnp, sgnm
  integer, allocatable, dimension(:) :: ip1_c, im1_c, ip2_c, im2_c

  integer :: left, right, up, down, a(2)
  
contains

#ifndef PARALLEL
  subroutine set_pbc(L)
#else
  subroutine set_pbc(L,cores)
    integer(i4), dimension(2), intent(in) :: cores
#endif
    integer(i4), dimension(2), intent(in) :: L

    allocate(sgnp(L(1)),sgnm(L(1)))

#ifndef PARALLEL    
    allocate(ip1(L(1)), im1(L(1)))
    allocate(ip2(L(2)), im2(L(2)))
    call set_periodic_bounds(ip1,im1,L(1))
    call set_periodic_bounds(ip2,im2,L(2))
#elif defined(PARALLEL)
    call set_periodic_bounds(ip1_c,im1_c,cores(1))
    call set_periodic_bounds(ip1_c,im1_c,cores(2))

    a = get_index_array(this_image(),cores)

    left  = get_index(im_core(a,1),cores)
    right = get_index(ip_core(a,1),cores)
    up    = get_index(im_core(a,2),cores)
    down  = get_index(ip_core(a,2),cores)

#endif

    sgnp = 1
    sgnm = 1

    sgnp(L(1)) = -1
    sgnm(1) = -1


    
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
    integer(i4), dimension(DIM), intent(in) :: x
    integer(i4) :: mu
    integer(i4), dimension(size(x)) :: ip

    ip = x

#ifndef PARALLEL
    select case(mu)
    case(1)
       ip(mu) = ip1(x(mu))
    case(2)
       ip(mu) = ip2(x(mu))
    end select
#else
    ip(mu) = x(mu) + 1
#endif
  end function ip

  function im(x, mu)
    integer(i4), dimension(DIM), intent(in) :: x
    integer(i4) :: mu
    integer(i4), dimension(size(x)) :: im
    
    im = x

#ifndef PARALLEL
    select case(mu)
    case(1)
       im(mu) = im1(x(mu))
    case(2)
       im(mu) = im2(x(mu))
    end select
#else
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
  
end module pbc
