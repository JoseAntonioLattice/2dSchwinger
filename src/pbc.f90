module pbc

  use iso_fortran_env, only : i4 => int32
  
  implicit none

  integer, allocatable, dimension(:) :: ip, im, sgnp,sgnm

contains

  subroutine set_pbc(L)
    
    integer(i4), intent(in) :: L

    allocate(ip(L), im(L))
    allocate(sgnp(L),sgnm(L))
    call initialize(L)
    
  end subroutine set_pbc

  subroutine initialize(L)

    integer(i4), intent(in) :: L
    integer(i4) :: i

    do i = 1, L
       ip(i) = i + 1
       im(i) = i - 1
       
    end do
    ip(L) = 1
    im(1) = L

    sgnp = 1
    sgnm = 1

    sgnp(L) = -1
    sgnm(1) = -1
    
  end subroutine initialize

  function ipf(vector, mu)

    integer(i4), dimension(:), intent(in) :: vector
    integer(i4) :: mu

    integer(i4), dimension(size(vector)) :: ipf
    
    ipf = vector

    ipf(mu) = ip(vector(mu))
    
  end function ipf

  function imf(vector, mu)

    integer(i4), dimension(:), intent(in) :: vector
    integer(i4) :: mu

    integer(i4), dimension(size(vector)) :: imf
    
    imf = vector

    imf(mu) = im(vector(mu))
    
  end function imf

  
end module pbc
