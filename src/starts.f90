#include "inc"
module starts
  use iso_fortran_env, only : dp => real64
  use constants, only : i, pi
  implicit none
  
contains
  
  subroutine hot_start(u)
    use parameters, only : Lx, Ly    
    complex(dp), intent(out), dimension(DIM) CODIM :: u
    real(dp), dimension(2,Lx,Ly) :: phi
    
    call random_number(phi)

    phi = 2*pi*phi
    u(:,1:Lx,1:Ly) = exp(i*phi); SYNC(u)
  end subroutine hot_start

  subroutine cold_start(u)
    complex(dp), intent(out), dimension(:,:,:) :: u
    u = 1.0_dp
  end subroutine cold_start

  
end module starts
