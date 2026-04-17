#include "inc"
module dirac
  use iso_fortran_env, only : dp => real64, i4 => int32
  use pbc, only : ip, im, sgnp, sgnm
  use constants, only : i
  implicit none
  
contains

  function D(phi,U)
    use parameters, only : Lx, Ly, m0
    complex(dp), dimension(DIM), intent(in) :: U, phi
    complex(dp), dimension(2,Lx,Ly) :: D
    integer(i4) :: x1, x2, mu
    integer(i4), dimension(2) :: x, xm1, xm2,xp1, xp2

    do x1 = 1, Lx
       do x2 = 1, Ly
          x = [x1,x2]
          xm1 = im(x,1)
          xm2 = im(x,2)
          xp1 = ip(x,1)
          xp2 = ip(x,2)
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
    use parameters, only : Lx, Ly, m0
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


end module dirac
