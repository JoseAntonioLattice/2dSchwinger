#include "inc"
module Gradient_flow
  use iso_fortran_env, only : dp => real64, i4 => int32
  use arrays, only : Up
  use u1
  implicit none
  
contains
  
  subroutine wilson_flow(U,beta)
    use parameters, only : Lx, Ly
    complex(dp), intent(inout) :: U(DIM)
    real(dp), intent(in) :: beta
    real(dp), parameter :: epsilon = 1.0E-3_dp
    integer, parameter :: Nt = 100
    integer :: x,y,t
    
    do t = 1, Nt
       do x = 1, Lx
          do y = 1, Ly
             Up(1,x,y) = exp(-epsilon*beta*real(U(1,x,y)*conjg(staples(U,[x,y],1))))*U(1,x,y)
             Up(2,x,y) = exp(-epsilon*beta*real(U(2,x,y)*conjg(staples(U,[x,y],2))))*U(2,x,y)
          end do
       end do
       !SYNC(Up)
       U = Up
    end do

  end subroutine wilson_flow

  
end module Gradient_flow
