#include "inc"
module observables
  use iso_fortran_env, only : dp => real64, i4 => int32
  use u1
  use constants, only :  pi
  implicit none

contains

  function DS(uold, unew, beta,stp)
    real(dp) :: DS,beta
    complex(dp) :: uold, unew, stp
    
    DS = -beta * real( (unew - uold) * conjg(stp) )
    
  end function DS
  

  function action(u)
    use parameters, only : Lx, Ly
    real(dp) :: action
    complex(dp), dimension(DIM), intent(in) :: u
    integer(i4) :: x,y

    action = 0.0_dp
    
    do x = 1, Lx
       do y = 1, Ly
          action = action + real(plaquette(u,[x,y]))
       end do
    end do

  end function action

  
  function action2(u,beta)
    use parameters, only : Lx, Ly
    real(dp) :: action2
    complex(dp), dimension(DIM), intent(in) :: u
    real(dp), intent(in) :: beta
    integer(i4) :: x,y

    action2 = 0.0_dp
    
    do x = 1, Lx
       do y = 1, Ly
          action2 = action2 + beta*(1.0_dp - real(plaquette(u,[x,y])))
       end do
    end do
    
  end function action2

  subroutine topological_charge_density(u,top)
    use parameters, only : Lx, Ly
    complex(dp), dimension(DIM), intent(in) :: u
    real(dp), intent(out) :: top(Lx,Ly) 
    integer(i4) :: x, y
    complex(dp) :: plq
    
    do x = 1, Lx
       do y = 1, Ly
          plq = plaquette(u,[x,y])
          top(x,y) = atan2(plq%im,plq%re)!/(2*pi)
       end do
    end do
       
  end subroutine topological_charge_density

  function slab_top_char(top,ix)
    real(dp), intent(in) :: top(:,:)
    integer(i4), intent(in) :: ix
    real(dp) :: slab_top_char
    complex(dp) :: plq

    slab_top_char = (sum(top(1:ix,:))/(2*pi))**2!*sum(top(ix+1:,:))
    
  end function slab_top_char


  function correlation_polyakov(U) result(corr_poly)
    use parameters, only : L
    complex(dp), intent(in) :: U(2,L(1),L(2))
    integer(i4) :: x, r, t
    complex(dp), dimension(0:L(1)/2-1) :: corr_poly
    complex(dp), dimension(L(1)) :: polyakov_loop

    do x = 1, L(1)
       polyakov_loop(x) = product(U(2,x,:))
    end do

    corr_poly = 0.0_dp
    do t = 0, L(1)/2 - 1
       do x = 1, L(1)
          r = mod(x-1+t,L(1))+1
          corr_poly(t) = corr_poly(t)+polyakov_loop(x)*conjg(polyakov_loop(r))
       end do
    end do
    corr_poly = corr_poly/L(1)
  end function correlation_polyakov

end module observables
