#include "inc"
module CG
  use iso_fortran_env, only : dp => real64, i4 => int32
  use arrays, only : cg_p, cg_Ddagp
  use dirac
  implicit none
contains

  function conjugate_gradient(phi,U) result(x) 
    use parameters, only : Lx, Ly, max_iter, tol
    complex(dp), dimension(DIM), intent(in) :: U, phi
    complex(dp), dimension(2,Lx,Ly) :: x
    complex(dp), dimension(2,Lx,Ly) :: r,Ap
    complex(dp) :: alpha, pAp
    real(dp)    :: phi_norm2, beta, rr, rr_new
    integer(i4) :: k

    
    cg_Ddagp(:,1:Lx,1:Ly) = Ddagger(phi,U)
    SYNC(cg_Ddagp) 
    Ap = D(cg_Ddagp,U)
    x = phi(:,1:Lx,1:Ly)
        
    r = x - Ap
    cg_p(:,1:Lx,1:Ly) = r; SYNC(cg_p)

    rr = real(sum(r*conjg(r)))
    phi_norm2 = real(sum(phi(:,1:Lx,1:Ly)*conjg(phi(:,1:Lx,1:Ly))))

#ifdef PARALLEL
    call co_sum(rr)
    call co_sum(phi_norm2)
#endif

    do k = 1, max_iter
       cg_Ddagp(:,1:Lx,1:Ly) = Ddagger(cg_p,U)
       SYNC(cg_Ddagp)
       Ap = D(cg_Ddagp,U)
       pAp = sum(cg_p(:,1:Lx,1:Ly)*conjg(Ap))
       
#ifdef PARALLEL
       call co_sum(pAp)
#endif

       alpha = rr/pAp
       x = x + alpha*cg_p(:,1:Lx,1:Ly)
       r = r - alpha*Ap
       
       rr_new = real(sum(r*conjg(r)))
#ifdef PARALLEL
       call co_sum(rr_new)
#endif
       if( rr_new < tol*phi_norm2 ) exit
       beta = rr_new/rr
       cg_p(:,1:Lx,1:Ly) = r + beta*cg_p(:,1:Lx,1:Ly)
       SYNC(cg_p)
       rr = rr_new
    end do
    
  end function conjugate_gradient

  
end module CG
