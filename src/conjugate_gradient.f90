module conjugate_gradient

  use iso_fortran_env, only : dp => real64
  implicit none
  
contains

  subroutine conjg_grad(A,b,x,tol)
    real(dp), dimension(:,:), intent(in) :: A
    real(dp), dimension(:), intent(in) :: b
    real(dp), dimension(:), intent(inout) :: x
    real(dp), intent(in) :: tol
    
    real(dp) :: err, alpha, beta
    real(dp), dimension(size(x)) :: r, r0, d, m

    integer :: iter
    
    r = b - matmul(A,x)
    d = r
    iter = 0
    err = dot_product(r,r)
    print*, iter, err
    do while(err >= tol)
       m = matmul(A,d)
       alpha = dot_product(r,r)/dot_product(d,m)
       r0 = r
       x = x + alpha*d
       r = r - alpha*m
       err = dot_product(r,r)
       
       beta = err/dot_product(r0,r0)
       
       d = r + beta*d
       iter = iter + 1
       print*, iter, err
       if(iter >= 100) stop "Method did not converge"
    end do
  end subroutine conjg_grad

  
end module conjugate_gradient
