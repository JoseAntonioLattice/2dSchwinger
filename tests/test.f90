program test

  use conjugate_gradient

  implicit  none

  integer, parameter :: n = 2
  real(dp) :: A(n,n), x(n), b(n)

  A = reshape([4.0_dp,1.0_dp,1.0_dp,3.0_dp],shape(A))
  b = [1.0_dp,2.0_dp]
  x = [2.0_dp,1._dp]

  call conjg_grad(A,b,x,10E-6_dp)

  print*, x
  
end program test
