program mass_eff

  implicit none
  integer :: un
  integer, parameter :: L = 16
  real(8), dimension(0:L) :: c, err_c
  integer :: i,id
  real(8) :: x0, pi_mass(L)
  
  open(newunit = un, file = "pion_correlator.dat", status = "old", action = "read")
  do i = 0, L
     read(un,*) id, c(i), err_c(i)
  end do
  close(un)
  open(newunit = un, file = "effective_mass_m=-0.05.dat")

  x0 = 1.2d0
  do i = 0, L-1
     call newton_rhapson(i,c,L,x0)
     print*, i, x0
     write(un,*) i, x0
     pi_mass(i+1) = x0
  end do
  write(un,"(2/)")
  print*, avr(pi_mass(4:11)), std_err(pi_mass(4:11))
  write(un,*) avr(pi_mass(4:11)), std_err(pi_mass(4:11))

  

contains

  subroutine newton_rhapson(i,c,L,x0)
   ! procedure(fun) :: f, df
    integer, intent(in) :: i, L
    real(8), intent(in) :: c(:)
    real(8), intent(inout) :: x0
    real(8) :: err, xold 
    real(8), parameter :: tol = 0.001d0
    integer :: k
    
    err = 10.0d0
    k = 0
    
    do while (k < 100 )
       xold = x0
       x0 = xold - f(i,c,L,xold)/df(i,c,L,xold)
       err = abs(x0-xold)
       
       k = k + 1
       !print*, "iteration", k, err
    end do
    
  end subroutine newton_rhapson
  
  
  function f(i,c,L,x)
    integer, intent(in) :: i, L
    real(8), intent(in) :: x, c(0:)
    real(8) :: f
    f = -c(i)/c(i+1) + cosh(0.5*x*(2*i-L))/cosh(0.5*x*(2*(i+1)-L))
  end function f

  function df(i,c,L,x)
    integer, intent(in) :: i,L
    real(8), intent(in) :: x, c(0:)
    real(8) :: df, h

    h = 0.001d0

    df = ( f(i,c,L,x+h) - f(i,c,L,x) )/h 
    !df=-(sinh(0.5*x*(2*i-L))*0.5*(2*i-L)*cosh(0.5*x*(2*(i+1)-L)) - &
    !     cosh(0.5*x*(2*i-L))*0.5*(2*(i+1)-L)*sinh(0.5*x*(2*(i+1)-L)) ) / &

    !cosh(0.5*x*(2*(i+1)-L))**2
  end function df

   function avr(x)
    real(8), intent(in) :: x(:)
    real(8) :: avr

    avr = sum(x)/size(x)
  end function avr

  
  function var(x)
    real(8), intent(in) :: x(:)
    real(8) :: var, avg

    avg = avr(x)
    var = sum((x-avg)**2)/(size(x)-1)
  end function var

  
  function std_err(x)
    real(8), intent(in) :: x(:)
    real(8) :: std_err

    std_err = sqrt(var(x)/size(x))
  end function std_err
  
end program mass_eff
