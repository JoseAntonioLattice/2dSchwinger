program main
  implicit none

  integer :: a

  a = this_image()

  call co_sum(a)
  print*, this_image(), a
  
  
end program main
