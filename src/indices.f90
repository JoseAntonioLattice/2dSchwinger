module indices

  implicit none

contains

  function get_index(vector,length_lattice)

    integer :: get_index

    integer, dimension(:), intent(in) :: vector
    integer :: dimension_space
    integer, intent(in) :: length_lattice(:)

    integer :: suma, prod
    integer :: i

    dimension_space = size(length_lattice)

    suma = vector(1)
    prod = 1
    if( dimension_space > 1)then
       do i = 2, dimension_space
          suma = suma + (vector(i) - 1) * product(length_lattice(1:i-1))
       end do
    end if

    get_index = suma


  end function get_index


  function get_index_array(idx,L) result(vector)

    integer, intent(in) :: idx
    integer :: d
    integer, intent(in) :: L(:)

    integer, dimension(size(L)) :: vector

    integer :: i, n, modx, suma, prod1, prod2

    d = size(L)
    
    modx = mod(idx,L(1))
    vector(1) = modx
    if(modx == 0) vector(1) = L(1)

    suma = vector(1)
    do i = 2, d
      modx = mod(idx,product(L(1:i)))
      if (i > 2) suma = suma + product(L(1:i-2))*(vector(i-1)-1)
      vector(i) = (modx - suma)/product(L(1:i-1)) + 1
      if(modx == 0) vector(i) = L(i)
    end do

  end function get_index_array


end module indices
