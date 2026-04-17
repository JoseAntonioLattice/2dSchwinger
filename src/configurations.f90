#include "inc"
module configurations
  use pbc
  use iso_fortran_env, only : dp => real64, i4 => int32
  use number2string
  use check_files_directories 
  use parameters, only : L, Lx, Ly,m0
   
  implicit none
contains

  subroutine save_configuration(U,beta)
     complex(dp), intent(in) :: U(DIM)
    real(dp), intent(in) :: beta
    integer(i4) :: un, ix,ex, iy, ey, a(2)
    complex(dp), allocatable :: U_global(:,:,:)CDIM2
    character(:), allocatable :: path,file_name
    character(100), dimension(6) :: directory_array

#ifndef PARALLEL
    allocate(U_global(2,L(1),L(2)))
    U_global = U
#else   
    allocate(U_global(2,L(1),L(2))[*])
    
    a = get_index_array(this_image(),cores)
    ix = L(1)/cores(1)*(a(1)-1)+1
    ex = L(1)/cores(1)*a(1)
    iy = L(2)/cores(2)*(a(2)-1)+1
    ey = L(2)/cores(2)*a(2)
    sync all
    U_global(:,ix:ex,iy:ey)[1] = U(:,1:Lx,1:Ly)
    sync all
    if( this_image()==1 ) then
#endif
       directory_array = [character(100):: "data","configurations", "m0="//real2str(m0,nint(log10(abs(m0)))+1,4),&
            "Lx="//trim(int2str(Lx)),"Lt="//trim(int2str(Ly)),"beta="//trim(real2str(beta,1,4))]
    
       call check_directory(directory_array,path)
       call numbered_files(path,"U",".bin",file_name)
       open(newunit = un, file = file_name, access = "sequential", form = "unformatted")
       !open(newunit = un, file = "U"//int2str(n)//".dat")
       write(un) U_global
       close(un)
       deallocate(U_global)
#ifdef PARALLEL
    endif
    sync all
#endif
    
  end subroutine save_configuration

  subroutine read_configuration(U,filename)
    complex(dp), intent(out) :: U(DIM)
    character(*), intent(in) :: filename 
    integer(i4) :: un, ix,ex, iy, ey, a(2)
    complex(dp), allocatable :: U_global(:,:,:)CDIM2

    allocate(U_global(2,L(1),L(2))CDIM)
#ifdef PARALLEL   
    a = get_index_array(this_image(),cores)
    ix = L(1)/cores(1)*(a(1)-1)+1
    ex = L(1)/cores(1)*a(1)
    iy = L(2)/cores(2)*(a(2)-1)+1
    ey = L(2)/cores(2)*a(2)
    
    if( this_image()==1 ) then
#endif
       open(newunit = un, file = filename, access = "sequential", form = "unformatted")  
       read(un) U_global
       close(un)
       
#ifdef PARALLEL
    endif
    call co_broadcast(U_global,source_image=1)
    U(:,1:Lx,1:Ly) = U_global(:,ix:ex,iy:ey)
#else
    U = U_global
#endif
    deallocate(U_global)
  end subroutine read_configuration



end module configurations
