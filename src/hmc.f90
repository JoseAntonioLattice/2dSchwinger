#include "inc"
module hmc
  use dirac
  use CG
  use pbc, only : ip, im, sgnp
  use u1
  use constants, only : pi, i
  use arrays, only : up, chi, psi, phi
  use iso_fortran_env, only : dp => real64, i4 => int32
  implicit none

contains
  
    subroutine hmc_1(U, beta)
    use parameters, only : m0, Lx, Ly, N => MD_steps, Time => trajectory_length
    complex(dp), intent(inout) CODIM :: U(DIM)
    real(dp), intent(in) :: beta
    real(dp), dimension(2,Lx,Ly) :: Forces, p, pnew
    integer(i4) :: k, x, y, mu
    real(dp) :: r, deltaT, DS, S0
    real(dp) :: DeltaH
    logical :: condition    

    
    deltat = Time/N
    call generate_pi(p)
    
    chi%im = 0.0_dp
    call generate_pi(chi(:,1:Lx,1:Ly)%re)
    SYNC(chi)
    phi(:,1:Lx,1:Ly) = D(chi,U)
    SYNC(phi)
    psi(:,1:Lx,1:Ly) = conjugate_gradient(phi,U)
    SYNC(psi)
    S0 = sum(conjg(phi(:,1:Lx,1:Ly))*psi(:,1:Lx,1:Ly))
    
    
    !! k = 0
    !U_0
    up = u
    !P_0
    pnew = p 
    
    !Compute F[U_0]
    call compute_forces(Forces,beta,up,psi,chi)
    
    ! Compute P_{1/2} = P_0 + 0.5*dt*F[U_0]
    pnew = pnew + 0.5*deltaT*Forces
 
    ! k = 1, n -1
    do k = 1, N - 1
       !U_k = exp(i*dt*P_{k-1/2})U_{k-1}
       up(:,1:Lx,1:Ly) = up(:,1:Lx,1:Ly) * exp(i*DeltaT*pnew)
       SYNC(Up)
       !compute F[U_k]
       psi(:,1:Lx,1:Ly)  = conjugate_gradient(phi,Up)
       SYNC(psi)
       chi(:,1:Lx,1:Ly)  = Ddagger(psi,Up)
       SYNC(chi)
       call compute_forces(Forces,beta,up,psi,chi)
      
       !P_{k+1/2} = P_{k-1/2} + dt*F[U_k]
       pnew = pnew + deltaT*Forces
       
    end do
    
    ! k = n
    !U_n = exp(i*dt*P_{n-1/2})U_{n-1}
    up(:,1:Lx,1:Ly)  = up(:,1:Lx,1:Ly) * exp(i*DeltaT*pnew)
    SYNC(Up)
    psi(:,1:Lx,1:Ly) = conjugate_gradient(phi,Up)
    SYNC(psi)
    chi(:,1:Lx,1:Ly) = Ddagger(psi,Up)
    SYNC(chi)
    
    !compute F[U_n]
    call compute_forces(Forces,beta,up,psi,chi)

    !P_n = P_{n-1/2} + 0.5*dt*F[U_n]
    pnew = pnew + 0.5*DeltaT*Forces

    ! Metropolis step
    DS = sum(conjg(phi(:,1:Lx,1:Ly))*psi(:,1:Lx,1:Ly)) - S0 
    DeltaH = DH(u,up,p,pnew,beta)+DS
    
#if defined(PARALLEL)
    sync all
    call co_sum(DeltaH,result_image = 1)
    if(this_image() == 1) then 
#endif
       call random_number(r)
       condition = (r <= exp(-DeltaH))
#if defined(PARALLEL)
    end if
    call co_broadcast(condition,source_image=1)
#endif
    if( condition ) u = up


  end subroutine hmc_1

!!$  subroutine hmc_jaime(U,beta,Ntime,time)
!!$    use parameters, only : Lx, Ly
!!$    integer(i4), intent(in) ::  Ntime
!!$    complex(dp), dimension(2,Lx,Ly), intent(inout) :: U
!!$    real(dp), intent(in) :: beta, time
!!$
!!$    complex(dp), dimension(2,Lx,Ly) :: unew
!!$    real(dp), dimension(2,Lx,Ly) :: p, pnew, force
!!$
!!$    real(dp) :: r,u1,u2, DeltaH, dt
!!$    integer(i4) :: x, y, mu, k
!!$
!!$    dt = time/NTime
!!$    call generate_pi(p)
!!$     
!!$    unew = u
!!$    pnew = p
!!$
!!$    unew = unew*exp(0.5*i*dt*pnew)
!!$    !call compute_forces(force,beta,unew,L)
!!$    
!!$    do k = 1, Ntime - 2
!!$       pnew = pnew + dt*force
!!$       unew = unew*exp(i*dt*pnew)
!!$       !call compute_forces(force,beta,unew,L)
!!$    end do
!!$    
!!$    pnew = pnew + dt*force
!!$    unew = unew*exp(0.5*i*dt*pnew)
!!$    
!!$    call random_number(r)
!!$    
!!$    DeltaH = DH(u,unew,p,pnew,beta)
!!$    if( r <= exp(-DeltaH)) u = unew
!!$    
!!$  end subroutine hmc_jaime
!!$
!!$  subroutine hmc_knechtli(u,beta,nsteps,time)
!!$    use parameters, only : Lx, Ly
!!$    integer(i4), intent(in) ::  nsteps
!!$    complex(dp), intent(inout), dimension(2,Lx,Ly) :: u
!!$    real(dp), intent(in) :: time, beta
!!$    complex(dp), dimension(2,Lx,Ly) :: unew
!!$    real(dp), dimension(2,Lx,Ly) :: p, pnew, force 
!!$    real(dp) :: r, DeltaH, dt
!!$    integer(i4) :: x, y, mu, k
!!$
!!$    dt = time/nsteps
!!$    
!!$    call generate_pi(p)
!!$    
!!$    unew = u
!!$    pnew = p
!!$
!!$    !call compute_forces(force,beta,u,L)
!!$    do k = 1, nsteps
!!$       !P_{k-1/2} = P_{k-1} +F_{k-1}        
!!$       pnew = pnew + 0.5*dt*Force
!!$       !U_k = exp(i*dt*P_{k-1/2})U_{k-1}
!!$       unew = unew(mu,x,y) * exp(dt*i*pnew)
!!$       !F_k
!!$       !call compute_forces(force,beta,unew,L)
!!$       ! P_k = P_{k-1/2} + dt*F_k
!!$       pnew = pnew + 0.5*dt*Force
!!$    end do
!!$
!!$    !Metropolis step
!!$    call random_number(r)
!!$    if( r <= exp(-DeltaH)) u = unew
!!$    
!!$  end subroutine hmc_knechtli

  
  function DH(U,Unew,P,Pnew,beta)
    use parameters, only : Lx, Ly
    real(dp) :: DH
    complex(dp), dimension(DIM), intent(in) :: U, Unew
    real(dp), dimension(:,:,:), intent(in) :: P, Pnew
    real(dp), intent(in) :: beta
    integer(i4) :: x, y,mu
    real(dp) :: DeltaS

    DH = 0.0_dp
    DeltaS = 0.0_dp
    do x = 1, Lx
       do y = 1, Ly
          DeltaS = DeltaS + real(plaquette(u,[x,y]) - plaquette(unew,[x,y]))
          do mu = 1, 2
             DH = DH + (pnew(mu,x,y))**2 - (p(mu,x,y))**2
          end do
       end do
    end do

    DeltaS = beta*DeltaS

    DH = 0.5*DH + DeltaS
    
  end function DH

  
  subroutine generate_pi(p)
    use parameters, only : Lx, Ly
    real(dp), intent(out), dimension(2,Lx,Ly) :: p
    real(dp), dimension(2,Lx,Ly) :: u1, u2
    integer(i4) :: x, y

    call random_number(u1)
    call random_number(u2)
    
    p = sqrt(-2*log(1.0_dp-u1))*cos(2*pi*u2) 
   
  end subroutine generate_pi

    subroutine compute_forces(forces,beta,u,psi,chi)
    use parameters, only : Lx, Ly
    complex(dp), intent(in), dimension(DIM) :: u , psi, chi
    real(dp), intent(out), dimension(2,Lx,Ly) :: Forces
    real(dp), intent(in) :: beta
    complex(dp) :: stp, h
    integer(i4) :: x(2), xp1(2), xp2(2), x1, x2
    
    do x1 = 1, Lx
       do x2 = 1, Ly
          x = [x1,x2]
          xp1 = ip(x,1)
          xp2 = ip(x,2)
          forces(1,x(1),x(2)) = aimag( &
               -beta*u(1,x(1),x(2))*conjg(staples(u,x,1)) + &
               U(1,x(1),x(2))*conjg(psi(1,x(1),x(2))-psi(2,x(1),x(2)))*sgnp(x1)*(chi(1,xp1(1),xp1(2)) - chi(2,xp1(1),xp1(2)))-&
               conjg(U(1,x(1),x(2))*sgnp(x1)*(psi(1,xp1(1),xp1(2))+psi(2,xp1(1),xp1(2))))*(chi(1,x(1),x(2))+ chi(2,x(1),x(2))) &
               )
          
          forces(2,x(1),x(2)) = aimag( &
               -beta*u(2,x(1),x(2))*conjg(staples(u,x,2)) + &
               U(2,x(1),x(2))*(conjg(psi(1,x(1),x(2))+i*psi(2,x(1),x(2))))*(chi(1,xp2(1),xp2(2)) + i*chi(2,xp2(1),xp2(2)))+&
               conjg(U(2,x(1),x(2))*(psi(1,xp2(1),xp2(2))-i*psi(2,xp2(1),xp2(2))))*(-chi(1,x(1),x(2))+ i*chi(2,x(1),x(2))) &
               )
       end do
    end do

  end subroutine compute_forces

  
end module hmc
