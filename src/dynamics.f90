#if defined(PARALLEL)
#define DIM 2,0:Lx+1,0:Ly+1
#define DIM3 :,:,:
#define CODIM ,codimension[*]
#define CODIM2 ,codimension[:]
#define CDIM [*]
#define CDIM2 [:]
#define ALLOC ,allocatable
#define DIM_CODIM_ALLOC dimension(:,:,:),allocatable,codimension[:] 
#define SYNC(u) call sync_lattice(u)
#elif !defined(PARALLEL)
#define DIM 2,Lx,Ly
#define DIM3 2,Lx,Ly
#define CODIM 
#define CDIM
#define CODIM2
#define CDIM2
#define ALLOC
#define DIM_CODIM_ALLOC dimension(:,:,:)
#define SYNC(u)
#endif

module dynamics

  use iso_fortran_env, only : dp => real64, i4 => int32
  use pbc
  
  implicit none

  real(dp), parameter :: pi = acos(-1.0_dp)
  complex(dp), parameter :: i = (0.0_dp, 1.0_dp)
  complex(dp), dimension(2,2,2) :: sigma

contains
 
  subroutine set_memory(u,L,beta,betai,betaf,nbeta, plqaction,top_char,slb_top_char,pion_correlator,n_measurements)
    use parameters, only : Lx, Ly
    complex(dp), dimension(:,:,:), allocatable CODIM2 :: u
    integer(i4), intent(in) :: L(2)
    real(dp), allocatable, dimension(:) :: beta
    real(dp), allocatable, dimension(:,:) CODIM2 :: slb_top_char, pion_correlator
    real(dp), intent(in) :: betai, betaf
    real(dp), allocatable, dimension(:) CODIM2 :: plqaction, top_char
    integer(i4), intent(in) :: nbeta, n_measurements
    real(dp), allocatable, dimension(:) :: beta_copy
    integer(i4) :: i_beta
#ifdef PARALLEL
    integer(i4), dimension(:), allocatable :: cores CDIM2
#endif
    
    Lx = L(1)
    Ly = L(2)
#ifdef PARALLEL
    allocate(cores(2)[*])
    if(this_image() == 1) then
       print*, "Enter cores arrays"
       read(*,*) cores
       print*, "user typed # cores", cores
       
    end if
    call co_broadcast(cores,source_image=1)
    sync all
    
    Lx = Lx/cores(1)
    Ly = Ly/cores(2)
    call set_pbc([Lx,Ly],cores)
    
#else
    call set_pbc([Lx,Ly])
#endif

#ifdef PARALLEL
    !print*, this_image(), left, right, up, down, left_down, left_up, right_down, right_up
    print*, this_image(), ",", im([1,1],1), im([1,1],2)
    print*, this_image(), ",", ip([Lx,Ly],1), ip([Lx,Ly],2)
#endif
    allocate(u(DIM)CDIM)
    allocate(beta(nbeta))
    
    do i_beta = 1, nbeta 
       beta(i_beta) = betai + (betaf - betai)/(nbeta-1)*(i_beta-1)
    end do
   
    sigma(1,:,:) = reshape([0.0_dp,1.0_dp,1.0_dp,0.0_dp],[2,2])
    sigma(2,:,:) = reshape([(0.0_dp,0.0_dp),i,-i,(0.0_dp,0.0_dp)],[2,2])


    allocate(plqaction(n_measurements)CDIM)

    allocate(top_char(N_measurements)CDIM)
    allocate(slb_top_char(L(1),n_measurements)CDIM)
    allocate(pion_correlator(L(1),n_measurements)CDIM)
  end subroutine set_memory
  
  subroutine initialization(u,plqaction,top_char,slb_top_char,pion_correlator,beta,N_thermalization,N_measurements, N_skip)
    use parameters, only : L, Lx, Ly
    complex(dp), dimension(DIM) CODIM, intent(inout) :: u
    real(dp), dimension(:) CODIM, intent(out) :: plqaction, top_char
    integer(i4), intent(in) :: N_thermalization, N_measurements, N_skip
    real(dp), dimension(:,:) CODIM, intent(out) :: slb_top_char, pion_correlator
    real(dp), intent(in) :: beta
    real(dp) :: top(Lx,Ly)
    integer(i4) :: i_skip, i_sweeps, j

   ! print*, This_image(), "inside initialization"
    call thermalization(u, N_thermalization, beta)

    !open(unit = 69, file = 'data/history_top.dat',status = 'unknown')
    !do i_sweeps = 1, N_measurements
    i_sweeps = 0
    do while(i_sweeps < N_measurements)
       do i_skip = 1, n_skip
          call sweeps(u,beta)
       end do
       
       !call topological_charge_density(u,top)


       
       !if(nint(abs(sum(top)/(2*pi))) == 2) then
          i_sweeps = i_sweeps + 1
          plqaction(i_sweeps) = action(u)
        !  top_char(i_sweeps) = sum(top)/(2*pi)
          
#ifdef PARALLEL
          sync all
          call co_sum(plqaction(i_sweeps),result_image = 1)
        !  call co_sum(top_char(i_sweeps),result_image = 1)
#endif
          !do j = 1, Lx
          !   slb_top_char(j,i_sweeps) = slab_top_char(top,j)
          !end do
          !pion_correlator(:,i_sweeps) = pion_propagator(U)/sqrt(1.0_dp*Ly)
          !write(69,*) top_char(i_sweeps)
          !flush(69)
       !end if
    end do
    
  end subroutine initialization
  
  subroutine thermalization(u,N_thermalization,beta)
    use parameters, only : Lx, Ly
    complex(dp), dimension(DIM) CODIM, intent(inout) :: u
    integer(i4) :: N_thermalization
    real(dp) :: beta

    integer(i4) :: i_sweeps

    do i_sweeps = 1, N_thermalization
       call sweeps(u,beta)
       !print*, this_image(), "Done sweep"
    end do
    
  end subroutine thermalization

  subroutine hot_start(u)
    use parameters, only : Lx, Ly    
    complex(dp), intent(out), dimension(DIM) CODIM :: u
    real(dp), dimension(2,Lx,Ly) :: phi
    
    call random_number(phi)

    phi = 2*pi*phi
    u(:,1:Lx,1:Ly) = exp(i*phi); SYNC(u)
  end subroutine hot_start

  subroutine cold_start(u)
    complex(dp), intent(out), dimension(:,:,:) :: u
    u = 1.0_dp
  end subroutine cold_start
  
  subroutine sweeps(u,beta)
    use parameters, only : Lx, Ly
    complex(dp), dimension(DIM) CODIM, intent(inout) :: u
    real(dp), intent(in) :: beta
    integer(i4) :: x,y,mu

    call hmc(U,beta)
       
  end subroutine sweeps


  subroutine hmc(U, beta)
    use parameters, only : m0, Lx, Ly, N => MD_steps, Time => trajectory_length
    complex(dp), intent(inout) CODIM :: U(DIM)
    real(dp), intent(in) :: beta
    real(dp), dimension(2,Lx,Ly) :: Forces, p, pnew
    complex(dp), dimension(DIM3) CODIM2 ALLOC :: Up, psi, chi, phi
    integer(i4) :: k, x, y, mu
    real(dp) :: r, deltaT, DS, S0
    real(dp) ALLOC CODIM2 :: DeltaH
    logical :: condition
    

#ifdef PARALLEL
    allocate(Up(DIM)[*])
    allocate(psi(DIM)[*])
    allocate(phi(DIM)[*])
    allocate(chi(DIM)[*])
    allocate(DeltaH[*])
#endif
    
    deltat = Time/N
    call generate_pi(p)
    
    chi%im = 0.0_dp
    call generate_pi(chi(:,1:Lx,1:Ly)%re); SYNC(chi)
    phi(:,1:Lx,1:Ly) = D(chi,U); SYNC(phi)
    print*, phi(1,0,1), phi(1,Lx,1)[left]
    psi(:,1:Lx,1:Ly) = conjugate_gradient(phi,U); SYNC(psi)
    S0 = sum(conjg(phi(:,1:Lx,1:Ly))*psi(:,1:Lx,1:Ly))
    
    !! k = 0
    !U_0
    up = u
    !P_0
    pnew = p 
    
    !Compute F[U_0]
    call compute_forces(Forces,beta,up,psi, chi)
    
    ! Compute P_{1/2} = P_0 + 0.5*dt*F[U_0]
    pnew = pnew + 0.5*deltaT*Forces
 
    ! k = 1, n -1
    do k = 1, N - 1
       !U_k = exp(i*dt*P_{k-1/2})U_{k-1}
       up(:,1:Lx,1:Ly) = up(:,1:Lx,1:Ly) * exp(i*DeltaT*pnew); SYNC(up)
       !compute F[U_k]
       psi(:,1:Lx,1:Ly) = conjugate_gradient(phi,Up); SYNC(psi)
       chi(:,1:Lx,1:Ly) = Ddagger(psi,Up); SYNC(chi)
       call compute_forces(Forces,beta,up,psi,chi)
      
       !P_{k+1/2} = P_{k-1/2} + dt*F[U_k]
       pnew = pnew + deltaT*Forces
       
    end do
    !print*, this_image(), "done MC dynamics"
    ! k = n
    !U_n = exp(i*dt*P_{n-1/2})U_{n-1}
    up(:,1:Lx,1:Ly)  = up(:,1:Lx,1:Ly) * exp(i*DeltaT*pnew); SYNC(up)
    psi(:,1:Lx,1:Ly) = conjugate_gradient(phi,Up); SYNC(psi)
    chi(:,1:Lx,1:Ly) = Ddagger(psi,Up); SYNC(chi)
    
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

#ifdef PARALLEL
    deallocate(Up)
    deallocate(psi)
    deallocate(phi)
    deallocate(chi)
    deallocate(DeltaH)
    !print*, this_image(), "inside hmc. Deallocation suceesful"
#endif

  end subroutine hmc

  subroutine hmc_jaime(U,beta,Ntime,time)
    use parameters, only : Lx, Ly
    integer(i4), intent(in) ::  Ntime
    complex(dp), dimension(2,Lx,Ly), intent(inout) :: U
    real(dp), intent(in) :: beta, time

    complex(dp), dimension(2,Lx,Ly) :: unew
    real(dp), dimension(2,Lx,Ly) :: p, pnew, force

    real(dp) :: r,u1,u2, DeltaH, dt
    integer(i4) :: x, y, mu, k


    dt = time/NTime
    call generate_pi(p)
     
    unew = u
    pnew = p

    unew = unew*exp(0.5*i*dt*pnew)
    !call compute_forces(force,beta,unew,L)
    
    do k = 1, Ntime - 2
       pnew = pnew + dt*force
       unew = unew*exp(i*dt*pnew)
       !call compute_forces(force,beta,unew,L)
    end do
    
    pnew = pnew + dt*force
    unew = unew*exp(0.5*i*dt*pnew)
    
    call random_number(r)
    
    DeltaH = DH(u,unew,p,pnew,beta)
    if( r <= exp(-DeltaH)) u = unew
    
  end subroutine hmc_jaime

  subroutine hmc_knechtli(u,beta,nsteps,time)
    use parameters, only : Lx, Ly
    integer(i4), intent(in) ::  nsteps
    complex(dp), intent(inout), dimension(2,Lx,Ly) :: u
    real(dp), intent(in) :: time, beta
    complex(dp), dimension(2,Lx,Ly) :: unew
    real(dp), dimension(2,Lx,Ly) :: p, pnew, force 
    real(dp) :: r, DeltaH, dt
    integer(i4) :: x, y, mu, k

    dt = time/nsteps
    
    call generate_pi(p)
    
    unew = u
    pnew = p

    !call compute_forces(force,beta,u,L)
    do k = 1, nsteps
       !P_{k-1/2} = P_{k-1} +F_{k-1}        
       pnew = pnew + 0.5*dt*Force
       !U_k = exp(i*dt*P_{k-1/2})U_{k-1}
       unew = unew(mu,x,y) * exp(dt*i*pnew)
       !F_k
       !call compute_forces(force,beta,unew,L)
       ! P_k = P_{k-1/2} + dt*F_k
       pnew = pnew + 0.5*dt*Force
    end do

    !Metropolis step
    call random_number(r)
    if( r <= exp(-DeltaH)) u = unew
    
  end subroutine hmc_knechtli

  
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
    real(dp) :: u1, u2, u3, u4
    integer(i4) :: x, y

    do x = 1, Lx
       do y = 1, Ly
          call random_number(u1)
          call random_number(u2)
          call random_number(u3)
          call random_number(u4)
          p(1,x,y) = sqrt(-2*log(1.0_dp-u1))*cos(2*pi*u2) 
          p(2,x,y) = sqrt(-2*log(1.0_dp-u3))*cos(2*pi*u4)
       end do
    end do
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
               -beta*u(1,x(1),x(2))*conjg(staples(u,x,1))  + &
               U(1,x(1),x(2))*conjg(psi(1,x(1),x(2))-psi(2,x(1),x(2)))*sgnp(x1)*(chi(1,xp1(1),xp1(2)) - chi(2,xp1(1),xp1(2)))-&
               conjg(U(1,x(1),x(2))*sgnp(x1)*(psi(1,xp1(1),xp1(2))+psi(2,xp1(1),xp1(2))))*(chi(1,x(1),x(2))+ chi(2,x(1),x(2))) &
               )
          
          forces(2,x(1),x(2)) = aimag( &
               -beta*u(2,x(1),x(2))*conjg(staples(u,x,2))  + &
               U(2,x(1),x(2))*(conjg(psi(1,x(1),x(2))+i*psi(2,x(1),x(2))))*(chi(1,xp2(1),xp2(2)) + i*chi(2,xp2(1),xp2(2)))+&
               conjg(U(2,x(1),x(2))*(psi(1,xp2(1),xp2(2))-i*psi(2,xp2(1),xp2(2))))*(-chi(1,x(1),x(2))+ i*chi(2,x(1),x(2))) &
               )
       end do
    end do

  end subroutine compute_forces

  function D(phi,U)
    use parameters, only : Lx, Ly, m0
    complex(dp), dimension(DIM), intent(in) :: U, phi
    complex(dp), dimension(2,Lx,Ly) :: D
    integer(i4) :: x1, x2, mu, alpha, beta
    integer(i4), dimension(2) :: x, xm1, xm2,xp1, xp2
      
    do x1 = 1, Lx
       do x2 = 1, Ly
          x = [x1,x2]
          xm1 = im(x,1)
          xm2 = im(x,2)
          xp1 = ip(x,1)
          xp2 = ip(x,2)
          D(1,x(1),x(2)) = (m0+2.0_dp) * phi(1,x(1),x(2)) - 0.5*( &
               U(1,x(1),x(2))*sgnp(x1) *(phi(1,xp1(1),xp1(2)) -  phi(2,xp1(1),xp1(2))) + &
               U(2,x(1),x(2))          *(phi(1,xp2(1),xp2(2)) +i*phi(2,xp2(1),xp2(2))) + &
               conjg(U(1,xm1(1),xm1(2)))*sgnm(x1)*(phi(1,xm1(1),xm1(2)) +  phi(2,xm1(1),xm1(2))) + &
               conjg(U(2,xm2(1),xm2(2)))         *(phi(1,xm2(1),xm2(2)) -i*phi(2,xm2(1),xm2(2))) &
               )
          D(2,x(1),x(2)) = (m0+2.0_dp) * phi(2,x(1),x(2)) - 0.5*( &
               U(1,x(1),x(2))*sgnp(x1) *(  -phi(1,xp1(1),xp1(2)) + phi(2,xp1(1),xp1(2))) + &
               U(2,x(1),x(2))          *(-i*phi(1,xp2(1),xp2(2)) + phi(2,xp2(1),xp2(2))) + &
               conjg(U(1,xm1(1),xm1(2)))*sgnm(x1)*(  phi(1,xm1(1),xm1(2)) + phi(2,xm1(1),xm1(2))) + &
               conjg(U(2,xm2(1),xm2(2)))         *(i*phi(1,xm2(1),xm2(2)) + phi(2,xm2(1),xm2(2))) &
               )
       end do
    end do
    
  end function D

  function Ddagger(phi,U)
    use parameters, only : Lx, Ly, m0
    complex(dp), dimension(DIM), intent(in) :: U, phi
    complex(dp), dimension(2,Lx,Ly) :: Ddagger
    integer(i4) :: x(2), xm1(2), xm2(2),xp1(2), xp2(2), alpha, beta, mu
    integer(i4) :: x1, x2
    
    do x1 = 1, Lx
       do x2 = 1, Ly
          x = [x1,x2]
          xm1 = im(x,1)
          xm2 = im(x,2)
          xp1 = ip(x,1)
          xp2 = ip(x,2)
          Ddagger(1,x(1),x(2)) = (m0+2.0_dp) * phi(1,x(1),x(2)) - 0.5*( &
               conjg(U(1,xm1(1),xm1(2)))*sgnm(x1)*(phi(1,xm1(1),xm1(2)) -  phi(2,xm1(1),xm1(2))) + &
               conjg(U(2,xm2(1),xm2(2)))         *(phi(1,xm2(1),xm2(2)) +i*phi(2,xm2(1),xm2(2))) + &
               U(1,x(1),x(2))*sgnp(x1)*(phi(1,xp1(1),xp1(2)) +  phi(2,xp1(1),xp1(2))) + &
               U(2,x(1),x(2))         *(phi(1,xp2(1),xp2(2)) -i*phi(2,xp2(1),xp2(2))) &
               )
          Ddagger(2,x(1),x(2)) = (m0+2.0_dp) * phi(2,x(1),x(2)) - 0.5*( &
               conjg(U(1,xm1(1),xm1(2)))*sgnm(x1)*(  -phi(1,xm1(1),xm1(2)) + phi(2,xm1(1),xm1(2))) + &
               conjg(U(2,xm2(1),xm2(2)))         *(-i*phi(1,xm2(1),xm2(2)) + phi(2,xm2(1),xm2(2))) + &
               U(1,x(1),x(2))*sgnp(x1)*(  phi(1,xp1(1),xp1(2)) + phi(2,xp1(1),xp1(2))) + &
               U(2,x(1),x(2))         *(i*phi(1,xp2(1),xp2(2)) + phi(2,xp2(1),xp2(2))) &
               )
       end do
    end do
        
  end function Ddagger
  
  function DDdagger(phi,U)
    use parameters, only : Lx, Ly
    complex(dp), dimension(DIM), intent(in) :: U, phi
    complex(dp), dimension(2,Lx,Ly) :: DDdagger

    DDdagger = D(Ddagger(phi,U),U)

  end function DDdagger
  
  function Dinv(phi,U)
    use parameters, only : Lx, Ly
    complex(dp), dimension(DIM), intent(in) :: phi, U
    complex(dp), dimension(2,Lx,Ly) :: Dinv, CG

    CG = conjugate_gradient(phi,U)
    Dinv = Ddagger(CG,U) 
    
  end function Dinv

  subroutine check_CG()
    use parameters, only : Lx, Ly
    complex(dp), dimension(DIM) :: U, phi,phi_p, psi,Dphi,chi,D1,D2
    real(dp), dimension(2,Lx,Ly) :: r

    
    
  end subroutine check_CG
 
  function conjugate_gradient(phi,U) result(x) 
    use parameters, only : Lx, Ly, max_iter, tol
    complex(dp), dimension(DIM), intent(in) :: U, phi
    complex(dp), dimension(2,Lx,Ly) :: x
    complex(dp), dimension(2,Lx,Ly) :: r,Ap
    complex(dp), dimension(DIM3) CODIM2 ALLOC :: p
    complex(dp) :: alpha, pAp
    real(dp)    :: phi_norm2, beta, rr, rr_new
    integer(i4) :: k

#ifdef PARALLEL
    allocate(p(2,0:Lx+1,0:Ly+1)[*])
#endif

    x = phi(:,1:Lx,1:Ly)
    Ap = DDdagger(phi,U)
    r = phi(:,1:Lx,1:Ly) - Ap
    p(:,1:Lx,1:Ly) = r; SYNC(p)

    rr = real(sum(r*conjg(r)))
    phi_norm2 = real(sum(phi(:,1:Lx,1:Ly)*conjg(phi(:,1:Lx,1:Ly))))

#ifdef PARALLEL
    call co_sum(rr)
    call co_sum(phi_norm2)
#endif

    do k = 1, max_iter
       
       Ap = DDdagger(p,U)
       pAp = sum(p(:,1:Lx,1:Ly)*conjg(Ap))
       
#ifdef PARALLEL
       call co_sum(pAp)
#endif

       alpha = rr/pAp
       
       x = x + alpha*p(:,1:Lx,1:Ly)
       r = r - alpha*Ap

       rr_new = real(sum(r*conjg(r)))
#ifdef PARALLEL
       call co_sum(rr_new)
#endif
       if( rr_new*phi_norm2 < tol ) exit
       beta = rr_new/rr
       p(:,1:Lx,1:Ly) = r + beta*p(:,1:Lx,1:Ly); SYNC(p)
       rr = rr_new
    end do
    
  end function conjugate_gradient

  function staples(u,x,mu)
    use parameters, only : Lx, Ly
    complex(dp) :: staples
    complex(dp), dimension(DIM), intent(in) :: u
    integer(i4), intent(in) :: x(2), mu
    integer(i4), dimension(2) :: x2, x3, x4, x5, x6

    integer(i4) :: nu
    
    if ( mu == 1 ) then
       nu = 2
    elseif( mu == 2)then
       nu = 1
    end if

    x2 = ip(x,nu)
    x3 = ip(x,mu)
    x4 = im(x,nu)
    x5 = x4
    x6 = im(x3,nu)
    
    staples = u(nu,x(1),x(2)) * u(mu,x2(1), x2(2)) * conjg( u(nu,x3(1), x3(2)) ) + &
         conjg( u(nu,x4(1),x4(2)) ) * u(mu,x5(1), x5(2)) * u(nu,x6(1), x6(2))
    
  end function staples
  
  function DS(uold, unew, beta,stp)

    real(dp) :: DS,beta
    complex(dp) :: uold, unew, stp

    DS = -beta * real( (unew - uold) * conjg(stp) )
    
  end function DS
  
  function plaquette(u,x)
    use parameters, only : Lx, Ly
    complex(dp) :: plaquette
    
    complex(dp), dimension(DIM), intent(in) :: u
    integer(i4), intent(in) :: x(2)
    integer(i4), dimension(2) :: x2, x3


    x2 = ip(x,1)
    x3 = ip(x,2)
    
    plaquette = U(1,x(1),x(2)) * U(2,x2(1),x2(2)) * &
          conjg(U(1,x3(1),x3(2))) * conjg(U(2,x(1),x(2)))
    
  end function plaquette

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

    !action2 = beta*(L**2 - action2)
    
  end function action2

  subroutine topological_charge_density(u,top)
    use parameters, only : Lx, Ly
    complex(dp), dimension(DIM), intent(in) :: u
    real(dp), intent(out) :: top(Lx,Ly) 
    integer(i4) :: i, j
    complex(dp) :: plq
    
    do i = 1, size(u(1,:,1))
       do j = 1, size(u(1,1,:))
          plq = plaquette(u,[i,j])
          top(i,j) = atan2(plq%im,plq%re)!/(2*pi)
       end do
    end do
       
  end subroutine topological_charge_density

  function slab_top_char(top,ix)
    real(dp), intent(in) :: top(:,:)
    integer(i4), intent(in) :: ix
    real(dp) :: slab_top_char
    integer(i4) :: i, j
    complex(dp) :: plq

    slab_top_char = (sum(top(1:ix,:))/(2*pi))**2!*sum(top(ix+1:,:))
    
  end function slab_top_char


#ifdef PARALLEL
  subroutine sync_lattice(u)
    use parameters, only : Lx, Ly
    complex(dp), dimension(DIM) :: u[*]

    sync all
    !Send edges
    u(:, Lx+1, 1:Ly)[left]  = u(:, 1 , 1:Ly)
    u(:, 0   , 1:Ly)[right] = u(:, Lx, 1:Ly)
    u(:, 1:Lx, Ly+1)[down]  = u(:, 1:Lx, 1)
    u(:, 1:Lx, 0   )[up]    = u(:, 1:Lx, Ly)

    sync all
    !Send corners
    u(:, Lx+1, Ly+1 )[left_down]  = u(:, 1 , 1 )
    u(:, Lx+1, 0    )[left_up]    = u(:, 1 , Ly)
    u(:, 0   , Ly+1 )[right_down] = u(:, Lx, 1 )
    u(:, 0   , 0    )[right_up]   = u(:, Lx, Ly)

    sync all
  end subroutine sync_lattice
#endif
  
end module dynamics
