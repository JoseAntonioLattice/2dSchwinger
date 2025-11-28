module dynamics

  use iso_fortran_env, only : dp => real64, i4 => int32
  use pbc
  
  implicit none

  real(dp), parameter :: pi = acos(-1.0_dp)
  complex(dp), parameter :: i = (0.0_dp, 1.0_dp)
  complex(dp), dimension(2,2,2) :: sigma

contains
 
  subroutine set_memory(u,L,beta,betai,betaf,nbeta, plqaction,top_char,slb_top_char,pion_correlator,n_measurements)

    complex(dp), allocatable, dimension(:,:,:) :: u
    integer(i4), intent(in) :: L(2)
    real(dp), allocatable, dimension(:) :: beta
    real(dp), allocatable, dimension(:,:) :: slb_top_char, pion_correlator
    real(dp), intent(in) :: betai, betaf
    real(dp), allocatable, dimension(:) :: plqaction, top_char
    integer(i4), intent(in) :: nbeta, n_measurements
    real(dp), allocatable, dimension(:) :: beta_copy
    integer(i4) :: i_beta
    
    call set_pbc(L(1))
    allocate(u(2,L(1),L(2)))
    allocate(beta(nbeta))

    do i_beta = 1, nbeta 
       beta(i_beta) = betai + (betaf - betai)/(nbeta-1)*(i_beta-1)
    end do
   
    !print*, sgnp
    !print*, sgnm
    sigma(1,:,:) = reshape([0.0_dp,1.0_dp,1.0_dp,0.0_dp],[2,2])
    sigma(2,:,:) = reshape([(0.0_dp,0.0_dp),i,-i,(0.0_dp,0.0_dp)],[2,2])
    allocate(plqaction(n_measurements),top_char(N_measurements))
    allocate(slb_top_char(L(1),n_measurements))
    allocate(pion_correlator(L(1),n_measurements))
  end subroutine set_memory
  
  subroutine initialization(u,plqaction,top_char,slb_top_char,pion_correlator,beta,N_thermalization,N_measurements, N_skip)
    use parameters, only : Lx, Ly
    complex(dp), dimension(:,:,:), intent(inout) :: u
    real(dp), dimension(:), intent(out) :: plqaction, top_char
    integer(i4), intent(in) :: N_thermalization, N_measurements, N_skip
    real(dp), dimension(:,:), intent(out) :: slb_top_char, pion_correlator
    real(dp), intent(in) :: beta
    real(dp) :: top(size(u(1,:,1)),size(u(1,1,:)))
    integer(i4) :: i_skip, i_sweeps, j

    call thermalization(u, N_thermalization, beta)

    open(unit = 69, file = 'data/history_top.dat',status = 'unknown')
    !do i_sweeps = 1, N_measurements
    i_sweeps = 0
    do while(i_sweeps < N_measurements)
       do i_skip = 1, n_skip
          call sweeps(u,beta)
       end do
       
       call topological_charge_density(u,top)
   
       !if(nint(abs(sum(top)/(2*pi))) == 2) then
          i_sweeps = i_sweeps + 1
          plqaction(i_sweeps) = action(u)
          top_char(i_sweeps) = sum(top)/(2*pi)
          do j = 1, Lx
             slb_top_char(j,i_sweeps) = slab_top_char(top,j)
          end do
          !pion_correlator(:,i_sweeps) = pion_propagator(U)/sqrt(1.0_dp*Ly)
          write(69,*) top_char(i_sweeps)
          flush(69)
       !end if
    end do
    
  end subroutine initialization
  
  subroutine thermalization(u,N_thermalization,beta)

    complex(dp), dimension(:,:,:), intent(inout) :: u
    integer(i4) :: N_thermalization
    real(dp) :: beta

    integer(i4) :: i_sweeps

    do i_sweeps = 1, N_thermalization
       call sweeps(u,beta)
    end do
    
  end subroutine thermalization

  subroutine hot_start(u)
    use parameters, only : Lx, Ly
    
    complex(dp), intent(out), dimension(2,Lx,Ly) :: u
    real(dp), dimension(2,Lx,Ly) :: phi

    call random_number(phi)

    phi = 2*pi*phi
    u = exp(i*phi)
    
  end subroutine hot_start

  subroutine cold_start(u)

    complex(dp), intent(out), dimension(:,:,:) :: u

    u = 1.0_dp
    
  end subroutine cold_start
  
  subroutine sweeps(u,beta)

    complex(dp), dimension(:,:,:), intent(inout) :: u
    real(dp), intent(in) :: beta
    integer(i4) :: x,y,mu
    
    call hmc(U,beta)
   
  end subroutine sweeps

  subroutine metropolis(u,x,mu,beta)

    complex(dp), dimension(:,:,:), intent(inout) :: u
    integer(i4), intent(in) :: x(2), mu
    real(dp), intent(in) :: beta
    
    complex(dp) :: u_new
    real(dp) :: phi
    real(dp) :: p, deltaS
    real(dp) :: r
    
    call random_number(phi)
    phi = 2*pi*phi
    u_new = exp(i*phi)

    deltaS = DS(u(mu,x(1),x(2)),u_new,beta,staples(u,x,mu))

    call random_number(r)
    p = min(1.0_dp,exp(-DeltaS))
    if ( r <= p )then
       u(mu,x(1),x(2)) = u_new
    end if
    
  end subroutine metropolis

  subroutine hmc(U, beta)
    use parameters, only : m0, Lx, Ly, N => MD_steps, Time => trajectory_length
    complex(dp), intent(inout) :: U(2,Lx,Ly)
    real(dp), intent(in) :: beta
    complex(dp), dimension(2,Lx,Ly) ::  Up
    real(dp), dimension(2,Lx,Ly) :: Forces, p, pnew
    complex(dp), dimension(2,Lx,Ly) :: psi, chi,phi
    integer(i4) :: k, x, y, mu
    real(dp) :: DeltaH,r, deltaT, DS, S0

    deltat = Time/N
    call generate_pi(p)
    chi%im = 0.0_dp
    call generate_pi(chi%re)

    phi = D(chi,U)
    psi = conjugate_gradient(phi,U)
    S0 = sum(conjg(phi)*psi)
    
    !! k = 0
    !U_0
    up = u
    !P_0
    pnew = p

    !Compute F[U_0]
    call compute_forces(Forces,beta,u,psi, chi)
    

    ! Compute P_{1/2} = P_0 + 0.5*dt*F[U_0]
    pnew = pnew + 0.5*deltaT*Forces
 
    ! k = 1, n -1
    do k = 1, N - 1
       !U_k = exp(i*dt*P_{k-1/2})U_{k-1}
       up = up * exp(i*DeltaT*pnew)
       
       !compute F[U_k]
       psi = conjugate_gradient(phi,Up)
       chi = Ddagger(psi,Up)
       call compute_forces(Forces,beta,up,psi,chi)
      
       !P_{k+1/2} = P_{k-1/2} + dt*F[U_k]
       pnew = pnew + deltaT*Forces
       
    end do 

    ! k = n
    !U_n = exp(i*dt*P_{n-1/2})U_{n-1}
    up = up * exp(i*DeltaT*pnew)

    psi = conjugate_gradient(phi,Up)
    chi = Ddagger(psi,Up)
    !compute F[U_n]
    call compute_forces(Forces,beta,up,psi,chi)

    
    !P_n = P_{n-1/2} + 0.5*dt*F[U_n]
    p = p + 0.5*DeltaT*Forces

    ! Metropolis step

    DS = sum(conjg(phi)*psi) - S0 
    DeltaH = DH(u,up,p,pnew,beta)+DS
    call random_number(r)
    if( r <= exp(-DeltaH) ) u = up
    
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
    real(dp) :: DH
    complex(dp), dimension(:,:,:), intent(in) :: U, Unew
    real(dp), dimension(:,:,:), intent(in) :: P, Pnew
    real(dp), intent(in) :: beta
    integer(i4) :: x, y,mu, Lx, Ly
    real(dp) :: DeltaS

    Lx = size(U(1,:,1))
    Ly = size(U(1,1,:))
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
    complex(dp), intent(in), dimension(2,Lx,Ly) :: u , psi, chi
    real(dp), intent(out), dimension(2,Lx,Ly) :: Forces
    real(dp), intent(in) :: beta
    complex(dp) :: stp, h
    integer(i4) :: x(2), xp1(2), xp2(2), ii, jj

    
    do ii = 1, Lx
       do jj = 1, Ly
          x = [ii,jj]
          xp1 = ipf(x,1)
          xp2 = ipf(x,2)
          forces(1,x(1),x(2)) = aimag( &
               -beta*u(1,x(1),x(2))*conjg(staples(u,x,1)) + &
               U(1,x(1),x(2))*conjg(psi(1,x(1),x(2))-psi(2,x(1),x(2)))*sgnp(ii)*(chi(1,xp1(1),xp1(2)) - chi(2,xp1(1),xp1(2)))-&
               conjg(U(1,x(1),x(2))*sgnp(ii)*(psi(1,xp1(1),xp1(2))+psi(2,xp1(1),xp1(2))))*(chi(1,x(1),x(2))+ chi(2,x(1),x(2))) &
               )
          
          forces(2,x(1),x(2)) = aimag( &
               -beta*u(2,x(1),x(2))*conjg(staples(u,x,2)) + &
               U(2,x(1),x(2))*(conjg(psi(1,x(1),x(2))+i*psi(2,x(1),x(2))))*(chi(1,xp2(1),xp2(2)) + i*chi(2,xp2(1),xp2(2)))+&
               conjg(U(2,x(1),x(2))*(psi(1,xp2(1),xp2(2))-i*psi(2,xp2(1),xp2(2))))*(-chi(1,x(1),x(2))+ i*chi(2,x(1),x(2))) &
               )
       end do
    end do


    
  end subroutine compute_forces

  function D(phi,U)
    use parameters, only : Lx, Ly, m0
    complex(dp), dimension(2,Lx,Ly), intent(in) :: U, phi
    complex(dp), dimension(2,Lx,Ly) :: D
    integer(i4) :: ii,jj,mu, alpha, beta
    integer(i4), dimension(2) :: x, xm1, xm2,xp1, xp2
    complex(dp) :: h

    D = (0.0_dp,0.0_dp)
    do ii = 1, Lx
       do jj = 1, Ly
          x = [ii,jj]
          xm1 = imf(x,1)
          xm2 = imf(x,2)
          xp1 = ipf(x,1)
          xp2 = ipf(x,2)
          D(1,x(1),x(2)) = (m0+2.0_dp) * phi(1,x(1),x(2)) - 0.5*( &
               U(1,x(1),x(2))*sgnp(ii)*(phi(1,xp1(1),xp1(2)) -  phi(2,xp1(1),xp1(2))) + &
               U(2,x(1),x(2))         *(phi(1,xp2(1),xp2(2)) +i*phi(2,xp2(1),xp2(2))) + &
               conjg(U(1,xm1(1),xm1(2)))*sgnm(ii)*(phi(1,xm1(1),xm1(2)) +  phi(2,xm1(1),xm1(2))) + &
               conjg(U(2,xm2(1),xm2(2)))         *(phi(1,xm2(1),xm2(2)) -i*phi(2,xm2(1),xm2(2))) &
               )
          D(2,x(1),x(2)) = (m0+2.0_dp) * phi(2,x(1),x(2)) - 0.5*( &
               U(1,x(1),x(2))*sgnp(ii)*(  -phi(1,xp1(1),xp1(2)) + phi(2,xp1(1),xp1(2))) + &
               U(2,x(1),x(2))         *(-i*phi(1,xp2(1),xp2(2)) + phi(2,xp2(1),xp2(2))) + &
               conjg(U(1,xm1(1),xm1(2)))*sgnm(ii)*(  phi(1,xm1(1),xm1(2)) + phi(2,xm1(1),xm1(2))) + &
               conjg(U(2,xm2(1),xm2(2)))         *(i*phi(1,xm2(1),xm2(2)) + phi(2,xm2(1),xm2(2))) &
               )
       end do
    end do
    
  end function D

  function Ddagger(phi,U)
    use parameters, only : Lx, Ly, m0
    complex(dp), dimension(2,Lx,Ly), intent(in) :: U, phi
    complex(dp), dimension(2,Lx,Ly) :: Ddagger
    integer(i4) :: x(2), xm1(2), xm2(2),xp1(2), xp2(2), alpha, beta, mu
    integer(i4) :: ii, jj
    
    do ii = 1, Lx
       do jj = 1, Ly
          x = [ii,jj]
          xm1 = imf(x,1)
          xm2 = imf(x,2)
          xp1 = ipf(x,1)
          xp2 = ipf(x,2)
          Ddagger(1,x(1),x(2)) = (m0+2.0_dp) * phi(1,x(1),x(2)) - 0.5*( &
               conjg(U(1,xm1(1),xm1(2)))*sgnm(ii)*(phi(1,xm1(1),xm1(2)) -  phi(2,xm1(1),xm1(2))) + &
               conjg(U(2,xm2(1),xm2(2)))         *(phi(1,xm2(1),xm2(2)) +i*phi(2,xm2(1),xm2(2))) + &
               U(1,x(1),x(2))*sgnp(ii)*(phi(1,xp1(1),xp1(2)) +  phi(2,xp1(1),xp1(2))) + &
               U(2,x(1),x(2))         *(phi(1,xp2(1),xp2(2)) -i*phi(2,xp2(1),xp2(2))) &
               )
          Ddagger(2,x(1),x(2)) = (m0+2.0_dp) * phi(2,x(1),x(2)) - 0.5*( &
               conjg(U(1,xm1(1),xm1(2)))*sgnm(ii)*(  -phi(1,xm1(1),xm1(2)) + phi(2,xm1(1),xm1(2))) + &
               conjg(U(2,xm2(1),xm2(2)))         *(-i*phi(1,xm2(1),xm2(2)) + phi(2,xm2(1),xm2(2))) + &
               U(1,x(1),x(2))*sgnp(ii)*(  phi(1,xp1(1),xp1(2)) + phi(2,xp1(1),xp1(2))) + &
               U(2,x(1),x(2))         *(i*phi(1,xp2(1),xp2(2)) + phi(2,xp2(1),xp2(2))) &
               )
       end do
    end do
    
    
  end function Ddagger
  
  function DDdagger(phi,U)
    use parameters, only : Lx, Ly
    complex(dp), dimension(2,Lx,Ly), intent(in) :: U, phi
    complex(dp), dimension(2,Lx,Ly) :: DDdagger

    DDdagger = D(Ddagger(phi,U),U)

  end function DDdagger
  
  function Dinv(phi,U)
    use parameters, only : Lx, Ly
    complex(dp), dimension(2,Lx,Ly), intent(in) :: phi, U
    complex(dp), dimension(2,Lx,Ly) :: Dinv

    Dinv = Ddagger(conjugate_gradient(phi,U),U) 
    
  end function Dinv

  subroutine check_CG()
    use parameters, only : Lx, Ly
    complex(dp), dimension(2,Lx,Ly) :: U, phi,phi_p, psi,Dphi,chi,D1,D2
    real(dp), dimension(2,Lx,Ly) :: r

    call random_number(r)

    U = exp(i*2*pi*r)
    chi%im = 0.0_dp
    call generate_pi(chi%re)

    phi = D(chi,U)

    !print*,"phi"
    !print*, phi
    !print*, "Dphi"
    !print*, D(phi,U)
    !print*, "Ddaggerphi"
    !print*, Ddagger(phi,U)
    !print*, "DDdaggerphi"
    !print*, DDdagger(phi,U)
    !psi = conjugate_gradient(phi,U)
    !phi_p = DDdagger(psi,U)

    D1 = conjugate_gradient(DDdagger(phi,U),U)!bcg(phi,U)
    D2 = DDdagger(conjugate_gradient(phi,U),U)!Dinv(phi,U)
    print*,phi(1,1,1)
    print*,D1(1,1,1)
    print*,D2(1,1,1)
  end subroutine check_CG
 
  function conjugate_gradient(phi,U) result(x)
    use parameters, only : Lx, Ly, max_iter, tol
    complex(dp), dimension(2,Lx,Ly), intent(in) :: U, phi
    complex(dp), dimension(2,Lx,Ly) :: x
    complex(dp), dimension(2,Lx,Ly) :: r,p,Ad
    complex(dp) :: alpha
    real(dp) :: err, phi_norm2, beta, rsq
    integer(i4) :: k

    x = phi
    Ad = DDdagger(x,U)
    r = phi - Ad
    p = r
    
    k = 0
    err = real(sum(r*conjg(r)))
    !print*, "err", err
    phi_norm2 = sqrt(real(sum(phi*conjg(phi))))
    !print*, "norm of phi", phi_norm2
    do while ( err > tol*phi_norm2 .and. k < max_iter)
       
       Ad = DDdagger(p,U)

       rsq = err
       alpha = rsq/sum(p*conjg(Ad))

       x = x + alpha*p
       r = r - alpha*Ad
       
       err = real(sum(r*conjg(r)))
       
       beta = err/rsq
       p = r + beta*p
       k = k + 1
       !print*,"k = ",k, "err old = ", rsq ,"err new = ", err
       if( k >= max_iter) stop "CG. Maximum number of iterations reached. No convergence."
    end do
    
  end function conjugate_gradient

  function bcg(phi,U) result(x)
    use parameters, only : Lx, Ly, max_iter, tol
    complex(dp), dimension(2,Lx,Ly), intent(in) :: U, phi
    complex(dp), dimension(2,Lx,Ly) :: x
    complex(dp), dimension(2,Lx,Ly) :: r,p,Ad,s,t,b,rtilde
    complex(dp) :: alpha, beta, rho_i, rho_i_2, omega
    real(dp) :: err, phi_norm, rsq
    integer(i4) :: it
    
    x = phi
    b = phi
    Ad = D(x,U)
    r = b-Ad
    rtilde = r
    phi_norm = sqrt(real(sum(phi*conjg(phi))))

    it = 0
    do while(it < max_iter)
       rho_i = sum(conjg(rtilde)*r)
       !if( rho_i <= 1e-10_dp ) stop "Method failed"
       if(it == 0) then
          p = r
       else
          beta = alpha*rho_i/(omega*rho_i_2)
          p = r + beta*(p-omega*Ad)
       end if

       Ad = D(p,U)
       alpha = rho_i/(sum(conjg(rtilde)*Ad))
       s = r - alpha*Ad
       err = sqrt(real(sum(conjg(s)*s)))
       if( err < tol*phi_norm )then
          x = x + alpha*p
          return
       end if
       t = D(s,U)
       omega = sum(conjg(s)*t)/real(sum(conjg(t)*t))
       r = s - omega*t
       x = x + alpha*p+omega*s
       rho_i_2 = rho_i
       it = it + 1
       if( it >= max_iter) stop "BCG. Maximum number of iterations reached. No convergence."
    end do
    
  end function bcg

  function pion_propagator(U)
    use parameters, only : Lx, Ly
    complex(dp), dimension(2,Lx,Ly), intent(in) :: U
    complex(dp), dimension(2,Lx,Ly) :: S1, S2, DinvS1, DinvS2
    real(dp), dimension(Lx) :: pion_propagator
    integer(i4) :: x, y

    S1 = (0.0_dp,0.0_dp)
    S2 = (0.0_dp,0.0_dp)

    S1(1,1,1) = (1.0_dp,0.0_dp)
    S2(2,1,1) = (1.0_dp,0.0_dp)
    
    DinvS1 = Dinv(S1,U)
    DinvS2 = Dinv(S2,U)

    pion_propagator = 0.0_dp
    do x = 1, Lx
       do y = 1, Ly
          pion_propagator(x) = pion_propagator(x) + &
               abs(DinvS1(1,x,y))**2 + &
               abs(DinvS1(2,x,y))**2 + &
               abs(DinvS2(1,x,y))**2 + &
               abs(DinvS2(2,x,y))**2 
       end do
    end do
    
  end function pion_propagator
  
  function staples(u,x,mu)

    complex(dp) :: staples
    complex(dp), dimension(:,:,:), intent(in) :: u
    integer(i4), intent(in) :: x(2), mu
    integer(i4), dimension(2) :: x2, x3, x4, x5, x6

    integer(i4) :: nu
    
    if ( mu == 1 ) then
       nu = 2
    elseif( mu == 2)then
       nu = 1
    end if

    x2 = ipf(x,nu)
    x3 = ipf(x,mu)
    x4 = imf(x,nu)
    x5 = x4
    x6 = imf(x3,nu)
    
    staples = u(nu,x(1),x(2)) * u(mu,x2(1), x2(2)) * conjg( u(nu,x3(1), x3(2)) ) + &
         conjg( u(nu,x4(1),x4(2)) ) * u(mu,x5(1), x5(2)) * u(nu,x6(1), x6(2))
    
  end function staples
  
  function DS(uold, unew, beta,stp)

    real(dp) :: DS,beta
    complex(dp) :: uold, unew, stp

    DS = -beta * real( (unew - uold) * conjg(stp) )
    
  end function DS
  
  function plaquette(u,x)

    complex(dp) :: plaquette
    
    complex(dp), dimension(:,:,:), intent(in) :: u
    integer(i4), intent(in) :: x(2)
    integer(i4), dimension(2) :: x2, x3


    x2 = ipf(x,1)
    x3 = ipf(x,2)
    
    plaquette = U(1,x(1),x(2)) * U(2,x2(1),x2(2)) * &
          conjg(U(1,x3(1),x3(2))) * conjg(U(2,x(1),x(2)))
    
  end function plaquette

  function action(u)
    use parameters, only : Lx, Ly
    real(dp) :: action
    complex(dp), dimension(2,Lx,Ly), intent(in) :: u
    integer(i4) :: x,y

    action = 0.0_dp
    
    do x = 1, Lx
       do y = 1, Ly
          action = action + real(plaquette(u,[x,y]))
       end do
    end do

    action = action / (Lx*Ly)
    
  end function action

  
  function action2(u,beta)

    use parameters, only : Lx, Ly
    real(dp) :: action2
    complex(dp), dimension(:,:,:), intent(in) :: u
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
    complex(dp), dimension(:,:,:), intent(in) :: u
    real(dp), intent(out) :: top(size(u(1,:,1)),size(u(1,1,:))) 
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

  
end module dynamics
