! this program calculates the dynamics of the spin-boson model
! using the dynamics fitering method and the RK4 propagator
! Authors: Tianchu Li and Qiang Shi.
!=========================================================
! 1. Main reference: Q. Shi, L.-P. Chen, G.-J. Nan, R.-X. Xu, and Y.-J. Yan, 
!    "Efficient hierarchical Liouville space propagator to quantum dissipative dynamics",
!    J. Chem. Phys. 130, 164518 (2009).

! 2. The filtered reduced density operator (RDO) and auxillary density 
!    operators (ADOs) are stored in a large array, operations are 
!    realized by using heap sort.

! 3. This is the spin-boson model, parameters: 
!    eps = (E_1-E_2)/2; vcoupling is the coupling constant

! 4. The harmonic bath uses the Debye-Drude spectral density, with parameters
!    J(\omega) = \frac{\eta \omega \omega_c}{\omega^2 + \omega_c^2}
!
! 5. passed using the intel fortran compiler ifort
!=========================================================

!                           DECLERATIONS
!                           ++++++++++++
module params
  implicit none
!=======================================================
! Physical constants
! ******************
  complex(8), parameter :: eye = (0.d0,1.d0)
  real(8), parameter    :: pi = 3.1415926535897932d0
  real(8), parameter    :: hbar = 1.d0

! model parameters
  real(8) :: beta = 1.d0, omega_c = 1.d0, eta = 0.5d0, &
             vcoupling = 1.d0, eps= 0.d0

! simulation parameters
  real(8), parameter    :: dt = 5.d-2
! cut-off parameters for dynamic filtering, check for convergence
  real(8), parameter    :: rcut1 = 1.d-7, rcut2 = 1.d-7, rcut3 = 1.d-9
  integer, parameter    :: nrmax = 3000000, neq =0, nsteps = 1000
! the largest truncating level :sum_nl for L-truncation, nb is basis size
  integer, parameter    :: sum_nl = 50
  integer, parameter    :: nb     = 50
! K = nk-1 is the number of Matsubara terms
  integer, parameter    :: nk = 3 
  integer, parameter    :: nlevel = nk
! number of system states
  integer, parameter    :: ndvr = 2

! variables for the HEOM, k0 is for the \delta function approximation,
! norm is the normalization factor, sdvr = diagonals of \simga_z
  real(8)  delta, k0, norm2fac, norm(nlevel)
  real(8)  hsys(ndvr,ndvr), sdvr(ndvr)

! for the exponential decomposition of the bath correlation function
! C(t) = \sum_k {[gval(k) + eye*gval_I(k)]*exp(-wval(k)*t)}
  real(8)    wval(nlevel), gval(nlevel)
  complex(8) gval_I(nlevel)

! data type to store the RDO/ADOs
! nl is the shifted index of the RDO/ADOs
! rval is the density matrix of the RDO/ADOs
  type  density                         
    integer nl(nlevel)           
    complex(8) rval(ndvr,ndvr)
  end type density

end module params

!=======================================================
! The main program
program main
  use params
  implicit none
  type(density) rhoall(nrmax)
  real(8)  t1
  integer i, istep, nrho, ios
  external delta_all

! initialize coefficents
  call setup_freqs_debye
!======================================
! -----------------------------------------------------------
! setup the initial RDO, which the first one in the array rhoall
  nrho = 1
  do i = 1, nlevel
      rhoall(1)%nl(i) = 1
  end do
  rhoall(1)%rval = 0.d0
  rhoall(1)%rval(1,1) = 1.d0

! -------------------------------------------------------------
! equilibrate on the donor state, set neq = 0 if not needed
  delta = 0.d0
  call setup_hsys_sb

  do istep = 1, neq
    t1 = (istep-1)*dt
    call rk4new(nrho,rhoall,dt,t1,delta_all)
    call prunewave(rhoall,nrho,rcut1)
  end do

! --------------------------------------------------------------
! propagation
  delta = vcoupling
  call setup_hsys_sb

  open(unit=23, file="./pop1.dat", iostat=ios, status="unknown", action="write")
  if ( ios /= 0 ) stop "Error opening file pop1.dat"
  open(unit=24, file="./pop2.dat", iostat=ios, status="unknown", action="write")
  if ( ios /= 0 ) stop "Error opening file pop2.dat"

  do istep = 1, nsteps
    if (mod(istep,100).eq.0) write(*,*) "istep = ", istep, "nrho =", nrho
    t1 = (istep-1)*dt
    write(23,*) t1, dreal(rhoall(1)%rval(1,1)), &
                                aimag(rhoall(1)%rval(1,1))
    write(24,*) t1, dreal(rhoall(1)%rval(2,2)),  &
                                aimag(rhoall(1)%rval(2,2))
    call rk4new(nrho,rhoall,dt,t1,delta_all)
    call prunewave(rhoall, nrho,rcut1)
  end do

  stop
end program

!==============================================
! setup the spin-boson Hamiltonian
subroutine setup_hsys_sb
  use params,only : eps,hsys,delta,sdvr
  implicit none

  hsys = 0.d0
  hsys(1,1) = eps;       hsys(2,2) = -eps
  hsys(1,2) = delta;     hsys(2,1) = delta
  sdvr(1) = 1.d0;        sdvr(2) = -1.d0

  return
end subroutine setup_hsys_sb


!=============================================
! the main subroutine for the RK4 propagation
subroutine rk4new(nwave,wave,dt,t_current,derivs)
  use params, only : density, nrmax
  implicit none
  integer, intent(in) :: nwave ! size of the input array
  type(density), intent(inout) :: wave(nrmax) !in : wave(time), out: wave(time+dt)
  real(8), intent(in) :: dt, t_current !time step and current time
  external derivs  !external subroutine to calculate \frac{d\rho}{dt}

  real(8)  dt2, dt6
  type(density) :: wave0(nrmax), dwave1(nrmax), dwave2(nrmax), dwave3(nrmax)
  integer  nwave0, ndwave1, ndwave2, ndwave3

  dt2=dt*0.5d0;  dt6=dt/6.d0

  !dwave1 = dwave/dt(t,wave(t)) ********
  call derivs(dwave1,wave,ndwave1,nwave,t_current)

  !wave0 = wave + dt2 * dwave1
  call addwave(wave0,wave,dwave1,nwave0,nwave,ndwave1,dt2)

  !dwave2 = dwave/dt(t+dt/2, rho(t)+dt/2*K1) *******
  call derivs(dwave2,wave0,ndwave2,nwave0,t_current+0.5d0*dt)

  !wave0 = wave + dt2 * dwave2
  call addwave(wave0,wave,dwave2,nwave0,nwave,ndwave2,dt2)

  !dwave3 = dwave/dt(t+dt/2, wave(t)+dt/2*K2) *******
  call derivs(dwave3,wave0,ndwave3,nwave0,t_current+0.5d0*dt)

  !wave0 = wave + dt * dwave3
  call addwave(wave0,wave,dwave3,nwave0,nwave,ndwave3,dt)

  !dwave3 = dwave2 + dwave3
  call addwave(dwave3,dwave3,dwave2,ndwave3,ndwave3,ndwave2,1.d0)

  !dwave2 = dwave/dt(t+dt, wave(t)+dt*K3)*********
  call derivs(dwave2,wave0,ndwave2,nwave0,t_current+dt)

  !wave = wave + dt6*(dwave1 + 2.d0*dwave3 + dwave2)
  call addwave(wave,wave,dwave1,nwave,nwave,ndwave1,dt6)
  call addwave(wave,wave,dwave2,nwave,nwave,ndwave2,dt6)
  call addwave(wave,wave,dwave3,nwave,nwave,ndwave3,2.d0*dt6)

  return
end subroutine rk4new

!=========================================
! calculates \frac{d\rho}{dt}, and store in an array
subroutine delta_all(drho,rho,ndrho,nrho,t1)
  use params
  implicit none

  type(density), intent(in)  :: rho(nrmax)
  type(density), intent(out) :: drho(nrmax)
  real(8),       intent(in)  :: t1
  integer,       intent(in)  :: nrho
  integer,       intent(out) :: ndrho

  complex(8) rhotmp(ndvr,ndvr), rhotmp_p(ndvr,ndvr), rhotmp_dc(ndvr,ndvr)
  real(8)    omegtmp
  integer    i, j, k, ii, jj

! ii labels the current size of the array
  ii = 0

  do i = 1, nrho
! The commutator and anticommutators with sigmaz, will be used later
    do j = 1, ndvr
      do k = 1 , ndvr
        rhotmp_p(j,k) = (sdvr(j) + sdvr(k)) * rho(i)%rval(j,k)
      end do
    end do

    do j = 1, ndvr
      do k = 1 , ndvr
        rhotmp(j,k) = (sdvr(j) - sdvr(k)) * rho(i)%rval(j,k)
      end do
    end do

! this is the Ishizaki-Tanimura delta-function term
    do j= 1, ndvr
      do k= 1, ndvr
        rhotmp_dc(j,k) = ((sdvr(j)**2+sdvr(k)**2)-2*sdvr(j)*sdvr(k)) &
                        * rho(i)%rval(j,k)
      end do
    end do

! diagonal terms
    omegtmp = 0.d0
    jj = 0.d0
    do j = 1, nlevel
      omegtmp = omegtmp + (rho(i)%nl(j)-1) * wval(j)
      jj = jj + (rho(i)%nl(j) - 1)  
    end do

    ii = ii+1
    if (ii.ge.nrmax) stop "too many elements!"
    drho(ii)%nl   = rho(i)%nl
    drho(ii)%rval = - omegtmp*rho(i)%rval       - eye*(matmul(hsys,rho(i)%rval) &
                    - matmul(rho(i)%rval,hsys)) - k0*rhotmp_dc

! hierachy terms
    do j = 1, nlevel
! from \rho_{n} to \rho_{n+1}
      if (rho(i)%nl(j).lt.nb .and. jj.lt.sum_nl)  then
        ii = ii + 1
        drho(ii)%nl    = rho(i)%nl
        drho(ii)%nl(j) = rho(i)%nl(j) + 1
        drho(ii)%rval  = gval_I(j)/norm(j)*sqrt(1.d0*rho(i)%nl(j))*rhotmp_p
        drho(ii)%rval  = drho(ii)%rval - eye*gval(j)/norm(j) &
                        *sqrt(1.d0*rho(i)%nl(j))*rhotmp
      end if

! from \rho_{n+1} to \rho_{n}
      if (rho(i)%nl(j).gt.1)  then
        ii = ii + 1
        drho(ii)%nl    = rho(i)%nl
        drho(ii)%nl(j) = rho(i)%nl(j) - 1
        drho(ii)%rval  = - eye*norm(j)*sqrt(rho(i)%nl(j)-1.d0)*rhotmp
      end if
    end do !j = 1,nlevel

  end do ! i = 1, nrho
  ndrho = ii
! sort and merge term with the same indices
  call mergewave(drho,ndrho)

  return
end subroutine delta_all

!===============================================
! add two different \rho arrays
subroutine addwave(wave1,wave,dwave,nwave1,nwave,ndwave,coef2)
  use params
  implicit none
  type(density) wave1(nrmax), wave(nrmax), dwave(nrmax)
  real(8), intent(in)  :: coef2
  integer nwave1, nwave, ndwave

  integer i, j
  logical eqcmp

  if (nwave+ndwave.ge.nrmax) stop "too many elements!"
  do i = 1, nwave
    wave1(i) = wave(i)
  end do
  do i = 1, ndwave
    wave1(nwave+i)%nl   = dwave(i)%nl
    wave1(nwave+i)%rval = dwave(i)%rval*coef2
  end do
  nwave1 = nwave + ndwave

! sort it using heap sort
  call hpsort4rho(nwave1,wave1)

! merge all same elements
  j = 1
  do i = 2, nwave1
    if (eqcmp(wave1(i),wave1(j))) then
      wave1(j)%rval = wave1(j)%rval + wave1(i)%rval
    else
      j = j + 1
      if (j.lt.i) wave1(j) = wave1(i)
    end if
  end do
  nwave1 = j

  call prunewave(wave1, nwave1,rcut2)

  return
end subroutine addwave

!================================================
! discard small ADOs for dynamic filtering
subroutine prunewave(wave1,nwave1, rsmall)
  use params
  implicit none
  type(density), intent(inout) :: wave1(nrmax)
  integer, intent(inout) :: nwave1
  real(8), intent(in) :: rsmall

  integer i, j

  ! prune small elements
  j = 0
  do i = 1, nwave1
    if (maxval(abs(wave1(i)%rval)) .gt. rsmall) then
      j = j + 1
      if (j.lt.i) wave1(j) = wave1(i)
    end if
  end do
  nwave1 = j

  return
end subroutine prunewave

!====================================================
! this subroutine sorts the array and 
! merges all terms with the same indices
subroutine mergewave(wave1,nwave1)
  use params
  implicit none
  type(density), intent(inout) :: wave1(nrmax)
  integer, intent(inout) :: nwave1
  logical eqcmp

  integer i, j
! sort the sum terms using heap sort
  call hpsort4rho(nwave1,wave1)

! merge all same elements
  j = 1
  do i = 2, nwave1
    if (eqcmp(wave1(i),wave1(j))) then
      wave1(j)%rval = wave1(j)%rval + wave1(i)%rval
    else
      j = j + 1
      if (j.lt.i) wave1(j) = wave1(i)
    end if
  end do
  nwave1 = j

  return
end subroutine mergewave

!===========================================================
! compare two ADO indices, returns true if they are the same
function eqcmp(r1,r2)
  use params
  implicit none
  type(density), intent(in) :: r1, r2
  logical eqcmp
  integer i

  eqcmp = .true.
  do i = 1, nlevel
    if (r1%nl(i) .ne. r2%nl(i)) then
       eqcmp = .false.
       exit
    end if
  end do

  return
end function

!====================================================
! compares two ADO indices, returns true if r1 < r2
! we arrange it such that (1,1,1) is the smallest
function ltcmp(r1,r2)
  use params
  implicit none
  type(density), intent(in) :: r1, r2
  logical ltcmp
  integer i

  ltcmp = .false.
  do i = 1, nlevel
    if (r1%nl(i) .lt. r2%nl(i)) then
      ltcmp = .true.
      exit
    else if (r1%nl(i) .gt. r2%nl(i)) then
      ltcmp = .false.
      exit
    end if
  end do

  return
end function

!==========================================
! the heap sort subroutine for the RDO/ADOs
! adapted from "Numerical Recipes in Fortran"
subroutine hpsort4rho(n,ra)
  use params
  implicit none
  integer n
  type(density) ra(n)
  integer i,ir,j,l
  type(density) rra
  logical ltcmp
  if (n.lt.2) return
  l=n/2+1
  ir=n

  do while( .true. )
    if(l.gt.1)then
      l=l-1
      rra=ra(l)
    else
      rra=ra(ir)
      ra(ir)=ra(1)
      ir=ir-1
      if(ir.eq.1)then
        ra(1)=rra
        return
      endif
    endif
    i=l
    j=l+l
    do while(j.le.ir)
      if(j.lt.ir)then
        if(ltcmp(ra(j),ra(j+1)))j=j+1
      endif
      if(ltcmp(rra,ra(j)))then
        ra(i)=ra(j)
        i=j
        j=j+j
      else
        j=ir+1
      endif
    end do
    ra(i)=rra
  end do

  return
end subroutine hpsort4rho

!==========================================
! obtain the parameters for the bath correlation function 
! decomposition, using the Matsubara frequencies
subroutine setup_freqs_debye
  use params
  implicit none
  real(8) omegtmp, a0coef
  integer i

  gval(1) = eta*omega_c/2.d0/tan(0.5d0*beta*omega_c)
  gval_I = 0.d0
  gval_I(1) = -eta*omega_c/2.d0
  wval(1) = omega_c
  do i = 1, nk-1 !note k starts from zero
    omegtmp = 2.d0*pi*i/beta
    gval(i+1) = - 2.d0*eta*omega_c/beta*omegtmp/(omega_c**2 - omegtmp**2)
    wval(i+1) = omegtmp
  end do

! this is the dephasing term
  k0 = (eta/beta-gval(1))/omega_c
  do i = 2, nk
    k0 = k0 - gval(i)/wval(i)
  end do

! normalization factor used for filtering
  a0coef = 0.d0
  do i = 1, nlevel
    norm(i) = sqrt(abs(gval(i)))
    a0coef   = a0coef + gval(i)
  end do
! an alternative normalization factor, not used here
  norm2fac = sqrt(abs(a0coef))

  return
end subroutine setup_freqs_debye
