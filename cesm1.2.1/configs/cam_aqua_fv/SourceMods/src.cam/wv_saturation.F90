module wv_saturation

!--------------------------------------------------------------------!
! Module Overview:                                                   !
!                                                                    !
! This module provides an interface to wv_sat_methods, providing     !
! saturation vapor pressure and related calculations to CAM.         !
!                                                                    !
! The original wv_saturation codes were introduced by J. J. Hack,    !
! February 1990. The code has been extensively rewritten since then, !
! including a total refactoring in Summer 2012.                      !
!                                                                    !
!--------------------------------------------------------------------!
! Methods:                                                           !
!                                                                    !
! Pure water/ice saturation vapor pressures are calculated on the    !
! fly, with the specific method determined by a runtime option.      !
! Mixed phase SVP is interpolated from the internal table, estbl,    !
! which is created during initialization.                            !
!                                                                    !
! The default method for calculating SVP is determined by a namelist !
! option, and used whenever svp_water/ice or qsat are called.        !
!                                                                    !
!--------------------------------------------------------------------!

use shr_kind_mod, only: r8 => shr_kind_r8
use physconst,    only: epsilo, &
                        latvap, &
                        latice, &
                        rh2o,   &
                        cpair,  &
                        cpliq,  &
                        tmelt,  &
                        h2otrip, &
                        rair,   &
                        cpwv

use wv_sat_methods, only: &
     svp_to_qsat => wv_sat_svp_to_qsat, &
     svp_to_qmmr => wv_sat_svp_to_qmmr

implicit none
private
save

! Public interfaces
! Namelist, initialization, finalization
public wv_sat_readnl
public wv_sat_init
public wv_sat_final

! Saturation vapor pressure calculations
public svp_water
public svp_ice
  
! Mixed phase (water + ice) saturation vapor pressure table lookup
public estblf

public svp_to_qsat
public svp_to_qmmr

! Subroutines that return both SVP and humidity
! Optional arguments do temperature derivatives
public qsat           ! Mixed phase
public qsat_water     ! SVP over water only
public qsat_ice       ! SVP over ice only

! Wet bulb temperature solver
public findsp_vc

! Unusual/non-default methods.
! Zhang-McFarlane scheme uses mass mixing ratio rather than
! specific humidity.
public qmmr

! Data

! This value is slightly high, but it seems to be the value for the
! steam point of water originally (and most frequently) used in the
! Goff & Gratch scheme.
real(r8), parameter :: tboil = 373.16_r8

! Table of saturation vapor pressure values (estbl) from tmin to
! tmax+1 Kelvin, in one degree increments.  ttrice defines the
! transition region, estbl contains a combination of ice & water
! values.
! Make these public parameters in case another module wants to see the
! extent of the table.
  real(r8), public, parameter :: tmin = 127.16_r8
!  real(r8), public, parameter :: tmax = 375.16_r8
  real(r8), public, parameter :: tmax = 500.16_r8

  real(r8), parameter :: ttrice = 20.00_r8  ! transition range from es over H2O to es over ice

  integer :: plenest                             ! length of estbl
  real(r8), allocatable :: estbl(:)              ! table values of saturation vapor pressure

  real(r8) :: omeps      ! 1.0_r8 - epsilo

  real(r8) :: c3         ! parameter used by findsp

  ! Set coefficients for polynomial approximation of difference
  ! between saturation vapor press over water and saturation pressure
  ! over ice for -ttrice < t < 0 (degrees C). NOTE: polynomial is
  ! valid in the range -40 < t < 0 (degrees C).
  real(r8) :: pcf(5) = (/ &
       5.04469588506e-01_r8, &
       -5.47288442819e+00_r8, &
       -3.67471858735e-01_r8, &
       -8.95963532403e-03_r8, &
       -7.78053686625e-05_r8 /)

!   --- Degree 6 approximation ---
!  real(r8) :: pcf(6) = (/ &
!       7.63285250063e-02, &
!       5.86048427932e+00, &
!       4.38660831780e-01, &
!       1.37898276415e-02, &
!       2.14444472424e-04, &
!       1.36639103771e-06 /)

contains

!---------------------------------------------------------------------
! ADMINISTRATIVE FUNCTIONS
!---------------------------------------------------------------------

subroutine wv_sat_readnl(nlfile)
  !------------------------------------------------------------------!
  ! Purpose:                                                         !
  !   Get runtime options for wv_saturation.                         !
  !------------------------------------------------------------------!

  use wv_sat_methods, only: wv_sat_get_scheme_idx, &
                            wv_sat_valid_idx, &
                            wv_sat_set_default

  use spmd_utils,      only: masterproc
  use namelist_utils,  only: find_group_name
  use units,           only: getunit, freeunit
  use mpishorthand
  use abortutils,      only: endrun

  character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input
   
  ! Local variables
  integer :: unitn, ierr

  character(len=32) :: wv_sat_scheme = "GoffGratch"

  character(len=*), parameter :: subname = 'wv_sat_readnl'

  namelist /wv_sat_nl/ wv_sat_scheme
  !-----------------------------------------------------------------------------

  if (masterproc) then
     unitn = getunit()
     open( unitn, file=trim(nlfile), status='old' )
     call find_group_name(unitn, 'wv_sat_nl', status=ierr)
     if (ierr == 0) then
        read(unitn, wv_sat_nl, iostat=ierr)
        if (ierr /= 0) then
           call endrun(subname // ':: ERROR reading namelist')
           return
        end if
     end if
     close(unitn)
     call freeunit(unitn)

  end if

#ifdef SPMD
  call mpibcast(wv_sat_scheme, len(wv_sat_scheme) , mpichar, 0, mpicom)
#endif

  if (.not. wv_sat_set_default(wv_sat_scheme)) then
     call endrun('wv_sat_readnl :: Invalid wv_sat_scheme.')
     return
  end if

end subroutine wv_sat_readnl

subroutine wv_sat_init
  !------------------------------------------------------------------!
  ! Purpose:                                                         !
  !   Initialize module (e.g. setting parameters, initializing the   !
  !   SVP lookup table).                                             !
  !------------------------------------------------------------------!

  use wv_sat_methods, only: wv_sat_methods_init, &
                            wv_sat_get_scheme_idx, &
                            wv_sat_valid_idx
  use spmd_utils,     only: masterproc
  use cam_logfile,    only: iulog
  use abortutils,     only: endrun
  use shr_assert_mod, only: shr_assert_in_domain
  use error_messages, only: handle_errmsg

  integer :: status

  ! For wv_sat_methods error reporting.
  character(len=256) :: errstring

  ! For generating internal SVP table.
  real(r8) :: t         ! Temperature
  integer  :: i         ! Increment counter

  ! Precalculated because so frequently used.
  omeps  = 1.0_r8 - epsilo

  ! Transition range method is only valid for transition temperatures at:
  ! -40 deg C < T < 0 deg C
  call shr_assert_in_domain(ttrice, ge=0._r8, le=40._r8, varname="ttrice",&
       msg="wv_sat_init: Invalid transition temperature range.")

! This parameter uses a hardcoded 287.04_r8?
  c3 = rair*(7.5_r8*log(10._r8))/cpair

! Init "methods" module containing actual SVP formulae.

  call wv_sat_methods_init(r8, tmelt, h2otrip, tboil, ttrice, &
       epsilo, errstring)

  call handle_errmsg(errstring, subname="wv_sat_methods_init")

  ! Add two to make the table slightly too big, just in case.
  plenest = ceiling(tmax-tmin) + 2

  ! Allocate SVP table.
  allocate(estbl(plenest), stat=status)
  if (status /= 0) then 
     call endrun('wv_sat_init :: ERROR allocating saturation vapor pressure table')
     return
  end if

  do i = 1, plenest
    estbl(i) = svp_trans(tmin + real(i-1,r8))
  end do

  if (masterproc) then
     write(iulog,*)' *** SATURATION VAPOR PRESSURE TABLE COMPLETED ***'
  end if

end subroutine wv_sat_init

subroutine wv_sat_final
  !------------------------------------------------------------------!
  ! Purpose:                                                         !
  !   Deallocate global variables in module.                         !
  !------------------------------------------------------------------!
  use abortutils,   only: endrun

  integer :: status

  if (allocated(estbl)) then

     deallocate(estbl, stat=status)

     if (status /= 0) then
        call endrun('wv_sat_final :: ERROR deallocating table')
        return
     end if

  end if

end subroutine wv_sat_final

!---------------------------------------------------------------------
! DEFAULT SVP FUNCTIONS
!---------------------------------------------------------------------

! Compute saturation vapor pressure over water
elemental function svp_water(t) result(es)

  use wv_sat_methods, only: &
       wv_sat_svp_water

  real(r8), intent(in) :: t ! Temperature (K)
  real(r8) :: es            ! SVP (Pa)

  es = wv_sat_svp_water(T)

end function svp_water

! Compute saturation vapor pressure over ice
elemental function svp_ice(t) result(es)

  use wv_sat_methods, only: &
       wv_sat_svp_ice

  real(r8), intent(in) :: t ! Temperature (K)
  real(r8) :: es            ! SVP (Pa)

  es = wv_sat_svp_ice(T)

end function svp_ice

! Compute saturation vapor pressure with an ice-water transition
elemental function svp_trans(t) result(es)

  use wv_sat_methods, only: &
       wv_sat_svp_trans

  real(r8), intent(in) :: t ! Temperature (K)
  real(r8) :: es            ! SVP (Pa)

  es = wv_sat_svp_trans(T)

end function svp_trans

!---------------------------------------------------------------------
! UTILITIES
!---------------------------------------------------------------------

! Does linear interpolation from nearest values found
! in the table (estbl).
elemental function estblf(t) result(es)

  real(r8), intent(in) :: t ! Temperature 
  real(r8) :: es            ! SVP (Pa)

  integer  :: i         ! Index for t in the table
  real(r8) :: t_tmp     ! intermediate temperature for es look-up

  real(r8) :: weight ! Weight for interpolation

  t_tmp = max(min(t,tmax)-tmin, 0._r8)   ! Number of table entries above tmin
  i = int(t_tmp) + 1                     ! Corresponding index.
  weight = t_tmp - aint(t_tmp, r8)       ! Fractional part of t_tmp (for interpolation).
  es = (1._r8 - weight)*estbl(i) + weight*estbl(i+1)

end function estblf

! Get enthalpy based only on temperature
! and specific humidity.
elemental function tq_enthalpy(t, q, hltalt) result(enthalpy)

  real(r8), intent(in) :: t      ! Temperature
  real(r8), intent(in) :: q      ! Specific humidity
  real(r8), intent(in) :: hltalt ! Modified hlat for T derivatives

  real(r8) :: enthalpy

  enthalpy = ((1-q)*cpair + q*cpwv) * t + hltalt * q
  
end function tq_enthalpy

!---------------------------------------------------------------------
! LATENT HEAT OF VAPORIZATION CORRECTIONS
!---------------------------------------------------------------------

elemental subroutine no_ip_hltalt(t, hltalt)
  !------------------------------------------------------------------!
  ! Purpose:                                                         !
  !   Calculate latent heat of vaporization of pure liquid water at  !
  !   a given temperature.                                           !
  !------------------------------------------------------------------!

  ! Inputs
  real(r8), intent(in) :: t        ! Temperature
  ! Outputs
  real(r8), intent(out) :: hltalt  ! Appropriately modified hlat

  hltalt = latvap

  ! Account for change of latvap with t above freezing where
  ! constant slope is given by -2369 j/(kg c) = cpv - cw
  if (t >= tmelt) then
     hltalt = hltalt - 2369.0_r8*(t-tmelt)
  end if

end subroutine no_ip_hltalt

elemental subroutine calc_hltalt(t, hltalt, tterm)
  !------------------------------------------------------------------!
  ! Purpose:                                                         !
  !   Calculate latent heat of vaporization of water at a given      !
  !   temperature, taking into account the ice phase if temperature  !
  !   is below freezing.                                             !
  !   Optional argument also calculates a term used to calculate     !
  !   d(es)/dT within the water-ice transition range.                !
  !------------------------------------------------------------------!

  ! Inputs
  real(r8), intent(in) :: t        ! Temperature
  ! Outputs
  real(r8), intent(out) :: hltalt  ! Appropriately modified hlat
  ! Term to account for d(es)/dT in transition region.
  real(r8), intent(out), optional :: tterm

  ! Local variables
  real(r8) :: tc      ! Temperature in degrees C
  real(r8) :: weight  ! Weight for es transition from water to ice
  ! Loop iterator
  integer :: i

  if (present(tterm)) tterm = 0.0_r8

  call no_ip_hltalt(t,hltalt)
  if (t < tmelt) then
     ! Weighting of hlat accounts for transition from water to ice.
     tc = t - tmelt

     if (tc >= -ttrice) then
        weight = -tc/ttrice

        ! polynomial expression approximates difference between es
        ! over water and es over ice from 0 to -ttrice (C) (max of
        ! ttrice is 40): required for accurate estimate of es
        ! derivative in transition range from ice to water
        if (present(tterm)) then
           do i = size(pcf), 1, -1
              tterm = pcf(i) + tc*tterm
           end do
           tterm = tterm/ttrice
        end if

     else
        weight = 1.0_r8
     end if

     hltalt = hltalt + weight*latice

  end if

end subroutine calc_hltalt

!---------------------------------------------------------------------
! OPTIONAL OUTPUTS
!---------------------------------------------------------------------

! Temperature derivative outputs, for qsat_*
elemental subroutine deriv_outputs(t, p, es, qs, hltalt, tterm, &
     gam, dqsdt)

  ! Inputs
  real(r8), intent(in) :: t      ! Temperature
  real(r8), intent(in) :: p      ! Pressure
  real(r8), intent(in) :: es     ! Saturation vapor pressure
  real(r8), intent(in) :: qs     ! Saturation specific humidity
  real(r8), intent(in) :: hltalt ! Modified latent heat
  real(r8), intent(in) :: tterm  ! Extra term for d(es)/dT in
                                 ! transition region.

  ! Outputs
  real(r8), intent(out), optional :: gam      ! (hltalt/cpair)*(d(qs)/dt)
  real(r8), intent(out), optional :: dqsdt    ! (d(qs)/dt)

  ! Local variables
  real(r8) :: desdt        ! d(es)/dt
  real(r8) :: dqsdt_loc    ! local copy of dqsdt
  real(r8) :: cpmix        ! cp of saturated env heat capacity

  if (qs == 1.0_r8) then
     dqsdt_loc = 0._r8
  else
     desdt = hltalt*es/(rh2o*t*t) + tterm
     dqsdt_loc = qs*p*desdt/(es*(p-omeps*es))
  end if

  cpmix = ((1-qs)*cpair + qs*cpwv)

  if (present(dqsdt)) dqsdt = dqsdt_loc
  if (present(gam))   gam   = dqsdt_loc * (hltalt/cpmix)

end subroutine deriv_outputs

!---------------------------------------------------------------------
! QSAT (SPECIFIC HUMIDITY) PROCEDURES
!---------------------------------------------------------------------

elemental subroutine qmmr(t, p, es, qm)
  !------------------------------------------------------------------!
  ! Purpose:                                                         !
  !     Provide saturation mass mixing ratio.                        !
  !                                                                  !
  ! Note that qmmr is a ratio over dry air, whereas qsat is a ratio  !
  ! over total mass. These differ mainly at low pressures, where     !
  ! SVP may exceed the actual pressure, and therefore qmmr can blow  !
  ! up.                                                              !
  !------------------------------------------------------------------!
  use wv_sat_methods, only: &
       wv_sat_svp_water


  ! Inputs
  real(r8), intent(in) :: t    ! Temperature
  real(r8), intent(in) :: p    ! Pressure
  ! Outputs
  real(r8), intent(out) :: es  ! Saturation vapor pressure
  real(r8), intent(out) :: qm  ! Saturation mass mixing ratio
                               ! (vapor mass over dry mass)

  es = wv_sat_svp_water(t)

  qm = svp_to_qmmr(es, p)

  ! Ensures returned es is consistent with limiters on qmmr.
  es = min(es, p)

end subroutine qmmr

elemental subroutine qsat(t, p, es, qs, gam, dqsdt, enthalpy)
  !------------------------------------------------------------------!
  ! Purpose:                                                         !
  !   Look up and return saturation vapor pressure from precomputed  !
  !   table, then calculate and return saturation specific humidity. !
  !   Optionally return various temperature derivatives or enthalpy  !
  !   at saturation.                                                 !
  !------------------------------------------------------------------!

  ! Inputs
  real(r8), intent(in) :: t    ! Temperature
  real(r8), intent(in) :: p    ! Pressure
  ! Outputs
  real(r8), intent(out) :: es  ! Saturation vapor pressure
  real(r8), intent(out) :: qs  ! Saturation specific humidity

  real(r8), intent(out), optional :: gam    ! (l/cpair)*(d(qs)/dt)
  real(r8), intent(out), optional :: dqsdt  ! (d(qs)/dt)
  real(r8), intent(out), optional :: enthalpy ! cpmix*t + hltalt*q

  ! Local variables
  real(r8) :: hltalt       ! Modified latent heat for T derivatives
  real(r8) :: tterm        ! Account for d(es)/dT in transition region

  es = estblf(t)

  qs = svp_to_qsat(es, p)

  ! Ensures returned es is consistent with limiters on qs.
  es = min(es, p)

  ! Calculate optional arguments.
  if (present(gam) .or. present(dqsdt) .or. present(enthalpy)) then

     ! "generalized" analytic expression for t derivative of es
     ! accurate to within 1 percent for 173.16 < t < 373.16
     call calc_hltalt(t, hltalt, tterm)

     if (present(enthalpy)) enthalpy = tq_enthalpy(t, qs, hltalt)

     call deriv_outputs(t, p, es, qs, hltalt, tterm, &
          gam=gam, dqsdt=dqsdt)

  end if

end subroutine qsat

elemental subroutine qsat_water(t, p, es, qs, gam, dqsdt, enthalpy)
  !------------------------------------------------------------------!
  ! Purpose:                                                         !
  !   Calculate SVP over water at a given temperature, and then      !
  !   calculate and return saturation specific humidity.             !
  !   Optionally return various temperature derivatives or enthalpy  !
  !   at saturation.                                                 !
  !------------------------------------------------------------------!

  use wv_sat_methods, only: wv_sat_qsat_water

  ! Inputs
  real(r8), intent(in) :: t    ! Temperature
  real(r8), intent(in) :: p    ! Pressure
  ! Outputs
  real(r8), intent(out) :: es  ! Saturation vapor pressure
  real(r8), intent(out) :: qs  ! Saturation specific humidity

  real(r8), intent(out), optional :: gam    ! (l/cpair)*(d(qs)/dt)
  real(r8), intent(out), optional :: dqsdt  ! (d(qs)/dt)
  real(r8), intent(out), optional :: enthalpy ! cpmix*t + hltalt*q

  ! Local variables
  real(r8) :: hltalt       ! Modified latent heat for T derivatives

  call wv_sat_qsat_water(t, p, es, qs)

  if (present(gam) .or. present(dqsdt) .or. present(enthalpy)) then

     ! "generalized" analytic expression for t derivative of es
     ! accurate to within 1 percent for 173.16 < t < 373.16
     call no_ip_hltalt(t, hltalt)

     if (present(enthalpy)) enthalpy = tq_enthalpy(t, qs, hltalt)

     ! For pure water/ice transition term is 0.
     call deriv_outputs(t, p, es, qs, hltalt, 0._r8, &
          gam=gam, dqsdt=dqsdt)

  end if

end subroutine qsat_water

elemental subroutine qsat_ice(t, p, es, qs, gam, dqsdt, enthalpy)
  !------------------------------------------------------------------!
  ! Purpose:                                                         !
  !   Calculate SVP over ice at a given temperature, and then        !
  !   calculate and return saturation specific humidity.             !
  !   Optionally return various temperature derivatives or enthalpy  !
  !   at saturation.                                                 !
  !------------------------------------------------------------------!

  use wv_sat_methods, only: wv_sat_qsat_ice

  ! Inputs
  real(r8), intent(in) :: t    ! Temperature
  real(r8), intent(in) :: p    ! Pressure
  ! Outputs
  real(r8), intent(out) :: es  ! Saturation vapor pressure
  real(r8), intent(out) :: qs  ! Saturation specific humidity

  real(r8), intent(out), optional :: gam    ! (l/cpair)*(d(qs)/dt)
  real(r8), intent(out), optional :: dqsdt  ! (d(qs)/dt)
  real(r8), intent(out), optional :: enthalpy ! cpair*t + hltalt*q

  ! Local variables
  real(r8) :: hltalt       ! Modified latent heat for T derivatives

  call wv_sat_qsat_ice(t, p, es, qs)

  if (present(gam) .or. present(dqsdt) .or. present(enthalpy)) then

     ! For pure ice, just add latent heats.
     hltalt = latvap + latice

     if (present(enthalpy)) enthalpy = tq_enthalpy(t, qs, hltalt)

     ! For pure water/ice transition term is 0.
     call deriv_outputs(t, p, es, qs, hltalt, 0._r8, &
          gam=gam, dqsdt=dqsdt)

  end if

end subroutine qsat_ice

!---------------------------------------------------------------------
! FINDSP (WET BULB TEMPERATURE) PROCEDURES
!---------------------------------------------------------------------

subroutine findsp_vc(q, t, p, use_ice, tsp, qsp)

  use cam_logfile,  only: iulog
  use abortutils,   only: endrun
  ! Wrapper for findsp which is 1D and handles the output status.
  ! Changing findsp to elemental restricted debugging output.
  ! If that output is needed again, it's preferable *not* to copy findsp,
  ! but to change the existing version.

  ! input arguments
  real(r8), intent(in) :: q(:)        ! water vapor (kg/kg)
  real(r8), intent(in) :: t(:)        ! temperature (K)
  real(r8), intent(in) :: p(:)        ! pressure    (Pa)
  logical,  intent(in) :: use_ice     ! flag to include ice phase in calculations

  ! output arguments
  real(r8), intent(out) :: tsp(:)     ! saturation temp (K)
  real(r8), intent(out) :: qsp(:)     ! saturation mixing ratio (kg/kg)

  integer :: status(size(q))   ! flag representing state of output
                               ! 0 => Successful convergence
                               ! 1 => No calculation done: pressure or specific
                               !      humidity not within usable range
                               ! 2 => Run failed to converge
                               ! 4 => Temperature fell below minimum
                               ! 8 => Enthalpy not conserved

  integer :: l_found(size(q)) !on which of the iterations the successful q was found

  real(r8) :: dq(size(q), 64)
  real(r8) :: dt(size(q), 64)
  real(r8) :: g(size(q), 64)
  real(r8) :: dgdt(size(q), 64)
  real(r8) :: dtfg(size(q))
  real(r8) :: enin(size(q))
  real(r8) :: enout(size(q), 64)
  real(r8) :: t1(size(q), 64)
  real(r8) :: q1(size(q), 64)
  real(r8) :: cpmix(size(q), 64)
  real(r8) :: gam(size(q), 64)
  real(r8) :: dcpqdt(size(q), 64)
  real(r8) :: hltalt(size(q), 64)


  integer :: n, i, l

  n = size(q)
  l=64 !max number of search iterations in findsp

  do i = 1,n
   call findsp(q(i), t(i), p(i), use_ice, tsp(i), qsp(i), status(i), l_found(i), dq(i,:), dt(i,:), g(i,:), dgdt(i,:), dtfg(i), enin(i), enout(i,:), t1(i,:), q1(i,:), cpmix(i,:), gam(i,:), dcpqdt(i,:), hltalt(i,:))
  end do

  ! Currently, only 2 and 8 seem to be treated as fatal errors.
  do i = 1,n
     if (status(i) == 2) then
        write(iulog,*) ' findsp not converging at i = ', i
        write(iulog,*) ' t, q, p ', t(i), q(i), p(i)
        write(iulog,*) ' tsp, qsp ', tsp(i), qsp(i)
        write(iulog,*) ' enin:', enin(i), 'dtfg:', dtfg(i)
        do l = 1,l_found(i)
            write(iulog,*) ' dt:', dt(i,l), 'dq:', dq(i,l), 'g:', g(i,l), 'dgdt:', dgdt(i,l)
            write(iulog,*) ' cpmix:', cpmix(i,l), "enout:", enout(i,l)
            write(iulog,*) ' gam:', gam(i,l), 'tsp*(cpwv-cpair)*(1/hltalt)', dcpqdt(i,l), "hltalt:", hltalt(i,l)
        end do
        call endrun ('wv_saturation::FINDSP -- not converging')
     else if (status(i) == 8) then
        write(iulog,*) ' the enthalpy is not conserved at i = ', i
        write(iulog,*) ' t, q, p ', t(i), q(i), p(i)
        write(iulog,*) ' tsp, qsp ', tsp(i), qsp(i)
        write(iulog,*) ' enin:', enin(i), 'dtfg:', dtfg(i)
        do l = 1,l_found(i)
            write(iulog,*) ' dt:', dt(i,l), 'dq:', dq(i,l), 'g:', g(i,l), 'dgdt:', dgdt(i,l)
            write(iulog,*) ' cpmix:', cpmix(i,l), "enout:", enout(i,l)
            write(iulog,*) ' gam:', gam(i,l), 'tsp*(cpwv-cpair)*(1/hltalt)', dcpqdt(i,l), "hltalt:", hltalt(i,l)
        end do
        call endrun ('wv_saturation::FINDSP -- enthalpy is not conserved')
     endif
  end do

end subroutine findsp_vc

subroutine findsp (q, t, p, use_ice, tsp, qsp, status, l_found, dq, dt, g, dgdt, dtfg, enin, enout, t1, q1, cpmix, gam, dcpqdt, hltalt)
! some changes to increase number of outputs to help debug traces
!----------------------------------------------------------------------- 
!
! Purpose: 
!     find the wet bulb temperature for a given t and q
!     in a longitude height section
!     wet bulb temp is the temperature and spec humidity that is 
!     just saturated and has the same enthalpy
!     if q > qs(t) then tsp > t and qsp = qs(tsp) < q
!     if q < qs(t) then tsp < t and qsp = qs(tsp) > q
!
! Method: 
! a Newton method is used
! first guess uses an algorithm provided by John Petch from the UKMO
! we exclude points where the physical situation is unrealistic
! e.g. where the temperature is outside the range of validity for the
!      saturation vapor pressure, or where the water vapor pressure
!      exceeds the ambient pressure, or the saturation specific humidity is 
!      unrealistic
! 
! Author: P. Rasch
! 
!-----------------------------------------------------------------------
!
!     input arguments
!
  use cam_logfile,  only: iulog


  real(r8), intent(in) :: q        ! water vapor (kg/kg)
  real(r8), intent(in) :: t        ! temperature (K)
  real(r8), intent(in) :: p        ! pressure    (Pa)
  logical,  intent(in) :: use_ice  ! flag to include ice phase in calculations
!
! output arguments
!
  real(r8), intent(out) :: tsp      ! saturation temp (K)
  real(r8), intent(out) :: qsp      ! saturation mixing ratio (kg/kg)
  integer,  intent(out) :: status   ! flag representing state of output
                                    ! 0 => Successful convergence
                                    ! 1 => No calculation done: pressure or specific
                                    !      humidity not within usable range
                                    ! 2 => Run failed to converge
                                    ! 4 => Temperature fell below minimum
                                    ! 8 => Enthalpy not conserved
  integer, intent(out) :: l_found   !additional output, stage of iteration at which we stop / converge
  real(r8), intent(out) :: dt(64)       ! additional output, change in guessed t at the last guess
  real(r8), intent(out) :: dq(64)       ! additional output, change in guessed q at the last guess
  real(r8), intent(out) :: g(64)        ! additional output, gam at last guess
  real(r8), intent(out) :: dgdt(64)     ! additional output, change in g with t at the last guess
  real(r8), intent(out) :: dtfg     ! additional output, change in t first guessed
  real(r8), intent(out) :: enin     ! additional output, input enthalpy
  real(r8), intent(out) :: enout(64)    ! additional output, final enthalpy
  real(r8), intent(out) :: cpmix(64)    ! additional output, final cpmix
  real(r8), intent(out) :: t1(64)       ! additional output, loop guess t
  real(r8), intent(out) :: q1(64)       ! additional output, loop guess q
  real(r8), intent(out) :: gam(64)      ! additional output, change in sat spec. hum. wrt temperature (times hltalt/cpair)
  real(r8), intent(out) :: dcpqdt(64)   ! additional output, change in enthalpy wrt T due to cpmix changes
  real(r8), intent(out) :: hltalt(64)   ! additional output, L



!
! local variables
!
  integer, parameter :: iter = 64   ! max number of times to iterate the calculation
  integer :: l                      ! iterator

  real(r8) es                   ! sat. vapor pressure
!   real(r8) gam                  ! change in sat spec. hum. wrt temperature (times hltalt/cpair)
!   real(r8) dgdt                 ! work variable
!   real(r8) g                    ! work variable
!   real(r8) hltalt               ! lat. heat. of vap.
  real(r8) qs                   ! spec. hum. of water vapor

! work variables
!   real(r8) t1, q1 !, dt, dq
  real(r8) qvd
  real(r8) r1b, c1, c2
  real(r8), parameter :: dttol = 1.e-4_r8 ! the relative temp error tolerance required to quit the iteration
  real(r8), parameter :: dqtol = 1.e-4_r8 ! the relative moisture error tolerance required to quit the iteration
!   real(r8) enin, enout
!   real(r8) cpmix ! composition adjusted cp

  dtfg = 0._r8
  enout(:) = 0._r8
  dt(:) = 0._r8
  dq(:) = 0._r8
  g(:) = 0._r8
  dgdt(:) = 0._r8
  cpmix(:) = 0._r8
  t1(:) = 0._r8
  q1(:) = 0._r8
  gam(:) = 0._r8
  hltalt(:) = 0._r8


  ! Saturation specific humidity at this temperature
  if (use_ice) then
     call qsat(t, p, es, qs)
  else
     call qsat_water(t, p, es, qs)
  end if

  ! make sure a meaningful calculation is possible
  if (p <= 5._r8*es .or. qs <= 0._r8 .or. qs >= 0.5_r8 &
       .or. t < tmin .or. t > tmax) then
     status = 1
     ! Keep initial parameters when conditions aren't suitable
     tsp = t
     qsp = q
     enin = 1._r8
     enout(1) = 1._r8
     l_found = 1
     return
  end if

  ! Prepare to iterate
  status = 2

  ! Get initial enthalpy
  if (use_ice) then
     call calc_hltalt(t,hltalt(1))
  else
     call no_ip_hltalt(t,hltalt(1))
  end if
  enin = tq_enthalpy(t, q, hltalt(1))
  enin = (1-q)*cpair*t + q*cpliq*t+q*hltalt(1)

  ! make a guess at the wet bulb temp using a UKMO algorithm (from J. Petch)
  c1 = hltalt(1)*c3
  c2 = (t + 36._r8)**2
  r1b = c2/(c2 + c1*qs)
  qvd = r1b * (q - qs)
  cpmix(1) = ((1-q)*cpair + q*cpwv)
  tsp = t + ((hltalt(1)/cpmix(1))*qvd)
  !dtfg = ((hltalt(1)/cpmix(1))*qvd)

  dtfg = ((hltalt(1)/cpmix(1))*qvd)
  dtfg = max(-15._r8, dtfg)
  dtfg = min(15._r8, dtfg)

  ! Generate qsp, gam, and enout from tsp.
  if (use_ice) then
     call qsat(tsp, p, es, qsp, gam=gam(1), enthalpy=enout(1))
     call calc_hltalt(tsp,hltalt)
  else
     call qsat_water(tsp, p, es, qsp, gam=gam(1), enthalpy=enout(1))
     call no_ip_hltalt(tsp,hltalt(1))
  end if

  enout(1) = (1-q)*cpair*tsp + q*cpliq*tsp + q*hltalt(1) + ((qsp-q)/(1-qsp))*(cpliq*tsp - cpliq*t + hltalt(1))

  ! iterate on first guess
  do l = 1, iter

     l_found = l
     g(l) = enin - enout(l)
     cpmix(l) = ((1-q)*cpair + q*cpwv)
     dgdt(l) = -cpmix(l) * (1 + gam(l) + tsp*(cpwv-cpair)*gam(l)*(1/hltalt(l)))  !changed to account for cp,m changing with Ts as well
     dgdt(l) = (1-q)*cpair + q*cpwv + ((qsp-q)/(1-qsp))*cpwv + ((1-q)/((1-qsp)**2))*gam(l)*(cpmix(l)/hltalt(l))*(cpliq*tsp - cpliq*t + hltalt(l))

     dcpqdt(l) = tsp*(cpwv-cpair)*(1/hltalt(l))

     ! New tsp
     t1(l) = tsp - g(l)/dgdt(l)
     dt(l) = abs(t1(l) - tsp)/t1(l)
     tsp = t1(l)

     ! bail out if past end of temperature range
     if ( tsp < tmin ) then
        tsp = tmin
        ! Get latent heat and set qsp to a value that preserves enthalpy.
        if (use_ice) then
           call calc_hltalt(tsp,hltalt(l))
        else
           call no_ip_hltalt(tsp,hltalt(l))
        end if
        qsp = (enin - cpair*tsp)/(hltalt(l) + tsp*(cpwv-cpair))  !TODO: UPDATE PRIORITY: LOW
        enout(l) = tq_enthalpy(tsp, qsp, hltalt(l))
        status = 4
        exit
     end if

     ! Re-generate qsp, gam, and enout from new tsp.
     if (use_ice) then
        call qsat(tsp, p, es, q1(l), gam=gam(l), enthalpy=enout(l))
     else
        call qsat_water(tsp, p, es, q1(l), gam=gam(l), enthalpy=enout(l))
     end if
     dq(l) = abs(q1(l) - qsp)/max(q1(l),1.e-12_r8)
     qsp = q1(l)

     if (qsp==1._r8) then !fatal error, q=1
         status = 2
         ! enout = (1-q)*cpair*tsp + q*cpliq*tsp + q*hltalt + ((qsp-q)/(1-qsp))*(cpliq*tsp - cpliq*t + hltalt)
         exit
     endif

     enout(l) = (1-q)*cpair*tsp + q*cpliq*tsp + q*hltalt(l) + ((qsp-q)/(1-qsp))*(cpliq*tsp - cpliq*t + hltalt(l)) !output correct enout


     ! if converged at this point, exclude it from more iterations
     if (dt(l) < dttol .and. dq(l) < dqtol) then
        status = 0
        exit
     endif
  end do

  ! Test for enthalpy conservation
  if (abs((enin-enout(l))/(enin+enout(l))) > 1.e-4_r8) status = 8

!   if (status==8) then
!       write(iulog,*) " "
!       write(iulog,*) "findsp not converging."
!       write(iulog,*) ' t, q, p ', t, q, p
!       write(iulog,*) ' tsp, qsp ', tsp(i), qsp(i)
!       write(iulog,*) ' dt:', dt(i), 'dq:', dq(i), 
!       write(iulog,*) ' g:', g(i), 'dgdt:', dgdt(i), 'dtfg:', dtfg(i)
!       write(iulog,*) ' enin:', enin, "enout:", enout
!       write(iulog,*) ' t1:', t1, 'q1:', q1
!       write(iulog,*) ' cpmix:', cpmix, 'gam:', gam, 'tsp*(cpwv-cpair)*(1/hltalt)', tsp*(cpwv-cpair)*(1/hltalt)
!   end if

end subroutine findsp

end module wv_saturation
