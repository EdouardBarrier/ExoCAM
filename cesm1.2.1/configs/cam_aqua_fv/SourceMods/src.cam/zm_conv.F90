module zm_conv

!---------------------------------------------------------------------------------
! Purpose:
!
! Interface from Zhang-McFarlane convection scheme, includes evaporation of convective 
! precip from the ZM scheme
!
! Apr 2006: RBN: Code added to perform a dilute ascent for closure of the CM mass flux
!                based on an entraining plume a la Raymond and Blythe (1992)
!
! Author: Byron Boville, from code in tphysbc
!
!cas 2011 implementing b.sanderson's new ientropy closure 
!
!eb 2023 for non-dilute atmospheres
!
!---------------------------------------------------------------------------------
  use shr_kind_mod,    only: r8 => shr_kind_r8
  use spmd_utils,      only: masterproc
  use ppgrid,          only: pcols, pver, pverp
  use cloud_fraction,  only: cldfrc_fice
  use physconst,       only: cpair, epsilo, gravit, latice, latvap, tmelt, rair, &
                             cpliq, rh2o, zvir !different cpwv
  use shr_const_mod,   only: shr_const_cpice, shr_const_mwwv, shr_const_mwdair
  use abortutils,      only: endrun
  use cam_logfile,     only: iulog

  implicit none

  save
  private                         ! Make default type private to the module
!
! PUBLIC: interfaces
!
  public zmconv_readnl            ! read zmconv_nl namelist
  public zm_convi                 ! ZM schemea
  public zm_convr                 ! ZM schemea
  public zm_conv_evap             ! evaporation of precip from ZM schemea
  public convtran                 ! convective transport
  public momtran                  ! convective momentum transport
  public entropy                  ! returns entropy for a given p, q, T
  public ientropy                 ! given the entropy, pressure, qt, returns temperature

!
! Private data
!
   real(r8), parameter :: unset_r8 = huge(1.0_r8)
   real(r8) :: zmconv_c0_lnd = unset_r8    
   real(r8) :: zmconv_c0_ocn = unset_r8    
   real(r8) :: zmconv_ke     = unset_r8    

   real(r8) rl         ! wg latent heat of vaporization.
   real(r8) cpres      ! specific heat at constant pressure in j/kg-degk.
   real(r8), parameter :: capelmt = 70._r8  ! threshold value for cape for deep convection.
   real(r8) :: ke           ! Tunable evaporation efficiency set from namelist input zmconv_ke
   real(r8) :: c0_lnd       ! set from namelist input zmconv_c0_lnd
   real(r8) :: c0_ocn       ! set from namelist input zmconv_c0_ocn
   real(r8) tau   ! convective time scale
   real(r8),parameter :: a = 21.656_r8
   real(r8),parameter :: b = 5418._r8
   real(r8),parameter :: c1 = 6.112_r8
   real(r8),parameter :: c2 = 17.67_r8
   real(r8),parameter :: c3 = 243.5_r8
   real(r8) :: tfreez
   real(r8) :: eps1
   real(r8),parameter :: cpwv = 1860._r8

   logical :: no_deep_pbl ! default = .false.
                          ! no_deep_pbl = .true. eliminates deep convection entirely within PBL 
   

!moved from moistconvection.F90
   real(r8) :: rgrav       ! reciprocal of grav
   real(r8) :: rgas        ! gas constant for dry air
   real(r8) :: grav        ! = gravit
   real(r8) :: cp          ! = cpres = cpair

   real(r8) :: mu_red      !reduced mean molar mass difference, Leconte 2017
   real(r8) :: cpice       !heat capacity of ice 1h
   real(r8) :: rlice       !latent heat of melting
   
   integer  limcnv       ! top interface level limit for convection

   real(r8),parameter ::  tiedke_add = 0.5_r8   

contains

subroutine zmconv_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'zmconv_readnl'

   namelist /zmconv_nl/ zmconv_c0_lnd, zmconv_c0_ocn, zmconv_ke
   !-----------------------------------------------------------------------------

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'zmconv_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, zmconv_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)

      ! set local variables
      c0_lnd = zmconv_c0_lnd
      c0_ocn = zmconv_c0_ocn
      ke = zmconv_ke

   end if

#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast(c0_lnd,            1, mpir8,  0, mpicom)
   call mpibcast(c0_ocn,            1, mpir8,  0, mpicom)
   call mpibcast(ke,                1, mpir8,  0, mpicom)
#endif

end subroutine zmconv_readnl


subroutine zm_convi(limcnv_in, no_deep_pbl_in)

   use dycore,       only: dycore_is, get_resolution

   integer, intent(in)           :: limcnv_in       ! top interface level limit for convection
   logical, intent(in), optional :: no_deep_pbl_in  ! no_deep_pbl = .true. eliminates ZM convection entirely within PBL 

   ! local variables
   character(len=32)   :: hgrid           ! horizontal grid specifier

   ! Initialization of ZM constants
   limcnv = limcnv_in
   tfreez = tmelt
   eps1   = epsilo
   rl     = latvap
   cpres  = cpair
   rgrav  = 1.0_r8/gravit
   rgas   = rair
   grav   = gravit
   cp     = cpres

   mu_red = (shr_const_mwwv - shr_const_mwdair)/shr_const_mwwv
   cpice = shr_const_cpice
   rlice = latice

   if ( present(no_deep_pbl_in) )  then
      no_deep_pbl = no_deep_pbl_in
   else
      no_deep_pbl = .false.
   endif

   ! tau=4800. were used in canadian climate center. however, in echam3 t42, 
   ! convection is too weak, thus adjusted to 2400.

   hgrid = get_resolution()
   tau = 3600._r8

   if ( masterproc ) then
      write(iulog,*) 'tuning parameters zm_convi: tau',tau
      write(iulog,*) 'tuning parameters zm_convi: c0_lnd',c0_lnd, ', c0_ocn', c0_ocn 
      write(iulog,*) 'tuning parameters zm_convi: ke',ke
      write(iulog,*) 'tuning parameters zm_convi: no_deep_pbl',no_deep_pbl
   endif

   if (masterproc) write(iulog,*)'**** ZM: DILUTE Buoyancy Calculation ****'

end subroutine zm_convi



subroutine zm_convr(lchnk   ,ncol    , &
   t       ,qh      ,prec    ,jctop_1 ,jctop_2 , &
   jcbot_1 ,jcbot_2 ,maxg_1  ,maxg_2  ,ideep   , &
   pblh    ,zm      ,geos    ,zi      ,qtnd    , &
   heat    ,pap     ,paph    ,dpp     , &
   delt    ,mcon    ,cme     ,cape_1  ,cape_2  , &
   tpert   ,dlf     ,pflx    ,zdu     ,rprd    , &
   mu      ,md      ,du      ,eu      ,ed      , &
   dp      ,dsubcld_1,   dsubcld_2    ,jt_1    ,jt_2 , &
   maxi_1  ,maxi_2  ,lengath ,ql      ,rliq    ,landfrac, sfcT, heat_2, qtnd_2, lcl_reached)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Main driver for zhang-mcfarlane convection scheme 
! 
! Method: 
! performs deep convective adjustment based on mass-flux closure
! algorithm.
! 
! Author:guang jun zhang, m.lazare, n.mcfarlane. CAM Contact: P. Rasch
!
! This is contributed code not fully standardized by the CAM core group.
! All variables have been typed, where most are identified in comments
! The current procedure will be reimplemented in a subsequent version
! of the CAM where it will include a more straightforward formulation
! and will make use of the standard CAM nomenclature
! 
!-----------------------------------------------------------------------
   use constituents, only: pcnst
   use phys_control, only: cam_physpkg_is

!
! ************************ index of variables **********************
!
!  wg * alpha    array of vertical differencing used (=1. for upstream).
!  w  * cape     convective available potential energy.
!  wg * capeg    gathered convective available potential energy.
!  c  * capelmt  threshold value for cape for deep convection.
!  ic  * cpres    specific heat at constant pressure in j/kg-degk.
!  i  * dpp      
!  ic  * delt     length of model time-step in seconds.
!  wg * dp       layer thickness in mbs (between upper/lower interface).
!  wg * dqdt     mixing ratio tendency at gathered points.
!  wg * dsdt     dry static energy ("temp") tendency at gathered points.
!  wg * dudt     u-wind tendency at gathered points.
!  wg * dvdt     v-wind tendency at gathered points.
!  wg * dsubcld  layer thickness in mbs between lcl and maxi.
!  ic  * grav     acceleration due to gravity in m/sec2.
!  wg * du       detrainment in updraft. specified in mid-layer
!  wg * ed       entrainment in downdraft.
!  wg * eu       entrainment in updraft.
!  wg * hmn      moist static energy.
!  wg * hsat     saturated moist static energy.
!  w  * ideep    holds position of gathered points vs longitude index.
!  ic  * pver     number of model levels.
!  wg * j0       detrainment initiation level index.
!  wg * jd       downdraft   initiation level index.
!  ic  * jlatpr   gaussian latitude index for printing grids (if needed).
!  wg * jt       top  level index of deep cumulus convection.
!  w  * lcl      base level index of deep cumulus convection.
!  wg * lclg     gathered values of lcl.
!  w  * lel      index of highest theoretical convective plume.
!  wg * lelg     gathered values of lel.
!  w  * lon      index of onset level for deep convection.
!  w  * maxi     index of level with largest moist static energy.
!  wg * maxg     gathered values of maxi.
!  wg * mb       cloud base mass flux.
!  wg * mc       net upward (scaled by mb) cloud mass flux.
!  wg * md       downward cloud mass flux (positive up).
!  wg * mu       upward   cloud mass flux (positive up). specified
!                at interface
!  ic  * msg      number of missing moisture levels at the top of model.
!  w  * p        grid slice of ambient mid-layer pressure in mbs.
!  i  * pblt     row of pbl top indices.
!  w  * pcpdh    scaled surface pressure.
!  w  * pf       grid slice of ambient interface pressure in mbs.
!  wg * pg       grid slice of gathered values of p.
!  w  * q        grid slice of mixing ratio.
!  wg * qd       grid slice of mixing ratio in downdraft.
!  wg * qg       grid slice of gathered values of q.
!  i/o * qh       grid slice of specific humidity.
!  w  * qh0      grid slice of initial specific humidity.
!  wg * qhat     grid slice of upper interface mixing ratio.
!  wg * ql       grid slice of cloud liquid water.
!  wg * qs       grid slice of saturation mixing ratio.
!  w  * qstp     grid slice of parcel temp. saturation mixing ratio.
!  wg * qstpg    grid slice of gathered values of qstp.
!  wg * qu       grid slice of mixing ratio in updraft.
!  ic  * rgas     dry air gas constant.
!  wg * rl       latent heat of vaporization.
!  w  * s        grid slice of scaled dry static energy (t+gz/cp).
!  wg * sd       grid slice of dry static energy in downdraft.
!  wg * sg       grid slice of gathered values of s.
!  wg * shat     grid slice of upper interface dry static energy.
!  wg * su       grid slice of dry static energy in updraft.
!  i/o * t       
!  o  * jctop    row of top-of-deep-convection indices passed out.
!  O  * jcbot    row of base of cloud indices passed out.
!  wg * tg       grid slice of gathered values of t.
!  w  * tl       row of parcel temperature at lcl.
!  wg * tlg      grid slice of gathered values of tl.
!  w  * tp       grid slice of parcel temperatures.
!  wg * tpg      grid slice of gathered values of tp.
!  i/o * u        grid slice of u-wind (real).
!  wg * ug       grid slice of gathered values of u.
!  i/o * utg      grid slice of u-wind tendency (real).
!  i/o * v        grid slice of v-wind (real).
!  w  * va       work array re-used by called subroutines.
!  wg * vg       grid slice of gathered values of v.
!  i/o * vtg      grid slice of v-wind tendency (real).
!  i  * w        grid slice of diagnosed large-scale vertical velocity.
!  w  * z        grid slice of ambient mid-layer height in metres.
!  w  * zf       grid slice of ambient interface height in metres.
!  wg * zfg      grid slice of gathered values of zf.
!  wg * zg       grid slice of gathered values of z.
!
!-----------------------------------------------------------------------
!
! multi-level i/o fields:
!  i      => input arrays.
!  i/o    => input/output arrays.
!  w      => work arrays.
!  wg     => work arrays operating only on gathered points.
!  ic     => input data constants.
!  c      => data constants pertaining to subroutine itself.
!
! input arguments
!
   integer, intent(in) :: lchnk                   ! chunk identifier
   integer, intent(in) :: ncol                    ! number of atmospheric columns

   real(r8), intent(in) :: t(pcols,pver)          ! grid slice of temperature at mid-layer.
   real(r8), intent(in) :: qh(pcols,pver,pcnst)   ! grid slice of specific humidity.
   real(r8), intent(in) :: pap(pcols,pver)     
   real(r8), intent(in) :: paph(pcols,pver+1)
   real(r8), intent(in) :: dpp(pcols,pver)        ! local sigma half-level thickness (i.e. dshj).
   real(r8), intent(in) :: zm(pcols,pver)
   real(r8), intent(in) :: geos(pcols)
   real(r8), intent(in) :: zi(pcols,pver+1)
   real(r8), intent(in) :: pblh(pcols)
   real(r8), intent(in) :: tpert(pcols)
   real(r8), intent(in) :: landfrac(pcols) ! RBN Landfrac
   real(r8), intent(in) :: sfcT(pcols)
!
! output arguments
!
   real(r8), intent(out) :: qtnd(pcols,pver)           ! specific humidity tendency (kg/kg/s)
   real(r8), intent(out) :: heat(pcols,pver)           ! heating rate (dry static energy tendency, W/kg)
   real(r8), intent(out) :: qtnd_2(pcols,pver)           ! specific humidity tendency (kg/kg/s)
   real(r8), intent(out) :: heat_2(pcols,pver)           ! heating rate (dry static energy tendency, W/kg)
   real(r8), intent(out) :: mcon(pcols,pverp)
   real(r8), intent(out) :: dlf(pcols,pver)    ! scattrd version of the detraining cld h2o tend
   real(r8), intent(out) :: pflx(pcols,pverp)  ! scattered precip flux at each level
   real(r8), intent(out) :: cme(pcols,pver)
   real(r8), intent(out) :: cape_1(pcols)        ! w  convective available potential energy.
   real(r8), intent(out) :: cape_2(pcols)  !need to make this an output
   real(r8), intent(out) :: zdu(pcols,pver)
   real(r8), intent(out) :: rprd(pcols,pver)     ! rain production rate
! move these vars from local storage to output so that convective
! transports can be done in outside of conv_cam.
   real(r8), intent(out) :: mu(pcols,pver)
   real(r8) mu_1(pcols,pver)
   real(r8) mu_2(pcols,pver)
   real(r8), intent(out) :: eu(pcols,pver)
   real(r8) eu_1(pcols,pver)
   real(r8) eu_2(pcols,pver)
   real(r8), intent(out) :: du(pcols,pver)
   real(r8) du_1(pcols,pver)
   real(r8) du_2(pcols,pver)
   real(r8), intent(out) :: md(pcols,pver)
   real(r8) md_1(pcols,pver)
   real(r8) md_2(pcols,pver)
   real(r8), intent(out) :: ed(pcols,pver)
   real(r8) ed_1(pcols,pver)
   real(r8) ed_2(pcols,pver)
   real(r8), intent(out) :: dp(pcols,pver)       ! wg layer thickness in mbs (between upper/lower interface).
   real(r8), intent(out) :: dsubcld_1(pcols)       ! wg layer thickness in mbs between lcl and maxi.
   real(r8), intent(out) :: dsubcld_2(pcols)       ! wg layer thickness in mbs between lcl and maxi.
   real(r8), intent(out) :: jctop_1(pcols)  ! o row of top-of-deep-convection indices passed out.
   real(r8), intent(out) :: jcbot_1(pcols)  ! o row of base of cloud indices passed out.
   real(r8), intent(out) :: jctop_2(pcols)  ! o row of top-of-deep-convection indices passed out.
   real(r8), intent(out) :: jcbot_2(pcols)  ! o row of base of cloud indices passed out.
   real(r8), intent(out) :: prec(pcols)
   real(r8), intent(out) :: rliq(pcols) ! reserved liquid (not yet in cldliq) for energy integrals
   logical,  intent(out) :: lcl_reached(pcols)  ! whether the lcl is reached in a given column


   real(r8) zs(pcols)
   real(r8) dlg(pcols,pver)    ! gathrd version of the detraining cld h2o tend
   real(r8) dlg_1(pcols,pver)    ! gathrd version of the detraining cld h2o tend
   real(r8) dlg_2(pcols,pver)    ! gathrd version of the detraining cld h2o tend
   real(r8) pflxg(pcols,pverp) ! gather precip flux at each level
   real(r8) pflxg_1(pcols,pverp) ! gather precip flux at each level
   real(r8) pflxg_2(pcols,pverp) ! gather precip flux at each level
   real(r8) cug(pcols,pver)    ! gathered condensation rate
   real(r8) cug_1(pcols,pver)    ! gathered condensation rate
   real(r8) cug_2(pcols,pver)    ! gathered condensation rate
   real(r8) evpg(pcols,pver)   ! gathered evap rate of rain in downdraft
   real(r8) evpg_1(pcols,pver)   ! gathered evap rate of rain in downdraft
   real(r8) evpg_2(pcols,pver)   ! gathered evap rate of rain in downdraft
   real(r8) mumax_1(pcols)
   real(r8) mumax_2(pcols)


   integer jt(pcols)                          ! wg top  level index of deep cumulus convection.
   integer jt_1(pcols)
   integer jt_2(pcols)
   integer maxg_1(pcols)                        ! wg gathered values of maxi.
   integer maxg_2(pcols)                      ! wg gathered values of mx_2
   integer, intent(out) ::  maxi_1(pcols)       ! w  index of level with largest moist static energy
   integer, intent(out) ::  maxi_2(pcols)       ! w  index of level with largest local moist static energy
   integer ideep(pcols)                       ! w holds position of gathered points vs longitude index.
   integer lengath
!     diagnostic field used by chem/wetdep codes
   real(r8) ql(pcols,pver)                    ! wg grid slice of cloud liquid water.
!
   real(r8) pblt(pcols)           ! i row of pbl top indices.




!
!-----------------------------------------------------------------------
!
! general work fields (local variables):
!
   real(r8) q(pcols,pver)              ! w  grid slice of mixing ratio.
   real(r8) p(pcols,pver)              ! w  grid slice of ambient mid-layer pressure in mbs.
   real(r8) z(pcols,pver)              ! w  grid slice of ambient mid-layer height in metres.
   real(r8) s(pcols,pver)              ! w  grid slice of scaled dry static energy (t+gz/cp).
   real(r8) tp(pcols,pver)             ! w  grid slice of parcel temperatures.
   real(r8) zf(pcols,pver+1)           ! w  grid slice of ambient interface height in metres.
   real(r8) pf(pcols,pver+1)           ! w  grid slice of ambient interface pressure in mbs.
   real(r8) qstp(pcols,pver)           ! w  grid slice of parcel temp. saturation mixing ratio.

   real(r8) tl(pcols)                  ! w  row of parcel temperature at lcl.

   integer lcl(pcols)                  ! w  base level index of deep cumulus convection.
   integer lcl_2(pcols)                ! w  base level index of deep cumulus convection.

   integer lel(pcols)                  ! w  index of highest theoretical convective plume.
   integer lon(pcols)                  ! w  index of onset level for deep convection.
   integer index(pcols)
   real(r8) precip

!
! gathered work fields:
!
   real(r8) qg(pcols,pver)             ! wg grid slice of gathered values of q.
   real(r8) tg(pcols,pver)             ! w  grid slice of temperature at interface.
   real(r8) pg(pcols,pver)             ! wg grid slice of gathered values of p.
   real(r8) zg(pcols,pver)             ! wg grid slice of gathered values of z.
   real(r8) sg(pcols,pver)             ! wg grid slice of gathered values of s.
   real(r8) tpg(pcols,pver)            ! wg grid slice of gathered values of tp.
   real(r8) zfg(pcols,pver+1)          ! wg grid slice of gathered values of zf.
   real(r8) pfg(pcols,pver+1)          ! wg grid slice of gathered values of pf.
   real(r8) qstpg(pcols,pver)          ! wg grid slice of gathered values of qstp.
   real(r8) qstpg_2(pcols,pver)
   real(r8) ug(pcols,pver)             ! wg grid slice of gathered values of u.
   real(r8) vg(pcols,pver)             ! wg grid slice of gathered values of v.
   real(r8) cmeg(pcols,pver)           ! wg grid slice of gathered values of net condensation
   real(r8) cmeg_1(pcols,pver)
   real(r8) cmeg_2(pcols,pver)

   real(r8) rprdg(pcols,pver)           ! wg gathered rain production rate
   real(r8) rprdg_1(pcols,pver)         ! wg gathered rain production rate
   real(r8) rprdg_2(pcols,pver)         ! wg gathered rain production rate

   real(r8) capeg(pcols)               ! wg gathered convective available potential energy.
   real(r8) capeg_2(pcols)
   real(r8) tlg(pcols)                 ! wg grid slice of gathered values of tl.
   real(r8) landfracg(pcols)           ! wg grid slice of landfrac  

   integer lclg(pcols)                 ! wg gathered values of lcl. (buoyan_dilute estimate)
   integer lclg_2(pcols)               ! wg gathered values of lcl.(buoyan_dilute estimate)
   integer lelg(pcols)
   integer lelg_2(pcols)
   real(r8) tpg_2(pcols,pver)
   logical lcl_reachedg(pcols) !wg gathered values of lcl reached
!
! work fields arising from gathered calculations.
!
   real(r8) dqdt(pcols,pver)           ! wg mixing ratio tendency at gathered points.
   real(r8) dqdt_1(pcols,pver)         ! wg mixing ratio tendency at gathered points.
   real(r8) dqdt_2(pcols,pver)         ! wg mixing ratio tendency at gathered points.

   real(r8) dsdt(pcols,pver)           ! wg dry static energy ("temp") tendency at gathered points.
   real(r8) dsdt_1(pcols,pver)         ! wg dry static energy ("temp") tendency at gathered points.
   real(r8) dsdt_2(pcols,pver)         ! wg dry static energy ("temp") tendency at gathered points.

!      real(r8) alpha(pcols,pver)      ! array of vertical differencing used (=1. for upstream).
   real(r8) sd(pcols,pver)             ! wg grid slice of dry static energy in downdraft.
   real(r8) sd_1(pcols,pver)             ! wg grid slice of dry static energy in downdraft.
   real(r8) sd_2(pcols,pver)             ! wg grid slice of dry static energy in downdraft.
   real(r8) qd(pcols,pver)             ! wg grid slice of mixing ratio in downdraft.
   real(r8) qd_1(pcols,pver)             ! wg grid slice of mixing ratio in downdraft.
   real(r8) qd_2(pcols,pver)             ! wg grid slice of mixing ratio in downdraft.
   real(r8) mc(pcols,pver)             ! wg net upward (scaled by mb) cloud mass flux.
   real(r8) mc_1(pcols,pver)             ! wg net upward (scaled by mb) cloud mass flux.
   real(r8) mc_2(pcols,pver)             ! wg net upward (scaled by mb) cloud mass flux.
   real(r8) qu(pcols,pver)             ! wg grid slice of mixing ratio in updraft.
   real(r8) qu_1(pcols,pver)             ! wg grid slice of mixing ratio in updraft.
   real(r8) qu_2(pcols,pver)             ! wg grid slice of mixing ratio in updraft.
   real(r8) su(pcols,pver)             ! wg grid slice of dry static energy in updraft.
   real(r8) su_1(pcols,pver)             ! wg grid slice of dry static energy in updraft.
   real(r8) su_2(pcols,pver)             ! wg grid slice of dry static energy in updraft.
   real(r8) qlg(pcols,pver)
   real(r8) qlg_1(pcols,pver)
   real(r8) qlg_2(pcols,pver)
   real(r8) qs(pcols,pver)             ! wg grid slice of saturation mixing ratio.
   real(r8) shat(pcols,pver)           ! wg grid slice of upper interface dry static energy.
   real(r8) qhat(pcols,pver)           ! wg grid slice of upper interface mixing ratio.
   real(r8) hmn(pcols,pver)            ! wg moist static energy.
   real(r8) hsat(pcols,pver)           ! wg saturated moist static energy.
   real(r8) dudt(pcols,pver)           ! wg u-wind tendency at gathered points.
   real(r8) dvdt(pcols,pver)           ! wg v-wind tendency at gathered points.

   real(r8) mb_1(pcols)                  ! wg cloud base mass flux.
   real(r8) mb_2(pcols)                  ! wg cloud base mass flux.


   integer jlcl_1(pcols)              ! wg  base level index of deep cumulus convection.
   integer jlcl_2(pcols)              ! wg  base level index of deep cumulus convection.
   integer j0_1(pcols)                 ! wg detrainment initiation level index.
   integer jd_1(pcols)                 ! wg downdraft initiation level index.
   integer j0_2(pcols)                 ! wg detrainment initiation level index.
   integer jd_2(pcols)                 ! wg downdraft initiation level index.

   real(r8) delt                     ! length of model time-step in seconds.

   integer i
   integer ii
   integer k
   integer msg                      !  ic number of missing moisture levels at the top of model.
   real(r8) qdifr
   real(r8) sdifr
   real(r8) cpmix
   integer y

   integer lel_2(pcols)
   real(r8) tp_2(pcols,pver)
   real(r8) qstp_2(pcols,pver)

   real(r8) max_mc
   real(r8) max_dsdt
   real(r8) max_dqdt
   integer max_mc_i
   integer max_mc_k
   integer max_dsdt_i
   integer max_dsdt_k
   integer max_dqdt_i
   integer max_dqdt_k

!
!--------------------------Data statements------------------------------
!
! Set internal variable "msg" (convection limit) to "limcnv-1"
!
   msg = limcnv - 1
!
! zero tracker variables for things going wrong
   max_mc = 0._r8
   max_dsdt = 0._r8
   max_dqdt = 0._r8
   max_mc_i = 0
   max_mc_k = 0
   max_dsdt_i = 0
   max_dsdt_k = 0
   max_dqdt_i = 0
   max_dqdt_k = 0
!
! initialize necessary arrays.
! zero out variables not used in cam
!
   qtnd(:,:) = 0._r8
   qtnd_2(:,:) = 0._r8
   heat(:,:) = 0._r8
   heat_2(:,:) = 0._r8
   mcon(:,:) = 0._r8
   rliq(:ncol)   = 0._r8
   lcl_reached(:pcols) = .false.
   lcl_reachedg(:pcols) = .false.
!
! initialize convective tendencies
!
   prec(:ncol) = 0._r8
   do k = 1,pver
      do i = 1,ncol
         dqdt(i,k)  = 0._r8
         dqdt_1(i,k)  = 0._r8
         dqdt_2(i,k)  = 0._r8
         dsdt(i,k)  = 0._r8
         dsdt_1(i,k)  = 0._r8
         dsdt_2(i,k)  = 0._r8
         dudt(i,k)  = 0._r8
         dvdt(i,k)  = 0._r8
         pflx(i,k)  = 0._r8
         pflxg(i,k) = 0._r8
         pflxg_1(i,k) = 0._r8
         pflxg_2(i,k) = 0._r8
         cme(i,k)   = 0._r8
         rprd(i,k)  = 0._r8
         zdu(i,k)   = 0._r8
         ql(i,k)    = 0._r8
         qlg(i,k)   = 0._r8
         qlg_1(i,k)   = 0._r8
         qlg_2(i,k)   = 0._r8
         dlf(i,k)   = 0._r8
         dlg(i,k)   = 0._r8
         dlg_1(i,k)   = 0._r8
         dlg_2(i,k)   = 0._r8
      end do
   end do
   do i = 1,ncol
      pflx(i,pverp) = 0
      pflxg(i,pverp) = 0
      pflxg_1(i,pverp) = 0
      pflxg_2(i,pverp) = 0

   end do
!
   do i = 1,ncol
      pblt(i) = pver
      jctop_1(i) = pver
      jcbot_1(i) = 1
      jctop_2(i) = pver
      jcbot_2(i) = 1
   end do
!
! calculate local pressure (mbs) and height (m) for both interface
! and mid-layer locations.
!
   do i = 1,ncol
      zs(i) = geos(i)*rgrav
      pf(i,pver+1) = paph(i,pver+1)*0.01_r8
      zf(i,pver+1) = zi(i,pver+1) + zs(i)
   end do
   do k = 1,pver
      do i = 1,ncol
         p(i,k) = pap(i,k)*0.01_r8
         pf(i,k) = paph(i,k)*0.01_r8
         z(i,k) = zm(i,k) + zs(i)
         zf(i,k) = zi(i,k) + zs(i)
      end do
   end do
!
   do k = pver - 1,msg + 1,-1
      do i = 1,ncol
         if (abs(z(i,k)-zs(i)-pblh(i)) < (zf(i,k)-zf(i,k+1))*0.5_r8) pblt(i) = k
      end do
   end do
!
! store incoming specific humidity field for subsequent calculation
! of precipitation (through change in storage).
! define dry static energy (normalized by cp).
!
   do k = 1,pver
      do i = 1,ncol
         q(i,k) = qh(i,k,1)
         cpmix = (1-q(i,k))*cpres + q(i,k)*cpwv
         s(i,k) = (cpmix/cpres)*t(i,k) + (grav/cpres)*z(i,k)
         tp(i,k)=0.0_r8
         tp_2(i,k) = 0.0_r8
         shat(i,k) = s(i,k)
         qhat(i,k) = q(i,k)
      end do
   end do

   do i = 1,ncol
      capeg(i) = 0._r8
      capeg_2(i) = 0._r8
      lclg(i) = 1
      lclg_2(i) = 1
      lelg(i) = pver
      maxg_1(i) = 1
      lelg_2(i) = pver
      maxg_2(i) = 1
      tlg(i) = 400._r8
      dsubcld_1(i) = 0._r8
      dsubcld_2(i) = 0._r8
      index(i) = 0
      ideep(i) = 0
      !write(iulog,*) "pcols:", pcols, "ncol:", ncol, "i:", i, "ideep(i):", ideep(i)
   end do

      !  Evaluate Tparcel, qsat(Tparcel), buoyancy and CAPE, 
      !     lcl, lel, parcel launch level at index maxi()=hmax

   call buoyan_dilute(lchnk   ,ncol    , &
               q       ,t       ,p       ,z       ,pf       , &
               tp      ,qstp    ,tl      ,rl      ,cape_1     , &
               pblt    ,lcl     ,lel     ,lon     ,maxi_1     , &
               maxi_2    ,rgas    ,grav    ,cpres   ,msg     , &
               tpert   ,tp_2    ,qstp_2  ,cape_2  ,lel_2,  lcl_2, sfcT)

!
! determine whether grid points will undergo some deep convection
! (ideep=1) or not (ideep=0), based on values of cape,lcl,lel
! (require cape.gt. 0 and lel<lcl as minimum conditions).
! check done using the maximum of the two convective region
!
   lengath = 0
   do i=1,ncol
      if (max(cape_1(i),cape_2(i)) > capelmt) then
         lengath = lengath + 1
         index(lengath) = i
      end if
   end do

!write(iulog,*) "lengath:", lengath
   if (lengath.eq.0) return
   do ii=1,lengath
      i=index(ii)
      ideep(ii)=i
      !if (lchnk==122) write(iulog,*) "i:", i, "ii:", ii, "ideep(i):", ideep(ii)
   end do

! write(iulog,*) "*** ZM_CONV: Finding convective regions. lchnk:", lchnk, "lengath: ", lengath, "ideep: ", ideep

!
! obtain gathered arrays necessary for ensuing calculations.
!
   do k = 1,pver
      do i = 1,lengath
         dp(i,k) = 0.01_r8*dpp(ideep(i),k)
         qg(i,k) = q(ideep(i),k)
         tg(i,k) = t(ideep(i),k)
         pg(i,k) = p(ideep(i),k)
         zg(i,k) = z(ideep(i),k)
         sg(i,k) = s(ideep(i),k)
         tpg(i,k) = tp(ideep(i),k)
         tpg_2(i,k) = tp_2(ideep(i),k)
         zfg(i,k) = zf(ideep(i),k)
         pfg(i,k) = pf(ideep(i),k)
         qstpg(i,k) = qstp(ideep(i),k)
         qstpg_2(i,k) = qstp_2(ideep(i),k)
         ug(i,k) = 0._r8
         vg(i,k) = 0._r8
      end do
   end do
!
   do i = 1,lengath
      zfg(i,pver+1) = zf(ideep(i),pver+1)
   end do
   do i = 1,lengath
      capeg(i) = cape_1(ideep(i))
      lclg(i) = lcl(ideep(i))
      lelg(i) = lel(ideep(i))
      maxg_1(i) = maxi_1(ideep(i))
      tlg(i) = tl(ideep(i))
      landfracg(i) = landfrac(ideep(i))
      capeg_2(i) = cape_2(ideep(i))
      lclg_2(i) = lcl_2(ideep(i))
      lelg_2(i) = lel_2(ideep(i))
      maxg_2(i) = maxi_2(ideep(i))
   end do

   ! capeg_2(:) = 0._r8
   ! lelg_2(:) = 2
   ! maxg_2(:) = 2
!
! calculate sub-cloud layer pressure "thickness" for use in
! closure and tendency routines.
!
   do k = msg + 1,pver
      do i = 1,lengath
         if (k >= maxg_1(i) .and. k<=(maxg_1(i)+2)) then
            dsubcld_1(i) = dsubcld_1(i) + dp(i,k)
         end if
         if (k >= maxg_2(i) .and. k<=(maxg_2(i)+2)) then
            dsubcld_2(i) = dsubcld_2(i) + dp(i,k)
         end if
      end do
   end do

! write(iulog,*) "lchnk:", lchnk,"maxg_1:", maxg_1, "dsubcld_1:", dsubcld_1
!
! define array of factors (alpha) which defines interfacial
! values, as well as interfacial values for (q,s) used in
! subsequent routines.
!
   do k = msg + 2,pver
      do i = 1,lengath
!            alpha(i,k) = 0.5
         sdifr = 0._r8
         qdifr = 0._r8
         if (sg(i,k) > 0._r8 .or. sg(i,k-1) > 0._r8) &
            sdifr = abs((sg(i,k)-sg(i,k-1))/max(sg(i,k-1),sg(i,k)))
         if (qg(i,k) > 0._r8 .or. qg(i,k-1) > 0._r8) &
            qdifr = abs((qg(i,k)-qg(i,k-1))/max(qg(i,k-1),qg(i,k)))
         if (sdifr > 1.E-6_r8) then
            shat(i,k) = log(sg(i,k-1)/sg(i,k))*sg(i,k-1)*sg(i,k)/(sg(i,k-1)-sg(i,k))
         else
            shat(i,k) = 0.5_r8* (sg(i,k)+sg(i,k-1))
         end if
         if (qdifr > 1.E-6_r8) then
            qhat(i,k) = log(qg(i,k-1)/qg(i,k))*qg(i,k-1)*qg(i,k)/(qg(i,k-1)-qg(i,k))
         else
            qhat(i,k) = 0.5_r8* (qg(i,k)+qg(i,k-1))
         end if
      end do
   end do
!
!
! obtain cloud properties for first convective region
!
   call cldprp(lchnk   , &
               qg      ,tg      ,ug      ,vg      ,pg      , &
               zg      ,sg      ,mu_1    ,eu_1    ,du_1    , &
               md_1    ,ed_1    ,sd_1    ,qd_1    ,mc_1    , &
               qu_1    ,su_1    ,zfg     ,qs      ,hmn     , &
               hsat    ,shat    ,qhat    ,qlg_1   ,pf      , &
               cmeg_1  ,maxg_1  ,lelg    ,jt_1    ,jlcl_1  , &
               maxg_1  ,j0_1    ,jd_1    ,rl      ,lengath , &
               rgas    ,grav    ,cpres   ,msg     ,capeg   , &
               pflxg_1 ,evpg_1  ,cug_1   ,rprdg_1 ,limcnv  ,landfracg, lcl_reachedg, lclg)
!
! convert detrainment from units of "1/m" to "1/mb".
!
   do k = msg + 1,pver
      do i = 1,lengath
         du_1   (i,k) = du_1   (i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
         eu_1   (i,k) = eu_1   (i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
         ed_1   (i,k) = ed_1   (i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
         cug_1  (i,k) = cug_1  (i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
         cmeg_1 (i,k) = cmeg_1 (i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
         rprdg_1(i,k) = rprdg_1(i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
         evpg_1 (i,k) = evpg_1 (i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
      end do
   end do
!
! obtain cloud properties for second convective region
!
   call cldprp(lchnk   , &
               qg      ,tg      ,ug      ,vg      ,pg      , &
               zg      ,sg      ,mu_2    ,eu_2    ,du_2    , &
               md_2    ,ed_2    ,sd_2    ,qd_2    ,mc_2    , &
               qu_2    ,su_2    ,zfg     ,qs      ,hmn     , &
               hsat    ,shat    ,qhat    ,qlg_2   ,pf      , &
               cmeg_2  ,maxg_2  ,lelg_2  ,jt_2    ,jlcl_2  , &
               maxg_2  ,j0_2    ,jd_2    ,rl      ,lengath , &
               rgas    ,grav    ,cpres   ,msg     ,capeg_2 , &
               pflxg_2 ,evpg_2  ,cug_2   ,rprdg_2 ,limcnv  ,landfracg, lcl_reachedg, lclg_2)
!
! convert detrainment from units of "1/m" to "1/mb".
!
   do k = msg + 1,pver
      do i = 1,lengath
         du_2   (i,k) = du_2   (i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
         eu_2   (i,k) = eu_2   (i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
         ed_2   (i,k) = ed_2   (i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
         cug_2  (i,k) = cug_2  (i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
         cmeg_2 (i,k) = cmeg_2 (i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
         rprdg_2(i,k) = rprdg_2(i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
         evpg_2 (i,k) = evpg_2 (i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
      end do
   end do

   ! mu_2(:,:) = 0._r8 !debug
   ! eu_2(:,:) = 0._r8 
   ! du_2(:,:) = 0._r8 
   ! md_2(:,:) = 0._r8 
   ! ed_2(:,:) = 0._r8 
   ! sd_2(:,:) = 0._r8 
   ! qd_2(:,:) = 0._r8 
   ! mc_2(:,:) = 0._r8 
   ! qu_2(:,:) = 0._r8 
   ! su_2(:,:) = 0._r8
   ! qlg_2(:,:) = 0._r8 
   ! cmeg_2(:,:) = 0._r8 
   ! pflxg_2(:,:) = 0._r8 
   ! evpg_2(:,:) = 0._r8 
   ! cug_2(:,:) = 0._r8 
   ! rprdg_2(:,:) = 0._r8
   ! maxg_2(:) = 2
   ! lelg_2(:) = 2
   ! jt_2(:) = 2
   ! jlcl_2(:) = 2
   ! j0_2(:) = 2
   ! jd_2(:) = 2



!call the closure on the two convective regions to find the base mass fluxes in each

   call closure(lchnk   ,&
                qg      ,tg      ,pg      ,zg      ,sg      , &
                tpg     ,qs      ,qu_1    ,su_1    ,mc_1    , &
                du_1    ,mu_1    ,md_1    ,qd_1    ,sd_1    , &
                qhat    ,shat    ,dp      ,qstpg   ,pfg     ,zfg     , &
                qlg_1   ,dsubcld_1 ,mb_1    ,capeg   ,tlg     , &
                lclg    ,lelg    ,jt_1    ,maxg_1  ,1       , &
                lengath ,rgas    ,grav    ,cpres   ,rl      , &
                msg     ,capelmt    )

   call closure(lchnk   ,&
                qg      ,tg      ,pg      ,zg      ,sg      , &
                tpg_2   ,qs      ,qu_2    ,su_2    ,mc_2    , &
                du_2    ,mu_2    ,md_2    ,qd_2    ,sd_2    , &
                qhat    ,shat    ,dp      ,qstpg_2 ,pfg     ,zfg     , &
                qlg_2   ,dsubcld_2 ,mb_2      ,capeg_2 ,tlg     , &
                lclg_2  ,lelg_2  ,jt_2    ,maxg_2  ,1       , &
                lengath ,rgas    ,grav    ,cpres   ,rl      , &
                msg     ,capelmt    )
!
! limit cloud base mass flux to theoretical upper bound.
!
   ! mb_2(:) = 0._r8 !debug

   do i=1,lengath
      mumax_1(i) = 0
      mumax_2(i) = 0
   end do
   do k=msg + 2,pver
      do i=1,lengath
         mumax_1(i) = max(mumax_1(i), mu_1(i,k)/dp(i,k))
         mumax_2(i) = max(mumax_2(i), mu_2(i,k)/dp(i,k))
      end do
   end do

   do i=1,lengath
      if (mumax_1(i) > 0._r8) then
         mb_1(i) = min(mb_1(i),0.5_r8/(delt*mumax_1(i)))
         ! if (0.5_r8/(delt*mumax_1(i)) < mb_1(i)) then
         !    write(iulog,*) "i:", i, "mb_1:", mb_1(i), "cap:",  0.5_r8/(delt*mumax_1(i))
         ! end if
      else
         mb_1(i) = 0._r8
      endif

      if (mumax_2(i) > 0._r8) then
         mb_2(i) = min(mb_2(i),0.5_r8/(delt*mumax_2(i)))
         ! if (0.5_r8/(delt*mumax_2(i)) < mb_2(i)) then
         !    write(iulog,*) "i:", i, "mb_2:", mb_2(i), "cap:",  0.5_r8/(delt*mumax_2(i))
         ! end if
      else
         mb_2(i) = 0._r8
      endif
   end do
   ! If no_deep_pbl = .true., don't allow convection entirely 
   ! within PBL (suggestion of Bjorn Stevens, 8-2000)

   if (no_deep_pbl) then
      do i=1,lengath
         if (zm(ideep(i),jt_1(i)) < pblh(ideep(i))) mb_1(i) = 0
         if (zm(ideep(i),jt_2(i)) < pblh(ideep(i))) mb_2(i) = 0
      end do
   end if

   !if (lchnk==248) write(iulog,*) "mb:", mb(9)

   do k=msg+1,pver
      do i=1,lengath
         mu_1   (i,k)  = mu_1   (i,k)*mb_1(i)
         md_1   (i,k)  = md_1   (i,k)*mb_1(i)
         mc_1   (i,k)  = mc_1   (i,k)*mb_1(i)
         du_1   (i,k)  = du_1   (i,k)*mb_1(i)
         eu_1   (i,k)  = eu_1   (i,k)*mb_1(i)
         ed_1   (i,k)  = ed_1   (i,k)*mb_1(i)
         cmeg_1 (i,k)  = cmeg_1 (i,k)*mb_1(i)
         rprdg_1(i,k)  = rprdg_1(i,k)*mb_1(i)
         cug_1  (i,k)  = cug_1  (i,k)*mb_1(i)
         evpg_1 (i,k)  = evpg_1 (i,k)*mb_1(i)
         
         mu_2   (i,k)  = mu_2   (i,k)*mb_2(i)
         md_2   (i,k)  = md_2   (i,k)*mb_2(i)
         mc_2   (i,k)  = mc_2   (i,k)*mb_2(i)
         du_2   (i,k)  = du_2   (i,k)*mb_2(i)
         eu_2   (i,k)  = eu_2   (i,k)*mb_2(i)
         ed_2   (i,k)  = ed_2   (i,k)*mb_2(i)
         cmeg_2 (i,k)  = cmeg_2 (i,k)*mb_2(i)
         rprdg_2(i,k)  = rprdg_2(i,k)*mb_2(i)
         cug_2  (i,k)  = cug_2  (i,k)*mb_2(i)
         evpg_2 (i,k)  = evpg_2 (i,k)*mb_2(i)

         pflxg_1(i,k+1)= pflxg_1(i,k+1)*mb_1(i)*100._r8/grav
         pflxg_2(i,k+1)= pflxg_2(i,k+1)*mb_2(i)*100._r8/grav
      end do
   end do
!
! compute temperature and moisture changes due to convection.
!
   call q1q2_pjr(lchnk   , &
                 dqdt_1  ,dsdt_1  ,qg      ,qs      ,qu_1    , &
                 su_1    ,du_1    ,qhat    ,shat    ,dp      , &
                 mu_1    ,md_1    ,sd_1    ,qd_1    ,qlg_1   , &
                 dsubcld_1 ,jt_1    ,maxg_1  ,1       ,lengath , &
                 cpres   ,rl      ,msg     ,          &
                 dlg_1   ,evpg_1  ,cug_1   )

   ! do k = msg + 1,pver
   !    do i = 1,lengath
   !       if (q(ideep(i),k) + 2._r8*delt*dqdt_1(i,k) > 0.7_r8) then
   !          write(iulog,*) "After q1q2_pjr_1. lchnk:", lchnk, "i:", i, "k:", k
   !          write(iulog,*) "New q:", (q(ideep(i),k) + 2._r8*delt*dqdt_1(i,k)),"dqdt:", dqdt_1(i,k)
   !          write(iulog,*) "qhat:", qhat(i,:)
   !          write(iulog,*) "evp:", evpg_1(i,:)
   !          write(iulog,*) "cug:", cug_1(i,:)
   !          write(iulog,*) "mu:", mu_1(i,:)
   !          write(iulog,*) "qu:", qu_1(i,:)
   !          write(iulog,*) "md:", md_1(i,:)
   !          write(iulog,*) "qd:", qd_1(i,:)
   !       end if
   !    end do
   ! end do


   call q1q2_pjr(lchnk   , &
                 dqdt_2  ,dsdt_2  ,qg      ,qs      ,qu_2    , &
                 su_2    ,du_2    ,qhat    ,shat    ,dp      , &
                 mu_2    ,md_2    ,sd_2    ,qd_2    ,qlg_2   , &
                 dsubcld_2 ,jt_2    ,maxg_2  ,1       ,lengath , &
                 cpres   ,rl      ,msg     ,          &
                 dlg_2   ,evpg_2  ,cug_2   )
!
! merge upper and lower convective band variables
!
   ! dqdt_2(:,:) = 0._r8 !debug
   ! dsdt_2(:,:) = 0._r8
   ! dlg_2(:,:) = 0._r8

   ! do k = msg + 1,pver
   !    do i = 1,lengath
   !       if (q(ideep(i),k) + 2._r8*delt*dqdt_2(i,k) > 0.7_r8) then
   !          write(iulog,*) "After q1q2_pjr_2. lchnk:", lchnk, "i:", i, "k:", k
   !          write(iulog,*) "New q:", (q(ideep(i),k) + 2._r8*delt*dqdt_2(i,k)),"dqdt:", dqdt_2(i,k)
   !          write(iulog,*) "qhat:", qhat(i,:)
   !          write(iulog,*) "evp:", evpg_2(i,:)
   !          write(iulog,*) "cug:", cug_2(i,:)
   !          write(iulog,*) "mu:", mu_2(i,:)
   !          write(iulog,*) "qu:", qu_2(i,:)
   !          write(iulog,*) "md:", md_2(i,:)
   !          write(iulog,*) "qd:", qd_2(i,:)
   !       end if
   !    end do
   ! end do

   do k = msg + 1,pver
         do i = 1,lengath
            dqdt(i,k) =  dqdt_1(i,k) + dqdt_2(i,k) !
            dsdt(i,k) = dsdt_1(i,k) + dsdt_2(i,k) ! 
            dlg(i,k) = dlg_1(i,k) + dlg_2(i,k) !
            pflxg(i,k) = pflxg_1(i,k) + pflxg_2(i,k) ! 
            cug(i,k) = cug_1(i,k) +cug_2(i,k) ! 
            cmeg(i,k) = cmeg_1(i,k) + cmeg_2(i,k) ! 
            rprdg(i,k) = rprdg_1(i,k) + rprdg_2(i,k) ! 
            evpg(i,k) = evpg_1(i,k) + evpg_2(i,k) ! 
            mu(i,k) = mu_1(i,k) + mu_2(i,k) ! 
            md(i,k) = md_1(i,k) + md_2(i,k) ! 
            mc(i,k) = mc_1(i,k) + mc_2(i,k) ! 
            eu(i,k) = eu_1(i,k) + eu_2(i,k) ! 
            du(i,k) = du_1(i,k) + du_2(i,k) ! 
            ed(i,k) = ed_1(i,k) + ed_2(i,k) ! 
            qu(i,k) = qu_1(i,k) + qu_2(i,k) ! 
            qd(i,k) = qd_1(i,k) + qd_2(i,k) ! 
            su(i,k) = su_1(i,k) + su_2(i,k) ! 
            sd(i,k) = sd_1(i,k) + sd_2(i,k) ! 
            qlg(i,k) = qlg_1(i,k) + qlg_2(i,k) ! 
         end do
      end do
!
! gather back temperature and mixing ratio.
!
   ! if (lchnk==122) write(iulog,*) "dqdt:", dqdt(1,:), "dsdt:", dsdt(1,:)
   ! if (maxval(dsdt(i,:))>0.01_r8) then


   do ii = msg + 1,pver
      do i = 1,lengath
         if (dsdt_1(i,ii) /=dsdt_1(i,ii)) then
            do k = maxg_1(i), jt_1(i), -1
               write(iulog,*) " "
               write(iulog,*) "*** ZM_CONV: q1q2_pjr. k:", k, "dsdt:", dsdt(i,k), "old s:", sg(i,k), "new s:", sg(i,k)+1800._r8*dsdt(i,k)
               write(iulog,*) "*** ZM_CONV: q1q2_pjr. dsdt:", dsdt(i,k), "dqdt:", dqdt(i,k), "dl:", dlg(i,k)
               write(iulog,*) "*** Going up column. mu:", mu(i,k), "md:", md(i,k), "dp:", dp(i,k)
               write(iulog,*) "*** Going up column, env properties. shat:", shat(i,k), "qhat:", qhat(i,k), "su:", su(i,k), "sd", sd(i,k)
               write(iulog,*) "*** Going up column. updraft properties. su", su(i,k), "qu:", qu(i,k)
               write(iulog,*) "*** Going up column. downdraft properties. sd", sd(i,k), "qd:", qd(i,k)
               write(iulog,*) "*** Going up column. emc:", evpg(i,k)-cug(i,k), "latent heat:", -rl/cpres*(evpg(i,k)-cug(i,k))
               write(iulog,*) "updraft s leaving level:", mu(i,k)*(su(i,k)-shat(i,k)), "downdraft s joining level:", md(i,k)*(sd(i,k)-shat(i,k))
               write(iulog,*) "updraft q leaving level:", mu(i,k)*(qu(i,k)-qhat(i,k)), "downdraft q joining level:", md(i,k)*(qd(i,k)-qhat(i,k))
            end do
         end if
         if (dsdt_2(i,ii) /=dsdt_2(i,ii)) then
            do k = maxg_2(i), jt_2(i), -1
               write(iulog,*) " "
               write(iulog,*) "*** ZM_CONV: q1q2_pjr. k:", k, "dsdt:", dsdt(i,k), "old s:", sg(i,k), "new s:", sg(i,k)+1800._r8*dsdt(i,k)
               write(iulog,*) "*** ZM_CONV: q1q2_pjr. dsdt:", dsdt(i,k), "dqdt:", dqdt(i,k), "dl:", dlg(i,k)
               write(iulog,*) "*** Going up column. mu:", mu(i,k), "md:", md(i,k), "dp:", dp(i,k)
               write(iulog,*) "*** Going up column, env properties. shat:", shat(i,k), "qhat:", qhat(i,k), "su:", su(i,k), "sd", sd(i,k)
               write(iulog,*) "*** Going up column. updraft properties. su", su(i,k), "qu:", qu(i,k)
               write(iulog,*) "*** Going up column. downdraft properties. sd", sd(i,k), "qd:", qd(i,k)
               write(iulog,*) "*** Going up column. emc:", evpg(i,k)-cug(i,k), "latent heat:", -rl/cpres*(evpg(i,k)-cug(i,k))
               write(iulog,*) "updraft s leaving level:", mu(i,k)*(su(i,k)-shat(i,k)), "downdraft s joining level:", md(i,k)*(sd(i,k)-shat(i,k))
               write(iulog,*) "updraft q leaving level:", mu(i,k)*(qu(i,k)-qhat(i,k)), "downdraft q joining level:", md(i,k)*(qd(i,k)-qhat(i,k))
            end do
         end if
      end do
   end do

   ! if (lchnk==60 .and. i==8) then
   ! do i = 1, lengath
   !    write(iulog,*) " "
   !    write(iulog,*) " "
   !    write(iulog,*) "*** ZM_CONV: q1q2_pjr. At lchnk:", lchnk, "i: ", i
   !    ! write(iulog,*) "*** ZM_CONV: q1q2_pjr. mx:", mx_1(i), "jt:", jt_1(i)
   !    do k = maxg_1(i), jt_1(i), -1
   !       write(iulog,*) " "
   !       write(iulog,*) "*** ZM_CONV: q1q2_pjr. k:", k, "dsdt:", dsdt(i,k), "old s:", sg(i,k), "new s:", sg(i,k)+1800._r8*dsdt(i,k)
   !       ! write(iulog,*) "*** ZM_CONV: q1q2_pjr. dsdt:", dsdt(10,k), "dqdt:", dqdt(10,k), "dl:", dlg(10,k)
   !       ! write(iulog,*) "*** Going up column. mu:", mu(10,k), "md:", md(10,k), "dp:", dp(10,k)
   !       ! write(iulog,*) "*** Going up column, env properties. shat:", shat(10,k), "qhat:", qhat(10,k), "su:", su(10,k), "sd", sd(10,k)
   !       ! write(iulog,*) "*** Going up column. updraft properties. su", su(10,k), "qu:", qu(10,k)
   !       ! write(iulog,*) "*** Going up column. downdraft properties. sd", sd(10,k), "qd:", qd(10,k)
   !       ! write(iulog,*) "*** Going up column. emc:", evpg(10,k)-cug(10,k), "latent heat:", -rl/cpres*(evpg(10,k)-cug(10,k))
   !       ! write(iulog,*) "updraft s leaving level:", mu(10,k)*(su(10,0)-shat(10,k)), "downdraft s joining level:", md(10,k)*(sd(10,k)-shat(10,k))
   !       ! write(iulog,*) "updraft q leaving level:", mu(10,k)*(qu(10,k)-qhat(10,k)), "downdraft q joining level:", md(10,k)*(qd(10,k)-qhat(10,k))
   !    end do
   !    ! write(iulog,*) " "
   ! ! end if
   ! end do
   ! do i = 1, lengath
   !    write(iulog,*) " "
   !    write(iulog,*) " "
   !    write(iulog,*) "*** ZM_CONV: q1q2_pjr. At lchnk:", lchnk, "i: ", i
   !    ! write(iulog,*) "*** ZM_CONV: q1q2_pjr. mx:", mx_2(i), "jt:", jt_2(i)
   !    do k = maxg_2(i), jt_2(i), -1
   !       write(iulog,*) " "
   !       write(iulog,*) "*** ZM_CONV: q1q2_pjr. k:", k, "dsdt:", dsdt(i,k), "old s:", sg(i,k), "new s:", sg(i,k)+1800._r8*dsdt(i,k)
   !    end do
   !    ! write(iulog,*) " "
   ! ! end if
   ! end do


   do k = msg + 1,pver
!DIR$ CONCURRENT
      do i = 1,lengath
!
! q is updated to compute net precip.
!
         q(ideep(i),k) = qh(ideep(i),k,1) + 2._r8*delt*dqdt(i,k)
         qtnd(ideep(i),k) = dqdt (i,k)
         qtnd_2(ideep(i),k) = dqdt (i,k)
         cme (ideep(i),k) = cmeg (i,k)
         rprd(ideep(i),k) = rprdg(i,k)
         zdu (ideep(i),k) = du   (i,k)
         mcon(ideep(i),k) = mc   (i,k)
         heat(ideep(i),k) = dsdt (i,k)*cpres
         heat_2(ideep(i),k) = dsdt (i,k)*cpres
         dlf (ideep(i),k) = dlg  (i,k)
         pflx(ideep(i),k) = pflxg(i,k)
         ql  (ideep(i),k) = qlg  (i,k)

         ! if (lchnk==56) then
         !    write(iulog,*) "lchnk:", lchnk, "i:", i, "k:", k, "heat:", dsdt (i,k)*cpres
         !    ! write(iulog,*) "mcon:", mc(i,k), "heat:", dsdt (i,k)*cpres, "qtnd:", dqdt(i,k)
         !    ! write(iulog,*) " "
         ! end if
         ! if (mc(i,k)>max_mc) then
         !    max_mc = mc(i,k)
         !    max_mc_i = i
         !    max_mc_k = k
         ! end if
         ! if (2._r8*delt*dqdt(i,k) > 0.1_r8) then
         !    ! write(iulog,*) "lchnk:", lchnk, "i:", i, "k:", k, "dq:",2._r8*delt*dqdt(i,k)
         !    write(iulog,*) "*** ZM_CONV: q1q2_pjr. At lchnk:", lchnk, "i: ", i, "k:", k
         !    write(iulog,*) "*** ZM_CONV: q1q2_pjr. mx_1:", maxg_1(i), "jt_1:", jt_1(i), "mx_2:", maxg_2(i), "jt_2:", jt_2(i)
         !    do y = maxg_1(i), jt_1(i), -1
         !       write(iulog,*) " "
         !       ! write(iulog,*) "*** ZM_CONV: q1q2_pjr. k:", k, "dsdt:", dsdt(i,k), "old s:", sg(i,k), "new s:", sg(i,k)+1800._r8*dsdt(i,k)
         !       write(iulog,*) "*** ZM_CONV: q1q2_pjr. k:", y, "dsdt:", dsdt(i,y), "dqdt:", dqdt(i,y), "dl:", dlg(i,y)
         !       write(iulog,*) "*** Going up column. mu:", mu(i,y), "md:", md(i,y), "dp:", dp(i,y)
         !       write(iulog,*) "*** Going up column, env properties. shat:", shat(i,y), "qhat:", qhat(i,y), "su:", su_1(i,y), "sd", sd(i,y)
         !       write(iulog,*) "*** Going up column. updraft properties. su", su_1(i,y), "qu:", qu_1(i,y)
         !       write(iulog,*) "*** Going up column. downdraft properties. sd", sd(i,y), "qd:", qd(i,y)
         !       write(iulog,*) "*** Going up column. emc:", evpg(i,y)-cug(i,y), "latent heat:", -rl/cpres*(evpg(i,y)-cug(i,y))
         !       write(iulog,*) "updraft s leaving level:", mu(i,y)*(su_1(i,y)-shat(i,y)), "downdraft s joining level:", md(i,y)*(sd(i,y)-shat(i,y))
         !       write(iulog,*) "updraft q leaving level:", mu(i,y)*(qu_1(i,y)-qhat(i,y)), "downdraft q joining level:", md(i,y)*(qd(i,y)-qhat(i,y))
         !    end do
         !    do y = maxg_2(i), jt_2(i), -1
         !       write(iulog,*) " "
         !       ! write(iulog,*) "*** ZM_CONV: q1q2_pjr. k:", k, "dsdt:", dsdt(i,k), "old s:", sg(i,k), "new s:", sg(i,k)+1800._r8*dsdt(i,k)
         !       write(iulog,*) "*** ZM_CONV: q1q2_pjr. k:", y, "dsdt:", dsdt(i,y), "dqdt:", dqdt(i,y), "dl:", dlg(i,y)
         !       write(iulog,*) "*** Going up column. mu:", mu(i,y), "md:", md(i,y), "dp:", dp(i,y)
         !       write(iulog,*) "*** Going up column, env properties. shat:", shat(i,y), "qhat:", qhat(i,y), "su:", su_2(i,y), "sd", sd(i,y)
         !       write(iulog,*) "*** Going up column. updraft properties. su", su_2(i,y), "qu:", qu_2(i,y)
         !       write(iulog,*) "*** Going up column. downdraft properties. sd", sd(i,y), "qd:", qd(i,y)
         !       write(iulog,*) "*** Going up column. emc:", evpg(i,y)-cug(i,y), "latent heat:", -rl/cpres*(evpg(i,y)-cug(i,y))
         !       write(iulog,*) "updraft s leaving level:", mu(i,y)*(su_2(i,y)-shat(i,y)), "downdraft s joining level:", md(i,y)*(sd(i,y)-shat(i,y))
         !       write(iulog,*) "updraft q leaving level:", mu(i,y)*(qu_2(i,y)-qhat(i,y)), "downdraft q joining level:", md(i,y)*(qd(i,y)-qhat(i,y))
         !    end do
         !    write(iulog,*) " "
         ! end if

         if (dsdt(i,k)*cpres < max_dsdt) then
            max_dsdt = dsdt(i,k)*cpres
            max_dsdt_i = i
            max_dsdt_k = k
         end if         
         if (dqdt(i,k)<max_dqdt) then
            max_dqdt = dqdt(i,k)
            max_dqdt_i = i
            max_dqdt_k = k
         end if

      end do
   end do

   ! do i=1,ncol
   !    if (i==1) then
   !       write(iulog,*) " "
   !       write(iulog,*) "lchnk:", lchnk, "i:", i, "ZM Deep Convection"
   !       do k=pver,25, -1
   !          write(iulog,*) "k:", k, "old s:", s(i,k), "new s:", s(i,k)+heat(i,k)*2._r8*delt/cpres, "change:", heat(i,k)*2._r8*delt/cpres
   !       end do
   !    end if
   ! end do

   ! do i=1,ncol
   ! if (i==1) then
   !    write(iulog,*) "*** ZM_CONV: q1q2_pjr. At lchnk:", lchnk, "i: ", i
   !    write(iulog,*) "*** ZM_CONV: q1q2_pjr. mx_1:", maxg_1(i), "jt_1:", jt_1(i), "mx_2:", maxg_2(i), "jt_2:", jt_2(i)
   !    do y = maxg_1(i), jt_1(i), -1
   !       write(iulog,*) " "
   !       ! write(iulog,*) "*** ZM_CONV: q1q2_pjr. k:", k, "dsdt:", dsdt(i,k), "old s:", sg(i,k), "new s:", sg(i,k)+1800._r8*dsdt(i,k)
   !       write(iulog,*) "*** ZM_CONV: q1q2_pjr. k:", y, "dsdt:", dsdt(i,y), "dqdt:", dqdt(i,y), "dl:", dlg(i,y)
   !       write(iulog,*) "*** Going up column. mu:", mu(i,y), "md:", md(i,y), "dp:", dp(i,y)
   !       write(iulog,*) "*** Going up column, env properties. shat:", shat(i,y), "qhat:", qhat(i,y), "su:", su_1(i,y), "sd", sd(i,y)
   !       write(iulog,*) "*** Going up column. updraft properties. su", su_1(i,y), "qu:", qu_1(i,y)
   !       write(iulog,*) "*** Going up column. downdraft properties. sd", sd(i,y), "qd:", qd(i,y)
   !       write(iulog,*) "*** Going up column. emc:", evpg(i,y)-cug(i,y), "latent heat:", -rl/cpres*(evpg(i,y)-cug(i,y))
   !       write(iulog,*) "updraft s leaving level:", mu(i,y)*(su_1(i,y)-shat(i,y)), "downdraft s joining level:", md(i,y)*(sd(i,y)-shat(i,y))
   !       write(iulog,*) "updraft q leaving level:", mu(i,y)*(qu_1(i,y)-qhat(i,y)), "downdraft q joining level:", md(i,y)*(qd(i,y)-qhat(i,y))
   !    end do
   !    do y = maxg_2(i), jt_2(i), -1
   !       write(iulog,*) " "
   !       ! write(iulog,*) "*** ZM_CONV: q1q2_pjr. k:", k, "dsdt:", dsdt(i,k), "old s:", sg(i,k), "new s:", sg(i,k)+1800._r8*dsdt(i,k)
   !       write(iulog,*) "*** ZM_CONV: q1q2_pjr. k:", y, "dsdt:", dsdt(i,y), "dqdt:", dqdt(i,y), "dl:", dlg(i,y)
   !       write(iulog,*) "*** Going up column. mu:", mu(i,y), "md:", md(i,y), "dp:", dp(i,y)
   !       write(iulog,*) "*** Going up column, env properties. shat:", shat(i,y), "qhat:", qhat(i,y), "su:", su_2(i,y), "sd", sd(i,y)
   !       write(iulog,*) "*** Going up column. updraft properties. su", su_2(i,y), "qu:", qu_2(i,y)
   !       write(iulog,*) "*** Going up column. downdraft properties. sd", sd(i,y), "qd:", qd(i,y)
   !       write(iulog,*) "*** Going up column. emc:", evpg(i,y)-cug(i,y), "latent heat:", -rl/cpres*(evpg(i,y)-cug(i,y))
   !       write(iulog,*) "updraft s leaving level:", mu(i,y)*(su_2(i,y)-shat(i,y)), "downdraft s joining level:", md(i,y)*(sd(i,y)-shat(i,y))
   !       write(iulog,*) "updraft q leaving level:", mu(i,y)*(qu_2(i,y)-qhat(i,y)), "downdraft q joining level:", md(i,y)*(qd(i,y)-qhat(i,y))
   !    end do
   ! end if
   ! end do

!if (lchnk==33) then
! write(iulog,*) "lchnk:", lchnk
! write(iulog,*) "max mc:", maxval(mc(1:lengath, :)), "at loc:", maxloc(mc(1:lengath, :))
! write(iulog,*) "max dsdt:", maxval(dsdt(1:lengath, :))*cpres, "at loc:", maxloc(dsdt(1:lengath, :))
! write(iulog,*) "min dsdt:", minval(dsdt(1:lengath, :))*cpres, "at loc:", minloc(dsdt(1:lengath, :))
! write(iulog,*) "max dqdt:", maxval(dqdt(1:lengath, :)), "at loc:", maxloc(dqdt(1:lengath, :))
! write(iulog,*) "min dqdt:", minval(dqdt(1:lengath, :)), "at loc:", minloc(dqdt(1:lengath, :))
! write(iulog,*) " "
!endif
! write(iulog,*) "max mc:", max_mc, "at i:", max_mc_i, "k:", max_mc_k
! write(iulog,*) "min dsdt:", max_dsdt, "at i:", max_dsdt_i, "k:", max_dsdt_k
! write(iulog,*) "min dqdt:", max_dqdt, "at i:", max_dqdt_i, "k:", max_dqdt_k
! write(iulog,*) " "
! if (maxval(heat_2)>100._r8)  write(iulog,*) "lchnk:", lchnk, "max dsdt:", maxval(heat_2), "at loc:", maxloc(heat_2)
! if (minval(heat_2)<-100._r8) write(iulog,*) "lchnk:", lchnk, "min dsdt:", minval(heat_2), "at loc:", minloc(heat_2) !, "dsdt:", heat_2(13,:)
! if (maxval(qtnd_2)>1._r8)  write(iulog,*) "lchnk:", lchnk, "max dqdt:", maxval(qtnd_2), "at loc:", maxloc(qtnd_2) !, "dqdt:", qtnd_2(13,:)
! if (minval(qtnd_2)<-1._r8) write(iulog,*) "lchnk:", lchnk, "min dqdt:", minval(qtnd_2), "at loc:", minloc(qtnd_2)
! write(iulog,*) " "
!
!DIR$ CONCURRENT
   do i = 1,lengath
      jctop_1(ideep(i)) = jt_1(i)
      jcbot_1(ideep(i)) = maxg_1(i)
      jctop_2(ideep(i)) = jt_2(i)
      jcbot_2(ideep(i)) = maxg_2(i)
!--bee
      pflx(ideep(i),pverp) = pflxg(i,pverp)
      lcl_reached(ideep(i)) = lcl_reachedg(i)
   end do

! Compute precip by integrating change in water vapor minus detrained cloud water
   do k = pver,msg + 1,-1
      do i = 1,ncol
         prec(i) = prec(i) - dpp(i,k)* (q(i,k)-qh(i,k,1)) - dpp(i,k)*dlf(i,k)*2*delt
      end do
   end do

! obtain final precipitation rate in m/s.
   do i = 1,ncol
      prec(i) = rgrav*max(prec(i),0._r8)/ (2._r8*delt)/1000._r8
   end do

! Compute reserved liquid (not yet in cldliq) for energy integrals.
! Treat rliq as flux out bottom, to be added back later.
   do k = 1, pver
      do i = 1, ncol
         rliq(i) = rliq(i) + dlf(i,k)*dpp(i,k)/gravit
      end do
   end do
   rliq(:ncol) = rliq(:ncol) /1000._r8

   capeg_2(:) = 0._r8
   lelg_2(:) = 2
   maxg_2(:) = 2

   return
end subroutine zm_convr

!===============================================================================
subroutine zm_conv_evap(ncol,lchnk, &
     t,pmid,pdel,q, &
     tend_s, tend_s_snwprd, tend_s_snwevmlt, tend_q, &
     prdprec, cldfrc, deltat,  &
     prec, snow, ntprprd, ntsnprd, flxprec, flxsnow )

!-----------------------------------------------------------------------
! Compute tendencies due to evaporation of rain from ZM scheme
!--
! Compute the total precipitation and snow fluxes at the surface.
! Add in the latent heat of fusion for snow formation and melt, since it not dealt with
! in the Zhang-MacFarlane parameterization.
! Evaporate some of the precip directly into the environment using a Sundqvist type algorithm
!-----------------------------------------------------------------------

    use wv_saturation,  only: qsat
    use phys_grid, only: get_rlat_all_p

!------------------------------Arguments--------------------------------
    integer,intent(in) :: ncol, lchnk             ! number of columns and chunk index
    real(r8),intent(in), dimension(pcols,pver) :: t          ! temperature (K)
    real(r8),intent(in), dimension(pcols,pver) :: pmid       ! midpoint pressure (Pa) 
    real(r8),intent(in), dimension(pcols,pver) :: pdel       ! layer thickness (Pa)
    real(r8),intent(in), dimension(pcols,pver) :: q          ! water vapor (kg/kg)
    real(r8),intent(inout), dimension(pcols,pver) :: tend_s     ! heating rate (J/kg/s)
    real(r8),intent(inout), dimension(pcols,pver) :: tend_q     ! water vapor tendency (kg/kg/s)
    real(r8),intent(out  ), dimension(pcols,pver) :: tend_s_snwprd ! Heating rate of snow production
    real(r8),intent(out  ), dimension(pcols,pver) :: tend_s_snwevmlt ! Heating rate of evap/melting of snow
    


    real(r8), intent(in   ) :: prdprec(pcols,pver)! precipitation production (kg/ks/s)
    real(r8), intent(in   ) :: cldfrc(pcols,pver) ! cloud fraction
    real(r8), intent(in   ) :: deltat             ! time step

    real(r8), intent(inout) :: prec(pcols)        ! Convective-scale preciptn rate
    real(r8), intent(out)   :: snow(pcols)        ! Convective-scale snowfall rate
!
!---------------------------Local storage-------------------------------

    real(r8) :: est    (pcols,pver)    ! Saturation vapor pressure
    real(r8) :: fice   (pcols,pver)    ! ice fraction in precip production
    real(r8) :: fsnow_conv(pcols,pver) ! snow fraction in precip production
    real(r8) :: qs     (pcols,pver)    ! saturation specific humidity
    real(r8),intent(out) :: flxprec(pcols,pverp)   ! Convective-scale flux of precip at interfaces (kg/m2/s)
    real(r8),intent(out) :: flxsnow(pcols,pverp)   ! Convective-scale flux of snow   at interfaces (kg/m2/s)
    real(r8),intent(out) :: ntprprd(pcols,pver)    ! net precip production in layer
    real(r8),intent(out) :: ntsnprd(pcols,pver)    ! net snow production in layer
    real(r8) :: work1                  ! temp variable (pjr)
    real(r8) :: work2                  ! temp variable (pjr)

    real(r8) :: evpvint(pcols)         ! vertical integral of evaporation
    real(r8) :: evpprec(pcols)         ! evaporation of precipitation (kg/kg/s)
    real(r8) :: evpsnow(pcols)         ! evaporation of snowfall (kg/kg/s)
    real(r8) :: snowmlt(pcols)         ! snow melt tendency in layer
    real(r8) :: flxsntm(pcols)         ! flux of snow into layer, after melting

    real(r8) :: evplimit               ! temp variable for evaporation limits
    real(r8) :: rlat(pcols)

    integer :: i,k                     ! longitude,level indices


!-----------------------------------------------------------------------

! convert input precip to kg/m2/s
    prec(:ncol) = prec(:ncol)*1000._r8

! determine saturation vapor pressure
    call qsat (t(1:ncol, 1:pver)    ,pmid(1:ncol, 1:pver),  &
               est(1:ncol, 1:pver)  ,qs(1:ncol, 1:pver))

! determine ice fraction in rain production (use cloud water parameterization fraction at present)
    call cldfrc_fice(ncol, t, fice, fsnow_conv)

! zero the flux integrals on the top boundary
    flxprec(:ncol,1) = 0._r8
    flxsnow(:ncol,1) = 0._r8
    evpvint(:ncol)   = 0._r8

    do k = 1, pver
       do i = 1, ncol

! Melt snow falling into layer, if necessary. 
          if (t(i,k) > tmelt) then
             flxsntm(i) = 0._r8
             snowmlt(i) = flxsnow(i,k) * gravit/ pdel(i,k)
          else
             flxsntm(i) = flxsnow(i,k)
             snowmlt(i) = 0._r8
          end if

! relative humidity depression must be > 0 for evaporation
          evplimit = max(1._r8 - q(i,k)/qs(i,k), 0._r8)

! total evaporation depends on flux in the top of the layer
! flux prec is the net production above layer minus evaporation into environmet
          evpprec(i) = ke * (1._r8 - cldfrc(i,k)) * evplimit * sqrt(flxprec(i,k))
!**********************************************************
!!          evpprec(i) = 0.    ! turn off evaporation for now
!**********************************************************

! Don't let evaporation supersaturate layer (approx). Layer may already be saturated.
! Currently does not include heating/cooling change to qsat
          evplimit   = max(0._r8, (qs(i,k)-q(i,k)) / deltat)

! Don't evaporate more than is falling into the layer - do not evaporate rain formed
! in this layer but if precip production is negative, remove from the available precip
! Negative precip production occurs because of evaporation in downdrafts.
!!$          evplimit   = flxprec(i,k) * gravit / pdel(i,k) + min(prdprec(i,k), 0.)
          evplimit   = min(evplimit, flxprec(i,k) * gravit / pdel(i,k))

! Total evaporation cannot exceed input precipitation
          evplimit   = min(evplimit, (prec(i) - evpvint(i)) * gravit / pdel(i,k))

          evpprec(i) = min(evplimit, evpprec(i))

! evaporation of snow depends on snow fraction of total precipitation in the top after melting
          if (flxprec(i,k) > 0._r8) then
!            evpsnow(i) = evpprec(i) * flxsntm(i) / flxprec(i,k)
!            prevent roundoff problems
             work1 = min(max(0._r8,flxsntm(i)/flxprec(i,k)),1._r8)
             evpsnow(i) = evpprec(i) * work1
          else
             evpsnow(i) = 0._r8
          end if

! vertically integrated evaporation
          evpvint(i) = evpvint(i) + evpprec(i) * pdel(i,k)/gravit

! net precip production is production - evaporation
          ntprprd(i,k) = prdprec(i,k) - evpprec(i)
! net snow production is precip production * ice fraction - evaporation - melting
!pjrworks ntsnprd(i,k) = prdprec(i,k)*fice(i,k) - evpsnow(i) - snowmlt(i)
!pjrwrks2 ntsnprd(i,k) = prdprec(i,k)*fsnow_conv(i,k) - evpsnow(i) - snowmlt(i)
! the small amount added to flxprec in the work1 expression has been increased from 
! 1e-36 to 8.64e-11 (1e-5 mm/day).  This causes the temperature based partitioning
! scheme to be used for small flxprec amounts.  This is to address error growth problems.
#ifdef PERGRO
          work1 = min(max(0._r8,flxsnow(i,k)/(flxprec(i,k)+8.64e-11_r8)),1._r8)
#else
          if (flxprec(i,k).gt.0._r8) then
             work1 = min(max(0._r8,flxsnow(i,k)/flxprec(i,k)),1._r8)
          else
             work1 = 0._r8
          endif
#endif
          work2 = max(fsnow_conv(i,k), work1)
          if (snowmlt(i).gt.0._r8) work2 = 0._r8
!         work2 = fsnow_conv(i,k)
          ntsnprd(i,k) = prdprec(i,k)*work2 - evpsnow(i) - snowmlt(i)
          tend_s_snwprd  (i,k) = prdprec(i,k)*work2*latice
          tend_s_snwevmlt(i,k) = - ( evpsnow(i) + snowmlt(i) )*latice

! precipitation fluxes
          flxprec(i,k+1) = flxprec(i,k) + ntprprd(i,k) * pdel(i,k)/gravit
          flxsnow(i,k+1) = flxsnow(i,k) + ntsnprd(i,k) * pdel(i,k)/gravit

! protect against rounding error
          flxprec(i,k+1) = max(flxprec(i,k+1), 0._r8)
          flxsnow(i,k+1) = max(flxsnow(i,k+1), 0._r8)
! more protection (pjr)
!         flxsnow(i,k+1) = min(flxsnow(i,k+1), flxprec(i,k+1))

! heating (cooling) and moistening due to evaporation 
! - latent heat of vaporization for precip production has already been accounted for
! - snow is contained in prec
          tend_s(i,k)   =-evpprec(i)*latvap + ntsnprd(i,k)*latice
          tend_q(i,k) = evpprec(i)
       end do
    end do

! set output precipitation rates (m/s)
    prec(:ncol) = flxprec(:ncol,pver+1) / 1000._r8
    snow(:ncol) = flxsnow(:ncol,pver+1) / 1000._r8

!**********************************************************
!!$    tend_s(:ncol,:)   = 0.      ! turn heating off
!**********************************************************

  end subroutine zm_conv_evap



subroutine convtran(lchnk   , &
                    doconvtran,q       ,ncnst   ,mu    ,md      , &
                    du      ,eu      ,ed      ,dp      , &
                    jt_1    ,mx_1    ,ideep   ,il1g    ,il2g    , &
                    nstep   ,fracis  ,dqdt    ,dpdry   ,jt_2    ,mx_2)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Convective transport of trace species
!
! Mixing ratios may be with respect to either dry or moist air
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: P. Rasch
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use constituents,    only: cnst_get_type_byind
   use ppgrid
   use abortutils, only: endrun

   implicit none
!-----------------------------------------------------------------------
!
! Input arguments
!
   integer, intent(in) :: lchnk                 ! chunk identifier
   integer, intent(in) :: ncnst                 ! number of tracers to transport
   logical, intent(in) :: doconvtran(ncnst)     ! flag for doing convective transport
   real(r8), intent(in) :: q(pcols,pver,ncnst)  ! Tracer array including moisture
   real(r8), intent(in) :: mu(pcols,pver)       ! Mass flux up
   real(r8), intent(in) :: md(pcols,pver)       ! Mass flux down
   real(r8), intent(in) :: du(pcols,pver)       ! Mass detraining from updraft
   real(r8), intent(in) :: eu(pcols,pver)       ! Mass entraining from updraft
   real(r8), intent(in) :: ed(pcols,pver)       ! Mass entraining from downdraft
   real(r8), intent(in) :: dp(pcols,pver)       ! Delta pressure between interfaces
   real(r8), intent(in) :: fracis(pcols,pver,ncnst) ! fraction of tracer that is insoluble

   integer, intent(in) :: jt_1(pcols)         ! Index of cloud top for each column - lower convective region
   integer, intent(in) :: jt_2(pcols)         ! Index of cloud top for each column - upper convective region
   integer, intent(in) :: mx_1(pcols)         ! Index of cloud top for each column - lower convective region
   integer, intent(in) :: mx_2(pcols)         ! Index of cloud top for each column - upper convective region
   integer, intent(in) :: ideep(pcols)      ! Gathering array
   integer, intent(in) :: il1g              ! Gathered min lon indices over which to operate
   integer, intent(in) :: il2g              ! Gathered max lon indices over which to operate
   integer, intent(in) :: nstep             ! Time step index

   real(r8), intent(in) :: dpdry(pcols,pver)       ! Delta pressure between interfaces


! input/output

   real(r8), intent(out) :: dqdt(pcols,pver,ncnst)  ! Tracer tendency array

!--------------------------Local Variables------------------------------

   integer i                 ! Work index
   integer k                 ! Work index
   integer kbm               ! Highest altitude index of cloud base
   integer kk                ! Work index
   integer kkp1              ! Work index
   integer km1               ! Work index
   integer kp1               ! Work index
   integer ktm               ! Highest altitude index of cloud top
   integer m                 ! Work index

   real(r8) cabv                 ! Mix ratio of constituent above
   real(r8) cbel                 ! Mix ratio of constituent below
   real(r8) cdifr                ! Normalized diff between cabv and cbel
   real(r8) chat(pcols,pver)     ! Mix ratio in env at interfaces
   real(r8) cond(pcols,pver)     ! Mix ratio in downdraft at interfaces
   real(r8) const(pcols,pver)    ! Gathered tracer array
   real(r8) fisg(pcols,pver)     ! gathered insoluble fraction of tracer
   real(r8) conu(pcols,pver)     ! Mix ratio in updraft at interfaces
   real(r8) dcondt(pcols,pver)   ! Gathered tend array
   real(r8) small                ! A small number
   real(r8) mbsth                ! Threshold for mass fluxes
   real(r8) mupdudp              ! A work variable
   real(r8) minc                 ! A work variable
   real(r8) maxc                 ! A work variable
   real(r8) fluxin               ! A work variable
   real(r8) fluxout              ! A work variable
   real(r8) netflux              ! A work variable

   real(r8) dutmp(pcols,pver)       ! Mass detraining from updraft
   real(r8) eutmp(pcols,pver)       ! Mass entraining from updraft
   real(r8) edtmp(pcols,pver)       ! Mass entraining from downdraft
   real(r8) dptmp(pcols,pver)    ! Delta pressure between interfaces
!-----------------------------------------------------------------------
!
   small = 1.e-36_r8
! mbsth is the threshold below which we treat the mass fluxes as zero (in mb/s)
   mbsth = 1.e-15_r8

! Find the highest level top and bottom levels of convection
   ktm = pver
   kbm = pver
   do i = il1g, il2g
      if (jt_2(i) == 2) then !left at default value, no upper convective region
         ktm = min(ktm,jt_1(i))
      else if (jt_1(i) == 2) then !no lower convective region
         ktm = min(ktm,jt_2(i))
      else
         ktm = min(ktm,min(jt_1(i),jt_2(i)))
      end if

      if (mx_2(i) == 2) then !left at default value, no upper convective region
         kbm = min(kbm,mx_1(i))
      else if (mx_1(i) == 2) then !no lower convective region (rare case?)
         kbm = min(kbm,mx_2(i))
      else
         kbm = min(kbm,min(mx_1(i),mx_2(i)))
      end if
   end do

! Loop ever each constituent
   do m = 2, ncnst
      if (doconvtran(m)) then

         if (cnst_get_type_byind(m).eq.'dry') then
            do k = 1,pver
               do i =il1g,il2g
                  dptmp(i,k) = dpdry(i,k)
                  dutmp(i,k) = du(i,k)*dp(i,k)/dpdry(i,k)
                  eutmp(i,k) = eu(i,k)*dp(i,k)/dpdry(i,k)
                  edtmp(i,k) = ed(i,k)*dp(i,k)/dpdry(i,k)
               end do
            end do
         else
            do k = 1,pver
               do i =il1g,il2g
                  dptmp(i,k) = dp(i,k)
                  dutmp(i,k) = du(i,k)
                  eutmp(i,k) = eu(i,k)
                  edtmp(i,k) = ed(i,k)
               end do
            end do
         endif
!        dptmp = dp

! Gather up the constituent and set tend to zero
         do k = 1,pver
            do i =il1g,il2g
               const(i,k) = q(ideep(i),k,m)
               fisg(i,k) = fracis(ideep(i),k,m)
            end do
         end do

! From now on work only with gathered data

! Interpolate environment tracer values to interfaces
         do k = 1,pver
            km1 = max(1,k-1)
            do i = il1g, il2g
               minc = min(const(i,km1),const(i,k))
               maxc = max(const(i,km1),const(i,k))
               if (minc < 0) then
                  cdifr = 0._r8
               else
                  cdifr = abs(const(i,k)-const(i,km1))/max(maxc,small)
               endif

! If the two layers differ significantly use a geometric averaging
! procedure
               if (cdifr > 1.E-6_r8) then
                  cabv = max(const(i,km1),maxc*1.e-12_r8)
                  cbel = max(const(i,k),maxc*1.e-12_r8)
                  chat(i,k) = log(cabv/cbel)/(cabv-cbel)*cabv*cbel

               else             ! Small diff, so just arithmetic mean
                  chat(i,k) = 0.5_r8* (const(i,k)+const(i,km1))
               end if

! Provisional up and down draft values
               conu(i,k) = chat(i,k)
               cond(i,k) = chat(i,k)

!              provisional tends
               dcondt(i,k) = 0._r8

            end do
         end do

! Do levels adjacent to top and bottom
         k = 2
         km1 = 1
         kk = pver
         do i = il1g,il2g
            mupdudp = mu(i,kk) + dutmp(i,kk)*dptmp(i,kk)
            if (mupdudp > mbsth) then
               conu(i,kk) = (+eutmp(i,kk)*fisg(i,kk)*const(i,kk)*dptmp(i,kk))/mupdudp
            endif
            if (md(i,k) < -mbsth) then
               cond(i,k) =  (-edtmp(i,km1)*fisg(i,km1)*const(i,km1)*dptmp(i,km1))/md(i,k)
            endif
         end do

! Updraft from bottom to top
         do kk = pver-1,1,-1
            kkp1 = min(pver,kk+1)
            do i = il1g,il2g
               mupdudp = mu(i,kk) + dutmp(i,kk)*dptmp(i,kk)
               if (mupdudp > mbsth) then
                  conu(i,kk) = (  mu(i,kkp1)*conu(i,kkp1)+eutmp(i,kk)*fisg(i,kk)* &
                                  const(i,kk)*dptmp(i,kk) )/mupdudp
               endif
            end do
         end do

! Downdraft from top to bottom
         do k = 3,pver
            km1 = max(1,k-1)
            do i = il1g,il2g
               if (md(i,k) < -mbsth) then
                  cond(i,k) =  (  md(i,km1)*cond(i,km1)-edtmp(i,km1)*fisg(i,km1)*const(i,km1) &
                                  *dptmp(i,km1) )/md(i,k)
               endif
            end do
         end do


         do k = ktm,pver
            km1 = max(1,k-1)
            kp1 = min(pver,k+1)
            do i = il1g,il2g

! version 1 hard to check for roundoff errors
!               dcondt(i,k) =
!     $                  +(+mu(i,kp1)* (conu(i,kp1)-chat(i,kp1))
!     $                    -mu(i,k)*   (conu(i,k)-chat(i,k))
!     $                    +md(i,kp1)* (cond(i,kp1)-chat(i,kp1))
!     $                    -md(i,k)*   (cond(i,k)-chat(i,k))
!     $                   )/dp(i,k)

! version 2 hard to limit fluxes
!               fluxin =  mu(i,kp1)*conu(i,kp1) + mu(i,k)*chat(i,k)
!     $                 -(md(i,k)  *cond(i,k)   + md(i,kp1)*chat(i,kp1))
!               fluxout = mu(i,k)*conu(i,k)     + mu(i,kp1)*chat(i,kp1)
!     $                 -(md(i,kp1)*cond(i,kp1) + md(i,k)*chat(i,k))

! version 3 limit fluxes outside convection to mass in appropriate layer
! these limiters are probably only safe for positive definite quantitities
! it assumes that mu and md already satify a courant number limit of 1
               fluxin =  mu(i,kp1)*conu(i,kp1)+ mu(i,k)*min(chat(i,k),const(i,km1)) &
                         -(md(i,k)  *cond(i,k) + md(i,kp1)*min(chat(i,kp1),const(i,kp1)))
               fluxout = mu(i,k)*conu(i,k) + mu(i,kp1)*min(chat(i,kp1),const(i,k)) &
                         -(md(i,kp1)*cond(i,kp1) + md(i,k)*min(chat(i,k),const(i,k)))

               netflux = fluxin - fluxout
               if (abs(netflux) < max(fluxin,fluxout)*1.e-12_r8) then
                  netflux = 0._r8
               endif
               dcondt(i,k) = netflux/dptmp(i,k)
            end do
         end do
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!DIR$ NOINTERCHANGE
         do k = kbm,pver
            km1 = max(1,k-1)
            do i = il1g,il2g
               if ((k == mx_1(i)) .or. (k == mx_2(i))) then

! version 1
!                  dcondt(i,k) = (1./dsubcld(i))*
!     $              (-mu(i,k)*(conu(i,k)-chat(i,k))
!     $               -md(i,k)*(cond(i,k)-chat(i,k))
!     $              )

! version 2
!                  fluxin =  mu(i,k)*chat(i,k) - md(i,k)*cond(i,k)
!                  fluxout = mu(i,k)*conu(i,k) - md(i,k)*chat(i,k)
! version 3
                  fluxin =  mu(i,k)*min(chat(i,k),const(i,km1)) - md(i,k)*cond(i,k)
                  fluxout = mu(i,k)*conu(i,k) - md(i,k)*min(chat(i,k),const(i,k))

                  netflux = fluxin - fluxout
                  if (abs(netflux) < max(fluxin,fluxout)*1.e-12_r8) then
                     netflux = 0._r8
                  endif
!                  dcondt(i,k) = netflux/dsubcld(i)
                  dcondt(i,k) = netflux/dptmp(i,k)

               else if (k > mx_1(i) .and. mx_2(i) ==2) then !only lower region
!                  dcondt(i,k) = dcondt(i,k-1)
                     dcondt(i,k) = 0._r8
               else if (k > mx_2(i) .and. mx_1(i) ==2) then !only upper region
!                  dcondt(i,k) = dcondt(i,k-1)
                     dcondt(i,k) = 0._r8
               else if ((k > mx_1(i)) .or. (k > mx_2(i) .and. k < jt_1(i)))  then! two convective regions, want to find points below mx_1 and between mx_2 and jt_1
                  dcondt(i,k) = 0._r8
               end if
            end do
         end do

! Initialize to zero everywhere, then scatter tendency back to full array
         dqdt(:,:,m) = 0._r8
         do k = 1,pver
            kp1 = min(pver,k+1)
!DIR$ CONCURRENT
            do i = il1g,il2g
               dqdt(ideep(i),k,m) = dcondt(i,k)
            end do
         end do

      end if      ! for doconvtran

   end do

   return
end subroutine convtran

!=========================================================================================

subroutine momtran(lchnk, ncol, &
                    domomtran,q       ,ncnst   ,mu      ,md    , &
                    du      ,eu      ,ed      ,dp      , &
                    jt_1    ,mx_1    ,jt_2    ,mx_2    ,ideep   ,il1g    ,il2g    , &
                    nstep   ,dqdt    ,pguall     ,pgdall, icwu, icwd, dt, seten    )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Convective transport of momentum
!
! Mixing ratios may be with respect to either dry or moist air
! 
! Method: 
! Based on the convtran subroutine by P. Rasch
! <Also include any applicable external references.> 
! 
! Author: J. Richter and P. Rasch
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use constituents,    only: cnst_get_type_byind
   use ppgrid
   use abortutils, only: endrun

   implicit none
!-----------------------------------------------------------------------
!
! Input arguments
!
   integer, intent(in) :: lchnk                 ! chunk identifier
   integer, intent(in) :: ncol                  ! number of atmospheric columns
   integer, intent(in) :: ncnst                 ! number of tracers to transport
   logical, intent(in) :: domomtran(ncnst)      ! flag for doing convective transport
   real(r8), intent(in) :: q(pcols,pver,ncnst)  ! Wind array
   real(r8), intent(in) :: mu(pcols,pver)       ! Mass flux up
   real(r8), intent(in) :: md(pcols,pver)       ! Mass flux down
   real(r8), intent(in) :: du(pcols,pver)       ! Mass detraining from updraft
   real(r8), intent(in) :: eu(pcols,pver)       ! Mass entraining from updraft
   real(r8), intent(in) :: ed(pcols,pver)       ! Mass entraining from downdraft
   real(r8), intent(in) :: dp(pcols,pver)       ! Delta pressure between interfaces
   real(r8), intent(in) :: dt                       !  time step in seconds : 2*delta_t

   integer, intent(in) :: jt_1(pcols)         ! Index of cloud top for each column - lower convective region
   integer, intent(in) :: jt_2(pcols)         ! Index of cloud top for each column - upper convective region
   integer, intent(in) :: mx_1(pcols)         ! Index of cloud top for each column - lower convective region
   integer, intent(in) :: mx_2(pcols)         ! Index of cloud top for each column - upper convective region
   integer, intent(in) :: ideep(pcols)      ! Gathering array
   integer, intent(in) :: il1g              ! Gathered min lon indices over which to operate
   integer, intent(in) :: il2g              ! Gathered max lon indices over which to operate
   integer, intent(in) :: nstep             ! Time step index



! input/output

   real(r8), intent(out) :: dqdt(pcols,pver,ncnst)  ! Tracer tendency array

!--------------------------Local Variables------------------------------

   integer i                 ! Work index
   integer k                 ! Work index
   integer kbm               ! Highest altitude index of cloud base
   integer kk                ! Work index
   integer kkp1              ! Work index
   integer kkm1              ! Work index
   integer km1               ! Work index
   integer kp1               ! Work index
   integer ktm               ! Highest altitude index of cloud top
   integer m                 ! Work index
   integer ii                 ! Work index

   real(r8) cabv                 ! Mix ratio of constituent above
   real(r8) cbel                 ! Mix ratio of constituent below
   real(r8) cdifr                ! Normalized diff between cabv and cbel
   real(r8) chat(pcols,pver)     ! Mix ratio in env at interfaces
   real(r8) cond(pcols,pver)     ! Mix ratio in downdraft at interfaces
   real(r8) const(pcols,pver)    ! Gathered wind array
   real(r8) conu(pcols,pver)     ! Mix ratio in updraft at interfaces
   real(r8) dcondt(pcols,pver)   ! Gathered tend array
   real(r8) small                ! A small number
   real(r8) mbsth                ! Threshold for mass fluxes
   real(r8) mupdudp              ! A work variable
   real(r8) minc                 ! A work variable
   real(r8) maxc                 ! A work variable
   real(r8) fluxin               ! A work variable
   real(r8) fluxout              ! A work variable
   real(r8) netflux              ! A work variable

   real(r8) momcu                ! constant for updraft pressure gradient term
   real(r8) momcd                ! constant for downdraft pressure gradient term
   real(r8) sum                  ! sum
   real(r8) sum2                  ! sum2
 
   real(r8) mududp(pcols,pver) ! working variable
   real(r8) mddudp(pcols,pver)     ! working variable

   real(r8) pgu(pcols,pver)      ! Pressure gradient term for updraft
   real(r8) pgd(pcols,pver)      ! Pressure gradient term for downdraft

   real(r8),intent(out) ::  pguall(pcols,pver,ncnst)      ! Apparent force from  updraft PG
   real(r8),intent(out) ::  pgdall(pcols,pver,ncnst)      ! Apparent force from  downdraft PG

   real(r8),intent(out) ::  icwu(pcols,pver,ncnst)      ! In-cloud winds in updraft
   real(r8),intent(out) ::  icwd(pcols,pver,ncnst)      ! In-cloud winds in downdraft

   real(r8),intent(out) ::  seten(pcols,pver) ! Dry static energy tendency
   real(r8)                 gseten(pcols,pver) ! Gathered dry static energy tendency

   real(r8)  mflux(pcols,pverp,ncnst)   ! Gathered momentum flux

   real(r8)  wind0(pcols,pver,ncnst)       !  gathered  wind before time step
   real(r8)  windf(pcols,pver,ncnst)       !  gathered  wind after time step
   real(r8) fkeb, fket, ketend_cons, ketend, utop, ubot, vtop, vbot, gset2
   

!-----------------------------------------------------------------------
!

! Initialize outgoing fields
   pguall(:,:,:)     = 0.0_r8
   pgdall(:,:,:)     = 0.0_r8
! Initialize in-cloud winds to environmental wind
   icwu(:ncol,:,:)       = q(:ncol,:,:)
   icwd(:ncol,:,:)       = q(:ncol,:,:)

! Initialize momentum flux and  final winds
   mflux(:,:,:)       = 0.0_r8
   wind0(:,:,:)         = 0.0_r8
   windf(:,:,:)         = 0.0_r8

! Initialize dry static energy

   seten(:,:)         = 0.0_r8
   gseten(:,:)         = 0.0_r8

! Define constants for parameterization

   momcu = 0.4_r8
   momcd = 0.4_r8

   small = 1.e-36_r8
! mbsth is the threshold below which we treat the mass fluxes as zero (in mb/s)
   mbsth = 1.e-15_r8

! Find the highest level top and bottom levels of convection
   ktm = pver
   kbm = pver
   do i = il1g, il2g
      if (jt_2(i) == 2) then !left at default value, no upper convective region
         ktm = min(ktm,jt_1(i))
      else if (jt_1(i) == 2) then !no lower convective region
         ktm = min(ktm,jt_2(i))
      else
         ktm = min(ktm,min(jt_1(i),jt_2(i)))
      end if

      if (mx_2(i) == 2) then !left at default value, no upper convective region
         kbm = min(kbm,mx_1(i))
      else if (mx_1(i) == 2) then !no lower convective region (rare case?)
         kbm = min(kbm,mx_2(i))
      else
         kbm = min(kbm,min(mx_1(i),mx_2(i)))
      end if
   end do

! Loop ever each wind component
   do m = 1, ncnst                    !start at m = 1 to transport momentum
      if (domomtran(m)) then

! Gather up the winds and set tend to zero
         do k = 1,pver
            do i =il1g,il2g
               const(i,k) = q(ideep(i),k,m)
                wind0(i,k,m) = const(i,k)
            end do
         end do


! From now on work only with gathered data

! Interpolate winds to interfaces

         do k = 1,pver
            km1 = max(1,k-1)
            do i = il1g, il2g

               ! use arithmetic mean
               chat(i,k) = 0.5_r8* (const(i,k)+const(i,km1))

! Provisional up and down draft values
               conu(i,k) = chat(i,k)
               cond(i,k) = chat(i,k)

!              provisional tends
               dcondt(i,k) = 0._r8

            end do
         end do


!
! Pressure Perturbation Term
! 

      !Top boundary:  assume mu is zero 

         k=1
         pgu(:il2g,k) = 0.0_r8
         pgd(:il2g,k) = 0.0_r8

         do k=2,pver-1
            km1 = max(1,k-1)
            kp1 = min(pver,k+1)
            do i = il1g,il2g
            
               !interior points

               mududp(i,k) =  ( mu(i,k) * (const(i,k)- const(i,km1))/dp(i,km1) &
                           +  mu(i,kp1) * (const(i,kp1) - const(i,k))/dp(i,k))

               pgu(i,k) = - momcu * 0.5_r8 * mududp(i,k)
                           

               mddudp(i,k) =  ( md(i,k) * (const(i,k)- const(i,km1))/dp(i,km1) &
                           +  md(i,kp1) * (const(i,kp1) - const(i,k))/dp(i,k))

               pgd(i,k) = - momcd * 0.5_r8 * mddudp(i,k)


            end do
         end do

       ! bottom boundary 
       k = pver
       km1 = max(1,k-1)
       do i=il1g,il2g

          mududp(i,k) =   mu(i,k) * (const(i,k)- const(i,km1))/dp(i,km1)
          pgu(i,k) = - momcu *  mududp(i,k)
          
          mddudp(i,k) =   md(i,k) * (const(i,k)- const(i,km1))/dp(i,km1) 

          pgd(i,k) = - momcd * mddudp(i,k)
          
       end do
       

!
! In-cloud velocity calculations
!

! Do levels adjacent to top and bottom
         k = 2
         km1 = 1
         kk = pver
         kkm1 = max(1,kk-1)
         do i = il1g,il2g
            mupdudp = mu(i,kk) + du(i,kk)*dp(i,kk)
            if (mupdudp > mbsth) then
                 
               conu(i,kk) = (+eu(i,kk)*const(i,kk)*dp(i,kk)+pgu(i,kk)*dp(i,kk))/mupdudp
            endif
            if (md(i,k) < -mbsth) then
               cond(i,k) =  (-ed(i,km1)*const(i,km1)*dp(i,km1))-pgd(i,km1)*dp(i,km1)/md(i,k)
            endif

                        
         end do



! Updraft from bottom to top
         do kk = pver-1,1,-1
            kkm1 = max(1,kk-1)
            kkp1 = min(pver,kk+1)
            do i = il1g,il2g
               mupdudp = mu(i,kk) + du(i,kk)*dp(i,kk)
               if (mupdudp > mbsth) then
            
                  conu(i,kk) = (  mu(i,kkp1)*conu(i,kkp1)+eu(i,kk)* &
                                  const(i,kk)*dp(i,kk)+pgu(i,kk)*dp(i,kk))/mupdudp
               endif
            end do

         end do


! Downdraft from top to bottom
         do k = 3,pver
            km1 = max(1,k-1)
            do i = il1g,il2g
               if (md(i,k) < -mbsth) then
                            
                  cond(i,k) =  (  md(i,km1)*cond(i,km1)-ed(i,km1)*const(i,km1) &
                                  *dp(i,km1)-pgd(i,km1)*dp(i,km1) )/md(i,k)

               endif
            end do
         end do


         sum = 0._r8
         sum2 = 0._r8


         do k = ktm,pver
            km1 = max(1,k-1)
            kp1 = min(pver,k+1)
            do i = il1g,il2g
               ii = ideep(i)
	
! version 1 hard to check for roundoff errors
               dcondt(i,k) =  &
                           +(mu(i,kp1)* (conu(i,kp1)-chat(i,kp1)) &
                           -mu(i,k)*   (conu(i,k)-chat(i,k))      &
                           +md(i,kp1)* (cond(i,kp1)-chat(i,kp1)) &
                           -md(i,k)*   (cond(i,k)-chat(i,k)) &
                          )/dp(i,k)

            end do
         end do

  ! dcont for bottom layer
          !
          !DIR$ NOINTERCHANGE
          do k = kbm,pver
             km1 = max(1,k-1)
             do i = il1g,il2g
               if (k == mx_1(i) .or. k == mx_2(i)) then
                   ! version 1
                   dcondt(i,k) = (1./dp(i,k))*   &  
                        (-mu(i,k)*(conu(i,k)-chat(i,k)) &
                        -md(i,k)*(cond(i,k)-chat(i,k)) &
                        )
                end if
             end do
          end do

! Initialize to zero everywhere, then scatter tendency back to full array
         dqdt(:,:,m) = 0._r8

         do k = 1,pver
            do i = il1g,il2g
               ii = ideep(i)
               dqdt(ii,k,m) = dcondt(i,k)
    ! Output apparent force on the mean flow from pressure gradient
               pguall(ii,k,m) = -pgu(i,k)
               pgdall(ii,k,m) = -pgd(i,k)
               icwu(ii,k,m)   =  conu(i,k)
               icwd(ii,k,m)   =  cond(i,k)
            end do
         end do

          ! Calculate momentum flux in units of mb*m/s2 

          do k = ktm,pver
             do i = il1g,il2g
                ii = ideep(i)
                mflux(i,k,m) = &
                     -mu(i,k)*   (conu(i,k)-chat(i,k))      &
                     -md(i,k)*   (cond(i,k)-chat(i,k))
             end do
          end do


          ! Calculate winds at the end of the time step 

          do k = ktm,pver
             do i = il1g,il2g
                ii = ideep(i)
                km1 = max(1,k-1)
                kp1 = k+1
                windf(i,k,m) = const(i,k)    -   (mflux(i,kp1,m) - mflux(i,k,m)) * dt /dp(i,k)

             end do
          end do

       end if      ! for domomtran
   end do

 ! Need to add an energy fix to account for the dissipation of kinetic energy
    ! Formulation follows from Boville and Bretherton (2003)
    ! formulation by PJR

    do k = ktm,pver
       km1 = max(1,k-1)
       kp1 = min(pver,k+1)
       do i = il1g,il2g

          ii = ideep(i)

          ! calculate the KE fluxes at top and bot of layer 
          ! based on a discrete approximation to b&b eq(35) F_KE = u*F_u + v*F_v at interface
          utop = (wind0(i,k,1)+wind0(i,km1,1))/2.
          vtop = (wind0(i,k,2)+wind0(i,km1,2))/2.
          ubot = (wind0(i,kp1,1)+wind0(i,k,1))/2.
          vbot = (wind0(i,kp1,2)+wind0(i,k,2))/2.
          fket = utop*mflux(i,k,1)   + vtop*mflux(i,k,2)    ! top of layer
          fkeb = ubot*mflux(i,k+1,1) + vbot*mflux(i,k+1,2)  ! bot of layer

          ! divergence of these fluxes should give a conservative redistribution of KE
          ketend_cons = (fket-fkeb)/dp(i,k)

          ! tendency in kinetic energy resulting from the momentum transport
          ketend = ((windf(i,k,1)**2 + windf(i,k,2)**2) - (wind0(i,k,1)**2 + wind0(i,k,2)**2))*0.5/dt

          ! the difference should be the dissipation
          gset2 = ketend_cons - ketend
          gseten(i,k) = gset2

       end do

    end do

    ! Scatter dry static energy to full array
    do k = 1,pver
       do i = il1g,il2g
          ii = ideep(i)
          seten(ii,k) = gseten(i,k)

       end do
    end do

   return
end subroutine momtran

!=========================================================================================

subroutine buoyan(lchnk   ,ncol    , &
                  q       ,t       ,p       ,z       ,pf      , &
                  tp      ,qstp    ,tl      ,rl      ,cape    , &
                  pblt    ,lcl     ,lel     ,lon     ,mx      , &
                  rd      ,grav    ,cp      ,msg     , &
                  tpert   )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! <Say what the routine does> 
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author:
! This is contributed code not fully standardized by the CCM core group.
! The documentation has been enhanced to the degree that we are able.
! Reviewed:          P. Rasch, April 1996
! 
!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------------
!
! input arguments
!
   integer, intent(in) :: lchnk                 ! chunk identifier
   integer, intent(in) :: ncol                  ! number of atmospheric columns

   real(r8), intent(in) :: q(pcols,pver)        ! spec. humidity
   real(r8), intent(in) :: t(pcols,pver)        ! temperature
   real(r8), intent(in) :: p(pcols,pver)        ! pressure
   real(r8), intent(in) :: z(pcols,pver)        ! height
   real(r8), intent(in) :: pf(pcols,pver+1)     ! pressure at interfaces
   real(r8), intent(in) :: pblt(pcols)          ! index of pbl depth
   real(r8), intent(in) :: tpert(pcols)         ! perturbation temperature by pbl processes

!
! output arguments
!
   real(r8), intent(out) :: tp(pcols,pver)       ! parcel temperature
   real(r8), intent(out) :: qstp(pcols,pver)     ! saturation mixing ratio of parcel
   real(r8), intent(out) :: tl(pcols)            ! parcel temperature at lcl
   real(r8), intent(out) :: cape(pcols)          ! convective aval. pot. energy.
   integer lcl(pcols)        !
   integer lel(pcols)        !
   integer lon(pcols)        ! level of onset of deep convection
   integer mx(pcols)         ! level of max moist static energy
!
!--------------------------Local Variables------------------------------
!
   real(r8) capeten(pcols,5)     ! provisional value of cape
   real(r8) tv(pcols,pver)       !
   real(r8) tpv(pcols,pver)      !
   real(r8) buoy(pcols,pver)

   real(r8) a1(pcols)
   real(r8) a2(pcols)
   real(r8) estp(pcols)
   real(r8) pl(pcols)
   real(r8) plexp(pcols)
   real(r8) hmax(pcols)
   real(r8) hmn(pcols)
   real(r8) y(pcols)

   logical plge600(pcols)
   integer knt(pcols)
   integer lelten(pcols,5)

   real(r8) cp
   real(r8) e
   real(r8) grav
   real(r8) zvirp1

   integer i
   integer k
   integer msg
   integer n

   real(r8) rd
   real(r8) rl
#ifdef PERGRO
   real(r8) rhd
#endif
!
!-----------------------------------------------------------------------
!
   zvirp1 = zvir + 1.0

   do n = 1,5
      do i = 1,ncol
         lelten(i,n) = pver
         capeten(i,n) = 0._r8
      end do
   end do
!
   do i = 1,ncol
      lon(i) = pver
      knt(i) = 0
      lel(i) = pver
      mx(i) = lon(i)
      cape(i) = 0._r8
      hmax(i) = 0._r8
   end do

   tp(:ncol,:) = t(:ncol,:)
   qstp(:ncol,:) = q(:ncol,:)

!!! RBN - Initialize tv and buoy for output.
!!! tv=tv : tpv=tpv : qstp=q : buoy=0.
   tv(:ncol,:) = t(:ncol,:) *(1._r8+zvirp1*q(:ncol,:))/ (1._r8+q(:ncol,:))
   tpv(:ncol,:) = tv(:ncol,:)
   buoy(:ncol,:) = 0._r8

!
! set "launching" level(mx) to be at maximum moist static energy.
! search for this level stops at planetary boundary layer top.
!
#ifdef PERGRO
   do k = pver,msg + 1,-1
      do i = 1,ncol
         hmn(i) = cp*t(i,k) + grav*z(i,k) + rl*q(i,k)
!
! Reset max moist static energy level when relative difference exceeds 1.e-4
!
         rhd = (hmn(i) - hmax(i))/(hmn(i) + hmax(i))
         if (k >= nint(pblt(i)) .and. k <= lon(i) .and. rhd > -1.e-4_r8) then
            hmax(i) = hmn(i)
            mx(i) = k
         end if
      end do
   end do
#else
   do k = pver,msg + 1,-1
      do i = 1,ncol
         hmn(i) = cp*t(i,k) + grav*z(i,k) + rl*q(i,k)
         if (k >= nint(pblt(i)) .and. k <= lon(i) .and. hmn(i) > hmax(i)) then
            hmax(i) = hmn(i)
            mx(i) = k
         end if
      end do
   end do
#endif
!
   do i = 1,ncol
      lcl(i) = mx(i)
      e = p(i,mx(i))*q(i,mx(i))/ (eps1+q(i,mx(i)))
      tl(i) = 2840._r8/ (3.5_r8*log(t(i,mx(i)))-log(e)-4.805_r8) + 55._r8
      if (tl(i) < t(i,mx(i))) then
         plexp(i) = (1._r8/ (0.2854_r8* (1._r8-0.28_r8*q(i,mx(i)))))
         pl(i) = p(i,mx(i))* (tl(i)/t(i,mx(i)))**plexp(i)
      else
         tl(i) = t(i,mx(i))
         pl(i) = p(i,mx(i))
      end if
   end do

!
! calculate lifting condensation level (lcl).
!
   do k = pver,msg + 2,-1
      do i = 1,ncol
         if (k <= mx(i) .and. (p(i,k) > pl(i) .and. p(i,k-1) <= pl(i))) then
            lcl(i) = k - 1
         end if
      end do
   end do
!
! if lcl is above the nominal level of non-divergence (600 mbs),
! no deep convection is permitted (ensuing calculations
! skipped and cape retains initialized value of zero).
!
   do i = 1,ncol
      plge600(i) = pl(i).ge.600._r8
   end do
!
! initialize parcel properties in sub-cloud layer below lcl.
!
   do k = pver,msg + 1,-1
      do i=1,ncol
         if (k > lcl(i) .and. k <= mx(i) .and. plge600(i)) then
            tv(i,k) = t(i,k)* (1._r8+zvirp1*q(i,k))/ (1._r8+q(i,k))
            qstp(i,k) = q(i,mx(i))
            tp(i,k) = t(i,mx(i))* (p(i,k)/p(i,mx(i)))**(0.2854_r8* (1._r8-0.28_r8*q(i,mx(i))))
!
! buoyancy is increased by 0.5 k as in tiedtke
!
!-jjh          tpv (i,k)=tp(i,k)*(1.+zvirp1*q(i,mx(i)))/
!-jjh     1                     (1.+q(i,mx(i)))
            tpv(i,k) = (tp(i,k)+tpert(i))*(1._r8+zvirp1*q(i,mx(i)))/ (1._r8+q(i,mx(i)))
            buoy(i,k) = tpv(i,k) - tv(i,k) + tiedke_add
         end if
      end do
   end do

!
! define parcel properties at lcl (i.e. level immediately above pl).
!
   do k = pver,msg + 1,-1
      do i=1,ncol
         if (k == lcl(i) .and. plge600(i)) then
            tv(i,k) = t(i,k)* (1._r8+zvirp1*q(i,k))/ (1._r8+q(i,k))
            qstp(i,k) = q(i,mx(i))
            tp(i,k) = tl(i)* (p(i,k)/pl(i))**(0.2854_r8* (1._r8-0.28_r8*qstp(i,k)))
!              estp(i)  =exp(a-b/tp(i,k))
! use of different formulas for est has about 1 g/kg difference
! in qs at t= 300k, and 0.02 g/kg at t=263k, with the formula
! above giving larger qs.
!
            estp(i) = c1*exp((c2* (tp(i,k)-tfreez))/((tp(i,k)-tfreez)+c3))
            call esat_buck(tp(i,k), estp(i))

            qstp(i,k) = eps1*estp(i)/ (p(i,k)-estp(i))
            a1(i) = cp / rl + qstp(i,k) * (1._r8+ qstp(i,k) / eps1) * rl * eps1 / &
                    (rd * tp(i,k) ** 2)
            a2(i) = .5_r8* (qstp(i,k)* (1._r8+2._r8/eps1*qstp(i,k))* &
                    (1._r8+qstp(i,k)/eps1)*eps1**2*rl*rl/ &
                    (rd**2*tp(i,k)**4)-qstp(i,k)* &
                    (1._r8+qstp(i,k)/eps1)*2._r8*eps1*rl/ &
                    (rd*tp(i,k)**3))
            a1(i) = 1._r8/a1(i)
            a2(i) = -a2(i)*a1(i)**3
            y(i) = q(i,mx(i)) - qstp(i,k)
            tp(i,k) = tp(i,k) + a1(i)*y(i) + a2(i)*y(i)**2
!          estp(i)  =exp(a-b/tp(i,k))
           estp(i) = c1*exp((c2* (tp(i,k)-tfreez))/ ((tp(i,k)-tfreez)+c3))
           call esat_buck(tp(i,k), estp(i))

            qstp(i,k) = eps1*estp(i) / (p(i,k)-estp(i))
!
! buoyancy is increased by 0.5 k in cape calculation.
! dec. 9, 1994
!-jjh          tpv(i,k) =tp(i,k)*(1.+zvirp1*qstp(i,k))/(1.+q(i,mx(i)))
!
            tpv(i,k) = (tp(i,k)+tpert(i))* (1._r8+zvirp1*qstp(i,k)) / (1._r8+q(i,mx(i)))
            buoy(i,k) = tpv(i,k) - tv(i,k) + tiedke_add
         end if
      end do
   end do
!
! main buoyancy calculation.
!
   do k = pver - 1,msg + 1,-1
      do i=1,ncol
         if (k < lcl(i) .and. plge600(i)) then
            tv(i,k) = t(i,k)* (1._r8+zvirp1*q(i,k))/ (1._r8+q(i,k))
            qstp(i,k) = qstp(i,k+1)
            tp(i,k) = tp(i,k+1)* (p(i,k)/p(i,k+1))**(0.2854_r8* (1._r8-0.28_r8*qstp(i,k)))
!          estp(i) = exp(a-b/tp(i,k))
           estp(i) = c1*exp((c2* (tp(i,k)-tfreez))/((tp(i,k)-tfreez)+c3))
           call esat_buck(tp(i,k), estp(i))

            qstp(i,k) = eps1*estp(i)/ (p(i,k)-estp(i))
            a1(i) = cp/rl + qstp(i,k)* (1._r8+qstp(i,k)/eps1)*rl*eps1/ (rd*tp(i,k)**2)
            a2(i) = .5_r8* (qstp(i,k)* (1._r8+2._r8/eps1*qstp(i,k))* &
                    (1._r8+qstp(i,k)/eps1)*eps1**2*rl*rl/ &
                    (rd**2*tp(i,k)**4)-qstp(i,k)* &
                    (1._r8+qstp(i,k)/eps1)*2._r8*eps1*rl/ &
                    (rd*tp(i,k)**3))
            a1(i) = 1._r8/a1(i)
            a2(i) = -a2(i)*a1(i)**3
            y(i) = qstp(i,k+1) - qstp(i,k)
            tp(i,k) = tp(i,k) + a1(i)*y(i) + a2(i)*y(i)**2
!          estp(i)  =exp(a-b/tp(i,k))
           estp(i) = c1*exp((c2* (tp(i,k)-tfreez))/ ((tp(i,k)-tfreez)+c3))
           call esat_buck(tp(i,k), estp(i))


            qstp(i,k) = eps1*estp(i)/ (p(i,k)-estp(i))
!-jjh          tpv(i,k) =tp(i,k)*(1.+zvirp1*qstp(i,k))/
!jt            (1.+q(i,mx(i)))
            tpv(i,k) = (tp(i,k)+tpert(i))* (1._r8+zvirp1*qstp(i,k))/(1._r8+q(i,mx(i)))
            buoy(i,k) = tpv(i,k) - tv(i,k) + tiedke_add
         end if
      end do
   end do

!
   do k = msg + 2,pver
      do i = 1,ncol
         if (k < lcl(i) .and. plge600(i)) then
            if (buoy(i,k+1) > 0._r8 .and. buoy(i,k) <= 0._r8) then
               knt(i) = min(5,knt(i) + 1)
               lelten(i,knt(i)) = k
            end if
         end if
      end do
   end do
!
! calculate convective available potential energy (cape).
!
   do n = 1,5
      do k = msg + 1,pver
         do i = 1,ncol
            if (plge600(i) .and. k <= mx(i) .and. k > lelten(i,n)) then
               capeten(i,n) = capeten(i,n) + rd*buoy(i,k)*log(pf(i,k+1)/pf(i,k))
            end if
         end do
      end do
   end do
!
! find maximum cape from all possible tentative capes from
! one sounding,
! and use it as the final cape, april 26, 1995
!
   do n = 1,5
      do i = 1,ncol
         if (capeten(i,n) > cape(i)) then
            cape(i) = capeten(i,n)
            lel(i) = lelten(i,n)
         end if
      end do
   end do
!
! put lower bound on cape for diagnostic purposes.
!
   do i = 1,ncol
      cape(i) = max(cape(i), 0._r8)
   end do
!
   return
end subroutine buoyan

subroutine cldprp(lchnk   , &
    q       ,t       ,u       ,v       ,p       , &
    z       ,s       ,mu      ,eu      ,du      , &
    md      ,ed      ,sd      ,qd      ,mc      , &
    qu      ,su      ,zf      ,qst     ,hmn     , &
    hsat    ,shat    ,qhat    ,ql      ,pf      , &
    cmeg    ,jb      ,lel     ,jt      ,jlcl    , &
    mx      ,j0      ,jd      ,rl      ,il2g    , &
    rd      ,grav    ,cp      ,msg     ,cape    , &
    pflx    ,evp     ,cu      ,rprd    ,limcnv  ,landfrac ,lcl_reached, lcl_est)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Actually calculates the properties of the updrafts and downdrafts 
! 
! Method: 
! may 09/91 - guang jun zhang, m.lazare, n.mcfarlane.
!             original version cldprop.
! 
! Author: See above, modified by P. Rasch
! This is contributed code not fully standardized by the CCM core group.
!
! this code is very much rougher than virtually anything else in the CCM
! there are debug statements left strewn about and code segments disabled
! these are to facilitate future development. We expect to release a
! cleaner code in a future release
!
! the documentation has been enhanced to the degree that we are able
!
!-----------------------------------------------------------------------
   use phys_grid, only: get_rlon_p, get_rlat_p

   implicit none

!------------------------------------------------------------------------------
!
! Input arguments
!
   integer, intent(in) :: lchnk                  ! chunk identifier

   real(r8), intent(in) :: q(pcols,pver)         ! spec. humidity of env
   real(r8), intent(in) :: t(pcols,pver)         ! temp of env
   real(r8), intent(in) :: p(pcols,pver)         ! pressure of env
   real(r8), intent(in) :: z(pcols,pver)         ! height of env
   real(r8), intent(in) :: s(pcols,pver)         ! normalized dry static energy of env
   real(r8), intent(in) :: zf(pcols,pverp)       ! height of interfaces
   real(r8), intent(in) :: pf(pcols,pverp)        !pressure height of interfaces
   real(r8), intent(in) :: u(pcols,pver)         ! zonal velocity of env
   real(r8), intent(in) :: v(pcols,pver)         ! merid. velocity of env

   real(r8), intent(in) :: landfrac(pcols)       ! RBN Landfrac
   real(r8), intent(in) :: cape(pcols)           ! CAPE as estimated in buoyan_dilute

   integer, intent(in) :: jb(pcols)              ! updraft base level
   integer, intent(in) :: lel(pcols)             ! updraft launch level
   integer, intent(out) :: jt(pcols)              ! updraft plume top
   integer, intent(out) :: jlcl(pcols)            ! updraft lifting cond level
   integer, intent(in) :: mx(pcols)              ! updraft base level (same is jb)
   integer, intent(out) :: j0(pcols)              ! level where updraft begins detraining
   integer, intent(out) :: jd(pcols)              ! level of downdraft
   integer, intent(in) :: limcnv                 ! convection limiting level
   integer, intent(in) :: il2g                   !CORE GROUP REMOVE
   integer, intent(in) :: msg                    ! missing moisture vals (always 0)
   real(r8), intent(in) :: rl                    ! latent heat of vap
   real(r8), intent(in) :: shat(pcols,pver)      ! interface values of dry stat energy
   real(r8), intent(in) :: qhat(pcols,pver)      ! interface values of q
!
! output
!
   real(r8), intent(out) :: rprd(pcols,pver)     ! rate of production of precip at that layer
   real(r8), intent(out) :: du(pcols,pver)       ! detrainement rate of updraft
   real(r8), intent(out) :: ed(pcols,pver)       ! entrainment rate of downdraft
   real(r8), intent(out) :: eu(pcols,pver)       ! entrainment rate of updraft
   real(r8), intent(out) :: hmn(pcols,pver)      ! moist stat energy of env
   real(r8), intent(out) :: hsat(pcols,pver)     ! sat moist stat energy of env
   real(r8), intent(out) :: mc(pcols,pver)       ! net mass flux
   real(r8), intent(out) :: md(pcols,pver)       ! downdraft mass flux
   real(r8), intent(out) :: mu(pcols,pver)       ! updraft mass flux
   real(r8), intent(out) :: pflx(pcols,pverp)    ! precipitation flux thru layer
   real(r8), intent(out) :: qd(pcols,pver)       ! spec humidity of downdraft
   real(r8), intent(out) :: ql(pcols,pver)       ! liq water of updraft
   real(r8), intent(out) :: qst(pcols,pver)      ! saturation spec humidity of env.
   real(r8), intent(out) :: qu(pcols,pver)       ! spec hum of updraft
   real(r8), intent(out) :: sd(pcols,pver)       ! normalized dry stat energy of downdraft
   real(r8), intent(out) :: su(pcols,pver)       ! normalized dry stat energy of updraft
   logical, intent(out) :: lcl_reached(pcols)    ! whether condensation has happened in a given column
   integer, intent(in) :: lcl_est(pcols)         ! estimate for lcl being reached from buoyan_dilute



   real(r8) rd                   ! gas constant for dry air
   real(r8) grav                 ! gravity
   real(r8) cp                   ! heat capacity of dry air

!
! Local workspace
!
   real(r8) gamma(pcols,pver)
   real(r8) dz(pcols,pver)
   real(r8) iprm(pcols,pver)
   real(r8) hu(pcols,pver)
   real(r8) hd(pcols,pver)
   real(r8) eps(pcols,pver)
   real(r8) f(pcols,pver)
   real(r8) k1(pcols,pver)
   real(r8) i2(pcols,pver)
   real(r8) ihat(pcols,pver)
   real(r8) i3(pcols,pver)
   real(r8) idag(pcols,pver)
   real(r8) i4(pcols,pver)
   real(r8) qsthat(pcols,pver)
   real(r8) hsthat(pcols,pver)
   real(r8) gamhat(pcols,pver)
   real(r8) cu(pcols,pver)
   real(r8) evp(pcols,pver)
   real(r8) cmeg(pcols,pver)
   real(r8) qds(pcols,pver)
   ! RBN For c0mask
   real(r8) c0mask(pcols)

   real(r8) hmin(pcols)
   real(r8) expdif(pcols)
   real(r8) expnum(pcols)
   real(r8) ftemp(pcols)
   real(r8) eps0(pcols)
   real(r8) rmue(pcols)
   real(r8) zuef(pcols)
   real(r8) zdef(pcols)
   real(r8) epsm(pcols)
   real(r8) ratmjb(pcols)
   real(r8) est(pcols)
   real(r8) totpcp(pcols)
   real(r8) totevp(pcols)
   real(r8) alfa(pcols)
   real(r8) ql1
   real(r8) estu
   real(r8) qstu

   real(r8) qv_in
   real(r8) ql_in
   real(r8) t_det


   real(r8) eta_min(pcols)                !entropy minimum in a single column
   real(r8) h_detrain(pcols,pver)         !moist static energy at a given height of a saturated parcel with the same buoyancy as the environment
   real(r8) s_detrain(pcols,pver)         !normalised dry static energy at a given height of a saturated parcel with the same buoyancy as the environment
   real(r8) q_detrain(pcols,pver)         !water vapour mixing ratio at a given height of a saturated parcel with the same buoyancy as the environment
   real(r8) eta_detrain(pcols,pver)       !entropy at a given height of a saturated parcel with the same buoyancy as the environment
   real(r8) eta_t_detrain(pcols,pver)     !entropy at a given height of a saturated parcel with the same buoyancy as the environment, with contribution for entropy of precipitation lost
   real(r8) eta(pcols,pver)               !entropy of environment at each height
   real(r8) eta_u(pcols,pver)             !entropy of updraft
   real(r8) eta_d(pcols,pver)             !entropy of downdraft
   real(r8) cpmix                         !holder for composition-weighted heat capacity
   real(r8) dqv                           !holder for change in moisture of downdraft as a result of saturation at every level
   real(r8) mass_scaling                  !holder for change in mass of plume due to precipitation
   real(r8) qsat                          !holder for saturation q when detraining
   real(r8) qtot                          !total q in a saturated rising plume
   real(r8) qv,qliq, qv_out, qd1
   real(r8) tu(pcols,pver)                !temperature of updraft (defined at interfaces)
   real(r8) that(pcols,pver)              !temperature of env at interfaces
   real(r8) tv(pcols,pver)                !environmental virtual temperature at interfaces
   real(r8) tpv(pcols,pver)               !updraft parcel virtual temperature (defined at interfaces)
   real(r8) tfg
   real(r8) tp                            !parcel temperature returned at various points
   real(r8) t_out                         !downdraft temperature returned at various points
   real(r8) mbm(pcols, pver)              !effective starting mass of updraft, adjusted for precipitation
   real(r8) mtm(pcols, pver)              !effective starting mass of downdraft, adjusted for evaporation
   real(r8) q_max(pcols)                  !maximum env q encountered in course of plume
   real(r8) q_min(pcols)                  !minimum env q encountered in course of plume
   real(r8) s_max(pcols)                  !maximum env s encountered in course of plume
   real(r8) s_min(pcols)                  !minimum env s encountered in course of plume   
   real(r8) ent_max(pcols)                !maximum env eta encountered in course of plume
   real(r8) ent_min(pcols)                !minimum env eta encountered in course of plume

   real(r8) etag(pcols,pver)
   real(r8) etag_detrain(pcols,pver)
   real(r8) etag_t_detrain(pcols,pver)
   real(r8) iprm_g(pcols,pver)
   real(r8) k1_g(pcols,pver)
   real(r8) i2_g(pcols,pver)
   real(r8) ihat_g(pcols,pver)
   real(r8) i3_g(pcols,pver)
   real(r8) idag_g(pcols,pver)
   real(r8) i4_g(pcols,pver)
   real(r8) expnum_g(pcols)
   real(r8) etag_min(pcols)
   integer j0_g(pcols)


   real(r8) small
   real(r8) mdt

   integer khighest
   integer klowest
   integer kount
   integer i,k
   integer rcall

   logical doit(pcols)
   logical done(pcols)
   logical plume_top_reached(pcols)
   logical dry_downdraft(pcols)
   real(r8) max_entrainment(pcols)

   real(r8) this_lat, this_lon
!
!------------------------------------------------------------------------------
!

   !write(iulog,*) "***ZM_CONV: Finding cloud properties. il2g: ", il2g

   do i = 1,il2g
      ftemp(i) = 0._r8
      expnum(i) = 0._r8
      expdif(i) = 0._r8
      c0mask(i)  = c0_ocn * (1._r8-landfrac(i)) +   c0_lnd * landfrac(i)
      max_entrainment(i) = 7.e-6_r8 * shr_const_mwdair/(1-mu_red*q(i,mx(i)))      ! Maximum entrainment rate, adjusted so = 2e-4 / m for Earth atmosphere
   end do
!
!jr Change from msg+1 to 1 to prevent blowup
!
   do k = 1,pver
      do i = 1,il2g
         dz(i,k) = zf(i,k) - zf(i,k+1)
      end do
   end do

!
! initialize many output and work variables to zero
!
   pflx(:il2g,1) = 0

   do k = 1,pver
      do i = 1,il2g
         k1(i,k) = 0._r8
         i2(i,k) = 0._r8
         i3(i,k) = 0._r8
         i4(i,k) = 0._r8
         mu(i,k) = 0._r8
         f(i,k) = 0._r8
         eps(i,k) = 0._r8
         eu(i,k) = 0._r8
         du(i,k) = 0._r8
         ql(i,k) = 0._r8
         cu(i,k) = 0._r8
         evp(i,k) = 0._r8
         cmeg(i,k) = 0._r8
         qds(i,k) = q(i,k)
         md(i,k) = 0._r8
         ed(i,k) = 0._r8
         sd(i,k) = s(i,k)
         qd(i,k) = q(i,k)
         mc(i,k) = 0._r8
         qu(i,k) = q(i,k)
         su(i,k) = s(i,k)
         tu(i,k) = t(i,k)
!        est(i)=exp(a-b/t(i,k))
         est(i) = c1*exp((c2* (t(i,k)-tfreez))/((t(i,k)-tfreez)+c3))
         call esat_buck(t(i,k), est(i))

         !++bee
         if ( p(i,k)-est(i) > 0._r8 ) then
            qst(i,k) = eps1*est(i)/ (p(i,k) - est(i)*(1-eps1))
         else
            qst(i,k) = 1.0_r8
         end if
         !--bee
         gamma(i,k) = qst(i,k)*(1._r8 + qst(i,k)/eps1)*eps1*rl/(rd*t(i,k)**2)*rl/cp
         hmn(i,k) = ((1-q(i,k))*cpres + q(i,k)*cpliq)*t(i,k) + grav*z(i,k) + rl*q(i,k)  !changed to account for different env cp
         hsat(i,k) = ((1-q(i,k))*cpres + q(i,k)*cpliq)*t(i,k) + grav*z(i,k) + rl*qst(i,k)  !changed to account for different env cp
         hu(i,k) = hmn(i,k)
         hd(i,k) = hmn(i,k)
         rprd(i,k) = 0._r8

         rcall=2
         Tfg=t(i,k)
         call itv(rcall,i,lchnk,p(i,k),tp,q(i,k),t(i,k),qsat,q(i,mx(i)),Tfg, k, lcl_est(i)) !assume no significant cloud liquid in either parcel or env
         h_detrain(i,k) = ((1-qsat)*cpres + qsat*cpliq)*tp + grav*z(i,k) + rl*qsat
         s_detrain(i,k) = (((1-qsat)*cpres + qsat*cpwv)*tp + grav*z(i,k))/cpres
         q_detrain(i,k) = qsat
         call entropy(tp,p(i,k),qsat, eta_detrain(i,k))
         eta_t_detrain(i,k) = eta_detrain(i,k) + (q(i,mx(i))-qsat)*(cpliq*log(tp/tfreez)-eta_detrain(i,k))  !rough measure of entropy contribution from existing precipitation
         call entropy(t(i,k),p(i,k),q(i,k), eta(i,k))

         ! if (i==1) then
         !    write(iulog,*) "i:", i, "k:", k, "s:", s(i,k)
         !    write(iulog,*) "eta:", eta(i,k), "eta_det:", eta_detrain(i,k), "eta_t_det:", eta_t_detrain(i,k)
         !    write(iulog,*) "qmx:", q(i,mx(i)), "q:", q(i,k), "q_det:", qsat, "t:", t(i,k), "t_det:", tp
         ! end if
         
         eta_u(i,k) = eta(i,k)
         eta_d(i,k) = eta(i,k)

         ! if (i==1) then
         !    write(iulog,*) " "
         !   write(iulog,*) "*** CLDPRP: detraining plume properties. local P:", p(i,k), "T: ", t(i,k), "Q: ", q(i,k), "k:", k
         !   write(iulog,*) "*** CLDPRP. T_detrain: ", tp, "Q_detrain: ", qsat, "Local h: ", hmn(i,k), "h_detrain: ", h_detrain(i,k)
         !   write(iulog,*) "*** CLDPRP. s_detrain: ", s_detrain(i,k), "local eta: ", eta(i,k), ", eta_detrain: ", eta_detrain(i,k), "modified eta_t_detrain: ", eta_t_detrain(i,k)
         ! end if

         if ((lcl_est(i) < mx(i) - 10) .and. (k > lcl_est(i) + 3)) then  !we are in a regime where we want to pick the gas entropy
            call entropy_gas(t(i,k),p(i,k),q(i,k),etag(i,k))
            call entropy_gas(tp,p(i,k), qsat, etag_detrain(i,k))
            etag_t_detrain(i,k) = etag_detrain(i,k) ! as no precipitation
         end if

      end do
   end do

!
!jr Set to zero things which make this routine blow up
!
   do k=1,msg
      do i=1,il2g
         rprd(i,k) = 0._r8
      end do
   end do
!
! interpolate the layer values of qst, hsat and gamma to
! layer interfaces, define and calculate that and consequently tv, tpv
!
   do i = 1,il2g
      hsthat(i,msg+1) = hsat(i,msg+1)
      qsthat(i,msg+1) = qst(i,msg+1)
      gamhat(i,msg+1) = gamma(i,msg+1)
      totpcp(i) = 0._r8
      totevp(i) = 0._r8
      dry_downdraft(i) = .false.
   end do
   do k = msg + 2,pver
      do i = 1,il2g
         if (abs(qst(i,k-1)-qst(i,k)) > 1.E-6_r8) then
            qsthat(i,k) = log(qst(i,k-1)/qst(i,k))*qst(i,k-1)*qst(i,k)/ (qst(i,k-1)-qst(i,k))
         else
            qsthat(i,k) = qst(i,k)
         end if
         hsthat(i,k) = cp*shat(i,k) + rl*qsthat(i,k)
         if (abs(gamma(i,k-1)-gamma(i,k)) > 1.E-6_r8) then
            gamhat(i,k) = log(gamma(i,k-1)/gamma(i,k))*gamma(i,k-1)*gamma(i,k)/ &
                           (gamma(i,k-1)-gamma(i,k))
         else
            gamhat(i,k) = gamma(i,k)
         end if
         cpmix = (1 - qhat(i,k))*cpres + qhat(i,k)*cpwv
         that(i,k) = (cpres*shat(i,k) - grav*zf(i,k))/cpmix  !env T at interface
         tv(i,k) = that(i,k) * (1 - mu_red*qhat(i,k))
         tpv(i,k) = tv(i,k)
      end do
   end do
!
! initialize cloud top to highest plume top.
!jr changed hard-wired 4 to limcnv+1 (not to exceed pver)
!
   jt(:) = pver
   do i = 1,il2g
      jt(i) = max(lel(i),limcnv+1)
      jt(i) = min(jt(i),pver)
      jd(i) = pver
      jlcl(i) = lel(i)
      eta_min(i) = 1.E6_r8
   end do
!
! find the level of minimum entropy, where detrainment starts
!

   do k = msg + 1,pver
      do i = 1,il2g
         if (eta(i,k) <= eta_min(i) .and. k >= jt(i) .and. k <= jb(i)) then
            eta_min(i) = eta(i,k)
            j0(i) = k
         end if

         if (etag(i,k) <= etag_min(i) .and. k >= jt(i) .and. k <= jb(i)) then
            etag_min(i) = etag(i,k)
            j0_g(i) = k
         end if        
      end do
   end do
   do i = 1,il2g
      j0(i) = min(j0(i),jb(i)-2)
      j0(i) = max(j0(i),jt(i)+2)
      !
      ! Fix from Guang Zhang to address out of bounds array reference
      !
      j0(i) = min(j0(i),pver)

      j0_g(i) = min(j0_g(i),jb(i)-2)
      j0_g(i) = max(j0_g(i),jt(i)+2)
      j0_g(i) = min(j0_g(i),pver)

      !now combine the two entropy minima within the bounds of relevance to find the true entropy minimum.
      !if there are two minima, pick the higher one (should be closer to actual cloud deck)
      if (((lcl_est(i) < mx(i) - 10) .and. (j0_g(i) > lcl_est(i) + 3)) .and. (j0(i) < lcl_est(i) + 3)) then
         j0(i) = min(j0(i), j0_g(i))
      else if ((lcl_est(i) < mx(i) - 10) .and. (j0_g(i) > lcl_est(i) + 3)) then  !only gas phase entropy interests us
         j0(i) = j0_g(i)
      end if !in remaining moist only case we leave j0 as is
   end do
!
! Initialize certain arrays inside cloud
!
   do k = msg + 1,pver
      do i = 1,il2g
         if (k >= jt(i) .and. k <= jb(i)) then
            hu(i,k) = hmn(i,mx(i)) + cp*tiedke_add   !could change to cpmix*tiedke_add, not going to be very impactful
            su(i,k) = s(i,mx(i)) + tiedke_add
         end if
      end do
   end do
!
! *********************************************************
! compute taylor series for approximate eps(z) below
! *********************************************************
!
   do k = pver - 1,msg + 1,-1
      do i = 1,il2g
         if (k < jb(i) .and. k >= jt(i)) then
            k1(i,k) = k1(i,k+1) + (eta(i,mx(i))-eta(i,k))*dz(i,k)
            ihat(i,k) = 0.5_r8* (k1(i,k+1)+k1(i,k))
            i2(i,k) = i2(i,k+1) + ihat(i,k)*dz(i,k)
            idag(i,k) = 0.5_r8* (i2(i,k+1)+i2(i,k))
            i3(i,k) = i3(i,k+1) + idag(i,k)*dz(i,k)
            iprm(i,k) = 0.5_r8* (i3(i,k+1)+i3(i,k))
            i4(i,k) = i4(i,k+1) + iprm(i,k)*dz(i,k)

            ! find coefficients for gas entropy expression, these may be used later
            k1_g(i,k) = k1_g(i,k+1) + (etag(i,mx(i))-etag(i,k))*dz(i,k)
            ihat_g(i,k) = 0.5_r8* (k1_g(i,k+1)+k1_g(i,k))
            i2_g(i,k) = i2_g(i,k+1) + ihat_g(i,k)*dz(i,k)
            idag_g(i,k) = 0.5_r8* (i2_g(i,k+1)+i2_g(i,k))
            i3_g(i,k) = i3_g(i,k+1) + idag_g(i,k)*dz(i,k)
            iprm_g(i,k) = 0.5_r8* (i3_g(i,k+1)+i3_g(i,k))
            i4_g(i,k) = i4_g(i,k+1) + iprm_g(i,k)*dz(i,k)

         end if
      end do
   end do
!
! re-initialize hmin array for ensuing calculation. Choosing not to edit 
! this as this is only an approximate condition for overall convection already and should remain valid
! actually - h not decreasing with height does not necessarily mean no convection, should probs change this
!
   do i = 1,il2g
      hmin(i) = 1.E9_r8
   end do
   do k = msg + 1,pver
      do i = 1,il2g
         if (k >= j0(i) .and. k <= jb(i) .and. hmn(i,k) <= hmin(i)) then
            hmin(i) = hmn(i,k)
            expdif(i) = hmn(i,mx(i)) - hmin(i)
         end if
      end do
   end do
!
! *********************************************************
! compute approximate eps(z) using above taylor series
! *********************************************************
!
   do k = msg + 2,pver
      do i = 1,il2g
         expnum(i) = 0._r8
         ftemp(i) = 0._r8
         expnum_g(i) = 0._r8

         if (k < jt(i) .or. k >= jb(i)) then
            k1(i,k) = 0._r8
            expnum(i) = 0._r8

            k1_g(i,k) = 0._r8  !
            expnum_g(i) = 0._r8
         else
            expnum(i) = (eta(i,mx(i)) - (eta_t_detrain(i,k-1)*(zf(i,k)-z(i,k)) + &
                     eta_t_detrain(i,k)* (z(i,k-1)-zf(i,k)))/(z(i,k-1)-z(i,k))) !/ 2._r8 !can divide by 2 here

            expnum_g(i) = (etag(i,mx(i)) - (etag_t_detrain(i,k-1)*(zf(i,k)-z(i,k)) + &
                     etag_t_detrain(i,k)* (z(i,k-1)-zf(i,k)))/(z(i,k-1)-z(i,k))) !/ 2._r8 !can divide by 2 here
         end if

         if ((expdif(i) > 100._r8 .and. expnum(i) > 0._r8) .and. &
            k1(i,k) > expnum(i)*dz(i,k)) then
            ftemp(i) = expnum(i)/k1(i,k)
            f(i,k) = ftemp(i) + i2(i,k)/k1(i,k)*ftemp(i)**2 + &
                  (2._r8*i2(i,k)**2-k1(i,k)*i3(i,k))/k1(i,k)**2* &
                  ftemp(i)**3 + (-5._r8*k1(i,k)*i2(i,k)*i3(i,k)+ &
                  5._r8*i2(i,k)**3+k1(i,k)**2*i4(i,k))/ &
                  k1(i,k)**3*ftemp(i)**4
            ! if (i==1) then
            !    write(iulog,*) " "
            !    write(iulog,*) "lchnk:", lchnk, "i:", i, "k:", k, "s:", s(i,k)
            !    write(iulog,*) "expnum:", expnum(i), "ftemp:", ftemp(i), "f:", f(i,k), "max entrainment:", max_entrainment(i)
            !    write(iulog,*) "base eta:", eta(i,mx(i)), "level eta:", eta(i,k), "detraining eta_t:", eta_t_detrain(i,k)
            !    write(iulog,*) "k1:", k1(i,k), "i2:", i2(i,k), "i3:", i3(i,k), "i4:", i4(i,k)
            !    write(iulog,*) "2nd coeff:", i2(i,k)/k1(i,k), "3rd coeff:", (2._r8*i2(i,k)**2-k1(i,k)*i3(i,k))/k1(i,k)**2, & 
            !    "4th coeff:", (-5._r8*k1(i,k)*i2(i,k)*i3(i,k)+ 5._r8*i2(i,k)**3+k1(i,k)**2*i4(i,k))/k1(i,k)**3
            !    write(iulog,*) "1st cont:", ftemp(i), "2nd cont:", i2(i,k)/k1(i,k)*ftemp(i)**2, "3rd cont:", (2._r8*i2(i,k)**2-k1(i,k)*i3(i,k))/k1(i,k)**2 *ftemp(i)**3, & 
            !    "4th cont:", (-5._r8*k1(i,k)*i2(i,k)*i3(i,k)+ 5._r8*i2(i,k)**3+k1(i,k)**2*i4(i,k))/k1(i,k)**3*ftemp(i)**4
            ! end if
            f(i,k) = max(f(i,k),0._r8)
            f(i,k) = min(f(i,k),max_entrainment(i))  !can divide by 2 here
         end if
         ! if (i==1) then
         !    write(iulog,*) " "
         !    write(iulog,*) "k:", k, "expdif:", expdif(i), "expnum:", expnum(i), "ftemp:", ftemp(i), "f:", f(i,k)
         !    write(iulog,*) "base eta:", eta(i,mx(i)), "level eta:", eta(i,k), "detraining eta_t:", eta_t_detrain(i,k)
         ! end if

         if ((lcl_est(i) < mx(i) - 10) .and. (k > lcl_est(i) + 3)) then  !we are in a regime where we want to pick the gas entropy
            ! if we are in the regime to use gas entropy then provide a new lambda
            if ((expdif(i) > 100._r8 .and. expnum_g(i) > 0._r8) .and. &
            k1_g(i,k) > expnum_g(i)*dz(i,k)) then
               ftemp(i) = expnum_g(i)/k1_g(i,k)
               f(i,k) = ftemp(i) + i2_g(i,k)/k1_g(i,k)*ftemp(i)**2 + &
                     (2._r8*i2_g(i,k)**2-k1_g(i,k)*i3_g(i,k))/k1_g(i,k)**2* &
                     ftemp(i)**3 + (-5._r8*k1_g(i,k)*i2_g(i,k)*i3_g(i,k)+ &
                     5._r8*i2_g(i,k)**3+k1_g(i,k)**2*i4_g(i,k))/ &
                     k1_g(i,k)**3*ftemp(i)**4
                  ! if (i==1) then
                  !    write(iulog,*) " "
                  !    write(iulog,*) "GAS ENTROPY: lchnk:", lchnk, "i:", i, "k:", k, "s:", s(i,k)
                  !    write(iulog,*) "expnum:", expnum_g(i), "ftemp:", ftemp(i), "f:", f(i,k), "max entrainment:", max_entrainment(i)
                  !    write(iulog,*) "base eta:", etag(i,mx(i)), "level eta:", etag(i,k), "detraining eta_t:", etag_t_detrain(i,k)
                  !    write(iulog,*) "k1:", k1_g(i,k), "i2:", i2_g(i,k), "i3:", i3_g(i,k), "i4:", i4_g(i,k)
                  !    write(iulog,*) "2nd coeff:", i2_g(i,k)/k1_g(i,k), "3rd coeff:", (2._r8*i2_g(i,k)**2-k1_g(i,k)*i3_g(i,k))/k1_g(i,k)**2, & 
                  !    "4th coeff:", (-5._r8*k1_g(i,k)*i2_g(i,k)*i3_g(i,k)+ 5._r8*i2_g(i,k)**3+k1_g(i,k)**2*i4_g(i,k))/k1_g(i,k)**3
                  !    write(iulog,*) "1st cont:", ftemp(i), "2nd cont:", i2_g(i,k)/k1_g(i,k)*ftemp(i)**2, "3rd cont:", (2._r8*i2_g(i,k)**2-k1_g(i,k)*i3_g(i,k))/k1_g(i,k)**2 *ftemp(i)**3, & 
                  !    "4th cont:", (-5._r8*k1_g(i,k)*i2_g(i,k)*i3_g(i,k)+ 5._r8*i2_g(i,k)**3+k1_g(i,k)**2*i4_g(i,k))/k1_g(i,k)**3*ftemp(i)**4
                  ! end if
               f(i,k) = max(f(i,k),0._r8)
               f(i,k) = min(f(i,k),max_entrainment(i))  !can divide by 2 here
            end if
         end if
      end do
   end do
   do i = 1,il2g
      if (j0(i) < jb(i)) then
         if (f(i,j0(i)) < 1.E-6_r8 .and. f(i,j0(i)+1) > f(i,j0(i))) j0(i) = j0(i) + 1
      end if
   end do
   do k = msg + 2,pver
      do i = 1,il2g
         if (k >= jt(i) .and. k <= j0(i)) then
            f(i,k) = max(f(i,k),f(i,k-1))
         end if
      end do
   end do
   do i = 1,il2g
      eps0(i) = f(i,j0(i))
      eps(i,jb(i)) = eps0(i)
   end do
!
! This is set to match the Rasch and Kristjansson paper
!
   do k = pver,msg + 1,-1
      do i = 1,il2g
         if (k >= j0(i) .and. k <= jb(i)) then
         eps(i,k) = f(i,j0(i))
         end if
      end do
   end do
   do k = pver,msg + 1,-1
      do i = 1,il2g
         if (k < j0(i) .and. k >= jt(i)) eps(i,k) = f(i,k)
      end do
   end do

   do i = 1,il2g
      if (cape(i) < capelmt) then
         eps0(i) = 0._r8 !convection disabled, CAPE not high enough to trigger it
      end if
   end do


!do i = 1,il2g
!   this_lat = get_rlat_p(lchnk,i)*57.296_r8
!   this_lon = get_rlon_p(lchnk,i)*57.296_r8
   !write(iulog,*) "*** ZM_CONV: Buoyan_dilute. lat: ", this_lat, "lon: ", this_lon
   !if ((this_lat < -45._r8 .and. this_lat > -55._r8) .and. (this_lon < 15._r8 .and. this_lon > 5._r8)) then !region around substellar point
!   if (cape(i) > 2500._r8) then
!      write(iulog,*) " "
!      write(iulog,*) " "
!      write(iulog,*) "*** ZM_CONV: CLDPRP. At lchnk: ", lchnk, "i: ", i, " lat: ", this_lat, "lon: ", this_lon
!      write(iulog,*) "*** ZM_CONV: CLDPRP. jb: ", jb(i), "j0: ", j0(i), "jt: ", jt(i)
!      write(iulog,*) " "
!      do k = pver,1,-1
!         write(iulog,*) "*** Going up column.  T: ", t(i,k), "P: ", p(i,k), "Q: ", q(i,k), "Z: ", z(i,k)
!         write(iulog,*) "*** Going up column. eta: ", eta(i,k), "eta_detrain: ", eta_detrain(i,k), "eta_t_detrain: ", eta_t_detrain(i,k)
!         write(iulog,*) "*** Going up column. f: ", f(i,k), "eps: ", eps(i,k)
!      end do
!   end if
!end do

!
   do i = 1,il2g
      if (eps0(i) > 0._r8) then
         lcl_reached(i) = .false. 
         plume_top_reached(i) = .false.
      end if
   end do

!
!NOW START THE MAIN UPDRAFT LOOP
!here and below mu, eu,du, md and ed are all normalized by mb
!consider khighest/klowest?
   !not considering changing cloud level if only one level thick or similar things
!

   do k = pver, msg+1, -1
      do i = 1,il2g
         if (eps0(i) > 0._r8) then
            if (k == jb(i)) then
               !set mass fluxes, and plume values
               mbm(i,k) = 1
               mu(i,k) = 1
               eu(i,k) = mu(i,k)/dz(i,k)
               hu(i,k) = hmn(i,k) + cp*tiedke_add
               su(i,k) = s(i,k) + tiedke_add
               qu(i,k) = q(i,k)
               call entropy(t(i,k)+tiedke_add, pf(i,k),q(i,k), eta_u(i,k))
               tu(i,k) = t(i,k)
               tpv(i,k) = tu(i,k) * (1 - mu_red*qu(i,k))
               !check to see if plume is already saturated (as we would expect)
               estu = c1*exp((c2* (tu(i,k)+tiedke_add -tfreez))/ ((tu(i,k)+tiedke_add -tfreez)+c3))
               call esat_buck(tu(i,k)+tiedke_add, estu)

               qstu = eps1*estu/ ((p(i,k)+p(i,k-1))/2._r8 - estu*(1-eps1))
               if (qu(i,k) >= qstu) then
                  jlcl(i) = k
                  lcl_reached(i) = .true.
               end if

               q_max(i) = q(i,k)
               q_min(i) = q(i,k)
               s_max(i) = s(i,k)
               s_min(i) = s(i,k)
               ent_max(i) = eta(i,k)
               ent_min(i) = eta(i,k)
               ! if (qu(i,k)>0.3_r8) then
               !    write(iulog,*) " "
               !    write(iulog,*) " "
               !    write(iulog,*) "lchnk:", lchnk, "i:", i, "k:", K
               !    write(iulog,*) "q:", q(i,k), "qu:", qu(i,k), "qhat:", qhat(i,k)
               !    write(iulog,*) "s:", s(i,k), "su:", su(i,k), "shat:", shat(i,k), "tiedke_add:", tiedke_add
               !    write(iulog,*) "pf:", pf(i,k), "eta_u:", eta_u(i,k)
               !    write(iulog,*) " "
               !    write(iulog,*) " "
               ! end if
            
            else if ((k < jb(i) .and. k > lel(i)) .and. (.not. lcl_reached(i) .and. .not. plume_top_reached(i))) then !starting from the bottom, no precipitation
               !masses calculated without any mass loss, and assuming that the lcl is reached higher up
               zuef(i) = zf(i,k) - zf(i,jb(i))
               rmue(i) = (mbm(i,k+1)/eps0(i))* (exp(eps(i,k+1)*zuef(i))-1._r8)/zuef(i)
               mu(i,k) = (mbm(i,k+1)/eps0(i))* (exp(eps(i,k  )*zuef(i))-1._r8)/zuef(i)
               eu(i,k) = (rmue(i)-mu(i,k+1))/dz(i,k)
               du(i,k) = (rmue(i)-mu(i,k))/dz(i,k)
               mbm(i,k) = mbm(i,k+1) !because no mass loss through precipitation

               !update max/min q, s, eta
               if (q_max(i)<q(i,k)) q_max(i)=q(i,k)
               if (q_min(i)>q(i,k)) q_min(i)=q(i,k)
               if (s_max(i)<s(i,k)) s_max(i)=s(i,k)
               if (s_min(i)>s(i,k)) s_min(i)=s(i,k)
               if (ent_max(i)<eta(i,k)) ent_max(i)=eta(i,k)
               if (ent_min(i)>eta(i,k)) ent_min(i)=eta(i,k)
               !calculate plume properties
               s_detrain(i,k) = min(s_detrain(i,k), su(i,k+1)) !warmer plumes are not the ones detraining
               h_detrain(i,k) = min(h_detrain(i,k), hu(i,k+1)) !warmer plumes are not the ones detraining

               hu(i,k) = mu(i,k+1)/mu(i,k)*hu(i,k+1) + &
                         dz(i,k)/mu(i,k)* (eu(i,k)*hmn(i,k) - du(i,k)*h_detrain(i,k))
               su(i,k) = mu(i,k+1)/mu(i,k)*su(i,k+1) + &
                         dz(i,k)/mu(i,k)* (eu(i,k)*s(i,k)   - du(i,k)*s_detrain(i,k)) 
               qu(i,k) = mu(i,k+1)/mu(i,k)*qu(i,k+1) + &
                         dz(i,k)/mu(i,k)* (eu(i,k)*q(i,k)   - du(i,k)*q_detrain(i,k))

               !impose checks on s, q in updraft so they stay at sensible values
               if (su(i,k) > s_max(i)) then !we are capping su at highest s encountered along the way
                  s_detrain(i,k) = eu(i,k)/du(i,k)*s(i,k) + mu(i,k+1)/(dz(i,k)*du(i,k))*su(i,k+1) - &
                         mu(i,k)/(dz(i,k)*du(i,k))*s_max(i)
                  su(i,k) = s_max(i)
               end if
               if (su(i,k) < s_min(i)) then !we are capping su at lowest s encountered along the way (unlikely, but you never know)
                  s_detrain(i,k) = eu(i,k)/du(i,k)*s(i,k) + mu(i,k+1)/(dz(i,k)*du(i,k))*su(i,k+1) - &
                         mu(i,k)/(dz(i,k)*du(i,k))*s_min(i)
                  su(i,k) = s_min(i)
               end if
               if (qu(i,k) > q_max(i)) then
                  q_detrain(i,k) = (eu(i,k)/du(i,k))*q(i,k) + mu(i,k+1)/(dz(i,k)*du(i,k))*qu(i,k+1) - &
                         mu(i,k)/(dz(i,k)*du(i,k))*q_max(i)
                  qu(i,k) = q_max(i)
               end if
               if (qu(i,k) < q_min(i)) then
                  q_detrain(i,k) = (eu(i,k)/du(i,k))*q(i,k) + mu(i,k+1)/(dz(i,k)*du(i,k))*qu(i,k+1) - &
                         mu(i,k)/(dz(i,k)*du(i,k))*q_min(i)
                  qu(i,k) = q_min(i)
               end if


               cpmix = (1 - q_detrain(i,k))*cpres + q_detrain(i,k)*cpwv
               t_det = (cpres/cpmix)*s_detrain(i,k) - grav*z(i,k)/cpmix !calculate temperature of detraining material to find entropy
               call entropy(t_det, p(i,k),q_detrain(i,k), eta_detrain(i,k)) !updraft entropy at layer top

               cpmix = (1 - qu(i,k))*cpres + qu(i,k)*cpwv
               tu(i,k) = (cpres/cpmix)*su(i,k) - grav*zf(i,k)/cpmix !need to find tu at layer top to calculate if LCL reached at some point in layer
               call entropy(tu(i,k), pf(i,k),qu(i,k), eta_u(i,k)) !updraft entropy at layer top

             
               !check for saturation
               estu = c1*exp((c2* (tu(i,k)-tfreez))/ ((tu(i,k)-tfreez)+c3))
               call esat_buck(tu(i,k), estu)

               qstu = eps1*estu/ (pf(i,k) - estu*(1-eps1)) !need to change from estu to e

               ! if (qu(i,k)<0._r8) then
               !    write(iulog,*) " "
               !    this_lat = get_rlat_p(lchnk,i)*57.296_r8
               !    this_lon = get_rlon_p(lchnk,i)*57.296_r8
               !    write(iulog,*) "*** ZM_CONV: CLDPRP. At lchnk:", lchnk, "i:", i, " lat:", this_lat, "lon:", this_lon, "k:", k
               !    write(iulog,*) "*** ZM_CONV: CLDPRP. jb: ", jb(i), "j0: ", j0(i), "jt: ", jt(i)
               !    write(iulog,*) "This level: T:", t(i,k), "P:", p(i,k), "Q:", q(i,k), "Z:", z(i,k), "h:", hmn(i,k), "s:", s(i,k), "eta:", eta(i,k)
               !    write(iulog,*) "This level: h_detrain:", h_detrain(i,k), "s_detrain:", s_detrain(i,k), "q_detrain:", q_detrain(i,k)
               !    write(iulog,*) "This level: mu:", mu(i,k), "eu:", eu(i,k), "du:", du(i,k), "dz:",dz(i,k)
               !    write(iulog,*) "Plume mass: mbm:", mbm(i,k), "eps0:", eps0(i), "eps:", eps(i,k), "zuef:", zuef(i), "rmue:", rmue(i)
               !    write(iulog,*) "This level: hu:", hu(i,k), "su:", su(i,k), "qu:", qu(i,k), "eta_u:", eta_u(i,k), "tu:", tu(i,k)
               !    write(iulog,*) "This level: cpmix: ", cpmix, "qstu:", qstu, "eps1:", eps1
               !    write(iulog,*) "Level below: mu:", mu(i,k+1), "hu", hu(i,k+1), "su:", su(i,k+1), "qu:", qu(i,k+1), "eta_u:", eta_u(i,k+1), "q:", q(i,k)
               !    write(iulog,*) " "
               ! end if
               
               if (qu(i,k) >= qstu) then

                  ! if (qu(i,k)>1._r8) then
                  !    write(iulog,*) " "
                  !    this_lat = get_rlat_p(lchnk,i)*57.296_r8
                  !    this_lon = get_rlon_p(lchnk,i)*57.296_r8
                  !    write(iulog,*) "*** ZM_CONV: CLDPRP. At lchnk: ", lchnk, "i: ", i, " lat: ", this_lat, "lon: ", this_lon, "k:", k
                  !    write(iulog,*) "*** This level: T:", t(i,k), "P:", p(i,k), "Pf:", pf(i,k), "Q:", q(i,k), "Z:", z(i,k), "Zf:", zf(i,k)
                  !    write(iulog,*) "*** This level: That:", that(i,k), "qhat:", qhat(i,k), "h:", hmn(i,k), "s:", s(i,k), "eta:", eta(i,k)
                  !    write(iulog,*) "*** This level: h_det:", h_detrain(i,k), "s_det:", s_detrain(i,k), "q_det:", q_detrain(i,k), "eta_det", eta_detrain(i,k)
                  !    write(iulog,*) "*** This level: mu:", mu(i,k), "eu:", eu(i,k), "du:", du(i,k), "dz:", dz(i,k)
                  !    write(iulog,*) "Level below: mu:", mu(i,k+1), "hu", hu(i,k+1), "su:", su(i,k+1), "qu:", qu(i,k+1), "eta_u:", eta_u(i,k+1)
                  !    write(iulog,*) "*** Dry convection estimates: hu:", hu(i,k), "su:", su(i,k), "qu:", qu(i,k), "eta_u:", eta_u(i,k), "tu:", tu(i,k), "qstu:", qstu, "estu:", estu
                  ! end if
                  

                  jlcl(i) = k
                  lcl_reached(i) = .true.
                  !have to do layer again with correct supersaturation expressions
                  !now need to invert entropy
                  eta_u(i,k) = mu(i,k+1)/mu(i,k)*eta_u(i,k+1) + &
                  dz(i,k)/mu(i,k)* (eu(i,k)*eta(i,k) - du(i,k)*eta_detrain(i,k))

                  !check on eta in updraft to make sure it stays at sensible values
                  if (eta_u(i,k) > ent_max(i)) then
                     eta_detrain(i,k) = (eu(i,k)/du(i,k))*eta(i,k) + mu(i,k+1)/(dz(i,k)*du(i,k))*eta_u(i,k+1) - &
                            mu(i,k)/(dz(i,k)*du(i,k))*ent_max(i)
                     eta_u(i,k) = ent_max(i)
                  end if
                  if (eta_u(i,k) < ent_min(i)) then
                     eta_detrain(i,k) = (eu(i,k)/du(i,k))*eta(i,k) + mu(i,k+1)/(dz(i,k)*du(i,k))*eta_u(i,k+1) - &
                            mu(i,k)/(dz(i,k)*du(i,k))*ent_min(i)
                     eta_u(i,k) = ent_min(i)
                  end if

                  qv_in = qu(i,k)
                  ql_in = mu(i,k+1)/mu(i,k)*ql(i,k+1) - &
                              dz(i,k)/mu(i,k)* du(i,k)*ql(i,k+1)
                  qtot = qv_in + ql_in !includes all liquid and vapour contributions to system
                  rcall=3._r8
                  Tfg=t(i,k)
                  call ientropy(rcall,i,lchnk,eta_u(i,k),pf(i,k),qtot,tu(i,k),qsat,Tfg,qv_out,qliq)

                  qu(i,k) = qv_out
                  cu(i,k) = (qliq - ql_in) * mu(i,k)/dz(i,k)
                  ql(i,k) = qliq/ (1._r8+dz(i,k)*c0mask(i)) !reduced from cloud water conversion to rainwater
                  totpcp(i) = totpcp(i) + dz(i,k)*(cu(i,k)-du(i,k)*ql(i,k+1))

!!$                  cpmix = (1 - qu(i,k))*cpres + qu(i,k)*cpwv
!!$                  su(i,k) = tu + grav*zf(i,k)/cpres       !normalised by cpd
!!$                  mass_scaling = (cpmix*tu + grav*zf(i,k))/cpres !!!!USED AS A DUMMY VARIABLE DELETE THIS LINE AFTER TESTING
!!$                  write(iulog,*) "*** Moist convection pre mass scaling: eta_u:", eta_u(i,k), "tu:", tu, "qu:", qu(i,k), "qtot:", qtot, "cu:", cu(i,k), "su using cpd:", su(i,k), "su using cpmix:", mass_scaling

                  rprd(i,k) = c0mask(i)*mu(i,k)*ql(i,k)
                  mass_scaling = 1 / (1 - ql(i,k)*c0mask(i)*dz(i,k))

                  ! if (qliq>1._r8) then
                  !    write(iulog,*) " "
                  !    this_lat = get_rlat_p(lchnk,i)*57.296_r8
                  !    this_lon = get_rlon_p(lchnk,i)*57.296_r8
                  !    write(iulog,*) "*** ZM_CONV: CLDPRP. At lchnk: ", lchnk, "i: ", i, " lat: ", this_lat, "lon: ", this_lon, "k:", k
                  !    write(iulog,*) "*** T:", t(i,k), "P:", p(i,k), "Q:", q(i,k), "Z:", z(i,k),  "s:", s(i,k), "eta:", eta(i,k)
                  !    write(iulog,*) "*** This level: That:", that(i,k), "Pf:", pf(i,k),"qhat:", qhat(i,k),"Zf:", zf(i,k)
                  !    write(iulog,*) "*** This level: h_det:", h_detrain(i,k), "s_det:", s_detrain(i,k), "q_det:", q_detrain(i,k), "eta_det", eta_detrain(i,k)
                  !    write(iulog,*) "*** This level: mu:", mu(i,k), "eu:", eu(i,k), "du:", du(i,k), "dz:", dz(i,k)
                  !    write(iulog,*) "qu:", qu(i,k), "qstu:", qstu, "qtot:", qtot, "tu:", tu(i,k), "eta_u:", eta_u(i,k)
                  !    write(iulog,*) "qliq:", qliq, "ql_in", ql_in, "ql:", ql(i,k), "rprd:", rprd(i,k), "Mstar:", mass_scaling
                  !    write(iulog,*) "Level below: mu:", mu(i,k+1), "hu", hu(i,k+1), "su:", su(i,k+1), "qu:", qu(i,k+1), "eta_u:", eta_u(i,k+1)!, "q:", q(i,k)
                  !    write(iulog,*) " "
                  ! end if

                  ql(i,k) = ql(i,k) * mass_scaling
                  qu(i,k) = qu(i,k) * mass_scaling  !implictly qd is also updated here
                  qtot = ql(i,k)+qu(i,k)
                  call entropy(tu(i,k),p(i,k),qtot, eta_u(i,k))            !adjusted for precipitation change
                  cpmix = (1 - qu(i,k))*cpres + qu(i,k)*cpwv
                  su(i,k) = (cpmix*tu(i,k) + grav*zf(i,k))/cpres       !normalised by cpd
                  hu(i,k) = ((1 - (qu(i,k)+ql(i,k)))*cpres + (qu(i,k)+ql(i,k))*cpliq)*tu(i,k) + grav*zf(i,k) + rl*qu(i,k)

!!$                  write(iulog,*) "*** Moist convection post mass scaling: mass_scaling:", mass_scaling ,"eta_u:", eta_u(i,k), "qtot:", qtot, "qu:", qu(i,k), "su:", su(i,k), "hu:", hu(i,k)
!!$                  write(iulog,*) "*** cpmix:", cpmix, "tu:", tu, "grav:", grav, "zf:", zf(i,k), "cpres:", cpres, "cp:", cp, "rl:", rl, "qu:", qu(i,k), "su:", su(i,k), "hu:", hu(i,k)
!!$                  write(iulog,*) " "
                  ! if (qu(i,k)<0._r8) then
                  !    write(iulog,*) " "
                  !    this_lat = get_rlat_p(lchnk,i)*57.296_r8
                  !    this_lon = get_rlon_p(lchnk,i)*57.296_r8
                  !    write(iulog,*) "*** T:", t(i,k), "P:", p(i,k), "Pf:", pf(i,k), "Q:", q(i,k), "Z:", z(i,k), "Zf:", zf(i,k)
                  !    write(iulog,*) "*** This level: That:", that(i,k), "qhat:", qhat(i,k), "h:", hmn(i,k), "s:", s(i,k), "eta:", eta(i,k)
                  !    write(iulog,*) "*** This level: h_det:", h_detrain(i,k), "s_det:", s_detrain(i,k), "q_det:", q_detrain(i,k), "eta_det", eta_detrain(i,k)
                  !    write(iulog,*) "*** This level: mu:", mu(i,k), "eu:", eu(i,k), "du:", du(i,k), "dz:", dz(i,k)
                  !    ! write(iulog,*) "*** Dry convection estimates: hu:", hu(i,k), "su:", su(i,k), "qu:", qu(i,k), "eta_u:", eta_u(i,k), "tu:", tu(i,k), "qstu:", qstu, "estu:", estu
                  !    write(iulog,*) "qu:", qu(i,k), "qstu:", qstu, "tu:", tu(i,k), "eta_u:", eta_u(i,k)
                  !    ! write(iulog,*) "qhat:", qhat(i,k), "that:", that(i,k)
                  !    write(iulog,*) "cu:", cu(i,k), "qliq:", qliq, "ql_in", ql_in, "ql:", ql(i,k) !"mu:", mu(i,k), "dz:", dz(i,k)
                  !    ! write(iulog,*) "i:", i, "k:", k, "pcp contribution:", dz(i,k)*(cu(i,k)-du(i,k)*ql(i,k+1)), "dz:", dz(i,k), "cu:", cu(i,k),"du:", du(i,k), "ql:", ql(i,k+1)
                  !    write(iulog,*) "Level below: mu:", mu(i,k+1), "hu", hu(i,k+1), "su:", su(i,k+1), "qu:", qu(i,k+1), "eta_u:", eta_u(i,k+1)!, "q:", q(i,k)
                  !    write(iulog,*) " "
                  ! end if
                  !update mu
                  mbm(i,k) = mbm(i,k+1) * (1 - ql(i,k)*c0mask(i)*dz(i,k))
                  mu(i,k) = (mbm(i,k)/eps0(i))* (exp(eps(i,k)*zuef(i))-1._r8)/zuef(i)
               end if

               !additional possible max or minimum plume entropy (entropy could change during condensation)
               if (ent_max(i)<eta_u(i,k)) ent_max(i)=eta_u(i,k)
               if (ent_min(i)>eta_u(i,k)) ent_min(i)=eta_u(i,k)

               !check for plume ending criteria
               tpv(i,k) = tu(i,k)*(1 - mu_red*qu(i,k))
    
               if ((((tv(i,k) > (tpv(i,k)+tiedke_add)) .and. (tv(i,k+1) < (tpv(i,k+1)+tiedke_add))) .or. mu(i,k)<0.01_r8) .and. k<=jb(i)-2) then !plume as a whole not buoyant or too small

                  !if (i==5 .and. lchnk==33) then
                  ! write(iulog,*) " "
                  ! this_lat = get_rlat_p(lchnk,i)*57.296_r8
                  ! this_lon = get_rlon_p(lchnk,i)*57.296_r8
                  ! write(iulog,*) "*** ZM_CONV: CLDPRP. At lchnk: ", lchnk, "i: ", i, " lat: ", this_lat, "lon: ", this_lon
                  ! write(iulog,*) "*** ZM_CONV: CLDPRP. lel:", lel(i), "jt: ", jt(i)
                  ! write(iulog,*) "*** ZM_CONV: plume detraining below or at LCL. tv: ", tv(i,k), "tpv: ", tpv(i,k), "mu:", mu(i,k), "k: ", k
                  ! !write(iulog,*) "*** ZM_CONV: plume detraining. tv: ", tv, "tpv: ", tpv, "k: ", k
                  ! !write(iulog,*) "*** ZM_CONV: plume detraining. tp: ", tu, "qp: ", qu(i,k)
                  ! !write(iulog,*) "*** ZM_CONV: plume detraining. t: ", t(i,k), "q: ", q(i,k)
                  ! write(iulog,*) "*** ZM_CONV: This level: That:", that(i,k), "Pf:", pf(i,k), "Qhat:", qhat(i,k), "hsthat:", hsthat(i,k), "shat:", shat(i,k), "eta:", eta(i,k)
                  ! !write(iulog,*) "*** ZM_CONV: This level: mu:", mu(i,k), "eu:", eu(i,k), "du:", du(i,k)
                  ! write(iulog,*) "*** ZM_CONV: This level: hu:", hu(i,k), "su:", su(i,k), "qu:", qu(i,k), "eta_u:", eta_u(i,k), "tu:", tu(i,k)
                  ! !write(iulog,*) "*** ZM_CONV: This level: h_det:", h_detrain(i,k), "s_det:", s_detrain(i,k), "q_det:", q_detrain(i,k), "eta_det", eta_detrain(i,k)
                  ! write(iulog,*) "*** ZM_CONV: Level below tv:", tv(i,k+1), "tpv:", tpv(i,k+1)
                  ! write(iulog,*) "*** ZM_CONV: Level below: That:", that(i,k+1), "Pf:", pf(i,k+1), "Qhat:", qhat(i,k+1), "hsthat:", hsthat(i,k+1), "shat:", shat(i,k+1), "eta:", eta(i,k+1)
                  ! write(iulog,*) "*** ZM_CONV: Level below: hu", hu(i,k+1), "su:", su(i,k+1), "qu:", qu(i,k+1), "eta_u:", eta_u(i,k+1), "tu:", tu(i,k+1)
                  ! write(iulog,*) " "
                  !end if

                  jt(i) = k ! I think - maybe leave at k, but plume is loses buoyancy somewhere between this level and the last
                  plume_top_reached(i) = .true.
                  lcl_reached(i)= .false.  ! if lcl_reached has been declared true at the detraining level, it has not really been reached
                  totpcp(i) = totpcp(i) - dz(i,k)*(cu(i,k)-du(i,k)*ql(i,k+1))  !remove any precipitation contribution
                  mbm(i,k) = 0
                  mu(i,k) = 0
                  eu(i,k) = 0
                  du(i,k) = mu(i,k+1)/dz(i,k)

                  !fix plume variables at that of level below (implicitly assuming no entrainment, no precipitation)
                  hu(i,k) = hu(i,k+1)
                  su(i,k) = su(i,k+1)
                  eta_u(i,k) = eta_u(i,k+1)
                  qu(i,k) = qu(i,k+1)
                  ql(i,k) = ql(i,k+1)
                  rprd(i,k) = rprd(i,k+1)
                  cu(i,k) = cu(i,k+1)

                  !set plume variables to 0
                  hu(i,k) = 0
                  su(i,k) = 0
                  eta_u(i,k) = 0
                  qu(i,k) = 0
                  ql(i,k) = 0
                  rprd(i,k) = 0
                  cu(i,k) = 0
                  
               end if

            else if ((k < jb(i) .and. k > lel(i)) .and. (lcl_reached(i) .and. .not. plume_top_reached(i))) then
               !masses
               zuef(i) = zf(i,k) - zf(i,jb(i))
               rmue(i) = (mbm(i,k+1)/eps0(i))* (exp(eps(i,k+1)*zuef(i))-1._r8)/zuef(i)
               mu(i,k) = (mbm(i,k+1)/eps0(i))* (exp(eps(i,k  )*zuef(i))-1._r8)/zuef(i)
               eu(i,k) = (rmue(i)-mu(i,k+1))/dz(i,k)
               du(i,k) = (rmue(i)-mu(i,k))/dz(i,k)

               !update max/min q, s, eta
               if (q_max(i)<q(i,k)) q_max(i)=q(i,k)
               if (q_min(i)>q(i,k)) q_min(i)=q(i,k)
               if (s_max(i)<s(i,k)) s_max(i)=s(i,k)
               if (s_min(i)>s(i,k)) s_min(i)=s(i,k)
               if (ent_max(i)<eta(i,k)) ent_max(i)=eta(i,k)
               if (ent_min(i)>eta(i,k)) ent_min(i)=eta(i,k)

               !now need to invert entropy
               eta_u(i,k) = mu(i,k+1)/mu(i,k)*eta_u(i,k+1) + &
                           dz(i,k)/mu(i,k)* (eu(i,k)*eta(i,k) - du(i,k)*eta_detrain(i,k))
               qv_in = mu(i,k+1)/mu(i,k)*qu(i,k+1) + &
                           dz(i,k)/mu(i,k)* (eu(i,k)*q(i,k) - du(i,k)*q_detrain(i,k))
               ql_in = mu(i,k+1)/mu(i,k)*ql(i,k+1) - &
                           dz(i,k)/mu(i,k)* du(i,k)*ql(i,k+1)

               !check on eta and qv_in in updraft to make sure it stays at sensible values
               if (eta_u(i,k) > ent_max(i)) then
                  eta_detrain(i,k) = (eu(i,k)/du(i,k))*eta(i,k) + mu(i,k+1)/(dz(i,k)*du(i,k))*eta_u(i,k+1) - &
                           mu(i,k)/(dz(i,k)*du(i,k))*ent_max(i)
                  eta_u(i,k) = ent_max(i)
               end if
               if (eta_u(i,k) < ent_min(i)) then
                  eta_detrain(i,k) = (eu(i,k)/du(i,k))*eta(i,k) + mu(i,k+1)/(dz(i,k)*du(i,k))*eta_u(i,k+1) - &
                           mu(i,k)/(dz(i,k)*du(i,k))*ent_min(i)
                  eta_u(i,k) = ent_min(i)
               end if

               if (qv_in > q_max(i)) then
                  q_detrain(i,k) = (eu(i,k)/du(i,k))*q(i,k) + mu(i,k+1)/(dz(i,k)*du(i,k))*qu(i,k+1) - &
                         mu(i,k)/(dz(i,k)*du(i,k))*q_max(i)
                  qv_in = q_max(i)
               end if
               if (qv_in < q_min(i)) then
                  q_detrain(i,k) = (eu(i,k)/du(i,k))*q(i,k) + mu(i,k+1)/(dz(i,k)*du(i,k))*qu(i,k+1) - &
                         mu(i,k)/(dz(i,k)*du(i,k))*q_min(i)
                  qu(i,k) = q_min(i)
               end if

               qtot = qv_in + ql_in !includes all liquid and vapour contributions to system

               rcall=3._r8
               Tfg=t(i,k)
               call ientropy(rcall,i,lchnk,eta_u(i,k),pf(i,k),qtot,tu(i,k),qsat,Tfg,qv_out,qliq)
   
               ! if (qliq>1._r8) then
               !    write(iulog,*) " "
               !    this_lat = get_rlat_p(lchnk,i)*57.296_r8
               !    this_lon = get_rlon_p(lchnk,i)*57.296_r8
               !    write(iulog,*) "*** ZM_CONV: Above LCL. At lchnk: ", lchnk, "i: ", i, " lat: ", this_lat, "lon: ", this_lon, "k:", k
               !    write(iulog,*) "*** T:", t(i,k), "P:", p(i,k), "Pf:", pf(i,k), "Q:", q(i,k), "Z:", z(i,k), "Zf:", zf(i,k)
               !    write(iulog,*) "*** This level: That:", that(i,k), "qhat:", qhat(i,k), "h:", hmn(i,k), "s:", s(i,k), "eta:", eta(i,k)
               !    write(iulog,*) "*** This level: h_det:", h_detrain(i,k), "s_det:", s_detrain(i,k), "q_det:", q_detrain(i,k), "eta_det", eta_detrain(i,k)
               !    write(iulog,*) "*** This level: mu:", mu(i,k), "eu:", eu(i,k), "du:", du(i,k), "dz:", dz(i,k)
               !    ! write(iulog,*) "*** Dry convection estimates: hu:", hu(i,k), "su:", su(i,k), "qu:", qu(i,k), "eta_u:", eta_u(i,k), "tu:", tu(i,k), "qstu:", qstu, "estu:", estu
               !    write(iulog,*) "qu:", qu(i,k), "qstu:", qstu, "tu:", tu(i,k), "eta_u:", eta_u(i,k)
               !    ! write(iulog,*) "qhat:", qhat(i,k), "that:", that(i,k)
               !    write(iulog,*) "cu:", cu(i,k), "qliq:", qliq, "ql_in", ql_in, "ql:", ql(i,k) !"mu:", mu(i,k), "dz:", dz(i,k)
               !    ! write(iulog,*) "i:", i, "k:", k, "pcp contribution:", dz(i,k)*(cu(i,k)-du(i,k)*ql(i,k+1)), "dz:", dz(i,k), "cu:", cu(i,k),"du:", du(i,k), "ql:", ql(i,k+1)
               !    write(iulog,*) "Level below: mu:", mu(i,k+1), "hu", hu(i,k+1), "su:", su(i,k+1), "qu:", qu(i,k+1), "eta_u:", eta_u(i,k+1)!, "q:", q(i,k)
               ! end if

               qu(i,k) = qv_out
               cu(i,k) = (qliq - ql_in) * mu(i,k)/dz(i,k)
               ql(i,k) = qliq/ (1._r8+dz(i,k)*c0mask(i)) !reduced from cloud water conversion to rainwater

               rprd(i,k) = c0mask(i)*mu(i,k)*ql(i,k)
   
               mass_scaling = 1 / (1 - ql(i,k)*c0mask(i)*dz(i,k))
               ql(i,k) = ql(i,k) * mass_scaling
               qu(i,k) = qu(i,k) * mass_scaling  !implictly qd is also updated here
               qtot = ql(i,k)+qu(i,k)
               call entropy(tu(i,k),pf(i,k),qtot, eta_u(i,k))            !adjusted for precipitation change
               cpmix = (1 - qu(i,k))*cpres + qu(i,k)*cpwv
               su(i,k) = (cpmix*tu(i,k) + grav*zf(i,k))/cpres       !normalised by cpd
               hu(i,k) = ((1 - (qu(i,k)+ql(i,k)))*cpres + (qu(i,k)+ql(i,k))*cpliq)*tu(i,k) + grav*zf(i,k) + rl*qu(i,k) !different cpmix because defined differently

               ! if (qu(i,k)<0._r8) then
               !    write(iulog,*) " "
               !    this_lat = get_rlat_p(lchnk,i)*57.296_r8
               !    this_lon = get_rlon_p(lchnk,i)*57.296_r8
               !    write(iulog,*) "*** ZM_CONV: Above LCL. At lchnk: ", lchnk, "i: ", i, " lat: ", this_lat, "lon: ", this_lon, "k:", k
               !    write(iulog,*) "*** T:", t(i,k), "P:", p(i,k), "Pf:", pf(i,k), "Q:", q(i,k), "Z:", z(i,k), "Zf:", zf(i,k)
               !    write(iulog,*) "*** This level: That:", that(i,k), "qhat:", qhat(i,k), "h:", hmn(i,k), "s:", s(i,k), "eta:", eta(i,k)
               !    write(iulog,*) "*** This level: h_det:", h_detrain(i,k), "s_det:", s_detrain(i,k), "q_det:", q_detrain(i,k), "eta_det", eta_detrain(i,k)
               !    write(iulog,*) "*** This level: mu:", mu(i,k), "eu:", eu(i,k), "du:", du(i,k), "dz:", dz(i,k)
               !    ! write(iulog,*) "*** Dry convection estimates: hu:", hu(i,k), "su:", su(i,k), "qu:", qu(i,k), "eta_u:", eta_u(i,k), "tu:", tu(i,k), "qstu:", qstu, "estu:", estu
               !    write(iulog,*) "qu:", qu(i,k), "qstu:", qstu, "tu:", tu(i,k), "eta_u:", eta_u(i,k)
               !    ! write(iulog,*) "qhat:", qhat(i,k), "that:", that(i,k)
               !    write(iulog,*) "cu:", cu(i,k), "qliq:", qliq, "ql_in", ql_in, "ql:", ql(i,k) !"mu:", mu(i,k), "dz:", dz(i,k)
               !    ! write(iulog,*) "i:", i, "k:", k, "pcp contribution:", dz(i,k)*(cu(i,k)-du(i,k)*ql(i,k+1)), "dz:", dz(i,k), "cu:", cu(i,k),"du:", du(i,k), "ql:", ql(i,k+1)
               !    write(iulog,*) "Level below: mu:", mu(i,k+1), "hu", hu(i,k+1), "su:", su(i,k+1), "qu:", qu(i,k+1), "eta_u:", eta_u(i,k+1)!, "q:", q(i,k)
               ! end if

               !update mu
               mbm(i,k) = mbm(i,k+1) * (1 - ql(i,k)*c0mask(i)*dz(i,k))
               mu(i,k) = (mbm(i,k)/eps0(i))* (exp(eps(i,k)*zuef(i))-1._r8)/zuef(i)
               

               !check for plume ending criteria
               tpv(i,k) = tu(i,k)*(1 - mu_red*qu(i,k))
               
               !additional possible max or minimum plume entropy
               if (ent_max(i)<eta_u(i,k)) ent_max(i)=eta_u(i,k)
               if (ent_min(i)>eta_u(i,k)) ent_min(i)=eta_u(i,k)

               if ((((tv(i,k) > (tpv(i,k)+tiedke_add)) .and. (tv(i,k+1) < (tpv(i,k+1)+tiedke_add))) .or. mu(i,k)<0.01_r8) .and. k<=jb(i)-2) then !plume as a whole not buoyant or too small

!!$                  write(iulog,*) " "
!!$                  this_lat = get_rlat_p(lchnk,i)*57.296_r8
!!$                  this_lon = get_rlon_p(lchnk,i)*57.296_r8
!!$                  write(iulog,*) "*** ZM_CONV: CLDPRP. At lchnk:", lchnk, "i:", i, " lat:", this_lat, "lon:", this_lon
!!$                  write(iulog,*) "*** ZM_CONV: CLDPRP. lel:", lel(i), "jt:", jt(i)
!!$                  write(iulog,*) "*** ZM_CONV: plume detraining above LCL. tv: ", tv(i,k), "tpv: ", tpv(i,k), "mu:", mu(i,k), "k: ", k
!!$                  !write(iulog,*) "*** ZM_CONV: plume detraining. tp: ", tu, "qp: ", qu(i,k)
!!$                  !write(iulog,*) "*** ZM_CONV: plume detraining. t: ", t(i,k), "q: ", q(i,k)
!!$                  write(iulog,*) "*** ZM_CONV: This level: That:", that(i,k), "Pf:", pf(i,k), "Qhat:", qhat(i,k), "hsthat:", hsthat(i,k), "shat:", shat(i,k), "eta:", eta(i,k)
!!$                  !write(iulog,*) "*** ZM_CONV: This level: mu:", mu(i,k), "eu:", eu(i,k), "du:", du(i,k)
!!$                  write(iulog,*) "*** ZM_CONV: This level: hu:", hu(i,k), "su:", su(i,k), "qu:", qu(i,k), "eta_u:", eta_u(i,k), "tu:", tu(i,k)
!!$                  !write(iulog,*) "*** ZM_CONV: This level: h_det:", h_detrain(i,k), "s_det:", s_detrain(i,k), "q_det:", q_detrain(i,k), "eta_det", eta_detrain(i,k)
!!$                  write(iulog,*) "*** ZM_CONV: Level below tv:", tv(i,k+1), "tpv:", tpv(i,k+1)
!!$                  write(iulog,*) "*** ZM_CONV: Level below: That:", that(i,k+1), "Pf:", pf(i,k+1), "Qhat:", qhat(i,k+1), "hsthat:", hsthat(i,k+1), "shat:", shat(i,k+1), "eta:", eta(i,k+1)
!!$                  write(iulog,*) "*** ZM_CONV: Level below: hu", hu(i,k+1), "su:", su(i,k+1), "qu:", qu(i,k+1), "eta_u:", eta_u(i,k+1), "tu:", tu(i,k+1)
!!$                  write(iulog,*) " "

                  jt(i) = k ! I think - maybe leave at k, but plume is loses buoyancy somewhere between this level and the last
                  plume_top_reached = .true.
                  !change mass fluxes
                  mbm(i,k) = 0
                  mu(i,k) = 0
                  eu(i,k) = 0
                  du(i,k) = mu(i,k+1)/dz(i,k)

                  !fix plume variables at that of level below
                  hu(i,k) = hu(i,k+1)
                  su(i,k) = su(i,k+1)
                  eta_u(i,k) = eta_u(i,k+1)
                  qu(i,k) = qu(i,k+1)
                  ql(i,k) = ql(i,k+1)
                  rprd(i,k) = rprd(i,k+1)
                  cu(i,k) = cu(i,k+1)

                  !set plume variables to 0
                  hu(i,k) = 0
                  su(i,k) = 0
                  eta_u(i,k) = 0
                  qu(i,k) = 0
                  ql(i,k) = 0
                  rprd(i,k) = 0
                  cu(i,k) = 0


               end if

               totpcp(i) = totpcp(i) + dz(i,k)*(cu(i,k)-du(i,k)*ql(i,k+1))

               ! if (lchnk==248 .and. i==9) then
               !    write(iulog,*) "lcl reached. cu:", cu(i,k), "qliq:", qliq, "ql_in", ql_in, "mu:", mu(i,k), "dz:", dz(i,k)
               !    write(iulog,*) "i:", i, "k:", k, "pcp contribution:", dz(i,k)*(cu(i,k)-du(i,k)*ql(i,k+1)), "dz:", dz(i,k), "cu:", cu(i,k),"du:", du(i,k), "ql:", ql(i,k+1)
               ! end if

            end if

            else if (k == lel(i)) then
               !change mass fluxes
               mbm(i,k) = 0
               mu(i,k) = 0
               eu(i,k) = 0
               du(i,k) = mu(i,k+1)/dz(i,k)

               !everything else can be left as 0 or what it is set at when initialising
         end if
      end do
   end do
!
! specify downdraft properties (no downdrafts if jd.ge.jb).
! scale down downward mass flux profile so that net flux
! (up-down) at cloud base in not negative.
!
   do i = 1,il2g
      if (.not. lcl_reached(i)) then 
         jd(i) = jb(i)
      end if
      ! dry_downdraft(i) = .false. !lcl is reached, downdrafts exist and start moistening
   end do

   do i = 1,il2g
      alfa(i) = 0.1_r8 ! in normal downdraft strength run alfa=0.2.  In test4 alfa=0.1
      jt(i) = min(jt(i),jb(i)-1)
      jd(i) = max(j0(i),jt(i)+1)
      jd(i) = min(jd(i),jb(i))
      hd(i,jd(i)) = hmn(i,jd(i)-1)
      if (jd(i) < jb(i) .and. eps0(i) > 0._r8) then
         epsm(i) = eps0(i)
         md(i,jd(i)) = -alfa(i)*epsm(i)/eps0(i)
      end if
   end do

   small = 1.e-20_r8
   do k = msg+2, pver
      do i = 1,il2g
         if (k == jd(i) .and. eps0(i) > 0._r8 .and. jd(i) < jb(i)) then
            !md(i,jd(i)) already set as 1
            eta_d(i,k) = eta(i,k-1)  !starting entropy at interface is that of level above

            !still have to find increase in saturation at the downdraft launch level, it starts off saturated
            rcall = 4._r8
            !if (i==1) rcall = 5._r8
            tfg = t(i,k)
            call ientropy_downdraft(rcall,i,lchnk,eta_d(i,k),pf(i,k),q(i,k) &
               ,t_out,qds(i,k),tfg, dry_downdraft(i))
            dqv = (qds(i,k) - q(i,k))/(1-qds(i,k)) !increase in the water mmr due to the evaporation
            mtm(i,k) = (1 + dqv)  !fixed so that starting mtm=1
            md(i,k) = mtm(i,k) * md(i,k)

            qd(i,k) = qds(i,k)
            evp(i,k-1) = -md(i,k)*(qd(i,k) - q(i,k)) / dz(i,k-1)

            call entropy(t_out,pf(i,k),qd(i,k), eta_d(i,k))  !accounting for new evaporation
            sd(i,k) = (((1-qd(i,k))*cpres + qd(i,k)*cpwv)*t_out + grav*zf(i,k))/cpres
            hd(i,k) = (((1-qd(i,k))*cpres + qd(i,k)*cpliq)*t_out + grav*zf(i,k) + rl*qd(i,k))

            ! if (i==1) then ! .and. lchnk==122) then
            !    write(iulog,*) " "
            !    this_lat = get_rlat_p(lchnk,i)*57.296_r8
            !    this_lon = get_rlon_p(lchnk,i)*57.296_r8
            !    write(iulog,*) "*** ZM_CONV: CLDPRP. At lchnk:", lchnk, "i:", i, " lat:", this_lat, "lon:", this_lon
            !    write(iulog,*) "*** ZM_CONV: CLDPRP. jb:", jb(i), "j0:", j0(i), "jt:", jt(i), "k:", k
            !    write(iulog,*) "This level: T:", t(i,k), "P:", p(i,k), "Q:", q(i,k), "Z:", z(i,k), "h:", hmn(i,k), "s:", s(i,k), "eta:", eta(i,k)
            !    write(iulog,*) "This level: md:", md(i,k), "ed:", ed(i,k), "mtm:", mtm(i,k-1)
            !    write(iulog,*) "This level: zdef:", zdef(i), "prefactor:", -alfa(i)*mtm(i,k-1)/(2._r8*eps0(i)), "exponential:", (exp(2._r8*epsm(i)*zdef(i))-1._r8)/zdef(i) 
            !    write(iulog,*) "i:", i, "k:", k, "evp contribution:", - dz(i,k-1)*ed(i,k)*q(i,k), "dz(k-1):", dz(i,k-1), "ed:", ed(i,k), "q:", q(i,k)
            !    write(iulog,*) " " 
            ! end if
         
         else if (((k > jd(i) .and. k <= jb(i)) .and. (eps0(i) > 0._r8 .and. jd(i) < jb(i))) .and. .not. dry_downdraft(i)) then
         ! else if ((k > jd(i) .and. k <= jb(i)) .and. eps0(i) > 0._r8 .and. jd(i) < jb(i)) then

            zdef(i) = zf(i,jd(i)) - zf(i,k)
            md(i,k) = -alfa(i)*mtm(i,k-1)/ (2._r8*eps0(i))*(exp(2._r8*epsm(i)*zdef(i))-1._r8)/zdef(i) 

               ! if (i==10 .and. lchnk==56) then
               !    write(iulog,*) " "
               !    this_lat = get_rlat_p(lchnk,i)*57.296_r8
               !    this_lon = get_rlon_p(lchnk,i)*57.296_r8
               !    write(iulog,*) "*** ZM_CONV: CLDPRP. At lchnk:", lchnk, "i:", i, " lat:", this_lat, "lon:", this_lon
               !    write(iulog,*) "*** ZM_CONV: CLDPRP. jb:", jb(i), "j0:", j0(i), "jt:", jt(i), "k:", k
               !    write(iulog,*) "This level: T:", t(i,k), "P:", p(i,k), "Q:", q(i,k), "Z:", z(i,k), "h:", hmn(i,k), "s:", s(i,k), "eta:", eta(i,k)
               !    write(iulog,*) "This level: md:", md(i,k), "ed:", ed(i,k), "mtm:", mtm(i,k-1)
               !    write(iulog,*) "This level: zdef:", zdef(i), "prefactor:", -alfa(i)*mtm(i,k-1)/(2._r8*eps0(i)), "exponential:", (exp(2._r8*epsm(i)*zdef(i))-1._r8)/zdef(i) 
               !    write(iulog,*) "Level below: mu:", mu(i,k+1), "hu", hu(i,k+1), "su:", su(i,k+1), "qu:", qu(i,k+1), "eta_u:", eta_u(i,k+1)
               !    !write(iulog,*) " "
               ! end if

            ed(i,k-1) = (md(i,k-1)-md(i,k))/dz(i,k-1)
            mdt = min(md(i,k),-small)

            eta_d(i,k) = (md(i,k-1)*eta_d(i,k-1) - dz(i,k-1)*ed(i,k-1)*eta(i,k-1))/mdt 
            !not including entropy contribution from liquid water about to be evaporated
            qd1 = (md(i,k-1)*qd(i,k-1) - dz(i,k-1)*ed(i,k-1)*q(i,k-1))/mdt !incoming water into the layer pre evaporation
            rcall = 4._r8
            !if (i==1) rcall = 5._r8
            tfg = t(i,k)
            call ientropy_downdraft(rcall,i,lchnk,eta_d(i,k),pf(i,k),qd1 &
               ,t_out,qds(i,k),tfg, dry_downdraft(i))
            dqv = (qds(i,k) - qd1)/(1-qds(i,k)) !increase in the water mmr due to the evaporation
            mtm(i,k) = mtm(i,k-1)*(1 + dqv)
            md(i,k) = md(i,k) * (1 + dqv)

            qd(i,k) = qds(i,k)
            evp(i,k-1) = -md(i,k)*(qd(i,k) - qd1) / dz(i,k-1)

            call entropy(t_out,pf(i,k),qd(i,k), eta_d(i,k))  !accounting for new evaporation
            sd(i,k) = (((1-qd(i,k))*cpres + qd(i,k)*cpwv)*t_out + grav*zf(i,k))/cpres
            hd(i,k) = (((1-qd(i,k))*cpres + qd(i,k)*cpliq)*t_out + grav*zf(i,k) + rl*qd(i,k))
            totevp(i) = totevp(i) - dz(i,k-1)*ed(i,k-1)*q(i,k-1) !works in concert with top / bottom water mass below, do evaporation happening as parcel moved from level above to this level

            ! if (qd(i,k)>0.8_r8) then
            !    write(iulog,*) " "
            !    this_lat = get_rlat_p(lchnk,i)*57.296_r8
            !    this_lon = get_rlon_p(lchnk,i)*57.296_r8
            !    write(iulog,*) "*** ZM_CONV: CLDPRP. At lchnk:", lchnk, "i:", i, " lat:", this_lat, "lon:", this_lon
            !    write(iulog,*) "*** ZM_CONV: CLDPRP. jb:", jb(i), "j0:", j0(i), "jt:", jt(i), "k:", k
            !    write(iulog,*) "Level above: T:", t(i,k-1), "P:", p(i,k-1), "Q:", q(i,k-1), "Z:", z(i,k-1), "h:", hmn(i,k-1), "s:", s(i,k-1), "eta:", eta(i,k-1)
            !    write(iulog,*) "Level above: md:", md(i,k-1), "qd:", qd(i,k-1), "eta_d:", eta_d(i,k-1)
            !    write(iulog,*) "This level: T:", t(i,k), "P:", p(i,k), "Q:", q(i,k), "Z:", z(i,k), "h:", hmn(i,k), "s:", s(i,k), "eta:", eta(i,k)
            !    write(iulog,*) "This level: pre-evap md:", md(i,k), "ed:", ed(i,k), "mtm:", mtm(i,k-1)
            !    write(iulog,*) "This level: zdef:", zdef(i), "prefactor:", -alfa(i)*mtm(i,k-1)/(2._r8*eps0(i)), "exponential:", (exp(2._r8*epsm(i)*zdef(i))-1._r8)/zdef(i) 
            !    !write(iulog,*) " "
            !    write(iulog,*) "Post evaporation: qd1:", qd1, "qds:", qds(i,k), "dqv:", dqv, "mtm:", mtm(i,k), "md:", md(i,k)
            !    write(iulog,*) "i:", i, "k:", k, "evp contribution:", - dz(i,k-1)*ed(i,k-1)*q(i,k-1), "dz:", dz(i,k-1), "ed:", ed(i,k-1), "q:", q(i,k-1)
            !    write(iulog,*) " "
            ! end if

            if (dry_downdraft(i)) then !redo layer assuming that the downdraft is dry from here
               zdef(i) = zf(i,jd(i)) - zf(i,k)
               md(i,k) = -alfa(i)*mtm(i,k-1)/ (2._r8*eps0(i))*(exp(2._r8*epsm(i)*zdef(i))-1._r8)/zdef(i) 
               mtm(i,k) = mtm(i,k-1)
               ed(i,k-1) = (md(i,k-1)-md(i,k))/dz(i,k-1)
               mdt = min(md(i,k),-small)

               eta_d(i,k) = (md(i,k-1)*eta_d(i,k-1) - dz(i,k-1)*ed(i,k-1)*eta(i,k-1))/mdt 
               !not including entropy contribution from liquid water about to be evaporated
               qd(i,k) = (md(i,k-1)*qd(i,k-1) - dz(i,k-1)*ed(i,k-1)*q(i,k-1))/mdt !incoming water into the layer pre evaporation
               !as we are adiabatically warming the downdraft, it will all stay in vapour phase
               rcall = 5._r8
               tfg = t(i,k)
               call ientropy(rcall,i,lchnk,eta_d(i,k),pf(i,k),qd(i,k) &
                  ,t_out,qds(i,k),tfg, qv_out, qliq) !rcall,icol,lchnk,s,p,qt,T,qsat,Tfg, qv,qc
                  
               evp(i,k-1) = 0._r8 !qd(i,k)=qd1, evaporation has happened in higher up layers
               !no new condensation or evaporation so no entropy change in theory, run it anyway
               call entropy(t_out,pf(i,k),qd(i,k), eta_d(i,k))  !not actually sure we track temp of downdraft - could be worth doing?
               sd(i,k) = (((1-qd(i,k))*cpres + qd(i,k)*cpwv)*t_out + grav*zf(i,k))/cpres
               hd(i,k) = (((1-qd(i,k))*cpres + qd(i,k)*cpliq)*t_out + grav*zf(i,k) + rl*qd(i,k))
               totevp(i) = totevp(i) - dz(i,k-1)*ed(i,k-1)*q(i,k-1) 
               !works in concert with top / bottom water mass below, do evaporation happening as parcel moved from level above to this level
               !means it can be non-zero even though we are not evaporating anything at this level
            end if
         else if (((k > jd(i) .and. k <= jb(i)) .and. (eps0(i) > 0._r8 .and. jd(i) < jb(i))) .and. (dry_downdraft(i))) then
               zdef(i) = zf(i,jd(i)) - zf(i,k)
               md(i,k) = -alfa(i)*mtm(i,k-1)/ (2._r8*eps0(i))*(exp(2._r8*epsm(i)*zdef(i))-1._r8)/zdef(i) 
               mtm(i,k) = mtm(i,k-1)
               ed(i,k-1) = (md(i,k-1)-md(i,k))/dz(i,k-1)
               mdt = min(md(i,k),-small)

               eta_d(i,k) = (md(i,k-1)*eta_d(i,k-1) - dz(i,k-1)*ed(i,k-1)*eta(i,k-1))/mdt 
               !not including entropy contribution from liquid water about to be evaporated
               qd(i,k) = (md(i,k-1)*qd(i,k-1) - dz(i,k-1)*ed(i,k-1)*q(i,k-1))/mdt !incoming water into the layer pre evaporation
               !as we are adiabatically warming the downdraft, it will all stay in vapour phase
               rcall = 5._r8
               tfg = t(i,k)
               call ientropy(rcall,i,lchnk,eta_d(i,k),pf(i,k),qd(i,k) &
                  ,t_out,qds(i,k),tfg, qv_out, qliq) !rcall,icol,lchnk,s,p,qt,T,qsat,Tfg, qv,qc
                  
               evp(i,k-1) = 0._r8 !qd(i,k)=qd1, evaporation has happened in higher up layers
               !no new condensation or evaporation so no entropy change in theory, run it anyway
               call entropy(t_out,pf(i,k),qd(i,k), eta_d(i,k))  !not actually sure we track temp of downdraft - could be worth doing?
               sd(i,k) = (((1-qd(i,k))*cpres + qd(i,k)*cpwv)*t_out + grav*zf(i,k))/cpres
               hd(i,k) = (((1-qd(i,k))*cpres + qd(i,k)*cpliq)*t_out + grav*zf(i,k) + rl*qd(i,k))
               totevp(i) = totevp(i) - dz(i,k-1)*ed(i,k-1)*q(i,k-1) 
               !works in concert with top / bottom water mass below, do evaporation happening as parcel moved from level above to this level
               !means it can be non-zero even though we are not evaporating anything at this level
         end if
      end do
   end do




   do k = msg + 1,pver
      do i = 1,il2g
         if ((k >= jt(i) .and. k <= jb(i)) .and. eps0(i) > 0._r8 .and. jd(i) < jb(i)) then
            ratmjb(i) = min(abs(mu(i,jb(i))/md(i,jb(i))),1._r8)
            md(i,k) = md(i,k)*ratmjb(i)    
            ed(i,k) = ed(i,k)*ratmjb(i)
            evp(i,k)=evp(i,k)*ratmjb(i)      
         end if
      end do
   end do
!
   do i = 1,il2g
      qd(i,jd(i)) = qds(i,jd(i))
      sd(i,jd(i)) = (hd(i,jd(i)) - rl*qd(i,jd(i)))/cpres
   end do

do i = 1,il2g
      totevp(i) = totevp(i) + md(i,jd(i))*qd(i,jd(i)) - md(i,jb(i))*qd(i,jb(i))

      ! if (lchnk==56 .and. i==10) then
      !    write(iulog,*) " "
      !    write(iulog,*) "new physics. lchnk:", lchnk, "i:", i, "totevp:", totevp(i), "totpcp:", totpcp(i)
      !    write(iulog,*) "q_jd:", md(i,jd(i))*qd(i,jd(i)), "q_jb:", md(i,jb(i))*qd(i,jb(i))
      ! end if
      totpcp(i) = max(totpcp(i),0._r8)
      totevp(i) = max(totevp(i),0._r8)
   end do
!
   do k = msg + 2,pver
      do i = 1,il2g
         if (totevp(i) > 0._r8 .and. totpcp(i) > 0._r8) then
            md(i,k)  = md (i,k)*min(1._r8, totpcp(i)/(totevp(i)+totpcp(i)))
            ed(i,k)  = ed (i,k)*min(1._r8, totpcp(i)/(totevp(i)+totpcp(i)))
            evp(i,k) = evp(i,k)*min(1._r8, totpcp(i)/(totevp(i)+totpcp(i)))
         else
            md(i,k) = 0._r8
            ed(i,k) = 0._r8
            evp(i,k) = 0._r8
         end if
            ! cmeg is the cloud water condensed - rain water evaporated
            ! rprd is the cloud water converted to rain - (rain evaporated)
            cmeg(i,k) = cu(i,k) - evp(i,k)
            rprd(i,k) = rprd(i,k)-evp(i,k)
      end do
   end do

! compute the net precipitation flux across interfaces
   pflx(:il2g,1) = 0._r8
   do k = 2,pverp
      do i = 1,il2g
         pflx(i,k) = pflx(i,k-1) + rprd(i,k-1)*dz(i,k-1)
      end do
   end do
!
   do k = msg + 1,pver
      do i = 1,il2g
         mc(i,k) = mu(i,k) + md(i,k)
      end do
   end do
!

return
end subroutine cldprp

subroutine closure(lchnk   , &
    q       ,t       ,p       ,z       ,s       , &
    tp      ,qs      ,qu      ,su      ,mc      , &
    du      ,mu      ,md      ,qd      ,sd      , &
    qhat    ,shat    ,dp      ,qstp    ,pf      ,zf      , &
    ql      ,dsubcld ,mb      ,cape    ,tl      , &
    lcl     ,lel     ,jt      ,mx      ,il1g    , &
    il2g    ,rd      ,grav    ,cp      ,rl      , &
    msg     ,capelmt )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! <Say what the routine does> 
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: G. Zhang and collaborators. CCM contact:P. Rasch
! This is contributed code not fully standardized by the CCM core group.
!
! this code is very much rougher than virtually anything else in the CCM
! We expect to release cleaner code in a future release
!
! the documentation has been enhanced to the degree that we are able
! 
!-----------------------------------------------------------------------
use dycore,    only: dycore_is, get_resolution
use phys_grid, only: get_rlon_p, get_rlat_p

implicit none

!
!-----------------------------Arguments---------------------------------
!
integer, intent(in) :: lchnk                 ! chunk identifier

real(r8), intent(inout) :: q(pcols,pver)        ! spec humidity
real(r8), intent(inout) :: t(pcols,pver)        ! temperature
real(r8), intent(inout) :: p(pcols,pver)        ! pressure (mb)
real(r8), intent(inout) :: mb(pcols)            ! cloud base mass flux
real(r8), intent(in) :: z(pcols,pver)        ! height (m)
real(r8), intent(in) :: s(pcols,pver)        ! normalized dry static energy
real(r8), intent(in) :: tp(pcols,pver)       ! parcel temp
real(r8), intent(in) :: qs(pcols,pver)       ! sat spec humidity
real(r8), intent(in) :: qu(pcols,pver)       ! updraft spec. humidity
real(r8), intent(in) :: su(pcols,pver)       ! normalized dry stat energy of updraft
real(r8), intent(in) :: mc(pcols,pver)       ! net convective mass flux
real(r8), intent(in) :: du(pcols,pver)       ! detrainment from updraft
real(r8), intent(in) :: mu(pcols,pver)       ! mass flux of updraft
real(r8), intent(in) :: md(pcols,pver)       ! mass flux of downdraft
real(r8), intent(in) :: qd(pcols,pver)       ! spec. humidity of downdraft
real(r8), intent(in) :: sd(pcols,pver)       ! dry static energy of downdraft
real(r8), intent(in) :: qhat(pcols,pver)     ! environment spec humidity at interfaces
real(r8), intent(in) :: shat(pcols,pver)     ! env. normalized dry static energy at intrfcs
real(r8), intent(in) :: dp(pcols,pver)       ! pressure thickness of layers
real(r8), intent(in) :: qstp(pcols,pver)     ! spec humidity of parcel
real(r8), intent(in) :: pf(pcols,pver+1)     ! pressure of interface levels
real(r8), intent(in) :: zf(pcols,pver+1)     ! height of interface levels
real(r8), intent(in) :: ql(pcols,pver)       ! liquid water mixing ratio

real(r8), intent(in) :: cape(pcols)          ! available pot. energy of column
real(r8), intent(in) :: tl(pcols)
real(r8), intent(in) :: dsubcld(pcols)       ! thickness of subcloud layer

integer, intent(in) :: lcl(pcols)        ! index of lcl
integer, intent(in) :: lel(pcols)        ! index of launch leve
integer, intent(in) :: jt(pcols)         ! top of updraft
integer, intent(in) :: mx(pcols)         ! base of updraft
!
!--------------------------Local variables------------------------------
!
real(r8) dtpdt(pcols,pver)
real(r8) dqpdt(pcols,pver)
real(r8) dqsdtp(pcols,pver)
real(r8) dtmdt(pcols,pver)
real(r8) dqmdt(pcols,pver)
real(r8) dboydt(pcols,pver)
real(r8) thetavp(pcols,pver)
real(r8) thetavm(pcols,pver)

real(r8) dtbdt(pcols),dqbdt(pcols),dtldt(pcols)
real(r8) beta
real(r8) capelmt
real(r8) cp
real(r8) dadt(pcols)
real(r8) debdt
real(r8) dltaa
real(r8) eb
real(r8) grav
real(r8) zvirp1

real(r8) cpmix  !bulk heat capacity of gas, including background gas and water vapour
real(r8) dt,dq  !temperature and moisture perturbations to see effect of these starting quantities on parcel buoyancy
real(r8) pl(pcols) !exact pressure level of LCL
real(r8) tl_dummy(pcols)  !exact temperature level of LCL
integer lcl_dummy(pcols) !exact level of LCL
logical leave_column(pcols)

real(r8) q1(pcols,pver)
real(r8) t1(pcols,pver)
real(r8) tpt1(pcols,pver)
real(r8) tpq1(pcols,pver)
real(r8) qpt1(pcols,pver)
real(r8) qpq1(pcols,pver)
real(r8) tpert(pcols,pver)
real(r8) tpv(pcols,pver)

integer i
integer il1g
integer il2g
integer k, kmin, kmax
integer msg

real(r8) rd
real(r8) rl

real(r8) this_lat, this_lon

! change of subcloud layer properties due to convection is
! related to cumulus updrafts and downdrafts.
! mc(z)=f(z)*mb, mub=betau*mb, mdb=betad*mb are used
! to define betau, betad and f(z).
! note that this implies all time derivatives are in effect
! time derivatives per unit cloud-base mass flux, i.e. they
! have units of 1/mb instead of 1/sec.
!
zvirp1 = zvir + 1.0


do i = il1g,il2g
   mb(i) = 0._r8
   eb = p(i,mx(i))*q(i,mx(i))/ (eps1+q(i,mx(i)))
   dtbdt(i) = (1._r8/dsubcld(i))* (mu(i,mx(i))*(shat(i,mx(i))-su(i,mx(i)))+ &
      md(i,mx(i))* (shat(i,mx(i))-sd(i,mx(i))))
   dqbdt(i) = (1._r8/dsubcld(i))* (mu(i,mx(i))*(qhat(i,mx(i))-qu(i,mx(i)))+ &
   md(i,mx(i))* (qhat(i,mx(i))-qd(i,mx(i))))
   debdt = eps1*p(i,mx(i))/ (eps1+q(i,mx(i)))**2*dqbdt(i)
   dtldt(i) = -2840._r8* (3.5_r8/t(i,mx(i))*dtbdt(i)-debdt/eb)/ &
   (3.5_r8*log(t(i,mx(i)))-log(eb)-4.805_r8)**2
end do
!
!   dtmdt and dqmdt are cumulus heating and drying.
!
do k = msg + 1,pver
   do i = il1g,il2g
      dtmdt(i,k) = 0._r8
      dqmdt(i,k) = 0._r8
      dboydt(i,k) = 0._r8
      dtpdt(i,k) = 0._r8
      dqpdt(i,k) = 0._r8
   end do
end do
!
do k = msg + 1,pver - 1
   do i = il1g,il2g
      if (k == jt(i)) then
         dtmdt(i,k) = (1._r8/dp(i,k))*(mu(i,k+1)* (su(i,k+1)-shat(i,k+1)- & !are there not values we can use from level k?
                  rl/cp*ql(i,k+1))+md(i,k+1)* (sd(i,k+1)-shat(i,k+1)))
         dqmdt(i,k) = (1._r8/dp(i,k))*(mu(i,k+1)* (qu(i,k+1)- &
                  qhat(i,k+1)+ql(i,k+1))+md(i,k+1)*(qd(i,k+1)-qhat(i,k+1)))
      end if
   end do
end do
!
beta = 0._r8
do k = msg + 1,pver - 1
   do i = il1g,il2g
      if (k > jt(i) .and. k < mx(i)) then
         dtmdt(i,k) = (mc(i,k)* (shat(i,k)-s(i,k))+mc(i,k+1)* (s(i,k)-shat(i,k+1)))/ &
                  dp(i,k) - rl/cp*du(i,k)*(beta*ql(i,k)+ (1-beta)*ql(i,k+1))
         !          dqmdt(i,k)=(mc(i,k)*(qhat(i,k)-q(i,k))
         !     1                +mc(i,k+1)*(q(i,k)-qhat(i,k+1)))/dp(i,k)
         !     2                +du(i,k)*(qs(i,k)-q(i,k))
         !     3                +du(i,k)*(beta*ql(i,k)+(1-beta)*ql(i,k+1))

         dqmdt(i,k) = (mu(i,k+1)* (qu(i,k+1)-qhat(i,k+1)+cp/rl* (su(i,k+1)-s(i,k)))- &
                  mu(i,k)* (qu(i,k)-qhat(i,k)+cp/rl*(su(i,k)-s(i,k)))+md(i,k+1)* &
                  (qd(i,k+1)-qhat(i,k+1)+cp/rl*(sd(i,k+1)-s(i,k)))-md(i,k)* &
                  (qd(i,k)-qhat(i,k)+cp/rl*(sd(i,k)-s(i,k))))/dp(i,k) + &
                  du(i,k)* (beta*ql(i,k)+(1-beta)*ql(i,k+1))
      end if
   end do
end do


!
dt = 0.1
dq = 0.0001
do i = il1g, il2g
   tpt1(i,:) = tp(i,:)  !parcel temperature resulting from increasing the local temperature
   qpt1(i,:) = qstp(i,:)  !parcel q from increasing local T
   tpq1(i,:) = tp(i,:)  !parcel T from increasing local q
   qpq1(i,:) = qstp(i,:)  !parcel q from increasing local q
   t1(i,:) = t(i,:) + dt !0.1 degree increase in temp should be small enough
   q1(i,:) = q(i,:) + dq !smaller change in q because it varies between 0 and 1, and saturation humidity is very low at low T anyway
   lcl_dummy(i) = mx(i)
   tl_dummy(i) = t(i,mx(i))
   pl(i) = p(i,mx(i))
   leave_column(i) = .false. !we want to work through all the columns in parcel_dilute
end do

call parcel_dilute(lchnk, il2g, msg, mx, p, t1, q, tpert, tpt1, tpv, qpt1, pl, tl_dummy, lcl_dummy, leave_column) !pcols as well
call parcel_dilute(lchnk, il2g, msg, mx, p, t, q1, tpert, tpq1, tpv, qpq1, pl, tl_dummy, lcl_dummy, leave_column)

do k = msg + 1,pver
   do i = il1g,il2g
      if (k >= jt(i) .and. k < mx(i)) then !going over all the plume
         dqpdt(i,k) = ((qpq1(i,k)-qstp(i,k))/dq * dqbdt(i)) + ((qpt1(i,k)-qstp(i,k))/dt * dtbdt(i))
         !dqpdqb * dqbdt + dqpdtb * dtbdt
         dtpdt(i,k) = ((tpt1(i,k)-tp(i,k))/dt * dtbdt(i)) + ((tpq1(i,k)-tp(i,k))/dq * dqbdt(i))
         !dtpdtb * dtbdt +dtpdqb * dqbdt

         dboydt(i,k) = (1-mu_red*qstp(i,k))*dtpdt(i,k) - mu_red*tp(i,k)*dqpdt(i,k) &
         - (1-mu_red*q(i,k))*dtmdt(i,k) - mu_red*t(i,k)*dqmdt(i,k)
      end if
   end do
end do

!
! buoyant energy change is set to 2/3*excess cape per 3 hours
!
dadt(il1g:il2g)  = 0._r8
kmin = minval(lel(il1g:il2g))
kmax = maxval(mx(il1g:il2g)) - 1
do k = kmin, kmax
   do i = il1g,il2g
      if ( k >= lel(i) .and. k <= mx(i) - 1) then
         dadt(i) = dadt(i) + rd * dboydt(i,k) * log(pf(i,k+1)/pf(i,k))
      endif
   end do
end do
do i = il1g,il2g
   dltaa = -1._r8* (cape(i)-capelmt)
   if (dadt(i) /= 0._r8) mb(i) = max(dltaa/tau/dadt(i),0._r8)
end do

! do i = il2g,il2g
!    this_lat = get_rlat_p(lchnk,i)*57.296_r8
!    this_lon = get_rlon_p(lchnk,i)*57.296_r8
!    ! if (lchnk==101 .and. i==1) then
!       !write(iulog,*) " "
!       write(iulog,*) " "
!       write(iulog,*) "*** ZM_CONV: NEW CLOSURE. At lchnk: ", lchnk, "i: ", i, " lat: ", this_lat, "lon: ", this_lon
!       write(iulog,*) "*** ZM_CONV: NEW CLOSURE. dadt:", dadt(i), "mb:", mb(i), "lel:", lel(i), "mx:", mx(i)
!       write(iulog,*) "*** ZM_CONV: NEW CLOSURE. dtbdt:", dtbdt(i), "dqbdt:", dqbdt(i)
!       write(iulog,*) " "
!       !do k = kmax+1,kmin,-1
!       ! do k = kmax+1, 47, -1
!       !    write(iulog,*) " "
!       !    write(iulog,*) "*** Going up column. dboydt: ", dboydt(i,k), "env change:", -(1-mu_red*q(i,k))*dtmdt(i,k) - mu_red*t(i,k)*dqmdt(i,k), "parcel change:", (1-mu_red*qstp(i,k))*dtpdt(i,k) - mu_red*tp(i,k)*dqpdt(i,k), "k:", k
!       !    write(iulog,*) "*** Going up column. dtmdt: ", dtmdt(i,k), "dqmdt: ", dqmdt(i,k), "k:", k
!       !    write(iulog,*) "*** Going up column. dtpdt:", dtpdt(i,k), "dqpdt:", dqpdt(i,k)
!       !    write(iulog,*) "*** mu:", mu(i,k), "md:", md(i,k), "mc:", mc(i,k), "dp:", dp(i,k), "du:", du(i,k)
!       !    write(iulog,*) "*** s:", s(i,k), "shat:", shat(i,k), "su:", su(i,k), "sd", sd(i,k)
!       !    write(iulog,*) "*** qhat:", qhat(i,k), "ql:", ql(i,k), "qu", qu(i,k), "qd:", qd(i,k)
!       !    write(iulog,*) "*** updraft into layer:", mu(i,k+1)*(su(i,k+1)-shat(i,k+1)), "updraft out of layer:", mu(i,k)*(su(i,k)-shat(i,k))
!       !    write(iulog,*) "*** downdraft into layer:", -md(i,k)*(sd(i,k)-shat(i,k)), "downdraft out of layer:", -md(i,k+1)*(sd(i,k+1)-shat(i,k+1))
!       !    write(iulog,*) "*** evaporating detraining liquid:", rl/cp*du(i,k)*ql(i,k+1)*dp(i,k)
!       !    write(iulog,*) "*** q change q result:", qpq1(i,k)-qstp(i,k), "t change q result:", qpt1(i,k)-qstp(i,k)
!       !    write(iulog,*) "*** t change t result:", tpt1(i,k)-tp(i,k), "q change t result:", tpq1(i,k)-tp(i,k)
!       ! end do
!    ! end if
! end do

return
end subroutine closure

subroutine q1q2_pjr(lchnk   , &
                    dqdt    ,dsdt    ,q       ,qs      ,qu      , &
                    su      ,du      ,qhat    ,shat    ,dp      , &
                    mu      ,md      ,sd      ,qd      ,ql      , &
                    dsubcld ,jt      ,mx      ,il1g    ,il2g    , &
                    cp      ,rl      ,msg     ,          &
                    dl      ,evp     ,cu      )


   use phys_grid, only: get_rlon_p, get_rlat_p
  
   implicit none

!----------------------------------------------------------------------- 
! 
! Purpose: 
! <Say what the routine does> 
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: phil rasch dec 19 1995
! 
!-----------------------------------------------------------------------


   real(r8), intent(in) :: cp

   integer, intent(in) :: lchnk             ! chunk identifier
   integer, intent(in) :: il1g
   integer, intent(in) :: il2g
   integer, intent(in) :: msg

   real(r8), intent(in) :: q(pcols,pver)
   real(r8), intent(in) :: qs(pcols,pver)
   real(r8), intent(in) :: qu(pcols,pver)
   real(r8), intent(in) :: su(pcols,pver)
   real(r8), intent(in) :: du(pcols,pver)
   real(r8), intent(in) :: qhat(pcols,pver)
   real(r8), intent(in) :: shat(pcols,pver)
   real(r8), intent(in) :: dp(pcols,pver)
   real(r8), intent(in) :: mu(pcols,pver)
   real(r8), intent(in) :: md(pcols,pver)
   real(r8), intent(in) :: sd(pcols,pver)
   real(r8), intent(in) :: qd(pcols,pver)
   real(r8), intent(in) :: ql(pcols,pver)
   real(r8), intent(in) :: evp(pcols,pver)
   real(r8), intent(in) :: cu(pcols,pver)
   real(r8), intent(in) :: dsubcld(pcols)

   real(r8),intent(out) :: dqdt(pcols,pver),dsdt(pcols,pver)
   real(r8),intent(out) :: dl(pcols,pver)
   integer kbm
   integer ktm
   integer jt(pcols)
   integer mx(pcols)
!
! work fields:
!
   integer i
   integer k

   real(r8) emc
   real(r8) rl

   real(r8) this_lat, this_lon
!-------------------------------------------------------------------
   do k = msg + 1,pver
      do i = il1g,il2g
         dsdt(i,k) = 0._r8
         dqdt(i,k) = 0._r8
         dl(i,k) = 0._r8
      end do
   end do
!
! find the highest level top and bottom levels of convection
!
   ktm = pver
   kbm = pver
   do i = il1g, il2g
      ktm = min(ktm,jt(i))
      kbm = min(kbm,mx(i))
   end do

   do k = ktm,pver-1
      do i = il1g,il2g
         emc = -cu (i,k)               &         ! condensation in updraft
               +evp(i,k)                         ! evaporating rain in downdraft

         dsdt(i,k) = -rl/cp*emc &
                     + (+mu(i,k+1)* (su(i,k+1)-shat(i,k+1)) &
                        -mu(i,k)*   (su(i,k)-shat(i,k)) &
                        +md(i,k+1)* (sd(i,k+1)-shat(i,k+1)) &
                        -md(i,k)*   (sd(i,k)-shat(i,k)) &
                       )/dp(i,k)

         dqdt(i,k) = emc + &
                    (+mu(i,k+1)* (qu(i,k+1)-qhat(i,k+1)) &
                     -mu(i,k)*   (qu(i,k)-qhat(i,k)) &
                     +md(i,k+1)* (qd(i,k+1)-qhat(i,k+1)) &
                     -md(i,k)*   (qd(i,k)-qhat(i,k)) &
                    )/dp(i,k)

         dl(i,k) = du(i,k)*ql(i,k+1)

      end do
   end do

!
!DIR$ NOINTERCHANGE!
   do k = kbm,pver
      do i = il1g,il2g
         if (k == mx(i)) then
            dsdt(i,k) = (1._r8/dsubcld(i))* &
                        (-mu(i,k)* (su(i,k)-shat(i,k)) &
                         -md(i,k)* (sd(i,k)-shat(i,k)) &
                        )
            dqdt(i,k) = (1._r8/dsubcld(i))* &
                        (-mu(i,k)*(qu(i,k)-qhat(i,k)) &
                         -md(i,k)*(qd(i,k)-qhat(i,k)) &
                        )
         else if (k > mx(i) .and. k <=mx(i)+2) then
            dsdt(i,k) = dsdt(i,k-1)
            dqdt(i,k) = dqdt(i,k-1)
         end if
      end do
   end do
   !
   !
! if (lchnk==56) write(iulog,*) "il1g:", il1g, "il2g:", il2g
! do i = il1g,il2g
!    ! if (lchnk==60) write(iulog,*) "i:", i
!    this_lat = get_rlat_p(lchnk,i)*57.296_r8
!    this_lon = get_rlon_p(lchnk,i)*57.296_r8
!    !if (maxval(dsdt(i,:))>0.01_r8) then
!    if (lchnk==39 .and. i==11) then
!       write(iulog,*) " "
!       write(iulog,*) " "
!       write(iulog,*) "*** ZM_CONV: q1q2_pjr. At lchnk:", lchnk, "i: ", i, " lat: ", this_lat, "lon: ", this_lon
!       write(iulog,*) "*** ZM_CONV: q1q2_pjr. mx:", mx(i), "jt:", jt(i)
!       do k = mx(i), jt(i), -1
!          write(iulog,*) " "
!          write(iulog,*) "*** ZM_CONV: q1q2_pjr. k:", k
!          write(iulog,*) "*** ZM_CONV: q1q2_pjr. dsdt:", dsdt(i,k), "dqdt:", dqdt(i,k), "dl:", dl(i,k)
!          write(iulog,*) "*** Going up column. mu:", mu(i,k), "md:", md(i,k), "dp:", dp(i,k)
!          write(iulog,*) "*** Going up column, env properties. shat:", shat(i,k), "qhat:", qhat(i,k), "su:", su(i,k), "sd", sd(i,k)
!          write(iulog,*) "*** Going up column. updraft properties. su", su(i,k), "qu:", qu(i,k)
!          write(iulog,*) "*** Going up column. downdraft properties. sd", sd(i,k), "qd:", qd(i,k)
!          write(iulog,*) "*** Going up column. emc:", evp(i,k)-cu(i,k), "latent heat:", -rl/cp*(evp(i,k)-cu(i,k))
!          write(iulog,*) "updraft s leaving level:", mu(i,k)*(su(i,k)-shat(i,k)), "downdraft s joining level:", md(i,k)*(sd(i,k)-shat(i,k))
!          write(iulog,*) "updraft q leaving level:", mu(i,k)*(qu(i,k)-qhat(i,k)), "downdraft q joining level:", md(i,k)*(qd(i,k)-qhat(i,k))
!       end do
!       write(iulog,*) " "
!    end if
! end do

   
   return
end subroutine q1q2_pjr

subroutine buoyan_dilute(lchnk   ,ncol    , &
                  q       ,t       ,p       ,z       ,pf      , &
                  tp      ,qstp    ,tl      ,rl      ,cape    , &
                  pblt    ,lcl     ,lel     ,lon     ,mx      , &
                  mx_2    ,rd      ,grav    ,cp      ,msg     , &
                  tpert   ,tp_2    ,qstp_2  ,cape_2  ,lel_2   ,lcl_2  ,sfcT)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Calculates CAPE the lifting condensation level and the convective top
! where buoyancy is first -ve.
! 
! Method: Calculates the parcel temperature based on a simple constant
! entraining plume model. CAPE is integrated from buoyancy.
! 09/09/04 - Simplest approach using an assumed entrainment rate for 
!            testing (dmpdp). 
! 08/04/05 - Swap to convert dmpdz to dmpdp  
!
! SCAM Logical Switches - DILUTE:RBN - Now Disabled 
! ---------------------
! switch(1) = .T. - Uses the dilute parcel calculation to obtain tendencies.
! switch(2) = .T. - Includes entropy/q changes due to condensate loss and freezing.
! switch(3) = .T. - Adds the PBL Tpert for the parcel temperature at all levels.
! 
! References:
! Raymond and Blythe (1992) JAS 
! 
! Author:
! Richard Neale - September 2004
! 
!-----------------------------------------------------------------------
   use phys_grid, only: get_rlon_p, get_rlat_p

   implicit none
!-----------------------------------------------------------------------
!
! input arguments
!
   integer, intent(in) :: lchnk                 ! chunk identifier
   integer, intent(in) :: ncol                  ! number of atmospheric columns

   real(r8), intent(in) :: q(pcols,pver)        ! spec. humidity
   real(r8), intent(in) :: t(pcols,pver)        ! temperature
   real(r8), intent(in) :: p(pcols,pver)        ! pressure
   real(r8), intent(in) :: z(pcols,pver)        ! height
   real(r8), intent(in) :: pf(pcols,pver+1)     ! pressure at interfaces
   real(r8), intent(in) :: pblt(pcols)          ! index of pbl depth
   real(r8), intent(in) :: tpert(pcols)         ! perturbation temperature by pbl processes
   real(r8), intent(in) :: sfcT(pcols)          ! surface temperature from cam_in

!
! output arguments
!
   real(r8), intent(out) :: tp(pcols,pver)       ! parcel temperature
   real(r8), intent(out) :: tp_2(pcols,pver)     ! parcel T, second region
   real(r8), intent(out) :: qstp(pcols,pver)     ! saturation mixing ratio of parcel (only above lcl, just q below).
   real(r8), intent(out) :: qstp_2(pcols,pver)   ! parcel q, second region
   real(r8), intent(out) :: tl(pcols)            ! parcel temperature at lcl
   real(r8), intent(out) :: cape(pcols)          ! convective aval. pot. energy.
   real(r8), intent(out) :: cape_2(pcols)        ! CAPE - second convective region
   integer, intent(out) :: lcl(pcols)        !
   integer, intent(out) :: lcl_2(pcols)        !
   integer, intent(out) :: lel(pcols)        !
   integer, intent(out) :: lel_2(pcols)
   integer, intent(out) :: lon(pcols)        ! level of onset of deep convection
   integer, intent(out) :: mx(pcols)         ! level of first instability to convection
   integer, intent(out) :: mx_2(pcols)       ! level of first instability to convection - second convection
   ! integer mx_sat(pcols)     ! level of instability to moist convection in saturated env
   ! integer mx_ssat(pcols)    ! level of instability to dry convection in subsaturated env
   ! integer mx_sat_2(pcols)   ! level of instability to moist convection in saturated env - higher region
   ! integer mx_ssat_2(pcols)  ! level of instability to dry convection in subsaturated env - higher region
!
!--------------------------Local Variables------------------------------
!
   ! real(r8) capeten(pcols,5)     ! provisional value of cape
   ! real(r8) capeten_2(pcols,5)   ! provisional value of cape, second region
   ! real(r8) tv(pcols,pver)       !
   ! real(r8) tpv(pcols,pver)      !
   ! real(r8) tpv_2(pcols,pver)
   ! real(r8) buoy(pcols,pver)
   ! real(r8) buoy_2(pcols,pver)

   ! real(r8) beta_mix_integral(pcols) !used to calculate virtual "potential" temperature
   real(r8) vir_pot_t(pcols,pver) !virtual "potential" temperature
   ! real(r8) vir_pot_t_max(pcols)  !maximum in each model column

   ! real(r8) grad_t
   ! real(r8) moist_adiabat
   ! real(r8) amu
   ! real(r8) gammas
   ! real(r8) L
   ! real(r8) beta_mix
   ! real(r8) cp_mix
   real(r8) conv_stab(pcols, pver)
   ! real(r8) conv_integral(pcols)
   ! real(r8) conv_integral_max(pcols)
   real(r8) tdifr  !fractional temperature difference between layers
   real(r8) that(pcols,pver)  !upper grid interface temperatures
   real(r8) eta(pcols,pver) !entropy, for diagnostic purposes

   ! real(r8) a1(pcols)
   ! real(r8) a2(pcols)
   ! real(r8) estp(pcols)
   real(r8) pl(pcols)
   ! real(r8) plexp(pcols)
   ! real(r8) hmax(pcols)
   ! real(r8) hmn(pcols)
   ! real(r8) y(pcols)

   ! logical plge600(pcols)
   logical mx1_fixed(pcols)
   logical mx2_fixed(pcols)  !upper convective region
   ! integer knt(pcols)
   ! integer lelten(pcols,5)
   integer max_mx1(pcols)
   ! integer knt_2(pcols)
   ! integer lelten_2(pcols,5)
   integer max_mx2(pcols)

   real(r8) cp
   real(r8) e
   real(r8) esat
   real(r8) qsat(pcols,pver)
   real(r8) grav
   real(r8) zvirp1
  
   integer i
   integer k
   integer msg
   integer n

   real(r8) rd
   real(r8) rl
#ifdef PERGRO
   real(r8) rhd
#endif

   real(r8) this_lat, this_lon

   !test values for entropy and ientropy
   real(r8) pdummy, tdummy, qdummy, tfguess
   real(r8) tfin, qsatfin, qvfin, qlfin
   real(r8) sdummy, s1, s2, s3, s4
   integer rcall
!
!-----------------------------------------------------------------------
!
   !test section for entropy and ientropy
   !pdummy = 1000._r8  !either mb or Pa
   !tdummy = 300._r8  !K
   !qdummy = 1.e-4_r8  !kg/kg
   !sdummy = entropy(tdummy, pdummy, qdummy)
   !tfguess = tdummy
   !rcall = 10._r8
   !call ientropy(rcall, 1, lchnk, sdummy, pdummy, qdummy, tfin,qsatfin, tfguess, qvfin, qlfin)
   !write(iulog,*) "*** (I)ENTROPY TEST. p: ", pdummy, "t: ", tdummy, "q: ", qdummy, "s: ", sdummy, "t_out: ", tfin, "qsat_out", qsatfin
   !s1 = entropy(tdummy-10._r8, pdummy, qdummy)
   !s2 = entropy(tdummy-5._r8, pdummy, qdummy)
   !s3 = entropy(tdummy+5._r8, pdummy, qdummy)
   !s4 = entropy(tdummy+10._r8, pdummy, qdummy)
   !write(iulog,*) "*** ENTROPY TEST. s1: ", s1, "s2: ", s2, "sdummy: ", sdummy, "s3: ", s3, "s4: ", s4

   zvirp1 = zvir + 1.0

   do n = 1,5
      do i = 1,ncol
         ! lelten(i,n) = pver
         ! capeten(i,n) = 0._r8
         ! lelten_2(i,n) = pver
         ! capeten_2(i,n) = 0._r8
      end do
   end do
!
   do i = 1,ncol
      lon(i) = pver
      ! knt(i) = 0
      ! knt_2(i) = 0
      lel(i) = 2
      lel_2(i) = 2
      ! mx_sat(i) = 2
      ! mx_sat_2(i) = 2
      ! mx_ssat(i) = 2
      ! mx_ssat_2(i) = 2
      mx(i) = 2  !all set to second-highest model level as default in case no suitable start point exists
      mx_2(i) = 2
      cape(i) = 0._r8
      cape_2(i) = 0._r8
      ! hmax(i) = 0._r8
      ! beta_mix_integral(i) = 0._r8
      ! vir_pot_t_max(i)=0._r8
      ! conv_integral_max(i)=0._r8
      max_mx1(i) = pver
      max_mx2(i) = pver
      mx1_fixed(i) = .false.  
      mx2_fixed(i) = .false.  
   end do
!
! Calculate some variables including upper interface temperatures
   do k = msg + 1,pver
      do i = 1,ncol
         tdifr = 0._r8
         tdifr = abs((t(i,k)-t(i,k-1))/max(t(i,k-1),t(i,k)))
         if (tdifr > 1.E-6_r8) then
            that(i,k) = log(t(i,k-1)/t(i,k))*t(i,k-1)*t(i,k)/(t(i,k-1)-t(i,k))
         else
            that(i,k) = 0.5_r8* (t(i,k)+t(i,k-1))
         end if
         call entropy(t(i,k), p(i,k), q(i,k), eta(i,k))
         ! tv(i,k) = t(i,k) * (1 - mu_red*q(i,k))
         ! tpv(i,k) = tv(i,k)
         ! buoy(i,k) = 0._r8
         ! vir_pot_t(i,k) = 0._r8
         ! ! conv_stab(i,k) = 0._r8
         ! tp(i,k) = t(i,k)
         ! qstp(i,k) = q(i,k)

         !find q_sat and whether it is close to q, this lets us decide whether to look for a max in C or theta_v
         call esat_buck(t(i,k), esat)
         e = q(i,k)*p(i,k)/(q(i,k)+eps1*(1-q(i,k)))
         qsat(i,k) = eps1*esat/(p(i,k)-e*(1-eps1))

      end do
   end do

!find vir_pot_t and conv_stab to let us find convective unstable start points
   call layer_stabilities(lchnk, ncol, msg, t, that, sfcT, q, p, pf, vir_pot_t, conv_stab)

   !work through possible start points to find start point and CAPE of first convective region
   do while (.not. all(mx1_fixed(:ncol))) 
      call find_mx_CAPE(lchnk, ncol, msg, t, q, qsat, p, pf, tpert, vir_pot_t, conv_stab, mx1_fixed, max_mx1, mx, lel, cape, tp, qstp)
      ! write(iulog,*) "mx1_fixed:", mx1_fixed, "max_mx1:", max_mx1
   end do

   ! do same but with second convective region
   do i=1,ncol
      max_mx2(i) = lel(i) - 1
   end do
   do while (.not. all(mx2_fixed(:ncol))) 
      call find_mx_CAPE(lchnk, ncol, msg, t, q, qsat, p, pf, tpert, vir_pot_t, conv_stab, mx2_fixed, max_mx2, mx_2, lel_2, cape_2, tp_2, qstp_2)
      ! write(iulog,*) "mx2_fixed:", mx2_fixed, "max_mx2:", max_mx2
   end do



   ! do i = 1,ncol
   !    this_lat = get_rlat_p(lchnk,i)*57.296_r8
   !    this_lon = get_rlon_p(lchnk,i)*57.296_r8 
   !    if (cape(i)>70._r8 .and. ((this_lat < -22._r8 .and. this_lat > -23._r8) .and. (this_lon > 154._r8 .and. this_lon < 161._r8))) then
   !       write(iulog,*) " "
   !       write(iulog,*) " "
   !       write(iulog,*) "*** ZM_CONV: Buoyan_dilute. At lchnk: ", lchnk, "i: ", i, " lat: ", this_lat, "lon: ", this_lon, "calling parcel_dilute."
   !       write(iulog,*) "*** ZM_CONV: Buoyan_dilute. mx: ", mx(i), "CAPE: ", cape(i), "lel:", lel(i), "pconvb:", p(i,mx(i)), "pconvt:", p(i,lel(i)) !"mx_2: ", mx_2(i), "CAPE_2: ", cape_2(i), "lel_2:", lel_2(i)
   !       ! write(iulog,*) "*** ZM_CONV: Buoyan_dilute. mx: ", mx(i), "mx_sat: ", mx_sat(i), "mx_ssat: ", mx_ssat(i), "CAPE: ", cape(i), "lel:", lel(i)
   !       ! write(iulog,*) "*** ZM_CONV: Buoyan_dilute. mx_2: ", mx_2(i), "mx_sat_2: ", mx_sat_2(i), "mx_ssat_2: ", mx_ssat_2(i), "CAPE_2: ", cape_2(i), "lel_2:", lel_2(i)
   !       ! write(iulog,*) "*** ZM_CONV: Parcel dilute. lcl: ", lcl(i), "plcl: ", pl(i), "tlcl: ", tl(i)
   !       write(iulog,*) " "
   !       ! do k = mx(i)+2,lel(i)-2,-1
   !       do k = pver, lel(i)-2, -1
   !          write(iulog,*) "*** Going up column, LCR.  T: ", t(i,k), "P: ", p(i,k), "Q: ", q(i,k), "k:", k  !"qsat: ", qsat(i,k), 
   !          write(iulog,*) "*** Finding mx. virtual potential temperature: ", vir_pot_t(i,k), "moist conv stability: ", conv_stab(i,k)
   !          ! write(iulog,*) "*** Finding mx. virtual potential temperature: ", vir_pot_t(i,k), "entropy: ", eta(i,k)
   !          ! write(iulog,*) "*** Parcel_dilute. env tv:", tv(i,k), "parcel T: ", tp(i,k), "parcel tv: ", tpv(i,k), "cape contribution:", rd*buoy(i,k)*log(pf(i,k+1)/pf(i,k))  !, "parcel q: ", qstp(i,k)
   !          ! write(iulog,*) "*** Going up column. dp:", p(i,k)-p(i,k+1), "dz:", z(i,k)-z(i,k+1)
   !          write(iulog,*) " "
   !       end do
   !       ! do k = mx_2(i)+2,lel_2(i)-2,-1
   !       !    write(iulog,*) "*** Going up column, UCR.  T: ", t(i,k), "P: ", p(i,k), "Q: ", q(i,k), "env tv:", tv(i,k), "k:", k  !"qsat: ", qsat(i,k), 
   !       !    !write(iulog,*) "*** Finding mx. virtual potential temperature: ", vir_pot_t(i,k), "moist conv stability: ", conv_stab(i,k)
   !       !    write(iulog,*) "*** Finding mx. virtual potential temperature: ", vir_pot_t(i,k), "entropy: ", eta(i,k)
   !       !    write(iulog,*) "*** Parcel_dilute. env tv:", tv(i,k), "parcel T: ", tp_2(i,k), "parcel tv: ", tpv_2(i,k), "cape contribution:", rd*buoy_2(i,k)*log(pf(i,k+1)/pf(i,k))  !, "parcel q: ", qstp(i,k)
   !       !    !write(iulog,*) "*** Going up column. dp:", p(i,k)-p(i,k+1), "dz:", z(i,k)-z(i,k+1)
   !       !    write(iulog,*) " "
   !       ! end do
   !    end if
   ! end do

   return
end subroutine buoyan_dilute

subroutine layer_stabilities(lchnk, ncol, msg, t, that, sfcT, q, p, pf, vir_pot_t, conv_stab)

   ! Routine to determine 1. the virtual potential temperature at each layer, 
   ! and 2. the stability to moist convection of each layer

   integer, intent(in) :: lchnk
   integer, intent(in) :: ncol
   integer, intent(in) :: msg         ! top level of allowed convection

   real(r8), intent(in), dimension(pcols,pver) :: p      !midpoint pressure
   real(r8), intent(in), dimension(pcols,pver) :: t      !midpoint temperature
   real(r8), intent(in), dimension(pcols,pver) :: q      !midpoint water mmr
   real(r8), intent(in), dimension(pcols,pver) :: that   !upper interface temperature
   real(r8), intent(in), dimension(pcols,pver+1) :: pf   !upper interface pressure
   real(r8), intent(in), dimension(pcols) :: sfcT        !surface temperature

   real(r8), intent(out), dimension(pcols,pver) :: vir_pot_t   !virtual potential temperature at midpoint
   real(r8), intent(out), dimension(pcols,pver) :: conv_stab   !moist convective stability across layer

   real(r8) cp_mix
   real(r8) amu
   real(r8) gammas
   real(r8) grad_t
   real(r8) moist_adiabat
   real(r8) beta_mix
   real(r8) beta_mix_integral(pcols)
   real(r8) L
   integer i, k

   do i=1,ncol
      beta_mix_integral(i) = 0._r8
   end do
   do k = pver, msg+1, -1
      do i = 1,ncol
         !calculate some general quantities
         cp_mix = ((1-q(i,k))*cpres + q(i,k)*cpwv)
         beta_mix = (rgas*(1 - mu_red*q(i,k))/cp_mix)
         L = rl - (cpliq - cpwv)*(t(i,k)-tfreez)  !not sure if this is needed or if we can just use rl
         amu = mu_red*q(i,k) / (1-mu_red*q(i,k))
         gammas = (1-mu_red*q(i,k)) * L/(t(i,k)*rh2o)

         !calculate stability if subsaturated
         if (k==pver) then
            beta_mix_integral(i) = beta_mix_integral(i) + beta_mix * log(pf(i,k+1)/p(i,k)) !need to start somewhere so start at very model bottom
         else
            beta_mix_integral(i) = beta_mix_integral(i) + beta_mix * log(p(i,k+1)/p(i,k))
         end if
         vir_pot_t(i,k) = t(i,k) * (1 - mu_red*q(i,k)) * exp(beta_mix_integral(i))
         !calculate stability if saturated
         if (k == pver) then !that does not do surface T so use sfcT
            grad_t = p(i,k)/t(i,k) * (that(i,k)-sfcT(i))/(pf(i,k)-pf(i,k+1))
         else     
            grad_t = p(i,k)/t(i,k) * (that(i,k)-that(i,k+1))/(pf(i,k)-pf(i,k+1))
         end if
         moist_adiabat = beta_mix * (1+ q(i,k)/(1-q(i,k)) * L/(t(i,k)*rgas)) / (1+ q(i,k)/(1-q(i,k)) * gammas*L/(cp_mix*t(i,k)))
         conv_stab(i,k) = (grad_t - moist_adiabat) * (amu*gammas - 1)

         ! if ((i==1 .and. lchnk==33) .and. (k<52 .and. k>28)) then
         !    write(iulog,*) " "
         !    ! if (k==pver) then
         !    !    write(iulog,*) "*** ZM_CONV: Ts:", sfcT(i)
         !    ! end if
         !    ! write(iulog,*) "*** ZM_CONV: Qsat check. T: ", t(i,k), "P: ", p(i,k), "Q: ", q(i,k), "eps1: ", eps1, "esat: ", esat, "qsat: ", qsat(i,k)
         !    ! write(iulog,*) "*** ZM_CONV: grad_t test. T", t(i,k), "P:", p(i,k), "That:", that(i,k), "Pf:", pf(i,k)
         !    ! write(iulog,*) "*** ZM_CONV: grad_t test. grad_t:", grad_t, "moist_adiabat:", moist_adiabat
         !    ! write(iulog,*) "*** ZM_CONV: moist adiabat. amu:", amu, "gammas:", gammas
         !    write(iulog,*) "*** ZM_CONV: Finding mx. k: ", k, "T(K): ", t(i,k), "Q(kg/kg): ", q(i,k)!, "Qsat(kg/kg): ", qsat(i,k)
         !    write(iulog,*) "*** ZM_CONV: Finding mx. virtual potential temperature: ", vir_pot_t(i,k), "moist conv stability: ", conv_stab(i,k)
         !    ! ! write(iulog,*) "*** VIR_POT_T:", vir_pot_t(i,k), "T:", t(i,k), "Q:", q(i,k), "P:", p(i,k)
         !    ! write(iulog,*) "*** VIR_POT_T. beta_mix:", beta_mix, "Pf(k+1):", pf(i,k+1), "Pf(k):", pf(i,k), "integral:", beta_mix_integral(i)
         !    ! write(iulog,*) "beta_mix:", beta_mix, "cpmix:", cp_mix, "rd:", rgas
         ! end if         
      end do
   end do

   return
end subroutine layer_stabilities

subroutine find_mx_CAPE(lchnk, ncol, msg, t, q, qsat, p, pf, tpert, vir_pot_t, conv_stab, mx_fixed, max_next_mx, mx, lel, CAPE, tp, qstp)

   !using calculated layer stabilities and a maximum layer index max_mx, look to see for a possible convective zone
   integer, intent(in) :: lchnk
   integer, intent(in) :: ncol
   integer, intent(in) :: msg         ! top level of allowed convection

   real(r8), intent(in), dimension(pcols,pver) :: p      !midpoint pressure
   real(r8), intent(in), dimension(pcols,pver) :: t      !midpoint temperature
   real(r8), intent(in), dimension(pcols,pver) :: q      !midpoint water mmr
   real(r8), intent(in), dimension(pcols,pver) :: qsat   !midpoint saturation water mmr
   real(r8), intent(in), dimension(pcols,pver+1) :: pf   !upper interface pressure
   real(r8), intent(in), dimension(pcols) :: tpert       !PBL pertubation temperature
   real(r8), intent(in), dimension(pcols,pver) :: vir_pot_t   !virtual potential temperature at midpoint
   real(r8), intent(in), dimension(pcols,pver) :: conv_stab   !moist convective stability across layer

   logical, intent(inout), dimension(pcols) :: mx_fixed     ! whether we have settled on a starting point for this convection
   integer, intent(inout), dimension(pcols) :: max_next_mx  ! what the maximum (lowest) conv starting we will consider is going to be
   integer, intent(inout), dimension(pcols) :: mx
   integer, intent(inout), dimension(pcols) :: lel
   real(r8), intent(inout), dimension(pcols) :: CAPE

   real(r8), intent(out), dimension(pcols,pver) :: tp !parcel temperature
   real(r8), intent(out), dimension(pcols,pver) :: qstp !parcel qsat

   integer mx_sat(pcols)
   integer mx_ssat(pcols)
   integer lcl(pcols)
   real(r8) tl(pcols)
   real(r8) pl(pcols)

   logical plume_definitely_ended(pcols)
   integer knt(pcols)
   integer lelten(pcols,5)
   real(r8) capeten(pcols,5)     ! provisional value of cape
   real(r8) tv(pcols,pver)       !
   real(r8) tpv(pcols,pver)      !
   real(r8) buoy(pcols,pver)

   integer i, k, n

   !initialise relevant variables
   do i=1,ncol
      mx_sat(i) = 2
      mx_ssat(i) = 2
      lcl(i) = 2
      plume_definitely_ended(i) = .false.
      knt(i) = 0
   end do
   do n=1,5
      do i=1,ncol
         lelten(i,n) = 2
         capeten(i,n) = 0._r8
      end do
   end do
   do k=msg,pver
      do i = 1,ncol
         if (.not. mx_fixed(i)) then !we have not yet fixed a convective start point, and so set (and reset) the affected tp, qstp to the layer wvalues
            tp(i,k) = t(i,k)
            qstp(i,k) = q(i,k)
         end if
      end do          
   end do


   ! Loop through columns again and find mx_sat and mx_ssat
   do k = msg+2, pver
      do i = 1, ncol
         if (.not. mx_fixed(i) .and. max_next_mx(i)>=k) then !this model column is one where we are looking for our upper/lower convective zone, and we are looking high enough in the column
            if ((conv_stab(i,k) < 0._r8  .and.  mx_sat(i) < k) .and. (q(i,k)/qsat(i,k) > 0.4_r8)) then  !take lowest level unstable to moist convection, regardless of whether system is saturated - need to remove subsaturation requirement
               mx_sat(i) = k
            end if

            if ((qsat(i,k)-q(i,k))/(qsat(i,k)+q(i,k)) > 1.e-3_r8) then !UNSATURATED. too stringent q requirement? maybe drop to 0.01?

               if (k==pver) then  !on bottom level so cannot check vir_pot_t on level below
                  if (vir_pot_t(i,k) > vir_pot_t(i,k-1)) then  !maximum in vir_pot_t
                     mx_ssat(i) = k
                  end if
               else  !not on bottom model level
                  if (vir_pot_t(i,k) > vir_pot_t(i,k-1) .and. vir_pot_t(i,k) > vir_pot_t(i,k+1)) then  !maximum in vir_pot_t
                     mx_ssat(i) = k
                  end if
               end if
            end if
         end if
      end do
   end do

   do i = 1,ncol ! Initialise LCL variables and find correct mx
      if (.not. mx_fixed(i)) then
         mx(i) = max(mx_sat(i),mx_ssat(i))
         lcl(i) = mx(i)
         tl(i) = t(i,mx(i))
         pl(i) = p(i,mx(i))
      end if
   end do
!
! main buoyancy calculation.
!
   call parcel_dilute(lchnk, ncol, msg, mx, p, t, q, tpert, tp, tpv, qstp, pl, tl, lcl, mx_fixed)
!
   do k = pver,msg + 1,-1
      do i=1,ncol
         if (.not. mx_fixed(i)) then
            if (k <= mx(i)) then   ! Define buoy from launch level to cloud top.
               tv(i,k) = t(i,k) * (1 - mu_red*q(i,k))
               buoy(i,k) = tpv(i,k) - tv(i,k) + tiedke_add  ! +0.5K or not?
            else
               qstp(i,k) = q(i,k)
               tp(i,k)   = t(i,k)            
               tpv(i,k)  = tv(i,k)
            endif
         end if
      end do
   end do
!

! decide where the plume top is by considering changes in buoyancy
   do k = pver, msg+2,-1 !working up column
      do i = 1,ncol
         if (.not. mx_fixed(i)) then
            if (k < mx(i)-2) then
               if ((buoy(i,k+1) > 0. .and. buoy(i,k) <= 0._r8) .and.  .not. plume_definitely_ended(i)) then
                  knt(i) = min(5,knt(i) + 1)
                  lelten(i,knt(i)) = k
                  if (buoy(i,k) < -0.5_r8) plume_definitely_ended(i) = .true.
               end if
               if (knt(i)>=2 .and. buoy(i,k) <= -0.5_r8) plume_definitely_ended(i) = .true. !only trigger once one plume top has already been found
            end if
         end if
      end do
   end do
!
! calculate convective available potential energy (cape).
!
   do n = 1,5
      do k = msg + 1,pver
         do i = 1,ncol
            if (.not. mx_fixed(i)) then
               if (k <= mx(i) .and. k > lelten(i,n)) then
                  capeten(i,n) = capeten(i,n) + rgas*buoy(i,k)*log(pf(i,k+1)/pf(i,k))
               end if
            end if
         end do
      end do
   end do
!
! find maximum cape from all possible tentative capes from one sounding, and use it as the final cape, april 26, 1995
!
   do n = 1,5
      do i = 1,ncol
         if (.not. mx_fixed(i)) then
            if (capeten(i,n) > cape(i)) then
               cape(i) = capeten(i,n)
               lel(i) = lelten(i,n)
            end if
         end if
      end do
   end do
!
! put lower bound on cape for diagnostic purposes.
!
   do i = 1,ncol
      cape(i) = max(cape(i), 0._r8)
   end do
!
! set variables for where to keep looking for convection
!
   do i = 1,ncol
      if (.not. mx_fixed(i)) then
         !new mx has been found with non-zero CAPE
         if (cape(i)>0._r8) then
            mx_fixed(i) = .true.
            max_next_mx(i) = lel(i)-1
         end if
         !if CAPE is negative - but mx was set - keep looking in case there's a better mx higher up
         if (cape(i)==0._r8 .and. mx(i)>2) then
            mx_fixed(i) = .false. !is already false but restate here for clarity
            max_next_mx(i) = mx(i)-1
         end if
         !if mx was not found / left at 2 - stop looking
         if (cape(i)==0._r8 .and. mx(i)==2) then
            mx_fixed(i) = .true.
            max_next_mx(i) = 2
         end if
      end if
   end do

   ! if (lchnk==33) then
   !    ! write(iulog,*) "conv_stab:", conv_stab(1,:), "vir_pot_t:", vir_pot_t(1,:)
   !    write(iulog,*) "mx_fixed:", mx_fixed, "max_next_mx:", max_next_mx, "cape:", cape, "mx:", mx
   ! end if
   !ideally - nothing at all in here would be affected if mx_fixed is already true - we already have our CAPE, mx, lel and shouldn't be looking for more
   !key question - does the fashion we do plume_definition_ended affect anything? should it?
   !in reality it should always be true once we have worked through the whole column, it only matters for allowing us to look a bit further for a new lel
   !instead - do we want to impose a new constraint on convection that it must find buoyancy within a given number of layers / a given fraction of a scale height?
   return
end subroutine find_mx_CAPE

subroutine parcel_dilute (lchnk, ncol, msg, klaunch, p, t, q, tpert, tp, tpv, qstp, pl, tl, lcl, leave_column)

! Routine  to determine 
!   1. Tp   - Parcel temperature
!   2. qstp - Saturated mixing ratio at the parcel temperature.

!--------------------
implicit none
!--------------------

integer, intent(in) :: lchnk
integer, intent(in) :: ncol
integer, intent(in) :: msg

integer, intent(in), dimension(pcols) :: klaunch(pcols)
real(r8), intent(in), dimension(pcols,pver) :: p
real(r8), intent(in), dimension(pcols,pver) :: t
real(r8), intent(in), dimension(pcols,pver) :: q
real(r8), intent(in), dimension(pcols) :: tpert ! PBL temperature perturbation.
logical, intent(in), dimension(pcols) :: leave_column !true: don't do column, false: do it 

real(r8), intent(inout), dimension(pcols,pver) :: tp    ! Parcel temp.
real(r8), intent(inout), dimension(pcols,pver) :: qstp  ! Parcel water vapour (sat value above lcl).
real(r8), intent(inout), dimension(pcols) :: tl         ! Actual temp of LCL.
real(r8), intent(inout), dimension(pcols) :: pl          ! Actual pressure of LCL. 

integer, intent(inout), dimension(pcols) :: lcl ! Lifting condesation level (first model level with saturation).

real(r8), intent(out), dimension(pcols,pver) :: tpv   ! Virtual temperature of parcel

!--------------------

! Have to be careful as s is also dry static energy.

! If we are to retain the fact that CAM loops over grid-points in the internal
! loop then we need to dimension sp,atp,mp,xsh2o with ncol.


real(r8) tmix(pcols,pver)        ! Temperature of the entraining parcel.
real(r8) qtmix(pcols,pver)       ! Total water of the entraining parcel.
real(r8) qsmix(pcols,pver)       ! Saturated mixing ratio at the tmix.
real(r8) smix(pcols,pver)        ! Entropy of the entraining parcel.
real(r8) xsh2o(pcols,pver)       ! Precipitate lost from parcel.
real(r8) ds_xsh2o(pcols,pver)    ! Entropy change due to loss of condensate.
real(r8) ds_freeze(pcols,pver)   ! Entropy change sue to freezing of precip.

real(r8) mp(pcols)    ! Parcel mass flux.
real(r8) qtp(pcols)   ! Parcel total water.
real(r8) sp(pcols)    ! Parcel entropy.

real(r8) sp0(pcols)    ! Parcel launch entropy.
real(r8) qtp0(pcols)   ! Parcel launch total water.
real(r8) mp0(pcols)    ! Parcel launch relative mass flux.
real(r8) p0(pcols)     ! Parcel launch pressure
real(r8) t0(pcols)     ! Parcel launch temperature
real(r8) tlexp(pcols)  ! dry adiabat for our column

real(r8) lwmax      ! Maximum condesate that can be held in cloud before rainout.
real(r8) dmpdp      ! Parcel fractional mass entrainment rate (/mb).
!real(r8) dmpdpc     ! In cloud parcel mass entrainment rate (/mb).
real(r8) dmpdz(pcols)      ! Parcel fractional mass entrainment rate (/m)
real(r8) dpdz,dzdp  ! Hydrstatic relation and inverse of.
real(r8) senv       ! Environmental entropy at each grid point.
real(r8) qtenv      ! Environmental total water "   "   ".
real(r8) penv       ! Environmental total pressure "   "   ".
real(r8) tenv       ! Environmental total temperature "   "   ".
real(r8) new_s      ! Hold value for entropy after condensation/freezing adjustments.
real(r8) new_q      ! Hold value for total water after condensation/freezing adjustments.
real(r8) dp         ! Layer thickness (center to center)
real(r8) tfguess    ! First guess for entropy inversion - crucial for efficiency!
real(r8) tscool     ! Super cooled temperature offset (in degC) (eg -35).

real(r8) qxsk, qxskp1               ! LCL excess water (k, k+1)
real(r8) dsdp, dqtdp, dqxsdp        ! LCL s, qt, p gradients (k, k+1)
real(r8) slcl,qtlcl,qslcl,tlcl      ! LCL s, qt, qs , T values.
real(r8) zvirp1
real(r8) esat, qsat, e              ! saturation vapour pressure and mass mixing ratio, current vapour pressure
real(r8) mass_scaling               ! amount qv, ql need to be scaled by to account for precipitation losses
real(r8) total_mass                 ! parcel mass including entrained component (but before accounting for precipitation losses)
real(r8) qv(pcols,pver)             ! water vapour of entraining parcel
real(r8) ql(pcols,pver)             ! condensed water of entraining parcel

integer rcall       ! Number of ientropy call for errors recording
integer i,k,ii   ! Loop counters.
logical lcl_reached !has the lcl been reached on the ascent


!======================================================================
!    SUMMARY
!
!  9/9/04 - Assumes parcel is initiated from level of maxh (klaunch)
!           and entrains at each level with a specified entrainment rate.
!
! 15/9/04 - Calculates lcl(i) based on k where qsmix is first < qtmix.          
!
!======================================================================
!
! Set some values that may be changed frequently.
!

zvirp1 = zvir + 1 !MWDAIR/MWWV (~0.3 for H2 dominated, 1.608 on earth)
dmpdz=-1.e-3_r8        ! Entrainment rate. (-ve for /m), default 10^-3
!dmpdpc = 3.e-2_r8   ! In cloud entrainment rate (/mb).
lwmax = 1.e-3_r8    ! Need to put formula in for this. Used to be fixed at 1.e-3
tscool = 0.0_r8   ! Temp at which water loading freezes in the cloud.

qtmix=0._r8
smix=0._r8

qtenv = 0._r8
senv = 0._r8
tenv = 0._r8
penv = 0._r8

qtp0 = 0._r8
sp0  = 0._r8
mp0 = 0._r8

qtp = 0._r8
sp = 0._r8
mp = 0._r8

new_q = 0._r8
new_s = 0._r8

lcl_reached= .false.

do i=1,ncol
   dmpdz(i)=-7.e-6_r8 * shr_const_mwdair/(1-mu_red*q(i,klaunch(i)))      ! Entrainment rate. (-ve for /m), default 10^-3
end do

do k=pver, msg+1, -1 !going up column
    do i=1,ncol
      if (.not. leave_column(i)) then
         if (k == klaunch(i)) then
            qtp0(i) = q(i,k)  !assume all in vapour form, don't include existing liquid water here
            call entropy(t(i,k), p(i,k), qtp0(i), sp0(i))
            p0(i)=p(i,k)
            t0(i)=t(i,k)
            mp0(i) = 1
            mp(i) = 1 !initialise mp
            smix(i,k) = sp0(i) !starting entropy
            qtmix(i,k) = qtp0(i) !total water content
            tmix(i,k) = t(i,k)
            tpv(i,k) = (1 - mu_red*qtmix(i,k))*(tmix(i,k)+tpert(i)) !virtual temperature is ultimately what we are really after for

            esat = c1*exp(c2*(t0(i)-tfreez)/(c3+t0(i)-tfreez))
            call esat_buck(t0(i), esat)

            e = q(i,k)*p(i,k)/(q(i,k)+eps1*(1-q(i,k)))
            qsmix(i,k)=eps1*esat/(p0(i)-e*(1-eps1))!do not need to run whole ientropy just to get this
            
            ! if (i==5 .and. lchnk==33) then
            !    write(iulog,*) " "
            !    write(iulog,*) " "
            !    write(iulog,*) "*** ZM_CONV: parcel_dilute. lchnk: ", lchnk, "i: ", i, "mx: ",k
            !    write(iulog,*) "*** ZM_CONV: launch level. tmix:", tmix(i,k), "smix:", smix(i,k), "qtmix:", qtmix(i,k), "qsmix:", qsmix(i,k)
            ! end if
            
         end if !just the launch level

         if ((k<klaunch(i)) ) then !all the other levels regardless of above/below LCL

            !find properties of entrained environment
            dp = p(i,k) - p(i,k+1)
            qtenv = 0.5_r8 * (q(i,k) + q(i,k+1))
            tenv = 0.5_r8 * (t(i,k) + t(i,k+1))
            penv = 0.5_r8 * (p(i,k) + p(i,k+1))
            call entropy(tenv, penv, qtenv, senv)

            !find entrainment rates
            dpdz = -(penv*grav)/(rgas*tenv)
            dzdp = 1._r8/dpdz
            dmpdp = dmpdz(i)*dzdp !/mb fractional entrainment
            !write(iulog,*) " FRACTIONAL ENTRAINMENT: ", dmpdz*dzdp

            total_mass = mp(i) - dmpdp*dp !as the latter quantity is negative, this is >mp(i)
            smix(i,k) = (smix(i,k+1)*mp(i) - senv * dmpdp*dp) /total_mass !parcel entropy before precipitation adjustments
            qtmix(i,k) = (qtmix(i,k+1)*mp(i) - qtenv * dmpdp*dp) /total_mass !parcel water before precipitation adjustments
            mp(i) = total_mass

            !smix(i,k) = smix(i,k+1) !!!TEST
            !qtmix(i,k) = qtmix(i,k+1) !!!TEST

            rcall = 1._r8
            tfguess = t(i,k)

            ! if (i==5 .and. lchnk==33 .and. k>(klaunch(i) - 5)) then
            !    write(iulog,*) " "
            !    write(iulog,*) "*** ZM_CONV: parcel_dilute. lchnk: ", lchnk, "i: ", i, "k: ",k
            !    write(iulog,*) "*** ZM_CONV: qtenv: ", qtenv, "tenv: ", tenv, "penv: ", penv, "senv:", senv
            !    write(iulog,*) '*** ZM_CONV: fractional entrainment: ', dmpdp, "dp: ", dp, "mp: ", mp(i)
            !    write(iulog,*) "*** ZM_CONV: previous entropy: ", smix(i,k+1), "previous q:", qtmix(i,k+1)
            !    write(iulog,*) "*** ZM_CONV: ientropy inputs: smix: ", smix(i,k), "qtmix: ", qtmix(i,k), "p:", p(i,k)
            ! end if

            
            call ientropy(rcall, i, lchnk, smix(i,k), p(i,k), qtmix(i,k), tmix(i,k),qsmix(i,k), tfguess, qv(i,k), ql(i,k))
            !smix(i,k+1) used to account for the change in entropy from the precipitation
            !do we need to do a different one for the initial LCL level - if entropy is constant then maybe we don't

            xsh2o(i,k) = max(0._r8, ql(i,k) - lwmax) !new liquid water in cloud
            mass_scaling = 1 / (1-xsh2o(i,k)) !amount by which the mass mixing ratios need to be increased by
            ql(i,k) = (ql(i,k) - xsh2o(i,k)) * mass_scaling
            qv(i,k) = qv(i,k) * mass_scaling
            qtmix(i,k) = qv(i,k) + ql(i,k) !update parcel water to take into account precipitation losses
            mp(i) = mp(i) * (1 - xsh2o(i,k))

            call entropy(tmix(i,k), p(i,k), qtmix(i,k), smix(i,k)) !find updated specific entropy with new mixing ratios - should be slightly higher
            tpv(i,k) = (1 - qtmix(i,k) + zvirp1*qv(i,k))*(tmix(i,k)+tpert(i)) !virtual temperature is ultimately what we are really after for
            
            ! if (i==5 .and. lchnk==33 .and. k>(klaunch(i) - 5)) then
            !    write(iulog,*) " "
            !    write(iulog,*) "*** ZM_CONV: parcel_dilute. lchnk: ", lchnk, "i: ", i, "k: ",k
            !    write(iulog,*) '*** ZM_CONV: parcel T: ', tmix(i,k), "updated qtmix:", qtmix(i,k), 'updated smix: ', smix(i,k), "tfguess:", tfguess
            !    write(iulog,*) '*** ZM_CONV: parcel Qsat: ', qsmix(i,k), 'mass_scaling: ',mass_scaling,'mp: ',mp(i), "tpv:", tpv(i,k)
            ! end if


            if ( qsmix(i,k)<qtmix(i,k) .and. qsmix(i,k+1)>qtmix(i,k+1) ) then
                  lcl(i) = k

                  dp = (p(i,k) - p(i,k+1)) !In -ve mb as p decreasing with height, difference between centre of layers
                  qxsk = qtmix(i,k) - qsmix(i,k) !supersaturation in this level
                  qxskp1 = qtmix(i,k+1) - qsmix(i,k+1) !undersaturation in level below
                  dqxsdp = (qxsk - qxskp1)/dp !change in supersaturation with pressure
                  pl(i) = p(i,k+1) - qxskp1/dqxsdp !exact pressure where supersaturation is reached
                  slcl=sp0(i) !assuming no entrainment - looking for max CAPE after all
                  qtlcl=qtp0(i) !assuming no entrainment
                  tlcl = tmix(i,k+1) * (pl(i)/p(i,k+1))**tlexp(i)
                  tl(i) = tlcl

                  !if (klaunch(i)==40) then
                  !   write(iulog,*) " "
                  !   write(iulog,*) "*** ZM_CONV: parcel_dilute. lchnk: ", lchnk, "i: ", i, "mx: ",k
                  !   write(iulog,*) '*** ZM_CONV: At LCL. k: ', k, "parcel Q: ", qtp0(i), "parcel Qsat: ", qsat
                  !end if
            end if

         end if
         tp(i,k)=tmix(i,k) !doesn't give values below parcel launch
         qstp(i,k)=qtmix(i,k)
      end if
   end do
end do
!now we have T and q (assume ql is small or precipitates), we can calculate virtual temperature and so the buoyancy at each level
!as in the buoyan_dilute calculation




return
end subroutine parcel_dilute


!-----------------------------------------------------------------------------------------
SUBROUTINE entropy(TK,p,qtot, eta)
!-----------------------------------------------------------------------------------------
!
! TK(K),p(mb),qtot(kg/kg)
! from Raymond and Blyth 1992
!
     real(r8), intent(in) :: p,qtot,TK
     real(r8), intent(out) :: eta
     real(r8) :: qv,qc,ql,qi,qsat,e,esat,L,Lice,eref,pref
     real(r8) :: T_min, T_max !maximum and minimum temperatures for water freezing vs condensing

pref = 1000.0_r8           ! mb
eref = 6.106_r8            ! sat p at tfreez (mb)

L = rl - (cpliq - cpwv)*(TK-tfreez)         ! T IN CENTIGRADE
Lice = rlice + (cpliq - cpice)*(Tk-tfreez)
T_min = tfreez - 10
T_max = tfreez

! Replace call to satmixutils.

esat = c1*exp(c2*(TK-tfreez)/(c3+TK-tfreez))       ! esat(T) in mb
call esat_buck(TK, esat)

e = qtot*p/(qtot+eps1*(1-qtot))                    ! Using e instead of esat here as well, may be justified
qsat=eps1*esat/(p-e*(1-eps1))                      ! Sat. mixing ratio (in kg/kg).

qv = min(qtot,qsat)                         ! Partition qtot into vapor part only.
qc = max(0._r8,(qtot-qsat))                     ! Total condensed water sum

!!$e = qv*p / (eps1*(1-qv) +qv)
!!$
!!$entropy = ((1-qtot)*cpres + qtot*cpliq)*log( TK/tfreez) - (1-qtot)*rgas*log( (p-e)/pref ) + &
!!$     L*qv/TK - qv*rh2o*log(e/p)

ql = qc * (Tk - T_min)/(T_max - T_min)       ! Liquid water section of condensed amount
qi = qc * (T_max - Tk)/(T_max - T_min)       ! Solid water section of condensed amount
ql = max(0._r8,ql)
ql = min(ql,qc)
qi = max(0._r8,qi)
qi = min(qi,qc)

e = qv*p / (eps1*(1-qv) +qv)
qsat=eps1*esat/(p-e*(1-eps1))                      ! Sat. mixing ratio (in kg/kg).

eta = ((1-qtot)*cpres + qtot*cpliq)*log( TK/tfreez) - (1-qtot)*rgas*log( (p-e)/pref ) + &
        L*qv/TK - Lice*qi/Tk - qv*rh2o*log(e/esat)
! 
end SUBROUTINE entropy

!-----------------------------------------------------------------------------------------
SUBROUTINE entropy_gas(TK,p,qtot, eta)
   !-----------------------------------------------------------------------------------------
   !
   ! TK(K),p(mb),qtot(kg/kg)
   ! from Raymond and Blyth 1992
   !
        real(r8), intent(in) :: p,qtot,TK
        real(r8), intent(out) :: eta
        real(r8) :: eref,pref, e
  
   pref = 1000.0_r8           ! mb
   eref = 6.106_r8            ! sat p at tfreez (mb)
   
   e = qtot*p / (eps1*(1-qtot) +qtot)
   eta = ((1-qtot)*cpres + qtot*cpwv)*log( TK/tfreez) - (1-qtot)*rgas*log( (p-e)/pref ) - qtot*rh2o*log(e/pref)
   ! 
end SUBROUTINE entropy_gas


!-----------------------------------------------------------------------------------------
SUBROUTINE esat_buck(T,esat)
   !-----------------------------------------------------------------------------------------
   !
   ! TK(K),esat (mb)
   ! from wikipedia (2024)
   !
        real(r8), intent(in) :: T
        real(r8), intent(out) :: esat

        esat = 6.1121_r8 * exp( (18.678_r8 - (T-273.15_r8)/234.5_r8) * ((T-273.15_r8))/(T-16.01_r8))
   ! 
end SUBROUTINE esat_buck

!
!-----------------------------------------------------------------------------------------
   SUBROUTINE ientropy (rcall,icol,lchnk,s,p,qt,T,qsat,Tfg, qv,qc)
!-----------------------------------------------------------------------------------------
!
! p(mb), Tfg/T(K), qt/qv(kg/kg), s(J/kg). 
! Inverts entropy, pressure and total water qt 
! for T and saturated vapor mixing ratio
! 

     use phys_grid, only: get_rlon_p, get_rlat_p

     integer, intent(in) :: icol, lchnk, rcall
     real(r8), intent(in)  :: s, p, Tfg, qt
     real(r8), intent(out) :: qsat, T, qv, qc
     real(r8) :: Ts,dTs,fs1,fs2,esat,t2,t1,ql,qi  
     real(r8) :: pref,eref,L,Lice,e,T_max,T_min
     real(r8) :: this_lat,this_lon
     real(r8) :: a,b,c,d,ebr,fa,fb,fc,pbr,qbr,rbr,sbr,tol1,xm
     integer :: NTRY,LOOPMAX,i
     real(r8) :: tol,EPS

EPS=3.e-9
LOOPMAX = 100                   !* max number of iteration loops 
NTRY = 50
! Values for entropy
pref = 1000.0_r8           ! mb ref pressure.
eref = 6.106_r8           ! sat p at tfreez (mb)

! Invert the entropy equation -- use Brent's method
! Brent, R. P. Ch. 3-4 in Algorithms for Minimization Without Derivatives. Englewood Cliffs, NJ: Prentice-Hall, 1973.

Ts = Tfg                  ! Better first guess based on Tprofile from conv.

t2 = Tfg+30			!high bracket
t1 = Tfg-30			!low bracket
if (t1 < 40._r8) t1 = 40._r8      !stops suggesting temperatures which lead to entropy being undefined



   L = rl - (cpliq - cpwv)*(t2-tfreez) 
   Lice = rlice + (cpliq - cpice)*(t2-tfreez)
   T_min = tfreez - 10
   T_max = tfreez

   esat = c1*exp(c2*(t2-tfreez)/(c3+t2-tfreez)) ! Bolton (eq. 10)
   call esat_buck(t2, esat)

   e = qt*p/(qt+eps1*(1-qt))                    ! Using e instead of esat here as well, may be justified
   qsat = eps1*esat/(p-e*(1-eps1))     
   qv = min(qt,qsat)
   qc = max(0._r8,(qt-qsat)) 
   
!!$   e = qv*p / (eps1*(1-qv) +qv)  ! Bolton (eq. 16)
!!$   fs2 = ((1-qt)*cpres + qt*cpliq)*log( t2/tfreez) - (1-qt)*rgas*log( (p-e)/pref ) + &
!!$        L*qv/t2 - qv*rh2o*log(e/p) - s

   ql = qc * (T - T_min)/(T_max - T_min)
   qi = qc * (T_max - T)/(T_max - T_min)
   ql = max(0._r8,ql)
   ql = min(ql,qc)
   qi = max(0._r8,qi)
   qi = min(qi,qc)

   e = qv*p / (eps1*(1-qv) +qv)  ! Bolton (eq. 16)
   qsat = eps1*esat/(p-e*(1-eps1))     
   fs2 = ((1-qt)*cpres + qt*cpliq)*log( t2/tfreez) - (1-qt)*rgas*log( (p-e)/pref ) + &
   L*qv/t2 - Lice*qi/t2 - qv*rh2o*log(e/esat) - s


   L = rl - (cpliq - cpwv)*(t1-tfreez)
   Lice = rlice + (cpliq - cpice)*(t1-tfreez)
   T_min = tfreez - 10
   T_max = tfreez

   esat = c1*exp(c2*(t1-tfreez)/(c3+t1-tfreez)) ! Bolton (eq. 10)
   call esat_buck(t1, esat)

   e = qt*p/(qt+eps1*(1-qt))                    ! Using e instead of esat here as well, may be justified
   qsat = eps1*esat/(p-e*(1-eps1))     
   qv = min(qt,qsat)
   qc = max(0._r8,(qt-qsat))              

!!$   e = qv*p / (eps1*(1-qv) +qv)  ! Bolton (eq. 16)
!!$   fs1 = ((1-qt)*cpres + qt*cpliq)*log( t1/tfreez) - (1-qt)*rgas*log( (p-e)/pref ) + &
!!$        L*qv/t1 - qv*rh2o*log(e/p) - s

   ql = qc * (t1 - T_min)/(T_max - T_min)
   qi = qc * (T_max - t1)/(T_max - T_min)
   ql = max(0._r8,ql)
   ql = min(ql,qc)
   qi = max(0._r8,qi)
   qi = min(qi,qc)

   e = qv*p / (eps1*(1-qv) +qv)  ! Bolton (eq. 16)
   fs1 = ((1-qt)*cpres + qt*cpliq)*log( t1/tfreez) - (1-qt)*rgas*log( (p-e)/pref ) + &
   L*qv/t1 - Lice*qi/t1 - qv*rh2o*log(e/esat) - s

   !if (rcall==1 .and. icol==1) then
   !   write(iulog,*) "*** IENTROPY. t1: ", t1, "qt: ", qt, "qv: ", qv, "ql: ", ql, "e: ", e
   !end if
   
a=t1
b=t2
fa=fs1
!cas 3/12: bug  fb should be assigned to fs2; not fb=fs1
fb=fs2   
c=b
fc=fb
tol=0.0001_r8

! if (rcall==1 .and. icol==5 .and. lchnk==33) then
!    write(iulog,*) "*** IENTROPY. t2: ", t2, "fs2: ", fs2, "t1: ", t1, "fs1: ", fs1
! end if

converge: do i=0, LOOPMAX
if ((fb > 0.0 .and. fc > 0.0) .or. (fb < 0.0 .and. fc < 0.0)) then
			c=a
			fc=fa
			d=b-a
			ebr=d
		end if
if (abs(fc) < abs(fb)) then
			a=b
			b=c
			c=a
			fa=fb
			fb=fc
			fc=fa
		end if
tol1=2.0_r8*EPS*abs(b)+0.5_r8*tol
xm=0.5_r8*(c-b)
if (abs(xm) <= tol1 .or. fb == 0.0) then
   Ts=b
   exit converge
end if
if (abs(ebr) >= tol1 .and. abs(fa) > abs(fb)) then
   sbr=fb/fa
   if (a == c) then
      pbr=2.0_r8*xm*sbr
      qbr=1.0_r8-sbr
   else
      qbr=fa/fc
      rbr=fb/fc
      pbr=sbr*(2.0_r8*xm*qbr*(qbr-rbr)-(b-a)*(rbr-1.0_r8))
      qbr=(qbr-1.0_r8)*(rbr-1.0_r8)*(sbr-1.0_r8)
   end if
   if (pbr > 0.0) qbr=-qbr
   pbr=abs(pbr)
   if (2.0_r8*pbr  <  min(3.0_r8*xm*qbr-abs(tol1*qbr),abs(ebr*qbr))) then
      ebr=d
      d=pbr/qbr
   else
      d=xm
      ebr=d
   end if
else
   d=xm
   ebr=d
end if
a=b
fa=fb
b=b+merge(d,sign(tol1,xm), abs(d) > tol1 )
t2=b

   L = rl - (cpliq - cpwv)*(t2-tfreez)
   Lice = rlice + (cpliq - cpice)*(t2-tfreez)
   T_min = tfreez - 10
   T_max = tfreez

   esat = c1*exp(c2*(t2-tfreez)/(c3+t2-tfreez)) ! Bolton (eq. 10)
   call esat_buck(t2, esat)

   e = qt*p / (eps1*(1-qt) +qt)  ! using e instead of esat, may be justified
   qsat = eps1*esat/(p-e*(1-eps1))     
   qv = min(qt,qsat)
   qc = max(0._r8,(qt-qsat))

!!$   e = qv*p / (eps1*(1-qv) +qv)  ! Bolton (eq. 16)
!!$   fs2 = ((1-qt)*cpres + qt*cpliq)*log( t2/tfreez) - (1-qt)*rgas*log( (p-e)/pref ) + &
!!$        L*qv/t2 - qv*rh2o*log(e/p) - s

   
   ql = qc * (t2 - T_min)/(T_max - T_min)
   qi = qc * (T_max - t2)/(T_max - T_min)
   ql = max(0._r8,ql)
   ql = min(ql,qc)
   qi = max(0._r8,qi)
   qi = min(qi,qc)
   e = qv*p / (eps1*(1-qv) +qv)  ! Bolton (eq. 16)
   qsat = eps1*esat/(p-e*(1-eps1))     
   fs2 = ((1-qt)*cpres + qt*cpliq)*log( t2/tfreez) - (1-qt)*rgas*log( (p-e)/pref ) + &
      L*qv/t2 - Lice*qi/t2 - qv*rh2o*log(e/esat) - s

fb=fs2

! if (rcall==1 .and. icol==5 .and. lchnk==33) then
!    write(iulog,*) "*** IENTROPY. new T: ", t2, "new fs2: ", fs2
!  end if


   if (i .eq. LOOPMAX - 1) then
      this_lat = get_rlat_p(lchnk, icol)*57.296_r8
      this_lon = get_rlon_p(lchnk, icol)*57.296_r8
      write(iulog,*) '*** ZM_CONV: IENTROPY: Failed and about to exit, info follows ****'
      write(iulog,100) 'ZM_CONV: IENTROPY. Details: call#,lchnk,icol= ',rcall,lchnk,icol, &
       ' lat: ',this_lat,' lon: ',this_lon, &
       ' P(mb)= ', p, ' Tfg(K)= ', Tfg, ' qt(g/kg) = ', 1000._r8*qt, &
       ' qsat(g/kg) = ', 1000._r8*qsat,', s(J/kg) = ',s
      call endrun('**** ZM_CONV IENTROPY: Tmix did not converge ****')
   end if
enddo converge

! Replace call to satmixutils.

esat = c1*exp(c2*(Ts-tfreez)/(c3+Ts-tfreez))   
call esat_buck(Ts, esat)

qsat=eps1*esat/(p-esat*(1-eps1))

qv = min(qt,qsat)                             !       /* check for saturation */
qc = max(0._r8,(qt-qsat))

T = Ts 

! if (rcall==9) then
!    write(iulog,*) " "
!    write(iulog,*) "ientropy. input eta:", s, "p:", p, "qt:", qt, "tfg:", Tfg
!    write(iulog,*) "ientropy. output T:", T, "qsat:", qsat, "qv:", qv, "qc:", qc
! end if

 100    format (A,I1,I4,I4,7(A,F6.2))

return
end SUBROUTINE ientropy

!
!-----------------------------------------------------------------------------------------
SUBROUTINE itv (rcall,icol,lchnk,p,tp,qe,te,qsat,qmx,Tfg, k, lcl_est)
   !-----------------------------------------------------------------------------------------
   !
   ! p(mb), Tp, Tfg/T(K), qt/qv(kg/kg). 
   ! Inverts pressure, local T, and total water qt 
   ! for T_detrain and saturated vapor mixing ratio
   ! neglects effect of mass-loading on condensates
   ! 
   
        use phys_grid, only: get_rlon_p, get_rlat_p
   
        integer, intent(in) :: icol, lchnk, rcall, k, lcl_est
        real(r8), intent(in)  :: te, p, Tfg, qe, qmx
        real(r8), intent(out) :: qsat, tp
        real(r8) :: Ts,fs1,fs2,esat,t2,t1     
        real(r8) :: pref,eref
        real(r8) :: this_lat,this_lon
        real(r8) :: a,b,c,d,ebr,fa,fb,fc,pbr,qbr,rbr,sbr,tol1,xm
        integer :: NTRY,LOOPMAX,i
        real(r8) :: tol,EPS
   
   EPS=3.e-8
   LOOPMAX = 100                   !* max number of iteration loops 
   NTRY = 50
   ! Values for entropy
   pref = 1000.0_r8           ! mb ref pressure.
   eref = 6.106_r8           ! sat p at tfreez (mb)
   
   ! Invert the entropy equation -- use Brent's method
   ! Brent, R. P. Ch. 3-4 in Algorithms for Minimization Without Derivatives. Englewood Cliffs, NJ: Prentice-Hall, 1973.
   
   Ts = Tfg                  ! Better first guess based on Tprofile from conv.
   
   t2 = Tfg+30			!high bracket
   t1 = Tfg-30			!low bracket
   
   
      esat = c1*exp(c2*(t2-tfreez)/(c3+t2-tfreez)) ! Bolton (eq. 10)
      call esat_buck(t2, esat)

      qsat = eps1*esat/(p-esat*(1-eps1))
      if ((k < lcl_est + 4) .and. (qmx/qsat > 0.4_r8)) then
         qsat = min(qsat, qmx + 0.01_r8)  !allow a small degree of enhanced q in detraining in very limited conditions
      else 
         qsat = min(qsat, qmx)
      end if
      fs2 = (t2 * (1 - mu_red*qsat)) - (te * (1 - mu_red*qe)) !buoyancy difference
   
   
      esat = c1*exp(c2*(t1-tfreez)/(c3+t1-tfreez)) ! Bolton (eq. 10)
      call esat_buck(t1, esat)

      qsat = eps1*esat/(p-esat*(1-eps1))     
      if ((k < lcl_est + 4) .and. (qmx/qsat > 0.4_r8)) then
         qsat = min(qsat, qmx + 0.01_r8)  !allow a small degree of enhanced q in detraining in very limited conditions
      else 
         qsat = min(qsat, qmx)
      end if
      fs1 = (t1 * (1 - mu_red*qsat)) - (te * (1 - mu_red*qe)) !buoyancy difference
      
   a=t1
   b=t2
   fa=fs1
   !cas 3/12: bug  fb should be assigned to fs2; not fb=fs1
   fb=fs2   
   c=b
   fc=fb
   tol=0.001_r8
   
   converge: do i=0, LOOPMAX
   if ((fb > 0.0 .and. fc > 0.0) .or. (fb < 0.0 .and. fc < 0.0)) then
            c=a
            fc=fa
            d=b-a
            ebr=d
         end if
   if (abs(fc) < abs(fb)) then
            a=b
            b=c
            c=a
            fa=fb
            fb=fc
            fc=fa
         end if
   tol1=2.0_r8*EPS*abs(b)+0.5_r8*tol
   xm=0.5_r8*(c-b)
   if (abs(xm) <= tol1 .or. fb == 0.0) then
      Ts=b
      exit converge
   end if
   if (abs(ebr) >= tol1 .and. abs(fa) > abs(fb)) then
      sbr=fb/fa
      if (a == c) then
         pbr=2.0_r8*xm*sbr
         qbr=1.0_r8-sbr
      else
         qbr=fa/fc
         rbr=fb/fc
         pbr=sbr*(2.0_r8*xm*qbr*(qbr-rbr)-(b-a)*(rbr-1.0_r8))
         qbr=(qbr-1.0_r8)*(rbr-1.0_r8)*(sbr-1.0_r8)
      end if
      if (pbr > 0.0) qbr=-qbr
      pbr=abs(pbr)
      if (2.0_r8*pbr  <  min(3.0_r8*xm*qbr-abs(tol1*qbr),abs(ebr*qbr))) then
         ebr=d
         d=pbr/qbr
      else
         d=xm
         ebr=d
      end if
   else
      d=xm
      ebr=d
   end if
   a=b
   fa=fb
   b=b+merge(d,sign(tol1,xm), abs(d) > tol1 )
   t2=b
   
      esat = c1*exp(c2*(t2-tfreez)/(c3+t2-tfreez)) ! Bolton (eq. 10)
      call esat_buck(t2, esat)

      qsat = eps1*esat/(p-esat*(1-eps1))     
      if ((k < lcl_est + 4) .and. (qmx/qsat > 0.4_r8)) then
         qsat = min(qsat, qmx + 0.01_r8)  !allow a small degree of enhanced q in detraining in very limited conditions
      else 
         qsat = min(qsat, qmx)
      end if      
      fs2 = (t2 * (1 - mu_red*qsat)) - (te * (1 - mu_red*qe)) !buoyancy difference
   
   fb=fs2
   
      if (i .eq. LOOPMAX - 1) then
         this_lat = get_rlat_p(lchnk, icol)*57.296_r8
         this_lon = get_rlon_p(lchnk, icol)*57.296_r8
         write(iulog,*) '*** ZM_CONV: ITV: Failed and about to exit, info follows ****'
         write(iulog,100) 'ZM_CONV: ITV. Details: call#,lchnk,icol= ',rcall,lchnk,icol, &
          ' lat: ',this_lat,' lon: ',this_lon, &
          ' P(mb)= ', p, ' Tfg(K)= ', Tfg, ' qe(g/kg) = ', 1000._r8*qe, &
          ' Te(K) = ', te
         call endrun('**** ZM_CONV ITV: Tmix did not converge ****')
      end if
   enddo converge
   
   ! Replace call to satmixutils.
   
   esat = c1*exp(c2*(Ts-tfreez)/(c3+Ts-tfreez))
   call esat_buck(Ts, esat)

   qsat=eps1*esat/(p-esat*(1-eps1))
   if ((k < lcl_est + 4) .and. (qmx/qsat > 0.4_r8)) then
      qsat = min(qsat, qmx + 0.01_r8)  !allow a small degree of enhanced q in detraining in very limited conditions
   else 
      qsat = min(qsat, qmx)
   end if   
   tp = Ts 
   
    100    format (A,I1,I4,I4,6(A,F6.2))
   
   return
   end SUBROUTINE itv

   !
!-----------------------------------------------------------------------------------------
   SUBROUTINE ientropy_downdraft (rcall,icol,lchnk,s,p,qt,T,qsat,Tfg, dry_dwndraft)
      !-----------------------------------------------------------------------------------------
      !
      ! p(mb), Tfg/T(K), qt/qv(kg/kg), s(J/kg). 
      ! Inverts entropy, pressure and total water qt 
      ! for T and saturated vapor mixing ratio
      ! Difference is that we force the parcel to be saturated, as the downdraft is
      ! So the useful quantities returned are the temperature and the added water amount
      ! 
      
           use phys_grid, only: get_rlon_p, get_rlat_p
      
           integer, intent(in) :: icol, lchnk, rcall
           real(r8), intent(in)  :: s, p, Tfg, qt
           real(r8), intent(out) :: qsat, T
           logical, intent(inout) :: dry_dwndraft
           real(r8) :: Ts,dTs,fs1,fs2,esat,t2,t1, dq, dqi, dql, T_max, T_min
           real(r8) :: pref,eref,L,Lice,e
           real(r8) :: this_lat,this_lon
           real(r8) :: a,b,c,d,ebr,fa,fb,fc,pbr,qbr,rbr,sbr,tol1,xm
           integer :: NTRY,LOOPMAX,i
           real(r8) :: tol,EPS
      
      !if (rcall==1 .and. icol==1) then
      !write(iulog,*) " "
      !write(iulog,*) " "
      !write(iulog,*) "*** IENTROPY_DOWNDRAFT. lchnk:", lchnk, "i:", icol
      !write(iulog,*) "*** IENTROPY_DOWNDRAFT. s:", s, "qt:", qt, "p:", p, "Tfg:", Tfg
      !end if

      EPS=3.e-8
      LOOPMAX = 100                   !* max number of iteration loops 
      NTRY = 50
      ! Values for entropy
      pref = 1000.0_r8           ! mb ref pressure.
      eref = 6.106_r8           ! sat p at tfreez (mb)
      T_max = tfreez
      T_min = tfreez-10
      
      ! Invert the entropy equation -- use Brent's method
      ! Brent, R. P. Ch. 3-4 in Algorithms for Minimization Without Derivatives. Englewood Cliffs, NJ: Prentice-Hall, 1973.
      
      Ts = Tfg                  ! Better first guess based on Tprofile from conv.
      
      t2 = Tfg+30			!high bracket
      t1 = Tfg-30			!low bracket
      
      
         L = rl - (cpliq - cpwv)*(t2-tfreez)
         Lice = rlice + (cpliq - cpice)*(t2-tfreez)

         esat = c1*exp(c2*(t2-tfreez)/(c3+t2-tfreez)) ! Bolton (eq. 10)
         call esat_buck(t2, esat)

         e = qt*p / (eps1*(1-qt) +qt)  ! Bolton (eq. 16)
         qsat = eps1*esat/(p-e*(1-eps1))
         qsat = min(0.99_r8,qsat)
         dq = (qsat - qt)/(1-qsat) !additional moisture needed to make plume saturated, should be positive
         
!!$         e = qsat*p / (eps1*(1-qsat) +qsat)  ! Bolton (eq. 16)
!!$         !pd = p-e
!!$         fs2 = ((1-qsat)*cpres + qsat*cpliq)*log( t2/tfreez) - (1-qsat)*rgas*log( (p-e)/pref ) + &
!!$              L*qsat/t2 - qsat*rh2o*log(e/p) - (s + dq*cpliq*log(t2/tfreez))/(1+dq)

         dql = dq * (t2 - T_min)/(T_max - T_min)
         dqi = dq * (T_max - t2)/(T_max - T_min)
         dql = max(0._r8,dql)
         dql = min(dql,dq)
         dqi = max(0._r8,dqi)
         dqi = min(dqi,dq)
         
         e = qsat*p / (eps1*(1-qsat) +qsat)  ! Bolton (eq. 16)
         !pd = p-e
         fs2 = ((1-qsat)*cpres + qsat*cpice)*log( t2/tfreez) - (1-qsat)*rgas*log( (p-e)/pref ) + &
         L*qsat/t2 - qsat*rh2o*log(e/esat) - (s + dq*cpliq*log(t2/tfreez) - dqi*Lice/t2)/(1+dq)
      
         if (rcall == 5._r8) write(iulog,*) "initial guesses. t2:", t2, "qsat2:", qsat, "dq:", dq, "fs2:", fs2
      
         L = rl - (cpliq - cpwv)*(t1-tfreez)
         Lice = rlice + (cpliq - cpice)*(t1-tfreez)

         esat = c1*exp(c2*(t1-tfreez)/(c3+t1-tfreez)) ! Bolton (eq. 10)
         call esat_buck(t1, esat)

         e = qt*p / (eps1*(1-qt) +qt)  ! Bolton (eq. 16)
         qsat = eps1*esat/(p-esat*(1-eps1))
         qsat = min(0.99_r8,qsat)
         dq = (qsat - qt)/(1-qsat) !additional moisture needed to make plume saturated, should be positive
         
!!$         e = qsat*p / (eps1*(1-qsat) +qsat)  ! Bolton (eq. 16)
!!$         fs1 = ((1-qsat)*cpres + qsat*cpliq)*log( t1/tfreez) - (1-qsat)*rgas*log( (p-e)/pref ) + &
!!$              L*qsat/t1 - qsat*rh2o*log(e/p) - (s + dq*cpliq*log( t1/tfreez))/(1+dq)

         dql = dq * (t2 - T_min)/(T_max - T_min)
         dqi = dq * (T_max - t2)/(T_max - T_min)
         dql = max(0._r8,dql)
         dql = min(dql,dq)
         dqi = max(0._r8,dqi)
         dqi = min(dqi,dq)
         
         e = qsat*p / (eps1*(1-qsat) +qsat)  ! Bolton (eq. 16)
         fs1 = ((1-qsat)*cpres + qsat*cpice)*log( t1/tfreez) - (1-qsat)*rgas*log( (p-e)/pref ) + &
         L*qsat/t1 - qsat*rh2o*log(e/esat) - (s + dq*cpliq*log( t1/tfreez) - dqi*Lice/t1)/(1+dq)

         if (rcall == 5._r8) write(iulog,*) "initial guesses. t1:", t1, "qsat2:", qsat, "dq:", dq, "fs2:", fs1
         
      a=t1
      b=t2
      fa=fs1
      !cas 3/12: bug  fb should be assigned to fs2; not fb=fs1
      fb=fs2   
      c=b
      fc=fb
      tol=0.001_r8
      
      converge: do i=0, LOOPMAX
      if ((fb > 0.0 .and. fc > 0.0) .or. (fb < 0.0 .and. fc < 0.0)) then
               c=a
               fc=fa
               d=b-a
               ebr=d
               !write(iulog,*) "fb:", fb, "fc:", fc, "b:", b, "c:", c, "d:", d, "ebr:", ebr
            end if
      if (abs(fc) < abs(fb)) then
               a=b
               b=c
               c=a
               fa=fb
               fb=fc
               fc=fa
               !write(iulog,*) "fb:", fb, "fc:", fc, "b:", b, "c:", c, "d:", d
            end if
      tol1=2.0_r8*EPS*abs(b)+0.5_r8*tol
      xm=0.5_r8*(c-b)
      !write(iulog,*) "xm:", xm, "tol1:", tol1
      if (abs(xm) <= tol1 .or. fb == 0.0) then
         !write(iulog,*) "Solution found. T:", b, "fs:", fb
         Ts=b
         exit converge
      end if
      if (abs(ebr) >= tol1 .and. abs(fa) > abs(fb)) then
         !write(iulog,*) "fb:", fb, "fa:", fa, "ebr:", ebr, "tol1:", tol1
         sbr=fb/fa
         if (a == c) then
            pbr=2.0_r8*xm*sbr
            qbr=1.0_r8-sbr
         else
            qbr=fa/fc
            rbr=fb/fc
            pbr=sbr*(2.0_r8*xm*qbr*(qbr-rbr)-(b-a)*(rbr-1.0_r8))
            qbr=(qbr-1.0_r8)*(rbr-1.0_r8)*(sbr-1.0_r8)
         end if
         !write(iulog,*) "sbr:", sbr, "pbr:", pbr, "qbr:", qbr
         if (pbr > 0.0) qbr=-qbr
         pbr=abs(pbr)
         if (2.0_r8*pbr  <  min(3.0_r8*xm*qbr-abs(tol1*qbr),abs(ebr*qbr))) then
            ebr=d
            d=pbr/qbr
         else
            d=xm
            ebr=d
         end if
         !write(iulog,*) "d:", d, "ebr:", ebr
      else
         d=xm
         ebr=d
      end if
      a=b
      fa=fb
      b=b+merge(d,sign(tol1,xm), abs(d) > tol1 )
      t2=b
      
      L = rl - (cpliq - cpwv)*(t2-tfreez)
      Lice = rlice + (cpliq - cpice)*(t2-tfreez)

      esat = c1*exp(c2*(t2-tfreez)/(c3+t2-tfreez)) ! Bolton (eq. 10)
      call esat_buck(t2, esat)

      e = qt*p / (eps1*(1-qt) +qt)  ! Bolton (eq. 16)
      qsat = eps1*esat/(p-esat*(1-eps1))
      qsat = min(0.99_r8, qsat)
      dq = (qsat - qt)/(1-qsat) !additional moisture needed to make plume saturated, should be positive
      
!!$      e = qsat*p / (eps1*(1-qsat) +qsat)  ! Bolton (eq. 16)
!!$      fs2 = ((1-qsat)*cpres + qsat*cpliq)*log( t2/tfreez) - (1-qsat)*rgas*log( (p-e)/pref ) + &
!!$           L*qsat/t2 - qsat*rh2o*log(e/p) - (s + dq*cpliq*log(t2/tfreez))/(1+dq)

      dql = dq * (t2 - T_min)/(T_max - T_min)
      dqi = dq * (T_max - t2)/(T_max - T_min)
      dql = max(0._r8,dql)
      dql = min(dql,dq)
      dqi = max(0._r8,dqi)
      dqi = min(dqi,dq)
      
      e = qsat*p / (eps1*(1-qsat) +qsat)  ! Bolton (eq. 16)
      fs2 = ((1-qsat)*cpres + qsat*cpice)*log( t2/tfreez) - (1-qsat)*rgas*log( (p-e)/pref ) + &
      L*qsat/t2 - qsat*rh2o*log(e/esat) - (s + dq*cpliq*log(t2/tfreez) - dqi*Lice/t2)/(1+dq)
!!$      
      fb=fs2

      if (rcall == 5._r8) write(iulog,*) "new guess. t2:", t2, "qsat2:", qsat, "dq:", dq, "fs2:", fs2
      
         if (i .eq. LOOPMAX - 1) then
            this_lat = get_rlat_p(lchnk, icol)*57.296_r8
            this_lon = get_rlon_p(lchnk, icol)*57.296_r8
            write(iulog,*) '*** ZM_CONV: IENTROPY_DOWNDRAFT: Failed and about to exit, info follows ****'
            write(iulog,100) 'ZM_CONV: IENTROPY_DOWNDRAFT. Details: call#,lchnk,icol= ',rcall,lchnk,icol, &
             ' lat: ',this_lat,' lon: ',this_lon, &
             ' P(mb)= ', p, ' Tfg(K)= ', Tfg, ' qt(g/kg) = ', 1000._r8*qt, &
             ' qsat(g/kg) = ', 1000._r8*qsat,', s(J/kg) = ',s
            call endrun('**** ZM_CONV IENTROPY_DOWNDRAFT: Tmix did not converge ****')
         end if
      enddo converge
      
      ! Replace call to satmixutils.
      
      esat = c1*exp(c2*(Ts-tfreez)/(c3+Ts-tfreez))
      call esat_buck(Ts, esat)

      qsat=eps1*esat/(p-esat*(1-eps1))
      if (qsat > 0.99_r8) then
         qsat = 0.99_r8
         dry_dwndraft = .true. !downdraft has to become subsaturated from here
      end if
      
      T = Ts 
      
       100    format (A,I1,I4,I4,7(A,F6.2))
      
      return
   end SUBROUTINE ientropy_downdraft

end module zm_conv

