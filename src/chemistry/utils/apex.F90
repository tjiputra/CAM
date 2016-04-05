module apex
!
! April, 2013: B. Foster (NCAR/HAO)
!
! This is a refactored version of the legacy apex code, originally written
! by Art Richmond and Roy Barnes, and others in the 1995-2000 timeframe.
! This new version is written in free-format fortran90. Subroutines and
! module data may be use-associated from this module. 
!
! Original reference for the legacy code:
!          Richmond, A. D., Ionospheric Electrodynamics Using Magnetic Apex
!          Coordinates, J. Geomag. Geoelectr., 47, 191-212, 1995. 
!
! This code should produce near-identical results as the legacy code, altho 
! the refactored version does not provide all subroutines and options available 
! in the old code, notably the ability to write and read-back an external file.
! 
! A typical calling sequence for a code calling this module is as follows:
!
! subroutine ggrid (legacy SUBROUTINE GGRID): 
!   Make a global lat,lon,alt grid for use in later calls (optional)
!
! subroutine apex_mka (legacy SUBROUTINE APXMKA): 
!   Make magnetic arrays x,y,z,v for use in later routines
!   (geographic lat,lon grid and altitudes are input) 
!   (This must be called before apex_mall and apex_q2g)
! 
! subroutine apex_mall (legacy ENTRY APXMALL): 
!   Calculate modified Apex coordinates and other magnetic field parameters
!   (usually called from lat,lon,alt nested loop)
! 
! subroutine apex_q2g (legacy ENTRY APXQ2G): 
!   Convert from quasi-dipole to geodetic coordinates
!   (usually called from lat,lon,alt nested loop)
!
  use shr_kind_mod,  only : r8 => shr_kind_r8
  use cam_logfile,   only : iulog
  use cam_abortutils,only : endrun

  implicit none

  private
  public :: apex_mka
  public :: apex_mall
  public :: apex_dypol
  public :: apex_subsol
  public :: apex_magloctm
 
  real(r8),parameter :: re = 6371.2_r8, eps = 1.e-5_r8

  real(r8),allocatable,save :: &
    xarray(:,:,:), & ! cos(quasi-dipole latitude)*cos(apex longitude)
    yarray(:,:,:), & ! cos(quasi-dipole latitude)*sin(apex longitude)
    zarray(:,:,:), & ! sin(quasi-dipole latitude)
    varray(:,:,:)    ! (VMP/VP)*((RE+ALT)/RE)**2
!
! This grid (geolat,geolon,geoalt is equivalent to gdlat,gdlon,gdalt, 
! as passed to apex_mka.
!
  integer :: nglat,nglon,ngalt
  real(r8),allocatable,save :: geolat(:), geolon(:), geoalt(:)    

  integer,parameter :: nmax=13
  integer,parameter :: ncoef = nmax*nmax + 2*nmax + 1 ! 196
  real(r8),dimension(ncoef) :: &
    gb, & ! Coefficients for magnetic field calculation
    gv    ! Coefficients for magnetic potential calculation
!
  real(r8) :: &
    rtd,  & ! radians to degrees
    dtr,  & ! degrees to radians
    pola    ! pole angle (deg); when geographic lat is poleward of pola,
            ! x,y,z,v arrays are forced to be constant (pola=89.995) 

  real(r8),parameter ::       & ! Formerly common /APXCON/
    req  = 6378.160_r8,       & ! Equatorial earth radius
    precise = 7.6e-11_r8,     & ! Precision factor
    glatlim = 89.9_r8,        & ! Limit above which gradients are recalculated
    xmiss = -32767._r8
!
! colat,elon,vp,ctp,stp were in two commons in legacy code:
! /APXDIPL/ and /DIPOLE/. Need to check if these need to be separated.
!
  real(r8) ::  & ! Formerly /APXDIPL/ and /DIPOLE/
    colat, & ! Geocentric colatitude of geomagnetic dipole north pole (deg)
    elon,  & ! East longitude of geomagnetic dipole north pole (deg)
    vp,    & ! Magnitude, in T.m, of dipole component of magnetic
             ! potential at geomagnetic pole and geocentric radius re
    ctp,stp
!
  real(r8) ::  & ! Formerly /FLDCOMD/
    bx,    & ! X comp. of field vector at the current tracing point (Gauss)
    by,    & ! Y comp. of field vector at the current tracing point (Gauss)
    bz,    & ! Z comp. of field vector at the current tracing point (Gauss)
    bb       ! Magnitude of field vector at the current tracing point (Gauss)

  real(r8) ::      & ! Formerly /APXIN/
    yapx(3,3)    ! Matrix of cartesian coordinates (loaded columnwise) 
!
! /ITRA/ was only in subs linapx and itrace, so can probably be removed from module data
!
  integer ::   & ! Formerly /ITRA/ 
    nstp         ! Step count. Incremented in sub linapx.
  real(r8)    ::   & 
    y(3),      & ! Array containing current tracing point cartesian coordinates.
    yp(3),     & ! Array containing previous tracing point cartesian coordinates.
    sgn,       & ! Determines direction of trace. Set in subprogram linapx
    ds           ! Step size (Km) Computed in subprogram linapx.

  real(r8) ::         & ! limits beyond which east-west gradients are computed 
    glatmn,glatmx   ! differently to avoid potential underflow (apex_mka)

  logical :: initialized = .false.
contains
!-----------------------------------------------------------------------
subroutine ggrid(nvert,glatmin,glatmax,glonmin,glonmax,altmin,altmax, &
                 gplat,gplon,gpalt,mxlat,mxlon,mxalt,nlat,nlon,nalt)
!
! Given desired range of geographic latitude, longitude and altitude, 
! choose an appropriate grid that can be used in subsequent calls to 
! subs apex_mka, apex_mall, apex_q2g.
!
! Input args:
  integer,intent(in) :: nvert,mxlat,mxlon,mxalt
  real(r8),intent(in) :: glatmin,glatmax,glonmin,glonmax,altmin,altmax
!
! Output args:
  integer,intent(out) :: nlat,nlon,nalt
  real(r8),intent(out) :: gplat(mxlat),gplon(mxlon),gpalt(mxalt)
!
! Local:
  real(r8) :: dlon,dlat,diht,dnv,glonmaxx,x
  integer :: nlatmin,nlatmax,nlonmin,nlonmax,naltmin,naltmax
  integer :: i,j,k,kk
  character(len=128) :: errmsg
!
! Check inputs:
  if (glatmin > glatmax) then
    write(errmsg,"('>>> ggrid: glatmin=',f9.2,' must be <= glatmax=',f9.2)") glatmin,glatmax
    write(iulog,*) errmsg
    call endrun( trim(errmsg) )
  endif
  if (glonmin > glonmax) then
    write(errmsg,"('>>> ggrid: glonmin=',f9.2,' must be <= glonmax=',f9.2)") glonmin,glonmax
    write(iulog,*) errmsg
    call endrun( trim(errmsg) )
  endif
  if (altmin > altmax) then
    write(errmsg,"('>>> ggrid: altmin=',f9.2,' must be <= altmax=',f9.2)") altmin,altmax
    write(iulog,*) errmsg
    call endrun( trim(errmsg) )
  endif
!
! Init outputs:
  nlat = 0 ; nlon = 0 ; nalt = 0
  gplat = 0._r8 ; gplon = 0._r8 ; gpalt = 0._r8
!
  dnv = dble(nvert)
  dlon = 360._r8 / (5._r8*dnv)
  dlat = 180._r8 / (3._r8*dnv)
  diht = 1._r8   / dnv

  nlatmin = max(int((glatmin+90._r8)/dlat),0)
  nlatmax = min(int((glatmax+90._r8)/dlat+1._r8),3*nvert)
  nlonmin = max(int((glonmin+180._r8)/dlon),0)
 
  glonmaxx = min(glonmax,glonmin+360._r8)
  nlonmax = min(int((glonmaxx+180._r8)/dlon+1._r8),10*nvert)
    
  x = re/(re+altmax)/diht-eps
  naltmin = max(x,1._r8)
  naltmin = min(naltmin,nvert-1)
  x = re/(re+altmin)/diht+eps
  i = x + 1._r8
  naltmax = min(i,nvert)

  nlat = nlatmax - nlatmin + 1
  nlon = nlonmax - nlonmin + 1
  nlon = min(nlon,5*nvert+1)
  nalt = naltmax - naltmin + 1

  do j=1,nlat
    gplat(j) = dlat*dble(nlatmin+j-1) - 90._r8
  enddo
  do i=1,nlon
    gplon(i) = dlon*dble(nlonmin+i-1) - 180._r8
  enddo
  do k=1,nalt
    kk = naltmax - k +1
    gpalt(k) = re*(dble(nvert-kk) - eps) / (dble(kk)+eps)
  enddo
  if (gplon(nlon-1) >= glonmax) nlon = nlon-1
  gpalt(1) = max(gpalt(1),0._r8)

  write(iulog,"('ggrid: nlat=',i4,' gplat=',/,(6f9.2))") nlat,gplat
  write(iulog,"('ggrid: nlon=',i4,' gplon=',/,(6f9.2))") nlon,gplon
  write(iulog,"('ggrid: nalt=',i4,' gpalt=',/,(6f9.2))") nalt,gpalt

end subroutine ggrid
!-----------------------------------------------------------------------
subroutine apex_mka(date,gplat,gplon,gpalt,nlat,nlon,nalt,ier)
!
! Given a 3d lat,lon,altitude grid, calculate x,y,z,v arrays in module
! data above. These arrays are used later for calculating quantities
! involving gradients of Apex coordinates, such as base vectors in the
! Modified-Apex and Quasi-Dipole systems.
!
! This defines module 3d data xarray,yarray,zarray,varray
!
! Input args:
  real(r8),intent(in) :: date              ! year and fraction
  integer, intent(in) :: nlat,nlon,nalt ! dimensions of 3d grid
  real(r8),intent(inout) :: gplat(nlat),gplon(nlon),gpalt(nalt)
!
! Output args:
  integer,intent(out) :: ier
!
! Local:
  integer :: i,j,k,kpol,istat
  real(r8) :: reqore,rqorm1,cp,ct,st,sp,stmcpm,stmspm,ctm
  real(r8) :: aht,alat,phia,bmag,xmag,ymag,zdown,vmp ! apex_sub output
  real(r8) :: vnor,rp,reqam1,a,slp,clp,phiar

  if (initialized) return
!
! Some parts of the legacy apex code use constants to set dtr,rtd,
! other parts use rtd=45./atan(1.), dtr=1./rtd. Differences are
! on the order of 1.e-18 to 1.e-14. Here, the atan method is used. 
!
!  rtd  = 5.72957795130823E1
!  dtr  = 1.745329251994330E-2
!
   rtd  = 45._r8/atan(1._r8)
   dtr  = 1._r8/rtd
!
! pola:
!   Pole angle (deg); when the geographic latitude is poleward of POLA, 
!   X,Y,Z,V are forced to be constant for all longitudes at each altitude.  
!   This makes POLA = 89.995
!
  pola = 90._r8-sqrt(precise)*rtd    ! Pole angle (deg)

  ier = 0
!
! Allocate 3d x,y,z,v arrays:
! These are not deallocated by this module. They can be deallocated
! by the calling program following the last call to the apex subs.
!
  if (.not.allocated(xarray)) then
    allocate(xarray(nlat,nlon,nalt),stat=istat)
    if (istat /= 0) call endrun( 'allocate xarray' )
    xarray = 0._r8
  endif
  if (.not.allocated(yarray)) then
    allocate(yarray(nlat,nlon,nalt),stat=istat)
    if (istat /= 0) call endrun( 'allocate yarray' )
    yarray = 0._r8
  endif
  if (.not.allocated(zarray)) then
    allocate(zarray(nlat,nlon,nalt),stat=istat)
    if (istat /= 0) call endrun( 'allocate zarray' )
    zarray = 0._r8
  endif
  if (.not.allocated(varray)) then
    allocate(varray(nlat,nlon,nalt),stat=istat)
    if (istat /= 0) call endrun( 'allocate varray' )
    varray = 0._r8
  endif
!
! Set geographic grid in module data for later reference:
! (these also are not deallocated by this module)
!
  nglon=nlon ; nglat=nlat ; ngalt=nalt
  allocate(geolat(nglat),stat=istat)
  allocate(geolon(nglon),stat=istat)
  allocate(geoalt(ngalt),stat=istat)
  geolat(:) = gplat(:)
  geolon(:) = gplon(:)
  geoalt(:) = gpalt(:)
!
! Set coefficients gb,gv (module data) for requested year:
!
  call cofrm(date)

! write(iulog,"('apex_mka after cofrm: ncoef=',i4,' gb=',/,(6f12.3))") ncoef,gb
! write(iulog,"('apex_mka after cofrm: ncoef=',i4,' gv=',/,(6f12.3))") ncoef,gv

  call apex_dypol(colat,elon,vp)

  ctp = cos(colat*dtr)
  stp = sin(colat*dtr)

  reqore = req/re
  rqorm1 = reqore-1._r8

  do j=1,nlat
    ct = sin(gplat(j)*dtr)
    st = cos(gplat(j)*dtr)
    kpol = 0
    if (abs(gplat(j)) > pola) kpol = 1
    do i=1,nlon
      if (kpol==1.and.i > 1) then
        xarray(j,i,:) = xarray(j,1,:)
        yarray(j,i,:) = yarray(j,1,:)
        zarray(j,i,:) = zarray(j,1,:)
        varray(j,i,:) = varray(j,1,:)
        cycle
      endif  
      cp = cos((gplon(i)-elon)*dtr)
      sp = sin((gplon(i)-elon)*dtr)
!
!  ctm   is pseudodipole component of z
! -ctm   is pseudodipole component of v
!  stmcpm is pseudodipole component of x
!  stmspm is pseudodipole component of y
!
      ctm = ctp*ct + stp*st*cp
      stmcpm = st*ctp*cp - ct*stp
      stmspm = st*sp
      do k=1,nalt
        call apex_sub(date,gplat(j),gplon(i),gpalt(k),&
          aht,alat,phia,bmag,xmag,ymag,zdown,vmp)

        vnor = vmp/vp
        rp = 1._r8 + gpalt(k)/re
        varray(j,i,k) = vnor*rp*rp + ctm
        reqam1 = req*(aht-1._r8)
        slp = sqrt(max(reqam1-gpalt(k),0._r8)/(reqam1+re))
!
! Reverse sign of slp in southern magnetic hemisphere
!
        if (zdown.lt.0._r8) slp = -slp
        clp = sqrt (rp/(reqore*aht-rqorm1))
        phiar = phia*dtr
        xarray(j,i,k) = clp*cos (phiar) - stmcpm
        yarray(j,i,k) = clp*sin (phiar) - stmspm
        zarray(j,i,k) = slp - ctm
      enddo ! k=1,nalt
    enddo ! i=1,nlon
  enddo ! j=1,nlat
!
! Establish for this grid polar latitude limits beyond which east-west
! gradients are computed differently to avoid potential underflow
! (glatmx,glatmn are in module data, glatlim is parameter constant)
!
  glatmx = max( glatlim,gplat(nlat-2))
  glatmn = min(-glatlim,gplat(2))

  initialized = .true.

end subroutine apex_mka
!-----------------------------------------------------------------------
subroutine apex_mall(glat,glon,alt,hr, b,bhat,bmag,si,alon,xlatm,vmp,w,&
  d,be3,sim,d1,d2,d3,e1,e2,e3,xlatqd,f,f1,f2,ier)
!
! Compute Modified Apex coordinates, quasi-dipole coordinates,
! base vectors and other parameters by interpolation from
! precalculated arrays. Subroutine apex_mka must be called
! before calling this subroutine.
!
! Args:
  real(r8),intent(in)  :: & ! Input
    glat             ,& ! Geographic (geodetic) latitude (deg)
    glon             ,& ! Geographic (geodetic) longitude (deg)
    alt              ,& ! Altitude (km)
    hr                  ! Reference altitude (km)

  real(r8),intent(out) :: & ! Output
    b(3)             ,& ! Magnetic field components (east, north, up), in nT    
    bhat(3)          ,& ! components (east, north, up) of unit vector along 
                        ! geomagnetic field direction
    bmag             ,& ! Magnitude of magnetic field (nT)
    si               ,& ! sin(i)
    alon             ,& ! Apex longitude = modified apex longitude = 
                        ! quasi-dipole longitude (deg)
    xlatm            ,& ! Modified Apex latitude (deg)
    vmp              ,& ! Magnetic potential (T.m)
    w                ,& ! W of Richmond reference above, in km**2 /nT (i.e., 10**15 m**2 /T)
    d                ,& ! D of Richmond reference above
    be3              ,& ! B_e3 of reference above (= Bmag/D), in nT
    sim              ,& ! sin(I_m) described in Richmond reference above
    xlatqd           ,& ! Quasi-dipole latitude (deg)
    f                ,& ! F described in ref above for quasi-dipole coordinates
    f1(2),f2(2)         ! Components (east, north) of base vectors
!
  real(r8),dimension(3),intent(out) :: d1,d2,d3,e1,e2,e3 ! Components of base vectors
  integer,intent(out) :: ier ! error return
!
! Local:
  real(r8) :: glonloc,cth,sth,glatx,clm,r3_2
  real(r8) :: fx,fy,fz,fv
  real(r8) :: dfxdth,dfydth,dfzdth,dfvdth, &
              dfxdln,dfydln,dfzdln,dfvdln, &
              dfxdh ,dfydh ,dfzdh ,dfvdh
  real(r8),dimension(3) :: gradx,grady,gradz,gradv, grclm,clmgrp,rgrlp
  real(r8) ::                    & ! dummies for polar calls to intrp
    fxdum,fydum,fzdum,fvdum,     &
    dmxdth,dmydth,dmzdth,dmvdth, &
    dmxdh,dmydh,dmzdh,dmvdh
!
! Init:
!
  ier = 0
  glonloc = glon

  call intrp (glat,glonloc,alt, geolat,geolon,geoalt,nglat,nglon,ngalt, &
             fx,fy,fz,fv,                                               &
             dfxdth,dfydth,dfzdth,dfvdth,                               &
             dfxdln,dfydln,dfzdln,dfvdln,                               &
             dfxdh ,dfydh ,dfzdh ,dfvdh, ier)

  if (ier /= 0) then
    call setmiss(xmiss,xlatm,alon,vmp,b,bmag,be3,sim,si,f,d,w, &
      bhat,d1,d2,d3,e1,e2,e3,f1,f2)
    write(iulog,"('apex_mall called setmiss: glat,glon,alt=',3f12.3)") &
      glat,glon,alt
    return
  endif

  call adpl(glat,glonloc,cth,sth,fx,fy,fz,fv, &
    dfxdth,dfydth,dfzdth,dfvdth,dfxdln,dfydln,dfzdln,dfvdln)

  call gradxyzv(alt,cth,sth, &
    dfxdth,dfydth,dfzdth,dfvdth,dfxdln,dfydln,dfzdln,dfvdln, &
    dfxdh,dfydh,dfzdh,dfvdh,gradx,grady,gradz,gradv)
!
! If the point is very close to either the North or South
! geographic pole, recompute the east-west gradients after
! stepping a small distance from the pole.
!
  if (glat > glatmx .or. glat < glatmn) then
    glatx = glatmx
    if (glat < 0._r8) glatx = glatmn

    call intrp (glatx,glonloc,alt, geolat,geolon,geoalt,nglat,nglon,ngalt, &
               fxdum,fydum,fzdum,fvdum,                                    &
               dmxdth,dmydth,dmzdth,dmvdth,dfxdln,dfydln,dfzdln,           &
               dfvdln,dmxdh,dmydh,dmzdh,dmvdh, ier)

    call adpl(glatx,glonloc,cth,sth,fxdum,fydum,fzdum,fvdum, &
      dmxdth,dmydth,dmzdth,dmvdth,dfxdln,dfydln,dfzdln,dfvdln)

    call grapxyzv(alt,cth,sth,dfxdln,dfydln,dfzdln,dfvdln, &
      gradx,grady,gradz,gradv)
  endif

  call gradlpv(hr,alt,fx,fy,fz,fv,gradx,grady,gradz,gradv, &
    xlatm,alon,vmp,grclm,clmgrp,xlatqd,rgrlp,b,clm,r3_2)

  call basevec(hr,xlatm,grclm,clmgrp,rgrlp,b,clm,r3_2, &
               bmag,sim,si,f,d,w,bhat,d1,d2,d3,e1,e2,e3,f1,f2)

  be3 = bmag/d
  ier = 0

end subroutine apex_mall
!-----------------------------------------------------------------------
subroutine apex_q2g(qdlat,qdlon,alt,gdlat,gdlon,ier)
!
! Convert from quasi-dipole to geodetic coordinates. This subroutine
! (input magnetic, output geodetic) is the functional inverse of 
! subroutine apex_mall (input geodetic, output magnetic). Sub apex_mka
! must be called before this routine.
!
! Args:
  real(r8),intent(in) ::  & ! inputs
    qdlat,                & ! quasi-dipole latitude (deg)
    qdlon,                & ! quasi-dipole longitude (deg)
    alt                     ! altitude (km)

  real(r8),intent(out) :: & ! outputs
    gdlat,            & ! geodetic latitude (deg)
    gdlon               ! geodetic longitude (deg)
  integer,intent(out) :: ier ! error return
!
! Local:
  real(r8) :: x0,y0,z0,xnorm,xdif,ydif,zdif,dist2,hgrd2e,hgrd2n,hgrd2,&
    angdist,distlon,glatx,cal,sal,coslm,slm,cad,sad,slp,clm2,slm2,&
    sad2,cal2,clp2,clp,dylon
  real(r8) :: ylat,ylon ! first guess output by gm2gc, input to intrp
  integer :: iter
  integer,parameter :: niter=20
  real(r8) ::                    & ! output of sub intrp
    fx,fy,fz,fv,                 & ! interpolated values of x,y,z,v
    dfxdth,dfydth,dfzdth,dfvdth, & ! derivatives of x,y,z,v wrt colatitude
    dfxdln,dfydln,dfzdln,dfvdln, & ! derivatives of x,y,z,v wrt longitude
    dfxdh ,dfydh ,dfzdh ,dfvdh     ! derivatives of x,y,z,v wrt altitude
  real(r8) ::                    & ! dummies for polar calls to intrp
    fxdum,fydum,fzdum,fvdum,     &
    dmxdth,dmydth,dmzdth,dmvdth, &
    dmxdh,dmydh,dmzdh,dmvdh
  real(r8) :: cth,sth  ! output of adpl
  character(len=5) :: edge

  ier = 0 ; gdlat = 0._r8 ; gdlon = 0._r8
!
! Determine quasi-cartesian coordinates on a unit sphere of the
! desired magnetic lat,lon in quasi-dipole coordinates.
!
  x0 = cos (qdlat*dtr) * cos (qdlon*dtr)
  y0 = cos (qdlat*dtr) * sin (qdlon*dtr)
  z0 = sin (qdlat*dtr)
!
! Initial guess:  use centered dipole, convert to geocentric coords
!
  call gm2gc (qdlat,qdlon,ylat,ylon)
!
! Iterate until (angular distance)**2 (units: radians) is within
! precise of location (qdlat,qdlon) on a unit sphere. 
! (precise is a parameter in module data)
!
  do iter=1,niter
!
! geolat,lon,alt and nglat,lon,alt are in module data (set by apex_mka)
!
    call intrp (ylat,ylon,alt, geolat,geolon,geoalt,nglat,nglon,ngalt, &
               fx,fy,fz,fv,                                            &
               dfxdth,dfydth,dfzdth,dfvdth,                            &
               dfxdln,dfydln,dfzdln,dfvdln,                            &
               dfxdh ,dfydh ,dfzdh ,dfvdh, ier)
    if (ier /= 0) then
      write(iulog,"('>>> apex_q2g error from intrp')")
      call endrun( 'qpex_q2g intrp'  )
    endif
!
!  Add-back of pseudodipole component to x,y,z,v and their derivatives.
!
    call adpl(ylat,ylon,cth,sth,fx,fy,fz,fv, &
      dfxdth,dfydth,dfzdth,dfvdth,dfxdln,dfydln,dfzdln,dfvdln)
    distlon = cos(ylat*dtr)

    if (ylat > glatmx .or. ylat < glatmn) then ! glatmx,glatmn are module data
      glatx = glatmx
      if (ylat.lt.0._r8) glatx = glatmn
      distlon = cos (glatx*dtr)
      call intrp (glatx,ylon,alt, geolat,geolon,geoalt,nglat,nglon,ngalt, &
                 fxdum,fydum,fzdum,fvdum,                                 &
                 dmxdth,dmydth,dmzdth,dmvdth,dfxdln,dfydln,dfzdln,        &
                 dfvdln,dmxdh,dmydh,dmzdh,dmvdh, ier)
      if (ier /= 0) then
        write(iulog,"('>>> apex_q2g error from polar intrp')")
        call endrun( 'qpex_q2g intrp' )
      endif

      call adpl(glatx,ylon,cth,sth,fxdum,fydum,fzdum,fvdum, &
        dmxdth,dmydth,dmzdth,dmvdth,dfxdln,dfydln,dfzdln,dfvdln)
    endif
!
! At this point, FX,FY,FZ are approximate quasi-cartesian
! coordinates on a unit sphere for the quasi-dipole coordinates
! corresponding to the geodetic coordinates YLAT, YLON.
! Normalize the vector length of (FX,FY,FZ) to unity using XNORM
! so that the resultant vector can be directly compared with the
! target vector (X0,Y0,Z0).
!
    xnorm = sqrt(fx*fx + fy*fy + fz*fz)
    xdif = fx/xnorm - x0
    ydif = fy/xnorm - y0
    zdif = fz/xnorm - z0
!
! dist2 = square of distance between normalized (fx,fy,fz) and x0,y0,z0.
!
    dist2 = xdif*xdif + ydif*ydif + zdif*zdif

    if (dist2 <= precise) then
      ier = 0
      gdlat = ylat
      gdlon = ylon
      return
    endif
!
! hgrd2* = one-half of east or north gradient of dist2 on unit sphere.
!
    hgrd2e =  (xdif*dfxdln + ydif*dfydln + zdif*dfzdln)/distlon
    hgrd2n = -(xdif*dfxdth + ydif*dfydth + zdif*dfzdth)
    hgrd2  = sqrt(hgrd2e*hgrd2e + hgrd2n*hgrd2n)
!
! angdist = magnitude of angular distance to be moved for new guess
!           of ylat, ylon.
!
    angdist = dist2/hgrd2
!
! Following spherical trigonometry moves ylat,ylon to new location,
! in direction of grad(dist2), by amount angdist.
!
    cal = -hgrd2n/hgrd2
    sal = -hgrd2e/hgrd2 
    coslm = cos(ylat*dtr)
    slm = sin(ylat*dtr)
    cad = cos(angdist) 
    sad = sin(angdist)
    slp = slm*cad + coslm*sad*cal

    clm2 = coslm*coslm
    slm2 = slm*slm
    sad2 = sad*sad
    cal2 = cal*caL
    clp2 = clm2 + slm2*sad2 - 2._r8*slm*cad*coslm*sad*cal -clm2*sad2*cal2
    clp = sqrt (max(0._r8,clp2))
    ylat = atan2(slp,clp)*rtd
!
! Restrict latitude iterations to stay within the interpolation grid
! limits, but let intrp find any longitude exceedence.  This is only
! an issue when the interpolation grid does not cover the entire
! magnetic pole region.
!
    ylat = min(ylat,geolat(nglat))
    ylat = max(ylat,geolat(1))
    dylon = atan2 (sad*sal,cad*coslm-sad*slm*cal)*rtd
    ylon = ylon + dylon
    if (ylon > geolon(nglon)) ylon = ylon - 360._r8
    if (ylon < geolon(1))     ylon = ylon + 360._r8

  enddo ! iter=1,niter

  write(iulog,"('>>> apex_q2g: ',i3,' iterations only reduced the angular')") niter
  write(iulog,"('              difference to ',f10.5,' degrees, where test criterion')") &
    sqrt(dist2)*rtd
  write(iulog,"('              is ',f10.5,' degrees.')") sqrt(precise)*rtd
  edge = '     '
  if (ylat == geolat(nglat)) edge = 'north'
  if (ylat == geolat(1))     edge = 'south'
  if (edge /= '     ') then
    write(iulog,"('Coordinates are on the ',a,' edge of the interpolation grid ')") edge
    write(iulog,"('and latitude is constrained to stay within grid limits when iterating.')") 
  endif
  ier = 1

end subroutine apex_q2g
!-----------------------------------------------------------------------
subroutine gradxyzv(alt,cth,sth, &
    dfxdth,dfydth,dfzdth,dfvdth,dfxdln,dfydln,dfzdln,dfvdln, &
    dfxdh,dfydh,dfzdh,dfvdh,gradx,grady,gradz,gradv)
!
! Calculates east,north,up components of gradients of x,y,z,v in
! geodetic coordinates.  All gradients are in inverse km.  Assumes
! flatness of 1/298.25 and equatorial radius (REQ) of 6378.16 km.
! 940803 A. D. Richmond
!
! Args:
  real(r8),intent(in) :: alt,cth,sth
  real(r8),dimension(3),intent(out) :: gradx,grady,gradz,gradv
  real(r8),intent(in) ::         &
    dfxdth,dfydth,dfzdth,dfvdth, &
    dfxdln,dfydln,dfzdln,dfvdln, &
    dfxdh,dfydh,dfzdh,dfvdh
!
! Local:
  real(r8) :: d,d2,rho,dddthod,drhodth,dzetdth,ddisdth

!
!          40680925. = req**2 (rounded off)
!          272340.   = req**2 * E2, where E2 = (2. - 1./298.25)/298.25
!                      is square of eccentricity of ellipsoid.
!
  d2 = 40680925.e0_r8 - 272340.e0_r8*cth*cth
  d = sqrt(d2)
  rho = sth*(alt + 40680925.e0_r8/d)
  dddthod = 272340.e0_r8*cth*sth/d2
  drhodth = alt*cth + (40680925.e0_r8/d)*(cth-sth*dddthod)
  dzetdth =-alt*sth - (40408585.e0_r8/d)*(sth+cth*dddthod)
  ddisdth = sqrt(drhodth*drhodth + dzetdth*dzetdth)

  gradx(1) = dfxdln/rho
  grady(1) = dfydln/rho
  gradz(1) = dfzdln/rho
  gradv(1) = dfvdln/rho

  gradx(2) = -dfxdth/ddisdth
  grady(2) = -dfydth/ddisdth
  gradz(2) = -dfzdth/ddisdth
  gradv(2) = -dfvdth/ddisdth

  gradx(3) = dfxdh
  grady(3) = dfydh
  gradz(3) = dfzdh
  gradv(3) = dfvdh

end subroutine gradxyzv
!-----------------------------------------------------------------------
subroutine grapxyzv(alt,cth,sth, &
  dfxdln,dfydln,dfzdln,dfvdln,gradx,grady,gradz,gradv)
!
! Calculates east component of gradient near pole.
!
! Args:
  real(r8),intent(in) :: alt,cth,sth
  real(r8),intent(in) :: dfxdln,dfydln,dfzdln,dfvdln
  real(r8),dimension(3),intent(out) :: gradx,grady,gradz,gradv
!
! Local:
  real(r8) :: d,d2,rho,dddthod,drhodth,dzetdth,ddisdth
!
!          40680925. = req**2 (rounded off)
!          272340.   = req**2 * E2, where E2 = (2. - 1./298.25)/298.25
!                      is square of eccentricity of ellipsoid.
!
  d2 = 40680925.e0_r8 - 272340.e0_r8*cth*cth
  d = sqrt(d2)
  rho = sth*(alt + 40680925.e0_r8/d)
  dddthod = 272340.e0_r8*cth*sth/d2
  drhodth = alt*cth + (40680925.e0_r8/d)*(cth-sth*dddthod)
  dzetdth =-alt*sth - (40408585.e0_r8/d)*(sth+cth*dddthod)
  ddisdth = sqrt(drhodth*drhodth + dzetdth*dzetdth)

  gradx(1) = dfxdln/rho
  grady(1) = dfydln/rho
  gradz(1) = dfzdln/rho
  gradv(1) = dfvdln/rho

end subroutine grapxyzv
!-----------------------------------------------------------------------
subroutine gradlpv(hr,alt,fx,fy,fz,fv,gradx,grady,gradz,gradv, &
  xlatm,xlonm,vmp,grclm,clmgrp,qdlat,rgrlp,b,clm,r3_2)
!
! Uses gradients of x,y,z,v to compute geomagnetic field and
! gradients of apex latitude, longitude.
!
! Args:
  real(r8),intent(in) :: & ! scalar inputs
    hr,                  & ! reference altitude (km)
    alt,                 & ! altitude (km)
    fx,fy,fz,fv            ! interpolated values of x,y,z,v, plus 
                           ! pseudodipole component
  real(r8),dimension(3),intent(in) :: & ! 3-component inputs
    gradx,grady,gradz,gradv ! interpolated gradients of x,y,z,v,
                            ! including pseudodipole components (east,north,up)
!
! Local:
  integer :: i
  real(r8) :: rr,r,rn,sqrror,cpm,spm,bo,rn2,x2py2,xnorm,xlp,slp,clp,grclp

  real(r8),intent(out) :: & ! scalar outputs
    xlatm,  & !  modified apex latitude (lambda_m), degrees
    xlonm,  & !  apex longitude (phi_a), degrees
    vmp,    & !  magnetic potential, in T.m.
    qdlat,  & !  quasi-dipole latitude, degrees
    clm,    & !  cos(lambda_m)
    r3_2      !  ((re + alt)/(re + hr))**(3/2)

  real(r8),dimension(3),intent(out) :: & ! 3-component outputs 
    grclm,   & ! grad(cos(lambda_m)), in km-1
    clmgrp,  & ! cos(lambda_m)*grad(phi_a), in km-1
    rgrlp,   & ! (re + alt)*grad(lambda')
    b          ! magnetic field, in nT

  xlatm=0._r8 ; xlonm=0._r8 ; vmp=0._r8 ; grclm=0._r8 ; clmgrp=0._r8 
  rgrlp = 0._r8 ; b=0._r8 ; clm=0._r8 ; r3_2=0._r8 ; qdlat=0._r8

  rr = re + hr
  r  = re + alt
  rn = r/re
  sqrror = sqrt(rr/r)
  r3_2 = 1._r8/sqrror/sqrror/sqrror
  xlonm = atan2(fy,fx)
  cpm = cos(xlonm)
  spm = sin(xlonm)
  xlonm = rtd*xlonm ! output
  bo = vp*1.e6_r8 ! vp is module data; 1.e6 converts T to nT and km-1 to m-1
  rn2 = rn*rn
  vmp = vp*fv/rn2   ! output
  b(1) = -bo*gradv(1)/rn2
  b(2) = -bo*gradv(2)/rn2
  b(3) = -bo*(gradv(3)-2._r8*fv/r)/rn2

  x2py2 = fx*fx + fy*fy
  xnorm = sqrt(x2py2 + fz*fz)
  xlp = atan2(fz,sqrt(x2py2))
  slp = sin(xlp)
  clp = cos(xlp)
  qdlat = xlp*rtd   ! output
  clm = sqrror*clp  ! output
  if (clm > 1._r8) then
    write(iulog,"('>>> gradlpv: hr=',f12.3,' alt=',f12.3)") hr,alt
    write(iulog,"('    Point lies below field line that peaks at reference height.')")
    call endrun( 'gradlpv' )
  endif
  xlatm = rtd*acos(clm)
!
!  If southern magnetic hemisphere, reverse sign of xlatm
!
  if (slp < 0._r8) xlatm = -xlatm
  do i=1,3  
    grclp = cpm*gradx(i) + spm*grady(i)
    rgrlp(i) = r*(clp*gradz(i) - slp*grclp)
    grclm(i) = sqrror*grclp
    clmgrp(i) = sqrror*(cpm*grady(i)-spm*gradx(i))
  enddo
  grclm(3) = grclm(3) - sqrror*clp/(2._r8*r)

end subroutine gradlpv
!-----------------------------------------------------------------------
subroutine basevec(hr,xlatm,grclm,clmgrp,rgrlp,b,clm,r3_2, &
                   bmag,sim,si,f,d,w,bhat,d1,d2,d3,e1,e2,e3,f1,f2)
!
! Computes base vectors and other parameters for apex coordinates.
! Vector components:  east, north, up
!
! Args:
  real(r8),intent(in) :: & ! scalar inputs
    hr,      & ! reference altitude
    xlatm,   & ! modified apex latitude (deg)
    clm,     & ! cos(lambda_m)
    r3_2       ! ((re + altitude)/(re + hr))**(3/2)

  real(r8),dimension(3),intent(in) :: & ! 3-component inputs
    grclm,   & ! grad(cos(lambda_m)), in km-1
    clmgrp,  & ! cos(lambda_m)*grad(phi_a), in km-1
    rgrlp,   & ! (re + altitude)*grad(lambda')
    b          ! ((re + altitude)/(re + hr))**(3/2)

  real(r8),intent(out) :: & ! scalar output
    bmag,    & ! magnitude of magnetic field, in nT
    sim,     & ! sin(I_m) of Richmond reference
    si,      & ! sin(I)
    f,       & ! F of Richmond reference
    d,       & ! D of Richmond reference
    w          ! W of Richmond reference

  real(r8),dimension(3),intent(out) :: & ! 3-component outputs
    bhat,             & ! unit vector along geomagnetic field direction
    d1,d2,d3,e1,e2,e3   ! base vectors of Richmond reference
  real(r8),dimension(2),intent(out) :: & ! 2-component outputs
    f1(2),f2(2)         ! base vectors of Richmond reference
!
! Local:
  integer :: i
  real(r8) :: rr,simoslm,d1db,d2db

  rr = re + hr
  simoslm = 2._r8/sqrt(4._r8 - 3._r8*clm*clm)
  sim = simoslm*sin(xlatm*dtr)
  bmag = sqrt(b(1)*b(1) + b(2)*b(2) + b(3)*b(3))
  d1db = 0._r8
  d2db = 0._r8
  do i=1,3
    bhat(i) = b(i)/bmag
    d1(i) = rr*clmgrp(i)
    d1db = d1db + d1(i)*bhat(i)
    d2(i) = rr*simoslm*grclm(i)
    d2db = d2db + d2(i)*bhat(i)
  enddo
!
! Ensure that d1,d2 are exactly perpendicular to B:
!
  do i=1,3
    d1(i) = d1(i) - d1db*bhat(i)
    d2(i) = d2(i) - d2db*bhat(i)
  enddo  
  e3(1) = d1(2)*d2(3) - d1(3)*d2(2)
  e3(2) = d1(3)*d2(1) - d1(1)*d2(3)
  e3(3) = d1(1)*d2(2) - d1(2)*d2(1) 
  d = bhat(1)*e3(1) + bhat(2)*e3(2) + bhat(3)*e3(3)
  do i=1,3
    d3(i) = bhat(i)/d
    e3(i) = bhat(i)*d ! Ensure that e3 lies along bhat.
  enddo
  e1(1) = d2(2)*d3(3) - d2(3)*d3(2)
  e1(2) = d2(3)*d3(1) - d2(1)*d3(3)
  e1(3) = d2(1)*d3(2) - d2(2)*d3(1)
  e2(1) = d3(2)*d1(3) - d3(3)*d1(2)
  e2(2) = d3(3)*d1(1) - d3(1)*d1(3)
  e2(3) = d3(1)*d1(2) - d3(2)*d1(1)
  w = rr*rr*clm*abs(sim)/(bmag*d)
  si = -bhat(3)
  f1(1) =  rgrlp(2)
  f1(2) = -rgrlp(1)
  f2(1) = -d1(2)*r3_2
  f2(2) =  d1(1)*r3_2
  f = f1(1)*f2(2) - f1(2)*f2(1)

end subroutine basevec
!-----------------------------------------------------------------------
subroutine apex_dypol(colat,elon,vp)
!
! Output args:
  real(r8),intent(out) :: &
    colat, & ! Geocentric colatitude of geomagnetic dipole north pole (deg)
    elon,  & ! East longitude of geomagnetic dipole north pole (deg)
    vp       ! Magnitude, in T.m, of dipole component of magnetic
             ! potential at geomagnetic pole and geocentric radius re
!
! Local:
  real(r8) :: gpl,ctp,stp
!
! Compute geographic colatitude and longitude of the north pole of
! earth centered dipole  
!
  gpl = sqrt( gb(2  )**2+ gb(3  )**2+ gb(4  )**2)
  ctp = gb(2  )/gpl
  stp = sqrt(1._r8 - ctp*ctp)

  colat = (acos(ctp))*rtd
  elon = atan2( gb(4  ), gb(3  ))*rtd
!           
! Compute magnitude of magnetic potential at pole, radius Re.
!      .2 = 2*(10**-4 T/gauss)*(1000 m/km) (2 comes through f0 in COFRM).
!
  vp = .2_r8*gpl*re

end subroutine apex_dypol
!-----------------------------------------------------------------------
subroutine apex_sub(date,dlat,dlon,alt,aht,alat,alon,bmag,xmag,ymag,zmag,vmp)
!
! Args:
  real(r8),intent(in) :: date
  real(r8),intent(inout) :: dlat,dlon,alt
  real(r8),intent(out) :: aht,alat,alon,bmag,xmag,ymag,zmag,vmp
!
! Local:
  real(r8) :: clatp,polon,vpol,x,y,z,xre,yre,zre
  integer :: iflag

  call cofrm(date)
  call apex_dypol(clatp,polon,vpol)
!
! colat,ctp,stp,elon,vp are in module data.
!
  colat = clatp
  ctp   = cos(clatp*dtr)
  stp   = sqrt(1._r8-ctp*ctp)

  elon  = polon
  vp    = vpol

  vmp = 0._r8
!
! Last 7 args of linapx are output:
!
  call linapx(dlat,dlon,alt, aht,alat,alon,xmag,ymag,zmag,bmag)

  xmag = xmag*1.e5_r8
  ymag = ymag*1.e5_r8
  zmag = zmag*1.e5_r8
  bmag = bmag*1.e5_r8
  call gd2cart (dlat,dlon,alt,x,y,z)
  iflag = 3
  xre = x/re ; yre = y/re ; zre = z/re
  call feldg(iflag,xre,yre,zre,bx,by,bz,vmp)

end subroutine apex_sub
!-----------------------------------------------------------------------
subroutine linapx(gdlat,glon,alt,aht,alat,alon,xmag,ymag,zmag,fmag)
!
! Input Args:
!
  real(r8),intent(inout) :: & ! These may be changed by convrt, depending on iflag
    gdlat,                  & ! latitude of starting point (deg)
    glon,                   & ! longitude of starting point (deg)
    alt                       ! height of starting point (km)
!
! Output Args:
!
  real(r8),intent(out) :: &
    aht,              & ! (Apex height+req)/req, where req is equatorial earth radius
    alat,             & ! Apex latitude (deg)
    alon,             & ! Apex longitude (deg)
    xmag,             & ! North component of magnetic field at starting point
    ymag,             & ! East component of magnetic field at starting point
    zmag,             & ! Down component of magnetic field at starting point
    fmag                ! Magnetic field magnitude at starting point
!
! Local:
!
  real(r8) :: gclat,r,singml,cgml2,rho,xlat,xlon,ht
  real(r8) :: bnrth,beast,bdown,babs,y1,y2,y3
  integer :: iflag,iapx
  integer,parameter :: maxs = 200
!
! Determine step size as a function of geomagnetic dipole
! coordinates of the starting point
!
  iflag = 2 ! gclat,r are returned
  call convrt(iflag,gdlat,alt,gclat,r,'linapx')

  singml = ctp*sin(gclat*dtr) + stp*cos(gclat*dtr)*cos((glon-elon)*dtr)
  cgml2 = max(0.25_r8,1._r8-singml*singml)
  ds = .06_r8*r/cgml2 - 370._r8 ! ds is in module data

  yapx = 0._r8 ! init (module data)
!
! Convert from geodetic to earth centered cartesian coordinates:
!
  call gd2cart(gdlat,glon,alt,y(1),y(2),y(3))
  nstp = 0
!
! Get magnetic field components to determine the direction for tracing field line:
!
  iflag = 1 
  call feldg(iflag,gdlat,glon,alt,xmag,ymag,zmag,fmag)

  sgn = sign(1._r8,-zmag)
!
! Use cartesian coordinates to get magnetic field components
! (from which gradients steer the tracing)
!
100 continue
  iflag = 2 ! module data bx,by,bz,bb are returned
  y1 = y(1)/re ; y2 = y(2)/re ; y3 = y(3)/re
  call feldg(iflag,y1,y2,y3,bx,by,bz,bb)
  nstp = nstp + 1
!
! Quit if too many steps.
!
  if (nstp >= maxs) then
    rho = sqrt(y(1)*y(1) + y(2)*y(2))
    iflag = 3 ! xlat and ht are returned
    call convrt(iflag,xlat,ht,rho,y(3),'linapx')
    xlon = rtd*atan2(y(2),y(1))
    iflag = 1
    call feldg(iflag,xlat,xlon,ht,bnrth,beast,bdown,babs)
    call dipapx(xlat,xlon,ht,bnrth,beast,bdown,aht,alon)
    alat = -sgn*rtd*acos(sqrt(1._r8/aht))
    return
  endif
!
! Find next point using adams algorithm after 7 points
!
  call itrace(iapx)
  if (iapx == 1) goto 100
!
! Maximum radius just passed.  Find apex coords
!
  call fndapx(alt,zmag,aht,alat,alon)

end subroutine linapx
!-----------------------------------------------------------------------
subroutine convrt(iflag,gdlat,alt,x1,x2,caller)
!
! Convert space point from geodetic to geocentric or vice versa.
!
! iflag = 1: Convert from geodetic to cylindrical
!   Input:  gdlat = Geodetic latitude (deg)
!           alt   = Altitude above reference ellipsoid (km)
!   Output: x1    = Distance from Earth's rotation axis (km)
!           x2    = Distance above (north of) Earth's equatorial plane (km)
!
! iflag = 2: Convert from geodetic to geocentric spherical
!   Input:  gdlat = Geodetic latitude (deg) 
!           alt   = Altitude above reference ellipsoid (km)
!   Output: x1    = Geocentric latitude (deg)
!           x2    = Geocentric distance (km)
!
! iflag = 3: Convert from cylindrical to geodetic
!   Input:  x1    = Distance from Earth's rotation axis (km)
!           x2    = Distance from Earth's equatorial plane (km)
!   Output: gdlat = Geodetic latitude (deg)
!           alt   = Altitude above reference ellipsoid (km)
!
! iflag = 4: Convert from geocentric spherical to geodetic
!   Input:  x1    = Geocentric latitude (deg)
!           x2    = Geocentric distance (km)
!   Output: gdlat = Geodetic latitude (deg)
!           alt   = Altitude above reference ellipsoid (km)
!
! Args:
  integer,intent(in) :: iflag
  real(r8),intent(inout) :: gdlat,alt
  real(r8),intent(inout) :: x1,x2
  character(len=*),intent(in) :: caller
!
! Local:
  real(r8) :: sinlat,coslat,d,z,rho,rkm,scl,gclat,ri,a2,a4,a6,a8,&
    ccl,s2cl,c2cl,s4cl,c4cl,s8cl,s6cl,dltcl,sgl
  real(r8),parameter ::                                                &
    fltnvrs = 298.25_r8                                              , &
    e2=(2._r8-1._r8/fltnvrs)/fltnvrs                                 , &
    e4=e2*e2, e6=e4*e2, e8=e4*e4                                     , &
    ome2req = (1._r8-e2)*req                                         , &
    A21 = (512._r8*E2 + 128._r8*E4 + 60._r8*E6 + 35._r8*E8)/1024._r8 , &
    A22 = (                                 E6 +        E8)/  32._r8 , &
    A23 = -3._r8*(                    4._r8*E6 +  3._r8*E8)/ 256._r8 , &
    A41 =    -(          64._r8*E4 + 48._r8*E6 + 35._r8*E8)/1024._r8 , &
    A42 =     (           4._r8*E4 +  2._r8*E6 +        E8)/  16._r8 , &
    A43 =                                        15._r8*E8 / 256._r8 , &
    A44 =                                              -E8 /  16._r8 , &
    A61 =  3._r8*(                    4._r8*E6 +  5._r8*E8)/1024._r8 , &
    A62 = -3._r8*(                          E6 +        E8)/  32._r8 , &
    A63 = 35._r8*(                    4._r8*E6 +  3._r8*E8)/ 768._r8 , &
    A81 =                                        -5._r8*E8 /2048._r8 , &
    A82 =                                        64._r8*E8 /2048._r8 , &
    A83 =                                      -252._r8*E8 /2048._r8 , &
    A84 =                                       320._r8*E8 /2048._r8

  if (iflag < 3) then ! geodetic to geocentric
!
! Compute rho,z
    sinlat = sin(gdlat*dtr)
    coslat = sqrt(1._r8-sinlat*sinlat)
    d      = sqrt(1._r8-e2*sinlat*sinlat)
    z      = (alt+ome2req/d)*sinlat
    rho    = (alt+req/d)*coslat
    x1 = rho
    x2 = z
    if (iflag == 1) return
!
! Compute gclat,rkm
    rkm   = sqrt(z*z+rho*rho)
    gclat = rtd*atan2(z,rho)
    x1 = gclat
    x2 = rkm
    return    ! iflag == 2
  endif ! iflag < 3
!
! Geocentric to geodetic:
  if (iflag == 3) then
    rho = x1
    z   = x2
    rkm = sqrt(z*z+rho*rho)
    scl = z/rkm
    gclat = asin(scl)*rtd
  elseif (iflag == 4) then
    gclat = x1
    rkm = x2
    scl = sin(gclat*dtr)
  else
    return
  endif
!
! iflag == 3 or 4:
!
  ri = req/rkm
  a2 = ri*(a21+ri*(a22+ri* a23))
  a4 = ri*(a41+ri*(a42+ri*(a43+ri*a44)))
  a6 = ri*(a61+ri*(a62+ri* a63))
  a8 = ri*(a81+ri*(a82+ri*(a83+ri*a84)))
  ccl = sqrt(1._r8-scl*scl)
  s2cl = 2._r8*scl*ccL
  c2cl = 2._r8*ccl*ccl-1._r8
  s4cl = 2._r8*s2cl*c2cl
  c4cl = 2._r8*c2cl*c2cl-1._r8
  s8cl = 2._r8*s4cl*c4cl
  s6cl = s2cl*c4cl+c2cl*s4cl
  dltcl = s2cl*a2+s4cl*a4+s6cl*a6+s8cl*a8
  gdlat = dltcl*rtd+gclat
  sgl = sin(gdlat*dtr)
  alt = rkm*cos(dltcl)-req*sqrt(1._r8-e2*sgl*sgl)

end subroutine convrt
!-----------------------------------------------------------------------
subroutine gd2cart(gdlat,glon,alt,x,y,z)
!
! Arg:
  real(r8),intent(inout) :: gdlat,alt,z
  real(r8),intent(in) :: glon
  real(r8),intent(out) :: x,y
!
! Local:
  real(r8) :: ang,rho
  integer :: iflag

  iflag = 1 ! Convert from geodetic to cylindrical (rho,z are output)
  call convrt(iflag,gdlat,alt,rho,z,'gd2cart')

  ang = glon*dtr
  x = rho*cos(ang)
  y = rho*sin(ang)

end subroutine gd2cart
!-----------------------------------------------------------------------
subroutine feldg(iflag,glat,glon,alt,bnrth,beast,bdown,babs)
!
! Compute the DGRF/IGRF field components at the point glat,glon,alt.
! cofrm must be called to establish coefficients for correct date
! prior to calling FELDG.
!
! iflag = 1:
!   Inputs:
!     glat = Latitude of point (deg)
!     glon = Longitude of point (deg)
!     alt  = Height of point (km)
!   Outputs:
!     bnrth = North component of field vector (Gauss)
!     beast = East component of field vector (Gauss)
!     bdown = Downward component of field vector (Gauss)
!     babs  = Magnitude of field vector (Gauss)  
!
! iflag = 2:
!   Inputs:
!     glat = x coordinate (in units of earth radii 6371.2 km)
!     glon = y coordinate (in units of earth radii 6371.2 km)
!     alt  = z coordinate (in units of earth radii 6371.2 km)
!   Outputs:
!     bnrth = x component of field vector (Gauss)
!     beast = y component of field vector (Gauss)
!     bdown = z component of field vector (Gauss)
!     babs  = Magnitude of field vector (Gauss)
!
! iflag = 3:
!   Inputs:
!     glat = x coordinate (in units of earth radii 6371.2 km)
!     glon = y coordinate (in units of earth radii 6371.2 km)
!     alt  = z coordinate (in units of earth radii 6371.2 km)
!   Outputs:
!     bnrth = Dummy variable
!     beast = Dummy variable
!     babs  = Legacy code had "Dummy variable" here, but its
!             set at the end if iflag==3.
!
! Args:
  integer,intent(in) :: iflag
  real(r8),intent(in)    :: glon
  real(r8),intent(out)   :: glat,alt
  real(r8),intent(out)   :: bnrth,beast,bdown,babs
!
! Local:
  integer :: i,is,ihmax,last,imax,mk,k,ih,m,il,ihm,ilm
  real(r8) :: rlat,ct,st,rlon,cp,sp,xxx,yyy,zzz,rq,f,x,y,z
  real(r8) :: xi(3),h(ncoef),g(ncoef)
  real(r8) :: s,t,bxxx,byyy,bzzz,brho

  if (iflag == 1) then
    is   = 1
    rlat = glat*dtr
    ct   = sin(rlat)
    st   = cos(rlat)
    rlon = glon*dtr
    cp   = cos(rlon)
    sp   = sin(rlon)
    call gd2cart(glat,glon,alt,xxx,yyy,zzz)
    xxx = xxx/re
    yyy = yyy/re
    zzz = zzz/re
  else
    is  = 2
    xxx = glat
    yyy = glon
    zzz = alt
  endif
  rq    = 1._r8/(xxx**2+yyy**2+zzz**2)
  xi(1) = xxx*rq
  xi(2) = yyy*rq
  xi(3) = zzz*rq
  ihmax = nmax*nmax+1
  last  = ihmax+nmax+nmax
  imax  = nmax+nmax-1
!
! Legacy code checks here to see if iflag or last call to cofrm have changed.
! For now, just do it anyway:
!
  if (iflag /= 3) then
    do i=1,last
      g(i) = gb(i) ! gb is module data from cofrm
    enddo
  else
    do i=1,last
      g(i) = gv(i) ! gv is module data from cofrm
    enddo
  endif

  do i=ihmax,last
    h(i) = g(i)
  enddo

  mk = 3
  if (imax == 1) mk = 1

  do k=1,mk,2
    i = imax
    ih = ihmax

100 continue
    il = ih-i
    f = 2._r8/dble(i-k+2)
    x = xi(1)*f
    y = xi(2)*f
    z = xi(3)*(f+f)

    i = i-2
    if (i < 1) then
      h(il) = g(il) + z*h(ih) + 2._r8*(x*h(ih+1)+y*h(ih+2))
    elseif (i == 1) then
      h(il+2) = g(il+2) + z*h(ih+2) + x*h(ih+4) - y*(h(ih+3)+h(ih))
      h(il+1) = g(il+1) + z*h(ih+1) + y*h(ih+4) + x*(h(ih+3)-h(ih))
      h(il)   = g(il)   + z*h(ih)   + 2._r8*(x*h(ih+1)+y*h(ih+2))
    else
      do m=3,i,2
        ihm = ih+m
        ilm = il+m
        h(ilm+1) = g(ilm+1)+ z*h(ihm+1) + x*(h(ihm+3)-h(ihm-1))- &
                   y*(h(ihm+2)+h(ihm-2))
        h(ilm)   = g(ilm)  + z*h(ihm)   + x*(h(ihm+2)-h(ihm-2))+ &
                   y*(h(ihm+3)+h(ihm-1))
      enddo
      h(il+2) = g(il+2) + z*h(ih+2) + x*h(ih+4) - y*(h(ih+3)+h(ih))
      h(il+1) = g(il+1) + z*h(ih+1) + y*h(ih+4) + x*(h(ih+3)-h(ih))
      h(il)   = g(il)   + z*h(ih)   + 2._r8*(x*h(ih+1)+y*h(ih+2))
    endif

    ih = il
    if (i >= k) goto 100
  enddo ! k=1,mk,2

  s = .5_r8*h(1)+2._r8*(h(2)*xi(3)+h(3)*xi(1)+h(4)*xi(2))
  t = (rq+rq)*sqrt(rq)
  bxxx = t*(h(3)-s*xxx)
  byyy = t*(h(4)-s*yyy)
  bzzz = t*(h(2)-s*zzz)
  babs = sqrt(bxxx**2+byyy**2+bzzz**2)
  if (is .eq. 1) then            ! (convert back to geodetic)
    beast = byyy*cp-bxxx*sp
    brho  = byyy*sp+bxxx*cp
    bnrth = bzzz*st-brho*ct
    bdown = -bzzz*ct-brho*st
  elseif (is .eq. 2) then        ! (leave in earth centered cartesian)
    bnrth = bxxx
    beast = byyy
    bdown = bzzz
  endif
!
! Magnetic potential computation makes use of the fact that the
! calculation of V is identical to that for r*Br, if coefficients
! in the latter calculation have been divided by (n+1) (coefficients
! GV).  Factor .1 converts km to m and gauss to tesla.
!
  if (iflag == 3) babs = (bxxx*xxx + byyy*yyy + bzzz*zzz)*re*.1_r8

end subroutine feldg
!-----------------------------------------------------------------------
subroutine dipapx(gdlat,gdlon,alt,bnorth,beast,bdown,a,alon)
!
! Compute a, alon from local magnetic field using dipole and spherical approx.
! Reference from legacy code: 940501 A. D. Richmond
!
! Input:
!   gdlat  = geodetic latitude, degrees
!   gdlon  = geodetic longitude, degrees
!   alt    = altitude, km
!   bnorth = geodetic northward magnetic field component (any units)
!   beast  = eastward magnetic field component
!   bdown  = geodetic downward magnetic field component
! Output:
!   a      = apex radius, 1 + h_A/R_eq
!   alon   = apex longitude, degrees
!     
! Algorithm: 
!   Use spherical coordinates.
!   Let GP be geographic pole.
!   Let GM be geomagnetic pole (colatitude COLAT, east longitude ELON).
!   Let G be point at GDLAT,GDLON.
!   Let E be point on sphere below apex of dipolar field line passing through G.
!   Let TD be dipole colatitude of point G, found by applying dipole formula
!     for dip angle to actual dip angle.
!   Let B be Pi plus local declination angle.  B is in the direction
!     from G to E.
!   Let TG be colatitude of G.
!   Let ANG be longitude angle from GM to G.
!   Let TE be colatitude of E.
!   Let TP be colatitude of GM.
!   Let A be longitude angle from G to E.
!   Let APANG = A + ANG
!   Let PA be geomagnetic longitude, i.e., Pi minus angle measured
!     counterclockwise from arc GM-E to arc GM-GP.
!   Let TF be arc length between GM and E.
!   Then, using notation C=cos, S=sin, COT=cot, spherical-trigonometry formulas
!     for the functions of the angles are as shown below.  Note: STFCPA,
!     STFSPA are sin(TF) times cos(PA), sin(PA), respectively.
!
  real(r8),intent(in)  :: gdlat,gdlon,alt,bnorth,beast,bdown
  real(r8),intent(out) :: a,alon
!
! Local:
  real(r8) :: bhor,std,ctd,sb,cb,ctg,stg,ang,sang,cang,cte,ste,sa,ca, &
              cottd,capang,sapang,stfcpa,stfspa,ha,r

  bhor = sqrt(bnorth*bnorth + beast*beast)
  if (bhor == 0._r8) then
    alon = 0._r8
    a = 1.e34_r8
    return
  endif

  cottd = bdown*.5_r8/bhor
  std = 1._r8/sqrt(1._r8+cottd*cottd)
  ctd = cottd*std
  sb = -beast/bhor
  cb = -bnorth/bhor
  ctg = sin(gdlat*dtr)
  stg = cos(gdlat*dtr)
  ang = (gdlon-elon)*dtr
  sang = sin(ang)
  cang = cos(ang)
  cte = ctg*std + stg*ctd*cb
  ste = sqrt(1._r8 - cte*cte)
  sa = sb*ctd/ste
  ca = (std*stg - ctd*ctg*cb)/ste
  capang = ca*cang - sa*sang
  sapang = ca*sang + sa*cang
  stfcpa = ste*ctp*capang - cte*stp
  stfspa = sapang*ste
  alon = atan2(stfspa,stfcpa)*rtd
  r = alt + re
  ha = alt + r*cottd*cottd
  a = 1._r8 + ha/req

end subroutine dipapx
!-----------------------------------------------------------------------
subroutine itrace(iapx)
  save
!
! Uses 4-point ADAMS formula after initialization.
! First 7 iterations advance point by 3 steps.
!
! y(3), yp(3), yapx(3,3), sgn and nstp are in module data
! yploc(3,4) is local
!
! Arg:
  integer,intent(out) :: iapx
!
! Local:
  integer :: i,j
  real(r8) :: yploc(3,4) ! local yp (i.e., not module data yp)
  real(r8) :: term,d2,d6,d12,d24,rc,rp

  iapx = 1
!
! Field line is defined by the following differential equations
! in cartesian coordinates.
! (yapx,yp,y are module data)
!
  yploc(1,4) = sgn*bx/bb 
  yploc(2,4) = sgn*by/bb 
  yploc(3,4) = sgn*bz/bb 

  if (nstp > 7) then
    do i=1,3
      yapx(i,1) = yapx(i,2)
      yapx(i,2) = y(i)
      yp(i) = y(i)
      term = 55._r8*yploc(i,4)-59._r8*yploc(i,3)+37._r8*yploc(i,2)-9._r8*yploc(i,1)
      y(i) = yp(i) + d24*term
      yapx(i,3) = y(i)
      do j=1,3
        yploc(i,j) = yploc(i,j+1)
      enddo
    enddo
    rc = rdus ( y(1),  y(2),  y(3))
    rp = rdus (yp(1), yp(2), yp(3))
    if (rc < rp) iapx=2
    return
  endif

  do i=1,3
    select case (nstp)
      case (1)
        d2  = ds/2._r8
        d6  = ds/6._r8
        d12 = ds/12._r8
        d24 = ds/24._r8
        yploc(i,1)= yploc(i,4)
        yp(i)     = y(i)
        yapx(i,1) = y(i)
        y(i) = yp(i) + ds*yploc(i,1)
      case (2)
        yploc(i,2) = yploc(i,4)
        y(i) = yp(i) + d2*(yploc(i,2)+yploc(i,1))
      case (3)
        y(i) = yp(i) + d6*(2._r8*yploc(i,4)+yploc(i,2)+3._r8*yploc(i,1))
      case (4)
        yploc(i,2) = yploc(i,4)
        yapx(i,2)  = y(i)
        yp(i)      = y(i)
        y(i)       = yp(i) + d2*(3._r8*yploc(i,2)-yploc(i,1))
      case (5)
        y(i) = yp(i) + d12*(5._r8*yploc(i,4)+8._r8*yploc(i,2)-yploc(i,1))
      case (6)
        yploc(i,3) = yploc(i,4)
        yp(i)      = y(i)
        yapx(i,3)  = y(i)
        y(i)       = yp(i) + d12*(23._r8*yploc(i,3)-16._r8*yploc(i,2)+5._r8*yploc(i,1))
      case (7)
        yapx(i,1) = yapx(i,2)
        yapx(i,2) = yapx(i,3)
        y(i) = yp(i) + d24*(9._r8*yploc(i,4)+19._r8*yploc(i,3)-5._r8*yploc(i,2)+yploc(i,1))
        yapx(i,3) = y(i)
      case default
        write(iulog,"('>>> itrace: unresolved case nstp=',i4)") nstp
        call endrun( 'itrace' )
    end select
  enddo
!
! Signal if apex passed:
!
  if (nstp == 6 .or. nstp == 7) then
    rc = rdus( yapx(1,3), yapx(2,3), yapx(3,3))
    rp = rdus( yapx(1,2), yapx(2,2), yapx(3,2))
    if (rc < rp) iapx=2
  endif

end subroutine itrace
!-----------------------------------------------------------------------
real(r8) function rdus(d,e,f)
  real(r8),intent(in) :: d,e,f
  rdus = sqrt(d**2 + e**2 + f**2)
end function rdus
!-----------------------------------------------------------------------
subroutine fndapx(alt,zmag,a,alat,alon)
!
! Find apex coords once tracing has signalled that the apex has been passed.
!
! Args:
  real(r8),intent(in) :: alt,zmag
  real(r8),intent(out) :: a,alat,alon
!
! Local:
  integer :: i,iflag_convrt, iflag_feldg
  real(r8) :: z(3),ht(3),yloc(3),gdlt,gdln,x,ydum,f,rho,xinter,rasq,xlon,ang,&
    cang,sang,r,cte,ste,stfcpa,stfspa
!
! Get geodetic field components.
!
  iflag_feldg = 1
  iflag_convrt = 3
  do i=1,3
    rho = sqrt(yapx(1,i)**2+yapx(2,i)**2)
    call convrt(iflag_convrt,gdlt,ht(i),rho,yapx(3,i),'fndapx')
    gdln = rtd*atan2(yapx(2,i),yapx(1,i))
    call feldg(iflag_feldg,gdlt,gdln,ht(i),x,ydum,z(i),f)
  enddo 
!
! Find cartesian coordinates at dip equator by interpolation
!
  do i=1,3
    call fint(z(1),z(2),z(3),yapx(i,1),yapx(i,2),yapx(i,3),0._r8,yloc(i))
  enddo
!
! Find apex height by interpolation
!
  call fint(z(1),z(2),z(3),ht(1),ht(2),ht(3),0._r8,xinter)
!
! Ensure that XINTER is not less than original starting altitude:
  xinter = max(alt,xinter)
  a = (req+xinter)/req
!
! Find apex coordinates , giving alat sign of dip at starting point.  
! Alon is the value of the geomagnetic longitude at the apex.
!
  if (a < 1._r8) then
    write(iulog,"('>>> fndapx: a=',e12.4,' < 1.')") a
    call endrun( 'fndapx' )
  endif

  rasq = rtd*acos(sqrt(1._r8/a))
  alat = sign(rasq,zmag)
!
! Algorithm for ALON:
!   Use spherical coordinates.
!   Let GP be geographic pole.
!   Let GM be geomagnetic pole (colatitude COLAT, east longitude ELON).
!   Let XLON be longitude of apex.
!   Let TE be colatitude of apex.
!   Let ANG be longitude angle from GM to apex.
!   Let TP be colatitude of GM.
!   Let TF be arc length between GM and apex.
!   Let PA = ALON be geomagnetic longitude, i.e., Pi minus angle measured 
!     counterclockwise from arc GM-apex to arc GM-GP.
!   Then, using notation C=cos, S=sin, spherical-trigonometry formulas
!     for the functions of the angles are as shown below.  Note: STFCPA,
!     STFSPA are sin(TF) times cos(PA), sin(PA), respectively. 
!
  xlon = atan2(yloc(2),yloc(1))
  ang  = xlon-elon*dtr
  cang = cos(ang)
  sang = sin(ang)
  r    = sqrt(yloc(1)**2+yloc(2)**2+yloc(3)**2)
  cte  = yloc(3)/r
  ste  = sqrt(1._r8-cte*cte)
  stfcpa = ste*ctp*cang - cte*stp
  stfspa = sang*ste
  alon = atan2(stfspa,stfcpa)*rtd

end subroutine fndapx
!-----------------------------------------------------------------------
subroutine fint(a1,a2,a3,a4,a5,a6,a7,result)
!
! Second degree interpolation
!
! Args:
  real(r8),intent(in) :: a1,a2,a3,a4,a5,a6,a7
  real(r8),intent(out) :: result

  result = ((a2-a3)*(a7-a2)*(a7-a3)*a4-(a1-a3)*(a7-a1)*(a7-a3)*a5+ &
    (a1-a2)*(a7-a1)*(a7-a2)*a6)/((a1-a2)*(a1-a3)*(a2-a3))

end subroutine fint
!-----------------------------------------------------------------------
subroutine gm2gc(gmlat,gmlon,gclat,gclon)
!
! Args:
  real(r8),intent(in)  :: gmlat,gmlon
  real(r8),intent(out) :: gclat,gclon
!
! Local:
  real(r8) :: ylat,ylon,stm,ctm,ctc

  stm = cos(gmlat*dtr)
  ctm = sin(gmlat*dtr)
  ctc = ctp*ctm - stp*stm*cos(gmlon*dtr) ! ctp,stp are module data
  ctc = min(ctc,1._r8)
  ctc = max(ctc,-1._r8)
  gclat = asin(ctc)*rtd
  gclon = atan2(stp*stm*sin(gmlon*dtr),ctm-ctp*ctc)
!
! elon is in module data, and was set by dypol (called from apex_mka)
!
  gclon = gclon*rtd + elon 
  if (gclon < -180._r8) gclon = gclon + 360._r8

end subroutine gm2gc
!-----------------------------------------------------------------------
subroutine intrp(glat,glon,alt, gplat,gplon,gpalt, nlat,nlon,nalt, &
                 fx,fy,fz,fv,                                      &
                 dfxdth,dfydth,dfzdth,dfvdth,                      &
                 dfxdln,dfydln,dfzdln,dfvdln,                      &
                 dfxdh ,dfydh ,dfzdh ,dfvdh, ier)
!
! Args:
!
  real(r8),intent(in)    :: glat,glon,alt
  integer,intent(in) :: nlat,nlon,nalt
  real(r8),intent(in)    :: gplat(nlat),gplon(nlon),gpalt(nalt)
  real(r8),intent(out)   ::          &
    fx,fy,fz,fv,                 &
    dfxdth,dfydth,dfzdth,dfvdth, &
    dfxdln,dfydln,dfzdln,dfvdln, &
    dfxdh ,dfydh ,dfzdh ,dfvdh
  integer,intent(out) :: ier
!
! Local:
!
  integer :: i,j,k,i0,j0,k0
  real(r8) :: glonloc,xi,dlon,yj,dlat,hti,diht,zk,fac,omfac
  real(r8) :: dfxdn,dfxde,dfxdd, &
          dfydn,dfyde,dfydd, &
          dfzdn,dfzde,dfzdd, &
          dfvdn,dfvde,dfvdd, &
          dmf,dmdfdn,dmdfde,dmdfdd

  ier = 0
  glonloc = glon
  if (glonloc < gplon(1))    glonloc = glonloc + 360._r8
  if (glonloc > gplon(nlon)) glonloc = glonloc - 360._r8
!
  i0 = 0
  do i=1,nlat-1
    if (glat >= gplat(i).and.glat <= gplat(i+1)) then
      i0 = i
      dlat = gplat(i+1)-gplat(i)
      xi = (glat - gplat(i)) / dlat
      exit 
    endif
  enddo
  if (i0==0) then
    write(iulog,"('>>> intrp: could not bracket glat=',f9.3,' in gplat=',/,(6f9.2))") &
      glat,gplat
    ier = 1
    return 
  endif

  j0 = 0
  do j=1,nlon-1
    if (glon >= gplon(j).and.glon <= gplon(j+1)) then
      j0 = j
      dlon = gplon(j+1)-gplon(j)
      yj = (glon - gplon(j)) / dlon
      exit 
    endif
  enddo
  if (j0==0) then
    write(iulog,"('>>> intrp: could not bracket glon=',f9.3,' in gplon=',/,(6f9.2))") &
      glon,gplon
    ier = 1
    return 
  endif

  k0 = 0
  do k=1,nalt-1
    if (alt >= gpalt(k).and.alt <= gpalt(k+1)) then
      k0 = k
      hti = re/(re+alt)
      diht = re/(re+gpalt(k+1)) - re/(re+gpalt(k))
      zk = (hti - re/(re+gpalt(k))) / diht
      exit 
    endif
  enddo
  if (k0==0) then
    write(iulog,"('>>> intrp: could not bracket alt=',f12.3,' in gpalt=',/,(6f12.2))") &
      alt,gpalt
    ier = 1
    return 
  endif

  call trilin(xarray(i0:i0+1,j0:j0+1,k0:k0+1),nlat,nlon,xi,yj,zk,fx,dfxdn,dfxde,dfxdd)
  dfxdth = -dfxdn*rtd/dlat
  dfxdln =  dfxde*rtd/dlon
  dfxdh  = -hti*hti*dfxdd/(re*diht)

  call trilin(yarray(i0:i0+1,j0:j0+1,k0:k0+1),nlat,nlon,xi,yj,zk,fy,dfydn,dfyde,dfydd)
  dfydth = -dfydn*rtd/dlat
  dfydln =  dfyde*rtd/dlon
  dfydh  = -hti*hti*dfydd/(re*diht)

  call trilin(zarray(i0:i0+1,j0:j0+1,k0:k0+1),nlat,nlon,xi,yj,zk,fz,dfzdn,dfzde,dfzdd)
  dfzdth = -dfzdn*rtd/dlat
  dfzdln =  dfzde*rtd/dlon
  dfzdh  = -hti*hti*dfzdd/(re*diht)

  call trilin(varray(i0:i0+1,j0:j0+1,k0:k0+1),nlat,nlon,xi,yj,zk,fv,dfvdn,dfvde,dfvdd)
  dfvdth = -dfvdn*rtd/dlat
  dfvdln =  dfvde*rtd/dlon
  dfvdh  = -hti*hti*dfvdd/(re*diht)

  if (nlat < 3) return
!
! Improve calculation of longitudinal derivatives near poles
!
  if (glat < dlat-90._r8) then
    fac = .5_r8*xi
    omfac = 1._r8 - fac
    xi = xi - 1._r8
    i0 = i0 + 1
    call trilin (xarray(i0:i0+1,j0:j0+1,k0:k0+1),nlat,nlon,xi,yj,zk,dmf,dmdfdn,dmdfde,dmdfdd)
    dfxdln = dfxdln*omfac + fac*dmdfde*rtd/dlon
    call trilin (yarray(i0:i0+1,j0:j0+1,k0:k0+1),nlat,nlon,xi,yj,zk,dmf,dmdfdn,dmdfde,dmdfdd)
    dfydln = dfydln*omfac + fac*dmdfde*rtd/dlon
    call trilin (varray(i0:i0+1,j0:j0+1,k0:k0+1),nlat,nlon,xi,yj,zk,dmf,dmdfdn,dmdfde,dmdfdd)
    dfvdln = dfvdln*omfac + fac*dmdfde*rtd/dlon
  endif

  if (glat > 90._r8-dlat) then
    fac = .5_r8*(1._r8-xi)
    omfac = 1._r8 - fac
    xi = xi + 1._r8
    i0 = i0 - 1
    call trilin (xarray(i0:i0+1,j0:j0+1,k0:k0+1),nlat,nlon,xi,yj,zk,dmf,dmdfdn,dmdfde,dmdfdd)
    dfxdln = dfxdln*omfac + fac*dmdfde*rtd/dlon
    call trilin (yarray(i0:i0+1,j0:j0+1,k0:k0+1),nlat,nlon,xi,yj,zk,dmf,dmdfdn,dmdfde,dmdfdd)
    dfydln = dfydln*omfac + fac*dmdfde*rtd/dlon
    call trilin (varray(i0:i0+1,j0:j0+1,k0:k0+1),nlat,nlon,xi,yj,zk,dmf,dmdfdn,dmdfde,dmdfdd)
    dfvdln = dfvdln*omfac + fac*dmdfde*rtd/dlon
  endif

end subroutine intrp
!-----------------------------------------------------------------------
subroutine trilin(u,nlat,nlon,xi,yj,zk,fu,dfudx,dfudy,dfudz)
!
! Args:
  integer,intent(in) :: &
    nlat,               & ! first dimension of u from calling routine
    nlon                  ! second dimension of u from calling routine
  real(r8),intent(in) :: &
    u(1:2,1:2,1:2),      & ! u(1,1,1) is address of lower corner of interpolation box
    xi,  & ! fractional distance across box in x direction
    yj,  & ! fractional distance across box in y direction
    zk     ! fractional distance across box in z direction
  real(r8),intent(out)   :: &
    fu,                 & ! interpolated value of u
    dfudx,              & ! interpolated derivative of u with respect to i (x direction)
    dfudy,              & ! interpolated derivative of u with respect to j (y direction)
    dfudz                 ! interpolated derivative of u with respect to k (z direction)
!
! Local:
  real(r8) :: omxi,omyj,omzk

! write(iulog,"('Enter trilin: xi,yj,zk=',3e12.4)") xi,yj,zk
! write(iulog,"('Enter trilin: u(1,1,1),u(1,2,1),u(1,1,2),u(1,2,2)=',4e12.4)") &
!                          u(1,1,1),u(1,2,1),u(1,1,2),u(1,2,2)
! write(iulog,"('Enter trilin: u(2,1,1),u(2,2,1),u(2,1,2),u(2,2,2)=',4e12.4)") &
!                          u(2,1,1),u(2,2,1),u(2,1,2),u(2,2,2)

  omxi = 1._r8 - xi
  omyj = 1._r8 - yj
  omzk = 1._r8 - zk

  fu = u(1,1,1)*omxi*omyj*omzk &
     + u(2,1,1)*xi*omyj*omzk   &
     + u(1,2,1)*omxi*yj*omzk   &
     + u(1,1,2)*omxi*omyj*zk   &
     + u(2,2,1)*xi*yj*omzk     &
     + u(2,1,2)*xi*omyj*zk     &
     + u(1,2,2)*omxi*yj*zk     &
     + u(2,2,2)*xi*yj*zk

  dfudx = (u(2,1,1)-u(1,1,1))*omyj*omzk &
        + (u(2,2,1)-u(1,2,1))*yj*omzk   &
        + (u(2,1,2)-u(1,1,2))*omyj*zk   &
        + (u(2,2,2)-u(1,2,2))*yj*zk
  dfudy = (u(1,2,1)-u(1,1,1))*omxi*omzk &
        + (u(2,2,1)-u(2,1,1))*xi*omzk   &
        + (u(1,2,2)-u(1,1,2))*omxi*zk   &
        + (u(2,2,2)-u(2,1,2))*xi*zk
  dfudz = (u(1,1,2)-u(1,1,1))*omxi*omyj &
        + (u(2,1,2)-u(2,1,1))*xi*omyj   &
        + (u(1,2,2)-u(1,2,1))*omxi*yj   &
        + (u(2,2,2)-u(2,2,1))*xi*yj

end subroutine trilin
!-----------------------------------------------------------------------
subroutine adpl(glat,glon,cth,sth,fx,fy,fz,fv, &
                dfxdth,dfydth,dfzdth,dfvdth,dfxdln,dfydln,dfzdln,dfvdln)
!
!  Add-back of pseudodipole component to x,y,z,v and their derivatives.
!
! Args:
  real(r8),intent(in)      :: glat,glon
  real(r8),intent(out)     :: cth,sth
  real(r8),intent(inout)   ::    &
     fx,fy,fz,fv,                &
    dfxdth,dfydth,dfzdth,dfvdth, &
    dfxdln,dfydln,dfzdln,dfvdln
!
! Local:
  real(r8) :: cph,sph,ctm

  cph = cos((glon-elon)*dtr)
  sph = sin((glon-elon)*dtr)
  cth = sin(glat*dtr)
  sth = cos(glat*dtr)
  ctm = ctp*cth + stp*sth*cph
  fx = fx + sth*ctp*cph - cth*stp
  fy = fy + sth*sph
  fz = fz + ctm
  fv = fv - ctm

  dfxdth = dfxdth + ctp*cth*cph + stp*sth
  dfydth = dfydth + cth*sph 
  dfzdth = dfzdth - ctp*sth + stp*cth*cph
  dfvdth = dfvdth + ctp*sth - stp*cth*cph

  dfxdln = dfxdln - ctp*sth*sph
  dfydln = dfydln + sth*cph
  dfzdln = dfzdln - stp*sth*sph
  dfvdln = dfvdln + stp*sth*sph

end subroutine adpl
!-----------------------------------------------------------------------
subroutine setmiss(xmiss,xlatm,alon,vmp,b,bmag,be3,sim,si,f,d,w, &
  bhat,d1,d2,d3,e1,e2,e3,f1,f2)
!
! Args:
  real(r8),intent(in)  :: xmiss
  real(r8),intent(out) :: xlatm,alon,vmp,bmag,be3,sim,si,f,d,w
  real(r8),dimension(3),intent(out) :: bhat,d1,d2,d3,e1,e2,e3,b
  real(r8),dimension(2),intent(out) :: f1,f2

  xlatm = xmiss
  alon  = xmiss
  vmp   = xmiss
  bmag  = xmiss
  be3   = xmiss
  sim   = xmiss
  si    = xmiss
  f     = xmiss
  d     = xmiss
  w     = xmiss
  bhat  = xmiss
  d1    = xmiss
  d2    = xmiss
  d3    = xmiss
  e1    = xmiss
  e2    = xmiss
  e3    = xmiss
  b     = xmiss
  f1    = xmiss
  f2    = xmiss

end subroutine setmiss
!-----------------------------------------------------------------------
subroutine cofrm(date)
  implicit none
!
! Input arg:
  real(r8),intent(in) :: date
!
! Local:
  integer :: m,n,i,l,ll,lm,nmx,nc,kmx,k,mm,nn
  real(r8) :: t,one,tc,r,f,f0
  integer,parameter :: n1=120, n2=195, isv=0
  integer,parameter :: &
    ncn1=19, & ! number of coefficients dimensioned n1
    ncn2=5     ! number of coefficients dimensioned n2
  integer,parameter :: ngh = n1*ncn1 + n2*ncn2 + 1 ! not sure why the extra +1
  real(r8),save :: g1(n1,ncn1), g2(n2,ncn2), gh(ngh)
  real(r8),parameter :: alt = 0._r8

  if (date < 1900._r8 .or. date > 2020._r8) then
    write(iulog,"('>>> cofrm: date=',f8.2,' Date must be >= 1900 and <= 2020')") date
    call endrun( 'cofrm' )
  endif
  if (date > 2015._r8) then
    write(iulog,"('>>> WARNING cofrm:')")
    write(iulog,"(/,'   This version of IGRF is intended for use up to ')")
    write(iulog,"('     2015. Values for ',f9.3,' will be computed but')") date
    write(iulog,"('     may be of reduced accuracy.',/)")
  endif

  g1(:,1) = (/                                               & ! 1900
   -31543._r8,-2298._r8, 5922._r8, -677._r8, 2905._r8,-1061._r8,  924._r8, 1121._r8, &
     1022._r8,-1469._r8, -330._r8, 1256._r8,    3._r8,  572._r8,  523._r8,  876._r8, &
      628._r8,  195._r8,  660._r8,  -69._r8, -361._r8, -210._r8,  134._r8,  -75._r8, &
     -184._r8,  328._r8, -210._r8,  264._r8,   53._r8,    5._r8,  -33._r8,  -86._r8, &
     -124._r8,  -16._r8,    3._r8,   63._r8,   61._r8,   -9._r8,  -11._r8,   83._r8, &
     -217._r8,    2._r8,  -58._r8,  -35._r8,   59._r8,   36._r8,  -90._r8,  -69._r8, &
       70._r8,  -55._r8,  -45._r8,    0._r8,  -13._r8,   34._r8,  -10._r8,  -41._r8, &
       -1._r8,  -21._r8,   28._r8,   18._r8,  -12._r8,    6._r8,  -22._r8,   11._r8, &
        8._r8,    8._r8,   -4._r8,  -14._r8,   -9._r8,    7._r8,    1._r8,  -13._r8, &
        2._r8,    5._r8,   -9._r8,   16._r8,    5._r8,   -5._r8,    8._r8,  -18._r8, &
        8._r8,   10._r8,  -20._r8,    1._r8,   14._r8,  -11._r8,    5._r8,   12._r8, &
       -3._r8,    1._r8,   -2._r8,   -2._r8,    8._r8,    2._r8,   10._r8,   -1._r8, &
       -2._r8,   -1._r8,    2._r8,   -3._r8,   -4._r8,    2._r8,    2._r8,    1._r8, &
       -5._r8,    2._r8,   -2._r8,    6._r8,    6._r8,   -4._r8,    4._r8,    0._r8, &
        0._r8,   -2._r8,    2._r8,    4._r8,    2._r8,    0._r8,    0._r8,   -6._r8  /) 

  g1(:,2) = (/                                               & ! 1905
   -31464._r8,-2298._r8, 5909._r8, -728._r8, 2928._r8,-1086._r8, 1041._r8, 1065._r8, &
     1037._r8,-1494._r8, -357._r8, 1239._r8,   34._r8,  635._r8,  480._r8,  880._r8, &
      643._r8,  203._r8,  653._r8,  -77._r8, -380._r8, -201._r8,  146._r8,  -65._r8, &
     -192._r8,  328._r8, -193._r8,  259._r8,   56._r8,   -1._r8,  -32._r8,  -93._r8, &
     -125._r8,  -26._r8,   11._r8,   62._r8,   60._r8,   -7._r8,  -11._r8,   86._r8, &
     -221._r8,    4._r8,  -57._r8,  -32._r8,   57._r8,   32._r8,  -92._r8,  -67._r8, &
       70._r8,  -54._r8,  -46._r8,    0._r8,  -14._r8,   33._r8,  -11._r8,  -41._r8, &
        0._r8,  -20._r8,   28._r8,   18._r8,  -12._r8,    6._r8,  -22._r8,   11._r8, &
        8._r8,    8._r8,   -4._r8,  -15._r8,   -9._r8,    7._r8,    1._r8,  -13._r8, &
        2._r8,    5._r8,   -8._r8,   16._r8,    5._r8,   -5._r8,    8._r8,  -18._r8, &
        8._r8,   10._r8,  -20._r8,    1._r8,   14._r8,  -11._r8,    5._r8,   12._r8, &
       -3._r8,    1._r8,   -2._r8,   -2._r8,    8._r8,    2._r8,   10._r8,    0._r8, &
       -2._r8,   -1._r8,    2._r8,   -3._r8,   -4._r8,    2._r8,    2._r8,    1._r8, &
       -5._r8,    2._r8,   -2._r8,    6._r8,    6._r8,   -4._r8,    4._r8,    0._r8, &
        0._r8,   -2._r8,    2._r8,    4._r8,    2._r8,    0._r8,    0._r8,   -6._r8  /)

  g1(:,3) = (/                                               & ! 1910
   -31354._r8,-2297._r8, 5898._r8, -769._r8, 2948._r8,-1128._r8, 1176._r8, 1000._r8, &
     1058._r8,-1524._r8, -389._r8, 1223._r8,   62._r8,  705._r8,  425._r8,  884._r8, &
      660._r8,  211._r8,  644._r8,  -90._r8, -400._r8, -189._r8,  160._r8,  -55._r8, &
     -201._r8,  327._r8, -172._r8,  253._r8,   57._r8,   -9._r8,  -33._r8, -102._r8, &
     -126._r8,  -38._r8,   21._r8,   62._r8,   58._r8,   -5._r8,  -11._r8,   89._r8, &
     -224._r8,    5._r8,  -54._r8,  -29._r8,   54._r8,   28._r8,  -95._r8,  -65._r8, &
       71._r8,  -54._r8,  -47._r8,    1._r8,  -14._r8,   32._r8,  -12._r8,  -40._r8, &
        1._r8,  -19._r8,   28._r8,   18._r8,  -13._r8,    6._r8,  -22._r8,   11._r8, &
        8._r8,    8._r8,   -4._r8,  -15._r8,   -9._r8,    6._r8,    1._r8,  -13._r8, &
        2._r8,    5._r8,   -8._r8,   16._r8,    5._r8,   -5._r8,    8._r8,  -18._r8, &
        8._r8,   10._r8,  -20._r8,    1._r8,   14._r8,  -11._r8,    5._r8,   12._r8, &
       -3._r8,    1._r8,   -2._r8,   -2._r8,    8._r8,    2._r8,   10._r8,    0._r8, &
       -2._r8,   -1._r8,    2._r8,   -3._r8,   -4._r8,    2._r8,    2._r8,    1._r8, &
       -5._r8,    2._r8,   -2._r8,    6._r8,    6._r8,   -4._r8,    4._r8,    0._r8, &
        0._r8,   -2._r8,    2._r8,    4._r8,    2._r8,    0._r8,    0._r8,   -6._r8  /)

  g1(:,4) = (/                                               & ! 1915
   -31212._r8,-2306._r8, 5875._r8, -802._r8, 2956._r8,-1191._r8, 1309._r8,  917._r8, &
     1084._r8,-1559._r8, -421._r8, 1212._r8,   84._r8,  778._r8,  360._r8,  887._r8, &
      678._r8,  218._r8,  631._r8, -109._r8, -416._r8, -173._r8,  178._r8,  -51._r8, &
     -211._r8,  327._r8, -148._r8,  245._r8,   58._r8,  -16._r8,  -34._r8, -111._r8, &
     -126._r8,  -51._r8,   32._r8,   61._r8,   57._r8,   -2._r8,  -10._r8,   93._r8, &
     -228._r8,    8._r8,  -51._r8,  -26._r8,   49._r8,   23._r8,  -98._r8,  -62._r8, &
       72._r8,  -54._r8,  -48._r8,    2._r8,  -14._r8,   31._r8,  -12._r8,  -38._r8, &
        2._r8,  -18._r8,   28._r8,   19._r8,  -15._r8,    6._r8,  -22._r8,   11._r8, &
        8._r8,    8._r8,   -4._r8,  -15._r8,   -9._r8,    6._r8,    2._r8,  -13._r8, &
        3._r8,    5._r8,   -8._r8,   16._r8,    6._r8,   -5._r8,    8._r8,  -18._r8, &
        8._r8,   10._r8,  -20._r8,    1._r8,   14._r8,  -11._r8,    5._r8,   12._r8, &
       -3._r8,    1._r8,   -2._r8,   -2._r8,    8._r8,    2._r8,   10._r8,    0._r8, &
       -2._r8,   -1._r8,    2._r8,   -3._r8,   -4._r8,    2._r8,    2._r8,    1._r8, &
       -5._r8,    2._r8,   -2._r8,    6._r8,    6._r8,   -4._r8,    4._r8,    0._r8, &
        0._r8,   -2._r8,    1._r8,    4._r8,    2._r8,    0._r8,    0._r8,   -6._r8  /)

  g1(:,5) = (/                                               & ! 1920
   -31060._r8,-2317._r8, 5845._r8, -839._r8, 2959._r8,-1259._r8, 1407._r8,  823._r8, &
     1111._r8,-1600._r8, -445._r8, 1205._r8,  103._r8,  839._r8,  293._r8,  889._r8, &
      695._r8,  220._r8,  616._r8, -134._r8, -424._r8, -153._r8,  199._r8,  -57._r8, &
     -221._r8,  326._r8, -122._r8,  236._r8,   58._r8,  -23._r8,  -38._r8, -119._r8, &
     -125._r8,  -62._r8,   43._r8,   61._r8,   55._r8,    0._r8,  -10._r8,   96._r8, &
     -233._r8,   11._r8,  -46._r8,  -22._r8,   44._r8,   18._r8, -101._r8,  -57._r8, &
       73._r8,  -54._r8,  -49._r8,    2._r8,  -14._r8,   29._r8,  -13._r8,  -37._r8, &
        4._r8,  -16._r8,   28._r8,   19._r8,  -16._r8,    6._r8,  -22._r8,   11._r8, &
        7._r8,    8._r8,   -3._r8,  -15._r8,   -9._r8,    6._r8,    2._r8,  -14._r8, &
        4._r8,    5._r8,   -7._r8,   17._r8,    6._r8,   -5._r8,    8._r8,  -19._r8, &
        8._r8,   10._r8,  -20._r8,    1._r8,   14._r8,  -11._r8,    5._r8,   12._r8, &
       -3._r8,    1._r8,   -2._r8,   -2._r8,    9._r8,    2._r8,   10._r8,    0._r8, &
       -2._r8,   -1._r8,    2._r8,   -3._r8,   -4._r8,    2._r8,    2._r8,    1._r8, &
       -5._r8,    2._r8,   -2._r8,    6._r8,    6._r8,   -4._r8,    4._r8,    0._r8, &
        0._r8,   -2._r8,    1._r8,    4._r8,    3._r8,    0._r8,    0._r8,   -6._r8  /)

  g1(:,6) = (/                                               & ! 1925
   -30926._r8,-2318._r8, 5817._r8, -893._r8, 2969._r8,-1334._r8, 1471._r8,  728._r8, &
     1140._r8,-1645._r8, -462._r8, 1202._r8,  119._r8,  881._r8,  229._r8,  891._r8, &
      711._r8,  216._r8,  601._r8, -163._r8, -426._r8, -130._r8,  217._r8,  -70._r8, &
     -230._r8,  326._r8,  -96._r8,  226._r8,   58._r8,  -28._r8,  -44._r8, -125._r8, &
     -122._r8,  -69._r8,   51._r8,   61._r8,   54._r8,    3._r8,   -9._r8,   99._r8, &
     -238._r8,   14._r8,  -40._r8,  -18._r8,   39._r8,   13._r8, -103._r8,  -52._r8, &
       73._r8,  -54._r8,  -50._r8,    3._r8,  -14._r8,   27._r8,  -14._r8,  -35._r8, &
        5._r8,  -14._r8,   29._r8,   19._r8,  -17._r8,    6._r8,  -21._r8,   11._r8, &
        7._r8,    8._r8,   -3._r8,  -15._r8,   -9._r8,    6._r8,    2._r8,  -14._r8, &
        4._r8,    5._r8,   -7._r8,   17._r8,    7._r8,   -5._r8,    8._r8,  -19._r8, &
        8._r8,   10._r8,  -20._r8,    1._r8,   14._r8,  -11._r8,    5._r8,   12._r8, &
       -3._r8,    1._r8,   -2._r8,   -2._r8,    9._r8,    2._r8,   10._r8,    0._r8, &
       -2._r8,   -1._r8,    2._r8,   -3._r8,   -4._r8,    2._r8,    2._r8,    1._r8, &
       -5._r8,    2._r8,   -2._r8,    6._r8,    6._r8,   -4._r8,    4._r8,    0._r8, &
        0._r8,   -2._r8,    1._r8,    4._r8,    3._r8,    0._r8,    0._r8,   -6._r8  /)

  g1(:,7) = (/                                               & ! 1930
   -30805._r8,-2316._r8, 5808._r8, -951._r8, 2980._r8,-1424._r8, 1517._r8,  644._r8, &
     1172._r8,-1692._r8, -480._r8, 1205._r8,  133._r8,  907._r8,  166._r8,  896._r8, &
      727._r8,  205._r8,  584._r8, -195._r8, -422._r8, -109._r8,  234._r8,  -90._r8, &
     -237._r8,  327._r8,  -72._r8,  218._r8,   60._r8,  -32._r8,  -53._r8, -131._r8, &
     -118._r8,  -74._r8,   58._r8,   60._r8,   53._r8,    4._r8,   -9._r8,  102._r8, &
     -242._r8,   19._r8,  -32._r8,  -16._r8,   32._r8,    8._r8, -104._r8,  -46._r8, &
       74._r8,  -54._r8,  -51._r8,    4._r8,  -15._r8,   25._r8,  -14._r8,  -34._r8, &
        6._r8,  -12._r8,   29._r8,   18._r8,  -18._r8,    6._r8,  -20._r8,   11._r8, &
        7._r8,    8._r8,   -3._r8,  -15._r8,   -9._r8,    5._r8,    2._r8,  -14._r8, &
        5._r8,    5._r8,   -6._r8,   18._r8,    8._r8,   -5._r8,    8._r8,  -19._r8, &
        8._r8,   10._r8,  -20._r8,    1._r8,   14._r8,  -12._r8,    5._r8,   12._r8, &
       -3._r8,    1._r8,   -2._r8,   -2._r8,    9._r8,    3._r8,   10._r8,    0._r8, &
       -2._r8,   -2._r8,    2._r8,   -3._r8,   -4._r8,    2._r8,    2._r8,    1._r8, &
       -5._r8,    2._r8,   -2._r8,    6._r8,    6._r8,   -4._r8,    4._r8,    0._r8, &
        0._r8,   -2._r8,    1._r8,    4._r8,    3._r8,    0._r8,    0._r8,   -6._r8  /)

  g1(:,8) = (/                                               & ! 1935
   -30715._r8,-2306._r8, 5812._r8,-1018._r8, 2984._r8,-1520._r8, 1550._r8,  586._r8, &
     1206._r8,-1740._r8, -494._r8, 1215._r8,  146._r8,  918._r8,  101._r8,  903._r8, &
      744._r8,  188._r8,  565._r8, -226._r8, -415._r8,  -90._r8,  249._r8, -114._r8, &
     -241._r8,  329._r8,  -51._r8,  211._r8,   64._r8,  -33._r8,  -64._r8, -136._r8, &
     -115._r8,  -76._r8,   64._r8,   59._r8,   53._r8,    4._r8,   -8._r8,  104._r8, &
     -246._r8,   25._r8,  -25._r8,  -15._r8,   25._r8,    4._r8, -106._r8,  -40._r8, &
       74._r8,  -53._r8,  -52._r8,    4._r8,  -17._r8,   23._r8,  -14._r8,  -33._r8, &
        7._r8,  -11._r8,   29._r8,   18._r8,  -19._r8,    6._r8,  -19._r8,   11._r8, &
        7._r8,    8._r8,   -3._r8,  -15._r8,   -9._r8,    5._r8,    1._r8,  -15._r8, &
        6._r8,    5._r8,   -6._r8,   18._r8,    8._r8,   -5._r8,    7._r8,  -19._r8, &
        8._r8,   10._r8,  -20._r8,    1._r8,   15._r8,  -12._r8,    5._r8,   11._r8, &
       -3._r8,    1._r8,   -3._r8,   -2._r8,    9._r8,    3._r8,   11._r8,    0._r8, &
       -2._r8,   -2._r8,    2._r8,   -3._r8,   -4._r8,    2._r8,    2._r8,    1._r8, &
       -5._r8,    2._r8,   -2._r8,    6._r8,    6._r8,   -4._r8,    4._r8,    0._r8, &
        0._r8,   -1._r8,    2._r8,    4._r8,    3._r8,    0._r8,    0._r8,   -6._r8  /)

  g1(:,9) = (/                                               & ! 1940
   -30654._r8,-2292._r8, 5821._r8,-1106._r8, 2981._r8,-1614._r8, 1566._r8,  528._r8, &
     1240._r8,-1790._r8, -499._r8, 1232._r8,  163._r8,  916._r8,   43._r8,  914._r8, &
      762._r8,  169._r8,  550._r8, -252._r8, -405._r8,  -72._r8,  265._r8, -141._r8, &
     -241._r8,  334._r8,  -33._r8,  208._r8,   71._r8,  -33._r8,  -75._r8, -141._r8, &
     -113._r8,  -76._r8,   69._r8,   57._r8,   54._r8,    4._r8,   -7._r8,  105._r8, &
     -249._r8,   33._r8,  -18._r8,  -15._r8,   18._r8,    0._r8, -107._r8,  -33._r8, &
       74._r8,  -53._r8,  -52._r8,    4._r8,  -18._r8,   20._r8,  -14._r8,  -31._r8, &
        7._r8,   -9._r8,   29._r8,   17._r8,  -20._r8,    5._r8,  -19._r8,   11._r8, &
        7._r8,    8._r8,   -3._r8,  -14._r8,  -10._r8,    5._r8,    1._r8,  -15._r8, &
        6._r8,    5._r8,   -5._r8,   19._r8,    9._r8,   -5._r8,    7._r8,  -19._r8, &
        8._r8,   10._r8,  -21._r8,    1._r8,   15._r8,  -12._r8,    5._r8,   11._r8, &
       -3._r8,    1._r8,   -3._r8,   -2._r8,    9._r8,    3._r8,   11._r8,    1._r8, &
       -2._r8,   -2._r8,    2._r8,   -3._r8,   -4._r8,    2._r8,    2._r8,    1._r8, &
       -5._r8,    2._r8,   -2._r8,    6._r8,    6._r8,   -4._r8,    4._r8,    0._r8, &
        0._r8,   -1._r8,    2._r8,    4._r8,    3._r8,    0._r8,    0._r8,   -6._r8  /)

  g1(:,10) = (/                                              & ! 1945
   -30594._r8,-2285._r8, 5810._r8,-1244._r8, 2990._r8,-1702._r8, 1578._r8,  477._r8, &
     1282._r8,-1834._r8, -499._r8, 1255._r8,  186._r8,  913._r8,  -11._r8,  944._r8, &
      776._r8,  144._r8,  544._r8, -276._r8, -421._r8,  -55._r8,  304._r8, -178._r8, &
     -253._r8,  346._r8,  -12._r8,  194._r8,   95._r8,  -20._r8,  -67._r8, -142._r8, &
     -119._r8,  -82._r8,   82._r8,   59._r8,   57._r8,    6._r8,    6._r8,  100._r8, &
     -246._r8,   16._r8,  -25._r8,   -9._r8,   21._r8,  -16._r8, -104._r8,  -39._r8, &
       70._r8,  -40._r8,  -45._r8,    0._r8,  -18._r8,    0._r8,    2._r8,  -29._r8, &
        6._r8,  -10._r8,   28._r8,   15._r8,  -17._r8,   29._r8,  -22._r8,   13._r8, &
        7._r8,   12._r8,   -8._r8,  -21._r8,   -5._r8,  -12._r8,    9._r8,   -7._r8, &
        7._r8,    2._r8,  -10._r8,   18._r8,    7._r8,    3._r8,    2._r8,  -11._r8, &
        5._r8,  -21._r8,  -27._r8,    1._r8,   17._r8,  -11._r8,   29._r8,    3._r8, &
       -9._r8,   16._r8,    4._r8,   -3._r8,    9._r8,   -4._r8,    6._r8,   -3._r8, &
        1._r8,   -4._r8,    8._r8,   -3._r8,   11._r8,    5._r8,    1._r8,    1._r8, &
        2._r8,  -20._r8,   -5._r8,   -1._r8,   -1._r8,   -6._r8,    8._r8,    6._r8, &
       -1._r8,   -4._r8,   -3._r8,   -2._r8,    5._r8,    0._r8,   -2._r8,   -2._r8  /)

  g1(:,11) = (/                                              & ! 1950
   -30554._r8,-2250._r8, 5815._r8,-1341._r8, 2998._r8,-1810._r8, 1576._r8,  381._r8, &
     1297._r8,-1889._r8, -476._r8, 1274._r8,  206._r8,  896._r8,  -46._r8,  954._r8, &
      792._r8,  136._r8,  528._r8, -278._r8, -408._r8,  -37._r8,  303._r8, -210._r8, &
     -240._r8,  349._r8,    3._r8,  211._r8,  103._r8,  -20._r8,  -87._r8, -147._r8, &
     -122._r8,  -76._r8,   80._r8,   54._r8,   57._r8,   -1._r8,    4._r8,   99._r8, &
     -247._r8,   33._r8,  -16._r8,  -12._r8,   12._r8,  -12._r8, -105._r8,  -30._r8, &
       65._r8,  -55._r8,  -35._r8,    2._r8,  -17._r8,    1._r8,    0._r8,  -40._r8, &
       10._r8,   -7._r8,   36._r8,    5._r8,  -18._r8,   19._r8,  -16._r8,   22._r8, &
       15._r8,    5._r8,   -4._r8,  -22._r8,   -1._r8,    0._r8,   11._r8,  -21._r8, &
       15._r8,   -8._r8,  -13._r8,   17._r8,    5._r8,   -4._r8,   -1._r8,  -17._r8, &
        3._r8,   -7._r8,  -24._r8,   -1._r8,   19._r8,  -25._r8,   12._r8,   10._r8, &
        2._r8,    5._r8,    2._r8,   -5._r8,    8._r8,   -2._r8,    8._r8,    3._r8, &
      -11._r8,    8._r8,   -7._r8,   -8._r8,    4._r8,   13._r8,   -1._r8,   -2._r8, &
       13._r8,  -10._r8,   -4._r8,    2._r8,    4._r8,   -3._r8,   12._r8,    6._r8, &
        3._r8,   -3._r8,    2._r8,    6._r8,   10._r8,   11._r8,    3._r8,    8._r8  /)

  g1(:,12) = (/                                              & ! 1955
   -30500._r8,-2215._r8, 5820._r8,-1440._r8, 3003._r8,-1898._r8, 1581._r8,  291._r8, &
     1302._r8,-1944._r8, -462._r8, 1288._r8,  216._r8,  882._r8,  -83._r8,  958._r8, &
      796._r8,  133._r8,  510._r8, -274._r8, -397._r8,  -23._r8,  290._r8, -230._r8, &
     -229._r8,  360._r8,   15._r8,  230._r8,  110._r8,  -23._r8,  -98._r8, -152._r8, &
     -121._r8,  -69._r8,   78._r8,   47._r8,   57._r8,   -9._r8,    3._r8,   96._r8, &
     -247._r8,   48._r8,   -8._r8,  -16._r8,    7._r8,  -12._r8, -107._r8,  -24._r8, &
       65._r8,  -56._r8,  -50._r8,    2._r8,  -24._r8,   10._r8,   -4._r8,  -32._r8, &
        8._r8,  -11._r8,   28._r8,    9._r8,  -20._r8,   18._r8,  -18._r8,   11._r8, &
        9._r8,   10._r8,   -6._r8,  -15._r8,  -14._r8,    5._r8,    6._r8,  -23._r8, &
       10._r8,    3._r8,   -7._r8,   23._r8,    6._r8,   -4._r8,    9._r8,  -13._r8, &
        4._r8,    9._r8,  -11._r8,   -4._r8,   12._r8,   -5._r8,    7._r8,    2._r8, &
        6._r8,    4._r8,   -2._r8,    1._r8,   10._r8,    2._r8,    7._r8,    2._r8, &
       -6._r8,    5._r8,    5._r8,   -3._r8,   -5._r8,   -4._r8,   -1._r8,    0._r8, &
        2._r8,   -8._r8,   -3._r8,   -2._r8,    7._r8,   -4._r8,    4._r8,    1._r8, &
       -2._r8,   -3._r8,    6._r8,    7._r8,   -2._r8,   -1._r8,    0._r8,   -3._r8  /)

  g1(:,13) = (/                                              & ! 1960
   -30421._r8,-2169._r8, 5791._r8,-1555._r8, 3002._r8,-1967._r8, 1590._r8,  206._r8, &
     1302._r8,-1992._r8, -414._r8, 1289._r8,  224._r8,  878._r8, -130._r8,  957._r8, &
      800._r8,  135._r8,  504._r8, -278._r8, -394._r8,    3._r8,  269._r8, -255._r8, &
     -222._r8,  362._r8,   16._r8,  242._r8,  125._r8,  -26._r8, -117._r8, -156._r8, &
     -114._r8,  -63._r8,   81._r8,   46._r8,   58._r8,  -10._r8,    1._r8,   99._r8, &
     -237._r8,   60._r8,   -1._r8,  -20._r8,   -2._r8,  -11._r8, -113._r8,  -17._r8, &
       67._r8,  -56._r8,  -55._r8,    5._r8,  -28._r8,   15._r8,   -6._r8,  -32._r8, &
        7._r8,   -7._r8,   23._r8,   17._r8,  -18._r8,    8._r8,  -17._r8,   15._r8, &
        6._r8,   11._r8,   -4._r8,  -14._r8,  -11._r8,    7._r8,    2._r8,  -18._r8, &
       10._r8,    4._r8,   -5._r8,   23._r8,   10._r8,    1._r8,    8._r8,  -20._r8, &
        4._r8,    6._r8,  -18._r8,    0._r8,   12._r8,   -9._r8,    2._r8,    1._r8, &
        0._r8,    4._r8,   -3._r8,   -1._r8,    9._r8,   -2._r8,    8._r8,    3._r8, &
        0._r8,   -1._r8,    5._r8,    1._r8,   -3._r8,    4._r8,    4._r8,    1._r8, &
        0._r8,    0._r8,   -1._r8,    2._r8,    4._r8,   -5._r8,    6._r8,    1._r8, &
        1._r8,   -1._r8,   -1._r8,    6._r8,    2._r8,    0._r8,    0._r8,   -7._r8  /)

  g1(:,14) = (/                                              & ! 1965
   -30334._r8,-2119._r8, 5776._r8,-1662._r8, 2997._r8,-2016._r8, 1594._r8,  114._r8, &
     1297._r8,-2038._r8, -404._r8, 1292._r8,  240._r8,  856._r8, -165._r8,  957._r8, &
      804._r8,  148._r8,  479._r8, -269._r8, -390._r8,   13._r8,  252._r8, -269._r8, &
     -219._r8,  358._r8,   19._r8,  254._r8,  128._r8,  -31._r8, -126._r8, -157._r8, &
      -97._r8,  -62._r8,   81._r8,   45._r8,   61._r8,  -11._r8,    8._r8,  100._r8, &
     -228._r8,   68._r8,    4._r8,  -32._r8,    1._r8,   -8._r8, -111._r8,   -7._r8, &
       75._r8,  -57._r8,  -61._r8,    4._r8,  -27._r8,   13._r8,   -2._r8,  -26._r8, &
        6._r8,   -6._r8,   26._r8,   13._r8,  -23._r8,    1._r8,  -12._r8,   13._r8, &
        5._r8,    7._r8,   -4._r8,  -12._r8,  -14._r8,    9._r8,    0._r8,  -16._r8, &
        8._r8,    4._r8,   -1._r8,   24._r8,   11._r8,   -3._r8,    4._r8,  -17._r8, &
        8._r8,   10._r8,  -22._r8,    2._r8,   15._r8,  -13._r8,    7._r8,   10._r8, &
       -4._r8,   -1._r8,   -5._r8,   -1._r8,   10._r8,    5._r8,   10._r8,    1._r8, &
       -4._r8,   -2._r8,    1._r8,   -2._r8,   -3._r8,    2._r8,    2._r8,    1._r8, &
       -5._r8,    2._r8,   -2._r8,    6._r8,    4._r8,   -4._r8,    4._r8,    0._r8, &
        0._r8,   -2._r8,    2._r8,    3._r8,    2._r8,    0._r8,    0._r8,   -6._r8  /)

  g1(:,15) = (/                                              & ! 1970
   -30220._r8,-2068._r8, 5737._r8,-1781._r8, 3000._r8,-2047._r8, 1611._r8,   25._r8, &
     1287._r8,-2091._r8, -366._r8, 1278._r8,  251._r8,  838._r8, -196._r8,  952._r8, &
      800._r8,  167._r8,  461._r8, -266._r8, -395._r8,   26._r8,  234._r8, -279._r8, &
     -216._r8,  359._r8,   26._r8,  262._r8,  139._r8,  -42._r8, -139._r8, -160._r8, &
      -91._r8,  -56._r8,   83._r8,   43._r8,   64._r8,  -12._r8,   15._r8,  100._r8, &
     -212._r8,   72._r8,    2._r8,  -37._r8,    3._r8,   -6._r8, -112._r8,    1._r8, &
       72._r8,  -57._r8,  -70._r8,    1._r8,  -27._r8,   14._r8,   -4._r8,  -22._r8, &
        8._r8,   -2._r8,   23._r8,   13._r8,  -23._r8,   -2._r8,  -11._r8,   14._r8, &
        6._r8,    7._r8,   -2._r8,  -15._r8,  -13._r8,    6._r8,   -3._r8,  -17._r8, &
        5._r8,    6._r8,    0._r8,   21._r8,   11._r8,   -6._r8,    3._r8,  -16._r8, &
        8._r8,   10._r8,  -21._r8,    2._r8,   16._r8,  -12._r8,    6._r8,   10._r8, &
       -4._r8,   -1._r8,   -5._r8,    0._r8,   10._r8,    3._r8,   11._r8,    1._r8, &
       -2._r8,   -1._r8,    1._r8,   -3._r8,   -3._r8,    1._r8,    2._r8,    1._r8, &
       -5._r8,    3._r8,   -1._r8,    4._r8,    6._r8,   -4._r8,    4._r8,    0._r8, &
        1._r8,   -1._r8,    0._r8,    3._r8,    3._r8,    1._r8,   -1._r8,   -4._r8  /)

  g1(:,16) = (/                                              & ! 1975
   -30100._r8,-2013._r8, 5675._r8,-1902._r8, 3010._r8,-2067._r8, 1632._r8,  -68._r8, &
     1276._r8,-2144._r8, -333._r8, 1260._r8,  262._r8,  830._r8, -223._r8,  946._r8, &
      791._r8,  191._r8,  438._r8, -265._r8, -405._r8,   39._r8,  216._r8, -288._r8, &
     -218._r8,  356._r8,   31._r8,  264._r8,  148._r8,  -59._r8, -152._r8, -159._r8, &
      -83._r8,  -49._r8,   88._r8,   45._r8,   66._r8,  -13._r8,   28._r8,   99._r8, &
     -198._r8,   75._r8,    1._r8,  -41._r8,    6._r8,   -4._r8, -111._r8,   11._r8, &
       71._r8,  -56._r8,  -77._r8,    1._r8,  -26._r8,   16._r8,   -5._r8,  -14._r8, &
       10._r8,    0._r8,   22._r8,   12._r8,  -23._r8,   -5._r8,  -12._r8,   14._r8, &
        6._r8,    6._r8,   -1._r8,  -16._r8,  -12._r8,    4._r8,   -8._r8,  -19._r8, &
        4._r8,    6._r8,    0._r8,   18._r8,   10._r8,  -10._r8,    1._r8,  -17._r8, &
        7._r8,   10._r8,  -21._r8,    2._r8,   16._r8,  -12._r8,    7._r8,   10._r8, &
       -4._r8,   -1._r8,   -5._r8,   -1._r8,   10._r8,    4._r8,   11._r8,    1._r8, &
       -3._r8,   -2._r8,    1._r8,   -3._r8,   -3._r8,    1._r8,    2._r8,    1._r8, &
       -5._r8,    3._r8,   -2._r8,    4._r8,    5._r8,   -4._r8,    4._r8,   -1._r8, &
        1._r8,   -1._r8,    0._r8,    3._r8,    3._r8,    1._r8,   -1._r8,   -5._r8  /)

  g1(:,17) = (/                                              & ! 1980
   -29992._r8,-1956._r8, 5604._r8,-1997._r8, 3027._r8,-2129._r8, 1663._r8, -200._r8, &
     1281._r8,-2180._r8, -336._r8, 1251._r8,  271._r8,  833._r8, -252._r8,  938._r8, &
      782._r8,  212._r8,  398._r8, -257._r8, -419._r8,   53._r8,  199._r8, -297._r8, &
     -218._r8,  357._r8,   46._r8,  261._r8,  150._r8,  -74._r8, -151._r8, -162._r8, &
      -78._r8,  -48._r8,   92._r8,   48._r8,   66._r8,  -15._r8,   42._r8,   93._r8, &
     -192._r8,   71._r8,    4._r8,  -43._r8,   14._r8,   -2._r8, -108._r8,   17._r8, &
       72._r8,  -59._r8,  -82._r8,    2._r8,  -27._r8,   21._r8,   -5._r8,  -12._r8, &
       16._r8,    1._r8,   18._r8,   11._r8,  -23._r8,   -2._r8,  -10._r8,   18._r8, &
        6._r8,    7._r8,    0._r8,  -18._r8,  -11._r8,    4._r8,   -7._r8,  -22._r8, &
        4._r8,    9._r8,    3._r8,   16._r8,    6._r8,  -13._r8,   -1._r8,  -15._r8, &
        5._r8,   10._r8,  -21._r8,    1._r8,   16._r8,  -12._r8,    9._r8,    9._r8, &
       -5._r8,   -3._r8,   -6._r8,   -1._r8,    9._r8,    7._r8,   10._r8,    2._r8, &
       -6._r8,   -5._r8,    2._r8,   -4._r8,   -4._r8,    1._r8,    2._r8,    0._r8, &
       -5._r8,    3._r8,   -2._r8,    6._r8,    5._r8,   -4._r8,    3._r8,    0._r8, &
        1._r8,   -1._r8,    2._r8,    4._r8,    3._r8,    0._r8,    0._r8,   -6._r8  /)

  g1(:,18) = (/                                              & ! 1985
   -29873._r8,-1905._r8, 5500._r8,-2072._r8, 3044._r8,-2197._r8, 1687._r8, -306._r8, &
     1296._r8,-2208._r8, -310._r8, 1247._r8,  284._r8,  829._r8, -297._r8,  936._r8, &
      780._r8,  232._r8,  361._r8, -249._r8, -424._r8,   69._r8,  170._r8, -297._r8, &
     -214._r8,  355._r8,   47._r8,  253._r8,  150._r8,  -93._r8, -154._r8, -164._r8, &
      -75._r8,  -46._r8,   95._r8,   53._r8,   65._r8,  -16._r8,   51._r8,   88._r8, &
     -185._r8,   69._r8,    4._r8,  -48._r8,   16._r8,   -1._r8, -102._r8,   21._r8, &
       74._r8,  -62._r8,  -83._r8,    3._r8,  -27._r8,   24._r8,   -2._r8,   -6._r8, &
       20._r8,    4._r8,   17._r8,   10._r8,  -23._r8,    0._r8,   -7._r8,   21._r8, &
        6._r8,    8._r8,    0._r8,  -19._r8,  -11._r8,    5._r8,   -9._r8,  -23._r8, &
        4._r8,   11._r8,    4._r8,   14._r8,    4._r8,  -15._r8,   -4._r8,  -11._r8, &
        5._r8,   10._r8,  -21._r8,    1._r8,   15._r8,  -12._r8,    9._r8,    9._r8, &
       -6._r8,   -3._r8,   -6._r8,   -1._r8,    9._r8,    7._r8,    9._r8,    1._r8, &
       -7._r8,   -5._r8,    2._r8,   -4._r8,   -4._r8,    1._r8,    3._r8,    0._r8, &
       -5._r8,    3._r8,   -2._r8,    6._r8,    5._r8,   -4._r8,    3._r8,    0._r8, &
        1._r8,   -1._r8,    2._r8,    4._r8,    3._r8,    0._r8,    0._r8,   -6._r8  /)

  g1(:,19) = (/                                              & ! 1990
   -29775._r8,-1848._r8, 5406._r8,-2131._r8, 3059._r8,-2279._r8, 1686._r8, -373._r8, &
     1314._r8,-2239._r8, -284._r8, 1248._r8,  293._r8,  802._r8, -352._r8,  939._r8, &
      780._r8,  247._r8,  325._r8, -240._r8, -423._r8,   84._r8,  141._r8, -299._r8, &
     -214._r8,  353._r8,   46._r8,  245._r8,  154._r8, -109._r8, -153._r8, -165._r8, &
      -69._r8,  -36._r8,   97._r8,   61._r8,   65._r8,  -16._r8,   59._r8,   82._r8, &
     -178._r8,   69._r8,    3._r8,  -52._r8,   18._r8,    1._r8,  -96._r8,   24._r8, &
       77._r8,  -64._r8,  -80._r8,    2._r8,  -26._r8,   26._r8,    0._r8,   -1._r8, &
       21._r8,    5._r8,   17._r8,    9._r8,  -23._r8,    0._r8,   -4._r8,   23._r8, &
        5._r8,   10._r8,   -1._r8,  -19._r8,  -10._r8,    6._r8,  -12._r8,  -22._r8, &
        3._r8,   12._r8,    4._r8,   12._r8,    2._r8,  -16._r8,   -6._r8,  -10._r8, &
        4._r8,    9._r8,  -20._r8,    1._r8,   15._r8,  -12._r8,   11._r8,    9._r8, &
       -7._r8,   -4._r8,   -7._r8,   -2._r8,    9._r8,    7._r8,    8._r8,    1._r8, &
       -7._r8,   -6._r8,    2._r8,   -3._r8,   -4._r8,    2._r8,    2._r8,    1._r8, &
       -5._r8,    3._r8,   -2._r8,    6._r8,    4._r8,   -4._r8,    3._r8,    0._r8, &
        1._r8,   -2._r8,    3._r8,    3._r8,    3._r8,   -1._r8,    0._r8,   -6._r8  /)

  g2(:,1) = (/                                               & ! 1995
   -29692._r8,-1784._r8, 5306._r8,-2200._r8, 3070._r8,-2366._r8, 1681._r8, -413._r8, &
     1335._r8,-2267._r8, -262._r8, 1249._r8,  302._r8,  759._r8, -427._r8,  940._r8, &
      780._r8,  262._r8,  290._r8, -236._r8, -418._r8,   97._r8,  122._r8, -306._r8, &
     -214._r8,  352._r8,   46._r8,  235._r8,  165._r8, -118._r8, -143._r8, -166._r8, &
      -55._r8,  -17._r8,  107._r8,   68._r8,   67._r8,  -17._r8,   68._r8,   72._r8, &
     -170._r8,   67._r8,   -1._r8,  -58._r8,   19._r8,    1._r8,  -93._r8,   36._r8, &
       77._r8,  -72._r8,  -69._r8,    1._r8,  -25._r8,   28._r8,    4._r8,    5._r8, &
       24._r8,    4._r8,   17._r8,    8._r8,  -24._r8,   -2._r8,   -6._r8,   25._r8, &
        6._r8,   11._r8,   -6._r8,  -21._r8,   -9._r8,    8._r8,  -14._r8,  -23._r8, &
        9._r8,   15._r8,    6._r8,   11._r8,   -5._r8,  -16._r8,   -7._r8,   -4._r8, &
        4._r8,    9._r8,  -20._r8,    3._r8,   15._r8,  -10._r8,   12._r8,    8._r8, &
       -6._r8,   -8._r8,   -8._r8,   -1._r8,    8._r8,   10._r8,    5._r8,   -2._r8, &
       -8._r8,   -8._r8,    3._r8,   -3._r8,   -6._r8,    1._r8,    2._r8,    0._r8, &
       -4._r8,    4._r8,   -1._r8,    5._r8,    4._r8,   -5._r8,    2._r8,   -1._r8, &
        2._r8,   -2._r8,    5._r8,    1._r8,    1._r8,   -2._r8,    0._r8,   -7._r8, &
        0._r8,    0._r8,    0._r8,    0._r8,    0._r8,    0._r8,    0._r8,    0._r8, &
        0._r8,    0._r8,    0._r8,    0._r8,    0._r8,    0._r8,    0._r8,    0._r8, &
        0._r8,    0._r8,    0._r8,    0._r8,    0._r8,    0._r8,    0._r8,    0._r8, &
        0._r8,    0._r8,    0._r8,    0._r8,    0._r8,    0._r8,    0._r8,    0._r8, &
        0._r8,    0._r8,    0._r8,    0._r8,    0._r8,    0._r8,    0._r8,    0._r8, &
        0._r8,    0._r8,    0._r8,    0._r8,    0._r8,    0._r8,    0._r8,    0._r8, &
        0._r8,    0._r8,    0._r8,    0._r8,    0._r8,    0._r8,    0._r8,    0._r8, &
        0._r8,    0._r8,    0._r8,    0._r8,    0._r8,    0._r8,    0._r8,    0._r8, &
        0._r8,    0._r8,    0._r8,    0._r8,    0._r8,    0._r8,    0._r8,    0._r8, &
        0._r8,    0._r8,    0._r8 /)

  g2(:,2) = (/                                               & ! 2000
   -29619.4_r8,-1728.2_r8, 5186.1_r8,-2267.7_r8, 3068.4_r8,-2481.6_r8, 1670.9_r8, &
     -458.0_r8, 1339.6_r8,-2288.0_r8, -227.6_r8, 1252.1_r8,  293.4_r8,  714.5_r8, &
     -491.1_r8,  932.3_r8,  786.8_r8,  272.6_r8,  250.0_r8, -231.9_r8, -403.0_r8, &
      119.8_r8,  111.3_r8, -303.8_r8, -218.8_r8,  351.4_r8,   43.8_r8,  222.3_r8, &
      171.9_r8, -130.4_r8, -133.1_r8, -168.6_r8,  -39.3_r8,  -12.9_r8,  106.3_r8, &
       72.3_r8,   68.2_r8,  -17.4_r8,   74.2_r8,   63.7_r8, -160.9_r8,   65.1_r8, &
       -5.9_r8,  -61.2_r8,   16.9_r8,    0.7_r8,  -90.4_r8,   43.8_r8,   79.0_r8, &
      -74.0_r8,  -64.6_r8,    0.0_r8,  -24.2_r8,   33.3_r8,    6.2_r8,    9.1_r8, &
       24.0_r8,    6.9_r8,   14.8_r8,    7.3_r8,  -25.4_r8,   -1.2_r8,   -5.8_r8, &
       24.4_r8,    6.6_r8,   11.9_r8,   -9.2_r8,  -21.5_r8,   -7.9_r8,    8.5_r8, &
      -16.6_r8,  -21.5_r8,    9.1_r8,   15.5_r8,    7.0_r8,    8.9_r8,   -7.9_r8, &
      -14.9_r8,   -7.0_r8,   -2.1_r8,    5.0_r8,    9.4_r8,  -19.7_r8,    3.0_r8, &
       13.4_r8,   -8.4_r8,   12.5_r8,    6.3_r8,   -6.2_r8,   -8.9_r8,   -8.4_r8, &
       -1.5_r8,    8.4_r8,    9.3_r8,    3.8_r8,   -4.3_r8,   -8.2_r8,   -8.2_r8, &
        4.8_r8,   -2.6_r8,   -6.0_r8,    1.7_r8,    1.7_r8,    0.0_r8,   -3.1_r8, &
        4.0_r8,   -0.5_r8,    4.9_r8,    3.7_r8,   -5.9_r8,    1.0_r8,   -1.2_r8, &
        2.0_r8,   -2.9_r8,    4.2_r8,    0.2_r8,    0.3_r8,   -2.2_r8,   -1.1_r8, &
       -7.4_r8,    2.7_r8,   -1.7_r8,    0.1_r8,   -1.9_r8,    1.3_r8,    1.5_r8, &
       -0.9_r8,   -0.1_r8,   -2.6_r8,    0.1_r8,    0.9_r8,   -0.7_r8,   -0.7_r8, &
        0.7_r8,   -2.8_r8,    1.7_r8,   -0.9_r8,    0.1_r8,   -1.2_r8,    1.2_r8, &
       -1.9_r8,    4.0_r8,   -0.9_r8,   -2.2_r8,   -0.3_r8,   -0.4_r8,    0.2_r8, &
        0.3_r8,    0.9_r8,    2.5_r8,   -0.2_r8,   -2.6_r8,    0.9_r8,    0.7_r8, &
       -0.5_r8,    0.3_r8,    0.3_r8,    0.0_r8,   -0.3_r8,    0.0_r8,   -0.4_r8, &
        0.3_r8,   -0.1_r8,   -0.9_r8,   -0.2_r8,   -0.4_r8,   -0.4_r8,    0.8_r8, &
       -0.2_r8,   -0.9_r8,   -0.9_r8,    0.3_r8,    0.2_r8,    0.1_r8,    1.8_r8, &
       -0.4_r8,   -0.4_r8,    1.3_r8,   -1.0_r8,   -0.4_r8,   -0.1_r8,    0.7_r8, &
        0.7_r8,   -0.4_r8,    0.3_r8,    0.3_r8,    0.6_r8,   -0.1_r8,    0.3_r8, &
        0.4_r8,   -0.2_r8,    0.0_r8,   -0.5_r8,    0.1_r8,   -0.9_r8 /)

  g2(:,3) = (/                                               & ! 2005
   -29554.63_r8,-1669.05_r8, 5077.99_r8,-2337.24_r8, 3047.69_r8,-2594.50_r8,   &
    1657.76_r8, -515.43_r8, 1336.30_r8,-2305.83_r8, -198.86_r8, 1246.39_r8,    &
     269.72_r8,  672.51_r8, -524.72_r8,  920.55_r8,  797.96_r8,  282.07_r8,    &
     210.65_r8, -225.23_r8, -379.86_r8,  145.15_r8,  100.00_r8, -305.36_r8,    &
    -227.00_r8,  354.41_r8,   42.72_r8,  208.95_r8,  180.25_r8, -136.54_r8,    &
    -123.45_r8, -168.05_r8,  -19.57_r8,  -13.55_r8,  103.85_r8,   73.60_r8,    &
      69.56_r8,  -20.33_r8,   76.74_r8,   54.75_r8, -151.34_r8,   63.63_r8,    &
     -14.58_r8,  -63.53_r8,   14.58_r8,    0.24_r8,  -86.36_r8,   50.94_r8,    &
      79.88_r8,  -74.46_r8,  -61.14_r8,   -1.65_r8,  -22.57_r8,   38.73_r8,    &
       6.82_r8,   12.30_r8,   25.35_r8,    9.37_r8,   10.93_r8,    5.42_r8,    &
     -26.32_r8,    1.94_r8,   -4.64_r8,   24.80_r8,    7.62_r8,   11.20_r8,    &
     -11.73_r8,  -20.88_r8,   -6.88_r8,    9.83_r8,  -18.11_r8,  -19.71_r8,    &
      10.17_r8,   16.22_r8,    9.36_r8,    7.61_r8,  -11.25_r8,  -12.76_r8,    &
      -4.87_r8,   -0.06_r8,    5.58_r8,    9.76_r8,  -20.11_r8,    3.58_r8,    &
      12.69_r8,   -6.94_r8,   12.67_r8,    5.01_r8,   -6.72_r8,  -10.76_r8,    &
      -8.16_r8,   -1.25_r8,    8.10_r8,    8.76_r8,    2.92_r8,   -6.66_r8,    &
      -7.73_r8,   -9.22_r8,    6.01_r8,   -2.17_r8,   -6.12_r8,    2.19_r8,    &
       1.42_r8,    0.10_r8,   -2.35_r8,    4.46_r8,   -0.15_r8,    4.76_r8,    &
       3.06_r8,   -6.58_r8,    0.29_r8,   -1.01_r8,    2.06_r8,   -3.47_r8,    &
       3.77_r8,   -0.86_r8,   -0.21_r8,   -2.31_r8,   -2.09_r8,   -7.93_r8,    &
       2.95_r8,   -1.60_r8,    0.26_r8,   -1.88_r8,    1.44_r8,    1.44_r8,    &
      -0.77_r8,   -0.31_r8,   -2.27_r8,    0.29_r8,    0.90_r8,   -0.79_r8,    &
      -0.58_r8,    0.53_r8,   -2.69_r8,    1.80_r8,   -1.08_r8,    0.16_r8,    &
      -1.58_r8,    0.96_r8,   -1.90_r8,    3.99_r8,   -1.39_r8,   -2.15_r8,    &
      -0.29_r8,   -0.55_r8,    0.21_r8,    0.23_r8,    0.89_r8,    2.38_r8,    &
      -0.38_r8,   -2.63_r8,    0.96_r8,    0.61_r8,   -0.30_r8,    0.40_r8,    &
       0.46_r8,    0.01_r8,   -0.35_r8,    0.02_r8,   -0.36_r8,    0.28_r8,    &
       0.08_r8,   -0.87_r8,   -0.49_r8,   -0.34_r8,   -0.08_r8,    0.88_r8,    &
      -0.16_r8,   -0.88_r8,   -0.76_r8,    0.30_r8,    0.33_r8,    0.28_r8,    &
       1.72_r8,   -0.43_r8,   -0.54_r8,    1.18_r8,   -1.07_r8,   -0.37_r8,    &
      -0.04_r8,    0.75_r8,    0.63_r8,   -0.26_r8,    0.21_r8,    0.35_r8,    &
       0.53_r8,   -0.05_r8,    0.38_r8,    0.41_r8,   -0.22_r8,   -0.10_r8,    &
      -0.57_r8,   -0.18_r8,   -0.82_r8 /)

  g2(:,4) = (/                                                & ! 2010
   -29496.5_r8,-1585.9_r8, 4945.1_r8,-2396.6_r8, 3026.0_r8,-2707.7_r8, 1668.6_r8,  &
     -575.4_r8, 1339.7_r8,-2326.3_r8, -160.5_r8, 1231.7_r8,  251.7_r8,  634.2_r8,  &
     -536.8_r8,  912.6_r8,  809.0_r8,  286.4_r8,  166.6_r8, -211.2_r8, -357.1_r8,  &
      164.4_r8,   89.7_r8, -309.2_r8, -231.1_r8,  357.2_r8,   44.7_r8,  200.3_r8,  &
      188.9_r8, -141.2_r8, -118.1_r8, -163.1_r8,    0.1_r8,   -7.7_r8,  100.9_r8,  &
       72.8_r8,   68.6_r8,  -20.8_r8,   76.0_r8,   44.2_r8, -141.4_r8,   61.5_r8,  &
      -22.9_r8,  -66.3_r8,   13.1_r8,    3.1_r8,  -77.9_r8,   54.9_r8,   80.4_r8,  &
      -75.0_r8,  -57.8_r8,   -4.7_r8,  -21.2_r8,   45.3_r8,    6.6_r8,   14.0_r8,  &
       24.9_r8,   10.4_r8,    7.0_r8,    1.6_r8,  -27.7_r8,    4.9_r8,   -3.4_r8,  &
       24.3_r8,    8.2_r8,   10.9_r8,  -14.5_r8,  -20.0_r8,   -5.7_r8,   11.9_r8,  &
      -19.3_r8,  -17.4_r8,   11.6_r8,   16.7_r8,   10.9_r8,    7.1_r8,  -14.1_r8,  &
      -10.8_r8,   -3.7_r8,    1.7_r8,    5.4_r8,    9.4_r8,  -20.5_r8,    3.4_r8,  &
       11.6_r8,   -5.3_r8,   12.8_r8,    3.1_r8,   -7.2_r8,  -12.4_r8,   -7.4_r8,  &
       -0.8_r8,    8.0_r8,    8.4_r8,    2.2_r8,   -8.4_r8,   -6.1_r8,  -10.1_r8,  &
        7.0_r8,   -2.0_r8,   -6.3_r8,    2.8_r8,    0.9_r8,   -0.1_r8,   -1.1_r8,  &
        4.7_r8,   -0.2_r8,    4.4_r8,    2.5_r8,   -7.2_r8,   -0.3_r8,   -1.0_r8,  &
        2.2_r8,   -4.0_r8,    3.1_r8,   -2.0_r8,   -1.0_r8,   -2.0_r8,   -2.8_r8,  &
       -8.3_r8,    3.0_r8,   -1.5_r8,    0.1_r8,   -2.1_r8,    1.7_r8,    1.6_r8,  &
       -0.6_r8,   -0.5_r8,   -1.8_r8,    0.5_r8,    0.9_r8,   -0.8_r8,   -0.4_r8,  &
        0.4_r8,   -2.5_r8,    1.8_r8,   -1.3_r8,    0.2_r8,   -2.1_r8,    0.8_r8,  &
       -1.9_r8,    3.8_r8,   -1.8_r8,   -2.1_r8,   -0.2_r8,   -0.8_r8,    0.3_r8,  &
        0.3_r8,    1.0_r8,    2.2_r8,   -0.7_r8,   -2.5_r8,    0.9_r8,    0.5_r8,  &
       -0.1_r8,    0.6_r8,    0.5_r8,    0.0_r8,   -0.4_r8,    0.1_r8,   -0.4_r8,  &
        0.3_r8,    0.2_r8,   -0.9_r8,   -0.8_r8,   -0.2_r8,    0.0_r8,    0.8_r8,  &
       -0.2_r8,   -0.9_r8,   -0.8_r8,    0.3_r8,    0.3_r8,    0.4_r8,    1.7_r8,  &
       -0.4_r8,   -0.6_r8,    1.1_r8,   -1.2_r8,   -0.3_r8,   -0.1_r8,    0.8_r8,  &
        0.5_r8,   -0.2_r8,    0.1_r8,    0.4_r8,    0.5_r8,    0.0_r8,    0.4_r8,  &
        0.4_r8,   -0.2_r8,   -0.3_r8,   -0.5_r8,   -0.3_r8,   -0.8_r8 /)

  g2(1:80,5) = (/                                            & ! 2012
      11.4_r8,   16.7_r8,  -28.8_r8,  -11.3_r8,   -3.9_r8,  -23.0_r8,    2.7_r8,  &
     -12.9_r8,    1.3_r8,   -3.9_r8,    8.6_r8,   -2.9_r8,   -2.9_r8,   -8.1_r8,  &
      -2.1_r8,   -1.4_r8,    2.0_r8,    0.4_r8,   -8.9_r8,    3.2_r8,    4.4_r8,  &
       3.6_r8,   -2.3_r8,   -0.8_r8,   -0.5_r8,    0.5_r8,    0.5_r8,   -1.5_r8,  &
       1.5_r8,   -0.7_r8,    0.9_r8,    1.3_r8,    3.7_r8,    1.4_r8,   -0.6_r8,  &
      -0.3_r8,   -0.3_r8,   -0.1_r8,   -0.3_r8,   -2.1_r8,    1.9_r8,   -0.4_r8,  &
      -1.6_r8,   -0.5_r8,   -0.2_r8,    0.8_r8,    1.8_r8,    0.5_r8,    0.2_r8,  &
      -0.1_r8,    0.6_r8,   -0.6_r8,    0.3_r8,    1.4_r8,   -0.2_r8,    0.3_r8,  &
      -0.1_r8,    0.1_r8,   -0.8_r8,   -0.8_r8,   -0.3_r8,    0.4_r8,    0.2_r8,  &
      -0.1_r8,    0.1_r8,    0.0_r8,   -0.5_r8,    0.2_r8,    0.3_r8,    0.5_r8,  &
      -0.3_r8,    0.4_r8,    0.3_r8,    0.1_r8,    0.2_r8,   -0.1_r8,   -0.5_r8,  &
       0.4_r8,    0.2_r8,    0.4_r8 /)

  g2(81:n2,5) = 0._r8
!
! Set gh from g1,g2:
!
  do n=1,ncn1
    i = (n-1)*n1
    gh(i+1:i+n1) = g1(:,n)
!   write(iulog,"('cofrm: n=',i3,' i+1:i+n1=',i4,':',i4)") n,i+1,i+n1
  enddo
  do n=1,ncn2
    i = n1*ncn1 + (n-1)*n2
    gh(i+1:i+n2) = g2(:,n)
!   write(iulog,"('cofrm: n=',i3,' i+1:i+n2=',i4,':',i4)") n,i+1,i+n2
  enddo
  gh(ngh) = 0._r8 ! not sure why gh is dimensioned with the extra element, so set it to 0. 
  
  if (date < 2010._r8) then
    t   = 0.2_r8*(date - 1900._r8)
    ll  = t
    one = ll
    t   = t - one
    if (date < 1995._r8) then
      nmx   = 10
      nc    = nmx*(nmx+2)
      ll    = nc*ll
      kmx   = (nmx+1)*(nmx+2)/2
    else
      nmx   = 13
      nc    = nmx*(nmx+2)
      ll    = 0.2_r8*(date - 1995.0_r8)
      ll    = 120*19 + nc*ll
      kmx   = (nmx+1)*(nmx+2)/2
    endif
    tc    = 1.0_r8 - t
    if (isv.eq.1) then
      tc = -0.2_r8
      t = 0.2_r8
    endif
  else ! date >= 2010
    t     = date - 2010.0_r8
    tc    = 1.0_r8
    if (isv.eq.1) then
      t = 1.0_r8
      tc = 0.0_r8
    end if
    ll    = 2865
    nmx   = 13
    nc    = nmx*(nmx+2)
    kmx   = (nmx+1)*(nmx+2)/2
  endif ! date < 2010
  r = alt
  l = 1
  m = 1
  n = 0
!
! Set outputs gb(ncoef) and gv(ncoef)
! These are module data above.
! 
  gb(1) = 0._r8
  gv(1) = 0._r8
  f0 = -1.e-5_r8
  do k=2,kmx
    if (n < m) then
      m = 0
      n = n+1
    endif ! n < m
    lm = ll + l
    if (m == 0) f0 = f0 * dble(n)/2._r8
    if (m == 0) f  = f0 / sqrt(2.0_r8)
    nn = n+1
    mm = 1      
      
    if (m /= 0) then
      f = f / sqrt(dble(n-m+1) / dble(n+m) )
      gb(l+1)  = (tc*gh(lm) + t*gh(lm+nc))* f
    else   
      gb(l+1)  = (tc*gh(lm) + t*gh(lm+nc))* f0
    endif  
    gv(l+1) = gb(l+1)/dble(nn)
    if (m /= 0) then
      gb(l+2)  = (tc*gh(lm+1) + t*gh(lm+nc+1))*f
      gv(l+2) = gb(l+2)/dble(nn)
      l = l+2
    else
      l = l+1
    endif
    m = m+1
  enddo

! write(iulog,"('cofrm: ncoef=',i4,' gb=',/,(6f12.3))") ncoef,gb
! write(iulog,"('cofrm: ncoef=',i4,' gv=',/,(6f12.3))") ncoef,gv

end subroutine cofrm
!-----------------------------------------------------------------------
subroutine apex_subsol(iyr,iday,ihr,imn,sec,sbsllat,sbsllon)
!
! Find subsolar geographic latitude and longitude given the
! date and time (Universal Time).
!     
! This is based on formulas in Astronomical Almanac for the
! year 1996, p.  C24. (U.S.  Government Printing Office,
! 1994).  According to the Almanac, results are good to at
! least 0.01 degree latitude and 0.025 degree longitude
! between years 1950 and 2050.  Accuracy for other years has
! not been tested although the algorithm has been designed to
! accept input dates from 1601 to 2100.  Every day is assumed
! to have exactly 86400 seconds; thus leap seconds that
! sometimes occur on June 30 and December 31 are ignored:
! their effect is below the accuracy threshold of the algorithm.
!     
! 961026 A. D. Richmond, NCAR
!
! Input Args:
  integer,intent(in) :: &
    iyr,   & ! Year (e.g., 1994). IYR must be in the range: 1601 to 2100.
    iday,  & ! Day number of year (e.g., IDAY = 32 for Feb 1)
    ihr,   & ! Hour of day    (e.g., 13 for 13:49)
    imn      ! Minute of hour (e.g., 49 for 13:49)
  real(r8),intent(in) :: sec ! Second and fraction after the hour/minute.
!
! Output Args:
  real(r8),intent(out) :: &
    sbsllat, & ! geographic latitude of subsolar point (degrees)
    sbsllon    ! geographic longitude of subsolar point (-180 to +180)
!
! Local:
  integer,parameter :: minyear=1601, maxyear = 2100
  real(r8),parameter :: & ! Use local params for compatability w/ legacy code,
                          ! but probably would be ok to use module data dtr,rtd
    d2r=0.0174532925199432957692369076847_r8, &
    r2d=57.2957795130823208767981548147_r8
  real(r8) :: yr,l0,g0,ut,df,lf,gf,l,g,grad,n,epsilon,epsrad,alpha,delta,&
    etdeg,aptime,lambda,lamrad,sinlam
  integer :: nleap,ncent,nrot

  sbsllat=0._r8 ; sbsllon=0._r8

  yr = iyr-2000
!
! nleap (final) = number of leap days from (2000 January 1) to (IYR January 1)
!                 (negative if iyr is before 1997)
  nleap = (iyr-1601)/4
  nleap = nleap - 99
  if (iyr <= 1900) then
    if (iyr < minyear) then
      write(iulog,*) 'subsolr invalid before ',minyear,': input year = ',iyr
      call endrun( 'subsolr' )
    endif
    ncent = (iyr-minyear)/100
    ncent = 3 - ncent
    nleap = nleap + ncent
  endif
  if (iyr > maxyear) then
    write(iulog,*) 'subsolr invalid after ',maxyear,':  input year = ',iyr
    call endrun( 'subsolr' )
  endif
!
! L0 = Mean longitude of Sun at 12 UT on January 1 of IYR:
!     L0 = 280.461 + .9856474*(365*(YR-NLEAP) + 366*NLEAP)
!          - (ARBITRARY INTEGER)*360.
!        = 280.461 + .9856474*(365*(YR-4*NLEAP) + (366+365*3)*NLEAP)
!          - (ARBITRARY INTEGER)*360.
!        = (280.461 - 360.) + (.9856474*365 - 360.)*(YR-4*NLEAP)
!          + (.9856474*(366+365*3) - 4*360.)*NLEAP,
!  where ARBITRARY INTEGER = YR+1.  This gives:
!
      l0 = -79.549_r8 + (-.238699_r8*(yr-4*nleap) + 3.08514e-2_r8*nleap)
!                 
! G0 = Mean anomaly at 12 UT on January 1 of IYR:
!     G0 = 357.528 + .9856003*(365*(YR-NLEAP) + 366*NLEAP)
!          - (ARBITRARY INTEGER)*360.
!        = 357.528 + .9856003*(365*(YR-4*NLEAP) + (366+365*3)*NLEAP) 
!          - (ARBITRARY INTEGER)*360.
!        = (357.528 - 360.) + (.9856003*365 - 360.)*(YR-4*NLEAP)
!          + (.9856003*(366+365*3) - 4*360.)*NLEAP,
!  where ARBITRARY INTEGER = YR+1.  This gives:
!
      g0 = -2.472_r8 + (-.2558905_r8*(yr-4*nleap) - 3.79617e-2_r8*nleap)
!     
! Universal time in seconds:
      ut = dble(ihr*3600 + imn*60) + sec
!
! Days (including fraction) since 12 UT on January 1 of IYR:
      df = (ut/86400._r8 - 1.5_r8) + iday
!
! Addition to Mean longitude of Sun since January 1 of IYR: 
      lf = .9856474_r8*df
! 
! Addition to Mean anomaly since January 1 of IYR:
      gf = .9856003_r8*df
! 
! Mean longitude of Sun:
      l = l0 + lf
! 
! Mean anomaly:
      g = g0 + gf
      grad = g*d2r
! 
! Ecliptic longitude:
      lambda = l + 1.915_r8*sin(grad) + .020_r8*sin(2._r8*grad)
      lamrad = lambda*d2r
      sinlam = sin(lamrad)
! 
! Days (including fraction) since 12 UT on January 1 of 2000:
      n = df + 365._r8*yr + dble(nleap)
! 
! Obliquity of ecliptic: 
      epsilon = 23.439_r8 - 4.e-7_r8*n
      epsrad = epsilon*d2r
! 
! Right ascension:
      alpha = atan2(cos(epsrad)*sinlam,cos(lamrad))*r2d
! 
! Declination:
      delta = asin(sin(epsrad)*sinlam)*r2d
! 
! Subsolar latitude (output argument):
      sbsllat = delta
! 
! Equation of time (degrees):
      etdeg = l - alpha
      nrot = nint(etdeg/360._r8)
      etdeg = etdeg - dble(360*nrot)
! 
! Apparent time (degrees):
! Earth rotates one degree every 240 s.
      aptime = ut/240._r8 + etdeg
!
! Subsolar longitude (output argument):
      sbsllon = 180._r8 - aptime
      nrot = nint(sbsllon/360._r8)
      sbsllon = sbsllon - dble(360*nrot)

end subroutine apex_subsol

!-----------------------------------------------------------------------
subroutine solgmlon(xlat,xlon,colat,elon,mlon)
!
! Compute geomagnetic longitude of the point with geocentric spherical
!  latitude and longitude of XLAT and XLON, respectively.
! 940719 A. D. Richmond, NCAR
!
! Algorithm:
!   Use spherical coordinates. 
!   Let GP be geographic pole. 
!   Let GM be geomagnetic pole (colatitude COLAT, east longitude ELON).
!   Let XLON be longitude of point P.
!   Let TE be colatitude of point P.
!   Let ANG be longitude angle from GM to P.
!   Let TP be colatitude of GM.
!   Let TF be arc length between GM and P.
!   Let PA = MLON be geomagnetic longitude, i.e., Pi minus angle measured
!     counterclockwise from arc GM-P to arc GM-GP. 
!   Then, using notation C=cos, S=sin, spherical-trigonometry formulas
!     for the functions of the angles are as shown below.  Note: STFCPA,
!     STFSPA are sin(TF) times cos(PA), sin(PA), respectively.
!
! Input Args:
  real(r8),intent(in)  :: xlat,xlon,colat,elon
! 
! Output Arg: Geomagnetic dipole longitude of the point (deg, -180. to 180.)
  real(r8),intent(out) :: mlon 
!
! Local:
  real(r8),parameter ::           &
    rtod=5.72957795130823e1_r8,  &
    dtor=1.745329251994330e-2_r8
  real(r8) :: ctp,stp,ang,cang,sang,cte,ste,stfcpa,stfspa

  ctp = cos(colat*dtor)
  stp = sqrt(1._r8 - ctp*ctp)
  ang = (xlon-elon)*dtor
  cang = cos(ang)
  sang = sin(ang)
  cte = sin(xlat*dtor)
  ste = sqrt(1._r8-cte*cte)
  stfcpa = ste*ctp*cang - cte*stp
  stfspa = sang*ste
  mlon = atan2(stfspa,stfcpa)*rtod

end subroutine solgmlon
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
subroutine apex_magloctm (alon,sbsllat,sbsllon,clatp,polon,mlt)
  !
  !-----------------------------------------------------------------------
  !  Computes magnetic local time from magnetic longitude, subsolar coordinates,
  !   and geomagnetic pole coordinates.
  !  950302 A. D. Richmond, NCAR
  !  Algorithm:  MLT is calculated from the difference of the apex longitude,
  !   alon, and the geomagnetic dipole longitude of the subsolar point.
  !
  !   Inputs:
  !    alon    = apex magnetic longitude of the point (deg)
  !    sbsllat = geographic latitude of subsolar point (degrees)
  !    sbsllon = geographic longitude of subsolar point (degrees)
  !    clatp   = Geocentric colatitude of geomagnetic dipole north pole (deg)
  !    polon   = East longitude of geomagnetic dipole north pole (deg)
  !
  !   Output:
  !    mlt (real) = magnetic local time for the apex longitude alon (hours)
  !
  !-----------------------------------------------------------------------
  !
  !------------------------------Arguments--------------------------------
  !
  REAL(r8) alon, sbsllat, sbsllon, clatp, polon, MLT
  !
  !---------------------------Local variables-----------------------------
  !
  real(r8) smlon
  !
  !-----------------------------------------------------------------------
  !
  call solgmlon (sbsllat,sbsllon,clatp,polon,smlon)
  mlt = (alon - smlon)/15.0_r8 + 12.0_r8
  if (mlt .ge. 24.0_r8) mlt = mlt - 24.0_r8
  if (mlt .lt.   0._r8) mlt = mlt + 24.0_r8
  return
end subroutine apex_magloctm

!-----------------------------------------------------------------------
end module apex
