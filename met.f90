!
!////////////////////////////////////////////////////////////////////////
!
!      met.f90
!      Created: 12 August 2016 14:39 
!      By: David Gillies  
!
!////////////////////////////////////////////////////////////////////////
!




MODULE DO3SE_Met




real, parameter :: DEG2RAD = 0.017453292519943295
real, parameter :: PAR_Wm2_to_photons = 4.57
real, parameter :: seaP = 101.325
real, parameter :: PARfrac = 0.45
  REAL, parameter :: UNDEF = -999.0
  INTEGER, parameter :: IUNDEF = -999
  
  interface is_undef
    module procedure is_undef_real
    module procedure is_undef_integer
  end interface is_undef
  public :: is_undef

  interface is_def
    module procedure is_def_real
    module procedure is_def_integer
  end interface is_def
  public :: is_def

CONTAINS

  subroutine calc_met_data(Ts_C, P, dd, hr, lat, lon, elev, albedo, &
                          CO2_method, CO2_constant, CO2, & ! ins
                          Rn, R, PAR, PPFD, VPD, RH, Idrctt, Idfuse, & ! inouts
                          esat, eact, sinB, VPD_dd) ! outs
   REAL, INTENT(in) :: Ts_C
   REAL, INTENT(inout) :: VPD
   REAL, INTENT(inout) :: RH
   REAL, INTENT(out) :: esat
   REAL, INTENT(out) :: eact
    REAL, INTENT(in) :: P        !< Atmospheric pressure (kPa)
    integer, intent(in) :: dd   !< Day of year (1--365)
    integer, intent(in) :: hr   !< Hour of day (0--23)
    REAL, INTENT(in) :: lat
    REAL, INTENT(in) :: lon
    REAL, INTENT(inout) :: elev
    REAL, INTENT(inout) :: albedo
    REAL, INTENT(inout) :: Rn
    REAL, INTENT(inout) :: R
    REAL, INTENT(inout) :: PAR 
    REAL, INTENT(inout) :: PPFD
    REAL, INTENT(inout) :: Idrctt
    REAL, INTENT(inout) :: Idfuse
    REAL, INTENT(out) :: sinB
   CHARACTER(len=16), INTENT(in) :: CO2_method
   REAL, INTENT(in) :: CO2_constant
   REAL, INTENT(inout) :: CO2
   REAL, INTENT(out) :: VPD_dd
   LOGICAL :: t_or_f
    ! Fixup met data

    call met_humidity(Ts_C, VPD, RH, esat, eact)
    call met_radiation(P, dd, hr, lat, lon, elev, albedo, Rn, R, PAR, PPFD, &
                                 Idrctt, Idfuse, sinB)

    call met_CO2(CO2_method, CO2_constant, CO2)

    ! Reset VPD sum at start of day
    
    t_or_f = start_of_day(hr)
    if (t_or_f) then
      VPD_dd = 0
    end if
    ! Only accumulate VPD sum during daylight hours
    t_or_f = is_daylight(R)
    if (t_or_f) then
      VPD_dd = VPD_dd + VPD
    end if
  end subroutine calc_met_data
  
  
  pure subroutine met_humidity(Ts_C, VPD, RH, esat, eact)
   REAL, INTENT(in) :: Ts_C
   REAL, INTENT(inout) :: VPD
   REAL, INTENT(inout) :: RH
   REAL, INTENT(out) :: esat
   REAL, INTENT(out) :: eact 

    esat = saturated_vapour_pressure(Ts_C)
    if (is_undef(VPD) .and. is_undef(RH)) then
      ! TODO: error if neither VPD or RH supplied
    else if (is_undef(RH)) then
      ! Calculate relative humidity from VPD
      eact = esat - VPD
      RH = eact / esat
    else if (is_undef(VPD)) then
      ! Calculate VPD from relative humidity
      eact = esat * RH
      VPD = esat - eact
    end if
  end subroutine met_humidity
  
  pure subroutine met_radiation(P, dd, hr, lat, lon, elev, albedo, Rn, R, PAR, PPFD, &
                                 Idrctt, Idfuse, sinB)
    REAL, INTENT(in) :: P        !< Atmospheric pressure (kPa)
    integer, intent(in) :: dd   !< Day of year (1--365)
    integer, intent(in) :: hr   !< Hour of day (0--23)
    REAL, INTENT(in) :: lat
    REAL, INTENT(in) :: lon
    REAL, INTENT(inout) :: elev
    REAL, INTENT(inout) :: albedo
    REAL, INTENT(inout) :: Rn
    REAL, INTENT(inout) :: R
    REAL, INTENT(inout) :: PAR 
    REAL, INTENT(inout) :: PPFD
    REAL, INTENT(inout) :: Idrctt
    REAL, INTENT(inout) :: Idfuse
    REAL, INTENT(out) :: sinB
   

    ! Solar elevation
    sinB = solar_elevation(lat, lon, dd, hr)
    ! Irradiance
    if (is_undef(PAR)) then
      ! Calculate PAR from some other available source
      if (is_def(Idrctt) .and. is_def(Idfuse)) then
        ! PAR is sum of direct and diffuse PAR
        PAR = Idrctt + Idfuse
      else if (is_def(PPFD)) then
        ! PAR from PPFD
        PAR = PPFD / PAR_Wm2_to_photons
      else if (is_def(R)) then
        ! Estimate PAR from global radiation
        PAR = R * PARfrac
      else
        ! TODO: error if no source of PAR is available
      end if
    end if
    if (is_undef(PPFD)) then
      ! If PPFD is missing, convert from PAR
      PPFD = PAR * PAR_Wm2_to_photons
    end if
    if (is_undef(Idrctt) .or. is_undef(Idfuse)) then
      ! If direct + diffuse are missing, estimate from PAR
      call PAR_direct_diffuse(PAR, sinB, P, Idrctt, Idfuse)
    end if
    if (is_undef(R)) then
      ! If global radiation is absent, estimate from PAR
      R = PAR / PARfrac
    end if

    ! Net radiation
    if (is_undef(Rn)) THEN
      Rn = net_radiation(lat, lon, elev, albedo, &
                            dd, hr, sinB, R, Ts_C, eact)
          
       
                     
    end if
  end subroutine met_radiation
  
  
  pure subroutine PAR_direct_diffuse(PAR, sinB, P, Idrctt, Idfuse)
    real, intent(in) :: PAR       !< Photosynthetically active radiation (W m-2)
    real, intent(in) :: sinB      !< sin() of solar elevation angle
    real, intent(in) :: P         !< Atmospheric pressure (kPa)

    real, intent(out) :: Idrctt   !< Direct PAR irradiance (W m-2)
    real, intent(out) :: Idfuse   !< Diffuse PAR irradiance (W m-2)

    real :: m, pPARdir, pPARdif, pPARtotal, ST, fPARdir, fPARdif

    if (sinB > 0.0) then
      m = 1.0 / sinB

      ! Potential direct PAR
      pPARdir = 600 * exp(-0.185 * (P / seaP) * m) * sinB
      ! Potential diffuse PAR
      pPARdif = 0.4 * (600 - pPARdir) * sinB
      ! Potential total PAR
      pPARtotal = pPARdir + pPARdif

      ! Sky transmissivity
      ST = max(0.21, min(0.9, PAR / pPARtotal))

      ! Direct and diffuse fractions
      fPARdir = (pPARdir / pPARtotal) * (1.0 - ((0.9 - ST) / 0.7)**(2.0/3.0))
      fPARdif = 1 - fPARdir

      ! Apply calculated direct and diffuse fractions to PARtotal
      Idrctt = fPARdir * PAR
      Idfuse = fPARdif * PAR
    else
      Idrctt = 0.0
      Idfuse = 0.0
    end if
  end subroutine PAR_direct_diffuse
  
  
  
  pure function solar_elevation(lat, lon, dd, hr) result(sinB)
    real, intent(in) :: lat     !< Latitude (degrees North)
    real, intent(in) :: lon     !< Longitude (degrees East)
    integer, intent(in) :: dd   !< Day of year (1--365)
    integer, intent(in) :: hr   !< Hour of day (0--23)

    real :: sinB      ! Output: sin() of solar elevation angle

    real :: t0_, dec, h

    t0_ = solar_noon(lon, dd)
    dec = solar_declination(dd)

    ! Hour-angle of the sun
    h = (15 * (hr - t0_)) * DEG2RAD

    ! sin() of solar elevation angle
    sinB = sin(lat * DEG2RAD)*sin(dec) + cos(lat * DEG2RAD)*cos(dec)*cos(h)
    ! TODO: should this line be removed? what effect does it have?  Does any 
    !       use of sinB happen when sinB < 0?
    sinB = max(0.0, sinB)
  end function solar_elevation

  !> Calculate solar declination (radians).
  pure real function solar_declination(dd)
    integer, intent(in) :: dd   !< Day of year

    solar_declination = (-23.4 * cos((360 * ((dd + 10) / 365.0))*DEG2RAD))*DEG2RAD
  end function solar_declination

  pure real function solar_noon(lon, dd)
    real, intent(in) :: lon     !< Longitude (degrees East)
    integer, intent(in) :: dd   !< Day of year

    real :: f, e, lonm, LC

    ! Solar noon correction for day of year
    f = (279.575 + (0.9856 * dd)) * DEG2RAD
    e = (-104.7*sin(f) + 596.2*sin(2*f) + 4.3*sin(3*f) - 12.7*sin(4*f) &
        - 429.3*cos(f) - 2.0*cos(2*f) + 19.3*cos(3*f)) / 3600
    ! Calculate the longitudinal meridian
    lonm = nint(lon / 15.0) * 15.0
    ! Solar noon, with day of year and longitudinal correction
    LC = (lon - lonm) / 15
    solar_noon = 12 - LC - e
  end function solar_noon
 
  pure elemental logical function is_def_real(x)
    real, intent(in) :: x

    is_def_real = x > UNDEF
  end function is_def_real

  pure elemental logical function is_def_integer(x)
    integer, intent(in) :: x

    is_def_integer = x /= IUNDEF
  end function is_def_integer


  pure elemental logical function is_undef_real(x)
    real, intent(in) :: x

    is_undef_real = x <= UNDEF
  end function is_undef_real

  pure elemental logical function is_undef_integer(x)
    integer, intent(in) :: x

    is_undef_integer = x == IUNDEF
  end function is_undef_integer
  
  
  pure real function net_radiation(lat, lon, elev, albedo, dd, hr, sinB, R, Ts_C, eact)
    real, intent(in) :: lat     !< Latitude (degrees North)
    real, intent(in) :: lon     !< Longitude (degrees East)
    real, intent(in) :: elev    !< Elevation (m above sea level)
    real, intent(in) :: albedo  !< Surface albedo (fraction)
    integer, intent(in) :: dd   !< Day of year (1--365)
    integer, intent(in) :: hr   !< Hour of day (0--23)
    real, intent(in) :: sinB    !< sin() of solar elevation angle
    real, intent(in) :: R       !< Global radiation (W m-2)
    real, intent(in) :: Ts_C    !< Surface air temperature (degrees C)
    real, intent(in) :: eact    !< Actual vapour pressure (kPa)

    real, parameter :: Gsc = 0.082          ! Solar constant (MJ/m^2/min)
    real, parameter :: SBC = 4.903e-9 / 24  ! Stephan Boltzman constant

    real :: lat_rad, R_MJ, t0_, h, h1, h2, dr, dec, Re, pR, Rnl, Rns

    if (sinB <= 0) then
      net_radiation = 0.0
    else
      ! Latitude in radians
      lat_rad = lat * DEG2RAD

      ! Convert global radiation W m-2 to MJ m-2 s-1
      R_MJ = R * 0.0036

      ! Hour-angle of the sun
      t0_ = solar_noon(lon, dd)
      h = (15 * (hr - t0_)) * DEG2RAD
      h1 = h - (PI/24)
      h2 = h + (PI/24)

      dr = 1 + (0.033 * cos(((2 * PI) / 365) * dd))
      dec = solar_declination(dd)
      ! External radiation (with fix to stop div by zero)
      ! TODO: fix this to be less hackish
      Re = max(0.00000000001, &
               ((12 * 60) / PI) * Gsc * dr * ((h2 - h1) * sin(lat_rad) * sin(dec) &
               + cos(lat_rad) * cos(dec) * (sin(h2) - sin(h1))))
      ! TODO: what was this for?
      !Re = max(0.0, ((12*60)/PI)*Gsc*dr*sinB)

      ! Calculate net longwave radiation
      pR = (0.75 + (2e-5 * elev)) * Re

      Rnl = max(0.0, (SBC*((Ts_C + T0)**4)) * (0.34-(0.14*sqrt(eact))) &
                     * ((1.35*(min(1.0, R_MJ/pR)))-0.35))
      Rns = (1 - albedo) * R_MJ

      net_radiation = max(0.0, Rns - Rnl)
    end if
  end function net_radiation
  
  subroutine met_CO2(CO2_method, CO2_constant, CO2)
   CHARACTER(len=16), INTENT(in) :: CO2_method
   REAL, INTENT(in) :: CO2_constant
   REAL, INTENT(inout) :: CO2
    select case (CO2_method)
    case ("constant")
      CO2 = CO2_constant
    case ("input")
      !ASSERT_DEFINED(CO2)
    case default
      !UNKNOWN_STRING(config%CO2_method)
    end select
  end subroutine met_CO2
  
  
  logical function start_of_day(hr)
    INTEGER, INTENT(in) :: hr
    start_of_day = hr == 0
  end function start_of_day
  
  
  logical function is_daylight(R)
   real, INTENT(in) :: R
    is_daylight = R > 50.0
  end function is_daylight
  
  pure real function saturated_vapour_pressure(Ts_C)
    real, intent(in) :: Ts_C    !< Surface air temperature (degrees C)

    saturated_vapour_pressure = 0.611 * exp(17.27 * Ts_C / (Ts_C + 237.3))
  end function saturated_vapour_pressure
  
END MODULE
