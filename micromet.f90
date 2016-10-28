!
!////////////////////////////////////////////////////////////////////////
!
!      micromet.f90
!      Created: 12 August 2016 15:59 
!      By: David Gillies  
!
!////////////////////////////////////////////////////////////////////////
!
MODULE DO3SE_Micromet
USE DO3SE_Met
  !> Canopy displacement (fraction of canopy height)
  real, parameter :: CANOPY_D = 0.7
  !> Canopy roughness length (fraction of canopy height)
  real, parameter :: CANOPY_Z0 = 0.1
  
CONTAINS

subroutine calc_micromet_data(Idrctt, Idfuse, sinB, &
                       cosA, LAI, &
                       PARsun, PARshade, &
                       OTC, h_u, z_u, u_met, u_umet, u_50, ustar, canopy_height, &
                       LAIsunfrac)
    
    real, intent(in) :: Idrctt        !< Direct PAR irradiance (W m-2)
    real, intent(in) :: Idfuse        !< Diffuse PAR irradiance (W m-2)
    real, intent(in) :: sinB          !< sin() of solar elevation angle
    real, intent(in) :: cosA          !< cos(A), A = mean leaf inclination (0.5 = 60 degrees)
    real, intent(in) :: LAI           !< Leaf area index (m2 m-2)

    real, intent(out) :: PARsun       !< PAR received by sunlit leaves (W m-2)
    real, intent(out) :: PARshade     !< PAR received by shaded leaves (W m-2)
    LOGICAL, INTENT(in) :: OTC
    REAL, INTENT(inout) :: h_u
    REAL, INTENT(INout) :: z_u
    REAL, INTENT(in) :: u_met
    REAL, INTENT(out) :: u_umet
    REAL, INTENT(out) :: u_50
    REAL, INTENT(out) :: ustar
    real, intent(in) :: canopy_height  
    REAL, INTENT(out) :: LAIsunfrac
    
    call PAR_sun_shade(Idrctt, Idfuse, sinB, &
                       cosA, LAI, &
                       PARsun, PARshade)
    ! TODO: multi-layer PAR

    call met_windspeed(OTC, h_u, z_u, u_met, u_umet, u_50, ustar, canopy_height) ! h is canopy height

    ! Estimate sunlit LAI fractions (used later in f_light)
    ! TODO: remove the per-land-cover component of this?
    LAIsunfrac = MLMC_sunlit_LAI(LAI, sinB)
  end subroutine calc_micromet_data
  
  
  pure subroutine PAR_sun_shade(Idrctt, Idfuse, sinB, cosA, LAI, PARsun, PARshade)
    real, intent(in) :: Idrctt        !< Direct PAR irradiance (W m-2)
    real, intent(in) :: Idfuse        !< Diffuse PAR irradiance (W m-2)
    real, intent(in) :: sinB          !< sin() of solar elevation angle
    real, intent(in) :: cosA          !< cos(A), A = mean leaf inclination (0.5 = 60 degrees)
    real, intent(in) :: LAI           !< Leaf area index (m2 m-2)

    real, intent(out) :: PARsun       !< PAR received by sunlit leaves (W m-2)
    real, intent(out) :: PARshade     !< PAR received by shaded leaves (W m-2)

    if (sinB > 0.0) then
      ! PAR flux densities evaluated using method of Norman (1982, p.79):
      ! "conceptually, 0.07 represents a scattering coefficient"
      PARshade = Idfuse * exp(-0.5 * LAI**0.8) + &
        0.07 * Idrctt * (1.1 - (0.1 * LAI)) * exp(-sinB)
      PARsun = Idrctt * 0.8 * (cosA / sinB) + PARshade
    else
      PARshade = 0.0
      PARsun = 0.0
    end if
  end subroutine PAR_sun_shade
  
  subroutine met_windspeed(OTC, h_u, z_u, u_met, u_umet, u_50, ustar, canopy_height)
    LOGICAL, INTENT(in) :: OTC
    REAL, INTENT(inout) :: h_u
    REAL, INTENT(INout) :: z_u
    REAL, INTENT(in) :: u_met
    REAL, INTENT(out) :: u_umet
    REAL, INTENT(out) :: u_50
    REAL, INTENT(out) :: ustar
    real, intent(in) :: canopy_height                       !< Canopy height

    real, parameter :: MIN_WINDSPEED = 0.01 ! Minimum windspeed value (m s-1)
    real, parameter :: MIN_USTAR = 0.0001   ! Minimum friction velocity value (m s-1)


    real :: u_d, u_z0, d, z0
    real :: ustar_ref

    if (OTC) then
      h_u = canopy_height
      z_u = canopy_height
    else
      if (is_def(h_u)) then
        ! nothing to do
      else
        h_u = canopy_height
      end if
    end if

    u_d = h_u * CANOPY_D
    u_z0 = h_u * CANOPY_Z0
    d = canopy_height * CANOPY_D
    z0 = canopy_height * CANOPY_Z0

    ustar_ref = ustar_from_velocity(max(MIN_WINDSPEED, u), (z_u - u_d), u_z0)
    u_50 = max(MIN_WINDSPEED, velocity_from_ustar(ustar_ref, (50 - u_d), u_z0))
    ustar = max(MIN_USTAR, ustar_from_velocity(u_50, (50 - d), z0))
    u_umet = max(MIN_WINDSPEED, velocity_from_ustar(ustar, (canopy_height - d), z0))
  end subroutine met_windspeed
  
  pure elemental function ustar_from_velocity(u, z, z0) result (ustar)
    real, intent(in) :: u   ! Velocity at height above boundary (m/s)
    real, intent(in) :: z   ! Height above boundary, e.g. z - d (m)
    real, intent(in) :: z0  ! Roughness length, height at which u=0 (m)
    real :: ustar           ! Output: friction velocity, ustar (m/s)

    real, parameter :: K = 0.41 ! von Karman's constant

    ustar = (u * K) / log(z / z0)
  end function ustar_from_velocity

  pure elemental function velocity_from_ustar(ustar, z, z0) result (u)
    real, intent(in) :: ustar   ! Friction velocity (m/s)
    real, intent(in) :: z       ! Height above boundary, e.g. z - d (m)
    real, intent(in) :: z0      ! Roughness length, height at which u=0 (m)
    real :: u                   ! Output: velocity (m/s)

    real, parameter :: K = 0.41     ! von Karman's constant

    u = (ustar / K) * log(z / z0)
  end function velocity_from_ustar

  
  pure function MLMC_sunlit_LAI(LAI, sinB) result (LAIsunfrac)
    real, intent(in) :: LAI  !< Leaf area index (m^2/m^2)
    real, intent(in) :: sinB                 !< sin() of solar elevation angle

    real :: LAIsunfrac   !< Output: fraction of LAI that is sunlit

    REAL :: sunLAI_acc
    real :: sunLAI_layer
    real :: LAI_layer
    integer :: iL

    sunLAI_acc = 0.0

      ! How much of "canopy so far" is sunlit?
    sunLAI_acc = sunlit_LAI(LAI, sinB)
      ! How much of that is in this layer?
      
      ! Fraction of LAI which is sunlit
    LAI_layer = LAI
      if (LAI_layer > 0.0) then
        LAIsunfrac = sunLAI_acc / LAI_layer
      else
        LAIsunfrac = 0.0
      end if

  end function MLMC_sunlit_LAI
  
  
  pure real function sunlit_LAI(LAI, sinB)
    real, intent(in) :: LAI     !< Leaf area index (m^2/m^2)
    real, intent(in) :: sinB    !< sin() of solar elevation angle

    if (LAI > 0.0 .and. sinB > 0.0) then
      sunlit_LAI = ((1 - exp(-0.5 * LAI / sinB)) * (2 * sinB))
    else
      sunlit_LAI = 0.0
    end if
  end function sunlit_LAI
  
  
END MODULE
