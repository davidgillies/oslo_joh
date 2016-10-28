!
!////////////////////////////////////////////////////////////////////////
!
!      ozone_deposition.f90
!      Created: 15 August 2016 15:30 
!      By: David Gillies  
!
!////////////////////////////////////////////////////////////////////////
!
MODULE DO3SE_Ozonedep
USE DO3SE_Met

  real, parameter :: CANOPY_D = 0.7
  real, parameter :: CANOPY_Z0 = 0.1
  real, parameter :: DIFF_O3 = 0.000015
  
CONTAINS

  subroutine calc_ozone_deposition(ustar, canopy_height, SAI, bulk_gsto, O3_constant, &
                                   O3_met, O3_50_met, O3_umet, h_O3_loc, z_O3_loc, OTC, &
                                   Ra_c, Ra, Rb, Rinc, Rext, Rsto, Rgs, Rsur, Rtotal, LAI, &
                                   O3_method)
                             
    !class(DO3SE_State_t), intent(inout) :: this
    REAL, INTENT(in) :: ustar
    REAL,INTENT(in)  :: canopy_height
    REAL, INTENT(IN) :: SAI
    REAL, INTENT(in) :: bulk_gsto
    real,intent(out)  :: Ra_c
    real,intent(out)  :: Ra
    real,intent(out)  :: Rb
    real,intent(out)  :: Rinc
    real,intent(out)  :: Rext
    real,intent(out)  :: Rsto
    real,intent(out)  :: Rgs
    real,intent(out)  :: Rsur
    REAL,INTENT(out)  :: Rtotal
    REAL, INTENT(in) :: O3_constant
    REAL,INTENT(inout) :: O3_met
    REAL,INTENT(inout) :: O3_50_met
    REAL, INTENT(inout) :: O3_umet
    REAL, INTENT(inout) :: h_O3_loc
    REAL, INTENT(inout) :: z_O3_loc
    logical, INTENT(in) :: OTC
    REAL, INTENT(inout) :: LAI
    CHARACTER(len=16), INTENT(in) :: O3_method
    Ra_c = UNDEF
    Ra = UNDEF
    Rb = UNDEF
    Rinc = UNDEF
    Rext = UNDEF
    Rsto = UNDEF
    Rgs = UNDEF
    Rsur = UNDEF
    ! Calculate resistance model for O3 over the target canopy

    Ra_c = Ra_simple(ustar, (canopy_height * (CANOPY_D + CANOPY_Z0)), &
                               50.0, canopy_height * CANOPY_D)
    Ra = Ra_simple(ustar, canopy_height, 50.0, canopy_height * CANOPY_D)
    Rb = Rb_func(ustar, DIFF_O3)
    

    Rinc = Rinc_func(SAI, canopy_height, ustar)
      !this%V%rmodel_O3%Rinc(iL) = Rinc_prototype(sum(this%MLMC(iL,:)SAI), this%V%met%ustar)
    Rext = Rext_func(SAI)
    Rsto = Rsto_func(bulk_gsto)

    Rgs = Rsoil
    Rsur = Rsur_func(Rb, Rsto, Rext, LAI, SAI)
    Rtotal = Rtotal_func(Rsur, Rinc, Rgs)

    ! Vd calculation duplicated here just for comparison purposes
    ! TODO: better way to expose Vd
    Vd = deposition_velocity(Ra_c, Rtotal)

    ! Calculate O3 concentration at 50m and canopy height
    call met_O3(O3_method, O3_constant, O3_met, O3_50_met, O3_umet, h_O3_loc, &
                    z_O3_loc, Ra, Ra_c, Rtotal, canopy_height, OTC)
    
    
    ! Calculate per-layer O3 concentrations
    if (OTC) then
      ! In OTC, assume uniform O3 driven by external factors.
      !this%ML(:)%met%O3 = this%ML(1)%met%O3
    else
      ! Use multi-layer resistance model for per-layer O3 concentration.
      !call multi_layer_O3(this%V%rmodel_O3, this%V%met%O3_50, this%ML(:)%met%O3)
    end if
  end subroutine calc_ozone_deposition
  
  
  !> Calculate aerodynamic resistance (Ra, s m-1) between two heights using a
  !! simple, neutral stability model.
  !!
  !! Must satisfy \f$z_2 \leq z_1\f$, \f$z_2 \gt d\f$ and \f$z_1 \gt d\f$.
  pure real function Ra_simple(ustar, z1, z2, d)
    real, intent(in) :: ustar   !< Friction velocity (m/s)
    real, intent(in) :: z1      !< Lower height (m)
    real, intent(in) :: z2      !< Upper height (m)
    real, intent(in) :: d       !< Zero displacement height (m)

    real, parameter :: K = 0.41 ! von Karman's constant

    Ra_simple = (1.0 / (ustar * K)) * log((z2 - d) / (z1 - d))
  end function Ra_simple

  !> Calculate quasi-laminar boundary layer resistance (Rb, s m-1) based on a
  !! given friction velocity and diffusivity.
  pure real function Rb_func(ustar, diff)
    real, intent(in) :: ustar   !< Friction velocity (m s-1)
    real, intent(in) :: diff    !< Molecular diffusivity in air (m2 s-1)

    real, parameter :: PR = 0.72    ! Prandtl number
    real, parameter :: K = 0.41     ! von Karman's constant
    real, parameter :: V = 0.000015 ! Kinematic viscosity of air at 20 C (m2 s-1)

    Rb_func = (2.0 / (K * ustar)) * (((V/diff)/PR)**(2.0/3.0))
  end function Rb_func
  
  !> Estimate in-canopy aerodynamic resistance (Rinc, s m-1).
  !!
  !! This is the older single-layer DO3SE method.
  !!
  !! TODO: to use in a multilayer model, what does h represent?  Height above 
  !!       ground, or thickness of layer?
  pure real function Rinc_func(SAI, h, ustar)
    real, intent(in) :: SAI   !< Stand area index (m2 m-2)
    real, intent(in) :: h     !< Vegetation height (m)
    real, intent(in) :: ustar !< Friction velocity (m s-1)

    real, parameter :: Rinc_b = 14    ! Rinc coefficient

    Rinc_func = Rinc_b * SAI * h/ustar
  end function Rinc_func
  
  
  pure real function Rext_func(SAI)
    real, intent(in) :: SAI   !< Stand area index (m2 m-2)

    real, parameter :: Rext_base = 2500

    Rext_func = Rext_base / SAI
  end function Rext_func
  
  !> Convert stomatal conductance to stomatal resistance (Rsto, s m-1).
  !!
  !! The maximum stomatal resistance is capped to prevent infinite values
  !! when the conductance is 0.
  elemental real function Rsto_func(Gsto)
    real, intent(in) :: Gsto    !< Stomatal conductance (mmol m-2 s-1)

    real, parameter :: MAX_RSTO = 100000

    ! (gsto in m s-1) = 41000 * (gsto in mmol m-2 s-1)
    ! (rsto in s m-1) = 1 / (gsto in m s-1)
    Rsto_func = min(MAX_RSTO, 41000.0 / Gsto)
  end function Rsto_func
  
  
  !> Calculate per-layer surface resistance - combined Rb, Rsto and Rext.
  ! TODO: per-layer Rb
  pure function Rsur_func(Rb, Rsto, Rext, LAI, SAI)

    real, intent(in) :: Rb
    REAL, intent(in) :: Rsto
    REAL, intent(in) :: Rext
    REAL, intent(in) :: LAI
    REAL, intent(in) :: SAI

    real :: Rsur_func
    integer :: i


      ! TODO: let infinities happen and propagate through this? 1/Inf = 0?
      if (LAI > 0) then
        ! LAI (and SAI) > 0, include Rsto and Rext components
        Rsur = Rb + 1/(1/Rsto + 1/Rext)
      else if (SAI > 0) then
        ! Only SAI, omit the Rsto component
        Rsur = Rb + Rext
      else
        ! No foliage, very high resistance!
        ! TODO: find a justification for this, probably based on Rsto
        ! TODO: have an "R_INF" constant?
        Rsur_func = 1000000
      end if

  end function Rsur_func
  
  !> Calculate multi-layer Rtotal - the total resistance for each layer and
  !! everything below that layer.
  pure function Rtotal_func(Rsur, Rinc, Rgs)

    REAL, intent(in) :: Rsur
    real, intent(in) :: Rinc
    real, intent(in) :: Rgs

    real :: Rtotal_func

    real :: tmp
    integer :: i

    tmp = Rgs

    tmp = 1/(1/Rsur + 1/(Rinc + tmp))

    Rtotal_func = tmp
  end function Rtotal_func
  

  pure real function deposition_velocity(Ra_c, Rtotal)
    !type(ResistanceModel_t), intent(in) :: rmodel
    REAL, INTENT(in) :: Ra_c
    REAL, INTENT(in) :: Rtotal
    deposition_velocity = 1.0 / (Ra_c + Rtotal)
  end function deposition_velocity
  

    
  
  !> Define an ozone concentration at the canopy
  subroutine met_O3(O3_method, O3_constant, O3_met, O3_50_met, O3_umet, h_O3_loc, &
                    z_O3_loc, Ra, Ra_c, Rtotal, h, OTC)
    !type(MetConfig_t), intent(in) :: config         !< Met configuration
    !type(MetData_t), intent(inout) :: met           !< Met data
    !type(MicroMetData_t), intent(inout) :: umet     !< Top layer micromet data
    !type(Location_t), intent(in) :: loc
    !type(ResistanceModel_t), intent(in) :: rmodel   !< Resistance model for O3 over target canopy
    real, intent(in) :: O3_constant
    REAL, INTENT(inout) :: O3_met
    REAL, INTENT(inout) :: O3_50_met
    REAL, INTENT(inout) :: O3_umet
    REAL, INTENT(inout) :: h_O3_loc
    REAL, INTENT(inout) :: z_O3_loc
    REAL, INTENT(inout) :: Ra
    REAL, INTENT(inout) :: Ra_c
    REAL, INTENT(inout) :: Rtotal
    CHARACTER(len=16), INTENT(in) :: O3_method
    logical, INTENT(in) :: OTC
    real, intent(in) :: h                           !< Canopy height

    real :: h_O3
    real :: O3_d, O3_z0
    real :: ustar_ref
    !type(ResistanceModel_t) :: rmodel_ref

    select case (O3_method)
    case ("constant")
      !ASSERT_DEFINED(config%O3_constant)
      O3_met = O3_constant
    case ("offset")
      !ASSERT_DEFINED(config%O3_constant)
      !ASSERT_DEFINED(met%O3)
      O3_met = max(0.0, O3_met + O3_constant)
    case ("input")
      !ASSERT_DEFINED(met%O3)
    case default
      !UNKNOWN_STRING(config%O3_method)
    end select

    if (OTC) then
      O3_50_met = O3_met
      O3_umet = O3_met
    else
      if (is_def(h_O3_loc)) THEN ! h_O3 defined in location
        h_O3 = h_O3_loc
      else
        h_O3 = h
      end if

      O3_d = h_O3 * CANOPY_D
      O3_z0 = h_O3 * CANOPY_Z0

      ! TODO: this is unnecessary if h_O3 = h?
      ustar_ref = ustar_from_velocity(u_50_met, 50.0 - O3_d, O3_z0)
      !rmodel_ref = rmodel
      ref_Ra_c = Ra_simple(ustar_ref, O3_d + O3_z0, 50.0, O3_d)
      ref_Ra = Ra_simple(ustar_ref, z_O3_loc, 50.0, O3_d)
      ref_Rb = Rb_func(ustar_ref, DIFF_O3)

      O3_50_met = O3_transfer_up(ref_Ra, ref_Ra_c, ref_Rb, O3_met)
      O3_umet = O3_transfer_down(Ra, Ra_c, Rtotal, O3_50_met)
    end if
  end subroutine met_O3
  
  !> Estimate friction velocity from windspeed velocity.
  pure elemental function ustar_from_velocity(u, z, z0) result (ustar)
    real, intent(in) :: u   ! Velocity at height above boundary (m/s)
    real, intent(in) :: z   ! Height above boundary, e.g. z - d (m)
    real, intent(in) :: z0  ! Roughness length, height at which u=0 (m)
    real :: ustar           ! Output: friction velocity, ustar (m/s)

    real, parameter :: K = 0.41 ! von Karman's constant

    ustar = (u * K) / log(z / z0)
  end function ustar_from_velocity
  
  !> Calculate aerodynamic resistance (Ra, s m-1) between two heights using a
  !! simple, neutral stability model.
  !!
  !! Must satisfy \f$z_2 \leq z_1\f$, \f$z_2 \gt d\f$ and \f$z_1 \gt d\f$.

  
  

  !> Scale O3 concentration up from the resistance model's reference height
  !! to 50m.
  pure real function O3_transfer_up(Ra, Ra_c, Rtotal, O3)
    !type(ResistanceModel_t), intent(in) :: rmodel
    real, intent(in) :: O3
    REAL,INTENT(in) :: Ra
    REAL,INTENT(in) :: Ra_c
    REAL,INTENT(in) :: Rtotal
    real :: Vd

    Vd = deposition_velocity(Ra_c, Rtotal)
    O3_transfer_up = O3 / (1.0 - (Ra * Vd))
  end function O3_transfer_up
  
  !> Scale O3 concentration down from 50m to the resistance model's reference
  !! height.
  pure real function O3_transfer_down(Ra, Ra_c, Rtotal, O3_50)
    !type(ResistanceModel_t), intent(in) :: rmodel
    real, intent(in) :: O3_50
    REAL,INTENT(in) :: Ra
    REAL,INTENT(in) :: Ra_c
    REAL,INTENT(in) :: Rtotal
    real :: Vd

    Vd = deposition_velocity(Ra_c, Rtotal)
    O3_transfer_down = O3_50 * (1.0 - (Ra * Vd))
  end function O3_transfer_down
END MODULE
