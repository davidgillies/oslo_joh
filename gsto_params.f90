!
!////////////////////////////////////////////////////////////////////////
!
!      gsto_params.f90
!      Created: 15 August 2016 09:44 
!      By: David Gillies  
!
!////////////////////////////////////////////////////////////////////////
!
MODULE DO3SE_Gsto_params
USE DO3SE_Phenology

    REAL :: K_z = 24.0      !< Coefficient for ozone damage (dimensionless)
    real, parameter :: PAR_Wm2_to_photons = 4.57
CONTAINS

  subroutine calc_gsto_parameters(f_phen_method, SGS, EGS, dd, f_phen_a, f_phen_b, &
                                  f_phen_c, f_phen_d, f_phen_e, f_phen_1, f_phen_2, &
                                  f_phen_3, f_phen_4, f_phen_limA, f_phen_limB, &
                                  leaf_f_phen_method, leaf_f_phen_a, leaf_f_phen_b, &
                                  leaf_f_phen_c,  leaf_f_phen_1, leaf_f_phen_2, &
                                  Astart, Aend, f_light_method, f_lightfac, & 
                                  PARsun, PARshade, LAIsunfrac, PAR, &
                                  f_temp_method, Ts_C, T_min, T_opt, T_max, fmin, &
                                  f_VPD_method, VPD, VPD_max, VPD_min, f_SW_method, &
                                  fSWP_exp_a, fSWP_exp_b, SWP, SWP_min, SWP_max, &
                                  ASW_FC, ASW, f_O3_method, POD_0, AOT_0, O3_method, &
                                  FO3_eff, V_cmax_25, J_max_25, phenology_method, &
                                  f_phen, leaf_f_phen, f_light, leaf_f_light, &
                                  f_temp, f_VPD, f_SW, f_O3)
    !class(DO3SE_State_t), intent(inout) :: this
    !type(LandCover_t), intent(in) :: LC
    !type(V_t), intent(in) :: V
    !type(ML_t), intent(in) :: ll
    !type(MLMC_t), intent(inout) :: xx
    CHARACTER(len=16), INTENT(in) :: f_phen_method
    integer, intent(in) :: SGS            !< Start of growing season (day of year)
    integer, intent(in) :: EGS            !< End of growing season (day of year)
    integer, intent(in) :: dd             !< Day of year
    REAL, INTENT(in) :: f_phen_a
    REAL, INTENT(in) :: f_phen_b
    REAL, INTENT(in) :: f_phen_c
    REAL, INTENT(in) :: f_phen_d
    REAL, INTENT(in) :: f_phen_e
    INTEGER, INTENT(in) :: f_phen_1
    INTEGER, INTENT(in) :: f_phen_2
    INTEGER, INTENT(in) :: f_phen_3
    INTEGER, INTENT(in) :: f_phen_4
    INTEGER, INTENT(in) :: f_phen_limA
    INTEGER, INTENT(in) :: f_phen_limB
    CHARACTER(len=16), INTENT(in) :: leaf_f_phen_method
    REAL, INTENT(in) :: leaf_f_phen_a       !< f_phen at Astart
    REAL, INTENT(in) :: leaf_f_phen_b       !< f_phen at mid-season peak
    REAL, INTENT(in) :: leaf_f_phen_c       !< f_phen at Aend
    INTEGER, INTENT(in) :: leaf_f_phen_1    !< Time from _a to _b (days)
    INTEGER, INTENT(in) :: leaf_f_phen_2    !< Time from _b to _c (days)
    integer, intent(in) :: Astart         !< Start of accumulation period (day of year)
    integer, intent(in) :: Aend           !< End of accumulation period (day of year)
    CHARACTER(len=16), INTENT(in) :: f_light_method
    real, intent(in) :: f_lightfac        !< Single leaf f_light coefficient
    real, intent(in) :: PARsun            !< PAR received by sunlit leaves (W m-2)
    real, intent(in) :: PARshade          !< PAR received by shaded leaves (W m-2)
    real, intent(in) :: LAIsunfrac        !< Fraction of canopy component that is sunlit
    real, intent(in) :: PAR               !< Photosynthetically active radiation (W m-2)    
    CHARACTER(len=16), INTENT(in) :: f_temp_method
    real, intent(in) :: Ts_C    !< Air temperature (degrees C)
    real, intent(in) :: T_min   !< Minimum temperature (degrees C)
    real, intent(in) :: T_opt   !< Optimum temperature (degrees C)
    real, intent(in) :: T_max   !< Maximum temperature (degrees C)
    real, intent(in) :: fmin    !< Minimum f_temp
    CHARACTER(len=16), INTENT(in) :: f_VPD_method
    real, intent(in) :: VPD       !< Vapour pressure deficit (kPa)
    real, intent(in) :: VPD_max   !< VPD for maximum gsto (kPa)
    real, intent(in) :: VPD_min   !< VPD for minimum gsto (kPa)
    CHARACTER(len=16), INTENT(in) :: f_SW_method
    real, intent(in) :: fSWP_exp_a
    real, intent(in) :: fSWP_exp_b
    real, intent(in) :: SWP       !< Soil water potential (MPa)
    real, intent(in) :: SWP_min   !< SWP for minimum gsto (MPa)
    real, intent(in) :: SWP_max   !< SWP for maximum gsto (MPa)
    real, intent(in) :: ASW_FC    !< Available soil water at field capacity (m)
    real, intent(in) :: ASW       !< Available soil water (m)
    CHARACTER(len=16), INTENT(in) :: f_O3_method
    REAL, INTENT(in) :: POD_0
    REAL, INTENT(in) :: AOT_0
    CHARACTER(len=16), INTENT(in) :: O3_method
    REAL, INTENT(in) :: FO3_eff
    REAL, INTENT(inout) :: V_cmax_25
    REAL, INTENT(inout) :: J_max_25
    CHARACTER(len=16), INTENT(in) :: phenology_method
    REAL, INTENT(out) :: f_phen
    REAL, INTENT(out) :: leaf_f_phen
    REAL, INTENT(out) :: f_light
    REAL, INTENT(out) :: leaf_f_light
    REAL, INTENT(out) :: f_VPD
    REAL, INTENT(out) :: f_SW
    REAL, INTENT(out) :: f_O3
    
    
    
    
    real :: delta_V_cmax, mult_V_cmax

    !associate (season => LC%season, &
    !           gc => LC%gsto, &
     !          pgc => LC%pn_gsto, &
     !          gp => xx%gsto_params)
      ! Initialise gsto parameters from configuration
    !gp = GstoParams_t(fmin=gc%fmin, gmax=gc%gmax, gmorph=gc%gmorph)

      ! Calculate f_phen
    select case (f_phen_method)
      case ("disabled")
        ! Nothing to do
      case ("simple day PLF")
        f_phen = f_phen_simple_PLF(SGS, EGS, dd, f_phen_a, f_phen_c, &
                                  f_phen_e, f_phen_1, f_phen_4)
      case ("complex day PLF")
        f_phen = f_phen_complex_PLF(SGS, EGS, dd, f_phen_a, f_phen_b, f_phen_c, &
                                  f_phen_d, f_phen_e, f_phen_1, f_phen_2, &
                                  f_phen_3, f_phen_4, f_phen_limA, f_phen_limB)
      case default
        !UNKNOWN_STRING(gc%f_phen_method)
      end select

      ! Calculate leaf_f_phen
      select case (leaf_f_phen_method)
      case ("disabled")
        ! Nothing to do
      case ("f_phen")
        leaf_f_phen = f_phen
      case ("day PLF")
        leaf_f_phen = leaf_f_phen_PLF(leaf_f_phen_a, leaf_f_phen_b, leaf_f_phen_c, &
                                leaf_f_phen_1, leaf_f_phen_2, Astart, Aend, dd)
      case default
        !UNKNOWN_STRING(gc%leaf_f_phen_method)
      end select

      select case (f_light_method)
      case ("disabled")
        ! Nothing to do
      case ("enabled")
        if (LAI > 0 .and. sinB > 0) then
          ! Calculate f_light and leaf_f_light
          ! TODO: attenuate PAR properly through the canopy
          f_light = f_light_func(f_lightfac, PARsun, PARshade, LAIsunfrac)
          ! TODO: "grassland multilayer" model used leaf_flight = Flightsun, i.e.
          !       leaf_f_light(gc%f_lightfac, ll%met%PARsun) - which version is right?
          leaf_f_light = leaf_f_light_func(f_lightfac, PAR)
        else
          f_light = 0.0
          leaf_f_light = 0.0
        end if
      case default
        !UNKNOWN_STRING(gc%f_light_method)
      end select

      ! Calculate f_temp
      select case (f_temp_method)
      case ("disabled")
        ! Nothing to do
      case ("default")
        f_temp = f_temp_func(Ts_C, T_min, T_opt, T_max, fmin)
      case ("square high")
        f_temp = f_temp_square_high(Ts_C, T_min, T_opt, T_max, fmin)
      case default
        !UNKNOWN_STRING(gc%f_temp_method)
      end select

      ! Calculate f_VPD
      select case (f_VPD_method)
      case ("disabled")
        ! Nothing to do
      case ("linear")
        f_VPD = f_VPD_linear(VPD, VPD_max, VPD_min, fmin)
      case ("log")
        f_VPD = f_VPD_log(VPD, fmin)
      case default
        !UNKNOWN_STRING(gc%f_VPD_method)
      end select

      ! Calculate f_SW
      select case (f_SW_method)
      case ("disabled")
        ! Nothing to do
      case ("fSWP exp")
        f_SW = f_SWP_exp(fSWP_exp_a, SWP_exp_b, fmin, SWP)
      case ("fSWP linear")
        f_SW = f_SWP_linear(SWP_min, SWP_max, fmin, SWP)
      ! TODO: implement LWP
      !case ("fLWP exp")
      !  gp%f_SW = f_SWP_exp(gc%fSWP_exp_a, gc%fSWP_exp_b, gc%fmin, V%SMD%LWP)
      case ("fPAW")
        f_SW = f_PAW(ASW_FC, fmin, ASW)
      case default
        !UNKNOWN_STRING(gc%f_SW_method)
      end select

      ! Calculate f_O3
      select case (f_O3_method)
      case ("disabled")
        ! Nothing to do
      case ("wheat")
        f_O3 = ((1+(POD_0/11.5)**10)**(-1))
      case ("potato")
        f_O3 = ((1+(AOT_0/40)**5)**(-1))
      case default
        !UNKNOWN_STRING(gc%f_O3_method)
      end select

      ! Photosynthetic gsto: V_cmax_25 and J_max_25
      !!select case (V_J_method)
      !case ("input")
        ! Nothing to do
      !case ("constant")
      ! gets it from config, already known
        !xx%V_cmax_25 = pgc%V_cmax_25
        !xx%J_max_25 = pgc%J_max_25
      !case default
        !UNKNOWN_STRING(pgc%V_J_method)
      !end select

      ! Photosynthetic gsto: O3 effect on V_cmax_25
      select case (O3_method)
      case ("disabled")
        ! Nothing to do
      case ("martin2000")
        ! Percentage reduction in V_cmax (converting FO3_eff from nmol to mmol)
        delta_V_cmax = K_z * (FO3_eff / 1e6)
        ! Convert to multiplier
        mult_V_cmax = 1.0 - (delta_V_cmax / 100)
        ! Reduce V_cmax and J_max
        V_cmax_25 = max(0.1 * V_cmax_25, min(V_cmax_25 * mult_V_cmax, V_cmax_25))
        J_max_25 = max(0.1 * J_max_25, min(J_max_25 * mult_V_cmax, J_max_25))
      case default
        !UNKNOWN_STRING(pgc%V_J_method)
      end select

      ! Photosynthetic gsto: phenology effect
      select case (phenology_method)
      case ("disabled")
        ! Nothing to do
      case ("leaf_f_phen")
        V_cmax_25 = V_cmax_25 * leaf_f_phen
        J_max_25 =  J_max_25 * leaf_f_phen
      case default
        !UNKNOWN_STRING(pgc%phenology_method)
      end select

  end subroutine calc_gsto_parameters

  real function f_phen_simple_PLF(SGS, EGS, dd, f_phen_a, f_phen_c, &
                                  f_phen_e, f_phen_1, f_phen_4)
    !type(GstoConfig_t), intent(in) :: gc  !< gsto parameters
    
    integer, intent(in) :: SGS            !< Start of growing season (day of year)
    integer, intent(in) :: EGS            !< End of growing season (day of year)
    integer, intent(in) :: dd             !< Day of year
    REAL, INTENT(in) :: f_phen_a
    REAL, INTENT(in) :: f_phen_c
    REAL, INTENT(in) :: f_phen_e
    INTEGER, INTENT(in) :: f_phen_1
    INTEGER, INTENT(in) :: f_phen_4
    
    
    real, dimension(2, 6) :: func
    integer :: dd_adj

    ! Build function
    func = reshape((/ real :: &
      SGS, 0.0, &
      SGS, f_phen_a, &
      (SGS + f_phen_1), f_phen_c, &
      (EGS - f_phen_4), f_phen_c, &
      EGS, f_phen_e, &
      EGS, 0.0 /), shape(func))
    ! Re-index everything to SGS = 0, wrapping dates before SGS to the end of the year
    call reindex(func(1,:), real(SGS), 365.0)
    !call assert(all(func(1,1:5) <= func(1,2:6)), "f_phen_simple_PLF: points not in order")
    dd_adj = dd - SGS
    if (dd_adj < 0) then
      dd_adj = dd_adj + 365
    end if
    ! Lookup value in PLF
    f_phen_simple_PLF = PLF_value(func, real(dd_adj))
  end function f_phen_simple_PLF


  
  real function f_phen_complex_PLF(SGS, EGS, dd, f_phen_a, f_phen_b, f_phen_c, &
                                  f_phen_d, f_phen_e, f_phen_1, f_phen_2, &
                                  f_phen_3, f_phen_4, f_phen_limA, f_phen_limB)
    !type(GstoConfig_t), intent(in) :: gc  !< gsto parameters
    integer, intent(in) :: SGS            !< Start of growing season (day of year)
    integer, intent(in) :: EGS            !< End of growing season (day of year)
    integer, intent(in) :: dd             !< Day of year
    REAL, INTENT(in) :: f_phen_a
    REAL, INTENT(in) :: f_phen_b
    REAL, INTENT(in) :: f_phen_c
    REAL, INTENT(in) :: f_phen_d
    REAL, INTENT(in) :: f_phen_e
    INTEGER, INTENT(in) :: f_phen_1
    INTEGER, INTENT(in) :: f_phen_2
    INTEGER, INTENT(in) :: f_phen_3
    INTEGER, INTENT(in) :: f_phen_4
    INTEGER, INTENT(in) :: f_phen_limA
    INTEGER, INTENT(in) :: f_phen_limB
    real, dimension(2, 10) :: func
    integer :: dd_adj

    ! Build function
    func = reshape((/ real :: &
      SGS, 0.0, &
      SGS, f_phen_a, &
      (SGS + f_phen_1), f_phen_b, &
      f_phen_limA, f_phen_b, &
      (f_phen_limA + f_phen_2), f_phen_c, &
      (f_phen_limB - f_phen_3), f_phen_c, &
      f_phen_limB, f_phen_d, &
      (EGS - f_phen_4), f_phen_d, &
      EGS, f_phen_e, &
      EGS, 0.0 /), shape(func))
    ! Re-index everything to SGS = 0, wrapping dates before SGS to the end of the year
    call reindex(func(1,:), real(SGS), 365.0)
    !call assert(all(func(1,1:9) <= func(1,2:10)), "f_phen_complex_PLF: points not in order")
    dd_adj = dd - SGS
    if (dd_adj < 0) then
      dd_adj = dd_adj + 365
    end if
    ! Lookup value in PLF
    f_phen_complex_PLF = PLF_value(func, real(dd_adj))
  end function f_phen_complex_PLF
  
  real function leaf_f_phen_PLF(leaf_f_phen_a, leaf_f_phen_b, leaf_f_phen_c, &
                                leaf_f_phen_1, leaf_f_phen_2, Astart, Aend, dd)
    !type(GstoConfig_t), intent(in) :: gc  !< gsto parameters
    REAL, INTENT(in) :: leaf_f_phen_a       !< f_phen at Astart
    REAL, INTENT(in) :: leaf_f_phen_b       !< f_phen at mid-season peak
    REAL, INTENT(in) :: leaf_f_phen_c       !< f_phen at Aend
    INTEGER, INTENT(in) :: leaf_f_phen_1    !< Time from _a to _b (days)
    INTEGER, INTENT(in) :: leaf_f_phen_2    !< Time from _b to _c (days)
    integer, intent(in) :: Astart         !< Start of accumulation period (day of year)
    integer, intent(in) :: Aend           !< End of accumulation period (day of year)
    integer, intent(in) :: dd             !< Day of year

    real, dimension(2, 6) :: func
    integer :: dd_adj

    ! Build function
    func = reshape((/ real :: &
      Astart, 0.0, &
      Astart, leaf_f_phen_a, &
      (Astart + leaf_f_phen_1), leaf_f_phen_b, &
      (Aend - leaf_f_phen_2), leaf_f_phen_b, &
      Aend, leaf_f_phen_c, &
      Aend, 0.0 /), shape(func))
    ! Re-index everything to SGS = 0, wrapping dates before Astart to the end of the year
    call reindex(func(1,:), real(Astart), 365.0)
    !call assert(all(func(1,1:5) <= func(1,2:6)), "leaf_f_phen_PLF: points not in order")
    dd_adj = dd - Astart
    if (dd_adj < 0) then
      dd_adj = dd_adj + 365
    end if
    ! Lookup value in PLF
    leaf_f_phen_PLF = PLF_value(func, real(dd_adj))
  end function leaf_f_phen_PLF
  
  
  pure real function f_light_func(f_lightfac, PARsun, PARshade, LAIsunfrac)
    real, intent(in) :: f_lightfac        !< Single leaf f_light coefficient
    real, intent(in) :: PARsun            !< PAR received by sunlit leaves (W m-2)
    real, intent(in) :: PARshade          !< PAR received by shaded leaves (W m-2)
    real, intent(in) :: LAIsunfrac        !< Fraction of canopy component that is sunlit

    real :: Flightsun, Flightshade

    ! TODO: does this need albedo?
    Flightsun = leaf_f_light_func(f_lightfac, PARsun)
    Flightshade = leaf_f_light_func(f_lightfac, PARshade)

    f_light_func = LAIsunfrac * Flightsun + (1.0 - LAIsunfrac) * Flightshade
  end function f_light_func
  
  pure real function leaf_f_light_func(f_lightfac, PAR)
    real, intent(in) :: f_lightfac        !< Single leaf f_light coefficient
    real, intent(in) :: PAR               !< Photosynthetically active radiation (W m-2)
    leaf_f_light_func = 1.0 - exp(-f_lightfac * (PAR * PAR_Wm2_to_photons))
  end function leaf_f_light_func
  
  ! pure real function leaf_f_light(f_lightfac, PAR)
  !  real, intent(in) :: f_lightfac        !< Single leaf f_light coefficient
  !  real, intent(in) :: PAR               !< Photosynthetically active radiation (W m-2)!
!
 !   leaf_f_light = 1.0 - exp(-f_lightfac * (PAR * PAR_Wm2_to_photons))
  !end function leaf_f_light
  
  pure real function f_temp_func(Ts_C, T_min, T_opt, T_max, fmin)
    real, intent(in) :: Ts_C    !< Air temperature (degrees C)
    real, intent(in) :: T_min   !< Minimum temperature (degrees C)
    real, intent(in) :: T_opt   !< Optimum temperature (degrees C)
    real, intent(in) :: T_max   !< Maximum temperature (degrees C)
    real, intent(in) :: fmin    !< Minimum f_temp

    real :: bt

    bt = (T_max - T_opt) / (T_opt - T_min)
    f_temp = ((Ts_C - T_min) / (T_opt - T_min)) * ((T_max - Ts_C) / (T_max - T_opt))**bt
    f_temp = max(fmin, min(1.0, f_temp))
  end function f_temp_func


  pure real function f_temp_square_high(Ts_C, T_min, T_opt, T_max, fmin)
    real, intent(in) :: Ts_C    !< Air temperature (degrees C)
    real, intent(in) :: T_min   !< Minimum temperature (degrees C)
    real, intent(in) :: T_opt   !< Optimum temperature (degrees C)
    real, intent(in) :: T_max   !< Maximum temperature (degrees C)
    real, intent(in) :: fmin    !< Minimum f_temp

    if (Ts_C >= T_max) then
      f_temp_square_high = 0.0
    else if (Ts_C >= T_opt) then
      f_temp_square_high = 1.0
    else
      f_temp_square_high = f_temp_func(Ts_C, T_min, T_opt, T_max, fmin)
    end if
  end function f_temp_square_high
  
  pure real function f_VPD_linear(VPD, VPD_max, VPD_min, fmin)
    real, intent(in) :: VPD       !< Vapour pressure deficit (kPa)
    real, intent(in) :: VPD_max   !< VPD for maximum gsto (kPa)
    real, intent(in) :: VPD_min   !< VPD for minimum gsto (kPa)
    real, intent(in) :: fmin      !< Minimum f_VPD

    real :: slope

    slope = (1.0 - fmin) / (VPD_min - VPD_max)
    f_VPD_linear = max(fmin, min(1.0, (fmin + (VPD_min - VPD)*slope)))
  end function f_VPD_linear
  
  pure real function f_VPD_log(VPD, fmin)
    real, intent(in) :: VPD       !< Vapour pressure deficit (kPa)
    real, intent(in) :: fmin      !< Minimum f_VPD

    f_VPD_log = max(fmin, min(1.0, (1.0 - 0.6*log(VPD))))
  end function f_VPD_log
  
  pure real function f_SWP_exp(a, b, fmin, SWP)
    real, intent(in) :: a
    real, intent(in) :: b
    real, intent(in) :: fmin      !< Minimum f_SWP
    real, intent(in) :: SWP       !< Soil water potential (MPa)

    f_SWP_exp = min(1.0, max(fmin, a * (-SWP)**b))
  end function f_SWP_exp
  
  pure real function f_SWP_linear(SWP_min, SWP_max, fmin, SWP)
    real, intent(in) :: SWP_min   !< SWP for minimum gsto (MPa)
    real, intent(in) :: SWP_max   !< SWP for maximum gsto (MPa)
    real, intent(in) :: fmin      !< Minimum f_SWP
    real, intent(in) :: SWP       !< Soil water potential (MPa)

    f_SWP_linear = fmin + (1-fmin) * ((SWP_min - SWP)/(SWP_min - SWP_max))
    f_SWP_linear = min(1.0, max(fmin, f_SWP_linear))
  end function f_SWP_linear
  
  pure real function f_PAW(ASW_FC, fmin, ASW)
    real, intent(in) :: ASW_FC    !< Available soil water at field capacity (m)
    real, intent(in) :: fmin      !< Minimum f_PAW
    real, intent(in) :: ASW       !< Available soil water (m)

    f_PAW = fmin + (1-fmin) * ((100 * (ASW/ASW_FC)) - ASW_MIN)/(ASW_MAX - ASW_MIN)
    f_PAW = min(1.0, max(fmin, f_PAW))
  end function f_PAW

END MODULE
