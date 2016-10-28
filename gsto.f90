!
!////////////////////////////////////////////////////////////////////////
!
!      gsto.f90
!      Created: 15 August 2016 13:40 
!      By: David Gillies  
!
!////////////////////////////////////////////////////////////////////////
!
MODULE DO3SE_Gsto
USE DO3SE_Photosynthesis
CONTAINS

                          
                          
  subroutine calc_gsto(method, VPD_crit, VPD_dd, leaf_gsto, gmax, &
                                   leaf_f_phen, f_O3, leaf_f_light, fmin, &
                                   f_temp, f_VPD, f_SW, mean_gsto, gmorph, &
                                   f_phen, f_light, VPD_max, VPD_min, &
                                   Tleaf_method, g_sto_0, m, V_cmax_25, J_max_25, &
                                   Lm, Tair_C, u, CO2, &
                                   PPFD, Rn, R, albedo, P, eact, &
                                   Tleaf_C, LAI, bulk_gsto, D_0_method, f_VPD_method, &
                                   Ts_C, Tleaf_balance_threshold, &
                                     Tleaf_adjustment_factor, &
                                     Tleaf_max_iterations, A_n, g_sv, g_bv)

    !class(DO3SE_State_t), intent(inout) :: this
    !type(Location_t), intent(in) :: loc
    !type(LandCover_t), intent(in) :: LC
    !type(V_t), intent(in) :: V
    !type(ML_t), intent(inout) :: ll
    !type(MLMC_t), intent(inout) :: xx
    CHARACTER(len=16), INTENT(in) :: method
    REAL, INTENT(in) :: gmax
    REAL, INTENT(in) :: gmorph
    REAL, INTENT(in) :: fmin
    REAL, INTENT(in) :: VPD_crit
    REAL, INTENT(in) :: VPD_dd
    REAL, INTENT(inout) :: leaf_gsto
    REAL, INTENT(inout) :: mean_gsto
    REAL, INTENT(inout) :: bulk_gsto
    REAL, INTENT(in) :: leaf_f_phen
    REAL, INTENT(in) :: f_O3
    REAL, INTENT(in) :: leaf_f_light
    REAL, INTENT(IN) :: f_temp
    REAL, INTENT(in) :: f_VPD
    REAL, INTENT(in) :: f_SW
    REAL, INTENT(IN) :: f_phen
    REAL, INTENT(in) :: f_light
    REAL, INTENT(in) :: V_cmax_25
    REAL, INTENT(in) :: J_max_25
    REAL, INTENT(in) :: eact
    REAL, INTENT(in) :: P
    REAL, INTENT(in) :: albedo
    REAL, INTENT(in) :: R
    REAL, INTENT(in) :: Rn
    REAL, INTENT(IN) :: PPFD
    REAL, INTENT(in) :: CO2
    REAL, INTENT(IN) :: u
    REAL, INTENT(in) :: Lm
    REAL, INTENT(in) :: LAI
    CHARACTER(len=16), INTENT(in) :: D_0_method
    CHARACTER(len=16), INTENT(in) :: f_VPD_method
    REAL, INTENT(in) :: VPD_max
    REAL, INTENT(in) :: VPD_min
    CHARACTER(len=16), INTENT(in) :: Tleaf_method
    REAL, INTENT(in) :: g_sto_0
    REAL,INTENT(IN) :: m
    REAL, INTENT(out) :: A_n        !< Output: Net CO2 assimilation (umol m-2 PLA s-1)
    real :: g_sto      !< Output: Stomatal conductance (mmol m-2 PLA s-1)
    REAL, INTENT(out) :: g_sv       !< Output: Stomatal conductance to water vapour
    REAL, INTENT(out) :: g_bv       !< Output: Boundary conductance to water vapour
    real, intent(inout) :: Tleaf_C  !< Leaf temperature (degrees C)
    !real, intent(inout) :: Tair_C
    REAL, INTENT(in) :: Ts_C
    real, INTENT(in) :: Tleaf_balance_threshold
    real, INTENT(in) :: Tleaf_adjustment_factor
    INTEGER, INTENT(in) :: Tleaf_max_iterations
    
    real :: D_0
    

    !associate (gc => LC%gsto, &
    !           pgc => LC%pn_gsto, &
    !           gp => xx%gsto_params)
      ! Calculate gsto
      select case (method)
      case ("multiplicative")
        leaf_gsto = apply_VPD_crit(VPD_crit, VPD_dd, leaf_gsto, gsto_leaf(gmax, &
                                   leaf_f_phen, f_O3, leaf_f_light, fmin, &
                                   f_temp, f_VPD, f_SW))
        mean_gsto = apply_VPD_crit(VPD_crit, VPD_dd, mean_gsto, gsto_mean(gmax, gmorph, &
                                   f_phen, f_light, fmin, f_temp, &
                                   f_VPD, f_SW))
      case ("photosynthesis")
        select case (D_0_method)
        case ("constant")
          !D_0 = pgc%D_0
        case ("f_VPD")
          select case (f_VPD_method)
          case ("linear")
            D_0 = inverse_f_VPD_linear(0.5, VPD_max, VPD_min, fmin)
          case ("log")
            D_0 = inverse_f_VPD_log(0.5, fmin)
          end select
        end select

        ! TODO: make this multi-layer aware
        ! TODO: Tleaf_C isn't multi-component aware, and will get overwritten
        !       repeatedly!  Fixing this should be able to make `ll` intent(in)
         call gsto_pn(Tleaf_method, g_sto_0, m, V_cmax_25, J_max_25, D_0, &
                     Lm, Ts_C, u, CO2, PPFD, &
                     Rn, R, albedo, P, eact, &
                    Tleaf_C, A_n, leaf_gsto, g_sv, g_bv, Tleaf_balance_threshold, &
                                     Tleaf_adjustment_factor, &
                                     Tleaf_max_iterations)
        mean_gsto = leaf_gsto
      case default
        !UNKNOWN_STRING(gc%method)
      end select
      ! Scale mean gsto up to bulk gsto
      bulk_gsto = mean_gsto * LAI

  end subroutine


  pure real function apply_VPD_crit(VPD_crit, VPD_dd, old_gsto, new_gsto)
    real, intent(in) :: VPD_crit    !< Accumulated VPD threshold (kPa)
    real, intent(in) :: VPD_dd      !< Accumulated VPD for the day (kPa)
    real, intent(in) :: old_gsto    !< Previous stomatal conductance
    real, intent(in) :: new_gsto    !< New stomatal conductance

    if (VPD_dd >= VPD_crit) then
      apply_VPD_crit = min(old_gsto, new_gsto)
    else
      apply_VPD_crit = new_gsto
    end if
  end function apply_VPD_crit
  
  !> Mean stomatal conductance of top leaf (mmol O3 m-2 PLA s-1).
  pure real function gsto_leaf(gmax, leaf_f_phen, f_O3, leaf_f_light, fmin, &
                               f_temp, f_VPD, f_SW)
    !type(GstoParams_t), intent(in) :: gp    !< Multiplicative gsto parameters
    REAL, INTENT(in) :: gmax
    REAL, INTENT(in) :: leaf_f_phen
    REAL, INTENT(in) :: f_O3
    REAL, INTENT(in) :: leaf_f_light
    REAL, INTENT(in) :: fmin
    REAL, INTENT(in) :: f_temp
    REAL, INTENT(in) :: f_VPD
    REAL, INTENT(in) :: f_SW
    gsto_leaf = gmax * min(leaf_f_phen, f_O3) * leaf_f_light * &
                max(fmin, f_temp * f_VPD * f_SW)
  end function gsto_leaf
  
  pure real function gsto_mean(gmax, gmorph, f_phen, f_light, fmin, f_temp, &
                               f_VPD, f_SW)
    !type(GstoParams_t), intent(in) :: gp    !< Multiplicative gsto parameters
    REAL, INTENT(in) :: gmax
    REAL, INTENT(in) :: gmorph
    REAL, INTENT(in) :: f_phen
    REAL, INTENT(in) :: f_light
    REAL, INTENT(in) :: fmin
    REAL, INTENT(in) :: f_temp
    REAL, INTENT(in) :: f_VPD
    REAL, INTENT(in) :: f_SW

    gsto_mean = gmax * gmorph * f_phen * f_light * &
                max(fmin, f_temp * f_VPD * f_SW)
  end function gsto_mean
  
  
  !> Calculate the VPD for a given f_VPD based on the linear VPD model.
  pure real function inverse_f_VPD_linear(f_VPD, VPD_max, VPD_min, fmin)
    real, intent(in) :: f_VPD
    real, intent(in) :: VPD_max   !< VPD for maximum gsto (kPa)
    real, intent(in) :: VPD_min   !< VPD for minimum gsto (kPa)
    real, intent(in) :: fmin      !< Minimum f_VPD

    real :: slope

    slope = (1.0 - fmin) / (VPD_min - VPD_max)
    inverse_f_VPD_linear = VPD_min - ((f_VPD - fmin) / slope)
  end function inverse_f_VPD_linear
 
  pure real function inverse_f_VPD_log(f_VPD, fmin)
    real, intent(in) :: f_VPD
    real, intent(in) :: fmin      !< Minimum f_VPD

    inverse_f_VPD_log = exp((1 - f_VPD) / 0.6)
  end function inverse_f_VPD_log
  
END MODULE
