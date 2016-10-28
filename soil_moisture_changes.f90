!
!////////////////////////////////////////////////////////////////////////
!
!      soil_moisture_changes.f90
!      Created: 16 August 2016 11:10 
!      By: David Gillies  
!
!////////////////////////////////////////////////////////////////////////
!
MODULE DO3SE_Soilmoisturechanges
USE DO3SE_Ozonedep
USE DO3SE_Met

  real, parameter :: DIFF_H2O = 0.000025
  real, parameter :: DRATIO = 0.663
  real, parameter :: T0 = 273.15

CONTAINS

  !> Calculate inputs to soil moisture, if necessary.  Should probably be called
  !! as the very last thing in the model loop.  Currently only implements the
  !! P-M water vapour stuff.
  subroutine calc_soil_moisture_changes(source, hr, Rtotal_O3, Rb_H2o, Rsto_H2o, Rinc_H2o, &
                                        Rsoil_H2o, Ra_H2o, Ei, Et, Es, Eat, &
                                  Ei_acc, Et_acc, Es_acc, Eat_acc, & 
                                   input, run_off, run_off_acc, effective, &
                                   intercepted_evaporated, evapotranspiration, Sn_diff, percolated, &
                                   LAI, root, run_off_fraction, ASW, SMD, precip_acc, &
                                   f_SW_method)
                                   
    !class(DO3SE_State_t), intent(inout) :: this
    CHARACTER(len=16), INTENT(in) :: source
    CHARACTER(len=16), INTENT(in) :: f_SW_method
    INTEGER, INTENT(in) :: hr
    REAL,INTENT(inout)  :: Rtotal_O3
    LOGICAL :: Es_blocked
    real, intent(out) :: Rb_H2o
    real, intent(out) :: Rsto_H2o
    real, intent(out) :: Rinc_H2o
    real, intent(out) :: Rsoil_H2o
    real, intent(out) :: Ra_H2o
    real, intent(inout) :: Ei
    real, intent(inout) :: Et
    real, intent(inout) :: Es
    real, intent(inout) :: Eat
    real, intent(inout) :: Ei_acc
    real, intent(inout) :: Et_acc
    real, intent(inout) :: Es_acc
    real, intent(inout) :: Eat_acc
    real, intent(inout) :: precip_acc
    REAL, INTENT(inout) :: input
    REAL, INTENT(inout) :: run_off
    REAL, INTENT(inout) :: run_off_acc
    REAL, INTENT(inout) :: effective
    REAL, INTENT(inout) :: intercepted_evaporated
    REAL, INTENT(inout) :: evapotranspiration
    REAL, INTENT(inout) :: Sn_diff
    REAL, INTENT(inout) :: percolated
    real, intent(inout) :: LAI     !< Leaf area index (m2 m-2)
    real, intent(inout) :: root    !< Root depth (m)
    real, intent(inout) :: run_off_fraction    !< Amount of precipitation that is lost as run-off
    real, intent(inout) :: ASW     !< Current available soil water, above PWP (m)
    real, intent(inout) :: SMD     !< Current soil moisture deficit, from FC (m)

    
    select case (source)
    case ("P-M")
      !
      ! Do hourly Penman-Monteith
      !

      ! Reset values at beginning of day
      if (start_of_day(hr)) then
            Ei_acc = 0.0
            Et_acc = 0.0
            Es_acc = 0.0
            Eat_acc = 0.0
            precip_acc = 0.0
      end if

      ! Adapt multi-layer O3 resistance model to single-layer H2O resistance
      !call init_ResistanceModel(this%V%rmodel_H2O, 1)
      Ra_H2O = Ra_O3
      Rb_H2O = Rb_func(ustar, DIFF_H2O)
      ! Combine multi-layer into single-layer
      Rinc_H2O = 1.0/(1.0/Rinc_O3)
      Rext_H2O = 1.0/(1.0/Rext_O3)
      Rsto_H2O = DRATIO * 1.0/(1.0/Rsto_O3)
      Rgs_H2O = Rgs_O3
      Rsur_H2O = Rsur_func(Rb_H2O, Rsto_H2O, Rext_H2O,LAI,SAI)
      Rtotal_O3 = Rtotal_func(Rsur_O3, Rinc_O3, Rgs_O3)
      ! Is soil evaporation blocked?
      ! TODO: this assumes the first land cover is the only one that matters
      !associate (gc => this%LCs(1)%gsto)
        select case (f_SW_method)
        case ("fSWP exp", "fLWP exp")
          Es_blocked = SWP < inverse_f_SWP_exp(fSWP_exp_a, fSWP_exp_b, 1.0)
        case ("fSWP linear")
          Es_blocked = SWP < SWP_max
        case ("fPAW")
          Es_blocked = ASW < inverse_f_PAW(ASW_FC, fmin, 1.0)
        case default
          Es_blocked = .true.
        end select
      !end associate

      ! Run hourly accumulations of precipitation and evapotranspiration

      precip_acc = precip_acc + precip/1000

      call penman_monteith_hourly(Rn*1000000, P*1000, Ts_C, &
                                  esat*1000, eact*1000, VPD*1000, Ra_H2O, &
                                  Rb_H2O, Rsto_H2O, Rsoil_H2O, Rinc_H2O, &
                                  LAI, logical(Es_blocked), Ei, Et, Es, Eat, &
                                  Ei_acc, Et_acc, Es_acc, Eat_acc)

      ! Calculate daily balance
      if (end_of_day(hr)) then
        call penman_monteith_daily(LAI, root, Ei_acc, Eat_acc, &
                                   run_off_fraction, ASW, SMD, & 
                                   recip_acc, input, run_off, run_off_acc, effective, &
                                   intercepted_evaporated, evapotranspiration, Sn_diff, percolated)
      end if
    end select
  end subroutine calc_soil_moisture_changes


  !> Calculate the SWP (MPa) for a given f_SWP based on the exponential
  !! relationship.
  !!
  !! Useful for determining SWP_max when using exponential f_SWP instead of
  !! linear f_SWP.
  pure real function inverse_f_SWP_exp(a, b, f_SWP)
    real, intent(in) :: a
    real, intent(in) :: b
    real, intent(in) :: f_SWP     !< Target f_SWP

    inverse_f_SWP_exp = -(f_SWP/a)**(1.0/b)
  end function inverse_f_SWP_exp
  
  !> Calculate the ASW (m) for a given f_PAW based on the f_PAW relationship.
  pure real function inverse_f_PAW(ASW_FC, fmin, f_PAW)
    real, intent(in) :: ASW_FC    !< Available soil water at field capacity (m)
    real, intent(in) :: fmin      !< Minimum f_PAW
    real, intent(in) :: f_PAW     !< Target f_PAW

    inverse_f_PAW = (ASW_MIN + ((f_PAW-fmin)/(1-fmin)) * (ASW_MAX-ASW_MIN)) * (ASW_FC/100)
  end function inverse_f_PAW
  
  subroutine penman_monteith_hourly(Rn, P, Ts_C, esat, eact, VPD, Ra_H2O, &
                                  Rb_H2O, Rsto_H2O, Rsoil_H2O, Rinc_H2O, LAI, Es_blocked, Ei, Et, Es, Eat, &
                                  Ei_acc, Et_acc, Es_acc, Eat_acc)
    real, intent(in) :: Rn              !< Net radiation (J)
    real, intent(in) :: P               !< Atmospheric pressure (Pa)
    real, intent(in) :: Ts_C            !< Air temperature (degrees C)
    real, intent(in) :: esat            !< Saturated vapour pressure (Pa)
    real, intent(in) :: eact            !< Actual vapour pressure (Pa)
    real, intent(in) :: VPD             !< Vapour pressure deficit (Pa)
    !type(ResistanceModel_t), intent(in) :: rm
    real, intent(in) :: LAI             !< Leaf area index (m2 m-2)
    logical, intent(in) :: Es_blocked   !< Is soil evaporation blocked?
    real, intent(in) :: Rb_H2O
    real, intent(in) :: Rsto_H2O
    real, intent(in) :: Rinc_H2O
    real, intent(in) :: Rsoil_H2O
    real, intent(in) :: Ra_H2O
    real, intent(inout) :: Ei
    real, intent(inout) :: Et
    real, intent(inout) :: Es
    real, intent(inout) :: Eat
    real, intent(inout) :: Ei_acc
    real, intent(inout) :: Et_acc
    real, intent(inout) :: Es_acc
    real, intent(inout) :: Eat_acc
    !type(PM_State_t), intent(inout) :: state

    real :: Tvir, delta, lambda, psychro, Pair, Cair, G

    real :: Et_1, Et_2, Ei_3, Et_3
    real :: t, Es_Rn, Es_G, Es_1, Es_2, Es_3
    real :: SW_a, SW_s, SW_c, C_canopy, C_soil

    ! This model (probably) makes some one-layer assumptions, so don't allow 
    ! multi-layer resistance model.
    !ASSERT(rm%nL == 1)


    !associate (Ra => rm%Ra, &
    !           Rb_H2O => rm%Rb, &
    !           Rinc => rm%Rinc(1), &
    !           Rsto_H2O => rm%Rsto(1), &
    !           Rsoil => rm%Rgs, &
    !           Ei => state%Ei, &
    !           Et => state%Et, &
    !           Es => state%Es, &
    !           Eat => state%Eat)

      Tvir = (Ts_c+T0)/(1-(0.378*(eact/P)))
      delta= ((4098*esat)/((Ts_C+237.3)**2))
      lambda = (2501000-(2361*Ts_C))
      psychro = 1628.6 * (P/lambda)
      Pair = (0.003486*(P/Tvir))
      Cair = (0.622*((lambda*psychro)/P))

      G = 0.1 * Rn
   
      Et_1 = (delta * (Rn - G)) / lambda
      Et_2 = 3600 * Pair * Cair * VPD / Rb_H2O / lambda

      Ei_3 = delta + psychro
      Ei = (Et_1 + Et_2) / Ei_3 / 1000

      Et_3 = delta + psychro * (1 + Rsto_H2O / Rb_H2O)
      Et = (Et_1 + Et_2) / Et_3 / 1000

      if (Es_blocked) then
        Es = 0
      else
        t = exp(-0.5 * LAI)
        Es_Rn = Rn * t
        Es_G = 0.1 * Es_Rn
        Es_1 = (delta * (Rn - G)) / lambda
        Es_2 = ((3600 * Pair * Cair * VPD) - (delta * Rinc_H2O * ((Rn - G) - (Es_Rn - Es_G)))) / (Rinc_H2O + Rb_H2O) / lambda
        Es_3 = delta + (psychro * (1.0 + (Rsoil_H2O / (Rb_H2O + Rinc_H2O))))
        Es = (Es_1 + Es_2) / Es_3 / 1000
      end if

      ! Calculate Eat from Et and Es (after Shuttleworth and Wallace, 1985)
      SW_a = (delta + psychro) * Rb_H2O
      SW_s = (delta + psychro) * Rinc_H20 + (psychro * Rsoil_H2O)
      SW_c = psychro * Rsto_H2O  ! Boundary layer resistance = 0
      C_canopy = 1 / (1 + ((SW_c * SW_a) / (SW_s * (SW_c + SW_a))))
      C_soil = 1 / (1 + ((SW_s * SW_a) / (SW_c * (SW_s + SW_a))))
      if (Es <= 0) then
        Eat = Et
      else
        Eat = (C_canopy * Et) + (C_soil * Es)
      end if

    !end associate

    ! Accumulate values

    Ei_acc = Ei_acc + Ei
    Et_acc = Et_acc + Et
    Es_acc = Es_acc + Es
    Eat_acc = Eat_acc + Eat

  end subroutine penman_monteith_hourly

  !> Is it currently the last hour of the day?
  logical function end_of_day(hr)
    !class(DO3SE_State_t), intent(in) :: this
   INTEGER, INTENT(in) :: hr
    end_of_day = hr == 23
  end function end_of_day

  pure subroutine penman_monteith_daily(Ei_acc, Eat_acc, LAI, root, run_off_fraction, ASW, SMD, &
                                        precip_acc, input, run_off, run_off_acc, effective, &
                                        intercepted_evaporated, evapotranspiration, Sn_diff, percolated)
  
    ! type(PM_State_t), intent(inout) :: state
    real, intent(inout) :: LAI     !< Leaf area index (m2 m-2)
    real, intent(in) :: root    !< Root depth (m)
    real, intent(in) :: run_off_fraction    !< Amount of precipitation that is lost as run-off
    real, intent(in) :: ASW     !< Current available soil water, above PWP (m)
    real, intent(in) :: SMD     !< Current soil moisture deficit, from FC (m)
    real, intent(inout) :: Ei_acc
    real, intent(inout) :: Eat_acc
    real, intent(inout) :: precip_acc
    REAL, INTENT(inout) :: input
    REAL, INTENT(inout) :: run_off
    REAL, INTENT(inout) :: run_off_acc
    REAL, INTENT(inout) :: effective
    REAL, INTENT(inout) :: intercepted_evaporated
    REAL, INTENT(inout) :: evapotranspiration
    REAL, INTENT(inout) :: Sn_diff
    REAL, INTENT(inout) :: percolated
    
    real :: max_ET, delta_SM

    ! Start with full amount of precipitation
    input = precip_acc

    ! Estimate loss to run-off
    run_off = run_off_fraction * input
    run_off_acc = run_off_acc + run_off
    ! Calculate "effective irrigation"
    effective = input - run_off

    ! Estimate loss of intercepted precipitation to evaporation.  Intercepted
    ! precipitation is estimated as 0.0001*LAI, which is therefore a limit on
    ! how much can be evaporated.
    intercepted_evaporated = min(effective, 0.0001 * LAI, Ei_acc)

    ! Can't lose water below PWP, so constrain evapotranspiration to ASW by
    ! restricting evapotranspiration.
    max_ET = ASW + effective - intercepted_evaporated
    evapotranspiration = min(max_ET, Eat_acc)

    ! Total balance = input - run_off - evaporated - evapotranspiration
    delta_SM = effective - intercepted_evaporated - evapotranspiration
    ! Converted to volumetric change using root depth.
    Sn_diff = delta_SM / root

    ! Amount that will go to deep percolation = remainder if water balance
    ! refills soil water, i.e. if it is greater than SMD.
    percolated = max(0.0, delta_SM - SMD)
    percolated_acc = percolated_acc + percolated
  end subroutine penman_monteith_daily
  
END MODULE
