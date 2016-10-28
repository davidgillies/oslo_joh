!
!////////////////////////////////////////////////////////////////////////
!
!      ozone_dose.f90
!      Created: 16 August 2016 09:01 
!      By: David Gillies  
!
!////////////////////////////////////////////////////////////////////////
!
MODULE DO3SE_Ozonedose
  real, parameter :: LEAF_G_O3 = 0.105
  real, parameter :: DT = 60*60
CONTAINS

  subroutine calc_ozone_dose(Lm, u, Ts_C, P, leaf_gsto, O3_met, POD_0, POD_Y, &
                             AOT_0, AOT_40, O3_method, FO3_eff, F_0, Y)
    !class(DO3SE_State_t), intent(inout) :: this
   REAL, INTENT(in) :: Lm
   REAL, INTENT(in) :: u
   REAL, INTENT(in) :: Ts_C
   REAL, INTENT(in) :: P
   REAL, INTENT(in) :: leaf_gsto
   REAL, INTENT(in) :: O3_met
   REAL, INTENT(inout) :: POD_0
   REAL, INTENT(inout) :: POD_Y
   REAL, INTENT(inout) :: AOT_0
   REAL, INTENT(inout) :: AOT_40
   CHARACTER(len=16), INTENT(in) :: O3_method
   REAL, INTENT(inout) :: FO3_eff
   REAL, INTENT(in) :: F_0
   REAL, INTENT(in) :: Y

    real :: O3_ppb_to_nmol

    ! TODO: move this to its own subroutine?
    !this%MLMC(:,:)%leaf_rmodel_O3 = LeafResistanceModel_t()
    !do iL = 1, this%nL
    !  do iLC = 1, this%nLC
    !  associate (LC => this%LCs(iLC), ll => this%ML(iL), xx => this%MLMC(iL,iLC))
    Rb_leaf = leaf_rb(leaf_gb(LEAF_G_O3, Lm, u))
    !  end associate
    !  end do
    !end do
    Rext_leaf = Rext_func(1.0)
    Rsto_leaf = Rsto_func(leaf_gsto)

    O3_ppb_to_nmol = O3_ppb_to_nmol_factor(Ts_C, P)
    !do iL = 1, this%nL
    !  do iLC = 1, this%nLC
    !  associate (LC => this%LCs(iLC), xx => this%MLMC(iL,iLC))
        ! TODO: fix reliance on MAX_RSTO...
    if (leaf_gsto > 0) then
          
       Fst = O3_met * O3_ppb_to_nmol * stomatal_flux_rate(Rsto_leaf, Rext_leaf, Rb_leaf)
    else
       Fst = 0
    end if

    POD_0 = POD_0 + ((Fst*DT)/1000000)
    POD_Y = POD_Y + ((max(0.0, Fst - Y)*DT)/1000000)

        ! Default OT0/OT40 to 0
    OT_0 = 0
    OT_40 = 0

        ! Only accumulate OT when global radiation > 50 W m-2
        if (is_daylight(R)) then
          ! Only accumulate OT0 when leaf_fphen > 0
          if (leaf_f_phen > 0) then
            OT_0 = O3_met / 1000
          end if

          ! Only accumulate OT40 when fphen > 0
          if (f_phen > 0) then
            OT_40 = max(0.0, O3_met - 40) / 1000
          end if
        end if

        ! Accumulate OT0/OT40
        AOT_0 = AOT_0 + OT_0
        AOT_40 = AOT_40 + OT_40

        if (O3_method == "martin2000") then
          ! Effective ozone dose for effect on photosynthesis
          FO3_eff = FO3_eff + (Fst - F_0) * DT
          FO3_eff = max(0.0, FO3_eff)
        end if
    !  end associate
    !  end do
    !end do
  end subroutine calc_ozone_dose

  elemental real function leaf_rb(gb)
    real, intent(in) :: gb    !< Leaf boundary layer conductance (mol m-2 s-1)

    ! gb / 41 : 'mol m-2 s-1' to 'm s-1'
    ! 1 / gb  : 'm s-1' to 's m-1'
    leaf_rb = 41 / gb
  end function leaf_rb

  elemental real function leaf_gb(G, Lm, u)
    real, intent(in) :: G     !< Leaf surface conductance (mol m-2 s-1)
    real, intent(in) :: Lm    !< Cross-wind leaf dimension (m)
    real, intent(in) :: u     !< Wind speed (m s-1)

    ! G * 2 : from single surface to PLA (both sides of leaf)
    leaf_gb = (G * 2) * sqrt(u / Lm)
  end function leaf_gb

  !> Estimate external plant cuticle resistance (Rext, s m-1).
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
  
  !> Calculate the conversion factor between parts per billion and nmol m-3 for
  !! O3 at a given temperature and pressure.
  pure real function O3_ppb_to_nmol_factor(Ts_c, P)
    real, intent(in) :: Ts_C      !< Air temperature (degrees C)
    real, intent(in) :: P         !< Atmospheric pressure (kPa)

    real, parameter :: M_O3 = 48.0      ! Molecular weight of O3 (g)

    real :: Vn

    ! Specific molar volume of an ideal gas at this temperature and pressure
    Vn = 8.314510 * ((Ts_C + T0) / P)
    ! Conversion to nmol m-3 (1 microgram O3 = 20.833 nmol m-3)
    O3_ppb_to_nmol_factor = (1.0/Vn) * M_O3 * 20.833
  end function O3_ppb_to_nmol_factor

  elemental real function stomatal_flux_rate(Rsto, Rext, Rb)
    !type(LeafResistanceModel_t), intent(in) :: leaf_rmodel
    REAL, INTENT(IN) :: Rsto
    REAL, INTENT(IN) :: Rext
    REAL, INTENT(IN) :: Rb
    real :: leaf_r

    leaf_r = 1.0 / ((1.0 / Rsto) + (1.0 / Rext))
    stomatal_flux_rate = (1.0/Rsto) * (leaf_r / (Rb + leaf_r))
  end function stomatal_flux_rate
  
  logical function is_daylight(R)
    !class(DO3SE_State_t), intent(in) :: this
    REAL, INTENT(in) :: R
    is_daylight = R > 50.0
  end function is_daylight
END MODULE
