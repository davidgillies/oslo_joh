!
!////////////////////////////////////////////////////////////////////////
!
!      main.f90
!      Created: 12 August 2016 13:44 
!      By: David Gillies  
!
!////////////////////////////////////////////////////////////////////////
!
PROGRAM do3se
USE DO3SE_Gsto_params
USE DO3SE_Phenology
USE DO3SE_Met
USE DO3SE_Micromet
USE DO3SE_Gsto
USE DO3SE_Ozonedep
USE DO3SE_Ozonedose
USE DO3SE_Soilmoisture
USE DO3SE_Soilmoisturechanges

IMPLICIT NONE

!  INPUTS FOR PHENOLOGY
   CHARACTER(len=16) :: height_method
   REAL :: height
   CHARACTER(len=16) :: LAI_method
   INTEGER :: SGS
   REAL :: LAI_a
   REAL :: LAI_b
   INTEGER :: LAI_1
   REAL :: LAI_c
   INTEGER :: EGS
   REAL :: LAI_d
   INTEGER :: LAI_2
   INTEGER :: dd
   CHARACTER(len=16) :: SAI_method
   CHARACTER(len=16) :: season_SAI_method
   REAL :: canopy_height
   REAL  :: LAI
   REAL :: SAI
! END PHENOLOGY INPUTS

! INPUTS FOR MET 
   REAL :: Ts_C
   REAL :: VPD
   REAL :: RH
   REAL :: esat
   REAL :: eact
    REAL :: P        !< Atmospheric pressure (kPa)
    ! integer, intent(in) :: dd   DEFINED ABOVE
    integer :: hr   !< Hour of day (0--23)
    REAL :: lat
    REAL :: lon
    REAL :: elev
    REAL :: albedo
    REAL :: Rn
    REAL :: R
    REAL :: PAR 
    REAL :: PPFD
    REAL :: Idrctt
    REAL :: Idfuse
    REAL :: sinB
   CHARACTER(len=16) :: CO2_method
   REAL :: CO2_constant
   REAL :: CO2
   REAL :: VPD_dd

! END INPUTS FOR MET

! INPUTS FOR MICROMET

    !real :: Idrctt        !< Direct PAR irradiance (W m-2)
    !real :: Idfuse        !< Diffuse PAR irradiance (W m-2)
    !real :: sinB          !< sin() of solar elevation angle
    real :: cosA          !< cos(A), A = mean leaf inclination (0.5 = 60 degrees)
    !real :: LAI           !< Leaf area index (m2 m-2)

    real :: PARsun       !< PAR received by sunlit leaves (W m-2)
    real :: PARshade     !< PAR received by shaded leaves (W m-2)
    LOGICAL :: OTC
    REAL :: h_u
    REAL :: z_u
    REAL :: u_met
    REAL :: u_umet
    REAL :: u_50
    REAL :: ustar
    !real :: canopy_height  
    REAL :: LAIsunfrac

! END INPUTS FOR MICROMET

! INPUTS FOR GSTO_PARAMS

                         
                              


    CHARACTER(len=16) :: f_phen_method
    !integer, intent(in) :: SGS            !< Start of growing season (day of year)
    !integer, intent(in) :: EGS            !< End of growing season (day of year)
    !integer, intent(in) :: dd             !< Day of year
    REAL :: f_phen_a
    REAL :: f_phen_b
    REAL :: f_phen_c
    REAL :: f_phen_d
    REAL :: f_phen_e
    INTEGER :: f_phen_1
    INTEGER :: f_phen_2
    INTEGER :: f_phen_3
    INTEGER :: f_phen_4
    INTEGER :: f_phen_limA
    INTEGER :: f_phen_limB
    CHARACTER(len=16) :: leaf_f_phen_method
    REAL :: leaf_f_phen_a       !< f_phen at Astart
    REAL :: leaf_f_phen_b       !< f_phen at mid-season peak
    REAL :: leaf_f_phen_c       !< f_phen at Aend
    INTEGER :: leaf_f_phen_1    !< Time from _a to _b (days)
    INTEGER :: leaf_f_phen_2    !< Time from _b to _c (days)
    integer :: Astart         !< Start of accumulation period (day of year)
    integer :: Aend           !< End of accumulation period (day of year)
    CHARACTER(len=16) :: f_light_method
    real :: f_lightfac        !< Single leaf f_light coefficient

    !real :: LAIsunfrac        !< Fraction of canopy component that is sunlit
    !real :: PAR               !< Photosynthetically active radiation (W m-2)    
    CHARACTER(len=16) :: f_temp_method
    !real, intent(in) :: Ts_C    !< Air temperature (degrees C)
    real :: T_min   !< Minimum temperature (degrees C)
    real :: T_opt   !< Optimum temperature (degrees C)
    real :: T_max   !< Maximum temperature (degrees C)
    real :: fmin    !< Minimum f_temp
    CHARACTER(len=16) :: f_VPD_method
    !real :: VPD       !< Vapour pressure deficit (kPa)
    real :: VPD_max   !< VPD for maximum gsto (kPa)
    real :: VPD_min   !< VPD for minimum gsto (kPa)
    CHARACTER(len=16) :: f_SW_method
    real :: fSWP_exp_a
    real :: fSWP_exp_b
    real :: SWP       !< Soil water potential (MPa)
    real :: SWP_min   !< SWP for minimum gsto (MPa)
    real :: SWP_max   !< SWP for maximum gsto (MPa)
    real :: ASW_FC    !< Available soil water at field capacity (m)
    real :: ASW       !< Available soil water (m)
    CHARACTER(len=16) :: f_O3_method
    REAL :: POD_0
    REAL :: AOT_0
    CHARACTER(len=16) :: O3_method
    CHARACTER(len=16) :: Ozone_dep_O3_method
    REAL :: FO3_eff
    REAL :: V_cmax_25
    REAL :: J_max_25
    CHARACTER(len=16) :: phenology_method
    REAL :: f_phen
    REAL :: leaf_f_phen
    REAL :: f_light
    REAL :: f_VPD
    REAL :: f_SW
    REAL :: f_O3
    REAL :: leaf_f_light
    REAL :: f_temp              

! END INPUTS FOR GSTO_PARAMS

! INPUTS FOR GSTO

    CHARACTER(len=16) :: method
    REAL :: gmax
    REAL :: gmorph
    !REAL :: fmin
    REAL :: VPD_crit
    !REAL :: VPD_dd
    REAL :: leaf_gsto
    REAL :: mean_gsto
    REAL :: bulk_gsto
    !REAL :: leaf_f_phen
    !REAL :: f_O3
    !REAL :: leaf_f_light
    !REAL :: f_temp
    !REAL :: f_VPD
    !REAL :: f_SW
    !REAL :: f_phen
    !REAL :: f_light
    !REAL :: V_cmax_25
    !REAL :: J_max_25
    !REAL :: eact
    !REAL :: P
    !REAL :: albedo
    !REAL :: R
    !REAL :: Rn
    !REAL :: PPFD
    !REAL :: CO2
    REAL :: u
    REAL :: Lm
    REAL :: Y
    !REAL :: LAI
    CHARACTER(len=16) :: D_0_method
    !CHARACTER(len=16) :: f_VPD_method
    !REAL :: VPD_max
    !REAL :: VPD_min
    CHARACTER(len=16) :: Tleaf_method
    REAL :: g_sto_0
    REAL :: m
    real :: A_n        !< Output: Net CO2 assimilation (umol m-2 PLA s-1)
    real :: g_sto      !< Output: Stomatal conductance (mmol m-2 PLA s-1)
    real :: g_sv       !< Output: Stomatal conductance to water vapour
    real :: g_bv       !< Output: Boundary conductance to water vapour
    real :: Tleaf_C  !< Leaf temperature (degrees C)
    !real :: Tair_C
    !REAL :: Ts_C
    real :: Tleaf_balance_threshold
    real :: Tleaf_adjustment_factor
    INTEGER :: Tleaf_max_iterations


! END INPUTS FOR GSTO

! INPUTS FOR OZONE DEPOSITION

    !REAL :: ustar
    !REAL  :: canopy_height
    !REAL :: SAI
    !REAL :: bulk_gsto
    real :: Ra_c
    real  :: Ra
    real  :: Rb
    real  :: Rinc
    real  :: Rext
    real  :: Rsto
    real  :: Rgs
    real  :: Rsur
    REAL  :: Rtotal
    REAL:: O3_constant
    REAL :: O3_met
    REAL :: O3_50_met
    REAL :: O3_umet
    REAL :: h_O3_loc
    REAL  :: z_O3_loc
    !logical :: OTC
    !REAL :: LAI
    !CHARACTER(len=16) :: O3_method

! END INPUTS FOR OZONE DEPOSITION

! INPUTS FOR OZONE DOSE

   !REAL :: Lm
   !REAL :: u
   !REAL :: Ts_C
   !REAL :: P
   !REAL :: leaf_gsto
   !REAL :: O3_met
   !REAL :: POD_0
   REAL :: POD_Y
   !REAL :: AOT_0
   REAL :: AOT_40
   !CHARACTER(len=16) :: O3_method
   !REAL :: FO3_eff
   REAL :: F_0
! END INPUTS FOR OZONE DOSE

! INPUTS FOR SOIL MOISTURE

    CHARACTER(len=16) :: source
    real :: SWP_AE
    real :: b
    REAL :: PWP
    REAL :: FC
    REAL :: root
    !REAL :: ASW
    REAL :: SMD
    !REAL :: SWP
    REAL :: Sn
    !INTEGER :: hr

! END INPUTS FOR SOIL MOISTURE

! INPUTS FOR SOIL MOISTURE CHANGES

    !CHARACTER(len=16)  :: source
    !CHARACTER(len=16)  :: f_SW_method
    !INTEGER  :: hr
    REAL  :: Rtotal_O3
    real :: Rb_H2O
    real :: Rsto_H2O
    real :: Rinc_H2O
    real :: Rsoil_H2O
    real :: Ra_H2O
    real :: Ei
    real :: Et
    real :: Es
    real :: Eat
    real :: Ei_acc
    real :: Et_acc
    real :: Es_acc
    real :: Eat_acc
    real :: precip_acc
    REAL :: input
    REAL :: run_off
    REAL :: run_off_acc
    REAL :: effective
    REAL :: intercepted_evaporated
    REAL :: evapotranspiration
    REAL :: Sn_diff
    REAL :: percolated
    !real :: LAI     !< Leaf area index (m2 m-2)
    !real :: root    !< Root depth (m)
    real :: run_off_fraction    !< Amount of precipitation that is lost as run-off
    !real  :: ASW     !< Current available soil water, above PWP (m)
    !real  :: SMD     !< Current soil moisture deficit, from FC (m)

! END INPUTS FOR SOIL MOISTURE CHANGES

   integer :: ios
   REAL :: precip



!!!!!!! CONFIG
height_method = 'constant' ! CONFIG
height = 1.0 ! CONFIG
LAI_method = "estimate total" ! CONFIG
SGS = 95 ! CONFIG
LAI_a = 0.0 ! CONFIG - used in PLF to get LAI.
LAI_b = 3.5 ! CONFIG
LAI_1 = 70 ! CONFIG
LAI_c = 3.5 ! CONFIG
EGS = 208 ! CONFIG
LAI_d = 0.0 ! CONFIG
LAI_2 = 22 ! CONFIG
dd = 135 ! CONFIG
SAI_method = "estimate total" ! CONFIG
season_SAI_method = "wheat" ! CONFIG
RH = -999.0 ! CONFIG I think
lat = 53.227 ! CONFIG
lon = -4.129 ! CONFIG
elev = 10 ! CONFIG
albedo = 0.2 ! CONFIG
CO2_method = 'constant' ! CONFIG
CO2_constant = 391.0 ! CONFIG
CO2 = -999.0 ! INPUT or calculated?
cosA = 0.5 ! set in configtypes, can overrule in config
OTC = .TRUE. ! CONFIG
h_u = 1.0 ! CONFIG
z_u = 1.0 ! CONFIG
f_phen_method = 'simple day PLF' ! CONFIG
f_phen_a = 0.2 ! CONFIG
f_phen_b = 1.0 ! CONFIG
f_phen_c = 1.0 ! CONFIG
f_phen_d = 1.0 ! CONFIG
f_phen_e = 0.2 ! CONFIG
f_phen_1 = 25 ! CONFIG
f_phen_2 = -999 ! CONFIG
f_phen_3 = -999 ! CONFIG
f_phen_4 = 34 ! CONFIG
f_phen_limA = -999 ! CONFIG
f_phen_limB = -999 ! CONFIG
leaf_f_phen_method = 'day PLF' ! CONFIG
leaf_f_phen_a = 1.0 ! CONFIG
leaf_f_phen_b = 0.7 ! CONFIG
leaf_f_phen_c = 0 ! CONFIG
leaf_f_phen_1 = 84 ! CONFIG
leaf_f_phen_2 = 9 ! CONFIG
Astart = 150 ! CONFIG
Aend = 208 ! CONFIG
f_light_method = "enabled" ! CONFIG
f_lightfac = 0.0105 ! CONFIG
!PARsun = 0.3 ! where did these come from?
!PARshade = 0.7
f_temp_method = "default" ! CONFIG
T_min = 12 ! CONFIG
T_opt = 26 ! CONFIG
T_max = 40 ! CONFIG
fmin = 0.01 ! CONFIG
f_VPD_method = 'disabled' ! CONFIG
VPD_max = 1.2 ! CONFIG
VPD_min = 3.2 ! CONFIG
f_SW_method = 'fSWP linear' ! CONFIG
fSWP_exp_a = -999.0 ! CONFIG only used for f_SW_method with exp.
fSWP_exp_b = -999.0 ! CONFIG
SWP_min = -0.56 ! CONFIG
SWP_max = -0.049 ! CONFIG
ASW_FC = -999.0 ! CONFIG
f_O3_method = 'disabled' ! CONFIG
Ozone_dep_O3_method = 'input' ! find out which method used where!!
O3_method = 'disabled' ! CONFIG
V_cmax_25 = 60 ! CONFIG OR INPUT
J_max_25 = 30 ! CONFIG OR INPUT
phenology_method = 'disabled' ! CONFIG
method = 'multiplicative' ! CONFIG
gmax = 383.0 ! CONFIG
gmorph = 1.0 ! CONFIG
VPD_crit = 8.35 ! CONFIG
Lm = 0.02 ! CONFIG
Y = 6.0 ! CONFIG
D_0_method = 'f_VPD' ! CONFIG
Tleaf_method = 'Nikolov' ! CONFIG
g_sto_0 = 30000 ! CONFIG
m = 16.83 ! CONFIG
Tleaf_balance_threshold = 0.0010000000474974513 ! CONFIG
Tleaf_adjustment_factor = 0.019999999552965164 ! CONFIG
Tleaf_max_iterations = 50 ! CONFIG
O3_constant = -999.0  ! CONFIG
O3_met = 38.5975 ! CONFIG
h_O3_loc = 1.0 ! CONFIG
z_O3_loc = 1.0 ! CONFIG
F_0 = 1.0 ! CONFIG
source = 'P-M' ! CONFIG
SWP_AE = -999.0 ! CONFIG
b = -999.0 ! CONFIG
PWP = -4.0 ! CONFIG
FC = -999.0 ! CONFIG
root = 0.5 ! CONFIG


!!!!!!! INPUT
Ts_C = 12.393576 ! INPUT
VPD = 0.037988 ! INPUT
hr = 0 ! INPUT
P = 101.13 ! INPUT
PAR = 1508.191906 ! INPUT or set in calc_met_data
PPFD = -999.0 ! INPUT or set in calc_met_data
Idrctt = -999.0 ! INPUT or set in calc_met_data
Idfuse = -999.0 ! INPUT or set in calc_met_data
R = -999.0 ! INPUT or set in calc_met_data
Rn = -999.0 ! INPUT or set in calc_met_data
u_met = 0.75 ! INPUT
u = 0.75
SWP = -999 ! INPUT
precip = 0.01


!!!!!! CALCULATED
LAIsunfrac = -999.0 ! calculated in calc_micromet_data
ASW = -999.0 
POD_0 = 0.0
AOT_0 = -999.0
FO3_eff = 0.0
leaf_gsto = -999.0
mean_gsto = -999.0
bulk_gsto = -999.0
Tleaf_C = Ts_C ! INPUT OR CALCULATED
O3_50_met = -999.0 ! CALCULATED??
O3_umet = -999.0 ! CALCULATED??
POD_Y = 0.0
AOT_40 = 0.0
SMD = -999.0 
Sn = -999.0 
Rtotal_O3 = Rtotal ! Is this right?
Ei = 5.0
Et = 3.0
Es = 5.0
Eat = 6.0
input = 34.0
run_off = 12.0
run_off_acc = 400.0
effective = 12.0
intercepted_evaporated = 1.0
evapotranspiration = 6.0
Sn_diff = 4.0
percolated = 3.0
run_off_fraction = 3.0
! ACCUMULATORS
Ei_acc  = 0.0
Et_acc = 0.0
Es_acc = 0.0
Eat_acc = 0.0
precip_acc = 0.0


DO

read (*, *, iostat=ios) &
    dd, &
    hr, &
    Ts_C, &
    P, &
    precip, &
    u, &
    O3_met, &
    R, &
    PAR, &
    PPFD, &
    VPD


    if (ios < 0) then
      exit
    end if
   

CALL calc_phenology(height_method, height, LAI_method, SGS, LAI_a, LAI_b, LAI_1, &
                          LAI_c, EGS, LAI_d, LAI_2, dd, SAI_method, &
                           season_SAI_method, canopy_height, LAI, SAI)
! END CALL PHENOLOGY

! CALL MET

CALL calc_met_data(Ts_C, P, dd, hr, lat, lon, elev, albedo, &
                          CO2_method, CO2_constant, CO2, & ! ins
                          Rn, R, PAR, PPFD, VPD, RH, Idrctt, Idfuse, & ! inouts
                          esat, eact, sinB, VPD_dd)

CALL calc_micromet_data(Idrctt, Idfuse, sinB, &
                       cosA, LAI, &
                       PARsun, PARshade, &
                       OTC, h_u, z_u, u_met, u_umet, u_50, ustar, canopy_height, &
                       LAIsunfrac)

! END CALL MICROMET

! CALL GSTO_PARAMS

CALL calc_gsto_parameters(f_phen_method, SGS, EGS, dd, f_phen_a, f_phen_b, &
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

! END CALL GSTO_PARAMS

! CALL GSTO

CALL calc_gsto(method, VPD_crit, VPD_dd, leaf_gsto, gmax, &
                                   leaf_f_phen, f_O3, leaf_f_light, fmin, &
                                   f_temp, f_VPD, f_SW, mean_gsto, gmorph, &
                                   f_phen, f_light, VPD_max, VPD_min, &
                                   Tleaf_method, g_sto_0, m, V_cmax_25, J_max_25, &
                                   Lm, Ts_C, u, CO2, &
                                   PPFD, Rn, R, albedo, P, eact, &
                                   Tleaf_C, LAI, bulk_gsto, D_0_method, f_VPD_method, &
                                   Ts_C, Tleaf_balance_threshold, &
                                     Tleaf_adjustment_factor, &
                                     Tleaf_max_iterations, A_n, g_sv, g_bv)

! END CALL GSTO

! CALL OZONE DEPOSITION

CALL calc_ozone_deposition(ustar, canopy_height, SAI, bulk_gsto, O3_constant, &
                                   O3_met, O3_50_met, O3_umet, h_O3_loc, z_O3_loc, OTC, &
                                   Ra_c, Ra, Rb, Rinc, Rext, Rsto, Rgs, Rsur, Rtotal, LAI, &
                                   Ozone_dep_O3_method)


! END CALL OZONE DEPOSITION

! CALL OZONE DOSE

CALL calc_ozone_dose(Lm, u, Ts_C, P, leaf_gsto, O3_met, POD_0, POD_Y, &
                             AOT_0, AOT_40, O3_method, FO3_eff, F_0, Y)

! END CALL OZONE DOSE

! CALL SOIL MOISTURE

CALL calc_soil_moisture(source, SWP, SWP_AE, b, ASW, SMD, PWP, FC, root, &
                                Sn, hr)

! END CALL SOIL MOISTURE

! CALL SOIL MOISTURE CHANGES

CALL calc_soil_moisture_changes(source, hr, Rtotal_O3, Rb_H2O, Rsto_H2O, Rinc_H2O, &
                                        Rsoil_H2O, Ra_H2O, Ei, Et, Es, Eat, &
                                  Ei_acc, Et_acc, Es_acc, Eat_acc, & 
                                   input, run_off, run_off_acc, effective, &
                                   intercepted_evaporated, evapotranspiration, Sn_diff, percolated, &
                                  LAI, root, run_off_fraction, ASW, SMD, precip_acc, &
                                   f_SW_method)


write (*, *) &
    dd, hr, Ts_C, &
    Tleaf_C, VPD, &
    u, precip, &
    P, O3_met, &
    CO2, R, &
    PAR, Rn, &
    LAI, SAI, &
    ustar, O3_50_met, &
    sinB, PPFD, &
    Idrctt, Idfuse, &
    esat, eact, &
    RH, PARsun, &
    PARshade, u, &
    O3_met, Tleaf_C, &
    fmin, &
    gmax, &
    gmorph, &
    f_phen, &
    leaf_f_phen, &
    f_light, &
    leaf_f_light, &
    f_temp, &
    f_VPD, &
    f_SW, &
    f_O3, &
    leaf_gsto, &
    mean_gsto, bulk_gsto, &
    Ra_c, & ! nL,
    Ra, Rb, &
    Rinc, Rext, &
    Rext, Rsto, &
    Rgs, Rsur, &
    Rtotal, &
    !Vd,  & ! leaf_Rb,
    !leaf_Rext, &
    !leaf_Rsto, &
    !Fst, &
    POD_0, &
    POD_y, &
    !OT_0, &
    !OT_40, &
    AOT_0, &
    AOT_40, &
    Ra_c, & ! nL,
    Ra, Rb, &
    Rinc, Rext, &
    Rsto, Rgs, &
    Rsur, Rtotal, &
    Sn_diff, input, &
    run_off, effective, &
    intercepted_evaporated, &
    evapotranspiration, &
    percolated, &
    Sn, SWP, &
    ASW, SMD, &
    A_n, V_cmax_25, &
    g_sv, g_bv

END DO


WRITE(*,*) "done"


END PROGRAM
