!
!////////////////////////////////////////////////////////////////////////
!
!      photosynthesis.f90
!      Created: 15 August 2016 14:09 
!      By: David Gillies  
!
!////////////////////////////////////////////////////////////////////////
!
MODULE DO3SE_Photosynthesis
USE DO3SE_Met
  real, parameter :: LEAF_G_H2O = 0.147
  real, parameter :: DRATIO = 0.663
  
CONTAINS
 
 



  pure subroutine gsto_pn(Tleaf_method, g_sto_0, m,    V_cmax_25, J_max_25, D_0, Lm, Tair_C, u, CO2, &
                          PPFD, Rn, R, albedo, P, eact, &
                          Tleaf_C, A_n, g_sto, g_sv, g_bv,  Tleaf_balance_threshold, &
                                     Tleaf_adjustment_factor, &
                                     Tleaf_max_iterations)
    CHARACTER(len=16), INTENT(in) :: Tleaf_method
    REAL, INTENT(in) :: g_sto_0
    REAL,INTENT(IN) :: m
    !type(pngstoconfig_t), intent(in) :: pgc   !< Photosynthesis gsto parameters
    real, intent(in) :: V_cmax_25   !< Maximum catalytic rate at 25 degrees (umol m-2 s-1)
    real, intent(in) :: J_max_25    !< Maximum rate of electron transport at 25 degrees (umol m-2 s-1)
    real, intent(in) :: D_0         !< "The VPD at which g_sto is reduced by a factor of 2" (kPa) (Leuning et al. 1998)
    real, intent(in) :: Lm          !< Leaf dimension (m)
    real, intent(in) :: Tair_C      !< Ambient air temperature (degrees C)
    real, intent(in) :: u           !< Wind speed (m/s)
    real, intent(in) :: CO2         !< CO2 concentration (ppm)
    real, intent(in) :: Rn          !< Net radiation (MJ m-2 h-1)
    real, intent(in) :: R           !< Global radiation (W m-2)
    real, intent(in) :: albedo      !< Surface albedo (fraction)
    real, intent(in) :: PPFD        !< PPFD (umol m-2 s-1)
    real, intent(in) :: P           !< Ambient air pressure (kPa)
    real, intent(in) :: eact        !< Ambient vapour pressure (kPa)

    real, intent(inout) :: Tleaf_C  !< Leaf temperature (degrees C)
    real, intent(out) :: A_n        !< Output: Net CO2 assimilation (umol m-2 PLA s-1)
    real, intent(inout) :: g_sto      !< Output: Stomatal conductance (mmol m-2 PLA s-1)
    real, intent(out) :: g_sv       !< Output: Stomatal conductance to water vapour
    real, intent(out) :: g_bv       !< Output: Boundary conductance to water vapour
    real, INTENT(in) :: Tleaf_balance_threshold
    real, INTENT(in) :: Tleaf_adjustment_factor
    INTEGER, INTENT(in) :: Tleaf_max_iterations

    integer, parameter :: MAX_ITERATIONS = 10
    integer, parameter :: MAX_TLEAF_ITERATIONS = 5

    integer :: i, j
    real :: R_ni

    ! aproximates the boundary layer conductance for forced convection
    ! (converted to umol m-2 s-1)
    g_bv = leaf_gb(LEAF_G_H2O, Lm, max(0.01, u)) * 1e6

    select case (Tleaf_method)
    case ("input")
      ! Tleaf_C supplied, just run photosynthesis
      call farquhar_photosynthesis(Tleaf_C, CO2, eact*1e3, PPFD, g_bv, &
                                   g_sto_0, m, V_cmax_25, J_max_25, D_0*1e3, &
                                   g_sv, A_n)
    case ("ambient")
      ! Use ambient air temperature as leaf temperature
      Tleaf_C = Tair_C
      call farquhar_photosynthesis(Tleaf_C, CO2, eact*1e3, PPFD, g_bv, &
                                   g_sto_0, m, V_cmax_25, J_max_25, D_0*1e3, &
                                   g_sv, A_n)
    case ("Nikolov")
      ! Estimate Tleaf_C according to Nikolov (1995), iteratively solved with
      ! A_n and g_sto, starting with Tleaf_C = Tair_C.
      Tleaf_C = Tair_C
      do i = 1, MAX_ITERATIONS
        call farquhar_photosynthesis(Tleaf_C, CO2, eact*1e3, PPFD, g_bv, &
                                     g_sto_0, m, V_cmax_25, J_max_25, D_0*1e3, &
                                     g_sv, A_n)
        ! TODO: leaf wetness status
        Tleaf_C = nikolov_leaf_temperature(Tair_C, P*1e3, Rn*Rn_MJ_to_W, eact*1e3, &
                                           g_sv, g_bv, 0.0)
      end do
    case ("EB")
      ! Estimate Tleaf_C according to Campbell & Norman (1998), iteratively
      ! solved with A_n and g_sto, starting with Tleaf_C = Tair_C.
      Tleaf_C = Tair_C
      call farquhar_photosynthesis(Tleaf_C, CO2, eact*1e3, PPFD, g_bv, &
                                   g_sto_0, m, V_cmax_25, J_max_25, D_0*1e3, &
                                   g_sv, A_n)
      do i = 1, MAX_ITERATIONS
        do j = 1, MAX_TLEAF_ITERATIONS
          R_ni = leaf_net_radiation(R, eact, Tair_C, Tleaf_C, albedo, 1.0)
          Tleaf_C = leaf_temp_EB(Lm, Tair_C, P, eact, R_ni, u, g_sv*1e-6)
        end do
        call farquhar_photosynthesis(Tleaf_C, CO2, eact*1e3, PPFD, g_bv, &
                                     g_sto_0, m, V_cmax_25, J_max_25, D_0*1e3, &
                                     g_sv, A_n)
      end do
    case ("de Boeck")
      ! Estimate Tleaf_C according to de Boeck (2012), iteratively
      ! solved with A_n and g_sto, starting with Tleaf_C = Tair_C.
      Tleaf_C = Tair_C
      call farquhar_photosynthesis(Tleaf_C, CO2, eact*1e3, PPFD, g_bv, &
                                   g_sto_0, m, V_cmax_25, J_max_25, D_0*1e3, &
                                   g_sv, A_n)
      do i = 1, MAX_ITERATIONS
        ! TODO: real "hypostomatous" setting, real "cloud cover" value?
        Tleaf_C = leaf_temp_de_Boeck(R, eact*1e3, Tair_C, Tleaf_C, P*1e3, &
                                     u, g_sv*1e-6, .true., Lm, albedo, 1.0, &
                                     Tleaf_balance_threshold, &
                                     Tleaf_adjustment_factor, &
                                     Tleaf_max_iterations)
        call farquhar_photosynthesis(Tleaf_C, CO2, eact*1e3, PPFD, g_bv, &
                                     g_sto_0, m, V_cmax_25, J_max_25, D_0*1e3, &
                                     g_sv, A_n)
      end do
    end select

    ! Convert g_sto from umol to mmol, and from H2O to O3
    g_sto = DRATIO * (max(0.0, g_sv) / 1000)
  end subroutine gsto_pn
  
  
  
  
  elemental real function leaf_gb(G, Lm, u)
    real, intent(in) :: G     !< Leaf surface conductance (mol m-2 s-1)
    real, intent(in) :: Lm    !< Cross-wind leaf dimension (m)
    real, intent(in) :: u     !< Wind speed (m s-1)

    ! G * 2 : from single surface to PLA (both sides of leaf)
    leaf_gb = (G * 2) * sqrt(u / Lm)
  end function leaf_gb
  
  
  pure subroutine farquhar_photosynthesis(Tleaf_C, c_a, e_a, Q, g_bl, &
                                          g_sto_0, m, V_cmax_25, J_max_25, D_0, &
                                          gsto_final, pngsto_An)
    real, intent(in) :: Tleaf_C       !< Leaf temperature (degrees C)
    real, intent(in) :: c_a           !< CO2 concentration (ppm)
    real, intent(in) :: e_a           !< Ambient vapour pressure (Pa)
    real, intent(in) :: Q             !< PPFD (umol/m^2/s)
    real, intent(in) :: g_bl          !< Boundary layer conductance to H2O vapour (micromol m-2 PLA s-1)
    real, intent(in) :: g_sto_0       !< Closed stomata conductance (umol/m^2/s)
    real, intent(in) :: m             !< Species-specific sensitivity to An (dimensionless)
    real, intent(in) :: V_cmax_25     !< Maximum catalytic rate at 25 degrees (umol/m^2/s)
    real, intent(in) :: J_max_25      !< Maximum rate of electron transport at 25 degrees (umol/m^2/s)
    real, intent(in) :: D_0           !< "The VPD at which g_sto is reduced by a factor of 2" (Pa) (Leuning et al. 1998)

    real, intent(out) :: gsto_final   !< Output: Raw photosynthesis-based stomatal conductance (umol m-2 s-1)
    real, intent(out) :: pngsto_An    !< Output: net CO2 assimilation (umol m-2 s-1)

    ! parameters considered (or defined) to be constant for all species
    real, parameter :: O_i = 210.0           !O2 concentration                   [mmol/mol]
    real, parameter :: E_K_C = 79430.0       !activation energy of K_C           [J/mol]            Medlyn2002
    real, parameter :: E_K_O = 36380.0       !activation energy of K_O           [J/mol]            Medlyn2002
    real, parameter :: E_R_d = 53000.0       !activation energy of R_d           [J/mol]            Leuning1995
    real, parameter :: E_Gamma_star = 37830.0 !activation energy for C-comp-point [J/mol]            Medlyn2002
    real, parameter :: K_C_25 = 404.9        !K.C at reference temperature 25    [micro mol/mol]    Medlyn2002
    real, parameter :: K_O_25 = 278.4        !K.O at reference temperature 25    [mmol/mol]         Medlyn2002
    real, parameter :: R_d_20 = 0.32         !R_d at reference temperature 20    [micro mol/(m^2*s)]Leuning1995
    real, parameter :: Gamma_star_25 = 42.75 !CO2 compensation point at T= 25    [micro mol/mol]    Medlyn2002
    real, parameter :: A_j_a = 4.0           !electron requirement for NADPH formation
    real, parameter :: A_j_b = 8.0           !electron requirement for ATP formation

    ! species spedific model parameters (that don't tend to have species specific
    ! values, others are supplied as arguments)
    real, parameter :: alpha = 0.3           !efficiency light energy conversion [mol electrons/mol photons]
    real, parameter :: Teta = 0.90           !shape of J~Q determining factor    []
    real, parameter :: H_a_jmax = 50300      !activation energy for J_max        [J/mol]
    real, parameter :: H_d_jmax = 152044     !deactivation energy for J_max      [J/mol]
    real, parameter :: H_a_vcmax = 73637     !activation energy for V_cmax       [J/mol]
    real, parameter :: H_d_vcmax = 149252    !deactivation energy for V_cmax     [J/mol]
    real, parameter :: S_V_vcmax = 486       !entropy terms                      [J/(mol*K)]
    real, parameter :: S_V_jmax = 495        !entropy terms                      [J/(mol*K)

    ! Converted inputs
    real :: Tleaf_K

    ! state variables
    real :: A_n                             !netto assimilation rate            [micro mol/(m^2*s)]
    real :: A_c                             !Rub. activity. lim. ass. rate      [micro mol/(m^2*s)]
    real :: A_j                             !electr. transp. lim. ass. rate     [micro mol/(m^2*s)]
    real :: A_p                             !triose phosphate utilisation lim. ass. rate   [micro mol/(m^2*s)]
    real :: Gamma_star                      !CO2 comp. point without day resp.  [micro mol/mol]
    real :: R_d                             !day respiration rate               [micro mol/(m^2*s)]
    real :: K_C                             !Michaelis constant CO2             [micro mol/mol]
    real :: K_O                             !Michaelis constant O2              [mmol/mol]
    real :: J                               !Rate of electron transport         [micro mol/(m^2*s)]
    real :: Gamma                           !CO2 compensation point             [micro mol/mol]
    real :: J_max                           !Max rate of electron transport     [micro mol/(m^2*s)]
    real :: V_cmax                          !Max catalytic rate of Rubisco      [micro mol/(m^2*s)]
    real :: e_sat_i                         !internal saturation vapour pressure[Pa]
    real :: g_sto                           !two sided stomatal conduct.,vapour [micro mol/(m^2s)]
    real :: h_s                             !relative humidity at leaf surface  [decimal fraction]
    real :: h_s_VPD                         !VPD at leaf surface                [Pa]
    real :: c_s                             !CO2 concentration at leaf surface  [micromol/mol]
    real :: c_i                             !CO2 concentration inside stomata   [micromol/mol]

    ! iteration parameters
    integer :: iterations                   !number of the iterations bofore convergence
    real :: c_i_sup                         !CO2 concentration inside stomata possible through supply
    integer :: k                            !loop parameters

    Tleaf_K = T0 + Tleaf_C

    ! Calculation of the model variables which are only
    ! dependend on environmental conditions:

    Gamma_star = temp_dep(Gamma_star_25, T0 + 25, E_Gamma_star, Tleaf_K)

    K_C = temp_dep(K_C_25, T0 + 25, E_K_C, Tleaf_K)

    K_O = temp_dep(K_O_25, T0 + 25, E_K_O, Tleaf_K)

    R_d = temp_dep(R_d_20, T0 + 20, E_R_d, Tleaf_K)

    J_max = temp_dep_inhibit(J_max_25, T0 + 25, H_a_jmax, H_d_jmax, S_V_jmax, Tleaf_K)

    V_cmax = temp_dep_inhibit(V_cmax_25, T0 + 25, H_a_vcmax, H_d_vcmax, S_V_vcmax, Tleaf_K)

    ! Electron transport rate
    J = (J_max + alpha*Q - sqrt((J_max + alpha*Q)**2 - 4*alpha*Q*Teta*J_max)) / (2*Teta)

    e_sat_i    = 1000 * saturated_vapour_pressure(Tleaf_C)

    Gamma        = (Gamma_star+(K_C*R_d*(1+(O_i/K_O))/V_cmax))/&
                   (1-(R_d/V_cmax))

    !The following loop guesses a start value for c_i and tests whether
    !it satisfies all the relevant restrictions. If not a new value for
    !c_i is tested:

    c_i         = 0.0

    ! gsto needs a starting point, so let's set it to g_sto_0
    g_sto = g_sto_0

    do k=1,50

      ! Rubisco activity limited assimilation rate
      A_c = V_cmax * ((c_i - Gamma_star) &
                      / (c_i + (K_C * (1 + (O_i / K_O)))))

      ! RuBP regeneration (electron transport) limited assimilation rate
      A_j = J * ((c_i - Gamma_star) &
                 / ((A_j_a * c_i) + (A_j_b * Gamma_star)))

      ! Triose phosphate utilisation limited assimilation rate
      A_p = 0.5 * V_cmax

      ! CO2 assimilation rate
      A_n = min(A_c, A_j, A_p) - R_d

      ! Surface CO2
      c_s = c_a - (A_n * (1.37/g_bl))

      ! Surface humidity
      h_s = (g_sto*e_sat_i + g_bl*e_a) / (e_sat_i * (g_sto + g_bl))
      ! Convert relative humidity to VPD
      h_s_VPD = e_sat_i - (e_sat_i * h_s)

      ! Stomatal conductance
      ! TODO: use humidity deficit version instead
      g_sto = g_sto_0 + m * (A_n / ((1 + (h_s_VPD/D_0)) * (c_s - Gamma)))*1e6

      ! CO2 supply
      c_i_sup = c_a - ((A_n*(1.6/g_sto + 1.37/g_bl))*1e6)

      !exits the loop when c_i calculated with both ways meet the convergence
      !criterium:

      iterations = k
      if (abs(c_i - c_i_sup) < 0.001) then
        exit
      end if

      !Guesses a new c_i as the mean of the first guess and c_i resulting from
      !the supply function:

      c_i      = c_i-(c_i-c_i_sup)/2

    end do

    ! Calculate final stomatal conductances
    gsto_final = g_sto
    pngsto_An = A_n
  end subroutine farquhar_photosynthesis
  
  
  pure function nikolov_leaf_temperature(Tair_C, Pr, R_i, e_a, g_sv, g_bv, Wstat) result(Tleaf_C)
    real, intent(in) :: Tair_C    !< Ambient air temperature (degrees C)
    real, intent(in) :: Pr        !< Air pressure (Pa)
    real, intent(in) :: R_i       !< Net radiation absorbed by leaf (W m-2)
    real, intent(in) :: e_a       !< Water vapour pressure in ambient air (Pa)
    real, intent(in) :: g_sv      !< Leaf stomatal conductance to water vapour (micromol m-2 s-1)
    real, intent(in) :: g_bv      !< Leaf boundary layer conductance to water (micromol m-2 s-1)
    real, intent(in) :: Wstat     !< Leaf wetness status (0 = dry, 1 = wet)

    real :: Tleaf_C  !< Output: leaf temperature (degrees C)

    !> Specific heat capacity of dry air at standard pressure and 20C (J kg-1 K-1)
    real, parameter :: C_P = 1010.0
    !> Triple-point temperature of water (degrees K)
    real, parameter :: T3_H2O = 273.16
    !> Leaf thermal emissivity
    real, parameter :: LTE = 0.975
    !> Stefan-Boltzmann constant (W m-2 K-4)
    real, parameter :: SBC = 5.670373e-8

    ! Constants for saturation vapour pressure curve:
    !   e_s(T) = a*T**4 + b*T**3 + c*T**2 + d*T + e
    real, parameter :: SVP_A = 5.82436e-4
    real, parameter :: SVP_B = 1.5842e-2
    real, parameter :: SVP_C = 1.55186
    real, parameter :: SVP_D = 44.513596
    real, parameter :: SVP_E = 607.919

    ! For deriving the quartic coefficients
    real :: Cfm, lambda, psychro, Tvir, rho, h_e, h_t, k, a, b, c, d
    ! For solving the quartic
    real :: y, E, sa, P, Q, Dscr, R, t1

    ! Conversion from micromol m-2 s-1 to m s-1
    Cfm = 8.3089764 * 1e-6 * ((Tair_C + T0) / Pr)
    ! Latent heat of vapourisation for water (J kg-1)
    lambda = (-0.0000614342*Tair_C**3 + 0.00158927*Tair_C**2 - 2.36418*Tair_C + 2500.79) * 1000
    ! Psychrometric parameter (Pa C-1)
    psychro = (C_P * Pr) / (0.622 * lambda)
    ! Virtual temperature for density calcualation (K)
    Tvir = (Tair_C + T0) / (1 - (0.378 * (e_a / Pr)))
    ! Density of air (kg m-3)
    rho = Pr / (287.058 * Tvir)

    if (Wstat <= 0) then
      h_e = ((rho * C_P) / psychro) * Cfm * g_bv
    else
      h_e = ((rho * C_P) / psychro) * Cfm * ((g_sv * g_bv) / (g_sv + g_bv))
    end if

    ! 0.924 = ratio of heat conductance to vapour conductance
    h_t = rho * C_P * Cfm * (0.924 * g_bv)

    ! Set up coefficients for quartic: T**4 + a*T**3 + b*T**2 + c*T + d = 0
    k = 1 / ((2 * LTE * SBC) + (SVP_A * h_e))
    a = ((8  * LTE * SBC * T3_H2O   ) + (SVP_B * h_e)) * k
    b = ((12 * LTE * SBC * T3_H2O**2) + (SVP_C * h_e)) * k
    c = ((8  * LTE * SBC * T3_H2O**3) + (SVP_D * h_e) + h_t) * k
    d = ((2  * LTE * SBC * T3_H2O**4) + (SVP_E * h_e) - (h_t * Tair_C) - (h_e * e_a) - R_i) * k

    ! Presumably the analytical solution for the quartic?
    y = a * c - 4 * d
    E = b**2
    sa = a**2
    P = (3 * y - E) / 9
    Q = (b * (2 * E - 9 * y) - 27 * (d * (4 * b - sa) - c**2)) / 54
    Dscr = sqrt(Q**2 + P**3)
    y = exp(log(Q + Dscr) / 3) - exp(log(Dscr - Q) / 3) + b / 3
    R = sqrt(0.25 * sa + y - b)
    t1 = 0.25 * (a * (4 * b - sa) - 8 * c) / R
    E = 0.5 * sa - b - y
    Tleaf_C = -0.25 * a - 0.5 * (R - sqrt(E - t1))
  end function nikolov_leaf_temperature
  
  
  pure function leaf_net_radiation(R, e_a, T_air, T_leaf, albedo, cloud_cover) result(Rn_out)
    real, intent(in) :: R             !< Global radiation (W m-2)
    real, intent(in) :: e_a           !< Vapour pressure (kPa)
    real, intent(in) :: T_air         !< Air temperature (degrees C)
    real, intent(in) :: T_leaf        !< Leaf temperature (degrees C)
    real, intent(in) :: albedo        !< Soil albedo (fraction)
    real, intent(in) :: cloud_cover   !< Cloud cover (fraction)

    real :: Rn      !< Output: leaf net (short- and long-wave) radiation (W m-2)

    ! Leaf short-wave absorptivity (de Boeck (2012); from Campbell & Norman (1998), p153, in text)
    real, parameter :: alpha_s_leaf = 0.5
    ! Leaf long-wave absorptivity (de Boeck, 2012)
    real, parameter :: alpha_l_leaf = 0.97
    ! Soil long-wave emissivity (de Boeck, 2012)
    real, parameter :: eps_l_soil = 0.945
    ! Leaf long-wave emissivity (de Boeck, 2012)
    real, parameter :: eps_l_leaf = 0.97

    real :: T_soil, eps_ac, eps_l_sky, R_s_in, R_l_in, R_l_out

    ! TODO: should this be something else?
    T_soil = T_air

    ! Clear sky emissivity (p163, eq. 10.10)
    eps_ac = 1.72 * (e_a / (T_air + T0))**(1.0/7.0)
    ! Sky long-wave emmisivity (p164, eq. 10.12)
    eps_l_sky = (1 - 0.84 * cloud_cover) * eps_ac + 0.84 * cloud_cover

    ! Net radiation calculations from de Boeck (2012):
    ! Absorbed short-wave (direct on top of leaf, reflected from soil on underside)
    R_s_in = (0.5 * R * alpha_s_leaf) + (0.5 * R * albedo * alpha_s_leaf)
    ! Absorbed long-wave (emitted from soil on underside, emitted from sky on top)
    R_l_in = (0.5 * alpha_l_leaf * eps_l_soil * SBC * (T_soil + T0)**4) + &
             (0.5 * alpha_l_leaf * eps_l_sky * SBC * (T_air + T0)**4)
    ! Emitted long-wave
    R_l_out = eps_l_leaf * SBC * (T_leaf + T0)**4
    Rn_out = R_s_in + R_l_in - R_l_out
  end function leaf_net_radiation
  

  pure function leaf_temp_EB(d, T_air, P, e_a, R_ni, u, g_vs) result(T_leaf)
    real, intent(in) :: d       !< Leaf characteristic dimension (m)
    real, intent(in) :: T_air   !< Air temperature (degrees C)
    real, intent(in) :: P       !< Air pressure (kPa)
    real, intent(in) :: e_a     !< Vapour pressure (kPa)
    real, intent(in) :: R_ni    !< Leaf net radiation (W m-2)
    real, intent(in) :: u       !< Wind speed (m s-1)
    real, intent(in) :: g_vs    !< Stomatal conductance to water vapour (mol m-2 s-1)

    real :: T_leaf    !< Output: leaf temperature (degrees C)

    ! Specific heat capacity of dry air (J mol-1 C-1; p279, table A1)
    real, parameter :: c_p = 29.3
    ! Latent heat of vapourisation of water (J mol-1; p37, in text)
    real, parameter :: lam = 44e3
    ! Emissivity of natural surface (fraction; p163, below eq. 10.9)
    real, parameter :: eps_s = 0.97

    real :: T_ak, g_r, g_Ha, g_Hr, g_va, g_v, a, b, c, e_s, VPD, delta, s, gam, gamstar

    T_ak = T_air + T0
    ! Radiative conductance to heat transfer (mol m-2 s-1; p188, eq. 12.7)
    g_r = 4 * eps_s * SBC * T_ak**3 / c_p
    ! Boundary layer (forced convection) conductance to heat transfer (mol m-2 s-1; p101, eq. 7.30)
    ! 1.4 enhancement factor for outdoor "turbulent wind" (p108, section 7.13)
    g_Ha = 1.4 * 0.135 * sqrt(u / d)
    ! Convective-radiative conductance to heat transfer (mol m-2 s-1; p188, eq. 12.9)
    g_Hr = g_Ha + g_r

    ! Boundary layer conductance to water vapour (mol m-2 s-1)
    ! 1.4 enhancement factor for outdoor "turbulent wind" (p108, section 7.13)
    g_va = 1.4 * 0.147 * sqrt(u / d)
    ! Whole-leaf (stomatal and boundary) conductance to water vapour (mol m-2 s-1; p224, eq. 14.2)
    ! (hypostomatus, adaxial stomatal conductance = 0)
    ! TODO: hypostomatus boolean switch
    g_v = leaf_conductance(g_vs, 0.0, g_va)

    ! Saturation vapour pressure coefficients (p41, below eq. 3.8)
    a = 0.611   ! (kPa)
    b = 17.502
    c = 240.97  ! (degrees C)
    ! Saturated vapour pressure (kPa; p41, eq. 3.8)
    e_s = a * exp(b * T_air / (T_air + c))
    ! Vapour deficit (kPa)
    VPD = e_s - e_a
    ! Slope of saturation vapour pressure function (p41, eq. 3.9)
    delta = b * c * e_s / (c + T_air)**2
    ! Slope of saturation mole fraction (p41, eq. 3.10)
    s = delta / P

    ! Thermodynamic psychrometer constant (C-1; p44, below eq. 3.16)
    gam = c_p / lam
    ! Apparent psychrometer constant (C-1; p225, below eq. 14.6)
    gamstar = gam * g_Hr / g_v
    ! Leaf temperature (degrees C; p225, eq. 14.6)
    T_leaf = T_air + (gamstar / (s + gamstar)) * (R_ni / (g_Hr * c_p) - VPD / (P * gamstar))
  contains
    pure real function leaf_conductance(gs_ab, gs_ad, ga)
      real, intent(in) :: gs_ab   !< Abaxial stomatal conductance (mol m-2 s-1)
      real, intent(in) :: gs_ad   !< Adaxial stomatal conductance (mol m-2 s-1)
      real, intent(in) :: ga      !< Boundary layer conductance (mol m-2 s-1)

      leaf_conductance = (0.5 * gs_ab * ga / (gs_ab + ga)) + (0.5 * gs_ad * ga / (gs_ad + ga))
    end function leaf_conductance
  end function leaf_temp_EB
  
  
  pure function leaf_temp_de_Boeck(R, e_a, T_air, initial_T_leaf, P, u, g_vs, hypostomatous, d, albedo, cloud_cover, &
                                   balance_threshold, adjustment_factor, max_iterations) result(T_leaf)
    real, intent(in) :: R                   !< Global radiation (W m-2)
    real, intent(in) :: e_a                 !< Vapour pressure (Pa)
    real, intent(in) :: T_air               !< Air temperature (degrees C)
    real, intent(in) :: initial_T_leaf      !< T_leaf starting value (degrees C)
    real, intent(in) :: P                   !< Air pressure (Pa)
    real, intent(in) :: u                   !< Wind speed (m s-1)
    real, intent(in) :: g_vs                !< Stomatal conductance to water vapour (mol m-2 s-1)
    logical, intent(in) :: hypostomatous    !< Are leaves hypostomatous?
    real, intent(in) :: d                   !< Leaf characteristic dimension (m)
    real, intent(in) :: albedo              !< Surface albedo (fraction)
    real, intent(in) :: cloud_cover         !< Cloud cover (fraction)
    real, intent(in) :: balance_threshold   !< Threshold within which to accept leaf energy balance
    real, intent(in) :: adjustment_factor   !< Multiplier for energy balance based T_leaf adjustments
    integer, intent(in) :: max_iterations   !< Maximum iteration count

    real :: T_leaf    !< Output: Leaf temperature (degrees C)

    ! Stefan-Boltzmann constant, W m-2 K-4 (Campbell & Norman (1998), p281, table A5)
    real, parameter :: SBC = 5.670373e-8
    ! Specific heat capacity of dry air, J mol-1 C-1 (Campbell & Norman (1998), p279, table A1)
    real, parameter :: c_p = 29.3

    ! Leaf short-wave absorptivity (de Boeck (2012); from Campbell & Norman (1998), p153, in text)
    real, parameter :: alpha_s_leaf = 0.5
    ! Leaf long-wave absorptivity (de Boeck, 2012)
    real, parameter :: alpha_l_leaf = 0.97
    ! Soil long-wave emissivity (de Boeck, 2012)
    real, parameter :: eps_l_soil = 0.945
    ! Leaf long-wave emissivity (de Boeck, 2012)
    real, parameter :: eps_l_leaf = 0.97

    real :: rho_s_soil, VPL, eps_ac, eps_l_sky, lam, T_soil, T_leaf_adjust, g_Ha, g_vb, g_v, e_s_Tleaf
    real :: R_s_in, R_l_in, R_l_out, H, lam_E, energy_balance
    integer :: i

    ! Soil short-wave reflectivity
    rho_s_soil = albedo
    ! Water vapour path length (de Boeck 2012)
    VPL = 46.5 * ((0.01 * e_a) / (T_air + T0))
    ! Clear sky emissivity (de Boeck 2012)
    eps_ac = 1 - (1 + VPL) * exp(-(1.2 + 3.0 * VPL)**0.5)
    ! Sky long-wave emmisivity (de Boeck (2012); from Campbell & Norman (1998), p164, eq. 10.12)
    eps_l_sky = (1 - 0.84 * cloud_cover) * eps_ac + 0.84 * cloud_cover
    ! Latent heat of vapourisation, J mol-1 (de Boeck (2012))
    lam = -42.575 * T_air + 44994

    ! Assume soil temperature is equal to air temperature
    ! TODO: use the banded estimate from de Boeck (2012)?
    T_soil = T_air

    ! Starting point
    T_leaf = initial_T_leaf
    T_leaf_adjust = 0.0

    do i = 1, max_iterations
      ! Apply adjustment (from previous iteration)
      T_leaf = T_leaf + T_leaf_adjust

      ! Boundary layer conductances (de Boeck (2012))
      ! XXX: Assume forced convection:
      g_Ha = 1.4 * 0.135 * sqrt(u / d)
      g_vb = 1.4 * 0.147 * sqrt(u / d)

      ! Total conductivity to water vapour
      if (hypostomatous) then
        g_v = combined_leaf_conductance(g_vs, 0.0, g_vb)
      else
        g_v = combined_leaf_conductance(g_vs, g_vs, g_vb)
      end if

      ! Saturated vapour pressure for leaf temperature, Pa
      e_s_Tleaf = saturated_vapour_pressure(T_leaf) * 1000  ! Converted from kPa to Pa

      ! Absorbed short-wave (direct on top of leaf, reflected from soil on underside) (de Boeck (2012))
      R_s_in = (0.5 * R * alpha_s_leaf) + (0.5 * R * rho_s_soil * alpha_s_leaf)

      ! Absorbed long-wave (emitted from soil on underside, emitted from sky on top) (de Boeck (2012))
      R_l_in = (0.5 * alpha_l_leaf * eps_l_soil * SBC * (T_soil + T0)**4) + &
               (0.5 * alpha_l_leaf * eps_l_sky * SBC * (T_air + T0)**4)

      ! Emitted long-wave (de Boeck (2012))
      R_l_out = eps_l_leaf * SBC * (T_leaf + T0)**4

      ! Sensible heat flux (de Boeck (2012))
      H = g_Ha * c_p * (T_leaf - T_air)

      ! Latent heat flux (de Boeck (2012))
      lam_E = lam * g_v * (e_s_Tleaf - e_a) / P

      ! Energy balance equation (de Boeck (2012))
      energy_balance = R_s_in + R_l_in - R_l_out - H - lam_E

      if (abs(energy_balance) <= balance_threshold) then
        exit
      else
        T_leaf_adjust = energy_balance * adjustment_factor
      end if
    end do
  contains
    pure function combined_leaf_conductance(gs_ab, gs_ad, ga) result(g)
      real, intent(in) :: gs_ab   !< Abaxial stomatal conductance (mol m-2 s-1)
      real, intent(in) :: gs_ad   !< Adaxial stomatal conductance (mol m-2 s-1)
      real, intent(in) :: ga      !< Boundary layer conductance (mol m-2 s-1)
      real :: g     !< Output: combined leaf conductance (mol m-2 s-1)

      g = (0.5 * gs_ab * ga / (gs_ab + ga)) + (0.5 * gs_ad * ga / (gs_ad + ga))
    end function combined_leaf_conductance
  end function leaf_temp_de_Boeck
  

  !> Calculate parameter from temperature dependence curve
  pure function temp_dep(P_ref, T_ref, H_a, T) result (P)
    real, intent(in) :: P_ref   ! Parameter value at T_ref
    real, intent(in) :: T_ref   ! Reference temperature (degrees K)
    real, intent(in) :: H_a     ! Activation energy (J/mol)
    real, intent(in) :: T       ! Temperature (degrees K)

    real :: P

    P = P_ref * exp((H_a * (T - T_ref)) / (T_ref * R * T))
  end function temp_dep
  
  
  !> Calculate parameter from temperature dependence curve with high
  !! temperature inhibition
  pure function temp_dep_inhibit(P_ref, T_ref, H_a, H_d, S, T) result (P)
    real, intent(in) :: P_ref   ! Parameter value at T_ref
    real, intent(in) :: T_ref   ! Reference temperature (degrees K)
    real, intent(in) :: H_a     ! Activation energy (J/mol)
    real, intent(in) :: H_d     ! Deactivation energy (J/mol)
    real, intent(in) :: S       ! Entropy term (J/(mol*K))
    real, intent(in) :: T       ! Temperature (degrees K)

    real :: P

    P = P_ref * exp((H_a * (T - T_ref)) / (T_ref * R * T)) &
        * (  (1 + exp((T_ref*S - H_d) / (T_ref*R))) &
           / (1 + exp((T*S - H_d) / (T*R))) )
  end function temp_dep_inhibit



END MODULE
