!
!////////////////////////////////////////////////////////////////////////
!
!      soil_moisture.f90
!      Created: 16 August 2016 10:22 
!      By: David Gillies  
!
!////////////////////////////////////////////////////////////////////////
!
MODULE DO3SE_Soilmoisture
USE DO3SE_Met

  TYPE, public :: SMDData_t
    REAL :: Sn = UNDEF        !< Soil water content (m3 m-3)
    REAL :: SWP = UNDEF       !< Soil water potential (MPa)
    REAL :: ASW = UNDEF       !< Available soil water, above PWP (m)
    REAL :: SMD = UNDEF       !< Soil moisture deficit, below FC (m)
  end type SMDData_t

CONTAINS

  subroutine calc_soil_moisture(source, SWP, SWP_AE, b, ASW, SMD, PWP, FC, root, &
                                Sn, hr)
    !class(DO3SE_State_t), intent(inout) :: this
    CHARACTER(len=16), INTENT(in) :: source
    real, intent(in) :: SWP_AE
    real, intent(in) :: b
    REAL, INTENT(in) :: PWP
    REAL, INTENT(in) :: FC
    REAL, INTENT(in) :: root
    REAL, INTENT(inout) :: ASW
    REAL, INTENT(inout) :: SMD
    REAL, INTENT(inout) :: SWP
    REAL, INTENT(inout) :: Sn
    INTEGER, INTENT(in) :: hr
    type(SMDData_t) :: data
    select case (source)
    case ("disabled")
      ! Nothing to do
    case ("input SWP")
      data = soil_moisture_from_SWP(SWP, SWP_AE, b, ASW, SMD, PWP, FC, root)
    case ("input SWC")
      data = soil_moisture_from_SWC(SWP, SWP_AE, b, ASW, SMD, PWP, FC, root, Sn)
    case ("P-M")
      if (start_of_day(hr)) then
        ! Apply change in soil moisture from previous day's P-M summary.
        data = soil_moisture_from_SWC(SWP, SWP_AE, b, ASW, SMD, PWP, FC, root, Sn + Sn_diff)
      end if
    case default
      !UNKNOWN_STRING(this%SMD_conf%source)
    end SELECT
    
    Sn = DATA%Sn
    SWP = DATA%SWP
    ASW = DATA%ASW
    SMD = DATA%SMD
    
  end subroutine calc_soil_moisture

  !> Fill soil moisture data from a soil water potential value.
  pure function soil_moisture_from_SWP(SWP, SWP_AE, b, ASW, SMD, PWP, FC, root) result(data)
    !type(SMDConfig_t), intent(in) :: config
    real, intent(in) :: SWP     !< Soil water potential (MPa)
    real, intent(in) :: SWP_AE
    real, intent(in) :: b
    REAL, INTENT(in) :: PWP
    REAL, INTENT(in) :: FC
    REAL, INTENT(in) :: root
    REAL, INTENT(IN) :: ASW
    REAL, INTENT(in) :: SMD
    type(SMDData_t) :: data

    !associate (soil => config%soil)
    data = soil_moisture_from_SWC(SWP, ASW, SMD, SWP_AE, b, PWP, FC, root, SWP_to_SWC(SWP_AE, b, SWP))
    !end associate
  end function soil_moisture_from_SWP
  
  !> Fill soil moisture data from a soil water content value.
  pure function  soil_moisture_from_SWC(SWP, ASW, SMD, SWP_AE, b, PWP, FC, root, Sn) RESULT (DATA)
    !type(SMDConfig_t), intent(in) :: config
    
    real, intent(in) :: Sn      !< Soil water content (m3 m-3)
    REAL, INTENT(in) :: SWP_AE
    REAL, INTENT(in) :: b
    REAL, INTENT(in) :: PWP
    REAL, INTENT(in) :: FC
    REAL, INTENT(in) :: root
    REAL, INTENT(IN) :: ASW
    REAL, INTENT(in) :: SMD
    real, intent(in) :: SWP     !< Soil water potential (MPa)
    type(SMDData_t) :: data

    real :: PWP_vol

    !associate (soil => config%soil)
      ! Convert PWP to volumetric content to use as a minimum soil water content
      PWP_vol = SWP_to_SWC(SWP_AE, b, PWP)

      ! Constrain soil water content to be between field capacity and PWP
      DATA%Sn = max(PWP_vol, min(FC, Sn))

      ! Calculate soil water potential (SWP)
      DATA%SWP = SWC_to_SWP(SWP_AE, b, DATA%Sn)

      ! Calculate available soil water (ASW)
      DATA%ASW = (DATA%Sn - PWP_vol) * root

      ! Calculate soil moisture deficit (SMD)
      DATA%SMD = (FC - DATA%Sn) * root
    !end associate
  end function  soil_moisture_from_SWC
  
  
  !> Use soil water release curve to convert from soil water potential (MPa) to
  !! soil water content (m3 m-3).
  pure real function SWP_to_SWC(SWP_AE, b, SWP)
    real, intent(in) :: SWP_AE    !< Water potential at air entry (MPa)
    real, intent(in) :: b         !< Texture dependent soil conductivity parameter
    real, intent(in) :: SWP       !< Soil water potential (MPa)

    real, parameter :: SWC_sat = 0.4 ! Saturated soil water content for soil water release curve

    SWP_to_SWC = 1.0 / (((SWP/SWP_AE)**(1.0/b)) / SWC_sat)
  end function SWP_to_SWC
  
  !> Use soil water release curve to convert from soil water content (m3 m-3) to
  !! soil water potential (MPa).
  pure real function SWC_to_SWP(SWP_AE, b, SWC)
    real, intent(in) :: SWP_AE    !< Water potential at air entry (MPa)
    real, intent(in) :: b         !< Texture dependent soil conductivity parameter
    real, intent(in) :: SWC       !< Soil water content (m3 m-3)

    real, parameter :: SWC_sat = 0.4 ! Saturated soil water content for soil water release curve

    SWC_to_SWP = SWP_AE * ((SWC_sat / SWC)**b)
  end function SWC_to_SWP
  

  
END MODULE
