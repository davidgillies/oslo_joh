!
!////////////////////////////////////////////////////////////////////////
!
!      phenology.f90
!      Created: 12 August 2016 13:45 
!      By: David Gillies  
!
!////////////////////////////////////////////////////////////////////////
!
MODULE DO3SE_Phenology


CONTAINS
subroutine calc_phenology(height_method, height, LAI_method, SGS, LAI_a, LAI_b, LAI_1, &
                          LAI_c, EGS, LAI_d, LAI_2, dd, SAI_method, &
                           season_SAI_method, canopy_height, LAI, SAI)
   CHARACTER(len=16), INTENT(IN) :: height_method
   REAL, INTENT(IN) :: height
   CHARACTER(len=16), INTENT(IN) :: LAI_method
   INTEGER, INTENT(IN) :: SGS
   REAL, INTENT(IN) :: LAI_a
   REAL, INTENT(IN) :: LAI_b
   INTEGER, INTENT(IN) :: LAI_1
   REAL, INTENT(IN) :: LAI_c
   INTEGER, INTENT(IN) :: EGS
   REAL, INTENT(IN) :: LAI_d
   INTEGER, INTENT(IN) :: LAI_2
   INTEGER, INTENT(IN) :: dd
   CHARACTER(len=16), INTENT(IN) :: SAI_method
   CHARACTER(len=16), INTENT(IN) :: season_SAI_method
   REAL, INTENT(out) :: canopy_height
   REAL, INTENT(out) :: LAI
   REAL, INTENT(out) :: SAI

   select case (height_method)
      case ("input")
        ! Nothing to do
      case ("constant")
        ! Use height of primary land cover
        canopy_height = height
      case default
        !UNKNOWN_STRING(height_method)
      end select

      ! Calculate height of top of each canopy
      !this%ML(:)%layer_height = this%V%canopy_height * this%LC_conf%layer_height(:this%nL)

      select case (LAI_method)
      case ("input")
        ! Nothing to do
      case ("estimate total")
        ! Use primary land cover's estimate of total LAI
        LAI = LAI_day_PLF(SGS, LAI_a, LAI_b, LAI_1, &
                          LAI_c, EGS, LAI_d, LAI_2, dd)
        ! Spread single LAI value to layers and LCs
      case default
        !UNKNOWN_STRING(LAI_method)
      end select

      select case (SAI_method)
      case ("input")
        ! Nothing to do
      case ("input total")
        ! Spread single SAI value to layers and LCs
        !this%MLMC(:,:)%SAI = this%MLMC(1,1)%SAI * this%fLAI(:,:)
      case ("estimate total")
        ! Use primary land cover's estimate of total SAI
        select case (season_SAI_method)
        case ("LAI")
          SAI = LAI
        case ("forest")
          SAI = LAI + 1.0
        case ("wheat")
          SAI = SAI_wheat(SGS, EGS, LAI_1, dd, LAI_day_PLF(SGS, LAI_a, LAI_b, LAI_1, &
                          LAI_c, EGS, LAI_d, LAI_2, dd))
        case default
          !UNKNOWN_STRING(SAI_method)
        end select
        ! Spread single SAI value to layers and LCs
        !this%MLMC(:,:)%SAI = this%MLMC(1,1)%SAI * this%fLAI(:,:)
      case default
        !UNKNOWN_STRING(SAI_method)
      end select

      ! Calculate the distribution of LAI between land covers
      !if LAI <= 0.0 then
      !  LC_dist = 1.0
      !else
      !  LC_dist = 1.0
      !end if

  end subroutine calc_phenology
  
  pure subroutine reindex(a, zero, wrap)
    real, intent(inout) :: a(1:)
    real, intent(in) :: zero
    real, intent(in), optional :: wrap

    a = a - zero
    if (present(wrap)) then
      where (a < 0.0) a = a + wrap
    end if
  end subroutine reindex
   
   real function LAI_day_PLF(SGS, LAI_a, LAI_b, LAI_1, &
                          LAI_c, EGS, LAI_d, LAI_2, dd)
   INTEGER, INTENT(IN) :: SGS
   REAL, INTENT(IN) :: LAI_a
   REAL, INTENT(IN) :: LAI_b
   INTEGER, INTENT(IN) :: LAI_1
   REAL, INTENT(IN) :: LAI_c
   INTEGER, INTENT(IN) :: EGS
   REAL, INTENT(IN) :: LAI_d
   INTEGER, INTENT(IN) :: LAI_2
   INTEGER, INTENT(IN) :: dd

    real, dimension(2, 5) :: func
    integer :: dd_adj

    ! Build function
    func = reshape((/ real :: &
      SGS, LAI_a, &
      (SGS + LAI_1), LAI_b, &
      (EGS - LAI_2), LAI_c, &
      EGS, LAI_d, &
      (SGS + 365), LAI_a /),  shape(func))
    ! Re-index everything to SGS = 0, wrapping dates before SGS to the end of the year
    call reindex(func(1,:), real(SGS), 365.0)
    !call assert(all(func(1,1:4) <= func(1,2:5)), "LAI_day_PLF: points not in order")
    dd_adj = dd - SGS
    if (dd_adj < 0) then
      dd_adj = dd_adj + 365
    end if
    ! Lookup value in PLF
    LAI_day_PLF = PLF_value(func, real(dd_adj))
  end function LAI_day_PLF
  
  real function SAI_wheat(SGS, EGS, LAI_1, dd, LAI)
   INTEGER, INTENT(IN) :: SGS
   INTEGER, INTENT(IN) :: LAI_1
   REAL, INTENT(IN) :: LAI
   INTEGER, INTENT(IN) :: EGS
    integer, intent(in) :: dd         !< Day of year        

    real, dimension(2, 4) :: func
    integer :: i, dd_adj
    
    SAI_wheat = UNDEF

    func = reshape((/ real :: &
      SGS, LAI, &
      (SGS + LAI_1), LAI + ((5.0/3.5) - 1) * LAI, &
      EGS + 1, LAI + 1.5, &
      SGS + 365, LAI /), shape(func))
    call reindex(func(1,:), real(SGS), 365.0)
    dd_adj = dd - SGS
    if (dd_adj < 0) then
      dd_adj = dd_adj + 365
    end if

    do i = 1, 4
      if (dd_adj < func(1,i)) then
        SAI_wheat = func(2,i)
        exit
      end if
    end do

    !call assert(is_def(SAI_wheat), "SAI_wheat: no result")
  end function SAI_wheat
  
  pure function PLF_value(points, x) result(y)
    real, dimension(:,:), intent(in) :: points
    real, intent(in) :: x
    real :: y

    ! TODO: sanity-check points: should be size 2 in first dimension, and
    !       x values should never decrease.

    integer :: i, n
    real :: ax, ay, bx, by

    n = ubound(points, 2)

    if (x < points(1,1)) then
      y = points(2,1)
    else if (x > points(1,n)) then
      y = points(2,n)
    else
      bx = points(1,1)
      by = points(2,1)

      do i = 2, n
        ax = bx
        ay = by
        bx = points(1,i)
        by = points(2,i)

        ! Skip zero-width pieces (this should be equivalent to an
        ! equality check, but checking floating point equality is evil
        ! and the compiler warns about it)
        if (abs(ax - bx) < epsilon(ax)) then
          cycle
        end if

        if (x <= bx) then
          y = ay + (by - ay) * ((x - ax) / (bx - ax))
          exit
        end if
      end do
    end if
  end function PLF_value

END MODULE
