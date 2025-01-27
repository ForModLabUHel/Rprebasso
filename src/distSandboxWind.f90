subroutine windrisk(siteInfoDist, spec, h, openedge, sitetype, tsum, tsincethin, &
  wrisk5dd1, wrisk5dd2, wrisk5dd3, wrisk0, wrisk5, wrisk)
  IMPLICIT NONE
  REAL (kind=8), intent(inout) ::  siteInfoDist(10) ! 5-year wind risk (suvanto output), pre-logit value, annual risk
  REAL (kind=8), intent(inout) ::  wrisk5, wrisk0, wrisk ! 5-year wind risk (suvanto output), pre-logit value, annual risk
  REAL (kind=8), intent(inout) :: wrisk5dd1, wrisk5dd2, wrisk5dd3 !5-year wind risk of each damage density class
  REAL (kind=8), intent(in) :: h ! input in m, converted to dm
  REAL (kind=8), intent(in) :: tsum ! effective temperature sums (degree days over 5°C); converted to 100 dd
  INTEGER, intent(in) :: spec !1 pine, 2 spruce, 3 other
  INTEGER, intent(in) :: openedge ! 0 = no open edge, 1 = open edge
  INTEGER, intent(in) :: sitetype ! prebas site types; reclassified to site fertility: 0 = infertile, 1 = fertile
  REAL (kind=8) :: wspeed ! localised 10 a max windspeed (m/s), Venäläinen 2017
  INTEGER, intent(in) :: tsincethin ! time since last thinning in years, categorised into 0-5, 6-10, >10 below
  INTEGER :: soiltype ! 0 = mineral, coarse; 1 = mineral, fine; 2 = organic
  INTEGER :: shallowsoil ! 1 = <30cm
  INTEGER :: sitefert ! site fertility, reclassified from sitetype (1:3 fertile / 1, 4:5 infertile/0)

  wspeed = siteInfoDist(1)
  !tsincethin = INT(siteInfoDist(2))
  soiltype = INT(siteInfoDist(3))
  shallowsoil = INT(siteInfoDist(4))

    IF (sitetype <= 3) sitefert = 1 !convert sitetypes to fert class
    IF (sitetype > 3) sitefert = 0

    wrisk0 = -14.690 + &
                  LOG(h*10)*1.661 + &
                  LOG(wspeed)*0.749 + &
                  tsum/100*0.096 + &!effective temp sum (100 degree days)
                  openedge * 0.310 + &
                  shallowsoil * 0.214 - &
                  sitefert * 0.425

    !categorical variables with more than two levels externalised:

    !!! time since thinning; reference: 0-5a
    if (tsincethin>5 .AND. tsincethin <= 10) wrisk0 = wrisk0 - 0.298
    if (tsincethin>10)  wrisk0 = wrisk0 - 0.844
    ! soiltype; reference: mineral, coarse
    if (soiltype == 1)  wrisk0 = wrisk0 - 0.356 !mineral, fine
    if (soiltype == 2)  wrisk0 = wrisk0 - 0.216!organic

    !!! species & spec/h interaction (ref: pine/1)
    IF (spec == 2) THEN
      wrisk0 = wrisk0 - 8.494
      wrisk0 = wrisk0 + LOG(h*10)*1.634
    ELSEIF (spec == 3) THEN
      wrisk0 = wrisk0 - 9.314
      wrisk0 = wrisk0 + LOG(h*10)*1.625
    ENDIF

    ! DAMAGE DENSITY (& logit transformation)
    ! spatial density of wind-disturbed NFI plots relative to all NFI plots
    !!! SEE SUPPLEMENTYRY INFO 1 of Suvanto et al. 2019
    ! for prediction, the model is run for each damage density class individually
    ! and a weighted average of the three is used; weights: damdens 0-2,2-3,>3: 0.905, 0.072, 0.023

    !damdens 0-2 (reference)
    wrisk5dd1 = EXP(wrisk0) / (1.0 + EXP(wrisk0)) ! logit transformation of reference
    !damdens 2-3
    wrisk0 = wrisk0 + 1.104
    wrisk5dd2 = EXP(wrisk0) / (1.0 + EXP(wrisk0))
    !damdens >3
    wrisk0 = wrisk0 + 1.898 - 1.104 !to avoid another variable
    wrisk5dd3 = EXP(wrisk0) / (1.0 + EXP(wrisk0))

    wrisk5 = wrisk5dd1 * 0.905 +  wrisk5dd2 * 0.072 + wrisk5dd3 * 0.023 !weighted average
    wrisk = wrisk5/5 ! annual risk

end subroutine


! LEGACY: version with individual inputs; new above: external inputs in wDistSiteInfo
subroutine windriskold(spec, h, tsincethin, wspeed, openedge, soiltype, shallowsoil, sitetype, tsum, &
  wrisk5dd1, wrisk5dd2, wrisk5dd3, wrisk0, wrisk5, wrisk)
  IMPLICIT NONE
  REAL (kind=8), intent(inout) ::  wrisk5, wrisk0, wrisk ! 5-year wind risk (suvanto output), pre-logit value, annual risk
  REAL (kind=8), intent(inout) :: wrisk5dd1, wrisk5dd2, wrisk5dd3 !5-year wind risk of each damage density class
  REAL (kind=8), intent(in) :: h ! input in m, converted to dm
  REAL (kind=8), intent(in) :: wspeed ! localised 10 a max windspeed (m/s), Venäläinen 2017
  REAL (kind=8), intent(in) :: tsum ! effective temperature sums (degree days over 5°C); converted to 100 dd
  INTEGER, intent(in) :: spec !1 pine, 2 spruce, 3 other
  INTEGER, intent(in) :: tsincethin ! time since last thinning in years, categorised into 0-5, 6-10, >10 below
  INTEGER, intent(in) :: openedge ! 0 = no open edge, 1 = open edge
  INTEGER, intent(in) :: soiltype ! 0 = mineral, coarse; 1 = mineral, fine; 2 = organic
  INTEGER, intent(in) :: shallowsoil ! 1 = <30cm
  INTEGER, intent(in) :: sitetype ! prebas site types; reclassified to site fertility: 0 = infertile, 1 = fertile
  INTEGER :: sitefert ! site fertility, reclassified from sitetype (1:3 fertile / 1, 4:5 infertile/0)

    IF (sitetype <= 3) sitefert = 1 !convert sitetypes to fert class
    IF (sitetype > 3) sitefert = 0

    wrisk0 = -14.690 + &
                  LOG(h*10)*1.661 + &
                  LOG(wspeed)*0.749 + &
                  tsum/100*0.096 + &!effective temp sum (100 degree days)
                  openedge * 0.310 + &
                  shallowsoil * 0.214 - &
                  sitefert * 0.425

    !categorical variables with more than two levels externalised:

    !!! time since thinning; reference: 0-5a
    if (tsincethin>5 .AND. tsincethin <= 10) wrisk0 = wrisk0 - 0.298
    if (tsincethin>10)  wrisk0 = wrisk0 - 0.844
    ! soiltype; reference: mineral, coarse
    if (soiltype == 1)  wrisk0 = wrisk0 - 0.356 !mineral, fine
    if (soiltype == 2)  wrisk0 = wrisk0 - 0.216!organic

    !!! species & spec/h interaction (ref: pine/1)
    IF (spec == 2) THEN
      wrisk0 = wrisk0 - 8.494
      wrisk0 = wrisk0 + LOG(h*10)*1.634
    ELSEIF (spec == 3) THEN
      wrisk0 = wrisk0 - 9.314
      wrisk0 = wrisk0 + LOG(h*10)*1.625
    ENDIF

    ! DAMAGE DENSITY (& logit transformation)
    ! spatial density of wind-disturbed NFI plots relative to all NFI plots
    !!! SEE SUPPLEMENTYRY INFO 1 of Suvanto et al. 2019
    ! for prediction, the model is run for each damage density class individually
    ! and a weighted average of the three is used; weights: damdens 0-2,2-3,>3: 0.905, 0.072, 0.023

    !damdens 0-2 (reference)
    wrisk5dd1 = EXP(wrisk0) / (1.0 + EXP(wrisk0)) ! logit transformation of reference
    !damdens 2-3
    wrisk0 = wrisk0 + 1.104
    wrisk5dd2 = EXP(wrisk0) / (1.0 + EXP(wrisk0))
    !damdens >3
    wrisk0 = wrisk0 + 1.898 - 1.104 !to avoid another variable
    wrisk5dd3 = EXP(wrisk0) / (1.0 + EXP(wrisk0))

    wrisk5 = wrisk5dd1 * 0.905 +  wrisk5dd2 * 0.072 + wrisk5dd3 * 0.023 !weighted average
    wrisk = wrisk5/5 ! annual risk

end subroutine





! subroutine prioritiseDistReact(siteOrder, outDist, nsites)
!   IMPLICIT NONE
!   REAL (kind=8), intent(inout) ::  outDist(nSites, 10)
!   INTEGER, intent(inout) ::  siteOrder(nSites)
!   INTEGER, intent(in) :: nSites !1 pine, 2 spruce, 3 o
!
!
! do i = 1, nSites
!   if (outDist(i, 7) == 1) then
!     siteOrder()
!
!
! end do
!
!
! end subroutine
!
! !
! subroutine move_element_to_front(vector, index, n)
!     implicit none
!     integer, intent(inout) :: vector(:)
!     integer, intent(in) :: index, n
!     integer :: temp, i
!
!     ! Check if index is valid
!     if (index < 1 .or. index > n) then
!         print *, "Error: Index out of bounds"
!         return
!     end if
!
!     ! Move the element to the front
!     temp = vector(index)
!     do i = index, 2, -1
!         vector(i) = vector(i-1)
!     end do
!     vector(1) = temp
! end subroutine move_element_to_front


subroutine move_element_to_front(siteorder, index, nsites, siteordertemp)
    implicit none
    integer, intent(inout) :: siteorder(nsites)
    integer, intent(in) :: index, nsites
    integer :: temp, i
    integer, intent(inout) :: siteordertemp(nsites)

    ! Check if index is valid
    if (index < 1 .or. index > nsites) then
        print *, "Error: Index out of bounds"
        return
    end if

    ! Move the element to the front
    siteordertemp(1) = siteorder(index) ! put focus site id to top
    do i = 1, index-1 ! shift all siteids in siteorder prior to index one down
        siteordertemp(i+1) = siteorder(i)
    end do
    siteordertemp((index+1):nsites) = siteorder((index+1):nsites) ! keep remaining siteorder as is
    !vector(1) = temp
end subroutine move_element_to_front


subroutine move_siteid_to_top(siteorder, siteid, nsites, siteordertemp, poutdist)
    implicit none
    integer, intent(inout) :: siteorder(nsites)
    integer, intent(in) :: siteid, nsites
    integer :: temp, i, index(1)
    integer, intent(inout) :: siteordertemp(nsites)
real (kind=8), intent(inout) :: poutdist(nsites, 10)
    !index = findloc(siteorder, siteid) !find location of siteid in question within siteorder
    siteordertemp = abs(siteorder-siteid)
    index = minloc(siteordertemp(:)) !find location of siteid in question within siteorder !! findloc only in fortran 2008 and later, workaround with abs/minloc

    ! Move the element to the front
    siteordertemp(1) = siteorder(index(1)) ! put focus site id to top
    do i = 1, index(1)-1 ! shift all siteids in siteorder prior to index one down
        siteordertemp(i+1) = siteorder(i)
    end do
    siteordertemp((index(1)+1):nsites) = siteorder((index(1)+1):nsites) ! keep remaining siteorder as is
end subroutine move_siteid_to_top

!
! subroutine move_siteid_to_top2(siteorder, siteid, nsites, siteordertemp, poutdist)
!     implicit none
!     integer, intent(inout) :: siteorder(nsites)
!     integer, intent(inout) :: nsites
!     integer :: temp, i, index(1), priosites(nsites), nsitesdist,  siteid
!     integer, intent(inout) :: siteordertemp(nsites)
! real (kind=8), intent(inout) :: poutdist(nsites, 10)
!
!   ! here: finding siteids of distcc-earmarked sites
!   ! for now limited to a single site
!
!     !nsitesdist = sum(poutdist(:,7))
!
!      priosites = findloc(poutdist(:,7), 1)
!
!     siteid = priosites(1)
!
!
!     ! below working if siteid is given
!     !index = findloc(siteorder, siteid) !find location of siteid in question within siteorder
!     siteordertemp = abs(siteorder-siteid)
!     index = minloc(siteordertemp(:)) !find location of siteid in question within siteorder !! findloc only in fortran 2008 and later, workaround with abs/minloc
!
!     ! Move the element to the front
!     siteordertemp(1) = siteorder(index(1)) ! put focus site id to top
!     do i = 1, index(1)-1 ! shift all siteids in siteorder prior to index one down
!         siteordertemp(i+1) = siteorder(i)
!     end do
!     siteordertemp((index(1)+1):nsites) = siteorder((index(1)+1):nsites) ! keep remaining siteorder as is
!
!     nsites = siteid
!
! end subroutine move_siteid_to_top2

subroutine find_row_indexes(matrix, n, row_indexes)
    implicit none
    integer, intent(in) :: matrix(n, 10)
    integer, intent(in) :: n
    integer, intent(inout) :: row_indexes(n)
    integer :: i, count

    count = 0
    do i = 1, n
        if (INT(matrix(i, 7)) == 1) then
            count = count + 1
            row_indexes(count) = i
        end if
    end do

    ! Fill remaining elements of row_indexes with 0
    do i = count + 1, n
        row_indexes(i) = 0
    end do
end subroutine find_row_indexes


subroutine prioDistInSO(outDist, nSites, maxYears, year, siteOrder)
! prioritise sites flagged for management reaction to wind disturbance in siteOrder (annualy randomised site order)
! puts sites flagged for mgmt reaction to disturbance in last year's outDist on top of site order
    implicit none
    integer, intent(in) :: nSites, year, maxYears
    integer, intent(inout) :: siteOrder(nSites,maxYears)! , siteOrderX(nSites)
    real (kind=8), intent(in) :: outDist(nSites, 10)
    integer :: sitexx, distpriositexx, ndistprio, temp, index(1), siteid, siteordertemp(nSites), priosites(nSites)
    ! fill priosites with siteids (outdist in 1:nsites order)
    ndistprio = 0

    do sitexx = 1, nSites
        if (outDist(sitexx, 8) == 1.) then
            ndistprio = ndistprio + 1
            priosites(ndistprio) = sitexx
        end if
    end do
    ! for each of these, put siteid on top and shift those siteids above in sitorder one down
    !write(1,*) ndistprio !tswrite
    ! siteOrder(2,year) = ndistprio !troubleshooting:
    do distpriositexx = 1, ndistprio
      siteid = priosites(distpriositexx)
          siteordertemp = abs(siteOrder(:,year)-siteid)
           index = minloc(siteordertemp(:)) !find location of siteid in question within siteorder !! findloc only in fortran 2008 and later, workaround with abs/minloc
          ! Move the element to the front
          siteordertemp(1) = siteOrder(index(1), year) ! put focus site id to top
          do sitexx = 1, index(1)-1 ! shift all siteids in siteorder prior to index one down
              siteordertemp(sitexx+1) = siteOrder(sitexx, year)
          end do
          siteordertemp((index(1)+1):nsites) = siteOrder((index(1)+1):nsites, year) ! keep remaining siteorder as is
          !siteOrder(1,year) = year+1

          siteOrder(:,year) = siteordertemp(:)
    end do
    ! to see if anything happening here is in output

      end subroutine prioDistInSO

!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!
! sample wind disturbance impact (relative damaged vol)
! based on pre-calculated probabilities (lognormal fitted to 2001 post-storm inventory data)


      subroutine sample_rdvol(rdvol_sampled)

        implicit none

        ! Declare variables
        real(8) :: cdf(1000, 2)        ! Matrix to hold the CDF (2 columns: probabilities and values)
        real(8):: rndm    ! Random number in the range [0, 1]
        real(8), intent(inout) :: rdvol_sampled    ! Sampled value based on the random number
        integer, parameter :: n = 1000 ! Number of rows in the CDF
        integer :: i             ! Loop counter

        ! Initialize the CDF matrix
        ! cdf(:, 1) = (/ 0.009364, 0.043472, 0.090473, 0.141734, 0.193028, 0.24238, 0.288928, 0.332368, 0.372677, &
        ! 0.409973, 0.444439, 0.476282, 0.505714, 0.532939, 0.55815, 0.581524, 0.603223, 0.623394, 0.64217,  &
        ! 0.659671, 0.676005, 0.69127, 0.705555, 0.718938, 0.731493, 0.743283, 0.754369, 0.764803, 0.774635,  &
        ! 0.783908, 0.792663, 0.800937, 0.808764, 0.816173, 0.823195, 0.829853, 0.836173, 0.842176, 0.847883,  &
        ! 0.853312, 0.858479, 0.863402, 0.868095, 0.872572, 0.876844, 0.880925, 0.884824, 0.888553, 0.89212,  &
        ! 0.895535, 0.898805, 0.901939, 0.904943, 0.907825, 0.910591, 0.913246, 0.915796, 0.918247, 0.920603,  &
        ! 0.922869, 0.925049, 0.927148, 0.929169, 0.931116, 0.932992, 0.9348, 0.936544, 0.938227, 0.93985,  &
        ! 0.941418, 0.942931, 0.944393, 0.945806, 0.947172, 0.948492, 0.949769, 0.951005, 0.952201, 0.953358,  &
        ! 0.954478, 0.955564, 0.956615, 0.957634, 0.958622, 0.959579, 0.960508, 0.961409, 0.962283, 0.963131,  &
        ! 0.963954, 0.964753, 0.965529, 0.966283, 0.967015, 0.967727, 0.968419, 0.969091, 0.969745, 0.970381,  &
        ! 0.971, 0.971601, 0.972187, 0.972757, 0.973311, 0.973852, 0.974378, 0.97489, 0.975389, 0.975875, 0.976349,  &
        ! 0.97681, 0.97726, 0.977699, 0.978127, 0.978544, 0.978951, 0.979348, 0.979736, 0.980114, 0.980482, 0.980843,  &
        ! 0.981194, 0.981537, 0.981872, 0.9822, 0.982519, 0.982832, 0.983137, 0.983435, 0.983726, 0.984011, 0.98429, &
        !  0.984562, 0.984828, 0.985088, 0.985343, 0.985592, 0.985835, 0.986074, 0.986307, 0.986535, 0.986758, 0.986977, &
        !   0.987191, 0.9874, 0.987605, 0.987806, 0.988002, 0.988195, 0.988384, 0.988568, 0.988749, 0.988927, 0.989101, &
        !    0.989271, 0.989438, 0.989601, 0.989762, 0.989919, 0.990073, 0.990365, 0.990654, 0.99094, 0.991223, 0.991503, &
        !     0.991781, 0.992056, 0.992329, 0.992599, 0.992867, 0.993132, 0.993395, 0.993656, 0.993914, 0.994171, 0.994425, &
        !      0.994677, 0.994927, 0.995176, 0.995422, 0.995666, 0.995909, 0.99615, 0.996389, 0.996626, 0.996861, 0.997095, &
        !       0.997328, 0.997558, 0.997787, 0.998015, 0.998241, 0.998466, 0.998689, 0.998911, 0.999131, 0.99935,  &
        !       0.999568, 0.999785, 1. /)  ! Cumulative probabilities
        ! cdf(:, 2) = (/ 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.055, 0.06, 0.065, &
        ! 0.07, 0.075, 0.08, 0.085, 0.09, 0.095, 0.1, 0.105, 0.11, 0.115, 0.12, 0.125, 0.13, 0.135, 0.14, &
        ! 0.145, 0.15, 0.155, 0.16, 0.165, 0.17, 0.175, 0.18, 0.185, 0.19, 0.195, 0.2, 0.205, 0.21, 0.215, &
        ! 0.22, 0.225, 0.23, 0.235, 0.24, 0.245, 0.25, 0.255, 0.26, 0.265, 0.27, 0.275, 0.28, 0.285, 0.29, &
        !  0.295, 0.3, 0.305, 0.31, 0.315, 0.32, 0.325, 0.33, 0.335, 0.34, 0.345, 0.35, 0.355, 0.36, 0.365, &
        !   0.37, 0.375, 0.38, 0.385, 0.39, 0.395, 0.4, 0.405, 0.41, 0.415, 0.42, 0.425, 0.43, 0.435, 0.44, &
        !   0.445, 0.45, 0.455, 0.46, 0.465, 0.47, 0.475, 0.48, 0.485, 0.49, 0.495, 0.5, 0.505, 0.51, 0.515, &
        !    0.52, 0.525, 0.53, 0.535, 0.54, 0.545, 0.55, 0.555, 0.56, 0.565, 0.57, 0.575, 0.58, 0.585, 0.59, &
        !     0.595, 0.6, 0.605, 0.61, 0.615, 0.62, 0.625, 0.63, 0.635, 0.64, 0.645, 0.65, 0.655, 0.66, 0.665, &
        !      0.67, 0.675, 0.68, 0.685, 0.69, 0.695, 0.7, 0.705, 0.71, 0.715, 0.72, 0.725, 0.73, 0.735, 0.74, &
        !       0.745, 0.75, 0.755, 0.76, 0.765, 0.77, 0.775, 0.78, 0.785, 0.79, 0.795, 0.8, 0.805, 0.81, 0.815, &
        !        0.82, 0.825, 0.83, 0.835, 0.84, 0.845, 0.85, 0.855, 0.86, 0.865, 0.87, 0.875, 0.88, 0.885, 0.89, &
        !         0.895, 0.9, 0.905, 0.91, 0.915, 0.92, 0.925, 0.93, 0.935, 0.94, 0.945, 0.95, 0.955, 0.96, 0.965, &
        !          0.97, 0.975, 0.98, 0.985, 0.99, 0.995, 1. /)  ! Corresponding values
        !


 cdf(:, 1) = (/ 0.000063, 0.000697, 0.002386, 0.005284, 0.009364, 0.014526, 0.020643, 0.027587, 0.035233, 0.043472, &
                0.052204, 0.061341, 0.070807, 0.080537, 0.090473, 0.100566, 0.110773, 0.121057, 0.131386, 0.141734, &
                0.152076, 0.162392, 0.172666, 0.182882, 0.193028, 0.203093, 0.213067, 0.222944, 0.232717, 0.24238, &
                0.25193, 0.261362, 0.270673, 0.279863, 0.288928, 0.297869, 0.306683, 0.315371, 0.323932, 0.332368, &
                0.340678, 0.348863, 0.356924, 0.364861, 0.372677, 0.380372, 0.387948, 0.395405, 0.402747, 0.409973, &
                0.417086, 0.424088, 0.430979, 0.437762, 0.444439, 0.45101, 0.457479, 0.463846, 0.470113, 0.476282, &
                0.482354, 0.488332, 0.494217, 0.50001, 0.505714, 0.511329, 0.516858, 0.522301, 0.527661, 0.532939, &
                0.538137, 0.543255, 0.548296, 0.553261, 0.55815, 0.562967, 0.567711, 0.572385, 0.576989, 0.581524, &
                0.585993, 0.590396, 0.594735, 0.59901, 0.603223, 0.607375, 0.611467, 0.6155, 0.619476, 0.623394, &
                0.627257, 0.631065, 0.634819, 0.63852, 0.64217, 0.645768, 0.649317, 0.652816, 0.656267, 0.659671, &
                0.663028, 0.666339, 0.669605, 0.672826, 0.676005, 0.67914, 0.682234, 0.685286, 0.688298, 0.69127, &
                0.694203, 0.697097, 0.699953, 0.702772, 0.705555, 0.708301, 0.711012, 0.713688, 0.71633, 0.718938, &
                0.721513, 0.724056, 0.726566, 0.729045, 0.731493, 0.73391, 0.736297, 0.738655, 0.740983, 0.743283, &
                0.745555, 0.747799, 0.750016, 0.752206, 0.754369, 0.756506, 0.758618, 0.760704, 0.762766, 0.764803, &
                0.766816, 0.768806, 0.770772, 0.772714, 0.774635, 0.776533, 0.778409, 0.780263, 0.782096, 0.783908, &
                0.785699, 0.78747, 0.789221, 0.790952, 0.792663, 0.794355, 0.796029, 0.797683, 0.799319, 0.800937, &
                0.802537, 0.80412, 0.805685, 0.807233, 0.808764, 0.810278, 0.811776, 0.813258, 0.814723, 0.816173, &
                0.817608, 0.819027, 0.820431, 0.82182, 0.823195, 0.824555, 0.8259, 0.827232, 0.828549, 0.829853, &
                0.831144, 0.832421, 0.833685, 0.834935, 0.836173, 0.837399, 0.838611, 0.839812, 0.841, 0.842176, &
                0.843341, 0.844494, 0.845635, 0.846765, 0.847883, 0.84899, 0.850087, 0.851172, 0.852247, 0.853312, &
                0.854366, 0.855409, 0.856443, 0.857466, 0.858479, 0.859483, 0.860477, 0.861462, 0.862437, 0.863402, &
                0.864359, 0.865306, 0.866245, 0.867174, 0.868095, 0.869007, 0.869911, 0.870806, 0.871693, 0.872572, &
                0.873442, 0.874305, 0.875159, 0.876006, 0.876844, 0.877675, 0.878499, 0.879315, 0.880124, 0.880925, &
                0.881719, 0.882506, 0.883286, 0.884058, 0.884824, 0.885583, 0.886336, 0.887081, 0.88782, 0.888553, &
                0.889279, 0.889999, 0.890712, 0.891419, 0.89212, 0.892815, 0.893504, 0.894187, 0.894864, 0.895535, &
                0.8962, 0.89686, 0.897514, 0.898162, 0.898805, 0.899443, 0.900075, 0.900701, 0.901323, 0.901939, &
                0.90255, 0.903156, 0.903757, 0.904353, 0.904943, 0.905529, 0.90611, 0.906687, 0.907258, 0.907825, &
                0.908387, 0.908945, 0.909498, 0.910047, 0.910591, 0.91113, 0.911666, 0.912197, 0.912723, 0.913246, &
                0.913764, 0.914278, 0.914788, 0.915294, 0.915796, 0.916294, 0.916788, 0.917278, 0.917765, 0.918247, &
                0.918726, 0.919201, 0.919672, 0.920139, 0.920603, 0.921063, 0.92152, 0.921973, 0.922423, 0.922869, &
                0.923312, 0.923751, 0.924187, 0.92462, 0.925049, 0.925476, 0.925898, 0.926318, 0.926735, 0.927148, &
                0.927558, 0.927966, 0.92837, 0.928771, 0.929169, 0.929564, 0.929956, 0.930346, 0.930732, 0.931116, &
                0.931496, 0.931874, 0.93225, 0.932622, 0.932992, 0.933359, 0.933723, 0.934085, 0.934444, 0.9348, &
                0.935154, 0.935505, 0.935854, 0.9362, 0.936544, 0.936886, 0.937224, 0.937561, 0.937895, 0.938227, &
                0.938556, 0.938883, 0.939208, 0.93953, 0.93985, 0.940168, 0.940484, 0.940797, 0.941109, 0.941418, &
                0.941725, 0.942029, 0.942332, 0.942633, 0.942931, 0.943228, 0.943522, 0.943815, 0.944105, 0.944393, &
                0.94468, 0.944964, 0.945247, 0.945528, 0.945806, 0.946083, 0.946358, 0.946631, 0.946902, 0.947172, &
                0.94744, 0.947705, 0.947969, 0.948232, 0.948492, 0.948751, 0.949008, 0.949264, 0.949517, 0.949769, &
                0.95002, 0.950268, 0.950516, 0.950761, 0.951005, 0.951247, 0.951488, 0.951727, 0.951964, 0.952201, &
                0.952435, 0.952668, 0.952899, 0.953129, 0.953358, 0.953585, 0.95381, 0.954035, 0.954257, 0.954478, &
                0.954698, 0.954917, 0.955134, 0.95535, 0.955564, 0.955777, 0.955988, 0.956199, 0.956408, 0.956615, &
                0.956822, 0.957027, 0.957231, 0.957433, 0.957634, 0.957834, 0.958033, 0.95823, 0.958427, 0.958622, &
                0.958816, 0.959008, 0.9592, 0.95939, 0.959579, 0.959767, 0.959954, 0.96014, 0.960324, 0.960508, &
                0.96069, 0.960872, 0.961052, 0.961231, 0.961409, 0.961586, 0.961761, 0.961936, 0.96211, 0.962283, &
                0.962454, 0.962625, 0.962794, 0.962963, 0.963131, 0.963297, 0.963463, 0.963627, 0.963791, 0.963954, &
                0.964116, 0.964276, 0.964436, 0.964595, 0.964753, 0.96491, 0.965066, 0.965221, 0.965376, 0.965529, &
                0.965682, 0.965833, 0.965984, 0.966134, 0.966283, 0.966431, 0.966579, 0.966725, 0.966871, 0.967015, &
                0.967159, 0.967303, 0.967445, 0.967586, 0.967727, 0.967867, 0.968006, 0.968145, 0.968282, 0.968419, &
                0.968555, 0.96869, 0.968825, 0.968958, 0.969091, 0.969224, 0.969355, 0.969486, 0.969616, 0.969745, &
                0.969874, 0.970002, 0.970129, 0.970255, 0.970381, 0.970506, 0.97063, 0.970754, 0.970877, 0.971, &
                0.971121, 0.971242, 0.971363, 0.971482, 0.971601, 0.97172, 0.971837, 0.971955, 0.972071, 0.972187, &
                0.972302, 0.972417, 0.972531, 0.972644, 0.972757, 0.972869, 0.97298, 0.973091, 0.973202, 0.973311, &
                0.973421, 0.973529, 0.973637, 0.973745, 0.973852, 0.973958, 0.974064, 0.974169, 0.974273, 0.974378, &
                0.974481, 0.974584, 0.974686, 0.974788, 0.97489, 0.974991, 0.975091, 0.975191, 0.97529, 0.975389, &
                0.975487, 0.975585, 0.975682, 0.975779, 0.975875, 0.975971, 0.976066, 0.976161, 0.976255, 0.976349, &
                0.976442, 0.976535, 0.976627, 0.976719, 0.97681, 0.976901, 0.976992, 0.977082, 0.977171, 0.97726, &
                0.977349, 0.977437, 0.977525, 0.977612, 0.977699, 0.977786, 0.977872, 0.977957, 0.978042, 0.978127, &
                0.978211, 0.978295, 0.978379, 0.978462, 0.978544, 0.978626, 0.978708, 0.97879, 0.978871, 0.978951, &
                0.979031, 0.979111, 0.979191, 0.97927, 0.979348, 0.979426, 0.979504, 0.979582, 0.979659, 0.979736, &
                0.979812, 0.979888, 0.979963, 0.980039, 0.980114, 0.980188, 0.980262, 0.980336, 0.980409, 0.980482, &
                0.980555, 0.980628, 0.9807, 0.980771, 0.980843, 0.980914, 0.980984, 0.981054, 0.981124, 0.981194, &
                0.981263, 0.981332, 0.981401, 0.981469, 0.981537, 0.981605, 0.981672, 0.981739, 0.981806, 0.981872, &
                0.981938, 0.982004, 0.98207, 0.982135, 0.9822, 0.982264, 0.982328, 0.982392, 0.982456, 0.982519, &
                0.982582, 0.982645, 0.982708, 0.98277, 0.982832, 0.982893, 0.982955, 0.983016, 0.983076, 0.983137, &
                0.983197, 0.983257, 0.983317, 0.983376, 0.983435, 0.983494, 0.983552, 0.983611, 0.983669, 0.983726, &
                0.983784, 0.983841, 0.983898, 0.983955, 0.984011, 0.984068, 0.984123, 0.984179, 0.984235, 0.98429, &
                0.984345, 0.984399, 0.984454, 0.984508, 0.984562, 0.984616, 0.984669, 0.984722, 0.984775, 0.984828, &
                0.984881, 0.984933, 0.984985, 0.985037, 0.985088, 0.98514, 0.985191, 0.985242, 0.985292, 0.985343, &
                0.985393, 0.985443, 0.985493, 0.985542, 0.985592, 0.985641, 0.98569, 0.985739, 0.985787, 0.985835, &
                0.985883, 0.985931, 0.985979, 0.986026, 0.986074, 0.986121, 0.986168, 0.986214, 0.986261, 0.986307, &
                0.986353, 0.986399, 0.986444, 0.98649, 0.986535, 0.98658, 0.986625, 0.986669, 0.986714, 0.986758, &
                0.986802, 0.986846, 0.98689, 0.986933, 0.986977, 0.98702, 0.987063, 0.987106, 0.987148, 0.987191, &
                0.987233, 0.987275, 0.987317, 0.987359, 0.9874, 0.987441, 0.987483, 0.987524, 0.987564, 0.987605, &
                0.987646, 0.987686, 0.987726, 0.987766, 0.987806, 0.987846, 0.987885, 0.987924, 0.987963, 0.988002, &
                0.988041, 0.98808, 0.988118, 0.988157, 0.988195, 0.988233, 0.988271, 0.988309, 0.988346, 0.988384, &
                0.988421, 0.988458, 0.988495, 0.988532, 0.988568, 0.988605, 0.988641, 0.988677, 0.988714, 0.988749, &
                0.988785, 0.988821, 0.988856, 0.988892, 0.988927, 0.988962, 0.988997, 0.989031, 0.989066, 0.989101, &
                0.989135, 0.989169, 0.989203, 0.989237, 0.989271, 0.989305, 0.989338, 0.989371, 0.989405, 0.989438, &
                0.989471, 0.989504, 0.989536, 0.989569, 0.989601, 0.989634, 0.989666, 0.989698, 0.98973, 0.989762, &
                0.989793, 0.989825, 0.989856, 0.989888, 0.989919, 0.98995, 0.989981, 0.990012, 0.990043, 0.990073, &
                0.990132, 0.99019, 0.990249, 0.990307, 0.990365, 0.990423, 0.990481, 0.990538, 0.990596, 0.990654, &
                0.990711, 0.990768, 0.990826, 0.990883, 0.99094, 0.990996, 0.991053, 0.99111, 0.991166, 0.991223, &
                0.991279, 0.991335, 0.991391, 0.991447, 0.991503, 0.991559, 0.991615, 0.99167, 0.991726, 0.991781, &
                0.991836, 0.991891, 0.991946, 0.992001, 0.992056, 0.992111, 0.992166, 0.99222, 0.992274, 0.992329, &
                0.992383, 0.992437, 0.992491, 0.992545, 0.992599, 0.992653, 0.992706, 0.99276, 0.992813, 0.992867, &
                0.99292, 0.992973, 0.993026, 0.993079, 0.993132, 0.993185, 0.993238, 0.99329, 0.993343, 0.993395, &
                0.993447, 0.9935, 0.993552, 0.993604, 0.993656, 0.993708, 0.99376, 0.993811, 0.993863, 0.993914, &
                0.993966, 0.994017, 0.994068, 0.99412, 0.994171, 0.994222, 0.994273, 0.994324, 0.994374, 0.994425, &
                0.994476, 0.994526, 0.994577, 0.994627, 0.994677, 0.994727, 0.994778, 0.994828, 0.994878, 0.994927, &
                0.994977, 0.995027, 0.995077, 0.995126, 0.995176, 0.995225, 0.995274, 0.995324, 0.995373, 0.995422, &
                0.995471, 0.99552, 0.995569, 0.995618, 0.995666, 0.995715, 0.995764, 0.995812, 0.995861, 0.995909, &
                0.995957, 0.996005, 0.996054, 0.996102, 0.99615, 0.996198, 0.996245, 0.996293, 0.996341, 0.996389, &
                0.996436, 0.996484, 0.996531, 0.996579, 0.996626, 0.996673, 0.99672, 0.996767, 0.996814, 0.996861, &
                0.996908, 0.996955, 0.997002, 0.997049, 0.997095, 0.997142, 0.997188, 0.997235, 0.997281, 0.997328, &
                0.997374, 0.99742, 0.997466, 0.997512, 0.997558, 0.997604, 0.99765, 0.997696, 0.997742, 0.997787, &
                0.997833, 0.997879, 0.997924, 0.99797, 0.998015, 0.99806, 0.998106, 0.998151, 0.998196, 0.998241, &
                0.998286, 0.998331, 0.998376, 0.998421, 0.998466, 0.998511, 0.998555, 0.9986, 0.998644, 0.998689, &
                0.998734, 0.998778, 0.998822, 0.998867, 0.998911, 0.998955, 0.998999, 0.999043, 0.999087, 0.999131, &
                0.999175, 0.999219, 0.999263, 0.999307, 0.99935, 0.999394, 0.999438, 0.999481, 0.999525, 0.999568, &
                0.999612, 0.999655, 0.999698, 0.999742, 0.999785, 0.999828, 0.999871, 0.999914, 0.999957, 1. /)  ! Cumulative probabilities
 cdf(:, 2) = (/ 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01, &
                0.011, 0.012, 0.013, 0.014, 0.015, 0.016, 0.017, 0.018, 0.019, 0.02, &
                0.021, 0.022, 0.023, 0.024, 0.025, 0.026, 0.027, 0.028, 0.029, 0.03, &
                0.031, 0.032, 0.033, 0.034, 0.035, 0.036, 0.037, 0.038, 0.039, 0.04, &
                0.041, 0.042, 0.043, 0.044, 0.045, 0.046, 0.047, 0.048, 0.049, 0.05, &
                0.051, 0.052, 0.053, 0.054, 0.055, 0.056, 0.057, 0.058, 0.059, 0.06, &
                0.061, 0.062, 0.063, 0.064, 0.065, 0.066, 0.067, 0.068, 0.069, 0.07, &
                0.071, 0.072, 0.073, 0.074, 0.075, 0.076, 0.077, 0.078, 0.079, 0.08, &
                0.081, 0.082, 0.083, 0.084, 0.085, 0.086, 0.087, 0.088, 0.089, 0.09, &
                0.091, 0.092, 0.093, 0.094, 0.095, 0.096, 0.097, 0.098, 0.099, 0.1, &
                0.101, 0.102, 0.103, 0.104, 0.105, 0.106, 0.107, 0.108, 0.109, 0.11, &
                0.111, 0.112, 0.113, 0.114, 0.115, 0.116, 0.117, 0.118, 0.119, 0.12, &
                0.121, 0.122, 0.123, 0.124, 0.125, 0.126, 0.127, 0.128, 0.129, 0.13, &
                0.131, 0.132, 0.133, 0.134, 0.135, 0.136, 0.137, 0.138, 0.139, 0.14, &
                0.141, 0.142, 0.143, 0.144, 0.145, 0.146, 0.147, 0.148, 0.149, 0.15, &
                0.151, 0.152, 0.153, 0.154, 0.155, 0.156, 0.157, 0.158, 0.159, 0.16, &
                0.161, 0.162, 0.163, 0.164, 0.165, 0.166, 0.167, 0.168, 0.169, 0.17, &
                0.171, 0.172, 0.173, 0.174, 0.175, 0.176, 0.177, 0.178, 0.179, 0.18, &
                0.181, 0.182, 0.183, 0.184, 0.185, 0.186, 0.187, 0.188, 0.189, 0.19, &
                0.191, 0.192, 0.193, 0.194, 0.195, 0.196, 0.197, 0.198, 0.199, 0.2, &
                0.201, 0.202, 0.203, 0.204, 0.205, 0.206, 0.207, 0.208, 0.209, 0.21, &
                0.211, 0.212, 0.213, 0.214, 0.215, 0.216, 0.217, 0.218, 0.219, 0.22, &
                0.221, 0.222, 0.223, 0.224, 0.225, 0.226, 0.227, 0.228, 0.229, 0.23, &
                0.231, 0.232, 0.233, 0.234, 0.235, 0.236, 0.237, 0.238, 0.239, 0.24, &
                0.241, 0.242, 0.243, 0.244, 0.245, 0.246, 0.247, 0.248, 0.249, 0.25, &
                0.251, 0.252, 0.253, 0.254, 0.255, 0.256, 0.257, 0.258, 0.259, 0.26, &
                0.261, 0.262, 0.263, 0.264, 0.265, 0.266, 0.267, 0.268, 0.269, 0.27, &
                0.271, 0.272, 0.273, 0.274, 0.275, 0.276, 0.277, 0.278, 0.279, 0.28, &
                0.281, 0.282, 0.283, 0.284, 0.285, 0.286, 0.287, 0.288, 0.289, 0.29, &
                0.291, 0.292, 0.293, 0.294, 0.295, 0.296, 0.297, 0.298, 0.299, 0.3, &
                0.301, 0.302, 0.303, 0.304, 0.305, 0.306, 0.307, 0.308, 0.309, 0.31, &
                0.311, 0.312, 0.313, 0.314, 0.315, 0.316, 0.317, 0.318, 0.319, 0.32, &
                0.321, 0.322, 0.323, 0.324, 0.325, 0.326, 0.327, 0.328, 0.329, 0.33, &
                0.331, 0.332, 0.333, 0.334, 0.335, 0.336, 0.337, 0.338, 0.339, 0.34, &
                0.341, 0.342, 0.343, 0.344, 0.345, 0.346, 0.347, 0.348, 0.349, 0.35, &
                0.351, 0.352, 0.353, 0.354, 0.355, 0.356, 0.357, 0.358, 0.359, 0.36, &
                0.361, 0.362, 0.363, 0.364, 0.365, 0.366, 0.367, 0.368, 0.369, 0.37, &
                0.371, 0.372, 0.373, 0.374, 0.375, 0.376, 0.377, 0.378, 0.379, 0.38, &
                0.381, 0.382, 0.383, 0.384, 0.385, 0.386, 0.387, 0.388, 0.389, 0.39, &
                0.391, 0.392, 0.393, 0.394, 0.395, 0.396, 0.397, 0.398, 0.399, 0.4, &
                0.401, 0.402, 0.403, 0.404, 0.405, 0.406, 0.407, 0.408, 0.409, 0.41, &
                0.411, 0.412, 0.413, 0.414, 0.415, 0.416, 0.417, 0.418, 0.419, 0.42, &
                0.421, 0.422, 0.423, 0.424, 0.425, 0.426, 0.427, 0.428, 0.429, 0.43, &
                0.431, 0.432, 0.433, 0.434, 0.435, 0.436, 0.437, 0.438, 0.439, 0.44, &
                0.441, 0.442, 0.443, 0.444, 0.445, 0.446, 0.447, 0.448, 0.449, 0.45, &
                0.451, 0.452, 0.453, 0.454, 0.455, 0.456, 0.457, 0.458, 0.459, 0.46, &
                0.461, 0.462, 0.463, 0.464, 0.465, 0.466, 0.467, 0.468, 0.469, 0.47, &
                0.471, 0.472, 0.473, 0.474, 0.475, 0.476, 0.477, 0.478, 0.479, 0.48, &
                0.481, 0.482, 0.483, 0.484, 0.485, 0.486, 0.487, 0.488, 0.489, 0.49, &
                0.491, 0.492, 0.493, 0.494, 0.495, 0.496, 0.497, 0.498, 0.499, 0.5, &
                0.501, 0.502, 0.503, 0.504, 0.505, 0.506, 0.507, 0.508, 0.509, 0.51, &
                0.511, 0.512, 0.513, 0.514, 0.515, 0.516, 0.517, 0.518, 0.519, 0.52, &
                0.521, 0.522, 0.523, 0.524, 0.525, 0.526, 0.527, 0.528, 0.529, 0.53, &
                0.531, 0.532, 0.533, 0.534, 0.535, 0.536, 0.537, 0.538, 0.539, 0.54, &
                0.541, 0.542, 0.543, 0.544, 0.545, 0.546, 0.547, 0.548, 0.549, 0.55, &
                0.551, 0.552, 0.553, 0.554, 0.555, 0.556, 0.557, 0.558, 0.559, 0.56, &
                0.561, 0.562, 0.563, 0.564, 0.565, 0.566, 0.567, 0.568, 0.569, 0.57, &
                0.571, 0.572, 0.573, 0.574, 0.575, 0.576, 0.577, 0.578, 0.579, 0.58, &
                0.581, 0.582, 0.583, 0.584, 0.585, 0.586, 0.587, 0.588, 0.589, 0.59, &
                0.591, 0.592, 0.593, 0.594, 0.595, 0.596, 0.597, 0.598, 0.599, 0.6, &
                0.601, 0.602, 0.603, 0.604, 0.605, 0.606, 0.607, 0.608, 0.609, 0.61, &
                0.611, 0.612, 0.613, 0.614, 0.615, 0.616, 0.617, 0.618, 0.619, 0.62, &
                0.621, 0.622, 0.623, 0.624, 0.625, 0.626, 0.627, 0.628, 0.629, 0.63, &
                0.631, 0.632, 0.633, 0.634, 0.635, 0.636, 0.637, 0.638, 0.639, 0.64, &
                0.641, 0.642, 0.643, 0.644, 0.645, 0.646, 0.647, 0.648, 0.649, 0.65, &
                0.651, 0.652, 0.653, 0.654, 0.655, 0.656, 0.657, 0.658, 0.659, 0.66, &
                0.661, 0.662, 0.663, 0.664, 0.665, 0.666, 0.667, 0.668, 0.669, 0.67, &
                0.671, 0.672, 0.673, 0.674, 0.675, 0.676, 0.677, 0.678, 0.679, 0.68, &
                0.681, 0.682, 0.683, 0.684, 0.685, 0.686, 0.687, 0.688, 0.689, 0.69, &
                0.691, 0.692, 0.693, 0.694, 0.695, 0.696, 0.697, 0.698, 0.699, 0.7, &
                0.701, 0.702, 0.703, 0.704, 0.705, 0.706, 0.707, 0.708, 0.709, 0.71, &
                0.711, 0.712, 0.713, 0.714, 0.715, 0.716, 0.717, 0.718, 0.719, 0.72, &
                0.721, 0.722, 0.723, 0.724, 0.725, 0.726, 0.727, 0.728, 0.729, 0.73, &
                0.731, 0.732, 0.733, 0.734, 0.735, 0.736, 0.737, 0.738, 0.739, 0.74, &
                0.741, 0.742, 0.743, 0.744, 0.745, 0.746, 0.747, 0.748, 0.749, 0.75, &
                0.751, 0.752, 0.753, 0.754, 0.755, 0.756, 0.757, 0.758, 0.759, 0.76, &
                0.761, 0.762, 0.763, 0.764, 0.765, 0.766, 0.767, 0.768, 0.769, 0.77, &
                0.771, 0.772, 0.773, 0.774, 0.775, 0.776, 0.777, 0.778, 0.779, 0.78, &
                0.781, 0.782, 0.783, 0.784, 0.785, 0.786, 0.787, 0.788, 0.789, 0.79, &
                0.791, 0.792, 0.793, 0.794, 0.795, 0.796, 0.797, 0.798, 0.799, 0.8, &
                0.801, 0.802, 0.803, 0.804, 0.805, 0.806, 0.807, 0.808, 0.809, 0.81, &
                0.811, 0.812, 0.813, 0.814, 0.815, 0.816, 0.817, 0.818, 0.819, 0.82, &
                0.821, 0.822, 0.823, 0.824, 0.825, 0.826, 0.827, 0.828, 0.829, 0.83, &
                0.831, 0.832, 0.833, 0.834, 0.835, 0.836, 0.837, 0.838, 0.839, 0.84, &
                0.841, 0.842, 0.843, 0.844, 0.845, 0.846, 0.847, 0.848, 0.849, 0.85, &
                0.851, 0.852, 0.853, 0.854, 0.855, 0.856, 0.857, 0.858, 0.859, 0.86, &
                0.861, 0.862, 0.863, 0.864, 0.865, 0.866, 0.867, 0.868, 0.869, 0.87, &
                0.871, 0.872, 0.873, 0.874, 0.875, 0.876, 0.877, 0.878, 0.879, 0.88, &
                0.881, 0.882, 0.883, 0.884, 0.885, 0.886, 0.887, 0.888, 0.889, 0.89, &
                0.891, 0.892, 0.893, 0.894, 0.895, 0.896, 0.897, 0.898, 0.899, 0.9, &
                0.901, 0.902, 0.903, 0.904, 0.905, 0.906, 0.907, 0.908, 0.909, 0.91, &
                0.911, 0.912, 0.913, 0.914, 0.915, 0.916, 0.917, 0.918, 0.919, 0.92, &
                0.921, 0.922, 0.923, 0.924, 0.925, 0.926, 0.927, 0.928, 0.929, 0.93, &
                0.931, 0.932, 0.933, 0.934, 0.935, 0.936, 0.937, 0.938, 0.939, 0.94, &
                0.941, 0.942, 0.943, 0.944, 0.945, 0.946, 0.947, 0.948, 0.949, 0.95, &
                0.951, 0.952, 0.953, 0.954, 0.955, 0.956, 0.957, 0.958, 0.959, 0.96, &
                0.961, 0.962, 0.963, 0.964, 0.965, 0.966, 0.967, 0.968, 0.969, 0.97, &
                0.971, 0.972, 0.973, 0.974, 0.975, 0.976, 0.977, 0.978, 0.979, 0.98, &
                0.981, 0.982, 0.983, 0.984, 0.985, 0.986, 0.987, 0.988, 0.989, 0.99, &
                0.991, 0.992, 0.993, 0.994, 0.995, 0.996, 0.997, 0.998, 0.999, 1. /)  ! Corresponding values





        ! Generate a random number in the range [0, 1]
        call random_seed()              ! Initialize the random number generator
        call random_number(rndm)

        ! Find the corresponding value in the CDF matrix
        rdvol_sampled = cdf(1, 2)  ! Default value in case of numerical errors
        do i = 1, n
          if (rndm <= cdf(i, 1)) then
            rdvol_sampled = cdf(i, 2)
            exit
          end if
        end do

        ! Output the result
        ! print *, "Random number: ", random_number
        ! print *, "Sampled value: ", sampled_value

      end subroutine sample_rdvol

