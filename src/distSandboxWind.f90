subroutine windrisk(siteInfoDist, spec, h, openedge, sitetype, tsum, tsincethin, &
  wrisk5dd1, wrisk5dd2, wrisk5dd3, wrisk0, wrisk5, wrisk)
  IMPLICIT NONE
  REAL (kind=8), intent(inout) ::  siteInfoDist(4) ! 5-year wind risk (suvanto output), pre-logit value, annual risk
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









