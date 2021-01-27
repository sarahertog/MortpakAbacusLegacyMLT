require(DemoTools)

lt_model_lq_single <-  function(
    Sex, # has to be specified always
    fitted_logquad = NULL, 
    q0_5 = NULL, 
    q0_1 = NULL, 
    q15_45 = NULL, 
    q15_35 = NULL, 
    e0 = NULL, 
    radix = 1e5, 
    tol = 1e-9, 
    maxit = 200, 
    axmethod = "pas", 
    a0rule = "ak", 
    IMR = NA, 
    region = "w", 
    mod = TRUE,
    SRB = 1.05,
    OAnew = 110,
    extrapLaw = "kannisto",
    extrapFit = 70:90) {
  
  # abridged LogQuad life table
  
  lt_abr <- lt_model_lq(Sex = Sex,
                        fitted_logquad = fitted_logquad, 
                        q0_5 = q0_5, 
                        q0_1 = q0_1, 
                        q15_45 = q15_45, 
                        q15_35 = q15_35, 
                        e0 = e0, 
                        radix = radix, 
                        tol = tol, 
                        maxit = maxit, 
                        axmethod = axmethod, 
                        a0rule = a0rule, 
                        IMR = IMR, 
                        region = region, 
                        mod = mod,
                        SRB = SRB)
  
  # graduate to single year of age truncate to open age 90+ bc LogQuad results at ages 95+ imply qx>1 and get automatically
  # overwritten by the lt_abridged function in lt_abridged2single, with undesirable results (discontinuity) for the mortality
  # pattern above age 90
  lt_sgl <- lt_abridged2single(nMx = lt_abr$lt$nMx, 
                               Age = lt_abr$lt$Age,
                               axmethod = axmethod,
                               a0rule = a0rule,
                               Sex = Sex,
                               region = region,
                               IMR = IMR,
                               mod = mod,
                               SRB = SRB, 
                               OAG = TRUE,
                               OAnew = 90,
                               extrapLaw = extrapLaw)
  
  # if user-specified OAnew > 90, then extend the the single age life table 
  if (OAnew > 90) {
  lt_sgl <- lt_single_mx(nMx = lt_sgl$nMx,
                         radix = radix,
                         a0rule = a0rule,
                         Sex = Sex,
                         region = region,
                         IMR = IMR, 
                         mod = mod,
                         SRB = SRB, 
                         OAG = TRUE,
                         OAnew = OAnew,
                         extrapLaw = extrapLaw,
                         extrapFit = extrapFit)
  }
  
  return(lt_sgl)
  
}