
#library(devtools)
#install_github("PPgp/HIV.LifeTables")
library(HIV.LifeTables)

library(DemoTools)


lt_model_hiv_single <- function (q5 = NULL,
                                 q1545 = NULL,
                                 e0 = NULL, 
                                 prev = NULL, 
                                 child.art = NULL, 
                                 adult.art = NULL, 
                                 hiv_region = 1, # 1 is Africa; 0 for all others
                                 opt = TRUE,  # vary intercept to approximate input 45q15
                                 recal = NULL, 
                                 axmethod = "un", 
                                 a0rule = "ak", 
                                 Sex = "m", 
                                 region = "w", 
                                 IMR = NA, 
                                 mod = TRUE, 
                                 SRB = 1.05, 
                                 OAnew = max(Age), 
                                 extrapLaw = "kannisto", 
                                 extrapFrom = NA, 
                                 extrapFit = NA) {
  
  stopifnot(Sex %in% c("m", "f"))
  hiv_sex <- ifelse(Sex == "f", 1, 0)
  
  q5_flag <- !is.null(q5)
  q1545_flag <- !is.null(q1545)
  e0_flag <- !is.null(e0)
  
  # indentify the valid inputs
  model_type <- NULL
  if (q5_flag & !q1545_flag & !e0_flag) { model_type = "5q0" }
  if (q5_flag & q1545_flag & !e0_flag) { model_type = "45q15" }
  if (!q5_flag & !q1545_flag & e0_flag) { model_type = "e0" }
  stopifnot(!is.null(model_type))
  
  ## Call HIV.LifeTables functions from David Sharrow's package
  # Abacus purported to be calling adult.art and child.art in addition to prevalence
  # but the function in Cran HIV.LifeTables do not reference these parameters so I have excluded them
  
  if (model_type == "5q0") { # child mortality only
    mx.out <- HIV.LifeTables::mortmod.5q0(child.mort = q5, 
                                          prev = prev, 
                                          child.art = child.art,
                                          adult.art = adult.art,
                                          region = hiv_region, 
                                          sex = hiv_sex, 
                                          opt = opt,
                                          recal = recal)
  }
  if (model_type == "45q15") { # child mortality and adult mortality
    mx.out <- HIV.LifeTables::mortmod.45q15(child.mort = q5, 
                                            adult.mort = q1545,
                                            prev = prev,  
                                            child.art = child.art,
                                            adult.art = adult.art,
                                            region = hiv_region, 
                                            sex = hiv_sex, 
                                            opt = opt,
                                            recal = recal)
  }
  if (model_type == "e0") { # e0 only
    mx.out <- HIV.LifeTables::mortmod.e0(e0 = e0, 
                                         prev = prev,  
                                         child.art = child.art,
                                         adult.art = adult.art,
                                         region = hiv_region, 
                                         sex = hiv_sex, 
                                         opt = opt,
                                         recal = recal)
  }
  
  
  # graduate the abridged mx pattern to single years of age
  Age <- c(0,1,seq(5,100,5))
  LT <- lt_abridged2single(nMx = mx.out,
                           Age = Age,
                           radix = 1e+05, 
                           axmethod = "un", 
                           a0rule = "ak", 
                           Sex = "m", 
                           region = "w", 
                           IMR = NA, 
                           mod = TRUE, 
                           SRB = 1.05, 
                           OAG = TRUE, 
                           OAnew = OAnew, 
                           extrapLaw = "kannisto", 
                           extrapFrom = extrapFrom, 
                           extrapFit = extrapFit)
  
  return(LT)
  
}
  
  

