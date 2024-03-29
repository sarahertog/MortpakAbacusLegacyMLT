# --------------------------------------------------------------------------------------------------------------------------
#      COMBIN is derived from MORTPAK software package and customized for Abacus
#               UNITED NATION SOFTWARE PACKAGE FOR MORTALITY MEASUREMENT   
  # Description of procedure
  # Get MLT by matching on q(15,n) to get set of q's
  # If Coale-Demeny models are used, use above q as a user defined model and use type = "user_defined" in BESTFT
  # For first 2 age groups convert to l1 and l5.  
  # If l1 is not available, use MLT and match on l5 to get l1.
  # When both l1 and l5 are available, use BESTFT first with 1q0 and 5q20 then with 4q1 and 5q20. Average
  # both results to get 5q5 and 5q10.  Note that 5q20 above is from MLT.
  # If l5 is not available, use BESTFT with 1q0 and 5q20 to get 4q1, 5q5 and 5q10. 
  # Note that early age groups come from second component fit.
# Sara Hertog modified to maintain method, but simplify syntax for ccmppWPP 2021 revision
# --------------------------------------------------------------------------------------------------------------------------

#' Estimate UN or Coale Demeney family model life tables, matching multiple input parameters
#' 
#' @details xxx
#' @param Sex Choose the sex of the population. 
#' #' The following options are available: \itemize{
#'   \item{\code{"b"}} -- Both sex; 
#'   \item{\code{"f"}} -- Females;
#'   \item{\code{"m"}} -- Males.
#'   }
#' @param type Choose the family of model life tables
#'  #' The following options are available: \itemize{
#'   \item{\code{"CD_East"}} -- Coale-Demeny East; 
#'   \item{\code{"CD_North"}} -- Coale-Demeny North;
#'   \item{\code{"CD_South"}} -- Coale-Demeny South;
#'   \item{\code{"CD_West"}} -- Coale-Demeny West;
#'   \item{\code{"UN_Chilean"}} -- UN Chilean;
#'   \item{\code{"UN_Far_Eastern"}} -- UN Far Eastern;
#'   \item{\code{"UN_General"}} -- UN General;
#'   \item{\code{"UN_Latin_American"}} -- UN Latin American;
#'   \item{\code{"UN_South_Asian"}} -- UN South Asian.
#'   }
#' @param q1 numeric. Probability of dying between birth and age 1.
#' @param q5 numeric. Probability of dying between birth and age 5.
#' @param indicator_adult character. Adult mortality life table indicator to match on
#'  #' The following options are available: \itemize{
#'   \item{\code{"35q15"}} -- Probability of dying between age 15 and age 50;
#'   \item{\code{"45q15"}} -- Probability of dying between age 15 and age 60.
#' @param value_adult numeric. Value of the adult mortality indicator to match on;
#' @inheritParams lt_abridged
#' @return data.frame. with two columns: age, giving abridged age groups, and mx
#' with model age specific mortality rates for abridged age groups
#' @importFrom stats uniroot MortCast DemoTools
#' @examples 
#' 
#'
#'lt <- lt_model_cdun_combin_single(type = "CD_West", Sex = "f",
#'                                  q1 = .04651, q5 = .06027,
#'                                  indicator_adult = "45q15", value_adult = .241,
#'                                  OAnew = 110)
#'



lt_model_cdun_combin_single <- function(type, 
                                        Sex, 
                                        q1 = NA, 
                                        q5 = NA, 
                                        indicator_adult, 
                                        value_adult,
                                        OAnew = 110,
                                        radix = 1e+05, 
                                        axmethod = "un", 
                                        a0rule = "ak", 
                                        region = "w", 
                                        IMR = NA, 
                                        mod = TRUE, 
                                        SRB = 1.05, 
                                        extrapLaw = "kannisto", 
                                        extrapFrom = OAnew-5, 
                                        extrapFit = seq(OAnew-30, OAnew-5, 5))   {
  
#############################
##############################
###############################
# New syntax for Abacus combin

# indicator_adult can be "35q15", "45q15"

  # first match on adult mortality indicator via lt_model_cdun_match function;
  lt_temp_adult_single <- lt_model_cdun_match_single(type = type,
                                                     Sex  = Sex,
                                                     indicator = indicator_adult,
                                                     value = value_adult,
                                                     OAnew = 130)
  
  # compute abridged lt that corresponds to single
  lt_temp_adult_abr <- lt_single2abridged(lx  = lt_temp_adult_single$lx,
                                          nLx = lt_temp_adult_single$nLx,
                                          ex  = lt_temp_adult_single$ex)
  
  
  # extract model prob of dying bw ages 20 and 25
  q5_20 <- lt_temp_adult_abr$nqx[lt_temp_adult_abr$Age==20] 
  
  # if type is CD family, then use above output 5qx pattern as the user-defined model pattern
  if (type %in% c("CD_West","CD_East","CD_North","CD_South")) {
    users_model_pattern <- lt_temp_adult_abr$nqx[1:18]
    type_combin <- "user_defined"
  } else {
    users_model_pattern <- NA
    type_combin <- type
  }
  
  # If q5 is not available, use BESTFT with 1q0 and 5q20 to get 4q1, 5q5 and 5q10. 
  if (!is.na(q1) & is.na(q5)) {
    # bestft on q1 and 5q20 to get 4q1, 5q5, and 5q10
    lt_temp_child <- lt_model_un_bestft(type = type_combin,
                                        Sex  = Sex,
                                        age_start_abridged = c(0,20),
                                        qx_abridged = c(q1,q5_20),
                                        user_pattern = users_model_pattern,
                                        lt_compute = FALSE)
    lt_temp_child <- lt_temp_child[lt_temp_child$bestft_components == 2,]
    
    # splice input q1, best_ft q for ages 1-4, 5-9, 10-14, and adult match for older ages
    qx_out <- c(q1,lt_temp_child$nqx[2:4],lt_temp_adult_abr$nqx[5:nrow(lt_temp_adult_abr)])
    
  }
  
  # If q1 is not available, use MLT and match on q5 to get q1.
  if (is.na(q1) & !is.na(q5)) {
    # match on q5 to get q1
    lt_temp_child <- lt_model_cdun_match_single(type = type,
                                                Sex  = Sex,
                                                indicator = "5q0",
                                                value = q5)
    q1 <- lt_temp_child$nqx[1]
  }
  
  # When both q1 and q5 are available, use BESTFT first with 1q0 and 5q20 then with 4q1 and 5q20. Average
  # both results to get 5q5 and 5q10.  Note that 5q20 above is from MLT.
  if (!is.na(q1) & !is.na(q5)) {
    l1 <- 100000*(1-q1)
    l5 <- 100000*(1-q5)
    q1_4 <- 1-(l5/l1)
    # bestft first with 1q0 and 5q20
    lt_out1 <- lt_model_un_bestft(type = type_combin,
                                  Sex = Sex,
                                  age_start_abridged = c(0,20),
                                  qx_abridged = c(q1,q5_20),
                                  user_pattern = users_model_pattern,
                                  lt_compute = FALSE)
    qx_out1 <- lt_out1[lt_out1$bestft_components == 2, "nqx"]
    # bestft second with 4q1 and 5q20
    lt_out2 <- lt_model_un_bestft(type = type_combin,
                                  Sex = Sex,
                                  age_start_abridged = c(1,20),
                                  qx_abridged = c(q1_4,q5_20),
                                  user_pattern = users_model_pattern,
                                  lt_compute = FALSE)
    qx_out2 <- lt_out2[lt_out2$bestft_components == 2, "nqx"]
    
    qx_out <- c(q1, q1_4, c((qx_out1+qx_out2)/2)[3:4], lt_temp_adult_abr$nqx[5:nrow(lt_temp_adult_abr)])
  }
  
lt_out_abr <- DemoTools::lt_abridged(nqx = qx_out,
                          Age = lt_temp_adult_abr$Age, 
                          Sex = Sex, 
                          axmethod = axmethod, 
                          a0rule = a0rule, 
                          OAnew = 130,
                          mod = mod,
                          extrapFrom = max(lt_temp_adult_abr$Age) - 5)

lt_out_sng <- DemoTools::lt_abridged2single(Age = lt_out_abr$Age,
                                 nMx = lt_out_abr$nMx,
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
                                 extrapFit = extrapFit,
                                 extrapFrom = extrapFrom)

 
return(lt_out_sng)

}
  
 