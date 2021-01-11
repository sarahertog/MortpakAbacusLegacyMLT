lt_model_cdun_combin <- function(type, 
                                 sex, 
                                 q1 = NA, 
                                 q5 = NA, 
                                 indicator_adult, 
                                 value_adult,
                                 axmethod = "un",
                                 a0rule = "cd",
                                 mod = FALSE,
                                 OAnew = 110)   {
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

#' Estimate UN or Coale Demeney family model life tables, matching multiple inptu parameters
#' 
#' @details xxx
#' @param sex Choose the sex of the population. 
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


  # sex = "f"
  # type = "CD_West"
  # q1 = .04651
  # q5 = .06027
  # indicator_adult <- "45q15"
  # value_adult = .241


#############################
##############################
###############################
# New syntax for Abacus combin

# indicator_adult can be "35q15", "45q15"

  # first match on adult mortality indicator via lt_model_cdun_match function;
  lt_temp_adult <- lt_model_cdun_match(type = type,
                                       sex  = sex,
                                       indicator = indicator_adult,
                                       value = value_adult)
  q5_20 <- lt_temp_adult$nqx[lt_temp_adult$Age==20] # model prob of dying bw ages 20 and 25
  
  # if type is CD family, then use above output qx pattern as the user-defined model pattern
  if (type %in% c("CD_West","CD_East","CD_North","CD_South")) {
    users_model_pattern <- lt_temp_adult$nqx[1:18]
    type_combin <- "user_defined"
  } else {
    users_model_pattern <- NA
    type_combin <- type
  }
  
  # If q5 is not available, use BESTFT with 1q0 and 5q20 to get 4q1, 5q5 and 5q10. 
  if (!is.na(q1) & is.na(q5)) {
    # bestft on q1 and 5q20 to get 4q1, 5q5, and 5q10
    lt_temp_child <- lt_model_un_bestft(type = type_combin,
                                        sex  = sex,
                                        age_start_abridged = c(0,20),
                                        qx_abridged = c(q1,q5_20),
                                        user_pattern = users_model_pattern)
    lt_temp_child <- lt_temp_child[lt_temp_child$bestft_components == 2,]
    
    # splice input q1, best_ft q for ages 1-4, 5-9, 10-14, and adult match for older ages
    qx_out <- c(q1,lt_temp_child$nqx[2:4],lt_temp_adult$nqx[5:24])
    
  }
  
  # If q1 is not available, use MLT and match on q5 to get q1.
  if (is.na(q1) & !is.na(q5)) {
    # match on q5 to get q1
    lt_temp_child <- lt_model_cdun_match(type = type,
                                         sex  = sex,
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
                                  sex = sex,
                                  age_start_abridged = c(0,20),
                                  qx_abridged = c(q1,q5_20),
                                  user_pattern = users_model_pattern,
                                  lt_compute = FALSE)
    qx_out1 <- lt_out1[lt_out1$bestft_components == 2, "nqx"]
    # bestft second with 4q1 and 5q20
    lt_out2 <- lt_model_un_bestft(type = type_combin,
                                  sex = sex,
                                  age_start_abridged = c(1,20),
                                  qx_abridged = c(q1_4,q5_20),
                                  user_pattern = users_model_pattern,
                                  lt_compute = FALSE)
    qx_out2 <- lt_out2[lt_out2$bestft_components == 2, "nqx"]
    
    qx_out <- c(q1, q1_4, c((qx_out1+qx_out2)/2)[3:4], lt_temp_adult$nqx[5:24])
  }
  
lt_out <- lt_abridged(Age = lt_temp_adult$Age, 
                      Sex = sex, 
                      nqx = qx_out, 
                      axmethod = axmethod, 
                      a0rule = a0rule, 
                      OAnew = OAnew,
                      mod = TRUE)

  
}
  
 