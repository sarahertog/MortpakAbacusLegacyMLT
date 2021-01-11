lt_model_un_bestft<-function(type, 
                             sex, 
                             age_start_abridged, 
                             qx_abridged, 
                             user_pattern=NA,
                             lt_compute = TRUE,
                             axmethod = "un",
                             a0rule = "cd",
                             mod = FALSE,
                             OAnew = 110) {
# --------------------------------------------------------------------------------------------------------------------------
#       BESTFT is derived from MORTPAK software package and customized for Abacus
#                  UNITED NATION SOFTWARE PACKAGE FOR MORTALITY MEASUREMENT    
# Sara Hertog modified to maintain method, but simplify syntax for ccmppWPP 2021 revision
# --------------------------------------------------------------------------------------------------------------------------

#' Identify best fit life table from UN model life table families fitting up to three components
#' 
#' @details xxx
#' @param sex Choose the sex of the population. 
#' #' The following options are available: \itemize{
#'   \item{\code{"f"}} -- Females;
#'   \item{\code{"m"}} -- Males.
#'   }
#' @param type Choose the family of model life tables
#'  #' The following options are available: \itemize{
#'   \item{\code{"UN_Chilean"}} -- UN Chilean;
#'   \item{\code{"UN_Far_Eastern"}} -- UN Far Eastern;
#'   \item{\code{"UN_General"}} -- UN General;
#'   \item{\code{"UN_Latin_American"}} -- UN Latin American;
#'   \item{\code{"UN_South_Asian"}} -- UN South Asian;
#'   \item{\code{"user_defined"}} -- user defined mortality pattern.
#'   }
#' @param age_start_abridged numeric. Vector of start ages of abridged age groups in input qx
#' @param qx_abridged numeric. Vector of probabilities of dying for abridged age groups to be fit
#' @param user_pattern numeric. Vector of qx to be used as model pattern when type=="user_defined"
#' 

## Do some checking/validating of input arguments

  sex <- tolower(substr(sex,1,1))
  sex_check <- sex %in% c("m","f")
  
  if (!sex_check) 
    stop("The sex argument must be either 'm' or 'f'.")
  
  type_check <- type %in% c("user_defined", "UN_Latin_American", "UN_Far_Eastern",
                            "UN_Chilean", "UN_South_Asian", "UN_General")
  
  if (!type_check) 
    stop("The type argument must be one of 'user_defined', 'UN_Latin_American', 'UN_Far_Eastern',
                            'UN_Chilean', 'UN_South_Asian', 'UN_General'.")
  
  # abridged ages
  age <- c(0,1,seq(5,80,5))
  age_check <- all(age_start_abridged %in% age)
  
  if (!age_check) 
    stop("The age argument must be starting ages for abridged age groups 0,1,5,10,15....")
  
  qx_check <- all(!is.na(qx_abridged) & qx_abridged != 0 )
  
  if (!qx_check) 
    stop("The qx_abridged argument must be all non_zero values.")
  
  if (type == "user_defined") {
    user_check <- length(user_pattern) == 18 # must be length 18 for this implementation
  }


  
  #  Empirical patterns for five UN life table families
  #  males at places 1:18; females at places 19:36
  
  if (type == "UN_Latin_American") {
    # UN-Latin American
    EMP  <- c(-1.12977,-1.49127,-2.13005,-2.40748,-2.21892,-2.01157, 
              -1.93591,-1.86961,-1.76133,-1.64220,-1.49651,-1.34160,-1.15720, 
              -.96945,-.74708,-.52259,-.29449,-.04031,-1.22452,-1.45667, 
              -2.13881,-2.46676,-2.31810,-2.14505,-2.03883,-1.93924,-1.83147, 
              -1.74288,-1.62385,-1.47924,-1.28721,-1.07443,-.83152,-.59239, 
              -.35970,-.08623)
  }
  if (type == "UN_Chilean") {
    # UN-Chilean
    EMP  <- c(-1.04722,-1.81992,-2.42430,-2.52487,-2.24491,-2.02821, 
              -1.90923,-1.78646,-1.66679,-1.52497,-1.37807,-1.21929,-1.03819, 
              -.84156,-.63201,-.42070,-.21110,+.01163,-1.12557,-1.82378, 
              -2.52319,-2.63933,-2.38847,-2.20417,-2.09701,-1.99128,-1.87930, 
              -1.75744,-1.61558,-1.45886,-1.26115,-1.05224,-.80346,-.58202, 
              -.35093,-.10587)
  }
  if (type == "UN_South_Asian") {
    # UN-South Asian
    EMP  <- c(-.97864,-1.24228,-2.01695,-2.44280,-2.35424,-2.27012, 
              -2.16833,-2.05942,-1.90053,-1.71213,-1.51120,-1.28493,-1.08192, 
              -.84671,-.62964,-.40229,-.19622,-.00129,-0.97055,-1.15424, 
              -1.93962,-2.36857,-2.19082,-2.09358,-2.04788,-1.95922,-1.87311, 
              -1.76095,-1.61425,-1.39012,-1.15515,-0.90816,-.68011,-.43231, 
              -.17489,0.05948)
  }
  if (type == "UN_Far_Eastern") {
    # UN-Far Eastern
    EMP  <- c(-1.53473,-2.15035,-2.61442,-2.66392,-2.42326,-2.23095, 
              -2.15279,-2.05765,-1.89129,-1.68244,-1.47626,-1.23020,-1.02801, 
              -.77148,-.54696,-.32996,-.11911,0.10572,-1.42596,-1.95200, 
              -2.55653,-2.68018,-2.33095,-2.15952,-2.03377,-1.94554,-1.82299, 
              -1.69084,-1.52189,-1.33505,-1.13791,-0.93765,-.72718,-.50916, 
              -.28389,-.01285)
  }
  if (type == "UN_General") {
    # UN-General
    EMP  <- c(-1.27638,-1.78957,-2.35607,-2.55527,-2.34263,-2.16193, 
              -2.09109,-2.00215,-1.86781,-1.70806,-1.52834,-1.33100,-1.12934, 
              -.91064,-.68454,-.45685,-.23002,0.00844,-1.35963,-1.77385, 
              -2.39574,-2.64549,-2.44766,-2.28991,-2.18850,-2.08535,-1.97231, 
              -1.84731,-1.69291,-1.50842,-1.30344,-1.08323,-.84402,-.59485, 
              -.34158,-.06493)
  }
  
  if (type == "user_defined") {
    # bring in user defined pattern
    EMP <- 0.50*log(user_pattern/(1.0-user_pattern))
    EMP <- rep(EMP,2)
  }
  
  if (sex == "m") {
    
    EMP <- EMP[1:18]
    
    VEC <- cbind(c(.23686,.36077,.33445,.30540,.28931,.28678,.27950,.28023,.26073,
                   .23626,.20794,.17804,.15136,.13217,.12243,.11457,.10445,.08878),
                 c(-.46007,-.68813,.06414,.12479,.24384,.10713,.06507,.03339,.02833,
                   .06473,.08705,.10620,.11305,.09467,.10809,.14738,.21037,.30918),
                 c(.09331,-.29269,-.47139,-.17403,.10715,.28842,.33620,.33692,.21354,
                   .15269,.06569,.00045,-.03731,-.10636,-.11214,-.22258,-.19631,-.38123))
  }
  if (sex == "f") {
    
    EMP <- EMP[19:36]
    
    VEC <- cbind(c(.18289,.31406,.31716,.30941,.32317,.32626,.30801,.29047,.25933,
                   .22187,.19241,.17244,.15729,.14282,.12711,.11815,.11591,.09772),
                 c(-.51009,-.52241,.08947,.03525,.03132,.07843,.06762,.00482,-.01409,
                   -.02178,.01870,.04427,.08201,.08061,.15756,.24236,.30138,.50530),
                 c(.23944,-.11117,.07566,.06268,-.26708,-.39053,-.28237,-.14277,-.05923,
                   .18909,.24773,.33679,.34121,.38290,.26731,.14442,.09697,-.13377))
  }
  


  # some parameters
  GAL1 <- VEC[,1]*VEC[,1]
  GAL2 <- VEC[,2]*VEC[,2]
  GAL3 <- VEC[,3]*VEC[,3]
  BAL1 <- VEC[,1]*VEC[,2]
  BAL2 <- VEC[,1]*VEC[,3]
  BAL3 <- VEC[,2]*VEC[,3]
  
  # initialize coefficients    
  GAMMA1 <- 0.0
  GAMMA2 <- 0.0
  GAMMA3 <- 0.0
  ALPHA1 <- 0.0
  ALPHA2 <- 0.0
  ALPHA3 <- 0.0
  BETA1 <- 0.0
  BETA2 <- 0.0
  BETA3 <- 0.0
  
  NCTR <- 0 # counter
  for (i in 1:length(age_start_abridged)){

    NCTR <- NCTR+1 # counter
    R      <- 0.5*log(qx_abridged[i]/(1.0-qx_abridged[i]))
    S      <- R-EMP[age == age_start_abridged[i]]
    SVEC1  <- S*VEC[age == age_start_abridged[i],1]
    SVEC2  <- S*VEC[age == age_start_abridged[i],2]
    SVEC3  <- S*VEC[age == age_start_abridged[i],3]
    GAMMA1 <- GAMMA1+GAL1[age == age_start_abridged[i]]
    GAMMA2 <- GAMMA2+GAL2[age == age_start_abridged[i]]
    GAMMA3 <- GAMMA3+GAL3[age == age_start_abridged[i]]
    ALPHA1 <- ALPHA1+SVEC1
    ALPHA2 <- ALPHA2+SVEC2
    ALPHA3 <- ALPHA3+SVEC3
    BETA1  <- BETA1+BAL1[age == age_start_abridged[i]]
    BETA2  <- BETA2+BAL2[age == age_start_abridged[i]]
    BETA3  <- BETA3+BAL3[age == age_start_abridged[i]]
        
  }

  lt_out <- list() # initialize component fits
  
  if(NCTR >= 3) NCTR <- 3
  NTEMP <- 4-NCTR
  for (i in NTEMP:3) {
    IC <- 4-i
    if(IC == 1 | IC == 2) GAMMA3 <- 1.0
    if(IC == 1 | IC == 2) BETA2 <- 0.0
    if(IC == 1 | IC == 2) BETA3 <- 0.0
    if(IC == 1) BETA1 <- 0.0
    if(IC == 1) GAMMA2 <- 1.0
    D  <- GAMMA1*GAMMA2*GAMMA3-GAMMA3*BETA1*BETA1-GAMMA2*BETA2*BETA2-GAMMA1*BETA3*BETA3+2.0*BETA1*BETA2*BETA3
    A1 <- ALPHA1*(GAMMA2*GAMMA3-BETA3*BETA3)+ALPHA2*(BETA2*BETA3-BETA1*GAMMA3)+ALPHA3*(BETA1*BETA3-BETA2*GAMMA2)
    A1 <- A1/D
    A2 <- ALPHA1*(BETA2*BETA3-BETA1*GAMMA3)+ALPHA2*(GAMMA1*GAMMA3-BETA2*BETA2)+ALPHA3*(BETA1*BETA2-BETA3*GAMMA1)
    A2 <- A2/D
    A3 <- ALPHA1*(BETA1*BETA3-BETA2*GAMMA2)+ALPHA2*(BETA1*BETA2-BETA3*GAMMA1)+ALPHA3*(GAMMA1*GAMMA2-BETA1*BETA1)
    A3 <- A3/D
    if (IC == 1) {
      A2 <- 0.0
      A3 <- 0.0
    }
    if (IC == 2) {
      A3 <- 0.0
    }
    CF <- EMP + A1*VEC[,1]+A2*VEC[,2]+A3*VEC[,3]
    CF <- exp(2.0*CF)/(1.0+exp(2.0*CF)) # nqx
    if (lt_compute == TRUE) {
    lt <- lt_abridged(Age = age, nqx = CF, Sex = sex, axmethod=axmethod, a0rule = a0rule, mod = mod, OAnew = OAnew)
    } else {
      lt <- data.frame(Age = age, nqx = CF)
      }
    lt$bestft_components <- IC
    lt_out[[i]] <- lt
    }

    out.data <- do.call(rbind,lt_out)

    return(out.data) 
  
}  

