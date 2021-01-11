lt_model_cdun_match <- function(type, 
                                sex, 
                                indicator, 
                                value, 
                                axmethod = "un",
                                a0rule = "cd",
                                mod = FALSE,
                                OAnew = 110)   {
# --------------------------------------------------------------------------------------------------------------------------
#      MATCH is derived from MORTPAK software package and customized for Abacus
#               UNITED NATION SOFTWARE PACKAGE FOR MORTALITY MEASUREMENT     
# Description of procedure
# This version is designed to use the same method as used in SLON and uses a look up table provided by the calling program
# Source of lookup table is MortCast MLTlookup
# Sara Hertog modified to maintain method, but simplify syntax for ccmppWPP 2021 revision
# --------------------------------------------------------------------------------------------------------------------------

#' Estimate UN or Coale Demeney family model life tables, matching on one input parameter
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
#' @param indicator character. Life table indicator to match on
#'  #' The following options are available: \itemize{
#'   \item{\code{"1q0"}} -- Probability of dying between birth and age 1; 
#'   \item{\code{"5q0"}} -- Probability of dying between birth and age 5; 
#'   \item{\code{"35q15"}} -- Probability of dying between age 15 and age 50;
#'   \item{\code{"45q15"}} -- Probability of dying between age 15 and age 60;
#'   \item{\code{"e0"}} -- Life expectancy at  birth;
#' @param values numeric. Values of the indicators to match on;
#' @inheritParams lt_abridged
#' @return data.frame. with two columns: age, giving abridged age groups, and mx
#' with model age specific mortality rates for abridged age groups
#' @importFrom stats uniroot MortCast DemoTools
#' @examples 
#' 


  # sex = "f"
  # type = "CD_West"
  # indicator <- "5q0"
  # value = .150

  require(MortCast)
  require(DemoTools)
  
  # parse MLT lookup table according to type and sex
  sexcode <- ifelse(substr(sex,1,1) == "m", 1, ifelse(substr(sex,1,1) == "f", 2, NA))
  MLTlookup <- MortCast::MLTlookup
  MLTlookup <- MLTlookup[MLTlookup$type == type & MLTlookup$sex == sexcode,]

    if (indicator == "e0") {
      mlts       <- MLTlookup
      mlts$level <- mlts$e0
    }
    
    if (indicator == "1q0") {
      # compute q1 levels for model life tables
      mlts       <- MLTlookup[MLTlookup$age == 1, c("e0", "lx")]
      mlts$level <- 1-(mlts$lx / 100000)
      mlts       <- merge(MLTlookup, mlts[, c("e0", "level")], by = "e0")
      
    }
  
    if (indicator == "5q0") {
      # compute q5 levels for model life tables
      mlts       <- MLTlookup[MLTlookup$age == 5, c("e0", "lx")]
      mlts$level <- 1-(mlts$lx / 100000)
      mlts       <- merge(MLTlookup, mlts[, c("e0", "level")], by = "e0")
      
    }
    
    if (indicator == "35q15") {
      # compute 35q15 levels for model life tables
      mlts       <- MLTlookup[MLTlookup$age %in% c(15,50), c("e0", "age", "lx")]
      mlts       <- reshape(mlts, direction = "wide", timevar = "age", 
                            idvar = "e0", sep = "_")
      mlts$level <- 1 - (mlts$lx_50 / mlts$lx_15)
      mlts       <- merge(MLTlookup, mlts[, c("e0", "level")], by = "e0")
    }
    
    if (indicator == "45q15") {
      # compute 45q15 levels for model life tables
      mlts       <- MLTlookup[MLTlookup$age %in% c(15,60), c("e0", "age", "lx")]
      mlts       <- reshape(mlts, direction = "wide", timevar = "age", 
                            idvar = "e0", sep = "_")
      mlts$level <- 1 - (mlts$lx_60 / mlts$lx_15)
      mlts       <- merge(MLTlookup, mlts[, c("e0", "level")], by = "e0")
    }
    
    # sort by level and age
    mlts         <- mlts[order(mlts$level, mlts$age),]
    
    # identify the model life tables with levels just below and above the value to match
    lvls   <- unique(mlts$level)
    iord   <- which.min(abs(value - lvls)) # identify closest level to value (could be higher or lower)
    lower  <- ifelse(lvls[iord] <= value, iord, iord-1)
    higher <- lower + 1
    
    # parse the two matched mlts
    mlt_match <- mlts[mlts$level %in% c(lvls[c(lower,higher)]), c("level","age","mx")]
    mlt_match <- reshape(mlt_match, direction = "wide", timevar = "level", 
                            idvar = "age", sep = "_")
    
    # interpolate log(mx) between the two matched mlts according to the position of the input value
    # relative to the two matched levels (replicates Abacus approach)
    pct_val <- (value - lvls[lower]) / (lvls[higher]-lvls[lower])
    mx_hat  <- exp((1.0-pct_val)*log(mlt_match[,2])+pct_val*log(mlt_match[,3]))  # m(x,n)
    
    # compute the life table
    lt_out <- DemoTools::lt_abridged(Age = mlt_match$age[1:27], # need to truncate series to OA<130
                                     nMx = mx_hat[1:27], 
                                     Sex = sex, 
                                     axmethod = axmethod, 
                                     a0rule = a0rule,
                                     OAG = FALSE,
                                     OAnew = OAnew,
                                     mod = mod)
  

  return(lt_out)
  
}


