## graduate UN and Coale-Demeny abridged model life tables to single year of age

# require(MortCast)
# 
# # parse MLT lookup table according to type and sex
# MLTlookup <- MortCast::MLTlookup
# MLTlookup$id <- paste(MLTlookup$type, MLTlookup$sex, MLTlookup$e0, sep = "-")
# 
# ids <- unique(MLTlookup$id)
# 
# 
# # for each abridged model life table, graduate to single year of age using ungroup package
# MLTlookup_sgl <- list()
# 
# for (i in 1:length(ids)) {
# 
#   mlt_abr <- MLTlookup[MLTlookup$id == ids[i],]
#   mlt_abr$dx <- c(-diff(mlt_abr$lx),mlt_abr$lx[nrow(mlt_abr)])
# 
#   Sex <- ifelse(mlt_abr$sex[1]==1, "m", "f")
# 
# 
#   lt_sgl <- lt_abridged2single(Age = mlt_abr$age,
#                                ndx = mlt_abr$dx,
#                                nLx = mlt_abr$Lx,
#                                Sex = Sex,
#                                a0rule = "ak",
#                                OAG = TRUE,
#                                OAnew = 130,
#                                extrapLaw = "kannisto")
# 
#   lt_sgl$id   <- ids[i]
#   lt_sgl$type <- mlt_abr$type[1]
#   lt_sgl$sex  <- Sex
#   lt_sgl$e0   <- mlt_abr$e0[1]
# 
#   MLTlookup_sgl[[i]] <- lt_sgl
#   rm(lt_sgl)
# 
# }
# 
# ## compile into a data frame
# MLTlookup_sgl <- do.call(rbind, MLTlookup_sgl)
# MLTlookup_sgl$index <- MLTlookup_sgl$e0
# 
# 
# ## check difference between e0 of original reference abridged lt and e0 of new single lt
# looke0 <- MLTlookup_sgl[MLTlookup_sgl$Age == 0,]
# 
# diff <- looke0$e0 - looke0$ex
# summary(diff)
# # > summary(diff)
# # Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# # -0.146166 -0.038588 -0.004869 -0.006832  0.024926  0.203691 
# 
# ## close enough
# 
# ##
# # save the mlt reference dataset
# 
# # index is the e0 level from original abridged MLT -- can use this later to relate single back to abridged
# MLTlookup_single <- MLTlookup_sgl[,c("type","sex","index","Age","AgeInt","nMx","nAx","nqx","lx",
#                                      "ndx","nLx","Sx","Tx","ex")]
# 
# save(MLTlookup_single, file = "data/MLTlookup_single.rda")

