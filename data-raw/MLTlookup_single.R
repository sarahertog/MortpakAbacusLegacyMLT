## graduate UN and Coale-Demeny abridged model life tables to single year of age

library(MortCast)
library(DemoTools)
library(tidyverse)
library(tictoc)

# parse MLT lookup table according to type and sex
MLTlookup <- MortCast::MLTlookup
MLTlookup$id <- paste(MLTlookup$type, MLTlookup$sex, MLTlookup$e0, sep = "-")

ids <- unique(MLTlookup$id)


# for each abridged model life table, graduate to single year of age using ungroup package
MLTlookup_sgl <- list()

tic()
for (i in 1:length(ids)) {

  mlt <- MLTlookup[MLTlookup$id == ids[i],]
  Sex <- ifelse(mlt$sex[1]==1, "m", "f")
  
  lt_sgl <- lt_abridged2single(Age = mlt$age[1:27],
                               nMx = mlt$mx[1:27],
                               Sex = Sex,
                               axmethod = "un",
                               a0rule = "ak",
                               OAG = TRUE,
                               OAnew = 130,
                               extrapLaw = "kannisto",
                               extrapFit = seq(110,125,5))
  
  lt_sgl$id   <- ids[i]
  lt_sgl$type <- mlt$type[1]
  lt_sgl$sex  <- Sex
  lt_sgl$e0_mlt   <- mlt$e0[1]

  MLTlookup_sgl[[i]] <- lt_sgl
  rm(lt_sgl)

}
toc()

## compile into a data frame
MLTlookup_sgl <- do.call(rbind, MLTlookup_sgl)
MLTlookup_sgl$index <- MLTlookup_sgl$e0_mlt


## check difference between e0 of original reference abridged lt and e0 of new single lt
looke0 <- MLTlookup_sgl[MLTlookup_sgl$Age == 0,]

diff1 <- looke0$e0_mlt - looke0$ex
summary(diff1)

# > summary(diff1)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# -0.289592 -0.083168 -0.038359 -0.049091 -0.005089  0.218819 
## not exact replication of original e0 value but close enough

MLTlookup_single <- MLTlookup_sgl[,c("type","sex","index","Age","AgeInt","nMx","nAx","nqx","lx",
                                     "ndx","nLx","Sx","Tx","ex")]

save(MLTlookup_single, file = "data/MLTlookup_single.rda")


load("data/MLTlookup_single.rda")

# Graduated e0 is slightly different from original abridged e0 that was the index level
# Check whether the mlt ranked by original abridged e0 index is the same as the rank by the new e0 of the graduated lt
mlt0 <- MLTlookup_single %>% 
  filter(Age == 0) %>% 
  arrange(type, sex, index) %>% 
  group_by(type, sex) %>% 
  mutate(rank_index = 1:length(index)) %>% 
  arrange(type, sex, ex) %>% 
  mutate(rank_e0 = 1:length(ex)) %>% 
  ungroup() %>% 
  mutate(rank_agree = rank_index==rank_e0)

table(mlt0$rank_agree)
# > table(mlt0$rank_agree)
# 
# TRUE 
# 630 

# rank is the same so we can use the new e0 as the index


MLTlookup$id <- paste(MLTlookup$type, MLTlookup$e0, sep = " - ")
MLTlookup$AgeInt <- age2int(MLTlookup$age)
MLTlookup_single$id <- paste(MLTlookup_single$type, MLTlookup_single$index, sep = " - ")

ids <- unique(MLTlookup_single$id)

pdf(file = "charts/Plot_single.pdf",
    width=8,height=10)

par(mfrow = c(3,2))

for (i in 1:length(ids)) {
  
  lt_abr <- MLTlookup %>% 
    filter(id == ids[i]) %>% 
    group_by(sex) %>% 
    mutate(AgeInt = age2int(age)) %>% 
    ungroup() 
  lt_new <- MLTlookup_single %>% 
    filter(id == ids[i])
    
  
  plot(c(0,130), c(min(log(lt_new$nMx)),0), type="n", xlab = "Age", ylab = "log(nMx)", 
       main = lt_new$id[1])
  lines(lt_new$Age[lt_new$sex == "m"] + (lt_new$AgeInt[lt_new$sex == "m"]/2),
        log(lt_new$nMx[lt_new$sex == "m"]),
        col = "blue",
        lty = 1)
  points(lt_abr$age[lt_abr$sex == 1] + (lt_abr$AgeInt[lt_abr$sex == 1]/2),
         log(lt_abr$mx[lt_abr$sex == 1]),
         col = "blue",
         pch = 19)
  lines(lt_new$Age[lt_new$sex == "f"] + (lt_new$AgeInt[lt_new$sex == "f"]/2),
        log(lt_new$nMx[lt_new$sex == "f"]),
        col = "red",
        lty = 1)
  points(lt_abr$age[lt_abr$sex == 2] + (lt_abr$AgeInt[lt_abr$sex == 2]/2),
         log(lt_abr$mx[lt_abr$sex == 2]),
         col = "red",
         pch = 19)
  legend("bottomright", c("males - abr mlt",
                          "males - new mlt",
                          "females - abr mlt",
                          "females - new mlt"),
         col = c("blue", "blue","red", "red"),
         lty = c(NA,1,NA,1),
         pch = c(19,NA,19,NA))
  
}
dev.off()



