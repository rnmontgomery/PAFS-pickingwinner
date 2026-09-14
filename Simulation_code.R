# Packages ----------------------------

library(MASS)
library(LaplacesDemon)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(see)
library(rsimsum)
library(ggpubr)


# Functions -----------------------------

# Function to draw from Normal Inverse-Wishart conjugate

draw_NIW <- function(data,  mu0, lambda0, v0, k0, ndraws){
  
  ybar <- colMeans(data)
  n_dat <- dim(data)[1]
  S <- cov(data)*(n_dat-1)
  
  # Posterior parameters
  mun <- (k0/(k0+n_dat))*mu0 + (n_dat/(k0+n_dat))*ybar
  kn <- k0 + n_dat
  vn <- v0 + n_dat
  lambda_n <- lambda0 + S + ((k0*n_dat)/(k0+n_dat))*(matrix(ybar-mu0)%*%t(matrix(ybar-mu0)))
  
  
  draws <- replicate(ndraws, {
    epsilon_draw <- rinvwishart(nu = vn, S = (lambda_n))
    mu_draw <- mvrnorm(n = 1, mu = mun, Sigma = epsilon_draw/kn)
    delta <- mu_draw/sqrt(diag(epsilon_draw))
    delta
  }, simplify = TRUE)
  
  return(t(draws))
  
}



# Run the following code to replicate simulation studies.
# The code takes several days to run, therefore
# The results was saved, and are available on github as combined_simresults.csv

set.seed(4387)


nvector = c(30,50,70,100)
endpoints = c(2,4,8)
correlation = c(0, 0.25)
meanvect = c("winner", "stagger")
effectlow = 0.25
effecthigh = 0.35
mu0v = c("pessimistic", "optimistic")
nsims = 5000
ndraws = 1000


  rows <- length(nvector)*length(endpoints)*length(correlation)*length(meanvect)*length(mu0v)*nsims
  emp <- 1
  
  outmat <- matrix(NA,rows,25)
  
  for (i in 1:length(nvector)){
    
    N <- nvector[i]
    
   for (j in 1:length(endpoints)){
    
    Nout <- endpoints[j]
    
   for (k in 1:length(correlation)){
    
    offdiag <- correlation[k]
    
   for (l in 1:length(meanvect)){
     
     if (meanvect[l] == "winner"){
       effect <- c(rep(effectlow, Nout-1), effecthigh)
     } else if (meanvect[l] == "stagger"){
       effect <- c( seq(effectlow, effecthigh, length.out = Nout))
     }
     
     for(m in 1:length(mu0v)){
       if (mu0v[m] == "pessimistic"){
         mu0 <- c(rep(0, Nout))
       } else if (mu0v[m] == "optimistic"){
         mu0 <- effect
       }

     sig <-  diag(Nout)
     sig[row(sig) != col(sig)] <- offdiag
     
     for (n in 1:nsims){
       mvndata <- mvrnorm(n = N, mu = c(effect), Sigma = sig)
       
       # Sample means
       ybars <- colMeans(mvndata)
       ybars
       sds <- apply(mvndata,2,sd)
       
       # Mean choice
       deltahat <- ybars/sds
       meanc <- which.max(deltahat)
       
       means_sampN <- power.t.test(n= NULL, delta = deltahat[meanc],sd = 1, power = 0.80, type = "two.sample" )$n # powering based on this
       pow <- power.t.test(n= means_sampN, delta = effect[meanc],sd = 1, power = NULL, type = "two.sample" )$power
       
       
       # Bayesian congjuate MvN
       #k0 <- N*k0p
       lambda0 = diag(Nout)
       v0 = Nout + 1
       ybar <- colMeans(mvndata)
       n_dat <- dim(mvndata)[1]
       S <- cov(mvndata)*(n_dat-1)
       
       
       # Draw from posterior .1
       draws1 <- draw_NIW(data = mvndata,  mu0, lambda0, v0, 0.1*N, ndraws = ndraws)
       
       max1 <- apply(draws1,1,which.max)
       bprob1 <- table(factor(max1, levels = 1:Nout))/ndraws
       deltaB1 <- colMeans(draws1)
       
       
       # Bayes choice .1
       bayesc1 <- as.numeric(names(which(bprob1 == max(bprob1))))
       if (length(bayesc1) >1){
         bayesc1 <- bayesc1[1]
       }
       deltahatB1 <- deltaB1[bayesc1]
       
       means_sampNB1 <- power.t.test(n= NULL, delta = deltahatB1,sd = 1, power = 0.80, type = "two.sample" )$n
       powB1 <- power.t.test(n= means_sampNB1, delta = effect[bayesc1],sd = 1, power = NULL, type = "two.sample" )$power
       
       
       # Draw from posterior .5
       draws5 <- draw_NIW(data = mvndata,  mu0, lambda0, v0, 0.5*N, ndraws = ndraws)
       
       max5 <- apply(draws5,1,which.max)
       bprob5 <- table(factor(max5, levels = 1:Nout))/ndraws
       deltaB5 <- colMeans(draws5)
       
       
       # Bayes choice .5
       bayesc5 <- as.numeric(names(which(bprob5 == max(bprob5))))
       if (length(bayesc5) >1){
         bayesc5 <- bayesc5[1]
       }
       deltahatB5 <- deltaB5[bayesc5]
       
       means_sampNB5 <- power.t.test(n= NULL, delta = deltahatB5,sd = 1, power = 0.80, type = "two.sample" )$n
       powB5 <- power.t.test(n= means_sampNB5, delta = effect[bayesc5],sd = 1, power = NULL, type = "two.sample" )$power
       
       
       # Draw from posterior 100
       draws100 <- draw_NIW(data = mvndata,  mu0, lambda0, v0, 1*N, ndraws = ndraws)
       
       max100 <- apply(draws100,1,which.max)
       bprob100 <- table(factor(max100, levels = 1:Nout))/ndraws
       deltaB100 <- colMeans(draws100)
       
       
       # Bayes choice 100
       bayesc100 <- as.numeric(names(which(bprob100 == max(bprob100))))
       if (length(bayesc100) >1){
         bayesc100 <- bayesc100[1]
       }
       deltahatB100 <- deltaB100[bayesc100]
       
       means_sampNB100 <- power.t.test(n= NULL, delta = deltahatB100,sd = 1, power = 0.80, type = "two.sample" )$n
       powB100 <- power.t.test(n= means_sampNB100, delta = effect[bayesc100],sd = 1, power = NULL, type = "two.sample" )$power
       
       
       # empty_rows <- which(rowSums(is.na(outmat)) == ncol(outmat))
       # emp <- empty_rows[1]
       if (emp %% 100 == 0){
         print(emp)
       }
       
       outmat[emp,1] <- N
       outmat[emp,2] <- Nout
       outmat[emp,3] <- offdiag
       outmat[emp,4] <- l
       outmat[emp,5] <- m
       
       # Cohen's D outcomes
       outmat[emp,6] <- meanc # mean choice
       outmat[emp,7] <- deltahat[meanc] # Estimated effect for choice
       outmat[emp,8] <- deltahat[meanc]/effect[meanc] # Bias
       outmat[emp,9] <- 0.8 - pow # power loss assuming powered for 80% trial given correct
       outmat[emp,10] <- ifelse(meanc == Nout,1,0) # Correct choice, yes/no
       
       
       # Bayesian 0.1
       outmat[emp,11] <- bayesc1 # bayes choice
       outmat[emp,12] <- deltahatB1# Estimated effect for choice
       outmat[emp,13] <- deltahatB1/effect[bayesc1] # Bias
       outmat[emp,14] <- 0.8 - powB1 # power loss assuming powered for 80% trial given correct
       outmat[emp,15] <- ifelse(bayesc1 == Nout,1,0) # Correct choice, yes/no
       
       
       
       # Bayesian 0.5
       outmat[emp,16] <- bayesc5 # bayes choice
       outmat[emp,17] <- deltahatB5# Estimated effect for choice
       outmat[emp,18] <- deltahatB5/effect[bayesc5] # Bias
       outmat[emp,19] <- 0.8 - powB5 # power loss assuming powered for 80% trial given correct
       outmat[emp,20] <- ifelse(bayesc5 == Nout,1,0) # Correct choice, yes/no
       
       
       # Bayesian 1.0
       outmat[emp,21] <- bayesc100 # bayes choice
       outmat[emp,22] <- deltahatB100# Estimated effect for choice
       outmat[emp,23] <- deltahatB100/effect[bayesc100] # Bias
       outmat[emp,24] <- 0.8 - powB100 # power loss assuming powered for 80% trial given correct
       outmat[emp,25] <- ifelse(bayesc100 == Nout,1,0) # Correct choice, yes/no
       
       emp <- emp+1
         
        }
       }
     }
   }
  }
}


  
# return(outmat)
# }
  


# Read the simulated results in and calculate estimates ------------------
  
  combined <- read.csv("--------YOURPATH-------------\\combined_simresults.csv", header = TRUE)
  
  
#combined <- read.csv("C:\\Users\\rmontgomery\\OneDrive - University of Kansas Medical Center\\Papers\\Pilot pick a winner\\Resubmission\\combined_simresults.csv", header = TRUE)

colnames(combined) <- c("N","Endpoints", "Correlation", "Meanvect", "mu0",
                      "meanchoice", "deltahat", "bias", "powerloss", "meancorrect",
                      "B01choice", "B01deltahat", "B01bias", "B01powerloss", "B01meancorrect",
                      "B50choice", "B50deltahat", "B50bias", "B50powerloss", "B50meancorrect",
                      "B100choice", "B100deltahat", "B100bias", "B100powerloss", "B100meancorrect")

# Meanvect: Winner, Stagger
# Mu0: Pessimistic, Correct

# Probabilities of correct choice ------

# Sample mean
combined %>% dplyr::group_by(N, Endpoints, Correlation, Meanvect, mu0) %>% dplyr::summarize(sample_correct= mean(meancorrect)) -> sample_correct
sample_correct$Meanvect <- factor(sample_correct$Meanvect ,levels = c(1,2), labels = c("winner","stagger"))
sample_correct$mu0 <- factor(sample_correct$mu0 ,levels = c(1,2), labels = c("pessimistic","correct"))

# Bayes 10
combined %>% dplyr::group_by(N, Endpoints, Correlation, Meanvect, mu0) %>% dplyr::summarize(Bayes_10correct= mean(B01meancorrect)) -> B01mean_correct
B01mean_correct$Meanvect <- factor(B01mean_correct$Meanvect ,levels = c(1,2), labels = c("winner","stagger"))
B01mean_correct$mu0 <- factor(B01mean_correct$mu0 ,levels = c(1,2), labels = c("pessimistic","correct"))

# bayes 50
combined %>% dplyr::group_by(N, Endpoints, Correlation, Meanvect, mu0) %>% dplyr::summarize(B50_meancorrect = mean(B50meancorrect )) -> B50mean_correct
B50mean_correct$Meanvect <- factor(B50mean_correct$Meanvect ,levels = c(1,2), labels = c("winner","stagger"))
B50mean_correct$mu0 <- factor(B50mean_correct$mu0 ,levels = c(1,2), labels = c("pessimistic","correct"))

# bayes 100
combined %>% dplyr::group_by(N, Endpoints, Correlation, Meanvect, mu0) %>% dplyr::summarize(B100_meancorrect = mean(B100meancorrect )) -> B100mean_correct
B100mean_correct$Meanvect <- factor(B100mean_correct$Meanvect ,levels = c(1,2), labels = c("winner","stagger"))
B100mean_correct$mu0 <- factor(B100mean_correct$mu0 ,levels = c(1,2), labels = c("pessimistic","correct"))




# Bias of choices -------------------

# Sample mean
combined %>% dplyr::group_by(N, Endpoints, Correlation, Meanvect, mu0) %>% dplyr::summarize(meanbias= mean(bias)) -> bias_sample
bias_sample$Meanvect <- factor(bias_sample$Meanvect ,levels = c(1,2), labels = c("winner","stagger"))
bias_sample$mu0 <- factor(bias_sample$mu0 ,levels = c(1,2), labels = c("pessimistic","correct"))

# Monte Carlo Standard error
combined %>% dplyr::group_by(N, Endpoints, Correlation, Meanvect, mu0) %>% dplyr::summarize(mmcse= sd(bias)/sqrt(n())) -> bias_samplese
bias_samplese$Meanvect <- factor(bias_samplese$Meanvect ,levels = c(1,2), labels = c("winner","stagger"))
bias_samplese$mu0 <- factor(bias_samplese$mu0 ,levels = c(1,2), labels = c("pessimistic","correct"))
max(bias_samplese$mmcse)

# Bias 10
combined %>% dplyr::group_by(N, Endpoints, Correlation, Meanvect, mu0) %>% dplyr::summarize(meanbias= mean(B01bias)) -> bias_delta_10
bias_delta_10$Meanvect <- factor(bias_delta_10$Meanvect ,levels = c(1,2), labels = c("winner","stagger"))
bias_delta_10$mu0 <- factor(bias_delta_10$mu0 ,levels = c(1,2), labels = c("pessimistic","correct"))

# Monte Carlo Standard error
combined %>% dplyr::group_by(N, Endpoints, Correlation, Meanvect, mu0) %>% dplyr::summarize(mmcse= sd(B01bias)/sqrt(n())) -> bias_delta_10se
bias_delta_10se$Meanvect <- factor(bias_delta_10se$Meanvect ,levels = c(1,2), labels = c("winner","stagger"))
bias_delta_10se$mu0 <- factor(bias_delta_10se$mu0 ,levels = c(1,2), labels = c("pessimistic","correct"))
max(bias_delta_10se$mmcse)

# Bias Bayes 50
combined %>% dplyr::group_by(N, Endpoints, Correlation, Meanvect, mu0) %>% dplyr::summarize(meanbias= mean(B50bias)) -> bias_delta_50
bias_delta_50$Meanvect <- factor(bias_delta_50$Meanvect ,levels = c(1,2), labels = c("winner","stagger"))
bias_delta_50$mu0 <- factor(bias_delta_50$mu0 ,levels = c(1,2), labels = c("pessimistic","correct"))

# Monte Carlo Standard error
combined %>% dplyr::group_by(N, Endpoints, Correlation, Meanvect, mu0) %>% dplyr::summarize(mmcse= sd(B50bias)/sqrt(n())) -> bias_delta_50se
bias_delta_50se$Meanvect <- factor(bias_delta_50se$Meanvect ,levels = c(1,2), labels = c("winner","stagger"))
bias_delta_50se$mu0 <- factor(bias_delta_50se$mu0 ,levels = c(1,2), labels = c("pessimistic","correct"))
max(bias_delta_50se$mmcse)

# Bias Bayes 100
combined %>% dplyr::group_by(N, Endpoints, Correlation, Meanvect, mu0) %>% dplyr::summarize(meanbias= mean(B100bias)) -> bias_delta_100
bias_delta_100$Meanvect <- factor(bias_delta_100$Meanvect ,levels = c(1,2), labels = c("winner","stagger"))
bias_delta_100$mu0 <- factor(bias_delta_100$mu0 ,levels = c(1,2), labels = c("pessimistic","correct"))

combined %>% dplyr::group_by(N, Endpoints, Correlation, Meanvect, mu0) %>% dplyr::summarize(mmcse= sd(B100bias)/sqrt(n())) -> bias_delta_100se
bias_delta_100se$Meanvect <- factor(bias_delta_100se$Meanvect ,levels = c(1,2), labels = c("winner","stagger"))
bias_delta_100se$mu0 <- factor(bias_delta_100se$mu0 ,levels = c(1,2), labels = c("pessimistic","correct"))
max(bias_delta_100se$mmcse)



# Power loss ------------------
# Sample Cohen's D
combined %>% dplyr::group_by(N, Endpoints, Correlation, Meanvect, mu0) %>% dplyr::summarize(powerloss= mean(powerloss)) -> pl_sample
pl_sample$Meanvect <- factor(pl_sample$Meanvect ,levels = c(1,2), labels = c("winner","stagger"))
pl_sample$mu0 <- factor(pl_sample$mu0 ,levels = c(1,2), labels = c("pessimistic","correct"))

#MCSE
combined %>% dplyr::group_by(N, Endpoints, Correlation, Meanvect, mu0) %>% dplyr::summarize(mcse= sd(powerloss)/sqrt(n()))  -> pl_samplese
pl_samplese$Meanvect <- factor(pl_samplese$Meanvect ,levels = c(1,2), labels = c("winner","stagger"))
pl_samplese$mu0 <- factor(pl_samplese$mu0 ,levels = c(1,2), labels = c("pessimistic","correct"))
max(pl_samplese$mcse)

# Bayes k = 0.1*N
combined %>% dplyr::group_by(N, Endpoints, Correlation, Meanvect, mu0) %>% dplyr::summarize(powerloss= mean(B01powerloss)) -> pl_B10
pl_B10$Meanvect <- factor(pl_B10$Meanvect ,levels = c(1,2), labels = c("winner","stagger"))
pl_B10$mu0 <- factor(pl_B10$mu0 ,levels = c(1,2), labels = c("pessimistic","correct"))

#MCSE
combined %>% dplyr::group_by(N, Endpoints, Correlation, Meanvect, mu0) %>% dplyr::summarize(mcse= sd(B01powerloss)/sqrt(n()))  -> pl_B10se
pl_B10se$Meanvect <- factor(pl_B10se$Meanvect ,levels = c(1,2), labels = c("winner","stagger"))
pl_B10se$mu0 <- factor(pl_B10se$mu0 ,levels = c(1,2), labels = c("pessimistic","correct"))
max(pl_B10se$mcse)

# Bayes k = 0.5*N
combined %>% dplyr::group_by(N, Endpoints, Correlation, Meanvect, mu0) %>% dplyr::summarize(powerloss= mean(B50powerloss)) -> pl_B50
pl_B50$Meanvect <- factor(pl_B50$Meanvect ,levels = c(1,2), labels = c("winner","stagger"))
pl_B50$mu0 <- factor(pl_B50$mu0 ,levels = c(1,2), labels = c("pessimistic","correct"))

#MCSE
combined %>% dplyr::group_by(N, Endpoints, Correlation, Meanvect, mu0) %>% dplyr::summarize(mcse= sd(B50powerloss)/sqrt(n()))  -> pl_B50se
pl_B50se$Meanvect <- factor(pl_B50se$Meanvect ,levels = c(1,2), labels = c("winner","stagger"))
pl_B50se$mu0 <- factor(pl_B50se$mu0 ,levels = c(1,2), labels = c("pessimistic","correct"))
max(pl_B50se$mcse)

# Bayes k = 1.0*N
combined %>% dplyr::group_by(N, Endpoints, Correlation, Meanvect, mu0) %>% dplyr::summarize(powerloss= mean(B100powerloss)) -> pl_B100
pl_B100$Meanvect <- factor(pl_B100$Meanvect ,levels = c(1,2), labels = c("winner","stagger"))
pl_B100$mu0 <- factor(pl_B100$mu0 ,levels = c(1,2), labels = c("pessimistic","correct"))

#MCSE
combined %>% dplyr::group_by(N, Endpoints, Correlation, Meanvect, mu0) %>% dplyr::summarize(mcse= sd(B100powerloss)/sqrt(n()))  -> pl_B100se
pl_B100se$Meanvect <- factor(pl_B100se$Meanvect ,levels = c(1,2), labels = c("winner","stagger"))
pl_B100se$mu0 <- factor(pl_B100se$mu0 ,levels = c(1,2), labels = c("pessimistic","correct"))
max(pl_B100se$mcse)

# Figure 1 for paper ----------------


# Two panel plot bias with N = 50
# Correlation = 0, Stagger effect

# Figure 1 a ----------------

combined50_2 <- combined[combined$N == 50 &combined$Endpoints == 2 & combined$Correlation == 0 & combined$Meanvect == 2 & combined$mu0 == 1,]

for (k in 1:dim(combined50_2)[1]){
  
  if(combined50_2[k,"meanchoice"] == 1){
    combined50_2[k,"meanchoice"] <- 2
  }else if(combined50_2[k,"meanchoice"] == 2){
    combined50_2[k,"meanchoice"] <- 1
  }
  }

combined50_2 |> group_by(meanchoice) |> summarize(mean = median(deltahat))

cutoffs <- data.frame(
  meanchoice = c(1, 2),
  cutoff = c(0.35, 0.25)
)

combined50_2 <- combined50_2 %>%
  group_by(meanchoice) %>%
  group_modify(~{
    d <- density(.x$deltahat)
    data.frame(x = d$x, y = d$y)
  }) %>%
  left_join(cutoffs, by = "meanchoice")


combined50_2$val <- ifelse(combined50_2$meanchoice == 1, 0.35, 0.25)

labs <- combined50_2 |>
  distinct(meanchoice, val) |>
  with(setNames(paste0("\u03b8 = ", val), meanchoice))

Fig1a <- ggplot(combined50_2, aes(x, y)) +
  geom_area(fill = "darkgrey") +

  geom_area(
    data = subset(combined50_2, x >= cutoff),
    fill = "orangered4"
  ) +
  facet_wrap(~meanchoice, nrow = 2, labeller = labeller(meanchoice = labs),
             scales = "free_x")+
  scale_x_continuous(
    limits = c(0, 1),      # shared range
    breaks = c(0, 0.25, 0.35,0.5, 0.70, 1.0)  # tick marks you want on each facet
  )+
  theme_bw()+
  ylab("Density estimates for winning endpoint") + xlab("Effect size") 



# Figure 1 b -----------------------
combined50_4 <- combined[combined$N == 50 &combined$Endpoints == 4 & combined$Correlation == 0 & combined$Meanvect == 2 & combined$mu0 == 1,]


for (k in 1:dim(combined50_4)[1]){
  
  if(combined50_4[k,"meanchoice"] == 4){
    combined50_4[k,"meanchoice"] <- 1
  }else if(combined50_4[k,"meanchoice"] == 3){
    combined50_4[k,"meanchoice"] <- 2
  }else if(combined50_4[k,"meanchoice"] == 2){
    combined50_4[k,"meanchoice"] <- 3
  }else if(combined50_4[k,"meanchoice"] == 1){
    combined50_4[k,"meanchoice"] <- 4
  }
}

combined50_4 |> group_by(meanchoice) |> summarize(mean = median(deltahat))


cutoffs <- data.frame(
  meanchoice = c(1, 2, 3, 4),
  cutoff = c(0.35,0.3167, 0.2834, 0.25)
)

combined50_4 <- combined50_4 %>%
  group_by(meanchoice) %>%
  group_modify(~{
    d <- density(.x$deltahat)
    data.frame(x = d$x, y = d$y)
  }) %>%
  left_join(cutoffs, by = "meanchoice")


combined50_4$val <- ifelse(combined50_4$meanchoice == 1, 0.35,
                           ifelse(combined50_4$meanchoice == 2, 0.3167,
                                  ifelse(combined50_4$meanchoice == 3,0.2834, 0.25 )))
                                  
                                  
labs <- combined50_4 |>
  distinct(meanchoice, val) |>
  with(setNames(paste0("\u03b8 = ", val), meanchoice))


facetspecific_breaks <- list(
  "1" = c(0, 0.25, 0.35, 0.50, 0.70, 1),
  "2" = c(0, 0.25, 0.3167, 0.50, 0.70, 1),
  "3" = c(0, 0.2834, 0.35, 0.50, 0.70, 1),
  "4" = c(0, 0.25, 0.35, 0.50, 0.70, 1)
)


Fig1b <-ggplot(combined50_4, aes(x, y)) +
  geom_area(fill = "darkgrey") +
  geom_area(
    data = subset(combined50_4, x >= cutoff),
    fill = "orangered4"
  ) +
  facet_wrap(~meanchoice, nrow = 4, labeller = labeller(meanchoice = labs),
             scales = "free_x") +
  scale_x_continuous(
    limits = c(0, 1),      
    breaks = c(0, 0.25, 0.35,0.5, 0.70, 1.0) )+
  theme_bw()+
  ylab("") + xlab("Effect size")





# Nested loop plot -------------------

# Correlation has little impact on the results, restricting to only 0

combined_corr0 <- combined[combined$Correlation == 0,]

stack_cd <- combined_corr0[ , c("N","Endpoints", "Correlation", "Meanvect", "mu0",
                                "meanchoice", "deltahat", "bias", "powerloss", "meancorrect")]
colnames(stack_cd) <- c("N","Endpoints", "Correlation", "Meanvect", "\u03BC",
                        "Chosen", "Estimate", "Bias ratio", "powerloss", "Correct")
stack_cd$Method <- "Cohen's D"

stack_B10 <- combined_corr0[ , c("N","Endpoints", "Correlation", "Meanvect", "mu0",
                                 "B01choice", "B01deltahat", "B01bias", "B01powerloss", "B01meancorrect")]
colnames(stack_B10) <- c("N","Endpoints", "Correlation", "Meanvect", "\u03BC",
                         "Chosen", "Estimate", "Bias ratio", "powerloss", "Correct")
stack_B10$Method <- "Bayes, k = 0.1"

stack_B50 <- combined_corr0[ , c("N","Endpoints", "Correlation", "Meanvect", "mu0",
                                 "B50choice", "B50deltahat", "B50bias", "B50powerloss", "B50meancorrect")]
colnames(stack_B50) <- c("N","Endpoints", "Correlation", "Meanvect", "\u03BC",
                         "Chosen", "Estimate", "Bias ratio", "powerloss", "Correct")
stack_B50$Method <- "Bayes, k = 0.5"

stack_B100 <- combined_corr0[ , c("N","Endpoints", "Correlation", "Meanvect", "mu0",
                                  "B100choice", "B100deltahat", "B100bias", "B100powerloss", "B100meancorrect")]
colnames(stack_B100) <- c("N","Endpoints", "Correlation", "Meanvect", "\u03BC",
                          "Chosen", "Estimate", "Bias ratio", "powerloss", "Correct")
stack_B100$Method <- "Bayes, k = 1"


# Stacking the data to use rsimsum functions
stacked <- rbind(stack_cd, stack_B10, stack_B50, stack_B100)


# Nested loop plot for proportion of correct choices----------

# Add labels for plot
stacked$μ <- factor(stacked$μ, levels = c(1,2), labels = c("Pessimistic Prior", "Correct Prior"))
stacked$Meanvect <- factor(stacked$Meanvect, levels = c(1,2), labels = c("Winner", "Stagger"))


nested_correct <- rsimsum::simsum(
  data = stacked, estvarname = "Correct", true = 1, se = NULL,
  methodvar = "Method", by = c("N","Endpoints",  "Meanvect", "\u03BC"),
  ref = "Cohen's D"
)
Nestedplot_correct <- autoplot(nested_correct, type = "nlp", stats = "thetamean", top = TRUE,
              target = 0, zoom = 1) 

Nestedplot_correct+ 
  ylab("Proportion of correct choice") + 
  theme_bw() 



Nestedplot_bias <- rsimsum::simsum(
  data = stacked, estvarname = "Bias ratio", true = 1, se = NULL,
  methodvar = "Method", by = c("N","Endpoints", "Meanvect", "\u03BC"),
  ref = "Cohen's D"
)
Nestedplot_bias <- autoplot(Nestedplot_bias, type = "nlp", stats = "thetamean") 

Nestedplot_bias+ 
  ylab("Relative Bias") + 
  theme_bw() 

