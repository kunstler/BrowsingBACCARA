##############
### FUNCTIONS FOR ANALYSIS

Compute_Surv_CI <- function(var_bin){
  print(length(var_bin))
  res <- glm(var_bin~ 1, data = data.frame(var_bin = var_bin),
             family = binomial(link = "logit"))
  preds <- predict(res,newdata = data.frame(var_bin= 1),
                   se.fit = TRUE, type = "link")
  critval <- 1.96 ## approx 95% CI
  upr <- preds$fit + (critval * preds$se.fit)
  lwr <- preds$fit - (critval * preds$se.fit)
  fit <- preds$fit
  fit2 <- res$family$linkinv(fit)
  upr2 <- res$family$linkinv(upr)
  lwr2 <- res$family$linkinv(lwr)
  return(c(mean = fit2, lwr = lwr2, upr = upr2))
}

AICmod1<-function(sp, bdd){
  aic1<-(glmer(B~1+(1|Site/Z.Classe),family = 'binomial',data=bdd[bdd$Species==sp,]))
  aic2<-(glmer(B~scale(Tm_hiv)+(1|Site/Z.Classe),family = 'binomial',data=bdd[bdd$Species==sp,]))
  aic3<-(glmer(B~scale(HT)+(1|Site/Z.Classe),family = 'binomial',data=bdd[bdd$Species==sp,]))
  aic4<-(glmer(B~scale(Tm_hiv)+scale(HT)+(1|Site/Z.Classe),family = 'binomial',data=bdd[bdd$Species==sp,]))
  aic5<-(glmer(B~scale(Tm_hiv)*scale(HT)+(1|Site/Z.Classe),family = 'binomial',data=bdd[bdd$Species==sp,]))
  res_list<-list(aic1,aic2,aic3,aic4,aic5)
  return(res_list)
}

### CI on browsing rate
fun_site_browsing_rate <- function(site_alt, df){
  tryCatch(Compute_Surv_CI(df[df$site_alt == site_alt,]$B),
           error =function(e) c(mean = NA, lwr = NA, upr = NA))
}


Brows_Prob_CI <- function(bdd){
  vec_site_alt <- unique(bdd$site_alt)
  list_res_ci <- vector("list")
  for(i in 1:4){
    list_res_ci[[i]] <- as.data.frame(t(sapply(vec_site_alt,
                                               fun_site_browsing_rate,
                                               bdd[bdd$Species == vec_sp[i], ])))
    list_res_ci[[i]]$sp <- vec_sp[i]
    list_res_ci[[i]]$site_alt <- vec_site_alt
    names(list_res_ci[[i]]) <- c("mean", "lwr", "upr", "sp", "site_alt")
  }
  df_ci <- as.data.frame(do.call("rbind", list_res_ci))
  df_ci$Tm_hiv<-bdd$Tm_hiv[match(df_ci$site_alt,bdd$site_alt)]
  
  return(df_ci)
}

Pred_CI <- function(res, newdat, CI_prob = 0.8){ # interval at 0.8 %
  z <- qnorm(CI_prob/2)
  
  mm<-model.matrix(~Tm_hiv_CR*HT_CR,newdat)
  y<-mm%*%fixef(res)
  pvar1 <- diag(mm %*% tcrossprod(vcov(res),mm))
  newdat2 <- data.frame(
    y=invlogit(y),
    plo = invlogit(y-1.96*sqrt(pvar1)),
    phi = invlogit(y+1.96*sqrt(pvar1))
  )
  newdat <-  cbind(newdat, newdat2)
  return(newdat)
}

Pred_Brows_Prob <- function(fglm, var_l, var_h, HT_val, df, sp, var = "Tm_hiv"){
  pred <- seq(from= var_l, to = var_h, length.out = 100)
  pred_CR <- (pred - mean(df[[var]]))/sd(df[[var]])
  HT <-  HT_val
  HT_CR <- (HT - mean(df$HT))/sd(df$HT)
  data_pred <- expand.grid(var= pred, HT = HT,
                           Site = df$Site[1000], Z_F = df$Z_F[1000], Species = sp)
  data_pred$var_CR <- (data_pred$var- mean(df[[var]]))/sd(df[[var]])
  data_pred$HT_CR <- (data_pred$HT - mean(df$HT))/sd(df$HT)
  names(data_pred)[c(1, 6)] <- c(var, paste0(var, "_CR"))
  df_out <- Pred_CI(fglm, data_pred)
  return(df_out)
}

