library(dplyr)
library (pROC)
library(ipw)
library(tidyr)
library(fabricatr)
library (survival)
library(ggfortify)
library(ggsci)
library(ggplot2)
library(ggfortify)
library(survminer)
library(ggpubr)
library(causalCmprsk)
library(randomForestSRC)
library(ggRandomForests)


df <- read.csv("Documents/FHS/FHS_Network_features_dementia_risk_dataforsurvivalanalysis_wave7.csv")

df$DEGREE_CENTRALITY = df$DEGREE_CENTRALITY*1000 
df$BETWEENNESS_CENTRALITY = df$BETWEENNESS_CENTRALITY*1000 
#df$Diversity_centered = df$Diversity - min((df$Diversity))
#diversity measure undefined
df <- within(df, CSD.D_loneliness_item[CSD.D_loneliness_item <= 2] <- 0)
df <- within(df, CSD.D_loneliness_item[CSD.D_loneliness_item >= 3] <- 1)

network_feature_list = list(df$DEGREE, df$DEGREE_CENTRALITY, df$CONSTRAINT, df$EIGENVECTOR, df$EFFECTIVE_SIZE, df$KATZ, df$BETWEENNESS_CENTRALITY, df$CLUSTCOEF, df$Diversity_centered, df$Diversity, df$SEX_SIMPSON, df$SEX_SHANNON, df$RELTYPE_SIMPSON, df$RELTYPE_SHANNON, df$ALTERTYPE_SIMPSON, df$ALTERTYPE_SHANNON)
network_feature_list = list('CLUSTCOEF', 'DEGREE','CONSTRAINT','EIGENVECTOR','EFFECTIVE_SIZE','KATZ','BETWEENNESS_CENTRALITY', 'SEX_SIMPSON', 'SEX_SHANNON', 'RELTYPE_SIMPSON', 'RELTYPE_SHANNON', 'ALTERTYPE_SIMPSON', 'ALTERTYPE_SHANNON')

for (network_feature in network_feature_list) {
  feature_name = network_feature
  network_feature = df[[network_feature]]

  # define groups splitting by the median
  network_feature_M <- split_quantile(x = network_feature, type = 2) 
  network_feature_M <- factor(network_feature_M, levels = c(1, 2))
  
  # define binary variable for IPW weight calculation (need to be numeric 0, 1)
  median_network_feature <- median(network_feature) 
  network_feature_level_IPW <- ifelse(network_feature>median_network_feature, 1, 
                                      0) # 1=high, 0=low
  median_network_feature <- median(network_feature) 
  network_feature_level_IPW <- ifelse(network_feature>median_network_feature, 1, 
                                      0) # 1=high, 0=low
  df$network_feature_level_IPW = network_feature_level_IPW
  
  df$Event = df$DEM_STATUS
  df$time.to.dementia = df$time.to.dementia/365
  df$time.to.event = df$time.to.dementia
  
  df_COX = df
  df_COX$network_feature = network_feature
  df_COX$network_feature_M = network_feature_M
  
  df_COX <- within(df_COX, Event[Event == 0 & Death_status_at_followup == 1] <- 0)
  df_COX <- within(df_COX, Event[Event == 0 & Death_status_at_followup == 0] <- 0)
  df_COX <- within(df_COX, Event[Event == 1 & Death_status_at_followup == 1] <- 1)
  df_COX <- within(df_COX, Event[Event == 1 & Death_status_at_followup == 0] <- 1)

  
  df_COX$p.score <- glm(network_feature_level_IPW ~ AGE7 + SEX 
                        + somecoll_education 
                        #+ SBP7 + DBP7
                        #+ BMI7  
                        #+ CURRSMK7 + CALC_LDL7 + HDL7 + TC7 + TRIG7
                        #+ APOE
                        #+ CURR_DIAB7 + CVD_HISTORY
                        #+ CES.D_score
                        , 
                        data = df_COX, 
                        family = "binomial")$fitted.values
  
  df_COX$ate.weights <- with(df_COX, network_feature_level_IPW * 1/p.score + (1-network_feature_level_IPW)* 1/(1-p.score))
  
  ## Report logistic model used for weight calculation
  a <- glm(network_feature_level_IPW ~ AGE7 + SEX 
           + somecoll_education 
           #+ SBP7 + DBP7
           #+ BMI7  
           #+ CURRSMK7 + CALC_LDL7 + HDL7 + TC7 + TRIG7
           #+ APOE
           #+ CURR_DIAB7 + CVD_HISTORY
           #+ CES.D_score
           ,
           data = df_COX, 
           family = "binomial")
  print(summary(a))

  fit <- coxph(Surv(time.to.event, Event) ~ AGE7 + SEX 
               + somecoll_education 
               #+ SBP7 + DBP7
               #+ BMI7  
               #+ CURRSMK7 + CALC_LDL7 + HDL7 + TC7 + TRIG7
               #+ APOE
               #+ CURR_DIAB7 + CVD_HISTORY
               #+ CES.D_score
               + network_feature, 
               data=df_COX, 
               weights = ate.weights)
  summary(fit)
  write.csv(summary(fit)$conf.int,paste(feature_name,'_confint_cox_weights.csv',sep=''))
  write.csv(summary(fit)$coefficients, paste(feature_name,'coef_cox_weights.csv',sep=''))  
  
  output <- tidy(fit)
  write.csv(output, paste(feature_name,'_cox_weights.csv',sep=''))  
  
  df_CIF = df
  df_CIF <- within(df_CIF, Event[Event == 1 & Death_status_at_followup == 1] <- 1)
  df_CIF <- within(df_CIF, Event[Event == 0 & Death_status_at_followup == 1] <- 2)
  df_CIF <- within(df_CIF, Event[Event == 1 & Death_status_at_followup == 0] <- 1)
  df_CIF <- within(df_CIF, Event[Event == 0 & Death_status_at_followup == 0] <- 0)
  
  covs.names = c("AGE7","SEX", "somecoll_education"#, 
                 #"BMI7", "SBP7", "DBP7",         
                 #"CURRSMK7", "CALC_LDL7", "HDL7", "TC7",
                 #"TRIG7", "APOE", 
                 #"CES.D_score",
                 #"CURR_DIAB7", "CVD_HISTORY"
  )

    # A small number of bootstrap replications (nbs.rep=80) can be chosen for illustration purposes.  
  df_CIF$network_feature_level_IPW = network_feature_level_IPW
  
  res.overlap1 <- fit.nonpar(df_CIF, X="time.to.event", E="Event", 
                             A="network_feature_level_IPW", # treatment arm (e.g. network feature groups in our case)
                             C=covs.names,
                             wtype="overlap", cens=0, conf.level=0.95, bs=TRUE, nbs.rep=80, seed=17, parallel = FALSE) 

  get.CIF.RMT_1 <- function(res) # this will only plot the below and above median
  {
    df_CIF.CIF.RMT <- rbind(data.frame(time=res$time, Event_TRT=1, ### Group No. 1 for plot
                                       # trt.1 = network feature, level high = 1
                                       # Ev= 1, outcome = 1 (dementia)
                                       CIF=res$trt.1[[paste("Ev=", 1, sep="")]]$CIF, 
                                       RMT=res$trt.1[[paste("Ev=", 1, sep="")]]$RMT,
                                       CumHaz=res$trt.1[[paste("Ev=", 1, sep="")]]$CumHaz,
                                       CIL.CIF=res$trt.1[[paste("Ev=", 1, sep="")]]$CIF.CI.L,
                                       CIU.CIF=res$trt.1[[paste("Ev=", 1, sep="")]]$CIF.CI.U,
                                       CIL.RMT=res$trt.1[[paste("Ev=", 1, sep="")]]$RMT.CI.L,
                                       CIU.RMT=res$trt.1[[paste("Ev=", 1, sep="")]]$RMT.CI.U,
                                       CIL.CumHaz=res$trt.1[[paste("Ev=", 1, sep="")]]$CumHaz.CI.L,
                                       CIU.CumHaz=res$trt.1[[paste("Ev=", 1, sep="")]]$CumHaz.CI.U),
                            data.frame(time=res$time, Event_TRT=2, ### Group No. 2 for plot
                                       # trt.0 = network feature, level low = 0
                                       # Ev= 1, outcome = 1 (dementia)
                                       CIF=res$trt.0[[paste("Ev=", 1, sep="")]]$CIF, 
                                       RMT=res$trt.0[[paste("Ev=", 1, sep="")]]$RMT,
                                       CumHaz=res$trt.0[[paste("Ev=", 1, sep="")]]$CumHaz,
                                       CIL.CIF=res$trt.0[[paste("Ev=", 1, sep="")]]$CIF.CI.L,
                                       CIU.CIF=res$trt.0[[paste("Ev=", 1, sep="")]]$CIF.CI.U,
                                       CIL.RMT=res$trt.0[[paste("Ev=", 1, sep="")]]$RMT.CI.L,
                                       CIU.RMT=res$trt.0[[paste("Ev=", 1, sep="")]]$RMT.CI.U,
                                       CIL.CumHaz=res$trt.0[[paste("Ev=", 1, sep="")]]$CumHaz.CI.L,
                                       CIU.CumHaz=res$trt.0[[paste("Ev=", 1, sep="")]]$CumHaz.CI.U)#,
    )
    df_CIF.CIF.RMT$Event_TRT <- factor(df_CIF.CIF.RMT$Event_TRT)
    levels(df_CIF.CIF.RMT$Event_TRT) <- c("Dementia network_feat high", "Dementia network_feat low")
    return(df_CIF.CIF.RMT)
  }
  df_CIF.CIF.RMT1 <- get.CIF.RMT_1(res.overlap1)
  combined <- rbind(df_CIF.CIF.RMT1)
  
  write.csv(combined, paste(feature_name,'_CIF.csv',sep=''))  
  

  p <- ggplot(combined, aes(x=time, y=CIF,
                            color=Event_TRT, 
                            fill=Event_TRT, 
                            shape=Event_TRT
  )) +
    geom_step(size=2.5) + ggtitle("Cumulative Incidence") + 
    geom_ribbon(aes(ymin=CIL.CIF, ymax=CIU.CIF), alpha=0.2, stat="identity") + # stat="stepribbon") +
    scale_fill_npg() +
    scale_color_npg() 
  #scale_colour_manual(values=c('1'="#000066",'2'="#663399",'3'="#339999"))
  #scale_x_continuous(breaks = round(seq(min(combined$time), max(combined$time), by = 5),5)) +# No breaks
  #xlim(0, 20)
  p <- p +  #geom_vline(xintercept=30, linetype="dashed") +
    xlab("time (years)") + ylab("Probability of event by time t") +
    theme(axis.text.x = element_text(face="bold", angle=45),
          axis.text.y = element_text(face="bold"), plot.title = element_text(hjust = 0.5))+
    theme(legend.position = c(0.2, 0.8),
          legend.background=element_rect(fill="transparent"),
          panel.grid.minor = element_line(size = .5,colour = "gray92"),
          panel.grid.major = element_line(size = .5,colour = "#C0C0C0"))
  # png(paste(feature_name,'_CIF_plot.png',sep='')) 
  p + theme_classic()
  
  ggsave(paste(feature_name,'_CIF_plot.png',sep=''), plot = p + theme_classic(), width = 30, height = 20, units = "cm")
  
  # dev.off()
  # Close the png file
     
  # plot adjusted models for each outcome/competing-event
  #plot.competing.risk(rf) 


  
}

