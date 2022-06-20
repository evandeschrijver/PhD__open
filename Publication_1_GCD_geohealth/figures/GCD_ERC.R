
pdf("C:/Users/ed19w187/Desktop/PhD/PhD/Publication_GCD/ERC/final/all3.pdf", width=6,height=6)
for(i in seq(length(namesreg_full))) {
  
  data <- fulldataset[[i]]
  
  ### GCD 1
  
  argvar <- list(x=data$Insitu,fun=varfun,
                 knots=quantile(data$Insitu,varper/100,na.rm=T),
                 Bound=range(data$Insitu,na.rm=T))
  if(!is.null(vardegree)) argvar$degree <- vardegree
  bvar <- do.call(onebasis,argvar)
  pred2 <- crosspred(bvar,coef=coefallinsitu[i,],
                     vcov=vcovallinsitu[[i]],
                     model.link="log",cen=mmt_list[i,1],by=0.1, from=quantile(data$Insitu,0.001), to=quantile(data$Insitu,0.999))
  res <- data.frame(cbind(pred2$predvar, pred2$allRRfit, pred2$allRRlow, pred2$allRRhigh))
  res$GCD <- "Insitu"
  names(res) <- c("temp", "rr", "rrlow", "rrhigh", "GCD")
  
  temp_GCD <- as.data.frame(data$Insitu)
  temp_GCD$GCD <- "Insitu"
  colnames(temp_GCD) <- c("Temp", "GCD")
  
  
  #### GCD 2
  argvar_GCD2 <- list(x=data$Haduk5,fun=varfun,
                      knots=quantile(data$Haduk5,varper/100,na.rm=T),
                      Bound=range(data$Haduk5,na.rm=T))
  if(!is.null(vardegree)) argvar_GCD2$degree <- vardegree
  bvar_GCD2 <- do.call(onebasis,argvar_GCD2)
  pred2_GCD2 <- crosspred(bvar_GCD2,coef=coefallhaduk5[i,],
                          vcov=vcovallhaduk5[[i]],
                          model.link="log",cen=mmt_list[i,2],by=0.1, from=quantile(data$Haduk5,0.001), to=quantile(data$Haduk5,0.999))
  res_GCD2 <- data.frame(cbind(pred2_GCD2$predvar, pred2_GCD2$allRRfit, pred2_GCD2$allRRlow, pred2_GCD2$allRRhigh))
  
  res_GCD2$GCD <- "Haduk5"
  names(res_GCD2) <- c("temp", "rr", "rrlow", "rrhigh", "GCD")
  
  temp_GCD2 <- as.data.frame(data$Haduk5)
  temp_GCD2$GCD <- "Haduk5"
  colnames(temp_GCD2) <- c("Temp", "GCD")
  
  
  
  #### GCD 3
  argvar_GCD3 <- list(x=data$ERA5,fun=varfun,
                      knots=quantile(data$ERA5,varper/100,na.rm=T),
                      Bound=range(data$ERA5,na.rm=T))
  if(!is.null(vardegree)) argvar_GCD3$degree <- vardegree
  bvar_GCD3 <- do.call(onebasis,argvar_GCD3)
  pred3_GCD3 <- crosspred(bvar_GCD3,coef=coefallERA5[i,],
                          vcov=vcovallERA5[[i]],
                          model.link="log",cen=mmt_list[i,4],by=0.1, from=quantile(data$ERA5,0.001), to=quantile(data$ERA5,0.999))
  res_GCD3 <- data.frame(cbind(pred3_GCD3$predvar, pred3_GCD3$allRRfit, pred3_GCD3$allRRlow, pred3_GCD3$allRRhigh))
  
  res_GCD3$GCD <- "ERA5"
  names(res_GCD3) <- c("temp", "rr", "rrlow", "rrhigh", "GCD")
  temp_GCD3 <- as.data.frame(data$ERA5)
  temp_GCD3$GCD <- "ERA5"
  colnames(temp_GCD3) <- c("Temp", "GCD")
  
  ######
  tempdistt <- rbind(temp_GCD,temp_GCD2, temp_GCD3)
  tempdistt$col[tempdistt$GCD == "Insitu"] <- col_mon
  tempdistt$col[tempdistt$GCD == "Haduk5"] <- col_haduk
  tempdistt$col[tempdistt$GCD == "ERA5"] <- col_ERA
  
  

  names(res_GCD2) <- c("temp", "rr", "rrlow", "rrhigh", "GCD")

 
  names(res_GCD3) <- c("temp", "rr", "rrlow", "rrhigh", "GCD")
  res <- rbind(res,res_GCD2, res_GCD3)
  
  res <- res %>% mutate(GCD = factor(res$GCD, levels=c( "ERA5", "Insitu", "Haduk5"))) 
  
  listerplot[[i]] <- ggplot(data=res, aes(x=temp, y=rr, group=GCD)) +
    geom_line(aes(colour=GCD), size=1) + scale_color_manual(values=c(col_ERA, col_mon, col_haduk)) +
    geom_ribbon(aes(ymin=rrlow, ymax=rrhigh), linetype=1, alpha=0.1)+
    ylim(c(0.7, 2.5))+ xlim(-8,30) +
    labs(y="Relative Risk", x="Temperature (?C)")+
    geom_vline(xintercept=quantile(data$Insitu,0.99, na.rm=T),
               linetype="dashed", color = col_mon) + 
    geom_vline(xintercept=quantile(data$Insitu,0.01, na.rm=T),
               linetype="dashed", color = col_mon) +
    
    geom_vline(xintercept=quantile(data$Haduk5,0.99, na.rm=T),
               linetype="dashed", color = col_haduk) + 
    geom_vline(xintercept=quantile(data$Haduk5,0.01, na.rm=T),
               linetype="dashed", color = col_haduk) +
    geom_vline(xintercept=quantile(data$ERA5,0.99, na.rm=T),
               linetype="dashed", color = col_ERA) + 
    geom_vline(xintercept=quantile(data$ERA5,0.01, na.rm=T),
               linetype="dashed", color = col_ERA) +
    
    geom_hline(yintercept=1,
               linetype="solid", color = "grey40") + 
    theme_minimal() +
    theme(legend.position = "none",
      axis.line = element_line(size =0.5),
      axis.text = element_text(size=10), 
      axis.title = element_text(size=10),
      axis.title.x = element_blank(),
      axis.ticks.x=element_blank(),
      axis.text.x = element_blank(),
      axis.line.x = element_blank(),
      plot.title = element_text(size=18, face="bold"),
      plot.subtitle = element_text(size=15, face="italic"),
      axis.ticks = element_line(size=0.5 ),
      axis.ticks.length = unit(4, "pt")) + theme(plot.margin = unit(c(1,1,0,1), "lines"))
  
  tempdistt$GCD <- factor(tempdistt$GCD, levels = c("ERA5", "Insitu", "Haduk5"))
  

  
  tempplot[[i]] <- tempdistt %>% ggplot(aes(x=Temp, group=GCD, fill=GCD)) + geom_density(alpha=0.4)+
    scale_fill_manual(values=c(col_ERA, col_mon, col_haduk), labels = c("ERA5", "Monitor", "Haduk")) +
    #scale_x_continuous(limits = c(min(res$temp), max(res$temp))) +
    labs(y="Density", x="Temperature (?C)")+ ylim(0,0.084) + xlim(-8,30) +
    
    geom_vline(xintercept=quantile(data$Insitu,0.99, na.rm=T),
               linetype="dashed", color = col_mon) + 
    geom_vline(xintercept=quantile(data$Insitu,0.01, na.rm=T),
               linetype="dashed", color = col_mon) +
    geom_vline(xintercept=quantile(data$Haduk5,0.99, na.rm=T),
               linetype="dashed", color = col_haduk) + 
    geom_vline(xintercept=quantile(data$Haduk5,0.01, na.rm=T),
               linetype="dashed", color = col_haduk) +
    geom_vline(xintercept=quantile(data$ERA5,0.99, na.rm=T),
               linetype="dashed", color = col_ERA) + 
    geom_vline(xintercept=quantile(data$ERA5,0.01, na.rm=T),
               linetype="dashed", color = col_ERA) +
    
    theme_minimal() + 
    theme(legend.position = "none",
      axis.line = element_line(size =0.5),
      axis.text = element_text(size=8), 
      axis.title = element_text(size=10),
      #legend.position="none",
      plot.title = element_text(size=18, face="bold"),
      plot.subtitle = element_text(size=15, face="italic"),
      axis.ticks = element_line(size=0.5 ),
      axis.ticks.length = unit(4, "pt")) + theme(plot.margin = unit(c(0,1,1,1), "lines"))
  
  
  
  finalplot_tot[[i]] <- egg::ggarrange(listerplot[[i]], tempplot[[i]], heights = c(1.25,1))
  
}
dev.off()



