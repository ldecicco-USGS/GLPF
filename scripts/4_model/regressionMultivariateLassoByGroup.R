# Next steps: 
# consider refitting with lm() and examining slopes by event
# consider automated way to keep track of slopes for individual events within group regressions
# -maybe use median regression for this.

# Heuristic overlap analysis: decision tree


library(glmnet)
library(dplyr)
library(RColorBrewer)
library(parallel)
library(doParallel)

#set up parallel cores for cv.glmnet
# Calculate the number of cores
no_cores <- detectCores() - 1

# Register  parallel backend
registerDoParallel(no_cores)


#setwd("D:/SRCData/Git/GLPF")
source("na.info.R")

df.orig <- readRDS("./cached_data/8_process_new_categories/rds/summary_noWW_noQA.rds")

#df.orig <- summaryDF
df <- df.orig
response <- c("lachno","bacHum")

df <- df[-which(is.na(df$lachno)),]

beginIV <- "Sag240_255"
endIV <- "rBS44_S45_BF"

begin <- which(names(df)==beginIV)
end <- which(names(df)==endIV)

IVs <- names(df)[begin:end]

na.info.list <- na.info(df[,-dim(df)[2]],first.col = beginIV)
rmRows <- unique(c(which(df$CAGRnumber %in% na.info.list$na.rows),
                   na.info.list$nan.rows,
                   na.info.list$inf.rows))
rmCols <- unique(which(names(df) %in% c(na.info.list$na.cols.partial,
                                        na.info.list$nan.cols,
                                        na.info.list$inf.cols)))
dfrmCols <- df[,-rmCols]
dfRmRows <- df[rmRows,]
df <- df[,-rmCols]

beginIV <- "Sag240_255"
endIV <- "rBS44_S45_BF"
begin <- which(names(df)==beginIV)
end <- which(names(df)==endIV)
IVs <- names(df)[begin:end]

groupFreq <- table(df$eventGroup2)

groups <- names(groupFreq)[which(groupFreq>21)]

mg.List <- list()
mg.cv.List <- list()

lambdaType <- "lambda.min"
filenm <- "MVLassoByGroupLminFull2.pdf"
pdf(filenm)
modelCoefList <- list()
for(i in 1:length(groups)){
  
  subdf <- df[which(df$eventGroup2==groups[i]),]
  IVs <- names(df)[begin:end]
  
  #subdf <- df
  foldID <- as.numeric(as.factor(subdf$eventNum))
  events <- unique(subdf$eventNum)
  
# #  Add events as separate dichotomous IVs
#   if(length(events)>1){
#     # eventDF <- data.frame(E1 = ifelse(subdf$eventNum == events[2],1,0))
#     for(j in 2:length(events)){
#       if(j==2)eventDF <- as.data.frame(ifelse(subdf$eventNum == events[j],1,0))
#       else eventDF <- cbind(eventDF,ifelse(subdf$eventNum == events[j],1,0))
#     }
#     names(eventDF) <- events[-1]
#     subdf <- cbind(subdf,eventDF)
#     IVs <- c(IVs,names(eventDF))
#   }
  
  y <- log10(as.matrix(subdf[,response]))
  x <- as.matrix(subdf[IVs])
  
  #If more than 2 events included in group, use event as fold ID, otherwise, use 5-fold XV
  if(length(unique(foldID))>2){
    mg.cv <- cv.glmnet(x=x, y=y,family="mgaussian",alpha=1,foldid = foldID,parallel = TRUE)
    mg <- glmnet(x=x, y=y,family="mgaussian", alpha=1)
  }else{
    mg.cv <- cv.glmnet(x=x, y=y,family="mgaussian",alpha=1,nfolds=5,parallel = TRUE)
    mg <- glmnet(x=x, y=y,family="mgaussian", alpha=1)
  }
  
  
  #Extract Coefficients from cv-determined model using lambda.1se
  if(lambdaType == "lambda.1se"){
    Coefficients <- coef(mg, s = mg.cv$lambda.1se)
    Active.Index <- which(Coefficients[[1]] != 0)
    Active.Coefficients <- Coefficients[[1]][Active.Index];Active.Coefficients
    Active.Coef.names <- row.names(Coefficients[[1]])[Active.Index];Active.Coef.names
  }
  
  
  #Extract Coefficients from cv-determined model using lambda.min
  if(lambdaType == "lambda.min"){
    Coefficients <- coef(mg, s = mg.cv$lambda.min)
    Active.Index <- which(Coefficients[[1]] != 0)
    Active.Coefficients <- Coefficients[[1]][Active.Index];Active.Coefficients
    Active.Coef.names <- row.names(Coefficients[[1]])[Active.Index];Active.Coef.names
  }
  
  modelCoefList[[i]] <- Active.Coef.names[-1]
  

  #Plot cross validated errors and other model results
  plot(mg.cv)
  
  predictions <- predict(mg.cv,newx=as.matrix(subdf[,IVs]),s=lambdaType,type = "response")
  
  plotpch <- 20
  colorOptions <- brewer.pal(9, "Set1")
  
  plotCol <- colorOptions[1:length(events)]
  names(plotCol) <- events
  plotcolors <- plotCol[subdf$eventNum]
  

  par(mfcol=c(2,1),mar=c(3,4,3,1),oma=c(0,2,0,4))
  #Plot Lachno
  plot(subdf[,response[1]],predictions[,1,1],col=plotcolors,pch=plotpch,log='x',xlab="",ylab="")
  mtext(response[1],line=1)
  mtext(paste(Active.Coef.names[-1],collapse=' + '),cex=0.7)
  mtext(groups[i],line=2,font=2)
  
  
  #Plot bacHum
  plot(subdf[,response[2]],predictions[,2,1],col=plotcolors,pch=plotpch,log='x',xlab="",ylab="")
  mtext("Predicted",side=2,line=-2,font=2,xpd=NA,outer=TRUE)
  mtext("Observed",side=1,line=2,font=2)
  mtext(response[2],line=1)
  
  legendNames <- names(plotCol)
  legend('bottomright',legend = legendNames,col=plotCol,pch=plotpch,inset = c(-0.15,0),bty = "n",xpd=NA)
 
  
  # calibrate Tobit regression and plot
  library(survival)

  IVs <- Active.Coef.names[-1]
  response <- response
  LOQ <- 225

  ## Compute survival coefficients for Lachno regression ##
  if(length(IVs) > 0){
  y <- Surv(log10(subdf[,response[1]]), subdf[,response[1]]>LOQ, type="left")
  #dfPredStd <- as.data.frame(scale(dfPred[,IVs]))
  form <- formula(paste('y ~',paste(IVs,collapse=' + ')))
  msurvStd <- survreg(form,data=subdf,dist='weibull')
  summary(msurvStd)
  predictions <- predict(msurvStd,newdata = subdf)
  
  par(mfcol=c(2,1),mar=c(3,4,3,1),oma=c(0,2,0,4))
  #Plot Lachno
  plot(subdf[,response[1]],predictions,col=plotcolors,pch=plotpch,log='x',xlab="",ylab="")
  mtext(paste(response[1],"Survival"),line=1)
  mtext(paste(Active.Coef.names[-1],collapse=' + '),cex=0.7)
  mtext(groups[i],line=2,font=2)
  abline(h=4,v=10000,col="blue",lty=2)
  
  
  ## Compute survival coefficients for Lachno regression ##
  y <- Surv(log10(subdf[,response[2]]), subdf[,response[2]]>LOQ, type="left")
  #dfPredStd <- as.data.frame(scale(dfPred[,IVs]))
  form <- formula(paste('y ~',paste(IVs,collapse=' + ')))
  msurvStd <- survreg(form,data=subdf,dist='weibull')
  summary(msurvStd)
  predictions <- predict(msurvStd,newdata = subdf)
  
  #Plot BacHum
  plot(subdf[,response[2]],predictions,col=plotcolors,pch=plotpch,log='x',xlab="",ylab="")
  mtext(paste(response[2],"Survival"),line=1)
  mtext(groups[i],line=2,font=2)
  abline(h=4,v=10000,col="blue",lty=2)
  }else{
    par(mfrow=c(1,1))
    plot(1:10,1:10,xaxt="n",yaxt="n",ylab="",xlab="",pch="")
  text(5,5,"No Lasso Variables")
  }
  
  mg.List[[i]] <- mg
  mg.cv.List[[i]] <- mg.cv
  
#----------------------------------------------------------------------------------
  # calibrate model with all data from individual group model
  
  subdf <- df
  y <- log10(as.matrix(subdf[,response]))
  x <- as.matrix(subdf[Active.Coef.names[-1]])
  foldIDs <- as.numeric(as.factor(subdf$eventNum))
  
  which(table(foldIDs)<10)

  if(length(modelCoefList[[i]]) > 1) {

    mg.cv <- cv.glmnet(x=x, y=y,family="mgaussian",alpha=1,foldid = foldIDs,parallel = TRUE)
    mg <- glmnet(x=x, y=y,family="mgaussian", alpha=1)

    


  #Extract Coefficients from cv-determined model using lambda.1se
  if(lambdaType == "lambda.1se"){
    Coefficients <- coef(mg, s = mg.cv$lambda.1se)
    Active.Index <- which(Coefficients[[1]] != 0)
    Active.Coefficients <- Coefficients[[1]][Active.Index];Active.Coefficients
    Active.Coef.names <- row.names(Coefficients[[1]])[Active.Index];Active.Coef.names
  }


  #Extract Coefficients from cv-determined model using lambda.min
  if(lambdaType == "lambda.min"){
    Coefficients <- coef(mg, s = mg.cv$lambda.min)
    Active.Index <- which(Coefficients[[1]] != 0)
    Active.Coefficients <- Coefficients[[1]][Active.Index];Active.Coefficients
    Active.Coef.names <- row.names(Coefficients[[1]])[Active.Index];Active.Coef.names
}
    #Plot cross validated errors and other model results
    plot(mg.cv)
    
    predictions <- predict(mg.cv,newx=x,s=lambdaType,type = "response")
    
    plotpch <- 20
    colorOptions <- brewer.pal(9, "Set1")
    
    plotCol <- colorOptions[1:length(events)]
    names(plotCol) <- events
    plotcolors <- "grey"
#      plotCol[subdf$eventNum]
    
    
    par(mfcol=c(2,1),mar=c(3,4,3,1),oma=c(0,2,0,4))
    #Plot Lachno
    plot(subdf[,response[1]],predictions[,1,1],col=plotcolors,pch=plotpch,log='x',xlab="",ylab="")
    mtext(response[1],line=1)
    mtext(paste(Active.Coef.names[-1],collapse=' + '),cex=0.7)
    mtext(groups[i],line=2,font=2)
    
    
    #Plot bacHum
    plot(subdf[,response[2]],predictions[,2,1],col=plotcolors,pch=plotpch,log='x',xlab="",ylab="")
    mtext("Predicted",side=2,line=-2,font=2,xpd=NA,outer=TRUE)
    mtext("Observed",side=1,line=2,font=2)
    mtext(response[2],line=1)
    
    legendNames <- names(plotCol)
    legend('bottomright',legend = legendNames,col=plotCol,pch=plotpch,inset = c(-0.15,0),bty = "n",xpd=NA)
    
}
  
}
dev.off()
shell.exec(filenm)

###------------------------------------------------------------------------
#Plot with all data by group and then by event

plotAll <- TRUE

lambdaType <- "lambda.min"
if(plotAll) {
  filenm <- "GroupLassoByEvent.pdf"
  pdf(filenm)
  subdf <- df
  y <- log10(as.matrix(subdf[,response]))
  x <- as.matrix(subdf[,IVs])
  
  events <- unique(subdf$eventNum)

  for(i in 1:length(groups)){
  
  mg <- mg.List[[i]]
  mg.cv <- mg.cv.List[[i]]
  
  

    #Extract Coefficients from cv-determined model using lambda.1se
    if(lambdaType == "lambda.1se"){
      Coefficients <- coef(mg, s = mg.cv$lambda.1se)
      Active.Index <- which(Coefficients[[1]] != 0)
      Active.Coefficients <- Coefficients[[1]][Active.Index];Active.Coefficients
      Active.Coef.names <- row.names(Coefficients[[1]])[Active.Index];Active.Coef.names
    }
    
    
    #Extract Coefficients from cv-determined model using lambda.min
    if(lambdaType == "lambda.min"){
      Coefficients <- coef(mg, s = mg.cv$lambda.min)
      Active.Index <- which(Coefficients[[1]] != 0)
      Active.Coefficients <- Coefficients[[1]][Active.Index];Active.Coefficients
      Active.Coef.names <- row.names(Coefficients[[1]])[Active.Index];Active.Coef.names
    }

    predictions <- predict(mg.cv,newx=x,s=lambdaType,type = "response")
    
    for(j in 1:length(events)){
      event <- events[j]
    plotpch <- 20
    # colorOptions <- brewer.pal(9, "Set1")
    # 
    # plotCol <- colorOptions[1:length(events)]
    # names(plotCol) <- events
    plotcolors <- "grey"
    eventcolor <- ifelse(subdf$eventNum==event,"blue",NA)
    #      plotCol[subdf$eventNum]
    
    ylim <- range(predictions[,1,1])
    ylim[1] <- ifelse(ylim[1] < 0,0,ylim[1])
    ylim[2] <- ifelse(ylim[2] > 8,8,ylim[2])
    
    par(mfcol=c(2,1),mar=c(3,4,3,1),oma=c(0,2,0,4))
    #Plot Lachno
    plot(subdf[,response[1]],predictions[,1,1],col=plotcolors,pch=plotpch,log='x',
         xlab="",ylab="",ylim=ylim)
    points(subdf[,response[1]],predictions[,1,1],col=eventcolor,pch=plotpch)
    mtext(response[1],line=1)
    mtext(paste(Active.Coef.names[-1],collapse=' + '),cex=0.7)
    mtext(paste(groups[i],";",event),line=2,font=2)
    
    
    #Plot bacHum
    ylim <- range(predictions[,2,1])
    ylim[1] <- ifelse(ylim[1] < 0,0,ylim[1])
    ylim[2] <- ifelse(ylim[2] > 8,8,ylim[2])
    
    plot(subdf[,response[2]],predictions[,2,1],col=plotcolors,pch=plotpch,log='x',
         xlab="",ylab="",ylim = ylim)
    points(subdf[,response[2]],predictions[,2,1],col=eventcolor,pch=plotpch)
    mtext("Predicted",side=2,line=-2,font=2,xpd=NA,outer=TRUE)
    mtext("Observed",side=1,line=2,font=2)
    mtext(response[2],line=1)
    
    legendNames <- names(plotCol)
    legend('bottomright',legend = legendNames,col=plotCol,pch=plotpch,inset = c(-0.15,0),bty = "n",xpd=NA)
    
    }
  }
}
dev.off()
shell.exec(filenm)
