library(USGSHydroOpt)
library(dplyr)
library(ggplot2)

cached.path <- "cached_data"
base.name <- "_noQA"
  
plotOpticalSummaries <- function(base.name, cached.path){
  
  dfall <- readRDS(file.path(cached.path,"5_process_filterData","rds",paste0("summary",base.name,".rds")))
  dfabs <- readRDS(file.path(cached.path,"5_process_filterData","rds",paste0("dfabs",base.name,".rds")))
  EEMs3D <- readRDS(file.path(cached.path,"7_process_summarize_optics","rds",paste0("EEMs3D",base.name,".rds")))
  
  ##############################################################################################
  ##############################################################################################
  # Plot EEMs
  a <- EEMs3D
  fieldIDs <- dfall[,"FieldID"]
  labIDs <- dfall[,"CAGRnumber"]
  Ex <- as.numeric(names(a[,1,1]))
  Em <- as.numeric(names(a[1,,1]))
  nlevels <- 50
  Peaks <- ex_ems
  peakCol <- "Peak"
  peakEx <- "ExCA"
  peakEm <- "EmCA"
  
  filenm <- file.path("cached_figures",paste0("EEMs",base.name,".pdf"))
  pdf(filenm)
  for (i in 1:dim(dfall)[1]){
    
    mat <- a[,,labIDs[i]]
    
    mainTitle <- paste("GLPF GLRI ",dfall[i,"USGSFieldID"],
                       dfall[i,"pdate"])
    
    EEMsPlot <- plotEEMs(mat=mat,Ex=Ex,Em=Em,
                         nlevels=nlevels,Peaks=Peaks,
                         peakCol=peakCol,
                         peakEx=peakEx,peakEm=peakEm,
                         mainTitle=mainTitle,titleSize = 1.0,
                         xlim = c(239,500), ylim = c(239,670))
  }
  dev.off()
  
  ##############################################################################################
  ##############################################################################################
  # Plot Abs
  
  fieldIDs <- dfall[,"USGSFieldID"]
  labIDs <- dfall[,"CAGRnumber"]
  
  filenm <- file.path("cached_figures",paste0("Abs",base.name,".pdf"))
  pdf(filenm)
  for (i in 1:dim(dfall)[1]){
    mat <- a[,,labIDs[i]]
    mainTitle <- paste("GLPF GLRI",dfall[i,"USGSFieldID"],dfall[i,"pdate"])
    
    plotAbs(dataAbs = dfabs,WaveCol = "nm",absCol = labIDs[i],main=mainTitle)
    
  }
  dev.off()
  
}

base.name <- "_noWW_noQA"
plotOpticalSummaries(base.name, cached.path)

base.name <- "_noQA"
plotOpticalSummaries(base.name, cached.path)

