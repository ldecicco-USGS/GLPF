library(dplyr)
library(data.table)

#Raw data folder:
raw.path <- "raw_data"
cached.path <- "cached_data"
cached.save <- "3_process_merge_bacteria"

mergeTrackingBact <- function(raw.path, cached.path, cached.save){
  df <- readRDS(file.path(cached.path,"0_download","tracking.rds" ))

  dfglpfSummary <- setDF(fread(file.path(raw.path,"optics","GLPF_Summary.csv")))
  load(file=file.path(raw.path,"PreCleaned",'GLRI01-04-16_mergedBact.RData'))
  load(file=file.path(raw.path,"PreCleaned",'GLRIWWMar162016summary.RData'))

  names(dfall)[names(dfall) == "USGS.NWIS.Station.ID..if.applicable."] <- "USGS_SiteID"
  dfall$USGS_SiteID <- zeroPad(dfall$USGS_SiteID, 8)
  dfall <- left_join(dfall, select(dfOpt, CAGRnumber, DOCResult, TDNResult), by="CAGRnumber")
  
  # "All GLRI samples were event samples, so just use that....Steve 8/26/16"
  # dfall == GLRI
  dfall$Comments[dfall$Comments == ""] <- "Event Grab"
  dfall$State <- gsub(' ','',dfall$State)
  dfall$pdate <- as.POSIXct(dfall$Start.date.time..mm.dd.yy.hh.mm.,format='%m/%d/%Y %H:%M',tz="GMT")

  GMTOffset <- ifelse(dfall$State=='WI',6,5)
  dfall$pdate <- dfall$pdate + GMTOffset*60*60
  names(dfall) <- gsub("\\.","",names(dfall))
  names(dfall) <- gsub('\\?','',names(dfall))
  
  dfall <- rename(dfall,
                  FilterA04µMUSGSMIBARL = FilterA04ÂµMUSGSMIBARL,
                  FilterB02µMUWMSFS = FilterB02ÂµMUWMSFS, 
                  enterocn100g = Entero,
                  Esp2cn100ml = espCN100mL,
                  lachnocn100g = Lachno,
                  bachumcn100g = BacHuman) %>%
    mutate(bachumcn100g = as.character(bachumcn100g),
           lachnocn100g = as.character(lachnocn100g),
           enterocn100g = as.character(enterocn100g),
           project = 'GLRI')
  
  df$MIBARLID <- as.character(df$MIBARLID)
  df$project <- 'GLPF'
  
  dfglpfDOC <- dfglpfSummary[,c('Grnumber','DOCResult','TDNResult')]
  dfglpfDOC <- dfglpfDOC[which(dfglpfDOC$Grnumber %in% df$CAGRnumber),]
  
  df <- left_join(df,
                  dfglpfDOC,
                  by = c('CAGRnumber'='Grnumber'))
  
  df <- full_join(df, dfall)

  df <- select(df, 
               -qpcrBacHumanResults,
               -qpcrEnterococcusResults,
               -qpcrLachno2Results) %>%
    filter(CAGRnumber != "")
  # Deal with zeros and BLDs, convert to numeric
  
  BacHumCensLevel <- 225
  LachnoCensLevel <- 225
  dogCensLevel <- 225
  racoonCensLevel <- 225
  entCensLevel <- 225
  eColiCensLevel <- 225
  espCensLevel <- 80
  ipaCensLevel <- 200

  noDetectIndicators <- c('BLD','BLQ','0','<LRL',"ND")

  splitCens <- function(df, base.name, column.name, noDetectIndicators, censor.level){
    df[[paste0(base.name,"_rk")]] <- ifelse(df[[column.name]] %in% noDetectIndicators, "<", "")
    df[[base.name]] <- suppressWarnings(ifelse(df[[paste0(base.name,"_rk")]] == "<", censor.level, as.numeric(df[[column.name]])))
    df <- select_(df,.dots=paste("-",column.name))
    return(df)
  }

  df <- splitCens(df, "bacHum", "bachumcn100g", noDetectIndicators, BacHumCensLevel)
  df <- splitCens(df, "lachno", "lachnocn100g", noDetectIndicators, LachnoCensLevel)
  df <- splitCens(df, "lachno3", "lachno3cn100ml", noDetectIndicators, LachnoCensLevel)
  df <- splitCens(df, "ent", "enterocn100g", noDetectIndicators, entCensLevel)
  df <- splitCens(df, "eColi", "ecolicn100g", noDetectIndicators, eColiCensLevel)
  df <- splitCens(df, "esp", "Esp2cn100ml", noDetectIndicators, espCensLevel)
  df <- splitCens(df, "ipaH", "ipaHcn100ml", noDetectIndicators, ipaCensLevel)
  df <- splitCens(df, "dog_marker", "dog", noDetectIndicators, dogCensLevel)
  df <- splitCens(df, "racoon_marker", "racoon", noDetectIndicators, racoonCensLevel)
  
  df$Site<-substr(df$FilterB02µMUWMSFS,1,nchar(df$FilterB02µMUWMSFS)-3)
  
  df$FieldID <- df$FilterB02µMUWMSFS
  
  df$hydroCondition <- "Base"
  df$hydroCondition[grep("event",df$Comments, ignore.case = TRUE)] <- "Event"

  df$eventNum <- 1

  df$eventNumInd <- 1
  
  eventNum <- ""
  for(state in c("MI","WI","NY")){
    subDF <- df[df$State == state,]
    subDF$Julian <- as.numeric(julian(subDF$pdate,origin=as.Date("1850-01-01")))
    for(i in 2:nrow(subDF)){
      behind <- subDF$Julian[i] - subDF$Julian[i-1]
      if(abs(behind) < 1){
        subDF$eventNumInd[i] <- subDF$eventNumInd[i-1] 
      } else {
        subDF$eventNumInd[i] <- subDF$eventNumInd[i-1] + 1
      }
    }
    df$eventNum[df$State == state] <- paste0(state,subDF$eventNumInd)
  }

  df <- df[!(df$CAGRnumber %in% c("x","N/A")),]
  
  dir.create(file.path(cached.path,cached.save), showWarnings = FALSE)
  
  saveRDS(df,file=file.path(cached.path,cached.save,'trackingBacteria.rds'))
  
}
  
mergeTrackingBact(raw.path, cached.path, cached.save)