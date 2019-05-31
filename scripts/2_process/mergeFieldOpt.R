library(dataRetrieval)
library(lubridate)
library(dplyr)

raw.path <- "raw_data"
cached.path <- "cached_data"
cached.save <- "4_process_merge_field_opt"

mergeFieldOpt <- function(raw.path, cached.path, cached.save){

  data.opt <- readRDS(file.path(cached.path,"3_process_merge_bacteria","trackingBacteria.rds"))
  
  data.opt <- mutate(data.opt,
                      USGSSTAID = zeroPad(ifelse(is.na(USGSSTAID), "", as.character(USGSSTAID)), padTo = 8),
                      USGSNWISStationIDifapplicable = zeroPad(ifelse(is.na(USGS_SiteID), "", as.character(USGS_SiteID)), padTo = 8),
                      startDateTime = pdate)#parse_date_time(Startdatetimemmddyyhhmm, c("%m/%d/%Y %H:%M:%S", "%m/%d/%Y %H:%M")),
                      
  data.opt <- mutate(data.opt,
                     endDateTime = parse_date_time(Enddatetimemmddyyhhmm, c("%m/%d/%Y %H:%M:%S", "%m/%d/%Y %H:%M")))

  ###############################################################
  # NY:
  ###############################################################
  data.ny <- readRDS(file.path(cached.path,"1_munge","dataNY.rds"))
  data.opt.ny <- filter(data.opt, State == "NY")
  
  data.ny$FieldID[data.ny$FieldID == "WWRCTQ-1" & data.ny$WT == 18] <- "WWRCTQ-002"
  data.ny$FieldID[data.ny$FieldID == "WWRCTE-1or2"] <- "WWRCTE-1"
  data.ny$FieldID[data.ny$FieldID == "WWRCTJ-1or2" & data.ny$SC == 1294] <- "WWRCTJ-1"
  data.ny$FieldID[data.ny$FieldID == "WWRCTJ-1or2"] <-"WWRCTJ-003"
  data.ny$FieldID[data.ny$FieldID == "WWRCTG-1" & data.ny$WT == 16.7] <- "WWRCTG-002"
  data.ny$FieldID[data.ny$FieldID == "WWRCTH-1" & data.ny$SC == 1368] <- "WWRCTH-002"
  data.ny$FieldID[data.ny$FieldID == "WWRCTF-1" & data.ny$SC == 1849] <- "WWRCTF-002"

  data.opt.ny$FieldID <- as.character(sapply(data.opt.ny$FieldID, function(x) {
    paste(strsplit(x,split = "-")[[1]][1],
          strsplit(x,split = "-")[[1]][2], sep="-")
    }))

  # Dates don't seem to be close
  data.ny.merge <- right_join(data.ny, data.opt.ny, 
                              by=c("SAMPLE_START_DT"="startDateTime",
                                   "FieldID"))
  data.ny.merge <- data.ny.merge[!is.na(data.ny.merge$CAGRnumber),]

  ###############################################################
  # MI:
  ###############################################################
  data.mi <- readRDS(file.path(cached.path,"1_munge","dataMI.rds"))
  data.opt.mi <- filter(data.opt, State == "MI")
  
  data.opt.mi$roundDate <- as.POSIXct(round(data.opt.mi$startDateTime,"min"))
  # data.opt.mi$FieldID_sm <- gsub("WW","", data.opt.mi$FieldID)
  # data.opt.mi$FieldID_sm <- sapply(strsplit(data.opt.mi$FieldID_sm, "-"), function(x) x[[1]][1])
  # 
  data.mi.merge <- right_join(select(data.mi,-FieldID), data.opt.mi, 
                             by=c("SAMPLE_START_DT"="roundDate",
                                  "SiteID" = "USGS_SiteID")) %>%
    filter(!is.na(CAGRnumber)) %>%
    mutate(startDateTime=SAMPLE_START_DT,
           UVch1 = NA,
           UVch2 = NA,
           UVch3 = NA) %>%
    rename(USGS_SiteID = SiteID)
  
  data.mi.merge <- data.mi.merge[, names(data.ny.merge)]
  
  ###############################################################
  # WI:
  ###############################################################
  data.wi <- readRDS(file.path(cached.path,"1_munge","dataWI.rds"))
  data.opt.wi <- filter(data.opt, State == "WI")

  data.wi.merge <- right_join(data.wi, data.opt.wi, by=c("FieldID"="FieldID")) %>%
    filter(!is.na(CAGRnumber))
  data.wi.merge <- data.wi.merge[!duplicated(data.wi.merge$CAGRnumber),]
  
  data.wi.merge <- data.wi.merge[, names(data.ny.merge)]
  
  ###############################################################
  # Full merge:
  ###############################################################
  
  data.merge <- bind_rows(data.mi.merge,data.wi.merge,data.ny.merge)
  
  cuts <- c(0,3,6,9)
  labs <- c("Winter", "Spring", "Summer", "Fall")
  months <- as.POSIXlt(data.merge$pdate)$mon
  months <- months + 1
  months[months == 12] <- 0
  
  data.merge$Season <- labs[findInterval(months,cuts)]
  names(data.merge) <- enc2utf8(names(data.merge))

  data.merge.check <- select(data.merge, CAGRnumber, State, FieldID, 
                             USGS_SiteID, VolumeFilteredL, Date,
                             FilterA04µMUSGSMIBARL,FilterB02µMUWMSFS,
                             Site,USGSNWISStationIDifapplicable,USGSSTAID,
                             SAMPLE_START_DT, endDateTime,eventNumInd,
                             Startdatetimemmddyyhhmm, Enddatetimemmddyyhhmm,
                             pdate, date, pedate,Time,Season,
                             MIBARLID,UWMFT,hydroCondition,
                             Comments, Comments2,VirusAutosampleorSewerGrab,IfTributarySewerOutfallManholeorDitch,
                             project,Matrix,SampleType,eventNum,
                             SampleType9regular2blank7replicate,
                             FilterAVolumemL,FilterBVolumemL,FilterAUSGSMIBARLVolumemL,FilterBUWMVolumemL,TotalAutoSamplerVolumeL)
  
  dir.create(file.path(cached.path,cached.save), showWarnings = FALSE)
  
  saveRDS(data.merge.check, file.path(cached.path,cached.save,"mergeCheck.rds"))
  
  merged.data <- data.merge[,names(data.merge)[!(names(data.merge) %in% names(data.merge.check))]]
  
  merged.data <- cbind(data.merge[,c("CAGRnumber","pdate","endDateTime","FieldID",
                                     "MIBARLID","UWMFT","eventNum",
                                     "Comments","hydroCondition","State","Season",
                                     "SampleType9regular2blank7replicate",
                                     "USGS_SiteID",
                                     "IfTributarySewerOutfallManholeorDitch",
                                     "VirusAutosampleorSewerGrab","project")], 
                       merged.data) %>%
    rename(SiteID=USGS_SiteID,
           pedate = endDateTime)
  
  dir.create(file.path(cached.path,cached.save), showWarnings = FALSE)
  
  saveRDS(merged.data, file.path(cached.path,cached.save,"mergedData.rds"))
  
  write.csv(merged.data, file.path(cached.path,cached.save,"mergedData.csv"), row.names = FALSE)
  write.csv(data.merge.check, file.path(cached.path,cached.save,"mergedDataSupplement.csv"), row.names = FALSE)
  
}

mergeFieldOpt(raw.optics, cached.path, cached.save)
