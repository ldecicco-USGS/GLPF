library(googlesheets)
library(dplyr)
library(lubridate)

cached.path <- "cached_data"
cached.save <- "0_download"
# Make sure Bison Connect is your default, and log out for good measure any personal gmail

getGoogleData <- function(cached.path, cached.save){
  
  
  my_sheets <- gs_ls()
  
  glpfTitle <- gs_title("GLPF sample tracking form.xlsx")
  # 
  GMToffset <- c(NA,5,5,6)
  
  for (i in 2:4){
    dfState <- gs_read(glpfTitle,ws=i,range=cell_cols("A:Z"))
    names(dfState) <- make.names(names(dfState))
    dfState <- filter(dfState,!is.na(Start.date.time..mm.dd.yy.hh.mm.))
    dfState$pdate <- as.POSIXct(dfState$Start.date.time..mm.dd.yy.hh.mm.,format='%m/%d/%Y %H:%M:%S',tz="GMT")
    dfState$pdate <- dfState$pdate + GMToffset[i]*60*60
    
    dfState$date <- as.Date(dfState$pdate)

    dfState$UWMFT <- suppressWarnings(as.integer(dfState$UWMFT))
    dfState$MIBARLID <- suppressWarnings(as.integer(dfState$MIBARLID))

    dfState$entero.cn.100g <- suppressWarnings(as.character(dfState$entero.cn.100g))
    dfState$bachum.cn.100g <- suppressWarnings(as.character(dfState$bachum.cn.100g))
    dfState$lachno.cn.100g <- suppressWarnings(as.character(dfState$lachno.cn.100g))
    dfState$lachno3.cn.100ml <- suppressWarnings(as.character(dfState$lachno3.cn.100ml))
    dfState$ecoli.cn.100g <- suppressWarnings(as.character(dfState$ecoli.cn.100g))
    dfState$dog.cn.100ml <- suppressWarnings(as.character(dfState$dog.cn.100ml))
    dfState$ecoli.cn.100g <- suppressWarnings(as.character(dfState$ecoli.cn.100g))
    dfState$racoon.cn.100ml..relative.. <- suppressWarnings(as.character(dfState$racoon.cn.100ml..relative..))
    
    dfState <- rename(dfState, 
                      racoon=racoon.cn.100ml..relative..,
                      dog = dog.cn.100ml)
    dfState$State <- gsub(' ','',dfState$State)
    
    names(dfState)[2] <- "USGS_SiteID"
    
    if(i == 2) {
      df <- dfState
    } else {
      df <- bind_rows(df,dfState)
    }
  }
  
  names(df) <- gsub('\\.','',names(df))
  names(df) <- gsub('\\?','',names(df))
  
  df <- df[!(is.na(df$CAGRnumber) | df$CAGRnumber %in% c("x")),]
  
  saveRDS(df, file.path(cached.path,cached.save,"tracking.rds"))

}

getGoogleData(cached.path, cached.save)

getGLRIsites <- function(cached.path, cached.save){
  
  glriTitle <- gs_title("GLRI Virus sample tracking form_V2_20150520.xlsx")
  
  sheets <- c("NY BacOptic","WI BacOptic")
  
  df <- data.frame()
  
  for(i in sheets){
    
    glriTracking <- gs_read(glriTitle,ws=i, range=cell_cols("A:Z"))
    names(glriTracking) <- make.names(names(glriTracking))
    
    glriTracking$pdate <- parse_date_time(glriTracking$Start.date.time..mm.dd.yy.hh.mm., orders = c("m/d/y H:M:S","m/d/y H:M"))
    
    glriTracking$pdate <- glriTracking$pdate + 5*60*60
    
    glriTracking$date <- as.Date(glriTracking$pdate)
    
    glriTracking$UWMFT <- suppressWarnings(as.integer(glriTracking$UWMFT))
    glriTracking$MIBARLID <- suppressWarnings(as.integer(glriTracking$MIBARLID))
  
    glriTracking$State <- gsub(' ','',glriTracking$State)
    
    names(glriTracking)[2] <- "USGS_SiteID"
    
    names(glriTracking) <- gsub('\\.','',names(glriTracking))
    names(glriTracking) <- gsub('\\?','',names(glriTracking))
    
    glriTracking <- glriTracking[!(is.na(glriTracking$CAGRnumber) | glriTracking$CAGRnumber %in% c("x")),]
  
    if(i == sheets[1]){
      df <- glriTracking
    } else {
      df <- dplyr::bind_rows(df, glriTracking)
    }
    
  }

  saveRDS(df, file.path(cached.path,cached.save,"tracking_GLRI.rds"))
  
  
}

getGLRIsites(cached.path, cached.save)
