library(dplyr)
library(ggplot2)
library(tidyr)

cached.path <- "cached_data"
base.name <- "_noWW_noQA"

summaryDF <- readRDS(file.path(cached.path,"8_process_new_categories","rds",paste0("summary",base.name,".rds")))

EEMs <- readRDS(file.path(cached.path,"7_process_summarize_optics","rds",paste0("EEMs3D",base.name,".rds")))

grNum <- summaryDF$CAGRnumber[1]

EEM <- as.data.frame(as.table(EEMs[,,grNum])) %>%
  mutate(Ex = as.numeric(as.character(Var1)),
         Em = as.numeric(as.character(Var2))) %>%
  select(-Var1, -Var2, Intensity = Freq) %>%
  filter(!is.na(Intensity)) 

intensity.limits = c(0,5)  

plot_variables_onEEMs <- function(variables, EEM, 
                                  intensity.limits = c(0,5)){
  
  text_list <- get_var_position(variables)
  
  colPal <- c(colors()[30], colors()[143],colors()[554])
  
  heatPage <- ggplot() +
    geom_tile(data = EEM, aes(x=Ex, y=Em, fill = Intensity)) +
    theme_bw() + 
    scale_fill_gradientn(colours=colPal,
                         na.value = "transparent",
                         limits= intensity.limits) +
    xlim(c(240,500)) +
    ylim(c(210, 600)) 
  
  if(nrow(text_list$point_text) > 0){
    heatPage <- heatPage +
      geom_text(data = text_list$point_text, 
                aes(x = ExCA, y = EmCA, label = Peak))
  }
  
  if(nrow(text_list$line_text) > 0){
    text_positions <- text_list$line_text %>%
      rowwise() %>%
      mutate(Ex_mean = mean(c(Ex1, Ex2), na.rm = TRUE),
             Em_mean = mean(c(Em1, Em2), na.rm = TRUE)) %>%
      data.frame()
    
    line_positions <- text_list$line_text %>%
      select(-Ex2) %>%
      gather(x, y, -Peak, -Ex1)
    
    heatPage <- heatPage +
      geom_line(data = line_positions,
                aes(x = Ex1, y = y, group = Peak), 
                color = "grey80", alpha = 0.3, size = 1.5) +
      geom_text(data = text_positions, 
                aes(x = Ex_mean, y = Em_mean, label = Peak)) 
  }
  
  if(nrow(text_list$box_text) > 0){
    text_positions <- text_list$box_text %>%
      rowwise() %>%
      mutate(Ex_mean = mean(c(Ex1, Ex2), na.rm = TRUE),
             Em_mean = mean(c(Em1, Em2), na.rm = TRUE)) %>%
      data.frame()

    heatPage <- heatPage +
      geom_rect(data = text_positions,
                aes(xmin = Ex1, xmax = Ex2,
                    ymin = Em1, ymax = Em2), 
                fill = "grey80", alpha = 0.3) +
      geom_text(data = text_positions, 
                aes(x = Ex_mean, y = Em_mean, label = Peak)) 
  } 
  
  print(heatPage)
  return(heatPage)
}

get_var_position <- function(variables){
  
  variables <- gsub("log","", variables)
  
  rat <- which(substring(variables, 1, 1) == "r")
  rB <- which(substring(variables, 1, 2) == "rB")
  rat <- rat[!(rat %in% rB)]
  
  rat_sig <- strsplit(gsub("r","",variables[rat]),"_")
  rat_sig <- t(sapply(rat_sig, function(x) x))
  

  rB_sig <- strsplit(gsub("rB","",variables[rB]),"_")
  rB_sig <- t(sapply(rB_sig, function(x) x))[,1:2]
  
  variables <- c(variables[-c(rat,rB)], c(rat_sig), c(rB_sig))
  variables <- unique(variables)
  
  Peaks <- USGSHydroOpt::ex_ems
  other_sigs <- read.csv(file.path("raw_data","opticalSummary","ex_ems_meansCA.csv"),stringsAsFactors = FALSE)
  other_sigs <- other_sigs[!(other_sigs$Peak %in% Peaks$Peak),]
  
  other_sig_points <- other_sigs[is.na(other_sigs$Ex2) & is.na(other_sigs$Em2),]
  other_sig_box <- other_sigs[!is.na(other_sigs$Ex2) & !is.na(other_sigs$Em2),]
  other_sig_lines <- other_sigs[is.na(other_sigs$Ex2) & !is.na(other_sigs$Em2),]
  
  point_text <- data.frame()
  line_text <- data.frame()
  ratio_text <- data.frame()
  box_text <- data.frame()
  
  if(any(variables %in% Peaks$Peak)){
    point_text <- Peaks[which(Peaks$Peak %in% 
                                variables[(variables %in% Peaks$Peak)]),]
  }
  
  if(any(variables %in% other_sig_points$Peak)){
    point_text_sub <- other_sig_points %>%
      filter(Peak %in% variables) %>%
      select(Peak, ExCA = Ex1, EmCA = Em1)
    
    point_text <- bind_rows(point_text, point_text_sub)
  }
  
  if(any(variables %in% other_sig_lines$Peak)){
    line_text <- other_sig_lines %>%
      filter(Peak %in% variables) %>%
      select(-Source)
  }

  if(length(grep("A\\d{3}", variables)) > 0){
    abs <- variables[grep("A\\d{3}", variables)]
    
    point_text_sub <- data.frame(Peak = abs,
                                 ExCA = as.numeric(gsub("A","",variables[grep("A\\d{3}", variables)])),
                                 EmCA = 250,
                                 stringsAsFactors = FALSE)
    
    point_text <- bind_rows(point_text, point_text_sub)
    
  }
  
  if(length(grep("Sag\\d{3}_\\d{3}", variables)) > 0){
    sags <- variables[grep("Sag\\d{3}_\\d{3}", variables)]
    
    waves <- strsplit(gsub("Sag","",variables[grep("Sag\\d{3}_\\d{3}", variables)]),"_")
    waves <- sapply(waves, function(x) x)
    
    wave_1 <- as.numeric(waves[1,])
    wave_2 <- as.numeric(waves[2,])
    mean_wave <- (wave_1 + wave_2)/2
    
    line_sub <- data.frame(
      Peak = sags, 
      Ex1 = wave_1, 
      Ex2 = wave_2,
      Em1 = 225,
      Em2 = 225,
      stringsAsFactors = FALSE)
    
    line_text <- bind_rows(line_text, line_sub)
    
  }
    
  if(any(variables %in% other_sig_box$Peak)){
    box_text <- other_sig_box %>%
      filter(Peak %in% variables) %>%
      select(-Source)
  }
  
  return(list(point_text=point_text,
              line_text=line_text,
              box_text = box_text,
              ratio_text = ratio_text))
}


variables <- c("Sag240_255",
               "rSn.5_Sn.4",
               "rSn.9_S1.50",
               "rBS44_S45_BF",
               "Mrange.25",
               "Sag290_350",
               "logH1",
               "A290",
               "W",
               "S1.50",
               "Sn.2",
               "T")

load("C:/Users/ldecicco/Downloads/modelCoefList.RData")

pdf("steves_variables.pdf", width = 11, height = 7.6)
for(i in modelCoefList){
  cat(i,"\n")
  browser()
  plot_variables_onEEMs(i, EEM )
}
dev.off()

pdf("test.pdf", width = 11, height = 7.6)
plot_variables_onEEMs(variables, EEM )
dev.off()



