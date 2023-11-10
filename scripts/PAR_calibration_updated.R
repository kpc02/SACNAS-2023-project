#install.packages("tidyverse")
#install.packages("readxl")
#install.packages("googlesheets4")
#install.packages("devtools")
#devtools::install_github("simonscmap/cmap4r/cmap4r")

library(tidyverse)
library(readxl)
library(cmap4r) 
library(googlesheets4)

# To set the API authorization key
set_authorization(cmap_key = "d24b5240-ea3a-45bd-8174-fdb7830d805f")


# To reset the authorization key
# set_authorization(reset = TRUE)

# To delete the key
# keyring::keyring_delete()




#--------------
# LOAD PAR data
#--------------
googledrive::drive_auth()
gs4_auth(token = googledrive::drive_token())
seaflow.meta <- read_sheet("https://docs.google.com/spreadsheets/d/1Tsi7OWIZWfCQJqLDpId2aG_i-8Cp-p63PYjjvDkOtH4")

#### this loads in the seaflow data from google sheets ####

read_sfl <- function(x){
  df <- read.csv(url(x), sep ="\t") #### this reads the url link u put in read_sfl and separates the string by a tab indentation ####

    #parse cruise name and serial number of instrument
    exp <- unlist(list(strsplit(sub(".sfl", "", basename(x)),"_"))) #### basename(x) gets rid of everything before cruise name and num of instrument, sub() gets rid of .sfl, strsplit creates strings by separating entire string using '_', bc of creation of new strings it made it into separate list and had to be listed and unlisted (??) ####
                                                                    #### 130 is version 2 of seaflow instrumnet ####
      if(length(exp) > 2) { cruise <- paste(exp[1],exp[2],sep="_")
      } else if(length(exp) ==2) cruise <- exp[1]
      print(cruise) #### if exp which is cruise name, number, and sometimes another num, is bigger than 2 then it will paste only the first 2 strings to avoid duplication in cruise. if its equal to 2 then just put first string in cruise (?) ####
      ship <- as.character(seaflow.meta[which(seaflow.meta$cruise == cruise),'Ship']) #### same thing as Boolean indexing in python (df[df[column == value]]); in seaflow meta df, finds cruise == cruise but its called ship in this case since cruise is for sfl ####
                                ##### cruise is TN397_130, so it finds cruisename within cruise and then finds matching ship ####
  df$cruise <- cruise
  df$ship <- ship
## $ refers to column like df.whatever
  return(df)
}

#### function called 'read_sfl' gets rid of specific cruise sfl link to leave behind cruise name and number ####  

# Merge all PAR data

## this is reading only one cruise and eventually i want to do it so it can read all the cruises
## sfl is a type of file that merges ship and seaflow data into one
## better to read sfl data from github than downloading all the cruise files onto computer

sfl <- read_sfl("https://raw.githubusercontent.com/seaflow-uw/seaflow-sfl/master/curated/TN397_130_130.sfl")

#### keep this below
#### when excluding these two lines of code, line 67 would not run because sfl was not defined so translate that
list.sfl <- list.files("~/Documents/Codes/seaflow-sfl/curated", pattern=".sfl", full.names=T) 
sfl <- do.call(rbind, lapply(list.sfl, function(x) read_sfl(x)))
#### list.sfl finds and lists any files with the format .sfl and makes it display the full name; since the files are from local directory (?), they need to pass through read_sfl function as well ####

## this code changes the date format from character to numeric datetime in UTC

sfl$DATE <- as.POSIXct(sfl$DATE, format = "%Y-%m-%dT%H:%M:%S", tz = "UTC")

## %>% is a more efficient way to write code and pipes through code to get last variable

sfl %>% group_by(cruise, DATE = cut(DATE, "1 hour")) %>%
        summarize(PAR = mean(PAR, na.rm=T)) %>% #### .agg is same as summarize ####
        ggplot(aes(as.POSIXct(DATE), PAR)) + 
        geom_line() +
        labs(x="time", y="Uncalibrated PAR") +
        theme_bw() + 
        facet_wrap(~ cruise, scales="free_x") 
#### plotting hourly par data ####

# PAR unit in microeinstein m^-2 s^-1, converted to einstein m^-2 day^-1 for compatibility with satellite data
# hes overriding the PAR column to convert to new unit
sfl$PAR <- 10^-6 * sfl$PAR * 24 * 60 * 60

# Calculate Daily mean
sfl_daily <- sfl %>% select(cruise, ship, DATE, LAT, LON, PAR) %>%
            group_by(cruise, ship, DATE = cut(DATE, breaks = "1 day")) %>%
            summarise(min_LAT = min(LAT, na.rm=T),
                      max_LAT = max(LAT, na.rm=T),
                      min_LON = min(LON, na.rm=T),
                      max_LON = max(LON, na.rm=T),
                      PAR = mean(PAR, na.rm=T)) %>%
            mutate(sat_PAR = 0) %>% ### adding sat par column and setting it to zero ####
            mutate(DATE = as.Date(DATE)) %>%
            arrange(DATE) ### sorts by date ###

#### using group by to get the daily mean of par based on cruise, ship, etc. and removes any NaN values found. gets minimum and max of lat and lon based on daily scale####


# clean up dataset
id <- which(!is.finite(sfl_daily$min_LAT))

if (length(id)!=0) {
  sfl_daily <- sfl_daily[-c(id),] #### if there is bad min lat values in id, then it will find everything but the bad values and put it in sfl_daily ####
}

#### tyring to find messed up min lat numbers (id is list of bad min lat) ####


#--------------------------
# Colocalize Modis_PAR data
#--------------------------

for(i in 1:nrow(sfl_daily)){

    print(sfl_daily$DATE[i])
    df <- get_spacetime(tableName = "tblModis_PAR",
                varName = "PAR",
                dt1 = as.character(sfl_daily$DATE[i]),
                dt2 = as.character(sfl_daily$DATE[i]),
                lat1 = sfl_daily$min_LAT[i],
                lat2 = sfl_daily$max_LAT[i],
                lon1 = sfl_daily$min_LON[i],
                lon2 = sfl_daily$max_LON[i],
                depth1 = 0,
                depth2 = 10)
    sfl_daily$sat_PAR[i] <- mean(df$PAR, na.rm=T)
}

#### getting cmap satellite par data and putting it into sat_par col in sfl_daily (putting avg satellite of specified location into daily ship data) ####

#----------------------------
# Calculate correction factor
#----------------------------

## troubleshoot and figure out which is what type of PAR data and mess around with plot to get the labels for the right PAR data
sfl_daily <- sfl_daily %>% group_by(cruise) %>%
              mutate(correction = median(sat_PAR / PAR, na.rm=T))

#### this creates a new column in sfl_daily (based on mutate function ?) and contains corrected median par values, that is grouped by the cruise (only 1 right now) ####

# Plot in-situ vs sat PAR
sfl_daily %>% ggplot(aes(DATE, PAR)) +
                geom_point(aes(color = "red"),size=2, alpha = 0.5, col = 2)  + 
                geom_point(aes(DATE, sat_PAR, colour = "black"), alpha = 0.5, col =1) +
                geom_point(aes(DATE, PAR * correction,colour = "green") , col = 3, alpha = 0.5) +
                labs(x="time", y="Uncalibrated PAR") +
                theme_bw() + 
                facet_wrap(~ cruise, scales="free_x")

#### creates a plot of date and par for different types of par data ####

##first geompoint is original PAR, second one is satellite and ship PAR (maybe), 
  #and third one is satellite and ship PAR calibrated (BIG MAYBE)


# Calculate correction factor for each cruise
calib <- sfl_daily %>% group_by(cruise, ship) %>%
                       summarize(time = mean(DATE),
                                correction = round(mean(correction),3)) %>%
                       arrange(time)
#### idk ####

## this is ... idk i forgot but she said not to worry about this so im commenting it out
#write_csv(calib, "~/Documents/Projects/CALIBRATION/PAR_calibration/par-calibration.csv")
#write_csv(calib, "~/Documents/Codes/popcycle/inst/par/par-calibration.csv")
