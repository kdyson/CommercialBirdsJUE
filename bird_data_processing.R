## This script takes care of bird data manipulation and processing for
## 'bird_data_analysis.R'.

## Data manipulation (this doc) --> analysis --> Rmarkdown/bookdown

## Written 10/02/2017. This code is based on splitting out the data manipulation
## and processing contained in previous files including 'bird_data_analysisX.R'
## and 'Office_Development_Bird_CommunitiesX.Rmd'.



## -- SETUP ----------------------------------------------------------------------

## Source helper R codes
source("../R_Scripts/triplet_fixer.R")
source("../Vegetation/vegetation_data_processing.R")
remove(dens.matrify.tree.sp.size, dens.tree.biomass, matrify.large.trees,
       matrify.shrub.zones, matrify.tree.sp.size #,
       #matrify.shrub.site, matrify.tree.sp.only, matrify.gc
)


setwd("../Birds")
options(scipen = 999)
options(tibble.width = Inf)


    ## Load necessary libraries
        library(tidyr)
        library(stringr)
        library(vegan)
        library(chron)
        library(dplyr)

    ## Set variables
        write.all <- FALSE



        ## -- ENVIRONMENTAL VARIABLES ----------------------------------------------------------------


        Dec_14 <-
            read.csv(file = "../../DataRepository/WeatherData/Hourly_WSU_W21/AWN_Dec_14.csv",
                     stringsAsFactors = F)
        Dec_15 <-
            read.csv(file = "../../DataRepository/WeatherData/Hourly_WSU_W21/AWN_Dec_15.csv",
                     stringsAsFactors = F)
        Jan_15 <-
            read.csv(file = "../../DataRepository/WeatherData/Hourly_WSU_W21/AWN_Jan_15.csv",
                     stringsAsFactors = F)
        Jan_16 <-
            read.csv(file = "../../DataRepository/WeatherData/Hourly_WSU_W21/AWN_Jan_16.csv",
                     stringsAsFactors = F)
        Feb_15 <-
            read.csv(file = "../../DataRepository/WeatherData/Hourly_WSU_W21/AWN_Feb_15.csv",
                     stringsAsFactors = F)
        Feb_16 <-
            read.csv(file = "../../DataRepository/WeatherData/Hourly_WSU_W21/AWN_Feb_16.csv",
                     stringsAsFactors = F)

        ## Create and join a weather data set based on hourly data from WSU_W21

        # Join all data sets
        Weather_Hour <-
            rbind.data.frame(Dec_14, Jan_15, Feb_15, Dec_15, Jan_16, Feb_16)
        # Get rid of empty column
        Weather_Hour <- Weather_Hour[-16]
        # Check n rows is right
        nrow(Weather_Hour) == nrow(Feb_16) + nrow(Feb_15) + nrow(Jan_16) + nrow(Jan_15) +
            nrow(Dec_15) + nrow(Dec_14)
        # Rename columns
        colnames(Weather_Hour) <-
            c(
                'Date',
                'Hour.PST',
                'Min.Air.Temp.F',
                'Avg.Air.Temp.F',
                'Max.Air.Temp.F',
                'Avg.Air.Temp.F.2nd',
                'Dew.Point.F',
                'Rel.Hum.pct',
                'Leaf.Wet.u',
                'Wind.Dir.and.Deg',
                'Speed.MPH',
                'Gust.MPH',
                'Soil.Temp.F.8in',
                'Tot.Precip.in',
                'Solar.Rad.W_m2'
            )



        # Properly format
        Weather_Hour$Date <- as.Date(Weather_Hour$Date, "%B %d, %Y")
        Weather_Hour$Hour.PST <-
            chron(times. = ifelse(
                Weather_Hour$Hour.PST == "24:00",
                as.character("00:00:00"),
                # This needs to be 00: not 24:
                paste(Weather_Hour$Hour.PST, ":00", sep = "")
            ))
        summary(Weather_Hour$Hour.PST) ## without ifelse, 24:00:00 is differently formatted than everything else...
        Weather_Hour$hours <- hours(Weather_Hour$Hour.PST)
        str(Weather_Hour) # should have 4362 observations; note raw downloaded data may have double data entries.


    # Create the native conifer density variable


        management.landscaping$tree.conifer.nat.dens <-
            dens.matrify.tree.sp.only$western.red.cedar +
            dens.matrify.tree.sp.only$western.hemlock +
            dens.matrify.tree.sp.only$douglas.fir
        management.landscaping$impervious.sqft <-
            dens.matrify.gc$impervious.sqft


        # Create the native shrub density variable
        native.shrubs <-
            shrub.native[which(shrub.native$shrub.origin == "native"), "shrub.species"]
        # native.shrubs.ns <- native.shrubs[-which(native.shrubs == "gaultheria.shallon")]

        management.landscaping$shrub.nat.dens <-
            rowSums(dens.matrify.shrub.site[, native.shrubs])
        management.landscaping$shrub.esr <-
            exp(diversity(dens.matrify.shrub.site[, native.shrubs]))


## -- BIRD OBSERVATIONS ------------------------------------------------------------------

    ## Pull in raw bird data and turn it into a useable table. need to calculate
    ## unknowns as portion of data: likely need to remove them because no way of
    ## telling if they are the same as other unknowns as there is with shrubs and
    ## trees.


        bird.data <- read.csv("../../DataRepository/BirdData/BirdData_all_07192016.csv",
                          stringsAsFactors = FALSE, header = FALSE)


    ## Fix column names (macro doesn't output column names; note colnames should be V1 etc due to header = FALSE)
        colnames(bird.data)
        colnames(bird.data) <-
            c(
                "site.name",
                "visit.number",
                "date",
                "cloud",
                "field.wind",
                "field.temp",
                "NWS.temp",
                "noise.level",
                "start",
                "end",
                "sweep.number",
                "on.site",
                "text.spp",
                "BBL",
                "A.vegetation.gleaning",
                "A.ground.foraging",
                "A.eating.berries",
                "A.eating.seeds",
                "A.sitting",
                "A.calling",
                "A.flitting",
                "detected.seen",
                "detected.call",
                "detected.other.noise",
                "detected.recorded.audio",
                "notes"
            )

        bird.data$date <-
            as.Date(bird.data$date, format = "%m/%d/%Y")
        bird.data$cloud <- as.factor(bird.data$cloud)
        bird.data$field.wind <- as.factor(bird.data$field.wind)
        bird.data$start <- chron(times. = bird.data$start)
        bird.data$end <- chron(times. = bird.data$end)
        # Note that there are some observations that did not occur during sweeps
        # and therefore do not have time data; these will show as na (170)



    ## If a bird was detected because it called, it should be recorded as
    ## activity:calling. Not sure why it wasn't consistently, other than saving
    ## time in the field by only checking one...
        bird.data$A.calling <-
            ifelse(
                bird.data$A.calling == "0" & bird.data$detected.call == "1",
                bird.data$detected.call,
                bird.data$A.calling
            )
        bird.data$start.hour <- hours(bird.data$start)


    # workaround in case time format is hh:mm instead of hh:mm:ss
    # bird.data$start <- chron(times. = paste(bird.data$start, ":00", sep=""))




    # Merge weather info with bird data.
        bird.data <- merge.data.frame(
            x = bird.data,
            y = Weather_Hour,
            by.x = c('date', 'start.hour'),
            by.y = c('Date', 'hours'),
            all.x = TRUE
        )

    ## There should be 3343 observations of birds on site; less 209 birds seen
    ## outside of official sweeps. Some of these two groups overlap so should have
    ## 3237 observations


        bird.data.onsite <-
            bird.data[which(bird.data$on.site == "yes" &
                                bird.data$sweep.number != "na"), ]
        summary(bird.data.onsite)
        colnames(bird.data.onsite)



    ## Call out where birds were observed foraging--Note that they must have been SEEN.


        bird.data.onsite$did.forage <-
            as.numeric(apply(bird.data.onsite[, 16:19], MARGIN = 1, FUN = max))

    ## Note, warning will pop up saying NAs introduced by coercion; this is
    ## because the no.birds have X's instead of numbers




## Incidence matrix

    ## Now, create the incidence matrix. This is for both years (or each year
    ## by itself). Collapses all data to give proportion of visits that species
    ## was observed per site.

        ## Calculate incidence
        incidence <-
            group_by(bird.data.onsite, site.name, text.spp) %>%
            summarise(
                visits.present = paste(as.character(unique(visit.number)), collapse = ", "),
                visits.2014 = sum(unique(visit.number) < 5),
                visits.2015 = sum(unique(visit.number) > 4),
                incidence = length(as.character(unique(visit.number))) /
                    8,
                incidence.2014 = visits.2014 / 4,
                incidence.2015 = visits.2015 / 4,
                ## Note: these just tell us if the bird was ever observed on site doing these actions.
                A.vegetation.gleaning = max(A.vegetation.gleaning),
                A.ground.foraging = max(A.ground.foraging),
                A.eating.berries = max(A.eating.berries),
                A.eating.seeds = max(A.eating.seeds),
                A.sitting = max(A.sitting),
                A.calling = max(A.calling),
                A.flitting = max(A.flitting),
                detected.seen = max(detected.seen),
                detected.call  = max(detected.call),
                detected.other.noise = max(detected.other.noise),
                detected.recorded.audio = max(detected.recorded.audio)
            )

        ## NOTE: as character must go outside unique, not inside!


        summary(incidence)      ## Mean incidence is slightly lower in 2015 compared with 2014.


## Presence/Absence
    ## Create Presence/Absence matrix for birds! Approach: Create two tables
    ## using summarise: one at the bird level and one at the site visit level.
    ## Then merge the two for a table describing presence absence of bird
    ## species for each site visit, combined with information about that visit.


        ## Create the bird level table
        bird.by.visit.PA <-
            group_by(bird.data.onsite, site.name, visit.number, text.spp) %>%
            summarise(
                BBL = max(BBL),
                sweeps.present = n(),
                sweeps.foraging = sum(did.forage),
                A.vegetation.gleaning = max(A.vegetation.gleaning),
                A.ground.foraging = max(A.ground.foraging),
                A.eating.berries = max(A.eating.berries),
                A.eating.seeds = max(A.eating.seeds),
                A.sitting = max(A.sitting),
                A.calling = max(A.calling),
                A.flitting = max(A.flitting),
                detected.seen = max(detected.seen),
                detected.call  = max(detected.call),
                detected.other.noise = max(detected.other.noise),
                detected.recorded.audio = max(detected.recorded.audio),
                notes = paste(unique(notes, na.rm = TRUE), collapse = ", ")
            )
        ## Note: here, vegetation gleaning etc. tell us if the bird was
        ## observed doing that activity on that visit!


## Need to add a few variables to describe aggregate foraging behavior, etc.


        ## When I know the bird was foraging (recorded)
        bird.by.visit.PA$B.foraging <-
            apply(X = bird.by.visit.PA[, 6:9], MARGIN = 1, FUN = max)
        ## When I know the bird wasn't foraging (saw bird and did not record foraging activity)
        bird.by.visit.PA$B.not.foraging <-
            ifelse(bird.by.visit.PA$B.foraging == "0" &
                       bird.by.visit.PA$detected.seen == "1",
                   "1",
                   "0")
        ## When I didn't see the bird, and therefore didn't know if it was foraging or not.
        bird.by.visit.PA$B.unknown.foraging <-
            ifelse(bird.by.visit.PA$B.foraging == "0" &
                       bird.by.visit.PA$detected.seen == "0",
                   "1",
                   "0")
        ## Bird engaging in another activity, such as singing.
        bird.by.visit.PA$B.other.activity <-
            apply(X = bird.by.visit.PA[, 10:12], MARGIN = 1, FUN = max)
        ## PA column for basis of PA matrify
        bird.by.visit.PA$PA <-
            ifelse(bird.by.visit.PA$sweeps.present > 0, 1, 0)

        str(bird.by.visit.PA)
        summary(bird.by.visit.PA)



    ## Now create another group based on bird.by.visit.PA grouped by site name and
    ## text.spp (like incidence is) in order to determine # times seen foraging.



        species.foraging <-
            group_by(bird.by.visit.PA, site.name, text.spp) %>%
            summarise(
                proportion.foraging = sum(sweeps.foraging > 0) / 8,
                proportion.ofseen.foraging = sum(sweeps.foraging > 0) /
                    sum(detected.seen == 1),
                string.sweeps.foraging = paste(sweeps.foraging, collapse = ", "),
                string.visits.seen = paste(as.character(unique(visit.number)),
                                           collapse = ", "),
                PA.foraging = ifelse(proportion.foraging > 0, 1, 0)
            )

        ## Since foraging requires sight detection, it might be worthwhile
    ## to convert proportion foraging to a 0/1 scale.



## Matrify bird data


        matrify.incidence <-
            ez.matrify(
                filename = incidence,
                species.name = "text.spp",
                site.name = "site.name",
                abundance = "incidence"
            )
        matrify.incidence <-
            matrify.incidence[, -which(
                colSums(matrify.incidence) == 0 |
                    colnames(matrify.incidence) %in%
                    c("sharp.shinned.hawk", "coopers.hawk", "unknown")
            )]

        sp.seen_2.clean <-
            colnames(matrify.incidence[, -which(colSums(matrify.incidence) < 0.250 |
                                                    colnames(matrify.incidence) == "no.birds")])
        matrify.incidence_2.clean <-
            matrify.incidence[, which(colnames(matrify.incidence) %in% sp.seen_2.clean)]


        matrify.foraging <-
            ez.matrify(
                filename = species.foraging,
                species.name = "text.spp",
                site.name = "site.name",
                abundance = "proportion.foraging"
            )
        matrify.foraging <-
            matrify.foraging[, -which(
                colSums(matrify.foraging) == 0 |
                    colnames(matrify.foraging) %in%
                    c("sharp.shinned.hawk", "coopers.hawk", "no.birds")
            )]


        bird.guilds <-
            read.csv("../../DataRepository/BirdData/birdguildinfo.csv",
                     stringsAsFactors = FALSE)


        ## species richness and sub-sets for the different guilds.

        # appear to be this shape: plot(c(1:10), 1 - c(1:10)^2)


        # Diet
        matrify.omnivore_2.clean <-
            matrify.incidence_2.clean[, which(colnames(matrify.incidence_2.clean) %in%
                                                  bird.guilds$species.name[which(bird.guilds$diet == "omnivore")])]
        matrify.grainivore_2.clean <-
            matrify.incidence_2.clean[, which(colnames(matrify.incidence_2.clean) %in%
                                                  bird.guilds$species.name[which(bird.guilds$diet == "grainivore")])]
        matrify.insectivore_2.clean <-
            matrify.incidence_2.clean[, which(colnames(matrify.incidence_2.clean) %in%
                                                  bird.guilds$species.name[which(bird.guilds$diet == "insectivore")])]
        # Foraging Substrate
        matrify.ground_2.clean <-
            matrify.incidence_2.clean[, which(colnames(matrify.incidence_2.clean) %in%
                                                  bird.guilds$species.name[which(bird.guilds$foraging.substrate ==
                                                                                     "ground")])]
        matrify.treshrb_2.clean <-
            matrify.incidence_2.clean[, which(colnames(matrify.incidence_2.clean) %in%
                                                  bird.guilds$species.name[which(bird.guilds$foraging.substrate ==
                                                                                     "trees.shrubs")])]
        matrify.generalist_2.clean <-
            matrify.incidence_2.clean[, which(colnames(matrify.incidence_2.clean) %in%
                                                  bird.guilds$species.name[which(bird.guilds$foraging.substrate ==
                                                                                     "generalist")])]

        # Forest Preference
        matrify.mixed_2.clean <-
            matrify.incidence_2.clean[, which(colnames(matrify.incidence_2.clean) %in%
                                                  bird.guilds$species.name[which(bird.guilds$forest.preference ==
                                                                                     "mixed")])]
        matrify.conifer_2.clean <-
            matrify.incidence_2.clean[, which(colnames(matrify.incidence_2.clean) %in%
                                                  bird.guilds$species.name[which(bird.guilds$forest.preference ==
                                                                                     "conifer")])]
        matrify.open_2.clean <-
            matrify.incidence_2.clean[, which(colnames(matrify.incidence_2.clean) %in%
                                                  bird.guilds$species.name[which(bird.guilds$forest.preference ==
                                                                                     "open.none")])]


        guild.eff.sp.rich_2.clean <-
            tibble(
                site.name = rownames(matrify.incidence_2.clean),
                omnivore.esr = exp(diversity(matrify.omnivore_2.clean)),
                grainivore.esr = exp(diversity(matrify.grainivore_2.clean)),
                insectivore.esr = exp(diversity(matrify.insectivore_2.clean)),

                ground.for.esr = exp(diversity(matrify.ground_2.clean)),
                treshrb.for.esr = exp(diversity(matrify.treshrb_2.clean)),
                generalist.for.esr = exp(diversity(matrify.generalist_2.clean)),

                mixed.esr = exp(diversity(matrify.mixed_2.clean)),
                conifer.esr = exp(diversity(matrify.conifer_2.clean)),
                open.esr = exp(diversity(matrify.open_2.clean))
            )










## -- SITE VARIABLES ----------------------------------------------------------------------

## Matrified tree and shrub data from Vegetation analysis included in
## source() from vegetation_data_processing.R

## Site covariates also included.

    vegetation.clusters <- read.csv("../Vegetation/vegetationClusters.csv")[-1]

## Create the visit level table with information about weather etc. per visit.


        bird.visit.info <-
            group_by(bird.data.onsite, site.name, visit.number) %>%
            summarise(
                site.name.visit = paste(unique(site.name), unique(visit.number), sep = "."),
                visit.date = max(date),
                noise.level = paste(unique(noise.level), collapse = ", "),
                loud.noise = max(ifelse(grepl(
                    noise.level, pattern = "loud"
                ), 1, 0)),
                visit.start = min(start),
                visit.end = max(end),
                number.visit.sweeps = max(as.numeric(sweep.number), na.rm = TRUE),
                visit.sweeps.length = as.numeric(number.visit.sweeps) * 20,
                cloud.all = paste(unique(cloud), collapse = ", "),
                fog = max(ifelse(cloud %in% c("fog", "light.fog"), 1 , 0), na.rm = TRUE),
                drizzle = max(ifelse(cloud == "drizzle", 1, 0), na.rm = TRUE),
                clouds.76.100 = max(ifelse(cloud == "76_100", 1, 0), na.rm = TRUE),
                clouds.0.15 = max(ifelse(cloud == "0_15", 1, 0), na.rm = TRUE),
                clouds.16.75 = max(ifelse(cloud %in% c("16_50", "51_71"), 1, 0), na.rm = TRUE),
                Min.Air.Temp.F.min = min(Min.Air.Temp.F),
                Avg.Air.Temp.F.mean = mean(Avg.Air.Temp.F),
                Avg.Air.Temp.F.median = median(Avg.Air.Temp.F),
                Max.Air.Temp.F.max = max(Max.Air.Temp.F),
                Wind.Direction = paste(unique(Wind.Dir.and.Deg), collapse = ", "),
                Speed.MPH.mean = mean(Speed.MPH),
                Speed.MPH.median = median(Speed.MPH),
                Gust.MPH.mean = mean(Gust.MPH),
                Gust.MPH.median = median(Gust.MPH),
                Tot.Precip.in.sum = sum(unique(Tot.Precip.in), na.rm = TRUE),
                Solar.Rad.W_m2.max = max(Solar.Rad.W_m2),
                Solar.Rad.W_m2.mean = mean(Solar.Rad.W_m2),
                Solar.Rad.W_m2.median = median(Solar.Rad.W_m2)
            )

        bird.visit.info$year = as.factor(ifelse(
            bird.visit.info$visit.number %in% c("1", "2", "3", "4"),
            "2014",
            "2015"
        ))

        str(bird.visit.info)
        summary(bird.visit.info)
        colnames(bird.visit.info)

        plot(bird.visit.info$Solar.Rad.W_m2.mean, pch = bird.visit.info$site.name)

## Create detection variables collapsed onto site to attach to
## incidence, based on the first collapse by visit...



        site.detection.info <- group_by(bird.visit.info, site.name) %>%
            summarise(
                dates = paste(unique(visit.date), collapse = ", "),
                proportion.loud.noise = mean(loud.noise),
                visit.start.mean = mean(visit.start),
                visit.end.mean = mean(visit.end),
                total.sweeps = sum(number.visit.sweeps),
                total.sweeps.length = sum(visit.sweeps.length),
                cloud.all = paste(cloud.all, collapse = ", "),
                fog.proportion = mean(fog),
                drizzle.proportion = mean(drizzle),
                clouds.76.100.proportion = mean(clouds.76.100),
                clouds.0.15.proportion = mean(clouds.0.15),
                clouds.16.75.proportion = mean(clouds.16.75),
                Min.Air.Temp.F.8v.mean = mean(Min.Air.Temp.F.min),
                Avg.Air.Temp.F.8v.mean = mean(Avg.Air.Temp.F.mean),
                Avg.Air.Temp.F.8v.median = median(Avg.Air.Temp.F.median),
                Max.Air.Temp.F.8v.mean = mean(Max.Air.Temp.F.max),
                Wind.Directions = paste(Wind.Direction, collapse = ", "),
                Speed.MPH.8v.mean = mean(Speed.MPH.mean),
                Speed.MPH.8v.median = median(Speed.MPH.median),
                Gust.MPH.8v.mean = mean(Gust.MPH.mean),
                Gust.MPH.8v.median = median(Gust.MPH.median),
                Tot.Precip.in.8v.sum = sum(Tot.Precip.in.sum, na.rm = TRUE),
                Solar.Rad.W_m2.8v.mean = mean(Solar.Rad.W_m2.mean),
                Solar.Rad.W_m2.8v.median = median(Solar.Rad.W_m2.median)
            )

        summary(site.detection.info)

        plot(
            site.detection.info$Solar.Rad.W_m2.8v.median,
            site.detection.info$Solar.Rad.W_m2.8v.mean
        )
        abline(a = 0, b = 1)

## Merge so that site visit info is attached to PA info.

    bird.by.visit.PA.info <- merge.data.frame(x = bird.by.visit.PA, y = bird.visit.info,
                                              by.x = c("site.name", "visit.number"),
                                              by.y = c("site.name", "visit.number"),
                                              all.x = TRUE)




## Now that we have a combined info matrix, clean up the global environment.
    remove(bird.data, bird.data.onsite, Dec_15, Dec_14, Feb_16, Feb_15, Jan_16, Jan_15, Weather_Hour
           )



## Add some shrub stats to the DF.








## -- DATA STANDARDIZTION -------------------------------------------------------

## Ocurs in bird_data_analysis:
    # univariate.sprich$exp.3_dead.wood <- univariate.sprich$dead.wood^(1/3)
    # univariate.sprich$exp3_age.2017 <- univariate.sprich$age.2017^3
    # univariate.sprich$exp2_height.m.median <- univariate.sprich$height.m.median^2









## -- WRITE ---------------------------------------------------------------------

    ## Export bird data!


    if (write.all == TRUE) {


        write.csv(incidence, file = "../Birds/bird_incidence.csv")
        write.csv(bird.by.visit.PA.info, file = "../Birds/birdbyvisitPAinfo.csv")
        write.csv(incidence, file = "../Birds/incidence.csv")
        write.csv(species.foraging, file = "../Birds/speciesforaging.csv")
        write.csv(site.detection.info, file = "../Birds/sitedetectioninfo.csv")

}
