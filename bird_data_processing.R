## This script takes care of bird data manipulation and processing for
## 'bird_data_analysis.R'.

## Data manipulation (this doc) --> analysis --> Rmarkdown/bookdown

## Written 10/02/2017. This code is based on splitting out the data manipulation
## and processing contained in previous files including 'bird_data_analysisX.R'
## and 'Office_Development_Bird_CommunitiesX.Rmd'.



## -- SETUP ----------------------------------------------------------------------

## Source helper R codes
source("../../../../RCode/R_Scripts/triplet_fixer.R")
source("VegetationData/vegetation_data_processing_Combined.R")

remove(
    count.tree.type,
    dens.matrify.tree.sp.size,
    dens.matrify.gc,
    DF.site.measures,
    matrify.large.trees,
    matrify.shrub.zones,
    matrify.tree.sp.size,
    remove.raw
)

    ## Make sure variables in same order (aka paranoia)
        management.landscaping <- management.landscaping[order(management.landscaping$Site),]
        sample.covariates <- sample.covariates[order(sample.covariates$SiteName),]
        
        dens.matrify.shrub.site <- dens.matrify.shrub.site[order(rownames(dens.matrify.shrub.site)),]
        dens.matrify.tree.sp.only <- dens.matrify.tree.sp.only[order(rownames(dens.matrify.tree.sp.only)),]
        matrify.gc <- matrify.gc[order(rownames(matrify.gc)),]
        matrify.shrub.site <- matrify.shrub.site[order(rownames(matrify.shrub.site)),]
        matrify.tree.sp.only <- matrify.tree.sp.only[order(rownames(matrify.tree.sp.only)),]
        
    # rename
        dens.shrub.site <- dens.matrify.shrub.site
        dens.tree.sp.only <- dens.matrify.tree.sp.only
        remove(dens.matrify.shrub.site, dens.matrify.tree.sp.only)

    ## Set options
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
            read.csv(file = "WeatherData/Hourly_WSU_W21/AWN_Dec_14.csv",
                     stringsAsFactors = F)
        Dec_15 <-
            read.csv(file = "WeatherData/Hourly_WSU_W21/AWN_Dec_15.csv",
                     stringsAsFactors = F)
        Jan_15 <-
            read.csv(file = "WeatherData/Hourly_WSU_W21/AWN_Jan_15.csv",
                     stringsAsFactors = F)
        Jan_16 <-
            read.csv(file = "WeatherData/Hourly_WSU_W21/AWN_Jan_16.csv",
                     stringsAsFactors = F)
        Feb_15 <-
            read.csv(file = "WeatherData/Hourly_WSU_W21/AWN_Feb_15.csv",
                     stringsAsFactors = F)
        Feb_16 <-
            read.csv(file = "WeatherData/Hourly_WSU_W21/AWN_Feb_16.csv",
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

    # Remove intermediate steps.
        remove(Dec_15, Dec_14, Feb_16, Feb_15, Jan_16, Jan_15
        )

## Create ground cover/vegetation covariate data tables
    # Ground cover
        dens.matrify.gc <- matrify.gc %>%
            mutate(site.name = rownames(matrify.gc),
                   dense.veg.pct = dense.veg/total.area,
                   dirt.litter.pct = dirt.litter/total.area,
                   grass.pct = grass/total.area,
                   gravel.pct = gravel/total.area,
                   ivy.pct = ivy/total.area,
                   mulch.pct = mulch/total.area,
                   water.pct = water/total.area,
                   impervious.pct = impervious.sqft/total.area,
                   pervious.pct = pervious.sqft/total.area)
        
        colnames(dens.matrify.gc)[colnames(dens.matrify.gc) == "total.area"] <- "total.area.sqft"
        remove(matrify.gc)
        
    # Vegetation
        management.landscaping[management.landscaping == "unknown"] <- "no.unkn"
        management.landscaping[management.landscaping == "no"] <- "no.unkn"

        
        dens.tree.sp.only.native <-
            dens.tree.sp.only[, colnames(dens.tree.sp.only) %in%
                                     tree.native$tree.species[tree.native$tree.origin == "native"]]
        
        matrify.tree.sp.only.native <-
            matrify.tree.sp.only[, colnames(matrify.tree.sp.only) %in%
                                     tree.native$tree.species[tree.native$tree.origin == "native"]]

        matrify.shrub.site.native <-
            matrify.shrub.site[, colnames(matrify.shrub.site) %in%
                                   shrub.native$shrub.scientific.update[shrub.native$shrub.origin == "native"]]

        dens.shrub.site.native <-
            dens.shrub.site[, colnames(dens.shrub.site) %in%
                                shrub.native$shrub.scientific.update[shrub.native$shrub.origin == "native"]]

        
        
        tree.metrics <- tibble(
            site.abbr = rownames(matrify.tree.sp.only),
            tree.abundance = rowSums(matrify.tree.sp.only),
            tree.nat.abundance = rowSums(matrify.tree.sp.only.native),
            tree.dens = rowSums(dens.tree.sp.only),
            tree.nat.dens = rowSums(dens.tree.sp.only.native),
            tree.nat.esr = exp(diversity(matrify.tree.sp.only.native)),
            conifer.nat.abundance = matrify.tree.sp.only$douglas.fir + 
                matrify.tree.sp.only$western.hemlock +
                matrify.tree.sp.only$western.red.cedar,
            conifer.nat.dens = dens.tree.sp.only$douglas.fir + 
                dens.tree.sp.only$western.hemlock + 
                dens.tree.sp.only$western.red.cedar
        )        
        
        
        shrub.metrics <- tibble(
            site.abbr = rownames(matrify.shrub.site),
            shrub.abundance = rowSums(matrify.shrub.site),
            shrub.nat.abundance = rowSums(matrify.shrub.site.native),
            shrub.dens = rowSums(dens.shrub.site),
            shrub.nat.dens = rowSums(dens.shrub.site.native),
            shrub.nat.esr = exp(diversity(matrify.shrub.site.native))
        )
        

    sample.covariates <- left_join(sample.covariates, dens.matrify.gc, c("SiteName" = "site.name")) %>%
                            left_join( . , management.landscaping, c("SiteName" = "Site")) %>%
                            left_join( . , tree.metrics, c("SiteName" = "site.abbr")) %>%
                            left_join( . , shrub.metrics, c("SiteName" = "site.abbr"))

    
        

## -- BIRD OBSERVATIONS ------------------------------------------------------------------

    ## Pull in raw bird data and turn it into a useable table. need to calculate
    ## unknowns as portion of data: likely need to remove them because no way of
    ## telling if they are the same as other unknowns as there is with shrubs and
    ## trees.
        
        
        bird.data <- read.csv("BirdData_all_07192016.csv",
                              stringsAsFactors = FALSE, na.strings = c("NA", "na"))[-1]
        
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

        bird.data.adjacent <- bird.data[bird.data$on.site == "adjacent", ]
        bird.data <- bird.data[bird.data$on.site == "yes" & !is.na(bird.data$sweep.number),] 
        # remove birds only seen on adjacent etc. properties for further analysis
        
        


    # Merge weather info with bird data.
        bird.data <- merge.data.frame(
            x = bird.data,
            y = Weather_Hour,
            by.x = c('date', 'start.hour'),
            by.y = c('Date', 'hours'),
            all.x = TRUE,
            sort = FALSE
        )

    ## There should be 3343 observations of birds on site; less 209 birds seen
    ## outside of official sweeps. Some of these two groups overlap so should have
    ## 3237 observations



    ## Call out where birds were observed foraging--Note that they must have been SEEN.


        bird.data$did.forage <-
            as.numeric(apply(bird.data[, 16:19], MARGIN = 1, FUN = max))

        ## Note, warning will pop up saying NAs introduced by coercion; this is
        ## because the no.birds have X's instead of numbers

    ## Create unique visit-sweep combo for sweep.number (otherwise it's repeated)
        
        bird.data$visit.sweep <- paste(bird.data$visit.number, bird.data$sweep.number, sep = "-")


## Incidence matrix

    ## This is for both years (or each year by itself). Collapses all data to
    ## give proportion of visits that species was observed per site.

        ## Calculate incidence
        incidence <-
            group_by(bird.data, site.name, text.spp) %>%
            summarise(
                visits.present = paste(as.character(unique(visit.number)), collapse = ", "),
                visits.2014 = sum(unique(visit.number) < 5),
                visits.2015 = sum(unique(visit.number) > 4),
                incidence = length(as.character(unique(visit.number))) /
                    8,
                incidence.2014 = visits.2014 / 4,
                incidence.2015 = visits.2015 / 4,
                sweeps.present = paste(as.character(unique(visit.sweep)), collapse = ", "),
                sweeps.present.count = length(as.character(unique(visit.sweep))),
                ## Note: these just tell us how often the bird was ever observed
                ## on site doing these actions. Number of sweeps isn't
                ## consistent between sites (by design)
                A.vegetation.gleaning = sum(A.vegetation.gleaning == "1"),
                A.ground.foraging = sum(A.ground.foraging == "1"),
                A.eating.berries = sum(A.eating.berries == "1"),
                A.eating.seeds = sum(A.eating.seeds == "1"),
                A.sitting = sum(A.sitting == "1"),
                A.calling = sum(A.calling == "1"),
                A.flitting = sum(A.flitting == "1"),
                detected.seen = sum(detected.seen == "1"),
                detected.call  = sum(detected.call == "1"),
                detected.other.noise = sum(detected.other.noise == "1"),
                detected.recorded.audio = sum(detected.recorded.audio == "1")
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
            group_by(bird.data, site.name, visit.number, text.spp) %>%
            summarise(
                BBL = max(BBL),
                sweeps.max = max(sweep.number),
                sweeps.present = n(),
                sweeps.foraging = sum(did.forage),
                A.vegetation.gleaning = sum(A.vegetation.gleaning == "1"),
                A.ground.foraging = sum(A.ground.foraging == "1"),
                A.eating.berries = sum(A.eating.berries == "1"),
                A.eating.seeds = sum(A.eating.seeds == "1"),
                A.sitting = sum(A.sitting == "1"),
                A.calling = sum(A.calling == "1"),
                A.flitting = sum(A.flitting == "1"),
                detected.seen = sum(detected.seen == "1"),
                detected.call  = sum(detected.call == "1"),
                detected.other.noise = sum(detected.other.noise == "1"),
                detected.recorded.audio = sum(detected.recorded.audio == "1"),
                notes = paste(unique(notes, na.rm = TRUE), collapse = ", ")
            )
        ## Note: here, vegetation gleaning etc. tell us if the bird was
        ## observed doing that activity on that visit!


## Need to add a few variables to describe aggregate foraging behavior, etc.

        A.foraging.cols <- c("A.vegetation.gleaning", "A.ground.foraging",
                             "A.eating.berries", "A.eating.seeds" )
        
        ## When I know the bird was foraging (recorded)
            bird.by.visit.PA$B.foraging <-
                apply(X = bird.by.visit.PA[, A.foraging.cols], MARGIN = 1, FUN = max)
            bird.by.visit.PA$B.foraging <- ifelse(bird.by.visit.PA$B.foraging > 0,
                                                  1,
                                                  0)
            
            
        ## When I know the bird wasn't foraging (saw bird and did not record foraging activity)
            bird.by.visit.PA$B.not.foraging <-
                ifelse(bird.by.visit.PA$B.foraging == 0 &
                           bird.by.visit.PA$detected.seen > 0,
                       1,
                       0)
            
            
        ## When I didn't see the bird, and therefore didn't know if it was foraging or not.
            bird.by.visit.PA$B.unknown.foraging <-
                ifelse(bird.by.visit.PA$B.foraging == 0 &
                           bird.by.visit.PA$detected.seen == 0,
                       1,
                       0)
            

        ## PA column for basis of PA matrify
            bird.by.visit.PA$PA <-
                ifelse(bird.by.visit.PA$sweeps.present > 0,
                       1,
                       0)
    
            str(bird.by.visit.PA)
            summary(bird.by.visit.PA)



    ## Now create another group based on bird.by.visit.PA grouped by site name and
    ## text.spp (like incidence is) in order to determine # times seen foraging.

        bird.foraging.incidence <-
            group_by(bird.by.visit.PA, site.name, text.spp) %>%
            summarise(
                proportion.foraging = sum(sweeps.foraging > 0) / 8,
                sweeps.foraging = paste(sweeps.foraging, collapse = ", "),
                visits.seen = paste(as.character(unique(visit.number)),
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
            matrify.incidence[ , 
                !(colnames(matrify.incidence) %in%
                    c("sharp.shinned.hawk", "coopers.hawk", "unknown")) ]

        sp.seen_2.clean <-
            colnames(matrify.incidence[, colSums(matrify.incidence) > 0.125 &
                                                    colnames(matrify.incidence) != "no.birds"] )
        matrify.incidence_2.clean <- matrify.incidence[, sp.seen_2.clean]


        matrify.foraging <-
            ez.matrify(
                filename = bird.foraging.incidence,
                species.name = "text.spp",
                site.name = "site.name",
                abundance = "proportion.foraging"
            )
        
        
        matrify.foraging <-
            matrify.foraging[, 
                colSums(matrify.foraging) > 0 &
                    !(colnames(matrify.foraging) %in%
                    c("sharp.shinned.hawk", "coopers.hawk", "no.birds")
            )]
        
        matrify.foraging_2 <- matrify.foraging[, colnames(matrify.foraging) %in% sp.seen_2.clean]


        bird.guilds <-
            read.csv("birdguildinfo.csv",
                     stringsAsFactors = FALSE)


        guild.eff.sp.rich_2.clean <-
            tibble(
                site.name = rownames(matrify.incidence_2.clean),
                omnivore.esr = exp(diversity( matrify.incidence_2.clean[, colnames(matrify.incidence_2.clean) %in%
                                                 bird.guilds$species.name[bird.guilds$diet == "omnivore"]]
                                                 )),
                grainivore.esr = exp(diversity(matrify.incidence_2.clean[, colnames(matrify.incidence_2.clean) %in%
                                                 bird.guilds$species.name[bird.guilds$diet == "grainivore"]]
                                               )),
                insectivore.esr = exp(diversity(matrify.incidence_2.clean[, colnames(matrify.incidence_2.clean) %in%
                                                 bird.guilds$species.name[bird.guilds$diet == "insectivore"]]
                                                )),
                
                ground.for.esr = exp(diversity( matrify.incidence_2.clean[, colnames(matrify.incidence_2.clean) %in%
                                                bird.guilds$species.name[bird.guilds$foraging.substrate == "ground"]]
                                                )),
                treshrb.for.esr = exp(diversity(matrify.incidence_2.clean[, colnames(matrify.incidence_2.clean) %in%
                                                bird.guilds$species.name[bird.guilds$foraging.substrate ==
                                                                                   "trees.shrubs"]]
                                                )),
                generalist.for.esr = exp(diversity(matrify.incidence_2.clean[, colnames(matrify.incidence_2.clean) %in%
                                                bird.guilds$species.name[bird.guilds$foraging.substrate == "generalist"]]
                                                )),

                mixed.esr = exp(diversity(matrify.incidence_2.clean[, colnames(matrify.incidence_2.clean) %in%
                                                bird.guilds$species.name[bird.guilds$forest.preference == "mixed"]]
                                          )),
                conifer.esr = exp(diversity(matrify.incidence_2.clean[, colnames(matrify.incidence_2.clean) %in%
                                                bird.guilds$species.name[bird.guilds$forest.preference == "conifer"]]
                                            )),
                open.esr = exp(diversity( matrify.incidence_2.clean[, colnames(matrify.incidence_2.clean) %in%
                                                bird.guilds$species.name[bird.guilds$forest.preference == "open.none"]]
                                          ))
            )









## -- SITE VARIABLES ----------------------------------------------------------------------

## Matrified tree and shrub data from Vegetation analysis included in
## source() from vegetation_data_processing.R

## Site covariates also included.

    vegetation.clusters <- read.csv("VegetationData/vegetationClusters.csv", stringsAsFactors = F)[-1]

## Create the visit level table with information about weather etc. per visit.


        bird.visit.info <-
            group_by(bird.data, site.name, visit.number) %>%
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
                clouds.16.75 = max(ifelse(cloud %in% c("16_50", "51_71"), 1, 0), na.rm = TRUE),
                clouds.0.15 = max(ifelse(cloud == "0_15", 1, 0), na.rm = TRUE),
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

## Create detection variables collapsed onto site to attach to incidence

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
                                              all.x = TRUE,
                                              sort = FALSE)


# And merge with sample variables:

        sample.covariates <- left_join( sample.covariates , vegetation.clusters[ , -6], c("SiteName" = "SiteName")) %>%
                                left_join( . , site.detection.info, c("SiteName" = "site.name"))
    
    

## Now that we have a combined info matrix, clean up the global environment.

    remove(Weather_Hour, A.foraging.cols)

    remove(dens.matrify.gc, management.landscaping, tree.metrics, shrub.metrics, site.detection.info)
    





## -- WRITE ---------------------------------------------------------------------

    ## Export bird data!


    if (write.all == TRUE) {


        write.csv(incidence, file = "bird_incidence.csv")
        write.csv(bird.by.visit.PA.info, file = "birdbyvisitPAinfo.csv")
        write.csv(bird.foraging.incidence, file = "bird_foraging_incidence.csv")
        write.csv(site.detection.info, file = "sitedetectioninfo.csv")

}
