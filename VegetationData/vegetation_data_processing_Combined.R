## This script takes care of vegetation data manipulation and processing for
## 'bird_data_processing.R' & 'bird_data_analysis.R'. Save for making it run
## properly in a subfolder it is the same code as published in Dyson, K., 2019.
## Vegetation communities on commercial developments are heterogenous and
## determined by development and landscaping decisions, not socioeconomics. PloS
## one, 14(9). Available online here:
## https://journals.plos.org/plosone/article?id=10.1371%2Fjournal.pone.0222069

## Data manipulation (this doc) --> analysis --> Rmarkdown/bookdown

## Written 7/12/2017. This code is based on splitting out the data manipulation
## and processing contained in previous files including 'veg_data_analysisVX.R'
## and 'OfficeDevelopmentVegetationvX.Rmd'.


## -- SETUP ----------------------------------------------------------------------

options(scipen = 999)
options(tibble.width = Inf)


    ## Source helper R codes
        source("../../../../RCode/R_Scripts/triplet_fixer.R")


    ## Load necessary libraries
        library(tidyr)
        library(dplyr)
        library(data.table)
        library(vegan)
        library(stringr)

    ## Set variables
        remove.raw <- TRUE
        write.all <- FALSE
        


## -- INDEPENDENT VARIABLES ------------------------------------------------------

    ## Load raw data sources

        population.covariates <- read.csv("VegetationData/popcovariateswtax.csv", stringsAsFactors = F)[-1]
        
        population.covariates$fixed_TotalApprLandVal <-
            ifelse(
                population.covariates$TotalApprLandVal == 0,
                yes = median(population.covariates$TotalApprLandVal),
                no = population.covariates$TotalApprLandVal
            )

        population.covariates$landrentperacre <-
            (population.covariates$fixed_TotalApprLandVal) /
            population.covariates$acres

        hist(population.covariates$fixed_TotalApprLandVal, 1000)


    ## This is commented out because there is a tax break or something
    ## causing the appraised value to be $1000 for some buildings...
    ## this metric is thus unreliable...

    # population.covariates$apprimps.netsqft <- population.covariates$TotalApprImpsVal/
    #population.covariates$TotalBldgNetSqFt

        sample.covariates <-
            population.covariates[population.covariates$in.sample == "Y" &
                                      population.covariates$SiteName != "FPX", ]
        

    # Independent variables (age, etc.)
        sample.covariates$age.2017 <- (2017 - sample.covariates$YrBuiltAvg)
        population.covariates$age.2017 <- (2017 - population.covariates$YrBuiltAvg)


## Management and landscaping variables:
        
    
    site.maintenance.raw <-
        read.csv("VegetationData/maintenance_data.csv")
    
    # Condense columns from y/n/etc to one column with 4 factors
        site.maintenance <- site.maintenance.raw[, which(
            colnames(site.maintenance.raw) %in%
                c("Site", "Snag.count", "Stump.count", "Log.count", "Notes")
        )]
        
        site.maintenance$dead.wood <-
            site.maintenance$Snag.count +
            site.maintenance$Stump.count +
            site.maintenance$Log.count
        
        ## For some reason you need to put the is.na test first, otherwise
        ## they just always stay NA. Also, technically I could turn this into
        ## a loop but ughhhh this is easy and works haha.

        site.maintenance$Cleanup <- ifelse(
            is.na(site.maintenance.raw$cleanup.unknown),
            "forest",
            ifelse(
                site.maintenance.raw$cleanup.n == 1,
                "no",
                ifelse(
                    site.maintenance.raw$cleanup.unknown == 1,
                    "unknown",
                    ifelse(site.maintenance.raw$cleanup.y == 1, "yes", "")
                )
            )
        )


        site.maintenance$Herbicide <- ifelse(
            is.na(site.maintenance.raw$herbicide.unknown),
            "forest",
            ifelse(
                site.maintenance.raw$herbicide.n == 1,
                "no",
                ifelse(
                    site.maintenance.raw$herbicide.unknown == 1,
                    "unknown",
                    ifelse(site.maintenance.raw$herbicide.y == 1, "yes", "")
                )
            )
        )


        site.maintenance$Fertilizer <- ifelse(
            is.na(site.maintenance.raw$fertilizer.unknown),
            "forest",
            ifelse(
                site.maintenance.raw$fertilizer.n == 1,
                "no",
                ifelse(
                    site.maintenance.raw$fertilizer.unknown == 1,
                    "unknown",
                    ifelse(site.maintenance.raw$fertilizer.y == 1, "yes", "")
                )
            )
        )


        site.maintenance$Insecticide <- ifelse(
            is.na(site.maintenance.raw$insecticide.unknown),
            "forest",
            ifelse(
                site.maintenance.raw$insecticide.n == 1,
                "no",
                ifelse(
                    site.maintenance.raw$insecticide.unknown == 1,
                    "unknown",
                    ifelse(site.maintenance.raw$insecticide.y == 1, "yes", "")
                )
            )
        )



        site.maintenance$Fungicide <- ifelse(
            is.na(site.maintenance.raw$fungicide.unknown),
            "forest",
            ifelse(
                site.maintenance.raw$fungicide.n == 1,
                "no",
                ifelse(
                    site.maintenance.raw$fungicide.unknown == 1,
                    "unknown",
                    ifelse(site.maintenance.raw$fungicide.y == 1, "yes", "")
                )
            )
        )


        site.maintenance$Irrigation <- ifelse(
            is.na(site.maintenance.raw$irrigation.unknown),
            "forest",
            ifelse(
                site.maintenance.raw$irrigation.n == 1,
                "no",
                ifelse(
                    site.maintenance.raw$irrigation.unknown == 1,
                    "unknown",
                    ifelse(site.maintenance.raw$irrigation.y == 1, "yes", "")
                )
            )
        )


        site.maintenance$Mulch <- ifelse(
            is.na(site.maintenance.raw$mulch.unknown),
            "forest",
            ifelse(
                site.maintenance.raw$mulch.n == 1,
                "no",
                ifelse(
                    site.maintenance.raw$mulch.unknown == 1,
                    "unknown",
                    ifelse(site.maintenance.raw$mulch.y == 1, "yes", "")
                )
            )
        )


        site.maintenance$Mushroom <- ifelse(
            is.na(site.maintenance.raw$mushroom.unknown),
            "forest",
            ifelse(
                site.maintenance.raw$mushroom.n == 1,
                "no",
                ifelse(
                    site.maintenance.raw$mushroom.unknown == 1,
                    "unknown",
                    ifelse(site.maintenance.raw$mushroom.y == 1, "yes", "")
                )
            )
        )


    
    
    DF.measurements <-
        read.csv("VegetationData/DF_measurements2.csv",
                 stringsAsFactors = FALSE)[-1]
    
    DF.measurements[which(DF.measurements$Height.m == "ND"), 3:7] <-
        0.0
    
        # An end row gets added for some reason...
    DF.measurements <- DF.measurements[-which(is.na(DF.measurements$Site)), ]
    
    DF.measurements[, 3:7] <-
        replace(DF.measurements[, 3:7], TRUE, lapply(DF.measurements[, 3:7], as.numeric))
    
    DF.site.measures <-
        group_by(DF.measurements[which(DF.measurements$Number %in%
                                           c("DF1", "DF2", "DF3", "DF4", "DF5")) , ],
                 Site) %>%
        summarise(
            height.m.mean = mean(Height.m, na.rm = TRUE),
            height.m.median = median(Height.m, na.rm = TRUE),
            CBH.m.mean = mean(CBH.m, na.rm = TRUE),
            CBH.m.median = median(CBH.m, na.rm = TRUE),
            DBH.in.mean = mean(DBH.in, na.rm = TRUE),
            DBH.in.median = median(DBH.in, na.rm = TRUE),
            stand.predate.development = max(as.character(native.stand.predate))
        ) %>% as.data.frame()


        # DF.site.measures <- rbind(DF.site.measures,
        #                           c("PS_F", rep(0, 7)),
        #                           c("WR_F", rep(0, 7)))

        # DF.site.measures[grepl(DF.site.measures$Site, pattern = "_F") , 2:7] <-
        #     NA

        DF.site.measures <- DF.site.measures[order(DF.site.measures$Site),]

        rownames(DF.site.measures) <-
            DF.site.measures$Site

        management.landscaping <- merge(DF.site.measures, site.maintenance,
                                        by = "Site", all.x = TRUE, sort = FALSE)




    if (remove.raw == TRUE) {
        remove(DF.measurements, site.maintenance.raw, site.maintenance)

    }








## -- TREE DATA ------------------------------------------------------------------

    ## Load raw data sources
        tree.data <- read.csv("VegetationData/treedata.csv",
                              stringsAsFactors = F)[-1]

        tree.native <- read.csv("VegetationData/tree_native.csv", stringsAsFactors = FALSE)
        
        
    ## Update naming etc.
        
        tree.data[tree.data$tree.species == "purple.beech", "tree.species"] <- c("beech")
        
        
        tree.native$tree.scientific.update <- ifelse(test = tree.native$tree.scientific.update == "",
                                                     yes =  tree.native$tree.scientific,
                                                     no = tree.native$tree.scientific.update)
        tree.native$tree.origin.update <- ifelse(test = tree.native$tree.origin.update == "",
                                                 yes =  tree.native$tree.origin,
                                                 no = tree.native$tree.origin.update)
        
        tree.data <- merge(tree.data[ , -7], tree.native,
                            by.x = "tree.species", by.y = "tree.species",
                            all.x = TRUE, sort = FALSE)
        
        
        

    ## Create count matrix

        # Need to sum within sites and tree categories for counts--This is a job
        # for dplyr.

        # Steps:
        # Figure out the unique tree/size combos
        # Count how many there are of each tree/size combo for each site
        # Create a table with: site, tree/size combo, and abundance

        count.tree.sp.size <- group_by(tree.data, site, tree.species, tree.size) %>%
                              summarise(
                                       tree.count = n(),
                                       tree.scientific = max(tree.scientific.update),
                                       tree.common = max(tree.common),
                                       tree.genus = max(genus),
                                       tree.type = max(tree.type)
                                       ) %>%
                              mutate(species.size = paste(tree.species,
                                                          tree.size, sep = ".") )



        count.tree.sp.only <- group_by(tree.data, site, tree.species) %>%
                              summarise(
                                        tree.count = n(),
                                        tree.latin = max(tree.latin),
                                        tree.scientific = max(tree.scientific.update),
                                        tree.common = max(tree.common),
                                        tree.genus = max(genus),
                                        tree.type = max(tree.type)
                                        )



        count.tree.type <- group_by(tree.data, site, tree.type) %>%
                           summarise(tree.type.count = n())


    ## Create species/site abundance matricies:



        matrify.tree.sp.size <- ez.matrify(filename = count.tree.sp.size,
                                           species.name = "species.size",
                                           site.name = 'site',
                                           abundance = "tree.count")

        matrify.tree.sp.only <- ez.matrify(filename = count.tree.sp.only,
                                           species.name = "tree.species",
                                           site.name = 'site',
                                           abundance = "tree.count")



    if (remove.raw == TRUE) {
        remove(tree.data, count.tree.sp.only, count.tree.sp.size)
    }





## Native classifications for trees:

    tree.native <- tree.native[!tree.native$tree.species == "purple.beech", ]
    tree.native$total.abundance <- colSums(matrify.tree.sp.only)

    

## -- SHRUB DATA -----------------------------------------------------------------

    ## Load raw data sources

    shrub.data <- read.csv("VegetationData/shrubdata.csv",
                               stringsAsFactors = F)[-1]

    shrub.native <- read.csv("VegetationData/shrub_native.csv",
                             stringsAsFactors = FALSE)

    
    ## Update naming etc.
    shrub.data[shrub.data$Species.Common == "spiraea", 5] <-
        "spiraea.japonica.gp"
    shrub.data[shrub.data$Species.Taxonomic == "spiraea.sp", 6] <-
        "spiraea.japonica.gp"
    
    shrub.data[shrub.data$Species.Common %in% c("dogwood", "tatarian.dogwood"), 5] <-
        "tatarian.dogwood.gp"
    shrub.data[shrub.data$Species.Taxonomic %in% c("cornus.sp", "cornus.alba"), 6] <-
        "cornus.alba.gp"
    
    shrub.native$shrub.scientific.update <- ifelse(
        test = shrub.native$shrub.scientific.update == "",
        yes =  shrub.native$shrub.scientific,
        no = shrub.native$shrub.scientific.update
    )
    shrub.native$shrub.origin.update <-
        ifelse(
            test = shrub.native$shrub.origin.update == "",
            yes =  shrub.native$shrub.origin,
            no = shrub.native$shrub.origin.update
        )
    
    shrub.data <- merge(shrub.data, shrub.native,
                        by.x = "Species.Taxonomic", by.y = "shrub.species",
                        all.x = TRUE, sort = FALSE)
    

    # don't need to do group_by for the shrubs since I summarized by site/zone
    # already--though do need to for sites only.


        count.shrub.site <- group_by(shrub.data, Site,
                                     Species.Taxonomic) %>%
                            summarise(
                                  shrub.count = sum(Count),
                                  shrub.common = max(Species.Common),
                                  shrub.common.pretty = max(shrub.common.pretty),
                                  shrub.scientific = max(shrub.scientific.update),
                                  shrub.genus = max(Genus)
                                  )



    ## need to get rid of ferns (and DNE/NOSHRUBS for site, since there aren't zones)
        count.shrub.site <- count.shrub.site[count.shrub.site$shrub.count < 999, ]

        count.shrub.zone <- shrub.data[!(shrub.data$Species.Common %in% c("ferns")), ]


    ## remove species with count = 0
        count.shrub.site <- count.shrub.site[count.shrub.site$shrub.count > 0, ]
        count.shrub.zone <- count.shrub.zone[count.shrub.zone$Count > 0, ]


    ## Create species/site matrices:

        matrify.shrub.zones <- ez.matrify(count.shrub.zone,
                                          species.name = "shrub.scientific",
                                          site.name = "Site.Standard.Group",
                                          abundance = "Count")

        ## Note: I considered combining prunus laurocerasus (cherry laurel) and
        ## var. zabeliana (Zabel laurel) because they're technically varieties
        ## of the same species. However, one is a large shrub (frequently used
        ## as 8'+ hedging) and the other is <2' tall. The habitat provided by
        ## the two is thus very different (birds regularly perch and nest in
        ## cherry laurel, where not so much in zabel...). Since the dogwoods are
        ## also split into tree and shrubs, I left it as is.

        matrify.shrub.zones$`Prunus laurocerasus` <-
            matrify.shrub.zones$`Prunus laurocerasus` +
            matrify.shrub.zones$`Prunus laurocerasus var. zabeliana`
        
        matrify.shrub.zones$`Prunus laurocerasus var. zabeliana` <- NULL

        matrify.shrub.site <- ez.matrify(count.shrub.site,
                                         species.name = "shrub.scientific",
                                         site.name = "Site",
                                         abundance = "shrub.count")
        
        matrify.shrub.site$`Prunus laurocerasus` <-
        matrify.shrub.site$`Prunus laurocerasus` +
        matrify.shrub.site$`Prunus laurocerasus var. zabeliana`
        matrify.shrub.site$`Prunus laurocerasus var. zabeliana` <- NULL


    if (remove.raw == TRUE) {
        remove(shrub.data, count.shrub.site, count.shrub.zone)
    }






    ## Native classifications for shrubs:
        # shrub.native <- shrub.native[shrub.native$shrub.scientific.update %in%
        #                                  colnames(matrify.shrub.site) , ]
                                          
      ## Need to make sure shrub.native is in the same order.
        shrub.native <- distinct(shrub.native, shrub.scientific.update, .keep_all = TRUE)
        shrub.native <- shrub.native[-which(shrub.native$shrub.scientific.update == 
                         "Prunus laurocerasus var. zabeliana"), ]
        
        
        shrub.native <- shrub.native[ order(shrub.native$shrub.scientific.update), ]
        matrify.shrub.site <- matrify.shrub.site[ , order(names(matrify.shrub.site))]
        
        shrub.native$shrub.scientific.update == colnames(matrify.shrub.site)
        
        
        shrub.native$total.abundance <- colSums(matrify.shrub.site)



## -- GROUND COVER DATA ----------------------------------------------------------

    ## I cannot find the original code that turned the area based gc data
    ## exported from GIS into the matrified form. NOTE that there was an error
    ## in the previous saved CSV, now called matrifygcdataBAD.csv. The following
    ## gives correct figures (checked against GIS).

    ground.cover <- read.csv("VegetationData/groundcover_area2.csv")

    count.gc <- group_by(ground.cover, Site, CovTyp) %>%
                summarise(
                    Area.sqft = sum(area.sqft)
                )

    matrify.gc <- ez.matrify(filename = count.gc, species.name = "CovTyp",
                             site.name = "Site", abundance = "Area.sqft")
    matrify.gc <- matrify.gc[ , -which(colnames(matrify.gc) == "forest")]

    matrify.gc <- merge(matrify.gc,
                        sample.covariates[ , colnames(sample.covariates) %in%
                                                       c("SiteName", "DissolvedArea")],
                        by.x = 0, by.y = "SiteName",
                        sort = FALSE,
                        all.x = TRUE)
    
            rownames(matrify.gc) <- matrify.gc$Row.names
            matrify.gc$Row.names <- NULL

        matrify.gc$impervious.sqft <- 
            matrify.gc$DissolvedArea - rowSums(matrify.gc[ , 1:7])
        matrify.gc$pervious.sqft <- 
            rowSums(matrify.gc[ , 1:7])

        colnames(matrify.gc) <- c("dense.veg",
                                  "dirt.litter",
                                  "grass",
                                  "gravel",
                                  "ivy",
                                  "mulch",
                                  "water",
                                  "total.area",
                                  "impervious.sqft",
                                  "pervious.sqft")
        
        matrify.gc <- matrify.gc[ , c(1:7, 9:10, 8)]



    if (remove.raw == TRUE) {
        remove(ground.cover, count.gc)
    }







## -- DATA STANDARDIZTION -------------------------------------------------------
 
    ## Wisconsin standardization doesn't make sense for my data: the between
    ## site effort is already standardized. However, looking at both relative (%
    ## of total trees/shrubs per site) and absolute (# per site) provides useful
    ## insight.

    # relativize by TOTAL area--it is the developer's decision on how much
    # imperviuos vs. pervious surface to create and a combination of the
    # developer and property owner on how to plant.

    # Note that these are number of shrubs or trees/acre

        if (length(setdiff(rownames(matrify.shrub.site),
                           sample.covariates$SiteName)) == 0 &
            rownames(matrify.shrub.site)[19] == sample.covariates$SiteName[19]) {
            dens.matrify.shrub.site <-
                matrify.shrub.site / sample.covariates$acres
            
        } else {
            dens.matrify.shrub.site <-
                merge(
                    matrify.shrub.site,
                    sample.covariates[, which(colnames(sample.covariates) %in%
                                                  c("SiteName", "acres"))],
                    by.x = 0,
                    by.y = "SiteName",
                    sort = FALSE,
                    all.x = TRUE
                )
            row.names(dens.matrify.shrub.site) <-
                dens.matrify.shrub.site$Row.names
            dens.matrify.shrub.site$Row.names <- NULL
            dens.matrify.shrub.site <-
                dens.matrify.shrub.site / dens.matrify.shrub.site$acres
            dens.matrify.shrub.site$acres <- NULL

        }
        
        
        
        if (length(setdiff(rownames(matrify.tree.sp.only),
                           sample.covariates$SiteName)) == 0 &
            rownames(matrify.tree.sp.only)[19] == sample.covariates$SiteName[19]) {
            dens.matrify.tree.sp.only <-
                matrify.tree.sp.only / sample.covariates$acres
            
        } else {
            dens.matrify.tree.sp.only <-
                merge(
                    matrify.tree.sp.only,
                    sample.covariates[, which(colnames(sample.covariates) %in%
                                                  c("SiteName", "acres"))],
                    by.x = 0,
                    by.y = "SiteName",
                    sort = FALSE,
                    all.x = TRUE
                )
            row.names(dens.matrify.tree.sp.only) <-
                dens.matrify.tree.sp.only$Row.names
            dens.matrify.tree.sp.only$Row.names <- NULL
            dens.matrify.tree.sp.only <-
                dens.matrify.tree.sp.only / dens.matrify.tree.sp.only$acres
            dens.matrify.tree.sp.only$acres <- NULL
            
        }
        
        dens.matrify.tree.sp.size <-
            matrify.tree.sp.size / sample.covariates$acres
        
        dens.matrify.gc <- matrify.gc / matrify.gc$total.area
        

    ## Transformations necessary for normality...

        ## Tree density:
            ## Quick test for multivariate normality

            # uniPlot((matrify.tree.sp.only[ ,25:35 ]^.5), type = "histogram")


            ## Third root seems to give more normal results... but ONLY for
            ## non-zeros. So many species get messed up, don't transform.

        ## Large trees only:
        
        matrify.large.trees <-
            matrify.tree.sp.size[, grepl(colnames(matrify.tree.sp.size), pattern = ".L")]
       

        ## Shrub density:
        ## 
            ## Quick test for multivariate normality

            # uniPlot((dens.matrify.shrub.site[ ,25:35 ]^.5), type = "histogram")

            ## Third root seems to give more normal results... but  again, ONLY
            ## for non-zeros. So many species get messed up, don't transform.



## -- WRITE ---------------------------------------------------------------------

        sample.covariates <- sample.covariates[ order(sample.covariates$SiteName), ]
        
        matrify.tree.sp.only <- matrify.tree.sp.only[ order(rownames(matrify.tree.sp.only)), ]
        matrify.tree.sp.size <- matrify.tree.sp.size[ order(rownames(matrify.tree.sp.size)), ]
        dens.matrify.tree.sp.only <- dens.matrify.tree.sp.only[ order(rownames(dens.matrify.tree.sp.only)), ]
        
        matrify.gc <- matrify.gc[ order(rownames(matrify.gc)), ]
        dens.matrify.gc <- dens.matrify.gc[ order(rownames(dens.matrify.gc)), ]
        
        matrify.shrub.site <- matrify.shrub.site[ order(rownames(matrify.shrub.site)), ]
        dens.matrify.shrub.site <- dens.matrify.shrub.site[ order(rownames(dens.matrify.shrub.site)), ]
        
    ## Optional: export matrified shrub and trees data; sample site covariates.

    if (write.all == TRUE) {

        write.csv(matrify.tree.sp.only, file = "matrify_tree_species.csv")
        write.csv(matrify.tree.sp.size, file = "matrify_tree_species_size.csv")
        write.csv(dens.matrify.tree.sp.only, file = "dens_matrify_tree_species.csv")

        write.csv(matrify.shrub.site, file = "matrify_shrub_bysite_comb.csv")
        write.csv(matrify.shrub.zones, file = "matrify_shrub_byzone_comb.csv")
        write.csv(dens.matrify.shrub.site, file = "dens_matrify_shrub_bysite_comb.csv")

        write.csv(matrify.gc, file = "matrify_gc.csv")
        write.csv(dens.matrify.gc, file = "dens_matrify_gc.csv")
        
        write.csv(DF.site.measures, file = "DF_site_measures.csv")
        write.csv(management.landscaping, file = "management_landscaping.csv")
        
        write.csv(sample.covariates, "sample_covariates.csv")
        write.csv(shrub.native, "shrub_native_updated.csv")
        write.csv(tree.native, "tree_native_updated.csv")

    }



