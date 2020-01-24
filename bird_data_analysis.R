## Bird analysis to support '03-birdanalysis.Rmd'.

## Data manipulation --> analysis (this doc) --> Rmarkdown/bookdown

## First written 10/06/2017. This analysis code is based on prior documents including
## 'bird_data_analysisX.R' and 'Office_Development_Bird_CommunitiesX.Rmd'.



## -- SETUP -------------------------------------------------------------

date <- format(Sys.Date(), "%m%d%Y")
perm <- 99999

## Load required libraries
    # library(labdsv)
    # library(chron)
    library(vegan)
    # library(cluster)
    # library(indicspecies)
    library(tidyr)
    library(data.table)
    # library(MVN)
    # library(dendextend)
    library(stringr)
    library(dplyr)
    # library(car)
    # library(AICcmodavg)


## Source data + helper R script

    source("bird_data_processing.R")
    source('../../../../RCode/R_Scripts/AICc_PERMANOVA.R')
    source('../../../../RCode/R_Scripts/repeat_multipatt.R')
    source("../../../../RCode/R_Scripts/group_PERMANOVA.R")
    source("../../../../RCode/R_Scripts/AICc_table_generation.R")

matrify.incidence_0.clean <- matrify.incidence[ , !(colnames(matrify.incidence) == "no.birds")]

vegetation.clusters <- vegetation.clusters[order(vegetation.clusters$SiteName), ]



## -- DESCRIPTIVE STATISTICS ---------------------------------------------

    # _2 suffix indicates that bird had to be seen at least twice

    descriptive.bird.bysite <-  tibble(
        site = rownames(matrify.incidence),
        site.area = sample.covariates$acres,
        bird.sp.richness = apply(matrify.incidence_0.clean,
                                 MARGIN = 1, specnumber),
        bird.sp.richness_2 = apply(matrify.incidence_2.clean,
                                   MARGIN = 1, specnumber),
        bird.Shannon = apply(matrify.incidence_0.clean,
                             MARGIN = 1, diversity),
        bird.Shannon_2 = apply(matrify.incidence_2.clean,
                               MARGIN = 1, diversity),
        bird.effective.sp = exp(bird.Shannon),
        bird.effective.sp_2 = exp(bird.Shannon_2),
        visits.no.birds = matrify.incidence$no.birds * 8,
        
        foraging.sp.richness = apply(matrify.foraging,
                                     MARGIN = 1, specnumber),
        foraging.sp.richness_2 = apply(matrify.foraging_2,
                                       MARGIN = 1, specnumber),
        foraging.Shannon = apply(matrify.foraging,
                                 MARGIN = 1, diversity),
        foraging.Shannon_2 = apply(matrify.foraging_2,
                                   MARGIN = 1, diversity),
        foraging.effective.sp = exp(foraging.Shannon),
        foraging.effective.sp_2 = exp(foraging.Shannon_2),
        
        conifer.nat.abundance = sample.covariates$conifer.nat.abundance,
        conifer.nat.dens = sample.covariates$conifer.nat.dens
    )

    temp <-
        summarise_all(descriptive.bird.bysite[, 2:17], tibble::lst(min, max, mean, sd, median)) 

        pretty.descriptive.birds.site <- tibble(
            Metric = colnames(descriptive.bird.bysite)[2:17],
            Minimum = unlist(temp[1:16]) %>% round(2),
            Maximum = unlist(temp[17:32]) %>% round(2),
            Mean = unlist(temp[33:48]) %>% round(2),
            `Standard Deviation` = unlist(temp[49:64]) %>% round(2),
            Median = unlist(temp[65:80]) %>% round(2)
        ) 

        remove(temp)

        pretty.descriptive.birds.site$Metric <- c("Site Area (acres)", "Species Richness",
                                                  "Sp. Richness (>1 sighting)",
                                                  "Shannon Entropy", "Shannon Entropy (>1 sighting)",
                                                  "Effective # Sp.",
                                                  "Effective # Sp. (>1 sighting)",
                                                  " No. Visits with No Birds in >0 Sweeps",
                                                  "Foraging Bird Sp. Richness",
                                                  "Foraging Bird Sp. Richness (>1 sighting)",
                                                  "Foraging Shannon Entropy", "Foraging Shannon (>1 sighting)",
                                                  "Foraging Effective # Sp.",
                                                  "Foraging Effective # Sp. (>1 sighting)",
                                                  "Native Conifer Abundance", "Native Conifer Density (#/acre)")




        descriptive.bird.bysp <- tibble(
            species.name = colnames(matrify.incidence),
            num.sites = apply(matrify.incidence, 2, specnumber),
            incidence.min = apply(matrify.incidence, 2, min),
            incidence.max = apply(matrify.incidence, 2, max),
            incidence.median = apply(matrify.incidence, 2, median))

        temp <- tibble(
            species.name = colnames(matrify.foraging),
            foraging.num.sites = apply(matrify.foraging, 2, specnumber),
            foraging.min = apply(matrify.foraging, 2, min),
            foraging.max = apply(matrify.foraging, 2, max),
            foraging.median = apply(matrify.foraging, 2, median)
        )

        descriptive.bird.bysp <- merge(descriptive.bird.bysp, temp, by = "species.name", all.x = TRUE)
        descriptive.bird.bysp[is.na(descriptive.bird.bysp)] <- 0

        remove(temp)

    bird.species.names <- tibble(common.name = unique(incidence$text.spp),
                               scientific.name = c("Corvus brachyrhynchos", "Turdus migratorius",
                                                   "Calypte anna", "Thryomanes bewickii",
                                                   "Poecile atricapillus", "Poecile rufescens",
                                                   "Junco hyemalis", "Regulus satrapa",
                                                   "No Birds", "Sitta canadensis", "Regulus calendula",
                                                   "Melospiza melodia",
                                                   "Pipilo maculatus", "Cyanocitta stelleri", "unknown",
                                                   "Dendroica coronata auduboni",
                                                   "Patagioenas fasciata", "Certhia americana",
                                                   "Vireo huttoni", "Colaptes auratus",
                                                   "Carduelis pinus", "Sphyrapicus ruber",
                                                   "Dendroica townsendi", "Ixoreus naevius",
                                                   "Dryocopus pileatus", "Psaltriparus minimus",
                                                   "Sturnus vulgaris", "Carpodacus mexicanus",
                                                   "Troglodytes pacificus", "Passerella iliaca",
                                                   "Loxia curvirostra", "Corvus corax",
                                                   "Charadrius vociferus", "Columba livia", "Accipiter striatus",
                                                   "Dryobates pubescens", "Accipiter cooperii",
                                                   "Carpodacus purpureus", "Spinus tristis",
                                                   "Passer domesticus"),
                               pretty.common.name = c("American Crow", "American Robin", "Anna's Hummingbird",
                                                      "Bewick's Wren",
                                                      "Black-capped Chickadee", "Chestnut-backed Chickadee",
                                                      "Dark-eyed Junco",
                                                      "Golden-crowned Kinglet", "No Birds", "Red-breasted Nuthatch",
                                                      "Ruby-crowned Kinglet", "Song Sparrow",
                                                      "Spotted Towhee", "Stellers Jay",
                                                      "Unknown sp.", "Audubon's Warbler",
                                                      "Band-tailed Pigeon", "Brown Creeper",
                                                      "Hutton's Vireo",
                                                      "Northern Flicker", "Pine Siskin",
                                                      "Red-breasted Sapsucker", "Townsend's Warbler",
                                                      "Varied Thrush", "Pileated Woodpecker",
                                                      "Bushtit", "European Starling",
                                                      "House Finch", "Pacific Wren", "Fox Sparrow",
                                                      "Red Crossbill", "Common Raven",
                                                      "Killdeer", "Rock Pigeon", "Sharp-shinned Hawk",
                                                      "Downy Woodpecker",
                                                      "Coopers Hawk", "Purple Finch",
                                                      "American Goldfinch", "House Sparrow"))

    pretty.descriptive.birds.species <- merge(descriptive.bird.bysp, bird.species.names, by.x = "species.name",
                                              by.y = "common.name", all.x = TRUE, sort = FALSE)
    pretty.descriptive.birds.species <- pretty.descriptive.birds.species[, c(11,10,2:9)]

    colnames(pretty.descriptive.birds.species) <- c("Common Name", "Scientific Name",
                                                    "Num Sites Seen", "Min Incidence", "Max Incidence",
                                                    "Median Incidence", "Num Sites Seen Foraging",
                                                    "Min Foraging Incidence",
                                                    "Max Foraging Incidence", "Median Foraging Incidence")

remove(descriptive.bird.bysp)


## -- CORRELATION INDEX / INDICATOR SPECIES ---------------------------------------------

    ## Point biserial correlation coeff.

        tree.cluster.pbiserial <-
            repeat.multipatt(
                matrix.name = matrify.incidence_0.clean,
                phi = FALSE,
                func.name = "r.g",
                cluster.name = vegetation.clusters$tree.cluster.name,
                plot.please = FALSE,
                quiet = FALSE,
                freq.cutoff = .5,
                stat.cutoff = .5
            )

        tree.cluster.pbiserial$clust.type <- "Tree"
        tree.cluster.pbiserial <-
            tree.cluster.pbiserial[order(tree.cluster.pbiserial$mean.stat, decreasing = TRUE), ]

        ## Cluster 1 associated with a whole host of 'forest' species. More
        ## things show up here, probably due to the fact that more information
        ## is held in the incidence matrix rather than the presence/absence
        ## matrix.

        shrub.cluster.pbiserial <-
            repeat.multipatt(
                matrix.name = matrify.incidence_0.clean,
                phi = FALSE,
                func.name = "r.g",
                cluster.name = vegetation.clusters$shrub.cluster.name,
                plot.please = FALSE,
                stat.cutoff = .5,
                freq.cutoff = .5,
                quiet = FALSE
            )

        shrub.cluster.pbiserial$clust.type <- "Shrub"
        shrub.cluster.pbiserial <-
            shrub.cluster.pbiserial[order(shrub.cluster.pbiserial$mean.stat,
                                          decreasing = TRUE), ]

        ## Appears to be some change in multipatt algorithm--AMCR & TOWA aren't
        ## showing up now. 
        
        ## Bird species show up for both groups now. American Crow associated
        ## with Cluster 1 (ornamental) and Northern FLicker and Townsends
        ## Warbler for Cluster 2 (Native). For birds that show up in both
        ## analyses, the results are similar (e.g. TOWA for shrub cluster is .51
        ## vs .54 with essentially identical frequency)

        ## The point biserial can also take negative values, which
        ## indicates that a species 'avoids' a specific area. With only two
        ## sites, this is the -value of the stat value

    pretty.pbiserial <- rbind(tree.cluster.pbiserial, shrub.cluster.pbiserial)
    pretty.pbiserial <- left_join(pretty.pbiserial, bird.guilds,  by = c("species" = "species.name"))
    pretty.pbiserial <- left_join(pretty.pbiserial, bird.species.names, by = c("species" = "common.name"))

    pretty.pbiserial$pval.range <- paste0(pretty.pbiserial$min.p.val, "-", pretty.pbiserial$max.p.val)
    pretty.pbiserial$foraging.substrate <- ifelse(pretty.pbiserial$foraging.substrate == "trees.shrubs",
                                                  "Trees & Shrubs",
                                                  "Ground")
    pretty.pbiserial$forest.preference <- ifelse(pretty.pbiserial$forest.preference == "conifer",
                                                 "Conifer",
                                                 "Open or No Preference")

    colnames(pretty.pbiserial) <- c("species", "Times Significant (out of 100)","Frequency Significant",
                                    "Mean Statistic", "group","Minimum p-value","Maximum p-value",
                                    "All p-values", "Community Type", "Vegetation Type", "Diet", "Foraging Substrate",
                                    "Forest Preference", "Scientific Name", "Bird Species", "p-value Range")



## Repeat for foraging birds.

    ## Point biserial correlation coeff.

        tree.cluster.pbiserial.foraging <-
            repeat.multipatt(
                matrix.name = matrify.foraging,
                phi = FALSE,
                func.name = "r.g",
                cluster.name = vegetation.clusters$tree.cluster.name,
                plot.please = FALSE,
                freq.cutoff = .5,
                stat.cutoff = .5,
                quiet = FALSE
            )
        tree.cluster.pbiserial.foraging$clust.type <- "Tree"
        tree.cluster.pbiserial.foraging <- tree.cluster.pbiserial.foraging[
                                        order(tree.cluster.pbiserial.foraging$groupname,
                                              tree.cluster.pbiserial.foraging$mean.stat,
                                              decreasing = TRUE),]

    ## Cluster 1 associated with a whole host of 'forest' species. Crows show up
    ## in the Ornamental clusters

        shrub.cluster.pbiserial.foraging <-
            repeat.multipatt(
                matrix.name = matrify.foraging,
                phi = FALSE,
                func.name = "r.g",
                cluster.name = vegetation.clusters$shrub.cluster.name,
                plot.please = FALSE,
                stat.cutoff = .5,
                freq.cutoff = .5,
                quiet = FALSE
            )
        # shrub.cluster.pbiserial.foraging$clust.type <- "Shrub"
        # shrub.cluster.pbiserial.foraging <- shrub.cluster.pbiserial.foraging[
        #                                         order(shrub.cluster.pbiserial.foraging$mean.stat,
        #                                                                            decreasing = TRUE),]

        ## Crow shows up again associated with Ornamental. Townsends shows up again
        ## for native, along with anna's hummingbird. They're territorial so not
        ## 100% sure how much is them claiming good habitat vs. nearby.
        
        ## Again the shrub results have changed. Not sure why.

    # Make a pretty table for foraging birds.

        # pretty.pbiserial.foraging <-
        #     rbind(tree.cluster.pbiserial.foraging,
        #           shrub.cluster.pbiserial.foraging)
        pretty.pbiserial.foraging <- tree.cluster.pbiserial.foraging
        pretty.pbiserial.foraging <-
            left_join(pretty.pbiserial.foraging,
                      bird.guilds,
                      by = c("species" = "species.name"))
        pretty.pbiserial.foraging <- left_join(pretty.pbiserial.foraging, bird.species.names, by = c("species" = "common.name"))
        
        pretty.pbiserial.foraging$pval.range <-
            paste0(pretty.pbiserial.foraging$min.p.val,
                   "-",
                   pretty.pbiserial.foraging$max.p.val)

        pretty.pbiserial.foraging$foraging.substrate <-
            ifelse(
                pretty.pbiserial.foraging$foraging.substrate == "trees.shrubs",

                "Trees & Shrubs",
                "Ground"
            )
        pretty.pbiserial.foraging$forest.preference <-
            ifelse(
                pretty.pbiserial.foraging$forest.preference == "conifer",
                "Conifer",
                "Open or No Preference"
            )
        pretty.pbiserial.foraging$diet <-
            ifelse(
                pretty.pbiserial.foraging$diet == "omnivore",
                "Omnivore",
                "Insectivore"
            )

        colnames(pretty.pbiserial.foraging) <-
            c("species", "Times Significant (out of 100)","Frequency Significant",
              "Mean Statistic", "group", "Minimum p-value", "Maximum p-value",
              "All p-values", "Community Type", "Vegetation Type", "Diet", "Foraging Substrate",
              "Forest Preference", "Scientific Name", "Bird Species", "p-value Range")








#
#
plot(sample.covariates$MedianHouseholdIncome_B19013e1, matrify.incidence$brown.creeper)
plot(sample.covariates$landrentperacre, matrify.incidence$brown.creeper)
plot(sample.covariates$impervious.sqft, matrify.incidence$brown.creeper)
plot(sample.covariates$conifer.nat.abundance, matrify.incidence$brown.creeper)
plot(sample.covariates$shrub.nat.dens, matrify.incidence$brown.creeper)


## -- REGRESSION PREPARATION ----------------------------------------------------------

        descriptive.bird.bysite$site
        sample.covariates$SiteName 
        vegetation.clusters$SiteName

        univariate.sprich <-
            tibble(site = descriptive.bird.bysite$site,
                   bird.effective.sp_2 = descriptive.bird.bysite$bird.effective.sp_2) %>%
                   left_join( . , sample.covariates, c("site" = "SiteName")) 

# Outliers
    outlier1 <- par(mfrow = c(4,3), mar = c(2,2,2,1))
    dotchart(univariate.sprich$bird.effective.sp_2, main = "Effective Sp. Richness (>2)")
    dotchart(univariate.sprich$acres, main = "Acres")
    dotchart(univariate.sprich$dissolved_parcel_500m_buffer_impervious_500m_mean,
             main = "500m Impervious Surface (%)")
    dotchart(univariate.sprich$Proportion_ForeignBorn_B99051e5, main = "% Foreign Born")
    dotchart(univariate.sprich$MedianHouseholdIncome_B19013e1, main = "Median Household Income ($)")
    dotchart(univariate.sprich$age.2017, main = "Site Age (baseline 2017)")
    dotchart(univariate.sprich$height.m.median, main = "Median DF Height (m)")
    dotchart(univariate.sprich$dead.wood, main = "All Dead Wood")
    dotchart(univariate.sprich$conifer.nat.abundance, main = "Native Conifer Density")
    dotchart(univariate.sprich$shrub.nat.dens, main = "Native Shrub Density")
    dotchart(univariate.sprich$impervious.sqft, main = "Site Impervious Surface (%)")
    (dev.off())

    # There are potential outliers--particularly Median DF Height (2 sites with
    # no DF) and all dead wood (two sites w/high values).

# Collinearity
    panel.hist <- function(x, ...)
    {
        usr <- par("usr"); on.exit(par(usr))
        par(usr = c(usr[1:2], 0, 1.5) )
        h <- hist(x, plot = FALSE)
        breaks <- h$breaks; nB <- length(breaks)
        y <- h$counts; y <- y/max(y)
        rect(breaks[-nB], 0, breaks[-1], y, col = "gray50", ...)
    }

    panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
    {
        usr <- par("usr"); on.exit(par(usr))
        par(usr = c(0, 1, 0, 1))
        r <- abs(cor(x, y))
        txt <- format(c(r, 0.123456789), digits = digits)[1]
        txt <- paste0(prefix, txt)
        if (missing(cex.cor)) cex.cor <- .9 + .5/strwidth(txt)
        text(0.5, 0.5, txt, cex = cex.cor + r)
    }

    pairs(
        select_if(univariate.sprich, is.numeric),
        panel = panel.smooth,
        upper.panel = panel.cor,
        diag.panel = panel.hist
    )

    # want histograms to be uniformly distributed (across a range of values),
    # correlations between explanatory variables low. Most are fine except for
    # impervious surface on site, which is highly correlated with dead wood and
    # native conifer density. Expect that median DF height is going to be
    # important.




## -- SPECIES RICHNESS REGRESSION ------------------------------------------------------------
    
    detection.vars <- c("proportion.loud.noise", "fog.proportion", "drizzle.proportion",
                        "clouds.76.100.proportion", "Avg.Air.Temp.F.8v.median", "Speed.MPH.8v.median",
                        "Gust.MPH.8v.median", "Tot.Precip.in.8v.sum", "Solar.Rad.W_m2.8v.median")
        
    neighborhood.vars <- c("dissolved_parcel_500m_buffer_impervious_500m_mean",
                           "dissolved_parcel_500m_intersections_NUMPOINTS",
                           "DistrictName",
                           "Proportion_ForeignBorn_B99051e5",
                           "MedianHouseholdIncome_B19013e1",
                           "acres",
                           "age.2017",
                           "TV_mean",
                           "SMV_mean",
                           "landrentperacre")
    
    management.vars <- c("impervious.pct",
                         "dead.wood",
                         "Herbicide",
                         "Fertilizer",
                         "Irrigation",
                         "Mulch",
                         "height.m.median",
                         "stand.predate.development",
                         "tree.nat.dens",
                         "tree.nat.esr",
                         "conifer.nat.dens",
                         "shrub.nat.dens",
                         "shrub.nat.esr")
    
## Site Detection Variables


    detection.uPERM.results <-
        group.univ.PERMANOVA(
            var.names = detection.vars,
            var.table = univariate.sprich,
            var.table.c = "univariate.sprich",
            num.control.vars = 0,
            species.vector = univariate.sprich$bird.effective.sp_2,
            species.vector.c = "univariate.sprich$bird.effective.sp_2",
            by.adonis2 = "terms"
        )

    ## Nothing is significant after, clouds before adjustment.

## Compare community groups.

    veg.group.vars <- c("tree.cluster","shrub.cluster")

    veg.uPERM.results <-
        group.univ.PERMANOVA(
            var.names = veg.group.vars,
            var.table = univariate.sprich,
            var.table.c = "univariate.sprich",
            num.control.vars = 1,
            control.vars = "clouds.76.100.proportion",
            species.vector = univariate.sprich$bird.effective.sp_2,
            species.vector.c = "univariate.sprich$bird.effective.sp_2",
            by.adonis2 = "terms"
        )


    # Tree cluster and shrub cluster are both significant before adjustment (particularly tree)

## Neighborhood Variables

    neighborhood.uPERM.results <-
        group.univ.PERMANOVA(
            var.names = neighborhood.vars,
            var.table = univariate.sprich,
            var.table.c = "univariate.sprich",
            num.control.vars = 1,
            control.vars = "clouds.76.100.proportion",
            species.vector = univariate.sprich$bird.effective.sp_2,
            species.vector.c = "univariate.sprich$bird.effective.sp_2",
            by.adonis2 = "terms"
        )
    
    # Nothing significant, even with intersections and district and land rent and SMV/TV means added


## Management and Landscape Variables


    mngmt.uPERM.results <-
        group.univ.PERMANOVA(
            var.names = management.vars,
            var.table = univariate.sprich,
            var.table.c = "univariate.sprich",
            num.control.vars = 1,
            control.vars = "clouds.76.100.proportion",
            species.vector = univariate.sprich$bird.effective.sp_2,
            species.vector.c = "univariate.sprich$bird.effective.sp_2",
            by.adonis2 = "terms"
        )

    ## native conifer density, native tree density, stand predate development,
    ## impervious, native shrub density and species richness, and median height
    ## are significant; most are median height and stands predating.


# Muliple var examination:
    sig.uPERM.results <- rbind(
        detection.uPERM.results[detection.uPERM.results$F.pval < 0.05,],
        neighborhood.uPERM.results[neighborhood.uPERM.results$F.pval < 0.05,],
        mngmt.uPERM.results[mngmt.uPERM.results$F.pval < 0.05,],
        veg.uPERM.results[veg.uPERM.results$F.pval < 0.05,]
    )
    
    multi.uPERM.results <-
        AICc.table.all(
            sig.vars = paste0("sample.covariates$", sig.uPERM.results$var.names[-1]),
            control.var.char = paste0("sample.covariates$", sig.uPERM.results$var.names[1]),
            matrix.char =  "univariate.sprich$bird.effective.sp_2",
            perm = perm,
            comb.incl = 1:3
        )

    # Check interaction effects for top models (none found)

    adonis2(
        univariate.sprich$bird.effective.sp_2 ~
            univariate.sprich$clouds.76.100.proportion *
            univariate.sprich$height.m.median )

    adonis2(
        univariate.sprich$bird.effective.sp_2 ~
            univariate.sprich$clouds.76.100.proportion *
            univariate.sprich$stand.predate.development  *
            univariate.sprich$height.m.median )



## Create a pretty table:
    
    
    pretty.all.uPERM.results <- rbind(
        detection.uPERM.results,
        neighborhood.uPERM.results,
        mngmt.uPERM.results,
        veg.uPERM.results
    )[, c(1, 3, 4:5, 8)]
    
    colnames(pretty.all.uPERM.results) <-
        c("Variable", "Variation Explained per df", "pseudo-_F_", "_p_-value", "AICc")
    
    pretty.all.uPERM.results$pretty.names <- c(
        "Median Average Air Temperature (F)",
        "__Overcast Visits (%; Control Var)__",
        "Drizzle Visits (%)",
        "Fog Visits (%)",
        "Median Wind Gusts (MPH)",
        "Loud Noise Visits (%)",
        "Median Solar Radiation (W/m^2)",
        "Median Wind Speed (MPH)",
        "Total Precipitation (in)",
        "Site Area (acres)",
        
        "Building Age (years)",
        "Mean Impervious within 500m (%)",
        "Major Intersections within 500m",
        "Town (Redmond/Bellevue)",
        "Assessed Value per Acre",
        "Median Household Income (USD)",
        "Foreign Born (%)",
        "Short and Medium Vegetation within 500m (% of area)",
        "Tall Vegetation within 500m (% of area)",
        
        "Native Conifer Density",
        "Dead Wood Abundance",
        "Fertilizer (Y/N)",
        "__Median Douglas Fir Height (m)__",
        "Herbicide (Y/N)",
        "__On Site Impervious (%)__",
        "Irrigation (Y/N)",
        "Mulch (Y/N)",
        "Native Shrub Density",
        "Native Shrub Effective Richness",
        "__Stand Predates Development__",
        "Native Tree Density",
        "Native Tree Effective Richness",
        "Shrub Community Type",
        "Tree Community Type"
    )
    
    ## TO HERE ######
    
    
    pretty.univ.results$`Delta AICc` <-
        as.numeric(pretty.univ.results$AICc) - min(as.numeric(pretty.univ.results$AICc))
    
    adonis2(univariate.sprich$bird.effective.sp_2 ~
                univariate.sprich$stand.predate.development +
                univariate.sprich$height.m.median + univariate.sprich$impervious.sqft)




    ## Compare competing models:

    univ.AICc.tab <-
        AICc.table.all(
            sig.vars = c(
                "univariate.sprich$tree.cluster.name",
                "univariate.sprich$shrub.cluster.name",
                "univariate.sprich$stand.predate.development",
                "univariate.sprich$impervious.sqft",
                "univariate.sprich$height.m.median"
            ),
            matrix.char = "univariate.sprich$bird.effective.sp_2",
            comb.incl = 1:3,
            control.var.char = "univariate.sprich$clouds.76.100.proportion",
            method = "eucli",
            extra.var = TRUE,
            extra.var.char = "univariate.sprich$MedianHouseholdIncome_B19013e1"
        )

    pretty.univ.AICc.tab <-
        univ.AICc.tab[c(which(univ.AICc.tab$`Delta AICc` <= 2), nrow(univ.AICc.tab)),]

    pretty.univ.AICc.tab$Model <-
        c(
            "Overcast and Median DF Height (m)",
            "Overcast, Stand Predating Development, and Median DF Height",
            "Overcast and Median HH Income ($)"
        )




    ## Look at different guild lm for relationships with significant variables (to report)

    guild.cols.totest.c <- c(
        "omnivore.esr",
        "grainivore.esr",
        "insectivore.esr",
        "ground.for.esr",
        "treshrb.for.esr",
        "mixed.esr",
        "conifer.esr",
        "open.esr"
    )
    guild.cols.totest.v <-
        as.data.frame(guild.eff.sp.rich_2.clean[, guild.cols.totest.c])

    guild.imperv.results <-
        data.frame(
            row.names = guild.cols.totest.c,
            corr.SP = rep(NA,
                          times = length(guild.cols.totest.c)),
            var.explnd = rep(NA),
            F.stat = rep(NA),
            prob.pval = rep(NA),
            AICc = rep(NA),
            adj.pval = rep(NA)
        )


    for (i in 1:length(guild.cols.totest.c)) {
        species.vect.c <-
            paste0("guild.eff.sp.rich_2.clean$", guild.cols.totest.c[i])

        PERM.temp <-
            group.univ.PERMANOVA(
                var.names = "impervious.sqft",
                var.table = univariate.sprich,
                var.table.c = "univariate.sprich",
                species.vector = guild.cols.totest.v[i],
                species.vector.c = species.vect.c,
                num.control.vars = 0
            )


        # add calculated values to table:
        guild.imperv.results$corr.SP[i] <-
            cor(univariate.sprich$impervious.sqft,
                guild.cols.totest.v[i],
                method = "sp")
        guild.imperv.results$var.explnd[i] <-
            PERM.temp$var.explnd
        guild.imperv.results$F.stat[i] <- PERM.temp$pseudo.F
        guild.imperv.results$prob.pval[i] <-
            PERM.temp$F.pval
        guild.imperv.results$AICc[i] <- PERM.temp$AICc[[1]]

    }

    guild.imperv.results$adj.pval <-
        p.adjust(p = guild.imperv.results$prob.pval,
                 method = "holm")


    # Stands predating development

            guild.spd.results <-
                data.frame(
                    row.names = guild.cols.totest.c,
                    # corr.SP = rep(NA,
                    #               times = length(guild.cols.totest.c)),
                    var.explnd = rep(NA,
                                     times = length(guild.cols.totest.c)),
                    F.stat = rep(NA),
                    prob.pval = rep(NA),
                    AICc = rep(NA),
                    adj.pval = rep(NA)
                )

            for (i in 1:length(guild.cols.totest.c)) {

                species.vect.c <- paste0("guild.eff.sp.rich_2.clean$",guild.cols.totest.c[i])

                PERM.temp <- group.univ.PERMANOVA(var.names = "stand.predate.development",
                                                  var.table = univariate.sprich,
                                                  var.table.c = "univariate.sprich",
                                                  species.vector = guild.cols.totest.v[i],
                                                  species.vector.c = species.vect.c,
                                                  num.control.vars = 0)


                # add calculated values to table:
                # guild.spd.results$corr.SP[i] <- icc(univariate.sprich$stand.predate.development,
                #                                        guild.cols.totest.v[i], method = "sp")
                guild.spd.results$var.explnd[i] <- PERM.temp$var.explnd
                guild.spd.results$F.stat[i] <- PERM.temp$pseudo.F
                guild.spd.results$prob.pval[i] <- PERM.temp$F.pval
                guild.spd.results$AICc[i] <- PERM.temp$AICc[[1]]

            }

            guild.spd.results$adj.pval <- p.adjust(p = guild.spd.results$prob.pval,
                                                   method = "holm")







        guild.DFheight.results <- data.frame(row.names = guild.cols.totest.c,
                                             corr.SP = rep(NA,
                                                           times = length(guild.cols.totest.c)),
                                             var.explnd = rep(NA),
                                             F.stat = rep(NA),
                                             prob.pval = rep(NA),
                                             AICc = rep(NA),
                                             adj.pval = rep(NA)
        )

        for (i in 1:length(guild.cols.totest.c)) {

            species.vect.c <- paste0("guild.eff.sp.rich_2.clean$",guild.cols.totest.c[i])

            PERM.temp <- group.univ.PERMANOVA(var.names = "height.m.median",
                                              var.table = univariate.sprich,
                                              var.table.c = "univariate.sprich",
                                              species.vector = guild.cols.totest.v[i],
                                              species.vector.c = species.vect.c,
                                              num.control.vars = 0)


            # add calculated values to table:
            guild.DFheight.results$corr.SP[i] <- cor(univariate.sprich$height.m.median,
                                                   guild.cols.totest.v[i], method = "sp")
            guild.DFheight.results$var.explnd[i] <- PERM.temp$var.explnd
            guild.DFheight.results$F.stat[i] <- PERM.temp$pseudo.F
            guild.DFheight.results$prob.pval[i] <- PERM.temp$F.pval
            guild.DFheight.results$AICc[i] <- PERM.temp$AICc[[1]]

        }

        guild.DFheight.results$adj.pval <- p.adjust(p = guild.DFheight.results$prob.pval,
                                                    method = "holm")



##  -- PERMANOVA -------------------------------------------------------------

    ## All bird incidence data

        ## Detection Variables


        PERMANOVA.detection <-
            group.PERMANOVA(
                var.names = detection.vars,
                var.table =  site.detection.info,
                var.table.c = "site.detection.info",
                species.table = matrify.incidence.clean,
                species.table.c = "matrify.incidence.clean",
                num.control.vars = 0,
                by.adonis2 = "terms"
            )
        ## clouds again significant before adjustments

        ## Group Membership


        PERMANOVA.veg.groups <-
            group.PERMANOVA(
                var.names = c(
                    "tree.cluster.name",
                    "shrub.cluster.name",
                    "vegetation_class"
                ),
                var.table =  all.covariates,
                var.table.c = "all.covariates",
                species.table = matrify.incidence.clean,
                species.table.c = "matrify.incidence.clean",
                num.control.vars = 1,
                control.vars = "clouds.76.100.proportion",
                by.adonis2 = "terms"
            )

        ## Tree cluster name is significant. beta dispersion is not!


        ## Socio-economic neighborhood variables...



        PERMANOVA.neighborhood <-
            group.PERMANOVA(
                var.names = neighborhood.vars,
                var.table =  all.covariates,
                var.table.c = "all.covariates",
                species.table = matrify.incidence.clean,
                species.table.c = "matrify.incidence.clean",
                num.control.vars = 1,
                control.vars = "clouds.76.100.proportion",
                by.adonis2 = "terms"
            )

        ## none significant before adjustment


        ## Landscaping and maintenance variables


        PERMANOVA.mangt.land <-
            group.PERMANOVA(
                var.names = mngmt.land.vars,
                var.table =  all.covariates,
                var.table.c = "all.covariates",
                species.table = matrify.incidence.clean,
                species.table.c = "matrify.incidence.clean",
                num.control.vars = 1,
                control.vars = "clouds.76.100.proportion",
                by.adonis2 = "terms"
            )

        ## dead wood, median height, impervious sufrace, shrub ESR,
        ## native shrub density, stands predating development, and
        ## native conifer density, are significant,. note that a number
        ## of the maintenance landscaping variables are not significant
        ## (location) but are significant in spread as with fungi.


        pretty.multi.results <- rbind(
            PERMANOVA.detection[order(rownames(PERMANOVA.detection)), ],
            PERMANOVA.neighborhood[order(rownames(PERMANOVA.neighborhood)), ],
            PERMANOVA.mangt.land[order(rownames(PERMANOVA.mangt.land)), ],
            PERMANOVA.veg.groups[order(rownames(PERMANOVA.veg.groups)), ]
        )[, c(1, 3, 4, 6)]
        colnames(pretty.multi.results) <-
            c("Variation Explained", "pseudo-_F_", "_p_-value", "Dispersion _p_-value")

        pretty.multi.results.sig <- pretty.multi.results[ pretty.multi.results$`_p_-value` < 0.05, ]


## AICc comparison:

            # Overall, six significant variables: tree group, median height, imp surface
            # stands predating and native conifer density. and shrub esr.


        sig.incidence.vars.PERMANOVA <-
            c(
                "all.covariates$tree.cluster.name",
                "all.covariates$stand.predate.development",
                "all.covariates$tree.conifer.nat.dens",
                "all.covariates$impervious.sqft",
                "all.covariates$dead.wood",
                "all.covariates$height.m.median",
                "all.covariates$shrub.esr"
            )


        incidence.AICc.tab <-
            AICc.table.all(
                sig.vars = sig.incidence.vars.PERMANOVA,
                matrix.char = "matrify.incidence.clean",
                comb.incl = 1:4,
                control.var.char = "all.covariates$clouds.76.100.proportion",
                extra.var = TRUE,
                extra.var.char = "all.covariates$MedianHouseholdIncome_B19013e1"
            )

## Create a pretty table: Model (name), Delta AICc, Pseudo-F, p-value


        pretty.incidence.AICc <-
            incidence.AICc.tab[c(which(incidence.AICc.tab$`Delta AICc` <= 2),
                                 1,
                                 nrow(incidence.AICc.tab)) ,]


        pretty.incidence.AICc$Model <-
            c(
                "Overcast and Native Conifer Density",
                "Overcast and Impervious (on site)",
                "Overcast and Median DF Height (m)",
                "Overcast, Tree Community Cluster, and Median DF Height",
                "Overcast, Native Conifer Density, and Median DF Height",
                "Overcast, Impervious, and Median DF Height",
                "Overcast, Dead Wood, and Median DF Height",
                "Overcast, Median DF Height, and Native Shrub Effective Richness",
                "Proportion of Visits Overcast",
                "Overcast and Median HH Income ($)"
            )






        end <- nrow(pretty.incidence.AICc)

        pretty.incidence.AICc <-
            pretty.incidence.AICc[, c(6, 2, 7, 8, 3:5)]

        pretty.incidence.AICc[-1] <-
            round(pretty.incidence.AICc[-1], 3)
        #
        #
        # pretty.incidence.AICc$`Delta AICc`[end] <-
        #     paste0("> ", pretty.incidence.AICc$`Delta AICc`[end])
        # pretty.incidence.AICc$`Pseudo-_F_`[end] <- paste0("< ",
        #                                                   pretty.incidence.AICc$`Pseudo-_F_`[end])
        # pretty.incidence.AICc$`Adjusted p-value` <-
        #     p.adjust(pretty.incidence.AICc$`p-value`,
        #              method = "holm")
        #
        #
        # pretty.incidence.AICc$`Adjusted p-value`[end] <-
        #     paste0("> ",
        #            pretty.incidence.AICc$`Adjusted p-value`[end])
        # pretty.incidence.AICc$`p-value`[end] <- paste0("> ",
        #                                                pretty.incidence.AICc$`p-value`[end])
        # pretty.incidence.AICc$`Relative Likelihood`[end] <-
        #     paste0("< ",
        #            pretty.incidence.AICc$`Relative Likelihood`[end])
        #






## Foraging incidence data

    ## Detection Variables


        PERMANOVA.foraging.detection <-
            group.PERMANOVA(
                var.names = detection.vars,
                var.table =  site.detection.info,
                var.table.c = "site.detection.info",
                species.table = matrify.foraging,
                species.table.c = "matrify.foraging",
                num.control.vars = 0,
                by.adonis2 = "terms"
            )
        ## none significant after adjustment (clouds significant before adjusting tho)


    ## Group Membership


        PERMANOVA.foraging.veg.groups <- group.PERMANOVA(
            var.names = c(
                "tree.cluster.name",
                "shrub.cluster.name",
                "vegetation_class"
            ),
            var.table =  all.covariates,
            var.table.c = "all.covariates",
            species.table = matrify.foraging,
            species.table.c = "matrify.foraging",
            num.control.vars = 1,
            control.vars = "clouds.76.100.proportion",
            by.adonis2 = "terms"
        )

        ## Tree cluster name is significant. beta dispersion is not!
        ## (somehow now veg class is not. idk.)

        ## Socio-economic neighborhood variables...


        PERMANOVA.foraging.neighborhood <-
            group.PERMANOVA(
                var.names = neighborhood.vars,
                var.table =  all.covariates,
                var.table.c = "all.covariates",
                species.table = matrify.foraging,
                species.table.c = "matrify.foraging",
                num.control.vars = 1,
                control.vars = "clouds.76.100.proportion",
                by.adonis2 = "terms"
            )
        ## none significant before or after adjustment. SMV is the only one that is close (.07)


        ## Landscaping and maintenance variables


        PERMANOVA.foraging.mangt.land <-
            group.PERMANOVA(
                var.names = mngmt.land.vars,
                var.table =  all.covariates,
                var.table.c = "all.covariates",
                species.table = matrify.foraging,
                species.table.c = "matrify.foraging",
                num.control.vars = 1,
                control.vars = "clouds.76.100.proportion",
                by.adonis2 = "terms"
            )

        ## dead wood, median height, stands predating and tree conifer density are
        ## all significant. So is shrub effective species richness and
        ## impervious surface. shrub density

        # only significant difference is that for foraging birds vegetation
        # class is also significant. Doesn't make sense to include both in
        # models though...


        ## AICc comparison (foraging):

        ## Overall, seven significant variables: tree group, veg group,
        ## median height, stands predating, native conifer density,
        ## impervious surface and shrub esr.

        sig.foraging.vars.PERMANOVA <- c(
            "all.covariates$tree.cluster.name",
            "all.covariates$stand.predate.development",
            "all.covariates$tree.conifer.nat.dens",
            "all.covariates$height.m.median",
            "all.covariates$impervious.sqft",
            "all.covariates$shrub.esr",
            "all.covariates$dead.wood",
            "all.covariates$shrub.nat.dens"
        )

        foraging.AICc.tab <-
            AICc.table.all(
                sig.vars = sig.foraging.vars.PERMANOVA,
                matrix.char = "matrify.foraging",
                comb.incl = 1:5,
                extra.var = TRUE,
                control.var.char = "all.covariates$clouds.76.100.proportion",
                extra.var.char = "all.covariates$MedianHouseholdIncome_B19013e1"
            )






        ## Create a pretty table: Model (name), Delta AICc, Pseudo-F, p-value

        pretty.foraging.AICc <-
            foraging.AICc.tab[c(which(foraging.AICc.tab$`Delta AICc` <= 2),
                                1,
                                (nrow(foraging.AICc.tab))),]

        pretty.foraging.AICc$Model <- c(
            "Overcast and Stand Predates Development",
            "Overcast and Native Conifer Density",
            "Overcast and Median DF Height (m)",
            "Overcast, Tree Community Cluster, and Median DF Height",
            "Overcast, Stand Predates Development, and Median DF Height",
            "Overcast, Median DF Height, and Native Shrub Effective Richness",

            "Overcast, Median DF Height, and Dead Wood Abundance",
            "Overcast",
            "Overcast and Socio-economic Variables (Median HH Income, $)"
        )



        pretty.foraging.AICc <-
            pretty.foraging.AICc[, c(6, 2, 7, 8, 3:5)]
        pretty.foraging.AICc[-1] <-
            round(pretty.foraging.AICc[-1], 3)

        # end <- nrow(pretty.foraging.AICc)
        #
        # pretty.foraging.AICc[-1] <-
        #     round(pretty.foraging.AICc[-1], 4)
        # pretty.foraging.AICc$`Delta AICc`[end] <-
        #     paste0("> ",
        #            pretty.foraging.AICc$`Delta AICc`[end])
        # pretty.foraging.AICc$`Pseudo-_F_`[end] <-
        #     paste0("< ",
        #            pretty.foraging.AICc$`Pseudo-_F_`[end])
        # pretty.foraging.AICc$`Adjusted p-value` <- p.adjust(pretty.foraging.AICc$`p-value`,
        #                                                             method = "holm")
        # pretty.foraging.AICc$`Adjusted p-value`[end] <-
        #     paste0("> ",
        #            pretty.foraging.AICc$`Adjusted p-value`[end])
        #
        #
        #
        # pretty.foraging.AICc$`p-value`[end] <- paste0("> ",
        #                                               pretty.foraging.AICc$`p-value`[end])
        # pretty.foraging.AICc$`Relative Likelihood`[end] <-
        #     paste0("< ",
        #            pretty.foraging.AICc$`Relative Likelihood`[end])
        #
        #





## -- NMDS VISUALIZATIONS ------------------------------------------------------------------

        # All birds

        incidence.NMDS <- metaMDS(
            comm = matrify.incidence.clean,
            distance = "bray",
            k = 2,
            trymax = 9999,
            autotransform = FALSE
        )      # initial

        for (i in 1:100) {
            incidence.NMDS <- metaMDS(
                comm = matrify.incidence.clean,
                distance = "bray",
                k = 2,
                trymax = 9999,
                previous.best = incidence.NMDS,
                autotransform = FALSE
            )           # subsequent runs
        }

        plot(incidence.NMDS, display = "sites")
        plot(goodness(incidence.NMDS))
        stressplot(incidence.NMDS)

        #incidence.NMDS$species


        # Foraging birds

        foraging.NMDS <- metaMDS(
            comm = matrify.foraging,
            distance = "bray",
            k = 2,
            trymax = 9999,
            autotransform = FALSE
        )      # initial
        for (i in 1:100) {
            foraging.NMDS <- metaMDS(
                comm = matrify.foraging,
                distance = "bray",
                k = 2,
                trymax = 9999,
                previous.best = foraging.NMDS,
                autotransform = FALSE
            )           # subsequent runs
        }

        plot(foraging.NMDS)
        plot(goodness(foraging.NMDS))
        stressplot(foraging.NMDS)

## -- TITAN Analysis ------------------------------------------------------------

        matrify.incidence.titan <-
            matrify.incidence[, colSums(matrify.incidence > 0) > 3]

        TITAN.inci.height <- titan(txa = matrify.incidence.titan,
                                   env = all.covariates$height.m.median,
                                   pur.cut = .75, rel.cut = .75, ncpus = 2)

        TITAN.inci.tcnd <- titan(txa = matrify.incidence.titan,
                                   env = all.covariates$tree.conifer.nat.dens,
                                   pur.cut = .75, rel.cut = .75, ncpus = 2)

        TITAN.inci.imp <- titan(txa = matrify.incidence.titan,
                                   env = all.covariates$impervious.sqft,
                                   pur.cut = .75, rel.cut = .75, ncpus = 2)

        TITAN.inci.clouds <- titan(txa = matrify.incidence.titan,
                                   env = all.covariates$clouds.76.100.proportion,
                                   pur.cut = .75, rel.cut = .75, ncpus = 2)

        TITAN.inci.sesr <- titan(txa = matrify.incidence.titan,
                                   env = all.covariates$shrub.esr,
                                   pur.cut = .75, rel.cut = .75, ncpus = 2)

        TITAN.inci.snd <- titan(txa = matrify.incidence.titan,
                                   env = all.covariates$shrub.nat.dens,
                                   pur.cut = .75, rel.cut = .75, ncpus = 2)

        TITAN.inci.deadw <- titan(txa = matrify.incidence.titan,
                                   env = all.covariates$dead.wood,
                                   pur.cut = .75, rel.cut = .75, ncpus = 2)



## -- Save everything... ----------------------------------------------------------

save.image(file = paste0("../Birds/bird_data_analysis_complete_", date, ".RData"))

