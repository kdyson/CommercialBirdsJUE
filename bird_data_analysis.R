## Bird analysis to support '03-birdanalysis.Rmd'.

## Data manipulation --> analysis (this doc) --> Rmarkdown/bookdown

## First written 10/06/2017. This analysis code is based on prior documents including
## 'bird_data_analysisX.R' and 'Office_Development_Bird_CommunitiesX.Rmd'.



## -- SETUP -------------------------------------------------------------

date <- format(Sys.Date(), "%m%d%Y")
perm <- 99999

## Load required libraries
    library(tibble)
    library(vegan)
    library(tidyr)
    library(data.table)
    library(stringr)
    library(dplyr)
    library(TITAN2)
library(sjPlot)


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
        descriptive.bird.bysp$incidence.range <- paste0(descriptive.bird.bysp$incidence.min,
                                                        "-", descriptive.bird.bysp$incidence.max)
        descriptive.bird.bysp$foraging.range <- paste0(descriptive.bird.bysp$foraging.min,
                                                        "-", descriptive.bird.bysp$foraging.max)
        
        
        pretty.guild <- bird.guilds
        pretty.guild$diet <- ifelse(pretty.guild$diet == "omnivore",
                                    "Omnivore",
                                    ifelse(pretty.guild$diet == "grainivore",
                                           "Grainivore",
                                           "Insectivore"))
        
        pretty.guild$foraging.substrate <- ifelse(
            pretty.guild$foraging.substrate == "ground",
            "Ground",
            ifelse(
                pretty.guild$foraging.substrate == "trees.shrubs",
                "Trees & Shrubs",
                "Generalist*"
            )
        )
        pretty.guild$forest.preference <- ifelse(
            pretty.guild$forest.preference == "conifer",
            "Conifer",
            ifelse(
                pretty.guild$forest.preference == "mixed",
                "Mixed",
                "Open or None"
            )
        )
        colnames(pretty.guild) <- c("Species Name", "Diet", "Foraging Substrate", "Forest Preference")

    bird.species.names <- tibble(common.name = unique(incidence$text.spp))
        bird.species.names$common.name <- bird.species.names$common.name[order(bird.species.names$common.name)]

        bird.species.names$`Scientific Name` <-
            c(
                "Corvus brachyrhynchos",
                "Spinus tristis",
                "Turdus migratorius",
                "Calypte anna",
                "Setophaga coronata auduboni",
                "Patagioenas fasciata",
                "Thryomanes bewickii",
                "Poecile atricapillus",
                "Certhia americana",
                "Psaltriparus minimus",
                "Poecile rufescens",
                "Corvus corax",
                "Accipiter cooperii",
                "Junco hyemalis",
                "Dryobates pubescens",
                "Sturnus vulgaris",
                "Passerella iliaca",
                "Regulus satrapa",
                "Carpodacus mexicanus",
                "Passer domesticus",
                "Vireo huttoni",
                "Charadrius vociferus",
                "No Birds",
                "Colaptes auratus",
                "Troglodytes pacificus",
                "Dryocopus pileatus",
                "Spinus pinus",
                "Haemorhous purpureus",
                "Sitta canadensis",
                "Sphyrapicus ruber",
                "Loxia curvirostra",
                "Columba livia",
                "Regulus calendula",
                "Accipiter striatus",
                "Melospiza melodia",
                "Pipilo maculatus",
                "Cyanocitta stelleri",
                "Setophaga townsendi",
                "Unknown spp.",
                "Ixoreus naevius"
            )
        
        bird.species.names$pretty.common.name <-
            c(
                "American Crow",
                "American Goldfinch",
                "American Robin",
                "Anna's Hummingbird",
                "Audubon's (Yellow-rumped) Warbler",
                "Band-tailed Pigeon",
                "Bewick's Wren",
                "Black-capped Chickadee",
                "Brown Creeper",
                "Bushtit",
                "Chestnut-backed Chickadee",
                "Common Raven",
                "Coopers Hawk",
                "Dark-eyed Junco",
                "Downy Woodpecker",
                "European Starling",
                "Fox Sparrow",
                "Golden-crowned Kinglet",
                "House Finch",
                "House Sparrow",
                "Hutton's Vireo",
                "Killdeer",
                "No Birds",
                "Northern Flicker",
                "Pacific Wren",
                "Pileated Woodpecker",
                "Pine Siskin",
                "Purple Finch",
                "Red-breasted Nuthatch",
                "Red-breasted Sapsucker",
                "Red Crossbill",
                "Rock Pigeon",
                "Ruby-crowned Kinglet",
                "Sharp-shinned Hawk",
                "Song Sparrow",
                "Spotted Towhee",
                "Stellers Jay",
                "Townsend's Warbler",
                "Unknown spp.",
                "Varied Thrush"
            )

    pretty.descriptive.birds.species <- merge(descriptive.bird.bysp, bird.species.names, by.x = "species.name",
                                              by.y = "common.name", all.x = TRUE, sort = FALSE) %>%
                                        merge(pretty.guild, by.x = "species.name", by.y = "Species Name",
                                              all.x = TRUE, sort = FALSE)
    pretty.descriptive.birds.species <- pretty.descriptive.birds.species[ , c("pretty.common.name", "Scientific Name",
                                                                              "Diet", "Foraging Substrate",
                                                                              "Forest Preference", "num.sites",
                                                                              "incidence.range", "incidence.median",
                                                                              "foraging.num.sites", "foraging.range",
                                                                              "foraging.median")]
    
    colnames(pretty.descriptive.birds.species) <-
        c(
            "Common Name",
            "Scientific Name",
            "Diet",
            "Foraging Substrate",
            "Forest Preference",
            "Num Sites Seen",
            "Incidence Range",
            "Median Incidence",
            "Num Sites Seen Foraging",
            "Foraging Incidence Range",
            "Median Foraging Incidence"
        )
    
    remove(descriptive.bird.bysp, pretty.guild)
    

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

    # pairs(
    #     select_if(univariate.sprich, is.numeric),
    #     panel = panel.smooth,
    #     upper.panel = panel.cor,
    #     diag.panel = panel.hist
    # )

    # want histograms to be uniformly distributed (across a range of values),
    # correlations between explanatory variables low. Most are fine except for
    # impervious surface on site, which is highly correlated with dead wood and
    # native conifer density. Expect that median DF height is going to be
    # important.

## -- Regression Variables ----------------------------------
    
    detection.vars <-
        c(
            "proportion.loud.noise",
            "fog.proportion",
            "drizzle.proportion",
            "clouds.76.100.proportion",
            "Avg.Air.Temp.F.8v.median",
            "Speed.MPH.8v.median",
            "Gust.MPH.8v.median",
            "Tot.Precip.in.8v.sum",
            "Solar.Rad.W_m2.8v.median"
        )
    
    veg.group.vars <- c("tree.cluster.name","shrub.cluster.name")
    
    neighborhood.vars <- c("dissolved_parcel_500m_buffer_impervious_500m_mean",
                           "dissolved_parcel_500m_intersections_NUMPOINTS",
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
                         
                         "tree.dens",
                         "tree.nat.dens",
                         "tree.nat.esr",
                         "tree.orn.dens",
                         "tree.orn.esr",
                         
                         "conifer.nat.dens",
                         
                         "shrub.dens",
                         "shrub.nat.dens",
                         "shrub.nat.esr",
                         "shrub.orn.dens",
                         "shrub.orn.esr"
                         )

    
    
## -- Eff. Sp. Richness PERMANOVA ------------------------------------------------------------
    
        
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

    plot(univariate.sprich$bird.effective.sp_2 ~ univariate.sprich$tree.orn.dens)
    plot(univariate.sprich$bird.effective.sp_2 ~ univariate.sprich$tree.nat.dens)
    



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
        "Assessed Value per Acre",
        "Median Household Income (USD)",
        "Foreign Born (%)",
        "Short and Medium Vegetation within 500m (% of area)",
        "Tall Vegetation within 500m (% of area)",
        
        "__Native Conifer Density__",
        "Dead Wood Abundance",
        "Fertilizer (Y/N)",
        "__Median Dominant Conifer Height (m)__",
        "Herbicide (Y/N)",
        "__On Site Impervious (%)__",
        "Irrigation (Y/N)",
        "Mulch (Y/N)",
        "Shrub Density",
        "__Native Shrub Density__",
        "__Native Shrub Effective Richness__",
        "Non-native Shrub Density",
        "Non-native Shrub Effective Richness",
        "__Stand Predates Development__",
        "Tree Density",
        "__Native Tree Density__",
        "Native Tree Effective Richness",
        "__Non-native Tree Density__",
        "Non-native Tree Effective Richness",
        "__Shrub Community Type__",
        "__Tree Community Type__"
    )
    

    
    pretty.all.uPERM.results$`Delta AICc` <-
        pretty.all.uPERM.results$AICc - min(pretty.all.uPERM.results$AICc)
    

    
    # Muliple var examination:
    sig.uPERM.results <- rbind(
        detection.uPERM.results[detection.uPERM.results$F.pval < 0.05,],
        neighborhood.uPERM.results[neighborhood.uPERM.results$F.pval < 0.05,],
        mngmt.uPERM.results[mngmt.uPERM.results$F.pval < 0.05,],
        veg.uPERM.results[veg.uPERM.results$F.pval < 0.05,]
    )
    sig.uPERM.results$delta.aic <- sig.uPERM.results$AIC.stat-
        min(sig.uPERM.results$AIC.stat)
    
    multi.uPERM.results <-
        AICc.table.all(
            sig.vars = paste0("sample.covariates$", sig.uPERM.results$var.names[-1]),
            control.var.char = paste0("sample.covariates$", sig.uPERM.results$var.names[1]),
            matrix.char =  "univariate.sprich$bird.effective.sp_2",
            perm = perm,
            comb.incl = 1:4
        )
    
    multi.uPERM.AICc.weights.byvar <- AICc.weights.byvar(
        sig.vars = paste0("sample.covariates$", sig.uPERM.results$var.names[-1]),
        multi.uPERM.results
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
    


    ## Make a table for 1+ variable models


    pretty.AICc.uPERM.tab <-
        multi.uPERM.results[c(1, which(multi.uPERM.results$`Delta AICc` <= 3))
                            , ]

    pretty.AICc.uPERM.tab$Model <-
        c(
            "Proportion of Visits Overcast",
            "Overcast + Median Dom. Conif. Height (m)",
            "Overcast + Median Dom. Conif. Height + Native Conifer Density",
            "Overcast + Median Dom. Conif. Height + Impervious on Site (%)",
            "Overcast + Median Dom. Conif. Height + Native Shrub Density",
            "Overcast + Median Dom. Conif. Height + Stand Predating Development",
            "Overcast + Median Dom. Conif. Height + Native Tree Density",
            "Overcast + Median Dom. Conif. Height + Tree Community Type"
            
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
        guild.eff.sp.rich_2.clean[, guild.cols.totest.c]

    
    
    
    guild.tests <- function(guild.var.name, guild.var.vector) {
        temp <- data.frame(
            row.names = guild.cols.totest.c,
            corr.SP = rep(NA,
                          times = length(guild.cols.totest.c)),
            var.explnd = rep(NA),
            F.stat = rep(NA),
            prob.pval = rep(NA),
            lin.R2 = rep(NA)
        )
        
        for (i in 1:length(guild.cols.totest.c)) {
            species.vect.c <-
                paste0("guild.eff.sp.rich_2.clean$", guild.cols.totest.c[i])
            
            PERM.temp <-
                group.univ.PERMANOVA(
                    var.names = guild.var.name,
                    var.table = univariate.sprich,
                    var.table.c = "univariate.sprich",
                    species.vector = guild.cols.totest.v[i],
                    species.vector.c = species.vect.c,
                    num.control.vars = 1, 
                    control.vars = "clouds.76.100.proportion"
                    
                )
            
            
            # add calculated values to table:
            if(is.numeric(guild.var.vector)){
                            temp$corr.SP[i] <-
                cor(
                    guild.var.vector,
                    guild.cols.totest.v[i],
                    method = "sp"
                )
            }
            else {temp$corr.SP[i] <- "not numeric"}
            
            temp$var.explnd[i] <- PERM.temp$var.explnd
            temp$F.stat[i] <- PERM.temp$pseudo.F
            temp$prob.pval[i] <- PERM.temp$F.pval

            model.lm <- lm(formula = guild.cols.totest.v[[i]] ~ guild.var.vector)
            temp$lin.R2[i] <- summary(model.lm)$r.squared
            
        }
        
        return(temp)
    }
    
    
    

    
## Guilds with median height of dominant trees
    
    guild.height.results <- guild.tests(guild.var.name = "height.m.median",
                                        guild.var.vector = univariate.sprich$height.m.median)
    

# Stands predating development

    guild.spd.results <- guild.tests("stand.predate.development",
                                     univariate.sprich$stand.predate.development)
    
    
# Less significant vars
    
    guild.imperv.results <- guild.tests("conifer.nat.dens",
                                        univariate.sprich$conifer.nat.dens)
    guild.snd.results <- guild.tests("shrub.nat.dens",
                                        univariate.sprich$shrub.nat.dens)
    guild.snesr.results <- guild.tests("shrub.nat.esr",
                                        univariate.sprich$shrub.nat.esr)
    guild.tnd.results <- guild.tests("tree.nat.dens",
                                        univariate.sprich$tree.nat.dens)
    guild.shrubcomm.results <- guild.tests("shrub.cluster.name",
                                        univariate.sprich$shrub.cluster.name)
    guild.treecomm.results <- guild.tests("tree.cluster.name",
                                        univariate.sprich$tree.cluster.name)
     
               








##  -- Incidence PERMANOVA ---------------------------------------------------------

    ## Detection Variables


        detection.mPERM.results <-
            group.PERMANOVA(
                var.names = detection.vars,
                var.table =  sample.covariates,
                var.table.c = "sample.covariates",
                species.table = matrify.incidence_0.clean,
                species.table.c = "matrify.incidence_0.clean",
                num.control.vars = 0,
                by.adonis2 = "terms"
            )
        ## clouds again significant before adjustments

    ## Group Membership

        veg.mPERM.results <-
            group.PERMANOVA(
                veg.group.vars,
                var.table =  sample.covariates,
                var.table.c = "sample.covariates",
                species.table = matrify.incidence_0.clean,
                species.table.c = "matrify.incidence_0.clean",
                num.control.vars = 1,
                control.vars = "clouds.76.100.proportion",
                by.adonis2 = "terms"
            )

        ## Tree community type is significant. beta dispersion is not!


    ## Socio-economic neighborhood variables...


        neighborhood.mPERM.results <-
            group.PERMANOVA(
                var.names = neighborhood.vars,
                var.table =  sample.covariates,
                var.table.c = "sample.covariates",
                species.table = matrify.incidence_0.clean,
                species.table.c = "matrify.incidence_0.clean",
                num.control.vars = 1,
                control.vars = "clouds.76.100.proportion",
                by.adonis2 = "terms"
            )

        ## none significant before adjustment


    ## Landscaping and maintenance variables


        mngmt.mPERM.results <-
            group.PERMANOVA(
                var.names = management.vars,
                var.table =  sample.covariates,
                var.table.c = "sample.covariates",
                species.table = matrify.incidence_0.clean,
                species.table.c = "matrify.incidence_0.clean",
                num.control.vars = 1,
                control.vars = "clouds.76.100.proportion",
                by.adonis2 = "terms"
            )

        ## dead wood, median height, impervious sufrace, shrub ESR,
        ## native shrub density, stands predating development, and
        ## native conifer density, etc. are significant,. note that a number
        ## of the maintenance landscaping variables are not significant
        ## (location) but are significant in spread as with fungi. 


        pretty.all.mPERM.results <- rbind(
            detection.mPERM.results,
            neighborhood.mPERM.results,
            mngmt.mPERM.results,
            veg.mPERM.results
        )[, c(1, 3, 4, 5, 10)]
        colnames(pretty.all.mPERM.results) <-
            c("Variable", "Variation Explained per df", "pseudo-_F_", "_p_-value", "AICc")
        
        sig.mPERM.results <- pretty.all.mPERM.results[ pretty.all.mPERM.results$`_p_-value` < 0.05, ]


## AICc comparison:

    # Overall, nine significant variables: native conifer density, dead
    # wood, median height, impervious on site, native shrub density,
    # native shrub esr, stands predating development, native tree
    # density, tree community


        multi.mPERM.results <-
            AICc.table.all(
                sig.vars = paste0("sample.covariates$",sig.mPERM.results$Variable[-1]),
                matrix.char = "matrify.incidence_0.clean",
                comb.incl = 1:4,
                control.var.char = "sample.covariates$clouds.76.100.proportion",
                perm = perm
            )
        
        multi.mPERM.AICc.weights.byvar <- AICc.weights.byvar(
            sig.vars = paste0("sample.covariates$",sig.mPERM.results$Variable[-1]),
            multi.mPERM.results
        )

## Create a pretty table: Model (name), Delta AICc, Pseudo-F, p-value


        pretty.AICc.mPERM.tab <-
            multi.mPERM.results[c(1, which(multi.mPERM.results$`Delta AICc` <= 2))
                                  ,]
        
        # pretty.AICc.mPERM.tab[ , select_if(is.numeric)] <-
        #     round(pretty.AICc.mPERM.tab[-1], 3)
        
        pretty.AICc.mPERM.tab$Model <-
            c(
                "Proportion of Visits Overcast",
                "Overcast + Native Conifer Density",
                "Overcast + Median Dom. Conif. Height (m)",
                "Overcast + Impervious on Site (%)",
                "Overcast + Native Conifer Density, and Median Dom. Conif. Height",
                "Overcast + Dead Wood + Median Dom. Conif. Height",
                "Overcast + Median Dom. Conif. Height + Impervious on Site",
                "Overcast + Median Dom. Conif. Height + Native Shrub ESR",
                "Overcast + Median Dom. Conif. Height + All Tree Density",
                "Overcast + Median Dom. Conif. Height + Native Tree Density",
                "Overcast + Median Dom. Conif. Height + Non-native Tree Density",
                "Overcast + Median Dom. Conif. Height + Tree Community Type"

            )





## Foraging incidence data

    ## Detection Variables

        detection.mPERM.foraging <-
            group.PERMANOVA(
                var.names = detection.vars,
                var.table =  sample.covariates,
                var.table.c = "sample.covariates",
                species.table = matrify.foraging,
                species.table.c = "matrify.foraging",
                num.control.vars = 0,
                by.adonis2 = "terms"
            )
        ## none significant after adjustment (clouds significant before adjusting tho)


    ## Group Membership

        veg.mPERM.foraging <- group.PERMANOVA(
            var.names = veg.group.vars,
            var.table =  sample.covariates,
            var.table.c = "sample.covariates",
            species.table = matrify.foraging,
            species.table.c = "matrify.foraging",
            num.control.vars = 1,
            control.vars = "clouds.76.100.proportion",
            by.adonis2 = "terms"
        )

        ## Tree cluster name is significant. beta dispersion is not!
        ## Shrubs not significant

    ## Socio-economic neighborhood variables...

        neighborhood.mPERM.foraging <-
            group.PERMANOVA(
                var.names = neighborhood.vars,
                var.table =  sample.covariates,
                var.table.c = "sample.covariates",
                species.table = matrify.foraging,
                species.table.c = "matrify.foraging",
                num.control.vars = 1,
                control.vars = "clouds.76.100.proportion",
                by.adonis2 = "terms"
            )
        ## none significant before or after adjustment. SMV is the only one that is close (.07)


    ## Landscaping and maintenance variables

        mngmt.mPERM.foraging <-
            group.PERMANOVA(
                var.names = management.vars,
                var.table =  sample.covariates,
                var.table.c = "sample.covariates",
                species.table = matrify.foraging,
                species.table.c = "matrify.foraging",
                num.control.vars = 1,
                control.vars = "clouds.76.100.proportion",
                by.adonis2 = "terms"
            )

        ## native conifer density, dead wood abundance, median dom conif height,
        ## impervious, native shurb density, native shrub esr, stands predating
        ## development, native tree density. Native tree esr borderline.

        ## Same as all birds. 

        
        
        pretty.all.mPERM.foraging <- rbind(
            detection.mPERM.foraging,
            neighborhood.mPERM.foraging,
            mngmt.mPERM.foraging,
            veg.mPERM.foraging
        )[, c(1, 3, 4, 5, 10)]
        colnames(pretty.all.mPERM.foraging) <-
            c("Variable", "Variation Explained per df", "pseudo-_F_", "_p_-value", "AICc")
        
        sig.mPERM.foraging <- pretty.all.mPERM.foraging[ pretty.all.mPERM.foraging$`_p_-value` < 0.05, ]
        

    ## AICc comparison (foraging):

        multi.mPERM.foraging <-
            AICc.table.all(
                sig.vars = paste0("sample.covariates$",sig.mPERM.foraging$Variable[-1]),
                matrix.char = "matrify.foraging",
                comb.incl = 1:5,
                extra.var = TRUE, 
                perm = perm,
                control.var.char = "sample.covariates$clouds.76.100.proportion",
                extra.var.char = "sample.covariates$MedianHouseholdIncome_B19013e1"
            )
        
        multi.mPERM.foraging.AICc.weights.byvar <- AICc.weights.byvar(
            sig.vars = paste0("sample.covariates$",sig.mPERM.foraging$Variable[-1]),
            multi.mPERM.foraging
        )


        ## Create a pretty table: Model (name), Delta AICc, Pseudo-F, p-value

        pretty.AICc.mPERM.foraging <-
            multi.mPERM.foraging[c(1, which(multi.mPERM.foraging$`Delta AICc` <= 2))
                                ,]

        
        pretty.AICc.mPERM.foraging$Model <- c(
            "Overcast",
            "Overcast + Native Conifer Density",
            "Overcast + Median Dom. Conif. Height (m)",
            "Overcast + Stand Predates Development",
            "Overcast + Non-native Tree Density",
            "Overcast + Median Dom. Conif. Height + Dead Wood Abundance",
            "Overcast + Median Dom. Conif. Height + Native Shrub Effective Richness",
            "Overcast + Median Dom. Conif. Height + Stand Predates Development",
            "Overcast + Median Dom. Conif. Height + Non-native Tree Density",
            "Overcast + Median Dom. Conif. Height + Tree Community Type"
        )






## -- NMDS VISUALIZATIONS ------------------------------------------------------------------

    # All birds

        incidence.NMDS <- metaMDS(
            comm = matrify.incidence_0.clean,
            distance = "bray",
            k = 2,
            trymax = 9999,
            autotransform = FALSE
        )      # initial

        for (i in 1:100) {
            incidence.NMDS <- metaMDS(
                comm = matrify.incidence_0.clean,
                distance = "bray",
                k = 2,
                trymax = 9999,
                previous.best = incidence.NMDS,
                autotransform = FALSE
            )           # subsequent runs
        }

        plot(incidence.NMDS)
        plot(goodness(incidence.NMDS))
        stressplot(incidence.NMDS)



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

        
## -- CORRELATION INDEX / INDICATOR SPECIES ---------------------------------------------
        
        ## Point biserial correlation coeff.
        
        tree.cluster.pbiserial <-
            repeat.multipatt(
                matrix.name = matrify.incidence_0.clean,
                phi = FALSE,
                func.name = "r.g",
                cluster.name = sample.covariates$tree.cluster.name,
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
                cluster.name = sample.covariates$shrub.cluster.name,
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
                cluster.name = sample.covariates$tree.cluster.name,
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
                cluster.name = sample.covariates$shrub.cluster.name,
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

        
        plot(sample.covariates$MedianHouseholdIncome_B19013e1, matrify.incidence$brown.creeper)
        plot(sample.covariates$landrentperacre, matrify.incidence$brown.creeper)
        plot(sample.covariates$impervious.sqft, matrify.incidence$brown.creeper)
        plot(sample.covariates$conifer.nat.abundance, matrify.incidence$brown.creeper)
        plot(sample.covariates$shrub.nat.dens, matrify.incidence$brown.creeper)
        

    # For stands of trees predating development:
        
        ## Point biserial correlation coeff.
        
        spd.cluster.pbiserial <-
            repeat.multipatt(
                matrix.name = matrify.incidence_0.clean,
                phi = FALSE,
                func.name = "r.g",
                cluster.name = sample.covariates$stand.predate.development,
                plot.please = FALSE,
                quiet = FALSE,
                freq.cutoff = .5,
                stat.cutoff = .5
            )
        
        spd.cluster.pbiserial$group <- ifelse(spd.cluster.pbiserial$group == "yes", "Yes, Stand Predates Development", "No")
        spd.cluster.pbiserial <-
            spd.cluster.pbiserial[order(spd.cluster.pbiserial$mean.stat, decreasing = TRUE), ]
        
        
        
        
        
        
## -- TITAN Analysis ------------------------------------------------------------

        matrify.incidence.titan <-
            matrify.incidence[, colSums(matrify.incidence > 0) > 3]
        
        TITAN.inci.clouds <- titan(txa = matrify.incidence.titan,
                                   env = sample.covariates$clouds.76.100.proportion,
                                   pur.cut = .75, rel.cut = .75, ncpus = 2)

        TITAN.inci.tcnd <- titan(txa = matrify.incidence.titan,
                                 env = sample.covariates$tree.nat.dens,
                                 pur.cut = .75, rel.cut = .75, ncpus = 2)

        TITAN.inci.deadw <- titan(txa = matrify.incidence.titan,
                                  env = sample.covariates$dead.wood,
                                  pur.cut = .75, rel.cut = .75, ncpus = 2)
                
        TITAN.inci.height <- titan(txa = matrify.incidence.titan,
                                   env = sample.covariates$height.m.median,
                                   pur.cut = .75, rel.cut = .75, ncpus = 2)

        TITAN.inci.imp <- titan(txa = matrify.incidence.titan,
                                   env = sample.covariates$impervious.sqft,
                                   pur.cut = .75, rel.cut = .75, ncpus = 2)
        
        TITAN.inci.snd <- titan(txa = matrify.incidence.titan,
                                env = sample.covariates$shrub.nat.dens,
                                pur.cut = .75, rel.cut = .75, ncpus = 2)
        
        TITAN.inci.sesr <- titan(txa = matrify.incidence.titan,
                                   env = sample.covariates$shrub.nat.esr,
                                   pur.cut = .75, rel.cut = .75, ncpus = 2)
        
        TITAN.inci.tnd <- titan(txa = matrify.incidence.titan,
                                env = sample.covariates$tree.nat.dens,
                                pur.cut = .75, rel.cut = .75, ncpus = 2)
        
        TITAN.inci.tod <- titan(txa = matrify.incidence.titan,
                                env = sample.covariates$tree.orn.dens,
                                pur.cut = .75, rel.cut = .75, ncpus = 2)
        
    
        plot_taxa(TITAN.inci.clouds, z2 = FALSE)
        plot_taxa(TITAN.inci.tcnd)
        plot_taxa(TITAN.inci.sesr, z1 = FALSE)
        plot_taxa(TITAN.inci.deadw)
        plot_taxa(TITAN.inci.height, z1 = FALSE)
        plot_taxa(TITAN.inci.imp)
        plot_taxa(TITAN.inci.snd)
        plot_taxa(TITAN.inci.tnd)
        plot_taxa(TITAN.inci.tod, z2 = FALSE)
        
        
## --- PLOT setup ------------------------------------------
        
        
        
        #sig.uPERM.results$var.names
        
        theme_set(theme_sjplot())
        
        # clouds only
        fit <-  lm(data = univariate.sprich, bird.effective.sp_2 ~ clouds.76.100.proportion)
        
        p.clouds <- plot_model(fit, type = "pred", terms = "clouds.76.100.proportion",
                               show.data = TRUE)
        
        # Impervious surface
        fit <-  lm(data = univariate.sprich, bird.effective.sp_2 ~ clouds.76.100.proportion + impervious.pct)
        
        p.impervious.pct <- plot_model(fit, type = "pred", terms = "impervious.pct", show.data = TRUE)

        
        # native shrub esr
        fit <-  lm(data = univariate.sprich, bird.effective.sp_2 ~ clouds.76.100.proportion + shrub.nat.esr)
        
        p.shrub.nat.esr <- plot_model(fit, type = "pred", terms = "shrub.nat.esr", show.data = TRUE)
        
        # native shrub density
        fit <-  lm(data = univariate.sprich, bird.effective.sp_2 ~ clouds.76.100.proportion + shrub.nat.dens)
        
        p.shrub.nat.dens <- plot_model(fit, type = "pred", terms = "shrub.nat.dens",
                                       show.data = TRUE)
        
        
        # tree cluster name
        fit <-  lm(data = univariate.sprich, bird.effective.sp_2 ~ clouds.76.100.proportion + tree.cluster.name)
        
        p.tree.cluster.name <- plot_model(fit, type = "pred", terms = "tree.cluster.name", show.data = TRUE)
        
        # shrub cluster name
        fit <-  lm(data = univariate.sprich, bird.effective.sp_2 ~ clouds.76.100.proportion + shrub.cluster.name)
        
        p.shrub.cluster.name <- plot_model(fit, type = "pred", terms = "shrub.cluster.name", show.data = TRUE)
        
        remove(fit)

## -- Save everything... ----------------------------------------------------------

save.image(file = paste0("bird_data_analysis_complete_", date, ".RData"))

