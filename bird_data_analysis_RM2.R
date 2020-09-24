## Bird analysis to support '03-birdanalysis.Rmd'.

## Data manipulation --> analysis (this doc) --> Rmarkdown/bookdown

## This analysis code is based on prior documents
## including 'bird_data_analysisX.R' and
## 'Office_Development_Bird_CommunitiesX.Rmd'.

## Updated for repeated measures for PERMANOVA 7/14
## Updated for colinearity 9/21 & fixing remaining units issues (e.g. acre --> hectare)


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
library(lme4)
library(MuMIn)


## Source data + helper R script

    source("bird_data_processing_RM.R")
    source('../../../../RCode/R_Scripts/AICc_PERMANOVA.R')
    source('../../../../RCode/R_Scripts/repeat_multipatt.R')
    source("../../../../RCode/R_Scripts/group_PERMANOVA.R")
    source("../../../../RCode/R_Scripts/group_lmer.R")
    source("../../../../RCode/R_Scripts/AICc_table_generation.R")

matrify.incidence.RM_0.clean <- matrify.incidence.RM[ , !(colnames(matrify.incidence.RM) == "no.birds")]
matrify.incidence_0.clean <- matrify.incidence[ , !(colnames(matrify.incidence) == "no.birds")]





## -- Descriptive Statistics  ---------------------------------------------

    # _2 suffix indicates that bird had to be seen at least twice
    # _0 is for the full data (incl only seen once)
    # clean indicates no 'no birds'

    descriptive.bird.bysite <-  tibble(
        site = rownames(matrify.incidence),
        bird.sp.richness = apply(matrify.incidence_0.clean,
                                 MARGIN = 1, specnumber),
        bird.sp.richness_2 = apply(matrify.incidence_2.clean,
                                   MARGIN = 1, specnumber),
        bird.sp.richness.2014_2 = apply(matrify.incidence.RM_2.clean[matrify.incidence.RM_2.clean$year==2014,1:31],
                                        MARGIN = 1, specnumber),
        bird.sp.richness.2015_2 = apply(matrify.incidence.RM_2.clean[matrify.incidence.RM_2.clean$year==2015,1:31],
                                        MARGIN = 1, specnumber),
        bird.Shannon = apply(matrify.incidence_0.clean,
                             MARGIN = 1, diversity),
        bird.Shannon_2 = apply(matrify.incidence_2.clean,
                               MARGIN = 1, diversity),
        bird.sp.Shannon.2014_2 = apply(matrify.incidence.RM_2.clean[matrify.incidence.RM_2.clean$year==2014,1:31],
                                        MARGIN = 1, diversity),
        bird.sp.Shannon.2015_2 = apply(matrify.incidence.RM_2.clean[matrify.incidence.RM_2.clean$year==2015,1:31],
                                        MARGIN = 1, diversity),
        bird.effective.sp = exp(bird.Shannon),
        bird.effective.sp_2 = exp(bird.Shannon_2),
        bird.effective.sp.2014_2 = exp(bird.sp.Shannon.2014_2),
        bird.effective.sp.2015_2 = exp(bird.sp.Shannon.2015_2),
        visits.no.birds = matrify.incidence$no.birds * 8,
        
        foraging.sp.richness = apply(matrify.foraging,
                                     MARGIN = 1, specnumber),
        foraging.sp.richness_2 = apply(matrify.foraging_2,
                                       MARGIN = 1, specnumber),
        foraging.sp.richness.2014_2 = apply(matrify.foraging.RM_2[matrify.foraging.RM_2$z.year==2014, 1:29],
                                        MARGIN = 1, specnumber),
        foraging.sp.richness.2015_2 = apply(matrify.foraging.RM_2[matrify.foraging.RM_2$z.year==2015,1:29],
                                        MARGIN = 1, specnumber),
        foraging.Shannon = apply(matrify.foraging,
                                 MARGIN = 1, diversity),
        foraging.Shannon_2 = apply(matrify.foraging_2,
                                   MARGIN = 1, diversity),
        foraging.sp.Shannon.2014_2 = apply(matrify.foraging.RM_2[matrify.foraging.RM_2$z.year==2014,1:29],
                                       MARGIN = 1, diversity),
        foraging.sp.Shannon.2015_2 = apply(matrify.foraging.RM_2[matrify.foraging.RM_2$z.year==2015,1:29],
                                       MARGIN = 1, diversity),
        
        foraging.effective.sp = exp(foraging.Shannon),
        foraging.effective.sp_2 = exp(foraging.Shannon_2),
        foraging.effective.sp.2014_2 = exp(foraging.sp.Shannon.2014_2),
        foraging.effective.sp.2015_2 = exp(foraging.sp.Shannon.2015_2),
        
        

    )

    temp <-
        summarise_all(descriptive.bird.bysite[, 2:26], tibble::lst(min, max, mean, sd, median)) 

        pretty.descriptive.birds.site <- tibble(
            Metric = colnames(descriptive.bird.bysite)[2:26],
            Minimum = unlist(temp[1:25]) %>% round(2),
            Maximum = unlist(temp[26:50]) %>% round(2),
            Mean = unlist(temp[51:75]) %>% round(2),
            `Standard Deviation` = unlist(temp[76:100]) %>% round(2),
            Median = unlist(temp[101:125]) %>% round(2)
        ) 

        remove(temp)

        pretty.descriptive.birds.site$Metric <- c(#"Site Area (hectares)",
                                                  "All yr Species Richness",
                                                  "All yr Sp. Richness (>1 sighting)",
                                                  "2014-2015 Sp. Richness (>1 sighting)",
                                                  "2015-2016 Sp. Richness (>1 sighting)",
                                                  "All yr Shannon Entropy",
                                                  "All yr Shannon Entropy (>1 sighting)",
                                                  "2014-2015 Shannon Entropy (>1 sighting)",
                                                  "2015-2016 Shannon Entropy (>1 sighting)",
                                                  "All yr Effective # Sp.",
                                                  "All yr Effective # Sp. (>1 sighting)",
                                                  "2014-2015 Effective # Sp. (>1 sighting)",
                                                  "2015-2016 Effective # Sp. (>1 sighting)",
                                                  
                                                  " No. Visits with No Birds in >0 Sweeps",
                                                  
                                                  "All yr Foraging Bird Sp. Richness",
                                                  "All yr Foraging Bird Sp. Richness (>1 sighting)",
                                                  "2014-2015 Foraging Bird Sp. Richness (>1 sighting)",
                                                  "2015-2016 Foraging Bird Sp. Richness (>1 sighting)",
                                                  "All yr Foraging Shannon Entropy",
                                                  "All yr Foraging Shannon (>1 sighting)",
                                                  "2014-2015 Foraging Shannon (>1 sighting)",
                                                  "2015-2016 Foraging Shannon (>1 sighting)",
                                                  "All yr Foraging Effective # Sp.",
                                                  "All yr Foraging Effective # Sp. (>1 sighting)",
                                                  "2014-2015 Foraging Effective # Sp. (>1 sighting)",
                                                  "2015-2016 Foraging Effective # Sp. (>1 sighting)"
                                                  
                                                  )




        descriptive.bird.bysp <- tibble(
            species.name = colnames(matrify.incidence),
            num.sites = apply(matrify.incidence, 2, specnumber),
            incidence.min = apply(matrify.incidence, 2, min),
            incidence.max = apply(matrify.incidence, 2, max),
            incidence.median = apply(matrify.incidence, 2, median))
        descriptive.bird.bysp.2014 <- tibble(
            species.name = colnames(matrify.incidence.2014[1:37]),
            num.sites = apply(matrify.incidence.2014[1:37], 2, specnumber),
            incidence.min = apply(matrify.incidence.2014[1:37], 2, min),
            incidence.max = apply(matrify.incidence.2014[1:37], 2, max),
            incidence.median = apply(matrify.incidence.2014[1:37], 2, median))
        
        descriptive.bird.bysp.2015 <- tibble(
            species.name = colnames(matrify.incidence.2015[1:37]),
            num.sites = apply(matrify.incidence.2015[1:37], 2, specnumber),
            incidence.min = apply(matrify.incidence.2015[1:37], 2, min),
            incidence.max = apply(matrify.incidence.2015[1:37], 2, max),
            incidence.median = apply(matrify.incidence.2015[1:37], 2, median))

        temp <- tibble(
            species.name = colnames(matrify.foraging),
            foraging.num.sites = apply(matrify.foraging, 2, specnumber),
            foraging.min = apply(matrify.foraging, 2, min),
            foraging.max = apply(matrify.foraging, 2, max),
            foraging.median = apply(matrify.foraging, 2, median)
        )
        temp.2014 <- tibble(
            species.name = colnames(matrify.foraging.2014[1:31]),
            foraging.num.sites = apply(matrify.foraging.2014[1:31], 2, specnumber),
            foraging.min = apply(matrify.foraging.2014[1:31], 2, min),
            foraging.max = apply(matrify.foraging.2014[1:31], 2, max),
            foraging.median = apply(matrify.foraging.2014[1:31], 2, median)
        )
        temp.2015 <- tibble(
            species.name = colnames(matrify.foraging.2015[1:31]),
            foraging.num.sites = apply(matrify.foraging.2015[1:31], 2, specnumber),
            foraging.min = apply(matrify.foraging.2015[1:31], 2, min),
            foraging.max = apply(matrify.foraging.2015[1:31], 2, max),
            foraging.median = apply(matrify.foraging.2015[1:31], 2, median)
        )

        descriptive.bird.bysp <- merge(descriptive.bird.bysp, temp, by = "species.name", all.x = TRUE)
        descriptive.bird.bysp[is.na(descriptive.bird.bysp)] <- 0
        
        descriptive.bird.bysp.2014 <- merge(descriptive.bird.bysp.2014, temp.2014, by = "species.name", all.x = TRUE)
        descriptive.bird.bysp.2014[is.na(descriptive.bird.bysp.2014)] <- 0
        
        descriptive.bird.bysp.2015 <- merge(descriptive.bird.bysp.2015, temp.2015, by = "species.name", all.x = TRUE)
        descriptive.bird.bysp.2015[is.na(descriptive.bird.bysp.2015)] <- 0

        remove(temp, temp.2014, temp.2015)
        descriptive.bird.bysp$incidence.range <- paste0(descriptive.bird.bysp$incidence.min,
                                                        "-", descriptive.bird.bysp$incidence.max)
        descriptive.bird.bysp$foraging.range <- paste0(descriptive.bird.bysp$foraging.min,
                                                        "-", descriptive.bird.bysp$foraging.max)
        descriptive.bird.bysp.2014$incidence.range <- paste0(descriptive.bird.bysp.2014$incidence.min,
                                                        "-", descriptive.bird.bysp.2014$incidence.max)
        descriptive.bird.bysp.2014$foraging.range <- paste0(descriptive.bird.bysp.2014$foraging.min,
                                                       "-", descriptive.bird.bysp.2014$foraging.max)
        descriptive.bird.bysp.2015$incidence.range <- paste0(descriptive.bird.bysp.2015$incidence.min,
                                                        "-", descriptive.bird.bysp.2015$incidence.max)
        descriptive.bird.bysp.2015$foraging.range <- paste0(descriptive.bird.bysp.2015$foraging.min,
                                                       "-", descriptive.bird.bysp.2015$foraging.max)
        
        
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

        
        
        
        
        
        
        
        bird.species.names <-
          tibble(common.name = unique(incidence$text.spp))
        bird.species.names$common.name <-
          bird.species.names$common.name[order(bird.species.names$common.name)]
        
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

        pretty.descriptive.birds.species <-
          merge(
            descriptive.bird.bysp,
            bird.species.names,
            by.x = "species.name",
            by.y = "common.name",
            all.x = TRUE,
            sort = FALSE
          ) %>%
          merge(
            pretty.guild,
            by.x = "species.name",
            by.y = "Species Name",
            all.x = TRUE,
            sort = FALSE
          )
        pretty.descriptive.birds.species <-
          pretty.descriptive.birds.species[, c("pretty.common.name",
                                               "Scientific Name",
                                               "Diet",
                                               "Foraging Substrate",
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
    
    remove(matrify.foraging.2014, matrify.foraging.2014_2,
           matrify.foraging.2015, matrify.foraging.2015_2,
           guild.eff.sp.rich.2014_2.clean, guild.eff.sp.rich.2015_2.clean, 
           matrify.incidence.2014, matrify.incidence.2014_2.clean,
           matrify.incidence.2015, matrify.incidence.2015_2.clean)

    descriptive.bird.bysite_RM <- tibble(
      site = rep(descriptive.bird.bysite$site, 2),
      year = c(rep(2014, 20), rep(2015, 20)),
      bird.effective.sp_RM = c(
        descriptive.bird.bysite$bird.effective.sp.2014_2,
        descriptive.bird.bysite$bird.effective.sp.2015_2
      )
    )
    
  ## Set up the tables for univariate comparison
    
    descriptive.bird.bysite_RM <- merge(
      descriptive.bird.bysite_RM,
      sample.covariates,
      by.x = c("site", "year"),
      by.y = c("SiteName", "year"),
      sort = F
    )
    
    guild.eff.sp.rich.RM_2.clean <-
      merge(
        guild.eff.sp.rich.RM_2.clean,
        sample.covariates,
        by.x = c("site.name", "year"),
        by.y = c("SiteName", "year"),
        sort = F
      )
    


    
    
    
## -- Testing GLMM implementation -----------------------------
    lme4::lmer(
      descriptive.bird.bysite_RM$bird.effective.sp_RM ~
        descriptive.bird.bysite_RM$year +
        (1 | descriptive.bird.bysite_RM$site)
    ) 

    
    summary(aov(
      bird.effective.sp_RM ~ clouds.76.100.proportion + stand.predate.development + Error(site),
      data = descriptive.bird.bysite_RM
    ))
    
    fit <- lmer(bird.effective.sp_RM ~ clouds.76.100.proportion + height.m.median + (1|site),
                data=descriptive.bird.bysite_RM)
    summary(fit)
    car::Anova(fit, test = "F")
    
    
    


## -- Regression Variables ----------------------------------
    
    detection.vars <-
        c(
            "proportion.loud.noise",
            "clouds.76.100.proportion",
            "Solar.Rad.W_m2.8v.median"
        )
    

    neighborhood.vars <- c("dissolved_parcel_500m_buffer_impervious_500m_mean",
                           "dissolved_parcel_500m_intersections_NUMPOINTS",
                           "Proportion_ForeignBorn_B99051e5",
                           "MedianHouseholdIncome_B19013e1",
                           "hectares",
                           "age.2017",
                           "TV_mean",
                           "SMV_mean",
                           "landrentperhect")
    
    management.vars <- c("impervious.pct",
                         "dead.wood",
                         "height.m.median",
                         "stand.predate.development",
                         
                         "tree.orn.dens",
                         "tree.orn.esr",
                         "conifer.nat.dens",
                         
                         "shrub.nat.dens",
                         "shrub.nat.esr",
                         "shrub.orn.dens",
                         "shrub.orn.esr"
                         )

## -- Outliers & Colin ----------------------------------------------------------
    
    
    # # Outliers
    #     outlier1 <- par(mfrow = c(4,3), mar = c(2,2,2,1))
    #     dotchart(univariate.sprich$bird.effective.sp_2, main = "Effective Sp. Richness (>2)")
    #     dotchart(univariate.sprich$hectares, main = "hectares")
    #     dotchart(univariate.sprich$dissolved_parcel_500m_buffer_impervious_500m_mean,
    #              main = "500m Impervious Surface (%)")
    #     dotchart(univariate.sprich$Proportion_ForeignBorn_B99051e5, main = "% Foreign Born")
    #     dotchart(univariate.sprich$MedianHouseholdIncome_B19013e1, main = "Median Household Income ($)")
    #     dotchart(univariate.sprich$age.2017, main = "Site Age (baseline 2017)")
    #     dotchart(univariate.sprich$height.m.median, main = "Median DF Height (m)")
    #     dotchart(univariate.sprich$dead.wood, main = "All Dead Wood")
    #     dotchart(univariate.sprich$conifer.nat.abundance, main = "Native Conifer Density")
    #     dotchart(univariate.sprich$shrub.nat.dens, main = "Native Shrub Density")
    #     dotchart(univariate.sprich$impervious.sqft, main = "Site Impervious Surface (%)")
    #     (dev.off())
    # 
    #     # There are potential outliers--particularly Median DF Height (2 sites with
    #     # no DF) and all dead wood (two sites w/high values).

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
      r <- abs(cor(x, y, method = "k"))
      txt <- format(c(r, 0.123456789), digits = digits)[1]
      txt <- paste0(prefix, txt)
      if (missing(cex.cor)) cex.cor <- .9 + .5/strwidth(txt)
      text(0.5, 0.5, txt, cex = cex.cor + r)
    }
    
    cor.temp <- select_if(descriptive.bird.bysite_RM[,
                                                     colnames(descriptive.bird.bysite_RM) %in%
                                                       c(detection.vars, neighborhood.vars, management.vars)],
                          is.numeric)
    
    colnames(cor.temp) <-
      c(
        "Impervious 500m",
        "Intersect 500m",
        "Foreign Born",
        "Median HH Income",
        "Tall Veg 500m",
        "S/M Veg 500m",
        "Building Age",
        "Site Imp.",
        "Med Conif Ht",
        "Dead Wood",
        "Non-nat Tree Dens",
        "Non-nat Tree ESR",
        "Nat Conif Dens",
        "Nat Shrb Dens",
        "Nat Shrb ESR",
        "Non-nat Shrb Dens",
        "Non-nat Shrb ESR",
        
        "Site Area",
        "Value per Hect.",

        "Loud Noise",
        "Overcast",
        "Solar Rad."
      )
    
    pairs(
      cor.temp,
      panel = panel.smooth,
      upper.panel = panel.cor
    )
    
    # want histograms to be uniformly distributed (across a range of values),
    # correlations between explanatory variables low. Most are fine except for
    # impervious surface on site, which is highly correlated with dead wood and
    # native conifer density. Expect that median DF height is going to be
    # important.
    
    
    
## -- Eff. Sp. Richness GLMM ------------------------------------------------------------
    
    
        
## Site Detection Variables


    detection.uLMER.results <-
        group.lmer(
            var.names = detection.vars,
            in.data = descriptive.bird.bysite_RM,
            random.subject.effect = "(1 | site)",
            dependent.var = "bird.effective.sp_RM"
        )

    ## Nothing is significant after, clouds before adjustment.

    
## Neighborhood Variables & Management and Landscape Variables


    neigh.mngmt.uLMER.results <-
        group.lmer(
            var.names = c(neighborhood.vars,management.vars),
            in.data = descriptive.bird.bysite_RM,
            random.subject.effect = "(1 | site)",
            dependent.var = "bird.effective.sp_RM"
        )

    ## most are median height and stands predating.

    
    ## S3 method for class 'AIClme'
    # modavgEffect(cand.set, modnames = NULL, newdata,
    #              second.ord = TRUE, nobs = NULL, uncond.se = "revised",
    #              conf.level = 0.95, ...) 
    # 
    
    
    
    
    
## Create a pretty table:
    
    
    pretty.all.uLMER.results <- rbind(
        detection.uLMER.results,
        neigh.mngmt.uLMER.results
    )
    
    colnames(pretty.all.uLMER.results) <-
        c("Variable", "Slope Estimate", "Std Error", "Mean Sq",
          "F statistic", "_p_-value", "Adjusted _p_-val",
          "AICc", "AICc REML", "log Likelihood")

    pretty.all.uLMER.results$`Delta AICc` <-
      pretty.all.uLMER.results$`AICc REML` - min(pretty.all.uLMER.results$`AICc REML`)
    
    pretty.vars <- tibble(
        coded.name = pretty.all.uLMER.results$Variable,
        `Variable Name` = c(
          "Loud Noise Visits (%)",
          "Overcast Visits (%; Control Var)",
          "Median Solar Radiation (W/m^2)",

          "Mean Impervious within 500m (%)",
          "Major Intersections within 500m",
          "Foreign Born (%)",
          "Median Household Income (USD 1'000s)",
          "Site Area (hectares)",
          "Building Age (years)",
          "Tall Vegetation within 500m (% of area)",
          "Short and Medium Vegetation within 500m (% of area)",
          "Assessed Value per Hectare (USD 1'000s)",
            
          "On Site Impervious (%)",
          "Dead Wood Abundance",
          "Median Dominant Conifer Height (m)",
          "Stand Predates Development",
          "Non-native Tree Density",
          "Non-native Tree Effective Richness",
          "Native Conifer Density",
          "Native Shrub Density",
          "Native Shrub Effective Richness",
          "Non-native Shrub Density",
          "Non-native Shrub Effective Richness"
        )
    )
    
    pretty.all.uLMER.results <- merge(pretty.all.uLMER.results,
                                      pretty.vars,
                                      by.x = "Variable",
                                      by.y = "coded.name", sort = F)

    
    sig.uLMER.results <-
      pretty.all.uLMER.results[pretty.all.uLMER.results$`_p_-value` < 0.05, ]
    
    



## -------------------- Guild tests ---------------------------------------------
    ## Look at different guild lm for relationships with significant variables (to report)

    guild.eff.sp.rich.RM_2.clean$stand.predate.developmentRL <- 
        relevel(as.factor(guild.eff.sp.rich.RM_2.clean$stand.predate.development), ref = "no.unkn")
    
    
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
        guild.eff.sp.rich.RM_2.clean[, guild.cols.totest.c]

    guild.tests <- function(guild.var.name, guild.var.vector, control.var = NULL, table.name = "guild.eff.sp.rich.RM_2.clean") {
      temp <- data.frame(
        row.names = guild.cols.totest.c,
        guild.names = guild.cols.totest.c,
        corr.SP = rep(NA, times = length(guild.cols.totest.c)),
        slope.estimate = rep(NA),
        std.error = rep(NA),
        Mean.Sq = rep(NA),
        F.stat = rep(NA),
        prob.pval = rep(NA),
        adj.pval = rep(NA),
        AICc = rep(NA),
        AIC.REML = rep(NA),
        logLik = rep(NA)
        
      )
        i = 1
        for (i in 1:length(guild.cols.totest.c)) {
            
            LMER.temp <-
                group.lmer(
                    var.names = guild.var.name,
                    in.data = guild.eff.sp.rich.RM_2.clean,
                    dependent.var = guild.cols.totest.v[i],
                    control.var = control.var,
                    random.subject.effect = "(1 | site.name)"
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
            else {temp$corr.SP[i] <- "Not Numeric"
            }
            
            
            temp$slope.estimate[i] <- LMER.temp$slope.estimate
            temp$std.error[i] <- LMER.temp$std.error
            
            temp$Mean.Sq[i] <- LMER.temp$Mean.Sq
            
            temp$F.stat[i] <- LMER.temp$F.stat
            temp$prob.pval[i] <- LMER.temp$prob.pval
            
            
            temp$AICc[i] <- LMER.temp$AICc
            temp$AIC.REML[i] <- LMER.temp$AIC.REML
            temp$logLik[i] <- LMER.temp$logLik
            
        }
        temp$adj.pval <- p.adjust(p = temp$prob.pval, method = "holm")
        
        temp$sig.stars <- ifelse(temp$adj.pval <0.001, "***",
                                 ifelse(temp$adj.pval < 0.01 & temp$adj.pval > 0.001, "**",
                                        ifelse(temp$adj.pval < 0.05 & temp$adj.pval > 0.01, "*", "")))

        if(is.numeric(guild.var.vector)){
          temp$pretty <- ifelse(
            temp$adj.pval < 0.05,
            yes = paste0(
              round(temp$slope.estimate, 3),
              " (SE: ",
              round(temp$std.error, 2),
              "; F stat:",
              round(temp$F.stat, 2),
              ")",
              temp$sig.stars
            ),
            no = "NS"
          )
        }
        else{
            temp$pretty <- ifelse(
                temp$adj.pval < 0.05,
                yes = paste0(
                    round(temp$slope.estimate,3),
                    " (SE: ",
                    round(temp$std.error, 2), "; F stat:",
                    round(temp$F.stat, 2),
                    ")",
                    temp$sig.stars
                ),
                no = "NS"
            )
        }
        
        
        return(temp)
    }
    
    
    

    
## Guilds with median height of dominant trees
    
    guild.height.results <-
        guild.tests(guild.var.name = "height.m.median",
                    guild.var.vector = guild.eff.sp.rich.RM_2.clean$height.m.median,
            )

# Stands predating development

    guild.spd.results <- guild.tests("stand.predate.developmentRL",
                                     guild.eff.sp.rich.RM_2.clean$stand.predate.developmentRL)

  # Less significant vars
    
    guild.tod.results <- guild.tests("tree.orn.dens",
                                     guild.eff.sp.rich.RM_2.clean$tree.orn.dens)
    
    guild.imperv.results <- guild.tests("impervious.pct",
                                        guild.eff.sp.rich.RM_2.clean$impervious.pct)

    guild.cnd.results <- guild.tests("conifer.nat.dens",
                                     guild.eff.sp.rich.RM_2.clean$conifer.nat.dens)
    guild.snd.results <- guild.tests("shrub.nat.dens",
                                     guild.eff.sp.rich.RM_2.clean$shrub.nat.dens)
    guild.snesr.results <- guild.tests("shrub.nat.esr",
                                       guild.eff.sp.rich.RM_2.clean$shrub.nat.esr)
    guild.smv.results <- guild.tests("SMV_mean",
                                     guild.eff.sp.rich.RM_2.clean$SMV_mean)
    guild.tv.results <- guild.tests("TV_mean",
                                     guild.eff.sp.rich.RM_2.clean$TV_mean)
    
    #Other veg variables
    guild.toesr.results <- guild.tests("tree.orn.esr",
                                    guild.eff.sp.rich.RM_2.clean$tree.orn.esr)
    
    guild.sod.results <- guild.tests("shrub.orn.dens",
                                    guild.eff.sp.rich.RM_2.clean$shrub.orn.dens)
    
    guild.soesr.results <- guild.tests("shrub.orn.esr",
                                     guild.eff.sp.rich.RM_2.clean$shrub.orn.esr)

    
  # Non-sig variables
    guild.yba.results <- guild.tests("YrBuiltAvg",
                                     guild.eff.sp.rich.RM_2.clean$YrBuiltAvg)
    guild.dw.results <- guild.tests("dead.wood",
                                     guild.eff.sp.rich.RM_2.clean$dead.wood)
    guild.lpa.results <- guild.tests("landrentperhect",
                                     guild.eff.sp.rich.RM_2.clean$landrentperhect)
    guild.ac.results <- guild.tests("hectares",
                                     guild.eff.sp.rich.RM_2.clean$hectares)
    guild.bimp.results <- guild.tests("dissolved_parcel_500m_buffer_impervious_500m_mean",
                                     guild.eff.sp.rich.RM_2.clean$dissolved_parcel_500m_buffer_impervious_500m_mean)
    guild.int.results <- guild.tests("dissolved_parcel_500m_intersections_NUMPOINTS",
                                      guild.eff.sp.rich.RM_2.clean$dissolved_parcel_500m_intersections_NUMPOINTS)
    guild.mhhi.results <- guild.tests("MedianHouseholdIncome_B19013e1",
                                     guild.eff.sp.rich.RM_2.clean$MedianHouseholdIncome_B19013e1)
    guild.pfb.results <- guild.tests("Proportion_ForeignBorn_B99051e5",
                                     guild.eff.sp.rich.RM_2.clean$Proportion_ForeignBorn_B99051e5)
    
    # Detection variables
    guild.lnv.results <- guild.tests("proportion.loud.noise",
                                     guild.eff.sp.rich.RM_2.clean$proportion.loud.noise)
    guild.cloud.results <- guild.tests("clouds.76.100.proportion",
                                     guild.eff.sp.rich.RM_2.clean$clouds.76.100.proportion)
    guild.sun.results <- guild.tests("Solar.Rad.W_m2.8v.median",
                                    guild.eff.sp.rich.RM_2.clean$Solar.Rad.W_m2.8v.median)

    
plot(guild.eff.sp.rich.RM_2.clean$insectivore.esr, guild.eff.sp.rich.RM_2.clean$tree.dens)

# note a lot of tree density comes from douglas fir. so this isn't surprising.

plot(guild.eff.sp.rich.RM_2.clean$conifer.esr, guild.eff.sp.rich.RM_2.clean$tree.orn.dens)


## ------------ Multivariate Regression Prep ----------------------------

plot(betadisper(vegdist(matrify.incidence.RM_0.clean[1:34]), 
           group = matrify.incidence.RM_0.clean$site))
plot(betadisper(vegdist(matrify.incidence.RM_0.clean[1:34]), 
                group = matrify.incidence.RM_0.clean$year))
    # sites differ in space (and dispersion some for VV and I1), but year does
    # not differ in space or dispersion.

matrify.incidence.centroids <- 
    betadisper(vegdist(matrify.incidence.RM_0.clean[1:34]),
           group = matrify.incidence.RM_0.clean$site,
           type = 'centroid',
           add = 'lingoes')$centroids

matrify.incidence.centroids <- as.data.frame(matrify.incidence.centroids,
                                             row.names =  matrify.incidence.RM_0.clean$site[1:20])


# Foraging birds only:

plot(betadisper(vegdist(matrify.foraging.RM[1:31]), 
                group = matrify.foraging.RM$z.site))
     
plot(betadisper(vegdist(matrify.foraging.RM[1:31]), 
                group = matrify.foraging.RM$z.year))

# same dispersion comments as full bird dataset

matrify.foraging.centroids <- 
    betadisper(vegdist(matrify.foraging.RM[1:31]),
               group = matrify.foraging.RM$z.site,
               type = 'centroid',
               add = 'lingoes')$centroids

matrify.foraging.centroids <-
    as.data.frame(matrify.foraging.centroids,
                  row.names =  matrify.foraging.RM$z.site[1:20])




##  -- Incidence PERMANOVA ---------------------------------------------------------

    ## Detection Variables


        detection.mPERM.results <-
            group.PERMANOVA(
                var.names = detection.vars,
                var.table =  sample.covariates.20,
                var.table.c = "sample.covariates.20",
                species.table = matrify.incidence.centroids,
                species.table.c = "matrify.incidence.centroids",
                method = "euclidean",
                by.adonis2 = "terms"
            )
        ## clouds significant before adjustments


    ## Socio-economic neighborhood variables...


        neighborhood.mPERM.results <-
            group.PERMANOVA(
                var.names = neighborhood.vars,
                var.table =  sample.covariates.20,
                var.table.c = "sample.covariates.20",
                species.table = matrify.incidence.centroids,
                species.table.c = "matrify.incidence.centroids",
                num.control.vars = 1,
                control.vars = "clouds.76.100.proportion",
                method = "euc",
                by.adonis2 = "terms"
            )

        ## none significant before adjustment


    ## Landscaping and maintenance variables


        mngmt.mPERM.results <-
            group.PERMANOVA(
                var.names = management.vars,
                var.table =  sample.covariates.20,
                var.table.c = "sample.covariates.20",
                species.table = matrify.incidence.centroids,
                species.table.c = "matrify.incidence.centroids",
                num.control.vars = 1,
                method = "euc",
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
            mngmt.mPERM.results
            
        )[, c(1, 3, 4, 5, 7, 10)]

        pretty.all.mPERM.results$delta.aicc <-
            pretty.all.mPERM.results$AIC.stat - min(pretty.all.mPERM.results$AIC.stat)
        
        colnames(pretty.all.mPERM.results) <-
            c("Variable", "Variation Explained per df", "pseudo-_F_", "_p_-value", "Dispersion _p_-value", "AICc", "Delta AICc")
        
        pretty.all.mPERM.results <- merge(pretty.all.mPERM.results,
                                          pretty.vars,
                                          by.x = "Variable",
                                          by.y = "coded.name")
        


## AICc comparison:

    # Overall, nine possibly significant variables: native conifer density, dead
    # wood, median height, impervious on site, native shrub density,
    # native shrub esr, stands predating development, native tree
    # density, tree community




## --- Foraging PERMANOVA ------------------------------------

## Foraging incidence data

    ## Detection Variables

        detection.mPERM.foraging <-
            group.PERMANOVA(
                var.names = detection.vars,
                var.table =  sample.covariates.20,
                var.table.c = "sample.covariates.20",
                species.table = matrify.foraging.centroids,
                species.table.c = "matrify.foraging.centroids",
                method = "euc",
                num.control.vars = 0,
                by.adonis2 = "terms"
            )
        ## none significant after adjustment (clouds significant before adjusting tho)



    ## Socio-economic neighborhood variables...

        neighborhood.mPERM.foraging <-
            group.PERMANOVA(
                var.names = neighborhood.vars,
                var.table =  sample.covariates.20,
                var.table.c = "sample.covariates.20",
                species.table = matrify.foraging.centroids,
                species.table.c = "matrify.foraging.centroids",
                num.control.vars = 1,
                control.vars = "clouds.76.100.proportion",
                by.adonis2 = "terms",
                method = "euc"
            )
        ## none significant before or after adjustment. SMV is the only one that is close (.07)


    ## Landscaping and maintenance variables

        mngmt.mPERM.foraging <-
            group.PERMANOVA(
                var.names = management.vars,
                var.table =  sample.covariates.20,
                var.table.c = "sample.covariates.20",
                species.table = matrify.foraging.centroids,
                species.table.c = "matrify.foraging.centroids",
                num.control.vars = 1,
                control.vars = "clouds.76.100.proportion",
                by.adonis2 = "terms",
                method = "euc"
            )

        ## native conifer density, dead wood abundance, median dom conif height,
        ## impervious, native shurb density, native shrub esr, stands predating
        ## development, native tree density. Native tree esr borderline.

        ## Same as all birds. 

        
        
        pretty.all.mPERM.foraging <- rbind(
            detection.mPERM.foraging,
            neighborhood.mPERM.foraging,
            mngmt.mPERM.foraging
        )[, c(1, 3, 4, 5, 7, 10)]
        
        pretty.all.mPERM.foraging$delta.AICc <- 
            pretty.all.mPERM.foraging$AIC.stat - min(pretty.all.mPERM.foraging$AIC.stat)
        
        colnames(pretty.all.mPERM.foraging) <-
          c("Variable", "Variation Explained per df", "pseudo-_F_", "_p_-value", "Dispersion _p_-value", "AICc", "Delta AICc")
        
        pretty.all.mPERM.foraging <- merge(pretty.all.mPERM.foraging,
                                           pretty.vars,
                                           by.x = "Variable",
                                           by.y = "coded.name")
        
        



## -- NMDS Visualisations ------------------------------------------------------------------

        # All birds
        
        
        
        pcoa.incidence <-
          cmdscale(vegdist(matrify.incidence.RM_0.clean[, 1:36], method = "bray"))
        plot(pcoa.incidence)
        
        sig.mPERM.results <-
          pretty.all.mPERM.results[pretty.all.mPERM.results$`Delta AICc` < 2,]
        
        
        # Foraging birds
        
        
        
        pcoa.foraging <-
          cmdscale(vegdist(matrify.foraging.RM[, 1:31], method = "bray"))
        plot(pcoa.foraging)
        
        sig.mPERM.foraging <-
          pretty.all.mPERM.foraging[pretty.all.mPERM.foraging$`Delta AICc` < 2,]
        
        
        
        
    incidence.envifit <- envfit(pcoa.incidence, 
               sample.covariates[, colnames(sample.covariates) %in% 
                                        pretty.all.mPERM.results$Variable])

    
    plot(incidence.NMDS)
    plot(envfit(pcoa, 
                sample.covariates.20[, colnames(sample.covariates.20) %in% 
                                       sig.all.mPERM.results$Variable]))
        
        
## -- CORRELATION INDEX / INDICATOR SPECIES ---------------------------------------------
        
    # total incidence used for these
    
    # For stands of trees predating development:
        
        ## Point biserial correlation coeff.
        
        spd.cluster.pbiserial <-
            repeat.multipatt(
                matrix.name = matrify.incidence_0.clean,
                phi = FALSE,
                func.name = "r.g",
                cluster.name = sample.covariates.20$stand.predate.development,
                plot.please = FALSE,
                quiet = FALSE,
                freq.cutoff = .5,
                stat.cutoff = .5
            )
        
        spd.cluster.pbiserial$group <- ifelse(spd.cluster.pbiserial$groupname == "yes", "Yes, Stand Predates Development", "No")
        spd.cluster.pbiserial <-
            spd.cluster.pbiserial[order(spd.cluster.pbiserial$mean.stat, decreasing = TRUE), ]
        
        spd.cluster.pbiserial <-
            left_join(spd.cluster.pbiserial,
                      bird.guilds,
                      by = c("species" = "species.name"))
        spd.cluster.pbiserial <- left_join(spd.cluster.pbiserial, bird.species.names, by = c("species" = "common.name"))
        
        spd.cluster.pbiserial$pval.range <-
            paste0(spd.cluster.pbiserial$min.p.val,
                   "-",
                   spd.cluster.pbiserial$max.p.val)
        
        spd.cluster.pbiserial$foraging.substrate <-
            ifelse(
                spd.cluster.pbiserial$foraging.substrate == "trees.shrubs",
                
                "Trees & Shrubs",
                "Ground"
            )
        spd.cluster.pbiserial$forest.preference <-
            ifelse(
                spd.cluster.pbiserial$forest.preference == "conifer",
                "Conifer",
                "Open or No Preference"
            )
        spd.cluster.pbiserial$diet <-
            ifelse(
                spd.cluster.pbiserial$diet == "omnivore",
                "Omnivore",
                "Insectivore"
            )
        
        spd.foraging.pbiserial <-
            repeat.multipatt(
                matrix.name = matrify.foraging,
                phi = FALSE,
                func.name = "r.g",
                cluster.name = sample.covariates.20$stand.predate.development,
                plot.please = FALSE,
                quiet = FALSE,
                freq.cutoff = .5,
                stat.cutoff = .5
            )
        
        spd.foraging.pbiserial$group <- ifelse(spd.foraging.pbiserial$groupname == "yes", "Yes, 1+ Stands Predate Development", "No Stands Predate Development")
        spd.foraging.pbiserial <-
            spd.foraging.pbiserial[order(spd.foraging.pbiserial$mean.stat, decreasing = TRUE), ]
        
        
        spd.foraging.pbiserial <-
            left_join(spd.foraging.pbiserial,
                      bird.guilds,
                      by = c("species" = "species.name"))
        spd.foraging.pbiserial <- left_join(spd.foraging.pbiserial, bird.species.names, by = c("species" = "common.name"))
        
        spd.foraging.pbiserial$pval.range <-
            paste0(spd.foraging.pbiserial$min.p.val,
                   "-",
                   spd.foraging.pbiserial$max.p.val)
        
        spd.foraging.pbiserial$foraging.substrate <-
            ifelse(
                spd.foraging.pbiserial$foraging.substrate == "trees.shrubs",
                
                "Trees & Shrubs",
                "Ground"
            )
        spd.foraging.pbiserial$forest.preference <-
            ifelse(
                spd.foraging.pbiserial$forest.preference == "conifer",
                "Conifer",
                "Open or No Preference"
            )
        spd.foraging.pbiserial$diet <-
            ifelse(
                spd.foraging.pbiserial$diet == "omnivore",
                "Omnivore",
                "Insectivore"
            )
        
        
        
## -- TITAN Analysis ------------------------------------------------------------

        matrify.incidence.titan <-
            matrify.incidence[, colSums(matrify.incidence > 0) > 3]
        
        TITAN.inci.clouds <- titan(txa = matrify.incidence.titan,
                                   env = sample.covariates.20$clouds.76.100.proportion,
                                   pur.cut = .75, rel.cut = .75, ncpus = 2)

        TITAN.inci.cnd <- titan(txa = matrify.incidence.titan,
                                 env = sample.covariates.20$conifer.nat.dens,
                                 pur.cut = .75, rel.cut = .75, ncpus = 2)

        TITAN.inci.deadw <- titan(txa = matrify.incidence.titan,
                                  env = sample.covariates.20$dead.wood,
                                  pur.cut = .75, rel.cut = .75, ncpus = 2)
                
        TITAN.inci.height <- titan(txa = matrify.incidence.titan,
                                   env = sample.covariates.20$height.m.median,
                                   pur.cut = .75, rel.cut = .75, ncpus = 2)

        TITAN.inci.imp <- titan(txa = matrify.incidence.titan,
                                   env = sample.covariates.20$impervious.pct,
                                   pur.cut = .75, rel.cut = .75, ncpus = 2)
        
        TITAN.inci.sesr <- titan(txa = matrify.incidence.titan,
                                   env = sample.covariates.20$shrub.nat.esr,
                                   pur.cut = .75, rel.cut = .75, ncpus = 2)
        
        TITAN.inci.snd <- titan(txa = matrify.incidence.titan,
                                env = sample.covariates.20$shrub.nat.dens,
                                pur.cut = .75, rel.cut = .75, ncpus = 2)
        
        TITAN.inci.tod <- titan(txa = matrify.incidence.titan,
                                 env = sample.covariates.20$tree.orn.dens,
                                 pur.cut = .75, rel.cut = .75, ncpus = 2)
        
        # plot_taxa(TITAN.inci.clouds, z2 = FALSE)
        # plot_taxa(TITAN.inci.tcnd)
        # plot_taxa(TITAN.inci.sesr, z1 = FALSE)
        # plot_taxa(TITAN.inci.deadw)
        # plot_taxa(TITAN.inci.height, z1 = FALSE)
        # plot_taxa(TITAN.inci.imp)
        # plot_taxa(TITAN.inci.snd)
        # plot_taxa(TITAN.inci.tnd)
        # plot_taxa(TITAN.inci.tod, z2 = FALSE)

    matrify.foraging.titan <-
        matrify.foraging[, colSums(matrify.foraging > 0) > 3]
        
        TITAN.fora.clouds <- titan(txa = matrify.foraging.titan,
                                   env = sample.covariates.20$clouds.76.100.proportion,
                                   pur.cut = .75, rel.cut = .75, ncpus = 2)
        
        TITAN.fora.cnd <- titan(txa = matrify.foraging.titan,
                                 env = sample.covariates.20$conifer.nat.dens,
                                 pur.cut = .75, rel.cut = .75, ncpus = 2)
        
        TITAN.fora.deadw <- titan(txa = matrify.foraging.titan,
                                  env = sample.covariates.20$dead.wood,
                                  pur.cut = .75, rel.cut = .75, ncpus = 2)
        
        TITAN.fora.height <- titan(txa = matrify.foraging.titan,
                                   env = sample.covariates.20$height.m.median,
                                   pur.cut = .75, rel.cut = .75, ncpus = 2)
        
        TITAN.fora.imp <- titan(txa = matrify.foraging.titan,
                                env = sample.covariates.20$impervious.pct,
                                pur.cut = .75, rel.cut = .75, ncpus = 2)
        
        TITAN.fora.sesr <- titan(txa = matrify.foraging.titan,
                                 env = sample.covariates.20$shrub.nat.esr,
                                 pur.cut = .75, rel.cut = .75, ncpus = 2)
        
        TITAN.fora.snd <- titan(txa = matrify.foraging.titan,
                                env = sample.covariates.20$shrub.nat.dens,
                                pur.cut = .75, rel.cut = .75, ncpus = 2)
        
        TITAN.fora.tod <- titan(txa = matrify.foraging.titan,
                                env = sample.covariates.20$tree.orn.dens,
                                pur.cut = .75, rel.cut = .75, ncpus = 2)

        
## --- PLOT setup ------------------------------------------
        
        
        library(sjPlot)
        fit <- lmer(bird.effective.sp_RM ~ height.m.median + (1|site),
                    data=descriptive.bird.bysite_RM)
        summary(fit)
        car::Anova(fit, test = "F")
        
        sjPlot::plot_model(fit)
        
        



## -- Save everything... ----------------------------------------------------------

save.image(file = paste0("bird_data_analysis_complete_", date, ".RData"))

