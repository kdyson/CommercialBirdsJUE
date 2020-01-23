# Script to anonymize bird data to prepare for publication. 

source("../../../../RCode/R_Scripts/anon_data.R")

# bird observations

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

bird.data <- anon.data(
    input.table = bird.data, 
    lookup.table = "../../DataAnalysis/site_lookup.csv",
    site.column = "site.name",
    match.column = "site.name",
    anon.column = "site.abbr",
    output.path = "BirdData_all_07192016.csv",
    write.to.disk = TRUE
)

