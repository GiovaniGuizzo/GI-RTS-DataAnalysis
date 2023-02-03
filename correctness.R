library(tidyverse)
library(GGally)
library(dplyr)

# Algorithms used in the experiments
algorithms <- c("gp", "ls")

# Programs used in the experiment
programs <- c("commons-codec",
              "commons-compress",
              "commons-csv",
              "commons-fileupload",
              "commons-imaging",
              "commons-text",
              "commons-validator",
              "gson",
              "jcodec",
              "jfreechart",
              "jgrapht",
              # "joda-time",
              "spatial4j")

# Strategies used in the experiment
strategies <- c("ekstazi", "starts", "random")

# Declare the treatment function
# It will be executed for each program, strategy, and run
treatmentFunction <- function(algorithm, runs){
  map_dfr(programs, function(program){
    # Each program has multiple strategies
    map_dfr(strategies, function(strategy){
      # Each strategy is ran for 20 independent runs
      map_dfr(1:runs, function(run){
        # First find the original algorithm results
        folderPath <- file.path(paste(algorithm, "results", sep = "-"), 
                                     program, 
                                     strategy)
        
        csvResults <- folderPath %>%
          file.path(paste(algorithm, "_result_", run, ".csv", sep = "")) %>%
          read_csv(col_types = cols(
            MethodName = col_character(),
            MethodIndex = col_integer(),
            Patch = col_character(),
            Compiled = col_logical(),
            AllTestsPassed = col_logical(),
            `TotalExecutionTime(ms)` = col_double(),
            Fitness = col_double(),
            FitnessImprovement = col_double(),
            TimeStamp = col_double()
          ))
        
        # Gets only the ones that passed
        csvResults <- csvResults %>%
          filter(Compiled == TRUE, AllTestsPassed == TRUE) %>%
          # Add extra info
          mutate(algorithm = algorithm,
                 program = program,
                 strategy = strategy,
                 run = run)
          # Rename columns to Original
        csvResults <- csvResults %>%
          rename(OriginalFitnessImprovement = FitnessImprovement)
        
        # Now gets results from patch analyser
        folderPatchPath <- file.path(paste(algorithm, "patch", "results", sep = "-"),
                                program,
                                strategy)

        csvPatchResults <- folderPatchPath %>%
          file.path(paste(algorithm, "_patch_result_", run, ".csv", sep = "")) %>%
          read_csv(col_types = cols(
            MethodName = col_character(),
            MethodIndex = col_integer(),
            Patch = col_character(),
            Compiled = col_logical(),
            AllTestsPassed = col_logical(),
            `TotalExecutionTime(ms)` = col_double(),
            Fitness = col_double(),
            FitnessImprovement = col_double(),
            TimeStamp = col_double(),
            NTests = col_double(),
            NPassed = col_double(),
            NFailed = col_double()
          ))

        # Add extra info
        csvPatchResults <- csvPatchResults %>%
          mutate(algorithm = algorithm,
                 program = program,
                 strategy = strategy,
                 run = run)

        # Maintain all the information from the patches, but also adds the
        # original results. Since the column FitnessImprovement is named
        # differently, it is kept in the new data frame
        left_join(csvPatchResults, csvResults,
                  by = c("Patch", "algorithm", "program", "strategy", "run"))
      })
    })
  })
}

# Reads everything
metrics <- map_dfr(algorithms, treatmentFunction, 20)

# Computes RS
addedRS <- metrics %>%
  filter(Patch != "|") %>%
  # Groups the patches per run
  group_by(algorithm, program, strategy, run) %>%
  # Finds the patch with the best original improvement
  slice_max(OriginalFitnessImprovement, n = 1, with_ties = FALSE) %>%
  # filter(OriginalFitnessImprovement > 0) %>%
  # Computes the relative safety of the patch. If a test case failed, then
  # RS < 1. If the patch passed with all test cases, then RS == 1.
  mutate(RS = NPassed / NTests)

# Save CSV file
write_csv(addedRS, "correctness-data.csv")

print(addedRS %>%
  group_by(algorithm, program, strategy, run) %>%
  summarise(mean_RS = format(RS, digits = 5),
            Sum_NFailed = sum(NFailed),
            Count_NFailed = length(NFailed[NFailed>0])) %>%
  filter(Sum_NFailed > 0), n = Inf)
