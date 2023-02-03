library(tidyverse)
library(GGally)
library(xtable)
library(pgirmess)
library(effsize)
library(hms)
library(lubridate)

# Algorithms used in the experiments
algorithms <- c("ls")

# Number of runs
runs <- 20

# Programs used in the experiment
programs <- c("commons-codec",
              "commons-compress",
              "commons-csv",
              "commons-fileupload",
              "commons-imaging",
              "commons-text",
              "commons-validator")

# Strategies used in the experiment
strategies <- c("none", "ekstazi", "starts")

# Reads result files
getCsvFileFunction <- function(algorithm, program, strategy, run){
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
  csvResults %>%
    filter(Compiled == TRUE, AllTestsPassed == TRUE) %>%
    # Add extra info
    mutate(algorithm = algorithm,
           program = program,
           strategy = strategy,
           run = run)
}

# Reads results for the none strategy
none_metrics <- map_dfr(algorithms, function(algorithm){
  map_dfr(programs, function(program){
    # Each program has multiple strategies
    map_dfr(strategies[1], function(strategy){
      # Each strategy is ran for 20 independent runs
      map_dfr(1:runs, function(run){
        getCsvFileFunction(algorithm, program, strategy, run) %>%
          filter(FitnessImprovement > 0)
      })
    })
  })
})

# Reads patch results
getCsvPatchFileFunction <- function(algorithm, program, strategy, run){
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
  csvPatchResults %>%
    mutate(algorithm = algorithm,
           program = program,
           strategy = strategy,
           run = run)
}

# Declare the treatment function
# It will be executed for each program, strategy, and run
treatmentFunction <- function(algorithm){
  map_dfr(programs, function(program){
    # Each program has multiple strategies
    map_dfr(strategies[2:3], function(strategy){
      # Each strategy is ran for 20 independent runs
      map_dfr(1:runs, function(run){
        # Now gets results from patch analyser
        getCsvPatchFileFunction(algorithm, program, strategy, run) %>%
          filter(Compiled == TRUE, AllTestsPassed == TRUE, FitnessImprovement > 0) %>%
          select(-NTests, -NPassed, -NFailed)
      })
    })
  })
}

# Reads everything
metrics <- map_dfr(algorithms, treatmentFunction)

metrics <- metrics %>%
  bind_rows(metrics, none_metrics)

# Modify names of strategies
metrics <- metrics %>%
  mutate(strategy = str_replace(strategy, "ekstazi", "Ekstazi"))
metrics <- metrics %>%
  mutate(strategy = str_replace(strategy, "none", "Gin"))
metrics <- metrics %>%
  mutate(strategy = str_replace(strategy, "starts", "STARTS"))

# Cleanup names
metrics <- metrics %>%
  mutate(program = str_replace_all(program, "commons-", ""))
programs <- str_replace_all(programs, "commons-", "")

# Order
metrics$strategy <- factor(metrics$strategy, levels = c("Gin", "Ekstazi", "STARTS"))

metrics <- metrics %>%
  filter(Patch != "|") %>%
  distinct(Patch, .keep_all = TRUE)

median_metrics <- metrics %>%
  group_by(program, strategy, run) %>%
  summarise(count_positive_patches = n())

# Compute Kruskal-Wallis
kruskalFunction <- function(sub_program){
  filtered_data <- median_metrics %>%
    filter(program == sub_program)
  
  print(sub_program)
  print(format(kruskal.test(count_positive_patches ~ strategy, data = filtered_data)$p.value, digits=4, nsmall = 2))
  print(kruskalmc(count_positive_patches ~ strategy, data = filtered_data))
}
mapply(kruskalFunction, programs)

median_metrics <- median_metrics %>%
  group_by(program, strategy) %>%
  summarise(median_count = median(count_positive_patches))

spreaded_data <- median_metrics %>%
  spread(strategy, median_count) %>%
  ungroup() %>%
  select(-program)

plotCD(spreaded_data, alpha = 0.05)

median_metrics %>%
  spread(strategy, median_count) %>%
  xtable() %>%
  print(include.rownames=FALSE, digits = 2)

median_metrics %>%
  group_by(strategy) %>%
  summarise(median_count = median(median_count)) %>%
  spread(strategy, median_count)

# Save CSV file
write_csv(median_metrics, "diverse.csv")
