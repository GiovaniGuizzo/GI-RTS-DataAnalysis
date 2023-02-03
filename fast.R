library(tidyverse)
library(GGally)
library(xtable)
library(pgirmess)
library(effsize)
library(hms)
library(lubridate)
library(scmamp)
library(ggplot2)
library(Rgraphviz)

# Algorithms used in the experiments
algorithms <- c("gp")

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

# Declare the treatment function
# It will be executed for each program, strategy, and run
treatmentFunction <- function(algorithm){
  map_dfr(programs, function(program){
    # Each program has multiple strategies
    map_dfr(strategies, function(strategy){
      # Each strategy is ran for 20 independent runs
      map_dfr(1:runs, function(run){
        # Now gets results from patch analyser
        getCsvFileFunction(algorithm, program, strategy, run) %>%
          filter(Compiled == TRUE, AllTestsPassed == TRUE)
      })
    })
  })
}

# Reads everything
metrics <- map_dfr(algorithms, treatmentFunction)

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

strategies <- c("Gin", "Ekstazi", "STARTS")

# Order
metrics$strategy <- factor(metrics$strategy, levels = strategies)

first_patch <- metrics %>%
  filter(Patch == "|")

metrics <- metrics %>%
  filter(Patch != "|", FitnessImprovement > 0)

metrics <- algorithms %>%
  map_dfr(function(subAlgorithm){
    programs %>% 
      map_dfr(function(subProgram){
        strategies %>%
          map_dfr(function(subStrategy){
            map_dfr(1:runs, function(subRun){
              first_patch <- filter(first_patch, program == subProgram,
                       algorithm == subAlgorithm,
                       strategy == subStrategy,
                       run == subRun)
              first_patch_timestamp <- first_patch$TimeStamp
              first_patch_time <- first_patch$`TotalExecutionTime(ms)`
              
              metrics %>%
                filter(program == subProgram,
                       algorithm == subAlgorithm,
                       strategy == subStrategy,
                       run == subRun) %>%
                slice_min(TimeStamp, n = 1) %>%
                mutate(time_to_find = TimeStamp - (first_patch_timestamp - first_patch_time))
            })
          })
      })
  })

median_metrics <- metrics %>%
  group_by(program, strategy, run) %>%
  summarise(avg_time_to_find = median(time_to_find/1000))

# Compute Kruskal-Wallis
kruskalFunction <- function(sub_program){
  filtered_data <- median_metrics %>%
    filter(program == sub_program)
  
  print(sub_program)
  print(format(kruskal.test(avg_time_to_find ~ strategy, data = filtered_data)$p.value, digits=4, nsmall = 2))
  print(kruskalmc(avg_time_to_find ~ strategy, data = filtered_data))
}
mapply(kruskalFunction, programs)

median_metrics <- metrics %>%
  group_by(program, strategy) %>%
  summarise(avg_time_to_find = median(time_to_find/1000))

spreaded_data <- median_metrics %>%
  spread(strategy, avg_time_to_find) %>%
  ungroup() %>%
  select(-program)

plotCD(spreaded_data, alpha = 0.05)

median_metrics %>%
  spread(strategy, avg_time_to_find) %>%
  xtable() %>%
  print(include.rownames=FALSE)

median_metrics %>%
  group_by(strategy) %>%
  summarise(avg_time_to_find = median(avg_time_to_find)) %>%
  spread(strategy, avg_time_to_find)

# Save CSV file
write_csv(median_metrics, "fast.csv")
