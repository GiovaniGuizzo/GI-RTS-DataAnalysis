library(tidyverse)
library(GGally)
library(xtable)
library(pgirmess)
library(effsize)
library(hms)
library(lubridate)

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
strategies <- c("none", "ekstazi", "starts", "random")

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
        getCsvFileFunction(algorithm, program, strategy, run)
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
    map_dfr(strategies[2:4], function(strategy){
      # Each strategy is ran for 20 independent runs
      map_dfr(1:runs, function(run){
        # Now gets results from patch analyser
        getCsvPatchFileFunction(algorithm, program, strategy, run) %>%
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
metrics <- metrics %>%
  mutate(strategy = str_replace(strategy, "random", "Random"))

# Cleanup names
metrics <- metrics %>%
  mutate(program = str_replace_all(program, "commons-", ""))
programs <- str_replace_all(programs, "commons-", "")

# Order
metrics$strategy <- factor(metrics$strategy, levels = c("Gin", "Ekstazi", "STARTS", "Random"))

# Compute best for each algorithm and program
avg_none <- metrics %>%
  filter(strategy == "Gin", Compiled == TRUE, AllTestsPassed == TRUE, Patch != "|") %>%
  group_by(algorithm, program, run) %>%
  summarise(best = max(FitnessImprovement)) %>%
  group_by(algorithm, program) %>%
  # Then average across all runs
  summarise(avg_best = median(best))

# Compute RIC
addedRIC <- algorithms %>%
  map_dfr(function(subAlgorithm){
    programs %>% 
      map_dfr(function(subProgram){
        avgBest <- filter(avg_none, program == subProgram, algorithm == subAlgorithm)$avg_best
        print(avgBest)
        metrics %>%
          filter(program == subProgram,
                 algorithm == subAlgorithm,
                 Compiled == TRUE,
                 AllTestsPassed == TRUE,
                 Patch != "|") %>%
          group_by(algorithm, program, strategy, run) %>%
          slice_max(FitnessImprovement, n = 1, with_ties = FALSE) %>%
          mutate(RIC = FitnessImprovement / avgBest)
      })
  })

addedRIC %>%
  group_by(program, strategy) %>%
  summarise(median = format(median(RIC), digits = 3)) %>%
  spread(strategy, median) %>%
  xtable()

addedRIC %>%
  group_by(program, strategy) %>%
  summarise(medain_RIC = format(median(RIC), digits = 3)) %>%
  spread(strategy, medain_RIC) %>%
  ungroup() %>%
  summarise(Gin = median(Gin),
            Ekstazi = median(Ekstazi),
            STARTS = median(STARTS),
            Random = median(Random))

kruskalFunction <- function(sub_program){
  filtered_data <- addedRIC %>%
    filter(program == sub_program)
  
  format(kruskal.test(RIC ~ strategy, data = filtered_data)$p.value, digits=4, nsmall = 2)
  # print(kruskalmc(RIC ~ strategy, data = filtered_data))
}
mapply(kruskalFunction, programs)

# Compute VDA
vdaFunction <- function(sub_program, groupA, groupB){
  filtered_data <- addedRIC %>%
    filter(program == sub_program)
  A<-filtered_data%>%filter(strategy == groupA)
  B<-filtered_data%>%filter(strategy == groupB)
  
  vda_result <- VD.A(A$RIC, B$RIC)
  c(format(vda_result$estimate, digits = 2), vda_result$magnitude)
}
mapply(vdaFunction, programs, "Gin", "Ekstazi")
mapply(vdaFunction, programs, "Gin", "STARTS")
mapply(vdaFunction, programs, "Ekstazi", "STARTS")

# Save CSV file
write_csv(addedRIC, "improvement-data.csv")
