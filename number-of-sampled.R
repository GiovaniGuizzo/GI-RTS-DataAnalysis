library(tidyverse)
library(GGally)
library(xtable)
library(pgirmess)
library(effsize)
library(hms)
library(lubridate)

# Algorithms used in the experiments
algorithms <- c("gp")

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
              "joda-time",
              "spatial4j")

# Strategies used in the experiment
strategies <- c("none", "ekstazi", "starts", "random")

runs <- 20

# Start loop with all programs and 20 independent runs
metrics <- map_dfr(algorithms, function(algorithm){
  map_dfr(programs, function(program){
    # Each program has multiple strategies
    map_dfr(strategies, function(strategy){
      map_dfr(1:runs, function(run){
        profileFolderPath <- file.path("profile-results", 
                                       program, 
                                       strategy)
        
        profileTiming <- profileFolderPath %>%
          file.path(paste("profile_result_", run,".csv", sep = "")) %>%
          read_csv(col_types = cols(
            Project = col_character(),
            MethodIndex = col_double(),
            Method = col_character(),
            Count = col_double(),
            Tests = col_character()
          )) %>%
          mutate(algorithm = algorithm,
                 program = program,
                 strategy = strategy,
                 run = run)
      })
    })
  })
})

metrics <- metrics %>%
  group_by(algorithm, strategy, program, run) %>%
  slice_max(Count, n = 1, with_ties = FALSE) %>%
  mutate(nTests = length(str_split(Tests, ",")[[1]]) - 1) %>%
  select(-Tests)

addedMedian <- metrics %>%
  group_by(algorithm, strategy, program) %>%
  summarise(mean_ntests = median(nTests))

addedPercentages <- map_dfr(algorithms, function(subAlgorithm){
  map_dfr(programs, function(subProgram){
        nTestsNone <- filter(addedMedian, algorithm == subAlgorithm, program == subProgram, strategy == "none")$mean_ntests
        metrics %>%
          filter(algorithm == subAlgorithm, program == subProgram) %>%
          mutate(percentage_tests = round((nTests / nTestsNone * 100), digits = 2))
      })
})

write_csv(addedPercentages, "number-of-sampled.csv")

addedPercentages <- addedPercentages %>%
  group_by(program, strategy) %>%
  summarise(median_value = median(percentage_tests))

addedPercentages <- addedPercentages %>%
  mutate(program = str_replace(program, "commons-", "")) %>%
  arrange(program)

print(addedPercentages %>%
  select(program, strategy, median_value) %>%
  spread(strategy, median_value) %>%
  rename(Ekstazi = ekstazi, STARTS = starts, Random = random) %>%
  select(Ekstazi, STARTS, Random) %>%
  xtable(), include.rownames = FALSE, booktabs = TRUE)

print(addedPercentages %>%
  group_by(strategy) %>%
  summarise(median_value = median(median_value)) %>%
  spread(strategy, median_value) %>%
  rename(Ekstazi = ekstazi, STARTS = starts, Random = random) %>%
  select(Ekstazi, STARTS, Random) %>%
  xtable(), include.rownames = FALSE, booktabs = TRUE)
  