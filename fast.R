library(tidyverse)
library(GGally)
library(xtable)
library(pgirmess)
library(effsize)
library(hms)
library(lubridate)
library(ggplot2)

# Algorithms used in the experiments
algorithms <- c("gp", "ls")

# Programs used in the experiment
programs <- c(
  "commons-codec",
  "commons-compress",
  "commons-csv",
  "commons-fileupload",
  "commons-imaging",
  "commons-text",
  "commons-validator",
  "gson",
  "jcodec",
  "jfreechart",
  "joda-time",
  "spatial4j"
)
strategies <- c("none", "ekstazi", "starts")

# Number of runs
runs <- 20

# Reads result files
getCsvFileFunction <- function(algorithm, program, strategy, run){
  
  folderPath <- file.path(paste(algorithm, "results", sep = "-"), 
                          program, 
                          strategy)
  csv_path <- folderPath %>%
    file.path(paste(algorithm, "_result_", run, ".csv", sep = ""))
  
  csvResults <- csv_path %>%
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
        folderPath <- file.path(paste(algorithm, "results", sep = "-"), 
                                program, 
                                strategy)
        csv_path <- folderPath %>%
          file.path(paste(algorithm, "_result_", run, ".csv", sep = ""))
        # Now gets results from patch analyser
        if(file.exists(csv_path)) {
          getCsvFileFunction(algorithm, program, strategy, run)
        }
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
metrics$index <- rownames(metrics)

original_programs <- metrics %>%
  filter(Patch == "|")

valid_patches <- metrics %>%
  filter(Patch != "|", FitnessImprovement > 0, Compiled == TRUE, AllTestsPassed == TRUE)

addedTimeToFind <- algorithms %>%
  map_dfr(function(subAlgorithm){
    programs %>% 
      map_dfr(function(subProgram){
        strategies %>%
          map_dfr(function(subStrategy){
            map_dfr(1:runs, function(subRun){
              print(paste(subAlgorithm, subProgram, subStrategy, subRun, sep = " | "))
              
              # Finds the first patch of the list, which is the original program
              original_program <- original_programs %>%
                filter(program == subProgram,
                       algorithm == subAlgorithm,
                       strategy == subStrategy,
                       run == subRun)
              # Saves the index
              original_program_index <- original_program$index
              
              # Finds the first valid and improving patch of the list that is not the original patch |
              first_valid_patch <- valid_patches %>%
                filter(program == subProgram,
                       algorithm == subAlgorithm,
                       strategy == subStrategy,
                       run == subRun) %>%
                slice_min(index, n = 1)
              
              if(nrow(first_valid_patch) > 0){
                # Saves the index
                first_valid_patch_index <- first_valid_patch$index
                # Compute the time to find as the sum of execution time of all previous patches
                first_valid_patch$time_to_find <- (metrics %>% 
                                                           slice(original_program_index:first_valid_patch_index) %>% 
                                                           summarise(sum_times = sum(`TotalExecutionTime(ms)`)))$sum_times
              }
              first_valid_patch
            })
          })
      })
  })

median_metrics <- addedTimeToFind %>%
  group_by(algorithm, program, strategy, run) %>%
  summarise(time_to_find = median(time_to_find/1000)) %>%
  mutate(algorithm = toupper(algorithm))

# Compute Kruskal-Wallis
# Compute pairwise difference
kruskal_result <- mapply(function(sub_algorithm) {
  mapply(function(sub_program) {
    filtered_data <- median_metrics %>%
      filter(program == sub_program, algorithm == toupper(sub_algorithm))
    print(sub_program)
    print(sub_algorithm)
    print(format(kruskal.test(time_to_find ~ strategy, data = filtered_data)$p.value, digits=4, nsmall = 2))
    print(kruskalmc(time_to_find ~ strategy, data = filtered_data))
  }, sort(programs))
}, sort(algorithms))

median_metrics <- median_metrics %>%
  group_by(algorithm, program, strategy) %>%
  summarise(avg_time_to_find = format(round(median(time_to_find), digits = 2), nsmall = 2))

formatted_table <- median_metrics %>%
  # Spread the strategies over columns
  pivot_wider(
    names_from = c(strategy),
    values_from = avg_time_to_find,
    names_sep = "+"
  )

# Find minimum value between approaches and add as column
formatted_table$min <-
  apply(
    formatted_table[3:5] %>% mutate_if(is.character, as.numeric) %>% replace(is.na(.), Inf),
    1,
    min
  )

print(
  formatted_table %>%
    # Round minimum value
    mutate(min = format(round(as.numeric(min), 2), digits = 2, nsmall = 2)) %>%
    # Add bold to highest value
    mutate(across(1:3, ~ ifelse(. == min, paste('\\textbf{', as.character(.), '}', sep = ""), as.character(.)))) %>%
    # Drop min column
    select(-min) %>%
    # Generate latex table
    xtable(),
  include.rownames = FALSE,
  booktabs = TRUE,
  sanitize.text.function = identity
)

median_metrics %>%
  group_by(algorithm, strategy) %>%
  summarise(avg_time_to_find = format(round(median(as.numeric(avg_time_to_find)), digits = 2), nsmall = 2)) %>%
  spread(strategy, avg_time_to_find)

median_metrics %>%
  group_by(strategy) %>%
  summarise(avg_time_to_find = format(round(median(as.numeric(avg_time_to_find)), digits = 2), nsmall = 2)) %>%
  spread(strategy, avg_time_to_find)

# Save CSV file
write_csv(median_metrics, "fast.csv")
