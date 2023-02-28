library(tidyverse)
library(GGally)
library(xtable)
library(pgirmess)
library(effsize)
library(hms)
library(lubridate)

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

# Read result files
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
  
  csvResults %>%
    # Add extra info
    mutate(algorithm = algorithm,
           program = program,
           strategy = strategy,
           run = run)
}

# Read results for the none strategy
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
    map_dfr(strategies[2:3], function(strategy){
      # Each strategy is ran for 20 independent runs
      map_dfr(1:runs, function(run){
        folderPath <- file.path(paste(algorithm, "results", sep = "-"), 
                                program, 
                                strategy)
        csv_path <- folderPath %>%
          file.path(paste(algorithm, "_result_", run, ".csv", sep = ""))
        if(file.exists(csv_path)) {
          # Now gets results from patch analyser
          getCsvPatchFileFunction(algorithm, program, strategy, run) %>%
            select(-NTests, -NPassed, -NFailed)
        }
      })
    })
  })
}

# Reads everything
metrics <- map_dfr(algorithms, treatmentFunction)

metrics <- bind_rows(metrics, none_metrics)

# Modify names of strategies
metrics <- metrics %>%
  mutate(strategy = str_replace(strategy, "ekstazi", "Ekstazi"))
metrics <- metrics %>%
  mutate(strategy = str_replace(strategy, "none", "GI"))
metrics <- metrics %>%
  mutate(strategy = str_replace(strategy, "starts", "STARTS"))

# Cleanup names
metrics <- metrics %>%
  mutate(program = str_replace_all(program, "commons-", ""))
programs <- str_replace_all(programs, "commons-", "")

# Order
metrics$strategy <- factor(metrics$strategy, levels = c("GI", "Ekstazi", "STARTS"))

valid_patches <- metrics %>%
  filter(Patch != "|", FitnessImprovement > 0, Compiled == TRUE, AllTestsPassed == TRUE) %>%
  distinct(Patch, .keep_all = TRUE)

median_metrics <- valid_patches %>%
  group_by(algorithm, program, strategy, run) %>%
  summarise(count_positive_patches = n())

# # Compute Kruskal-Wallis
# kruskalFunction <- function(sub_program){
#   filtered_data <- median_metrics %>%
#     filter(program == sub_program)
#   
#   print(sub_program)
#   print(format(kruskal.test(count_positive_patches ~ strategy, data = filtered_data)$p.value, digits=4, nsmall = 2))
#   print(kruskalmc(count_positive_patches ~ strategy, data = filtered_data))
# }
# mapply(kruskalFunction, programs)

median_metrics <- median_metrics %>%
  group_by(algorithm, program, strategy) %>%
  summarise(median_count = format(round(median(count_positive_patches), digits = 1), nsmall = 1))

formatted_table <- median_metrics %>%
  # Spread the strategies over columns
  pivot_wider(
    names_from = c(strategy),
    values_from = median_count,
    names_sep = "+"
  )

# Find minimum value between approaches and add as column
formatted_table$max <-
  apply(
    formatted_table[3:5] %>% mutate_if(is.character, as.numeric) %>% replace(is.na(.), 0),
    1,
    max
  )

print(
  formatted_table %>%
    # Round maximum value
    mutate(max = format(round(as.numeric(max), 1), digits = 1, nsmall = 1)) %>%
    # Add bold to highest value
    mutate(across(1:3, ~ ifelse(. == max, paste('\\textbf{', as.character(.), '}', sep = ""), as.character(.)))) %>%
    # Drop maximum column
    select(-max) %>%
    # Generate latex table
    xtable(),
  include.rownames = FALSE,
  booktabs = TRUE,
  sanitize.text.function = identity
)

median_metrics %>%
  group_by(algorithm, strategy) %>%
  summarise(median_count = format(round(median(as.numeric(median_count)), digits = 1), nsmall = 1)) %>%
  spread(strategy, median_count)

median_metrics %>%
  group_by(strategy) %>%
  summarise(median_count = format(round(median(as.numeric(median_count)), digits = 1), nsmall = 1)) %>%
  spread(strategy, median_count)

# Save CSV file
write_csv(median_metrics, "diverse.csv")
