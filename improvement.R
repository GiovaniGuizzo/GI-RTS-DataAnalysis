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
strategies <- c("none", "ekstazi", "starts", "random")

# Number of runs
runs <- 20

# Reads result files
getCsvFileFunction <- function(algorithm, program, strategy, run) {
  folderPath <- file.path(paste(algorithm, "results", sep = "-"),
                          program,
                          strategy)
  csvResults <- folderPath %>%
    file.path(paste(algorithm, "_result_", run, ".csv", sep = "")) %>%
    read_csv(
      col_types = cols(
        MethodName = col_character(),
        MethodIndex = col_integer(),
        Patch = col_character(),
        Compiled = col_logical(),
        AllTestsPassed = col_logical(),
        `TotalExecutionTime(ms)` = col_double(),
        Fitness = col_double(),
        FitnessImprovement = col_double(),
        TimeStamp = col_double()
      )
    )
  
  # Gets only the ones that passed
  csvResults %>%
    filter(Compiled == TRUE, AllTestsPassed == TRUE) %>%
    # Add extra info
    mutate(
      algorithm = algorithm,
      program = program,
      strategy = strategy,
      run = run
    )
}

# Reads results for the none strategy
none_metrics <- map_dfr(algorithms, function(algorithm) {
  map_dfr(programs, function(program) {
    # Each program has multiple strategies
    map_dfr(strategies[1], function(strategy) {
      # Each strategy is ran for 20 independent runs
      map_dfr(1:runs, function(run) {
        getCsvFileFunction(algorithm, program, strategy, run)
      })
    })
  })
})

# Reads patch results
getCsvPatchFileFunction <-
  function(algorithm, program, strategy, run) {
    folderPatchPath <-
      file.path(paste(algorithm, "patch", "results", sep = "-"),
                program,
                strategy)
    
    csvPatchResults <- folderPatchPath %>%
      file.path(paste(algorithm, "_patch_result_", run, ".csv", sep = "")) %>%
      read_csv(
        col_types = cols(
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
        )
      )
    
    # Add extra info
    csvPatchResults %>%
      mutate(
        algorithm = algorithm,
        program = program,
        strategy = strategy,
        run = run
      )
  }

# Declare the treatment function
# It will be executed for each program, strategy, and run
treatmentFunction <- function(algorithm, upToWhichStrategy) {
  map_dfr(programs, function(program) {
    # Each program has multiple strategies
    map_dfr(strategies[2:upToWhichStrategy], function(strategy) {
      # Each strategy is ran for 20 independent runs
      map_dfr(1:runs, function(run) {
        # Now gets results from patch analyser
        getCsvPatchFileFunction(algorithm, program, strategy, run) %>%
          select(-NTests,-NPassed,-NFailed)
      })
    })
  })
}

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
  "spatial4j"
)

# Reads everything
metrics <- map_dfr(algorithms, treatmentFunction, 4)

programs <- c("joda-time")
strategies <- c("none", "ekstazi", "random")

metrics <-
  bind_rows(metrics, map_dfr(algorithms, treatmentFunction, 3))

metrics <- bind_rows(metrics, none_metrics)

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
strategies <- c("none", "ekstazi", "starts", "random")

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
metrics$strategy <-
  factor(metrics$strategy,
         levels = c("Gin", "Ekstazi", "STARTS", "Random"))

# Compute best for each algorithm and program
avg_none <- metrics %>%
  filter(strategy == "Gin",
         Compiled == TRUE,
         AllTestsPassed == TRUE,
         Patch != "|") %>%
  group_by(algorithm, program, run) %>%
  summarise(best = max(FitnessImprovement)) %>%
  group_by(algorithm, program) %>%
  # Then average across all runs
  summarise(avg_best = median(best))

# Compute RIC
addedRIC <- algorithms %>%
  map_dfr(function(subAlgorithm) {
    programs %>%
      map_dfr(function(subProgram) {
        avgBest <-
          filter(avg_none, program == subProgram, algorithm == subAlgorithm)$avg_best
        metrics %>%
          filter(
            program == subProgram,
            algorithm == subAlgorithm,
            Compiled == TRUE,
            AllTestsPassed == TRUE,
            Patch != "|"
          ) %>%
          group_by(algorithm, program, strategy, run) %>%
          slice_max(FitnessImprovement, n = 1, with_ties = FALSE) %>%
          mutate(RIC = FitnessImprovement / avgBest)
      })
  })

# Transform algorithm name to upper case
addedRIC <- addedRIC %>%
  mutate(algorithm = toupper(algorithm))

formatted_table <- addedRIC %>%
  group_by(algorithm, program, strategy) %>%
  # Find the median RIC for each approach
  summarise(median = format(median(RIC), digits = 2, nsmall = 2)) %>%
  # Spread the strategies over columns
  pivot_wider(
    names_from = c(strategy),
    values_from = median,
    names_sep = "+"
  )

# addedRIC %>%
#   group_by(algorithm, program, strategy) %>%
#   summarise(median_RIC = format(median(RIC), digits = 3)) %>%
#   pivot_wider(names_from = c(algorithm, strategy), values_from = median_RIC, names_sep = "+") %>%
#   ungroup() %>%
#   mutate(across(-program, as.numeric)) %>%
#   summarise(across(-program, median, na.rm = TRUE))

pvalue_table <- algorithms %>%
  sort() %>%
  map_dfr(function(sub_algorithm) {
    programs %>%
      sort() %>%
      map_dfr(function(sub_program) {
        # Filter the data to match each algorithm-program tuple
        filtered_data <- addedRIC %>%
          filter(program == sub_program, algorithm == toupper(sub_algorithm))
        # Get p-value
        kruskal_result = format(
          kruskal.test(RIC ~ strategy, data = filtered_data)$p.value,
          digits = 4,
          nsmall = 2
        )
        # Generate small data frame with p-value
        data.frame(
          algorithm = sub_algorithm,
          program = sub_program,
          pvalue = format(round(as.numeric(kruskal_result), 3), nsmall = 3)
        )
      })
  }) %>%
  # Transform algorithm name to upper case
  mutate(algorithm = toupper(algorithm))

# Concatenate p-value as column
formatted_table <-
  left_join(formatted_table, pvalue_table, by = c("program", "algorithm"))

# Find maximum value between approaches and add as column
formatted_table$max <-
  apply(
    formatted_table[3:6] %>% mutate_if(is.character, as.numeric) %>% replace(is.na(.), 0),
    1,
    max
  )

print(
  formatted_table %>%
    # Round maximum value
    mutate(max = format(round(as.numeric(
      max
    ), 2), nsmall = 2)) %>%
    # Add bold to highest value
    mutate(across(1:4, ~ ifelse(
      . == max, paste('\\textbf{', ., '}', sep = ""), .
    ))) %>%
    # Add gray cell to values when p-value > 0.05
    mutate(across(1:4, ~ replace(
      ., pvalue > 0.05, paste('\\cellcolor{black!25}', .)
    ))) %>%
    # Fix 0 values to < 0.001
    mutate(pvalue = ifelse(pvalue == "0.000", "< 0.001", pvalue)) %>%
    # Bold significant p-values
    mutate(pvalue = ifelse(
      pvalue < 0.05, paste('\\textbf{', pvalue, '}', sep = ""), pvalue
    )) %>%
    # Drop max column
    select(-max) %>%
    # Generate latex table
    xtable(),
  include.rownames = FALSE,
  booktabs = TRUE,
  sanitize.text.function = identity
)

# Compute pairwise difference
mapply(function(sub_algorithm) {
  mapply(function(sub_program) {
    filtered_data <- addedRIC %>%
      filter(program == sub_program, algorithm == toupper(sub_algorithm))
    print(sub_program)
    print(sub_algorithm)
    print(kruskalmc(RIC ~ strategy, data = filtered_data))
  }, sort(programs))
}, sort(algorithms))

# Magnitude abbreviations
magnitudes = c("N", "S", "M", "L")

# Generate VDA table
vda_table <- algorithms %>%
  sort() %>%
  map_dfr(function(sub_algorithm) {
    programs %>%
      sort() %>%
      map_dfr(function(sub_program) {
        # Assign pairwise of groups
        list(c("Gin", "Ekstazi"),
             c("Gin", "STARTS"),
             c("Ekstazi", "STARTS")) %>%
          map_dfr(function(groups) {
            # Filter for the specific case
            filtered_data <- addedRIC %>%
              filter(program == sub_program,
                     algorithm == toupper(sub_algorithm))
            # Get groups
            groupA <- groups[1]
            groupB <- groups[2]
            # Get only the results of specific groups
            A <- filtered_data %>% filter(strategy == groupA)
            B <- filtered_data %>% filter(strategy == groupB)
            # Compute VDA
            vda_result <- VD.A(A$RIC, B$RIC)
            magnitude <- vda_result$magnitude
            estimate <- vda_result$estimate
            # Generate small dataframe with results
            data.frame(
              algorithm = sub_algorithm,
              program = sub_program,
              comparison = paste(groupA, groupB, sep = "/"),
              result = paste(format(round(estimate, 2), nsmall = 2), " (", magnitudes[magnitude], ")", sep = "")
            )
          })
      })
  }) %>%
  # Replace Gin with GI
  mutate(comparison = str_replace(comparison, "Gin", "GI")) %>%
  # Spread comparison names as columns in the table
  pivot_wider(names_from = comparison, values_from = result)

# print VDA table
print(
  vda_table %>%
    # Transform algorithm name to upper case
    mutate(algorithm = toupper(algorithm)) %>%
    # Bold for Large differences
    mutate(across(3:5, ~ ifelse(
      str_ends(., "L\\)"), paste('\\textbf{', ., '}', sep = ""), .
    ))) %>%
    # Generate latex table
    xtable(),
  include.rownames = FALSE,
  booktabs = TRUE,
  sanitize.text.function = identity
)

# Save CSV file
write_csv(addedRIC, "improvement-data.csv")
