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
strategies <- c("none", "ekstazi", "starts", "random")

# Number of runs
runs <- 20

metrics <- read.csv(file = "improvement-data.csv") %>%
  filter(strategy != "Random")
programs <- str_replace_all(programs, "commons-", "")

# Order
metrics$strategy <- factor(metrics$strategy, levels = c("Gin", "Ekstazi", "STARTS"))

metrics <- metrics %>%
  filter(Patch != "|")

median_metrics <- metrics %>%
  group_by(algorithm, program, strategy, run) %>%
  summarise(best = max(FitnessImprovement/1000))

# Compute Kruskal-Wallis
# Compute pairwise difference
kruskal_result <- mapply(function(sub_algorithm) {
  mapply(function(sub_program) {
    filtered_data <- median_metrics %>%
      filter(program == sub_program, algorithm == toupper(sub_algorithm))
    print(sub_program)
    print(sub_algorithm)
    print(format(kruskal.test(best ~ strategy, data = filtered_data)$p.value, digits=4, nsmall = 2))
    print(kruskalmc(best ~ strategy, data = filtered_data))
  }, sort(programs))
}, sort(algorithms))

median_metrics <- median_metrics %>%
  group_by(algorithm, program, strategy) %>%
  summarise(avg_best = round(median(best), digits = 2))

formatted_table <- median_metrics %>%
  # Spread the strategies over columns
  pivot_wider(
    names_from = c(strategy),
    values_from = avg_best,
    names_sep = "+"
  )

# Find maximum value between approaches and add as column
formatted_table$max <-
  apply(
    formatted_table[3:5] %>% mutate_if(is.character, as.numeric) %>% replace(is.na(.), 0),
    1,
    max
  )

print(
  formatted_table %>%
    # Round maximum value
    mutate(max = format(round(as.numeric(max), 2), nsmall = 2)) %>%
    # Add bold to highest value
    mutate(across(1:3, ~ ifelse(. == max, paste('\\textbf{', as.character(.), '}', sep = ""), as.character(.)))) %>%
    # Drop max column
    select(-max) %>%
    # Generate latex table
    xtable(),
  include.rownames = FALSE,
  booktabs = TRUE,
  sanitize.text.function = identity
)

median_metrics %>%
  group_by(algorithm, strategy) %>%
  summarise(avg_best = round(median(avg_best), digits = 2)) %>%
  spread(strategy, avg_best)

median_metrics %>%
  group_by(strategy) %>%
  summarise(avg_best = round(median(avg_best), digits = 2)) %>%
  spread(strategy, avg_best)

# Save CSV file
write_csv(median_metrics, "perfect.csv")
