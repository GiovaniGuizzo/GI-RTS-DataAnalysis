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
              "commons-validator")

# Strategies used in the experiment
strategies <- c("none", "ekstazi", "starts", "random")

runs <- 20

percentages <- read.csv("number-of-sampled.csv", header = TRUE)

percentages <- percentages %>%
  mutate(strategy = str_replace(strategy, "ekstazi", "Ekstazi"))
percentages <- percentages %>%
  mutate(strategy = str_replace(strategy, "none", "Gin"))
percentages <- percentages %>%
  mutate(strategy = str_replace(strategy, "starts", "STARTS"))
percentages <- percentages %>%
  mutate(strategy = str_replace(strategy, "random", "Random"))

percentages <- percentages %>%
  mutate(program = str_replace_all(program, "commons-", ""))
programs <- str_replace_all(programs, "commons-", "")

rc <- read.csv("timing-data.csv", header = TRUE)

joined <- left_join(percentages, rc, by = c("program", "strategy", "run"))

cor.test(joined$percentage_tests, joined$RC, method = "spearman")
