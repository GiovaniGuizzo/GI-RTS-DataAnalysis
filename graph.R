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
programs <- c("commons-codec",
              "commons-compress",
              "commons-csv",
              "commons-fileupload",
              "commons-imaging",
              "commons-text",
              "commons-validator")

# Strategies used in the experiment
strategies <- c("none", "ekstazi", "starts", "random")

# Declare the treatment function
# It will be executed for each program, strategy, and run
treatmentFunction <- function(algorithm, runs){
  map_dfr(programs, function(program){
    # Each program has multiple strategies
    map_dfr(strategies, function(strategy){
      # Each strategy is ran for 20 independent runs
      map_dfr(1:runs, function(run){
        profileFolderPath <- file.path("profile-results", 
                                       program, 
                                       strategy)
        profileTiming <- profileFolderPath %>%
          file.path(paste("profile_timing_", run, ".csv", sep = "")) %>%
          read_csv(col_names = c("profile_time"),
                   col_types = cols(profile_time = col_double()))
        
        folderPath <- file.path(paste(algorithm, "results", sep = "-"), 
                                program, 
                                strategy)
        
        csvTiming <- folderPath %>%
          file.path(paste(algorithm, "_timing_", run, ".csv", sep = "")) %>%
          read_csv(col_names = c("exec_time"),
                   col_types = cols(exec_time = col_double()))
        
        bind_cols(profileTiming, csvTiming) %>%
          mutate(algorithm = algorithm,
                 program = program,
                 strategy = strategy,
                 run = run)
      })
    })
  })
}

# Start loop with all programs and 20 independent runs
metrics <- map_dfr(algorithms, treatmentFunction, 20) %>%
  mutate(total_time = profile_time + exec_time)
# Cleanup names
metrics <- metrics %>%
  mutate(program = str_replace_all(program, "commons-", ""))
programs <- str_replace_all(programs, "commons-", "")
# Change names
metrics <- metrics %>%
  mutate(strategy = str_replace(strategy, "ekstazi", "Ekstazi"))
metrics <- metrics %>%
  mutate(strategy = str_replace(strategy, "none", "GI"))
metrics <- metrics %>%
  mutate(strategy = str_replace(strategy, "starts", "STARTS"))
metrics <- metrics %>%
  mutate(strategy = str_replace(strategy, "random", "Random"))

# Order
metrics$strategy <- factor(metrics$strategy, levels = c("GI", "Ekstazi", "STARTS", "Random"))

treatedMetrics <- metrics %>%
  group_by(program, strategy) %>%
  summarise(`Overhead` = as.numeric(dmilliseconds(sum(profile_time)), "hours"),
            `Optimisation Time` = as.numeric(dmilliseconds(sum(exec_time)), "hours"),
            `Total Time` = as.numeric(dmilliseconds(sum(total_time)), "hours"))

treatedMetrics <- treatedMetrics %>%
  gather(key = "Phase", value = "Value", `Overhead`, `Optimisation Time`)

plotCumulative <- treatedMetrics %>%
  ggplot(aes(x=strategy)) + 
  geom_bar(aes(fill=Phase, y=Value), position="stack", stat="identity") +
  geom_text(aes(label=round(`Total Time`, digits = 2), y=`Total Time`), position=position_dodge(width=0.9), vjust=-0.25) +
  theme_bw() +
  theme(panel.grid.minor=element_blank()) +
  theme(panel.grid.major.x=element_blank()) +
  theme(strip.background = element_rect(fill = "white")) +
  theme(text = element_text(size = 18)) +
  theme(legend.title = element_blank()) +
  facet_wrap(~program, nrow = 1) +
  labs(x = "Strategy", y = "Cumulative execution time (h)") + 
  theme(axis.title.y = element_text(margin = margin(r = 10))) +
  theme(axis.title.x = element_text(margin = margin(t = 10))) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  scale_y_continuous(breaks = seq(0,110,10)) +
  scale_fill_brewer(palette="Spectral", aesthetics = "colour")

ggsave("C:/Users/giova/Projects/images/cumulative.png",
       plot = plotCumulative,
       device = png(),
       width = 700,
       height = 200,
       units = "mm",
       scale = 0.6,
       dpi = "retina")
