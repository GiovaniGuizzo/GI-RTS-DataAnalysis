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

# Start loop with all programs and 20 independent runs
metrics <- map_dfr(algorithms, treatmentFunction, 20) %>%
  mutate(total_time = profile_time + exec_time)

programs <- c("joda-time")
strategies <- c("none", "ekstazi", "random")

# Start loop with all programs and 20 independent runs
metrics_joda <- map_dfr(algorithms, treatmentFunction, 20) %>%
  mutate(total_time = profile_time + exec_time)

metrics <- bind_rows(metrics, metrics_joda)

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

# Cleanup names
metrics <- metrics %>%
  mutate(program = str_replace_all(program, "commons-", ""))
programs <- str_replace_all(programs, "commons-", "")

# Compute average for each algorithm and program
avg_none <- metrics %>%
  filter(strategy == "none") %>%
  group_by(algorithm, program) %>%
  summarise(avg_time = median(total_time))

# Compute RC
addedRC <- algorithms %>%
  map_dfr(function(subAlgorithm){
    programs %>% 
      map_dfr(function(subProgram){
        averageTime <- filter(avg_none, program == subProgram, algorithm == subAlgorithm)$avg_time
        metrics %>%
          filter(program == subProgram, algorithm == subAlgorithm) %>%
          mutate(RC = total_time / averageTime)
      })
  })

# Transform algorithm name to upper case
addedRC <- addedRC %>%
  mutate(algorithm = toupper(algorithm))

addedRC <- addedRC[addedRC$strategy != "none" & addedRC$strategy != "random", ]

addedRC <- addedRC %>%
  mutate(strategy = str_replace(strategy, "ekstazi", "Ekstazi"))
addedRC <- addedRC %>%
  mutate(strategy = str_replace(strategy, "none", "GI"))
addedRC <- addedRC %>%
  mutate(strategy = str_replace(strategy, "starts", "STARTS"))
addedRC <- addedRC %>%
  mutate(strategy = str_replace(strategy, "random", "Random"))
# Order
addedRC$strategy <- factor(addedRC$strategy, levels = c("GI", "Ekstazi", "STARTS", "Random"))

formatted_table <- addedRC %>%
  group_by(algorithm, program, strategy)

print(formatted_table)

ylim1 = boxplot.stats(formatted_table$RC)$stats[c(1, 5)]
# Create the boxplot
plot_by_algorithm <- ggplot(formatted_table, aes(x = strategy, y = RC, fill = algorithm)) +
  geom_boxplot() +
  geom_hline(yintercept = 1, color="red", linetype="dashed") +
  labs(x = "Strategy", y = "RC") +
  facet_wrap(~program, nrow = 2) +
  theme_bw() +
  scale_y_continuous(breaks = seq(0,2.6,0.2)) +
  coord_cartesian(ylim = ylim1 * 1.3) +
  theme(panel.grid.minor=element_blank()) +
  theme(panel.grid.major.x=element_blank()) +
  theme(strip.background = element_rect(fill = "white")) +
  theme(text = element_text(size = 18)) +
  theme(axis.title.y = element_text(margin = margin(r = 10))) +
  theme(axis.title.x = element_text(margin = margin(t = 10))) +
  guides(fill=guide_legend(title="Algorithm"))

# Save plot
dir.create("images", showWarnings = FALSE)
ggsave("images/by_algorithm.png",
       plot = plot_by_algorithm,
       device = png(),
       width = 700,
       height = 400,
       units = "mm",
       scale = 0.6,
       dpi = "retina")
# # Remove Gin without RTS from the graph
# prunedRC <- addedRC %>%
#   filter(strategy != "GI")
# 
# # Modify names of strategies
# prunedRC <- prunedRC %>%
#   mutate(algorithm = str_replace(algorithm, "gp", "GP"))
# prunedRC <- prunedRC %>%
#   mutate(algorithm = str_replace(algorithm, "ls", "LS"))
# # Order
# prunedRC$algorithm <- factor(prunedRC$algorithm, levels = c("GP", "LS"))
# 
# # Compute lower and upper whiskers to crop weird outliers
# ylim1 = boxplot.stats(prunedRC$RC)$stats[c(1, 5)]
# 
# # Generate plot
# plotRC <- prunedRC %>%
#   ggplot(aes(x = algorithm, y = RC, fill = strategy)) +
#   geom_boxplot() +
#   geom_hline(yintercept = 1, color="red", linetype="dashed") +
#   labs(x = "Strategy", y = "Relative Cost (RC)")+
#   scale_y_continuous(breaks = seq(0,2.2,0.2)) +
#   coord_cartesian(ylim = ylim1*1.05) +
#   theme_bw() +
#   theme(panel.grid.minor=element_blank()) +
#   theme(panel.grid.major.x=element_blank()) +
#   theme(strip.background = element_rect(fill = "white")) +
#   theme(text = element_text(size = 18)) +
#   theme(axis.title.y = element_text(margin = margin(r = 10))) +
#   theme(axis.title.x = element_text(margin = margin(t = 10))) +
#   # theme(axis.text.x = element_text(angle = 45, hjust=1)) +
#   facet_wrap(~program, nrow = 2)
# print(plotRC)
# 
# # Save plot
# dir.create("images", showWarnings = FALSE)
# ggsave("images/RC.png",
#        plot = plotRC,
#        device = png(),
#        width = 700,
#        height = 400,
#        units = "mm",
#        scale = 0.6,
#        dpi = "retina")
# 
# treatedMetrics <- metrics %>%
#   group_by(algorithm, program, strategy) %>%
#   summarise(`Overhead` = as.numeric(dmilliseconds(sum(profile_time)), "hours"),
#             `Optimisation Time` = as.numeric(dmilliseconds(sum(exec_time)), "hours"))
# 
# treatedMetrics <- treatedMetrics %>%
#   pivot_longer(c(`Overhead`, `Optimisation Time`), names_to = "Phase", values_to = "Value") %>%
#   mutate(Phase = paste(str_to_upper(algorithm), Phase, sep = " ")) %>%
#   filter(!str_starts(Phase, "GP Overhead")) %>%
#   mutate(Phase = str_replace(Phase, "LS Overhead", "Overhead"))
# 
# # Modify names of strategies
# treatedMetrics <- treatedMetrics %>%
#   mutate(strategy = str_replace(strategy, "ekstazi", "GI+Ekstazi"))
# treatedMetrics <- treatedMetrics %>%
#   mutate(strategy = str_replace(strategy, "none", "GI"))
# treatedMetrics <- treatedMetrics %>%
#   mutate(strategy = str_replace(strategy, "starts", "GI+STARTS"))
# treatedMetrics <- treatedMetrics %>%
#   mutate(strategy = str_replace(strategy, "random", "GI+Random"))
# # Order
# treatedMetrics$strategy <- factor(treatedMetrics$strategy, levels = c("GI", "GI+Ekstazi", "GI+STARTS", "GI+Random"))
# 
# treatedMetrics <- treatedMetrics %>%
#   group_by(program, strategy) %>%
#   mutate(total = sum(`Value`))
# 
# plotCumulative <- treatedMetrics %>%
#   ggplot(aes(x=strategy)) + 
#   geom_bar(aes(fill=Phase, y=Value), position="stack", stat="identity") +
#   geom_text(aes(label=round(total, digits = 2), y=total), position=position_dodge(width=0.9), vjust=-0.25) +
#   theme_bw() +
#   theme(panel.grid.minor=element_blank()) +
#   theme(panel.grid.major.x=element_blank()) +
#   theme(strip.background = element_rect(fill = "white")) +
#   theme(text = element_text(size = 18)) +
#   theme(legend.title = element_blank()) +
#   facet_wrap(~program, nrow = 1) +
#   labs(x = "Strategy", y = "Cumulative execution time (h)") + 
#   theme(axis.title.y = element_text(margin = margin(r = 10))) +
#   theme(axis.title.x = element_text(margin = margin(t = 10))) +
#   theme(axis.text.x = element_text(angle = 45, hjust=1)) +
#   scale_y_continuous(breaks = seq(0,150,10)) +
#   scale_fill_brewer(palette="Spectral", aesthetics = "colour") +
#   facet_wrap(~program, nrow = 2)
# 
# # Save plot
# ggsave("images/cumulative.png",
#        plot = plotCumulative,
#        device = png(),
#        width = 700,
#        height = 400,
#        units = "mm",
#        scale = 0.6,
#        dpi = "retina")
# 
# treatedMetrics %>%
#   group_by(strategy) %>%
#   summarise(total_weeks = dhours(sum(total)), duration = as.numeric(dhours(sum(total)), "hous"))
