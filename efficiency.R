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
# Modify names of strategies
addedRC <- addedRC %>%
  mutate(strategy = str_replace(strategy, "ekstazi", "Ekstazi"))
addedRC <- addedRC %>%
  mutate(strategy = str_replace(strategy, "none", "Gin"))
addedRC <- addedRC %>%
  mutate(strategy = str_replace(strategy, "starts", "STARTS"))
addedRC <- addedRC %>%
  mutate(strategy = str_replace(strategy, "random", "Random"))
# Order
addedRC$strategy <- factor(addedRC$strategy, levels = c("Gin", "Ekstazi", "STARTS", "Random"))

# Save CSV file
write_csv(addedRC, "timing-data.csv")

# Get median RC of each strategy
summarisedRC <- addedRC %>%
  group_by(algorithm, program, strategy) %>%
  summarise(median_RC = median(RC)) %>%
  spread(strategy, median_RC)

# Print LaTeX table
summarisedRC %>%
  xtable()

# Get the overall median RC across programs
addedRC %>%
  group_by(algorithm, program, strategy) %>%
  summarise(median_RC = median(RC)) %>%
  spread(strategy, median_RC)  %>%
  ungroup() %>%
  summarise(Gin = median(Gin),
            Ekstazi = median(Ekstazi),
            STARTS = median(STARTS),
            Random = median(Random))


# Get total cost in hms
print("#### TOTAL COST ####")
addedRC %>%
  group_by(algorithm, strategy) %>%
  summarise(sum_cost = as_hms(sum(milliseconds(total_time))))

# Compute Kruskal-Wallis
kruskalFunction <- function(sub_program, sub_algorithm){
  filtered_data <- addedRC %>%
    filter(program == sub_program, algorithm == sub_algorithm)
  
  print(sub_program)
  print(sub_algorithm)
  print(format(kruskal.test(RC ~ strategy, data = filtered_data)$p.value, digits=4, nsmall = 2))
  print(kruskalmc(RC ~ strategy, data = filtered_data))
}
for (algorithm in algorithms) {
  mapply(kruskalFunction, programs, algorithm)
}

# Compute VDA
vdaFunction <- function(sub_program, sub_algorithm, groupA, groupB){
  filtered_data <- addedRC %>%
    filter(program == sub_program, algorithm == sub_algorithm)
  A<-filtered_data%>%filter(strategy == groupA)
  B<-filtered_data%>%filter(strategy == groupB)
  
  vda_result <- VD.A(A$RC, B$RC)
  print(c(format(vda_result$estimate, digits = 2), vda_result$magnitude))
}
for (algorithm in algorithms) {
  mapply(vdaFunction, programs, algorithm, "Gin", "Ekstazi")
  mapply(vdaFunction, programs, algorithm, "Gin", "STARTS")
  mapply(vdaFunction, programs, algorithm, "Ekstazi", "STARTS")
}
# Remove Gin without RTS from the graph
prunedRC <- addedRC %>%
  filter(strategy != "Gin")

# Modify names of strategies
prunedRC <- prunedRC %>%
  mutate(algorithm = str_replace(algorithm, "gp", "GP"))
prunedRC <- prunedRC %>%
  mutate(algorithm = str_replace(algorithm, "ls", "LS"))
# Order
prunedRC$algorithm <- factor(prunedRC$algorithm, levels = c("GP", "LS"))

# Compute lower and upper whiskers to crop weird outliers
ylim1 = boxplot.stats(prunedRC$RC)$stats[c(1, 5)]

# Generate plot
plotRC <- prunedRC %>%
  ggplot(aes(x = algorithm, y = RC, fill = strategy)) +
  geom_boxplot() +
  geom_hline(yintercept = 1, color="red", linetype="dashed") +
  labs(x = "Strategy", y = "Relative Cost (RC)")+
  scale_y_continuous(breaks = seq(0,2.2,0.2)) +
  coord_cartesian(ylim = ylim1*1.05) +
  theme_bw() +
  theme(panel.grid.minor=element_blank()) +
  theme(panel.grid.major.x=element_blank()) +
  theme(strip.background = element_rect(fill = "white")) +
  theme(text = element_text(size = 18)) +
  theme(axis.title.y = element_text(margin = margin(r = 10))) +
  theme(axis.title.x = element_text(margin = margin(t = 10))) +
  # theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  facet_wrap(~program, nrow = 1)
print(plotRC)

# Save plot
ggsave("C:/Users/giova/Projects/images/RC.png",
       plot = plotRC,
       device = png(),
       width = 700,
       height = 200,
       units = "mm",
       scale = 0.6,
       dpi = "retina")

