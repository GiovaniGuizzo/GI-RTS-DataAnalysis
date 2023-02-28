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
# Modify names of strategies
addedRC <- addedRC %>%
  mutate(strategy = str_replace(strategy, "ekstazi", "GI+Ekstazi"))
addedRC <- addedRC %>%
  mutate(strategy = str_replace(strategy, "none", "GI"))
addedRC <- addedRC %>%
  mutate(strategy = str_replace(strategy, "starts", "GI+STARTS"))
addedRC <- addedRC %>%
  mutate(strategy = str_replace(strategy, "random", "GI+Random"))
# Order
addedRC$strategy <- factor(addedRC$strategy, levels = c("GI", "GI+Ekstazi", "GI+STARTS", "GI+Random"))

# Save CSV file
write_csv(addedRC, "timing-data.csv")

# Transform algorithm name to upper case
addedRC <- addedRC %>%
  mutate(algorithm = toupper(algorithm))

formatted_table <- addedRC %>%
  group_by(algorithm, program, strategy) %>%
  # Find the median RC for each approach
  summarise(median = format(median(RC), digits = 2, nsmall = 2)) %>%
  # Spread the strategies over columns
  pivot_wider(
    names_from = c(strategy),
    values_from = median,
    names_sep = "+"
  )

addedRC %>%
  group_by(algorithm, program, strategy) %>%
  summarise(median_RC = format(median(RC), digits = 2, nsmall = 2)) %>%
  pivot_wider(names_from = c(algorithm, strategy), values_from = median_RC, names_sep = "+") %>%
  ungroup() %>%
  mutate(across(-program, as.numeric)) %>%
  summarise(across(-program, median, na.rm = TRUE))

addedRC %>%
  group_by(algorithm, program, strategy) %>%
  summarise(median_RC = format(median(RC), digits = 2, nsmall = 2)) %>%
  pivot_wider(names_from = c(strategy), values_from = median_RC, names_sep = "+") %>%
  ungroup() %>%
  mutate(across(c(-program, -algorithm), as.numeric)) %>%
  summarise(across(c(-program, -algorithm), median, na.rm = TRUE))

pvalue_table <- algorithms %>%
  sort() %>%
  map_dfr(function(sub_algorithm) {
    programs %>%
      sort() %>%
      map_dfr(function(sub_program) {
        # Filter the data to match each algorithm-program tuple
        filtered_data <- addedRC %>%
          filter(program == sub_program, algorithm == toupper(sub_algorithm))
        # Get p-value
        kruskal_result = format(
          kruskal.test(RC ~ strategy, data = filtered_data)$p.value,
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

# Find minimum value between approaches and add as column
formatted_table$min <-
  apply(
    formatted_table[3:6] %>% mutate_if(is.character, as.numeric) %>% replace(is.na(.), 0),
    1,
    min
  )

print(
  formatted_table %>%
    # Round minimum value
    mutate(min = format(round(as.numeric(
      min
    ), 2), nsmall = 2)) %>%
    # Add bold to highest value
    mutate(across(1:4, ~ ifelse(
      . == min, paste('\\textbf{', ., '}', sep = ""), .
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
    # Drop min column
    select(-min) %>%
    # Generate latex table
    xtable(),
  include.rownames = FALSE,
  booktabs = TRUE,
  sanitize.text.function = identity
)

# Get total cost in hms
print("#### TOTAL COST ####")
addedRC %>%
  group_by(algorithm, strategy) %>%
  summarise(sum_cost = as_hms(sum(milliseconds(total_time))))

# Compute pairwise difference
pairwise_result <- mapply(function(sub_algorithm) {
  mapply(function(sub_program) {
    filtered_data <- addedRC %>%
      filter(program == sub_program, algorithm == toupper(sub_algorithm))
    print(sub_program)
    print(sub_algorithm)
    print(kruskalmc(RC ~ strategy, data = filtered_data))
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
        list(c("GI", "GI+Ekstazi"),
             c("GI", "GI+STARTS"),
             c("GI+Ekstazi", "GI+STARTS")) %>%
          map_dfr(function(groups) {
            # Filter for the specific case
            filtered_data <- addedRC %>%
              filter(program == sub_program,
                     algorithm == toupper(sub_algorithm))
            # Get groups
            groupA <- groups[1]
            groupB <- groups[2]
            # Get only the results of specific groups
            A <- filtered_data %>% filter(strategy == groupA)
            B <- filtered_data %>% filter(strategy == groupB)
            # Compute VDA
            vda_result <- VD.A(A$RC, B$RC)
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

# Remove Gin without RTS from the graph
prunedRC <- addedRC %>%
  filter(strategy != "GI")

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
  facet_wrap(~program, nrow = 2)
print(plotRC)

# Save plot
dir.create("images", showWarnings = FALSE)
ggsave("images/RC.png",
       plot = plotRC,
       device = png(),
       width = 700,
       height = 400,
       units = "mm",
       scale = 0.6,
       dpi = "retina")

treatedMetrics <- metrics %>%
  group_by(algorithm, program, strategy) %>%
  summarise(`Overhead` = as.numeric(dmilliseconds(sum(profile_time)), "hours"),
            `Optimisation Time` = as.numeric(dmilliseconds(sum(exec_time)), "hours"))

treatedMetrics <- treatedMetrics %>%
  pivot_longer(c(`Overhead`, `Optimisation Time`), names_to = "Phase", values_to = "Value") %>%
  mutate(Phase = paste(str_to_upper(algorithm), Phase, sep = " ")) %>%
  filter(!str_starts(Phase, "GP Overhead")) %>%
  mutate(Phase = str_replace(Phase, "LS Overhead", "Overhead"))

# Modify names of strategies
treatedMetrics <- treatedMetrics %>%
  mutate(strategy = str_replace(strategy, "ekstazi", "GI+Ekstazi"))
treatedMetrics <- treatedMetrics %>%
  mutate(strategy = str_replace(strategy, "none", "GI"))
treatedMetrics <- treatedMetrics %>%
  mutate(strategy = str_replace(strategy, "starts", "GI+STARTS"))
treatedMetrics <- treatedMetrics %>%
  mutate(strategy = str_replace(strategy, "random", "GI+Random"))
# Order
treatedMetrics$strategy <- factor(treatedMetrics$strategy, levels = c("GI", "GI+Ekstazi", "GI+STARTS", "GI+Random"))

treatedMetrics <- treatedMetrics %>%
  group_by(program, strategy) %>%
  mutate(total = sum(`Value`))

plotCumulative <- treatedMetrics %>%
  ggplot(aes(x=strategy)) + 
  geom_bar(aes(fill=Phase, y=Value), position="stack", stat="identity") +
  geom_text(aes(label=round(total, digits = 2), y=total), position=position_dodge(width=0.9), vjust=-0.25) +
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
  scale_y_continuous(breaks = seq(0,150,10)) +
  scale_fill_brewer(palette="Spectral", aesthetics = "colour") +
  facet_wrap(~program, nrow = 2)

# Save plot
ggsave("images/cumulative.png",
       plot = plotCumulative,
       device = png(),
       width = 700,
       height = 400,
       units = "mm",
       scale = 0.6,
       dpi = "retina")

treatedMetrics %>%
  group_by(strategy) %>%
  summarise(total_weeks = dhours(sum(total)), duration = as.numeric(dhours(sum(total)), "hous"))
