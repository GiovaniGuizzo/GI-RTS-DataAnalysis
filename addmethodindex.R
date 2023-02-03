library(tidyverse)
library(GGally)

# Algorithms used in the experiments
algorithms <- c("gp", "ls")

# Programs used in the experiment
programs <- c("commons-codec",
              "commons-compress",
              "commons-fileupload",
              "commons-imaging")

# Strategies used in the experiment
strategies <- c("none", "ekstazi", "starts", "random")

for (algorithm in algorithms) {
  for (program in programs) {
    for (strategy in strategies) {
      for (run in 1:20) {
        profileFolderPath <- file.path("profile-results", 
                                       program, 
                                       strategy)
        realIndex <- profileFolderPath %>%
          file.path("profile_result_1.csv") %>%
          read_csv(col_types = cols(
            Project = col_character(),
            MethodIndex = col_double(),
            Method = col_character(),
            Count = col_double(),
            Tests = col_character()
          ))
        
        realIndex <- realIndex$MethodIndex[1]
        
        folderPath <- file.path(paste(algorithm, "results", sep = "-"), 
                                program, 
                                strategy)
        
        csvPath <- folderPath %>%
          file.path(paste(algorithm, "_result_", run, ".csv", sep = ""))
        
        csvResults <- csvPath %>%
          read_csv(col_types = cols(
            MethodName = col_character(),
            Patch = col_character(),
            Compiled = col_logical(),
            AllTestsPassed = col_logical(),
            `TotalExecutionTime(ms)` = col_double(),
            Fitness = col_double(),
            FitnessImprovement = col_double(),
            TimeStamp = col_double()
          )) %>%
          mutate(MethodIndex = realIndex)
        
        write_csv(csvResults, csvPath)
      }
    }
  }
}