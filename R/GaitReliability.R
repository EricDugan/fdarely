# -------------------------
# LOAD LIBRARIES
# -------------------------
library(here)
library(tidyverse)
library(tcltk)


# -------------------------
# LOAD THE DATA
# -------------------------

# Set the default search path using the here() function
default_path <- here("GaitData")

# Allow user to select multiple data files using a file picker with the default search path
data_paths <- tk_choose.files(default = default_path, caption = "Select data files", multi = TRUE)

# Read in the datasets from the selected paths
data_list <- lapply(data_paths, function(path) {
  data = read.csv(path, sep="\t", header=FALSE)
  
  # Dynamically determine the number of columns
  num_columns <- ncol(data)
  
  # Combine the 2nd and 5th rows of columns to create new column names
  new_col_names <- paste(data[2, 2:num_columns], data[5, 2:num_columns], sep="_")

  # Apply the function to the new column names
  new_col_names <- make_unique(new_col_names)
  
  # Set the new column names
  colnames(data)[2:num_columns] <- new_col_names
  
  # Remove the rows that were used to create the column names
  data <- data[-c(1:5), ]
  
  # Reset row names after removal
  rownames(data) <- NULL
  
  # Extract the first column as new column names
  new_headers <- data[, 1]
  
  # Drop the first column and transpose the rest
  data_transposed <- as.data.frame(t(data[,-1]))
  
  # Set the new column names
  colnames(data_transposed) <- new_headers
  
  # Create a new column "Variable" with current row names
  data_transposed$Variable <- rownames(data_transposed)
  
  # Reset the row names
  rownames(data_transposed) <- NULL
  
  # Reorder columns to ensure "Variable" is the first column
  data_transposed <- data_transposed %>%
    select(Variable, everything())
  
  # Extract file name without extension
  file_name <- tools::file_path_sans_ext(basename(path))
  
  # Add the file name as the first column
  data_transposed$FileName <- file_name
  
  # Reorder columns to ensure "FileName" and "Variable" are the first columns
  data_transposed <- data_transposed %>%
    select(FileName, Variable, everything())
  
  # Split the FileName column into multiple columns
  data_transposed <- data_transposed %>%
    separate(FileName, into = c("ControlSubject", "PT", "DataSet", "Session", "Type", "Extension"), sep = "_")
  
  # Drop the 'Extension' column as it's not needed
  data_transposed$Extension <- NULL

  return(data_transposed)
})


# Check the first few rows of the first dataset (as an example)
head(data_list[[1]])


# Combine all data frames in data_list into one data frame
combined_data <- bind_rows(data_list)

# View the first few rows of the combined data frame
head(combined_data)



# Filter rows that match the desired variables
filtered_data <- combined_data %>%
  filter(grepl("^LANKLE_ANGLE_X_", Variable))

# Determine the start and end columns for numeric conversion
start_col <- which(colnames(filtered_data) == "1")
end_col <- start_col + 100

# Convert the appropriate columns to numeric
filtered_data[, start_col:end_col] <- lapply(filtered_data[, start_col:end_col], as.numeric)

# Group by PT and Session and calculate the average for each of the columns 1 through 101
average_data <- filtered_data %>%
  ungroup() %>%
  group_by(PT, Session) %>%
  summarise(across(`1`:`101`, mean, .names = "avg_{.col}"), .groups = 'drop')

long_data <- average_data %>%
  pivot_longer(cols = starts_with("avg_"),
               names_to = "Timepoint",
               values_to = "Value")


ggplot(long_data, aes(x = as.numeric(str_replace(Timepoint, "avg_", "")), y = Value, color = interaction(PT, Session), group = interaction(PT, Session))) +
  geom_line() +
  labs(title = "Average Data per PT per Session",
       x = "% Gait Cycle",
       y = "Degrees",
       color = "Group (PT_Session)") +
  theme_minimal()


ggplot(long_data, aes(x = as.numeric(str_replace(Timepoint, "avg_", "")), 
                      y = Value, 
                      color = PT, 
                      linetype = Session, 
                      group = interaction(PT, Session))) +
  geom_line() +
  labs(title = "Average Data per PT per Session",
       x = "% Gait Cycle",
       y = "Degrees") +
  scale_color_manual(values = c("blue", "red")) + # Modify colors if needed
  theme_minimal()


average_data_num = average_data

# Create a mapping of characters to numeric codes
mapping <- c(AS = 1, SH = 2)
mapping_Session

# Convert PT values using the mapping
average_data_num$PT <- mapping[average_data_num$PT]





library("writexl")
write_xlsx(filtered_data,"C:\\Users\\eldugan\\Desktop\\GitHub Repos\\fdarely\\GaitData\\LeftAnkle.xlsx")