# Load required packages
library(tidyverse)
library(readr)
library(lubridate)
library(sf)

# Import entire eBird data from tab-separated text file
ebird_data <- read_tsv(
  "ebd_FJ_smp_relApr-2026.txt",
  col_types = cols(),
  na = c("", "NA", "null", "NULL", "NA", "N/A", "#N/A", "null", "NA", "Not Applicable")
)

# Define the geographic bounds (same as your image)
xmin <- 178
xmax <- 180
ymin <- -19
ymax <- -17

# Create a new dataframe with the columns we need
ebird_selected <- ebird_data %>%
  select(
    `COMMON NAME`,
    `SCIENTIFIC NAME`,
    `SUBSPECIES COMMON NAME`,
    `SUBSPECIES SCIENTIFIC NAME`,
    `OBSERVATION COUNT`,
    `BREEDING CODE`,
    `LOCALITY`,
    `LATITUDE`,
    `LONGITUDE`,
    `OBSERVATION DATE`,
    `TIME OBSERVATIONS STARTED`,
    `OBSERVER ID`
  )

# Clean the OBSERVATION COUNT column before converting to numeric
# First, convert to character and remove any non-numeric characters
ebird_selected <- ebird_selected %>%
  mutate(
    `OBSERVATION COUNT` = as.character(`OBSERVATION COUNT`),
    `OBSERVATION COUNT` = str_remove_all(`OBSERVATION COUNT`, "[^0-9.]"),  # Remove non-numeric except decimal points
    `OBSERVATION COUNT` = str_remove_all(`OBSERVATION COUNT`, "^\\.+"),  # Remove leading decimal points
    `OBSERVATION COUNT` = str_remove_all(`OBSERVATION COUNT`, "\\.$"),   # Remove trailing decimal points
    `OBSERVATION COUNT` = if_else(`OBSERVATION COUNT` == "", NA_character_, `OBSERVATION COUNT`)  # Replace empty strings with NA
  ) %>%
  # Convert to numeric - this should now work without warnings
  mutate(`OBSERVATION COUNT` = as.numeric(`OBSERVATION COUNT`))

# Check how many values were converted to NA
cat("=== OBSERVATION COUNT Conversion Summary ===\n")
cat("Total observations:", nrow(ebird_selected), "\n")
cat("NA values in OBSERVATION COUNT:", sum(is.na(ebird_selected$`OBSERVATION COUNT`)), "\n")

# Clean and filter the selected data
ebird_clean <- ebird_selected %>%
  # Remove rows with missing essential information
  filter(!is.na(`OBSERVATION DATE`) & 
           !is.na(`COMMON NAME`) & 
           !is.na(`LATITUDE`) & 
           !is.na(`LONGITUDE`) &
           !is.na(`OBSERVATION COUNT`)) %>%
  
  # Convert date and time to proper formats
  mutate(
    observation_datetime = paste(`OBSERVATION DATE`, `TIME OBSERVATIONS STARTED`),
    year = year(`OBSERVATION DATE`),
    month = month(`OBSERVATION DATE`),
    day = day(`OBSERVATION DATE`),
    hour = hour(`TIME OBSERVATIONS STARTED`),
    minute = minute(`TIME OBSERVATIONS STARTED`)
  ) %>%
  
  # Clean up names (remove extra spaces, etc.)
  mutate(`COMMON NAME` = trimws(`COMMON NAME`),
         `SCIENTIFIC NAME` = trimws(`SCIENTIFIC NAME`),
         `LOCALITY` = trimws(`LOCALITY`),
         `OBSERVER ID` = trimws(`OBSERVER ID`)) %>%
  
  # Remove observations with zero counts
  filter(`OBSERVATION COUNT` > 0) %>%
  
  # Remove duplicate observations (based on location, date, species)
  distinct(`LOCALITY`, `OBSERVATION DATE`, `COMMON NAME`, .keep_all = TRUE) %>%
  # Rename columns to tidy format (lowercase, underscore separated)
  rename(
    common_name = `COMMON NAME`,
    scientific_name = `SCIENTIFIC NAME`,
    subspecies_common_name = `SUBSPECIES COMMON NAME`,
    subspecies_scientific_name = `SUBSPECIES SCIENTIFIC NAME`,
    observation_count = `OBSERVATION COUNT`,
    breeding_code = `BREEDING CODE`,
    locality = `LOCALITY`,
    lat = `LATITUDE`,
    lon = `LONGITUDE`,
    observation_date = `OBSERVATION DATE`,
    time_observations_started = `TIME OBSERVATIONS STARTED`,
    observer_id = `OBSERVER ID`
  )

# Subset by geographic bounds FIRST
cat("=== Subsetting by Geographic Bounds ===\n")
ebird_bounded <- ebird_clean %>%
  filter(between(lon, xmin, xmax),
         between(lat, ymin, ymax))

cat("Observations after geographic filtering:", nrow(ebird_bounded), "\n")

# Convert to sf for spatial operations
ebird_sf <- ebird_bounded %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

# Exclude terrestrial observations but KEEP coordinates as regular columns
ebird_final <- ebird_sf %>%
  # Extract coordinates and convert back to regular dataframe
  st_coordinates() %>%
  as.data.frame() %>%
  # Rename coordinate columns properly
  setNames(c("lon", "lat")) %>%
  # Add row numbers to match original data
  mutate(row_num = 1:n()) %>%
  # Join with original data (excluding the geometry column)
  left_join(ebird_bounded %>% st_drop_geometry() %>% mutate(row_num = 1:n()), 
            by = "row_num") %>%
  select(-row_num) %>%  # Remove temporary row number column
  # Remove duplicate coordinate columns (keep only the ones from st_coordinates)
  select(-lat.y, -lon.y) %>%
  # Filter out unwanted bird species
  filter(!str_detect(common_name, 
                     regex("sandpiper|/|curlew|jaeger|koel|falcon|gull|ruff|knot|crake|sanderling|magpie|lori|Brown/Black|sp.|plover|whimbrel|egret|godwit|turnstone|heron|kingfisher|owl|cuckoo|hawk|robin|warbler|shrikebill|sparrow|tattler|fowl|thicket|fantail|parrot|dove|myna|starling|triller|avadavat|monarch|harrier|flycatcher|honeyeater|pigeon|swiftlet|Whistler|duck|swallow|myzomela|eye|lory|bulbul|woodswallow|thrush|lapwing|rail", 
                           ignore_case = TRUE))&
  !str_detect(observer_id, "obsr793933"))

cat("Observations after exclusion:", nrow(ebird_final), "\n")
cat("observations removed:", nrow(ebird_bounded) - nrow(ebird_final), "\n")

# Save final filtered dataframe with all variables
write_csv(ebird_final, "ebird_fiji_final_allspecies.csv")


# Summary of final filtered data
cat("\n=== Final Filtered Data Summary ===\n")
cat("Original rows:", nrow(ebird_data), "\n")
cat("Selected rows:", nrow(ebird_selected), "\n")
cat("Cleaned rows:", nrow(ebird_clean), "\n")
cat("After geographic filtering:", nrow(ebird_bounded), "\n")
cat("Final filtered rows:", nrow(ebird_final), "\n")
cat("Species found:", length(unique(ebird_final$common_name)), "\n")
cat("Observers:", length(unique(ebird_final$observer_id)), "\n")

# Basic exploration on final filtered data
cat("\n=== Most Common Species (Final Data) ===\n")
species_counts <- ebird_final %>%
  group_by(common_name, scientific_name) %>%
  summarize(
    Total_Observations = n(),
    Individuals = sum(observation_count, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(Total_Observations))

print(species_counts %>% head(10))
# Save species summary with common name and scientific name
write_csv(species_counts, "ebird_species_list.csv")

# Save species summary
write_csv(species_counts, "ebird_species_counts_final.csv")

# Geographic coverage of final filtered data
cat("\n=== Geographic Coverage (Final Data) ===\n")
geo_summary <- ebird_final %>%
  summarize(
    min_lat = min(lat.x),
    max_lat = max(lat.x),
    min_lon = min(lon.x),
    max_lon = max(lon.x),
    unique_localities = n_distinct(locality),
    .groups = "drop"
  )

print(geo_summary)


# Load required packages
library(ggplot2)
library(sf)
library(rnaturalearth)

# Get Fiji map data
fiji_map <- ne_countries(country = "Fiji", scale = 50, returnclass = "sf")

# Convert ebird data to sf object for plotting
plot_data_sf <- st_as_sf(
  ebird_final,
  coords = c("lon.x", "lat.x"),
  crs = 4326
)

# Create the plot using ggplot2
ebird_plot <- ggplot() +
  # Add Fiji map as background
  geom_sf(data = fiji_map, fill = "lightgray", color = "gray50") +
  # Add eBird observations as points
  geom_sf(data = plot_data_sf, 
          color = "#698E7CFF", 
          size = 2,
          alpha = 0.7) +
  # Set coordinate limits
  coord_sf(xlim = c(xmin, xmax), 
           ylim = c(ymin, ymax),
           expand = FALSE) +
  # Add title and labels
  labs(title = "eBird Observations - Final Filtered Area",
       x = "Longitude",
       y = "Latitude") +
  # Clean up the theme
  theme_minimal() +
  theme(panel.grid.major = element_line(color = "lightgray", linetype = "dotted"),
        panel.background = element_rect(fill = "white"),
        plot.title = element_text(hjust = 0.5, face = "bold"))

# Print the plot
print(ebird_plot)

# Save the plot as a PNG file
ggsave("ebird_observations_plot_final.png", 
       plot = plot, 
       width = 8, 
       height = 6, 
       dpi = 300,
       units = "in")


ebird_all_seabirds <- ebird_clean %>%
  # Filter out unwanted bird species
  filter(!str_detect(common_name, 
                     regex("silk|swamp|wax|munia|mallard|sandpiper|/|curlew|jaeger|koel|falcon|gull|ruff|knot|crake|sanderling|magpie|lori|Brown/Black|sp.|plover|whimbrel|egret|godwit|turnstone|heron|kingfisher|owl|cuckoo|hawk|robin|warbler|shrikebill|sparrow|tattler|fowl|thicket|fantail|parrot|dove|myna|starling|triller|avadavat|monarch|harrier|flycatcher|honeyeater|pigeon|swiftlet|Whistler|duck|swallow|myzomela|eye|lory|bulbul|woodswallow|thrush|lapwing|rail", 
                           ignore_case = TRUE))&
           !str_detect(observer_id, "obsr793933"))

unique(ebird_all_seabirds$common_name)
cat("Observations after exclusion:", nrow(ebird_final), "\n")
cat("observations removed:", nrow(ebird_bounded) - nrow(ebird_final), "\n")
