# ── Packages ────────────────────────────────────────────────

library(lme4)
library(unmarked)
library(vegan)
library(tidyverse)
library(readxl)
library(lubridate)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
# Right after your library block:
select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate
summarise <- dplyr::summarise
group_by <- dplyr::group_by
left_join <- dplyr::left_join
distinct <- dplyr::distinct
n_distinct <- dplyr::n_distinct
arrange <- dplyr::arrange
rowid_to_column <- tibble::rowid_to_column
pivot_longer <- tidyr::pivot_longer
pivot_wider <- tidyr::pivot_wider
# ── Read data ───────────────────────────────────────────────
survey_data <- read.csv("survey_data.csv")
all_spp_survey <- read.csv("species_list.csv")
ebird <- read.csv("ebird_fiji_final.csv")

# ── Process survey data ─────────────────────────────────────
survey_data <- survey_data %>%
  rowid_to_column("site_id") %>%
  mutate(site_id = as.character(site_id))

site_info <- survey_data %>%
  select(site_id, Lat, Long, Date, Time) %>%
  distinct()

survey_processed <- survey_data %>%
  pivot_longer(
    cols = c(COPE, TAPE, KEPE, WTSH, BUSH, SOSH, FFSH, BLND, BRND, RFBO, BRBO, LEFR, COTE, GBTE, WTTE, WFSP, UKN_Petrel, UKN_Shear, UNKN),
    names_to = "species_code",
    values_to = "observation_count",
    values_drop_na = TRUE
  ) %>%
  mutate(
    common_name = case_when(
      species_code == "COPE" ~ "Collared Petrel",
      species_code == "TAPE" ~ "Tahiti Petrel",
      species_code == "KEPE" ~ "Kermadec Petrel",
      species_code == "WTSH" ~ "Wedge-tailed Shearwater",
      species_code == "BUSH" ~ "Buller's Shearwater",
      species_code == "SOSH" ~ "Sooty Shearwater",
      species_code == "FFSH" ~ "Flesh-footed Shearwater",
      species_code == "BLND" ~ "Black Noddy",
      species_code == "BRND" ~ "Brown Noddy",
      species_code == "RFBO" ~ "Red-footed Booby",
      species_code == "BRBO" ~ "Brown Booby",
      species_code == "LEFR" ~ "Lesser Frigatebird",
      species_code == "COTE" ~ "Bridled Tern",
      species_code == "GBTE" ~ "Grey-backed Tern",
      species_code == "WTTE" ~ "White Tern",
      species_code == "WFSP" ~ "White-faced Storm-Petrel",
      species_code == "UKN_Petrel" ~ "Petrel sp",
      species_code == "UKN_Shear" ~ "Shearwater sp",
      species_code == "UNKN" ~ "Unknown sp",
      TRUE ~ NA_character_
    ),
    scientific_name = case_when(
      species_code == "COPE" ~ "Pterodroma brevicaudata",
      species_code == "TAPE" ~ "Pseudobulweria rostrata",
      species_code == "KEPE" ~ "Pterodroma neglecta",
      species_code == "WTSH" ~ "Ardenna pacifica",
      species_code == "BUSH" ~ "Ardenna tenuirostris",
      species_code == "SOSH" ~ "Ardenna grisea",
      species_code == "FFSH" ~ "Puffinus carneipes",
      species_code == "BLND" ~ "Anous minutus",
      species_code == "BRND" ~ "Anous stolidus",
      species_code == "RFBO" ~ "Sula sula",
      species_code == "BRBO" ~ "Sula leucogaster",
      species_code == "LEFR" ~ "Fregata ariel",
      species_code == "COTE" ~ "Sterna anaethetus",
      species_code == "GBTE" ~ "Onychoprion lunatus",
      species_code == "WTTE" ~ "Gygis alba",
      species_code == "WFSP" ~ "Fregetta grallaria",
      species_code == "UKN_Petrel" ~ "Petrel spp.",
      species_code == "UKN_Shear" ~ "Shearwater spp.",
      species_code == "UNKN" ~ "Unknown spp.",
      TRUE ~ NA_character_
    ),
    lat = -Lat,
    lon = Long
  ) %>%
  mutate(
    year = as.integer(format(as.Date(Date, format = "%d/%m/%Y"), "%Y")),
    month = as.integer(format(as.Date(Date, format = "%d/%m/%Y"), "%m")),
    day = as.integer(format(as.Date(Date, format = "%d/%m/%Y"), "%d")),
    hour = as.integer(substr(Time, 1, 2)),
    minute = as.integer(substr(Time, 4, 5))
  ) %>%
  select(-Lat, -Long, -Date, -Time) %>%
  left_join(site_info, by = "site_id") %>%
  select(-ends_with(".x")) %>%
  filter(!common_name %in% c("Petrel sp", "Shearwater sp", "Unknown sp")) %>%
  filter(observation_count > 0)

# ── Species-level metrics ───────────────────────────────────
species_metrics <- survey_processed %>%
  group_by(species_code, common_name, scientific_name) %>%
  summarise(
    total_abundance = sum(observation_count),
    mean_abundance = mean(observation_count),
    max_abundance = max(observation_count),
    n_surveys = n(),
    n_sites = n_distinct(site_id),
    .groups = "drop"
  ) %>%
  mutate(encounter_frequency = if_else(n_surveys > 0, n_sites / n_surveys, 0)) %>%
  arrange(desc(total_abundance))

# ── Site-level richness & diversity ─────────────────────────
site_richness <- survey_processed %>%
  group_by(site_id) %>%
  summarise(
    species_richness = n_distinct(species_code),
    total_abundance = sum(observation_count),
    n_observations = n(),
    .groups = "drop"
  )

pa_data <- survey_processed %>%
  group_by(site_id, species_code) %>%
  summarise(present = any(observation_count > 0), .groups = "drop")

pa_matrix <- pa_data %>%
  pivot_wider(id_cols = site_id, names_from = species_code,
              values_from = present, values_fill = FALSE)

pa_matrix_numeric <- as.matrix(select(pa_matrix, -site_id)) * 1

diversity_metrics <- data.frame(
  site_id = rep(pa_matrix$site_id, 4),
  Lat = rep(site_info$Lat[match(pa_matrix$site_id, site_info$site_id)], 4),
  Long = rep(site_info$Long[match(pa_matrix$site_id, site_info$site_id)], 4),
  shannon = c(diversity(pa_matrix_numeric, index = "shannon")),
  simpson = c(diversity(pa_matrix_numeric, index = "simpson")),
  richness = c(specnumber(pa_matrix_numeric)),
  pielou = c(diversity(pa_matrix_numeric, index = "shannon") / log(specnumber(pa_matrix_numeric)))
)

# ── Overall diversity statistics ────────────────────────────
total_species_richness <- length(unique(survey_processed$species_code))

shannon_mean <- mean(diversity_metrics$shannon, na.rm = TRUE)
shannon_sd   <- sd(diversity_metrics$shannon, na.rm = TRUE)
shannon_se   <- shannon_sd / sqrt(sum(!is.na(diversity_metrics$shannon)))

simpson_mean <- mean(diversity_metrics$simpson, na.rm = TRUE)
simpson_sd   <- sd(diversity_metrics$simpson, na.rm = TRUE)
simpson_se   <- simpson_sd / sqrt(sum(!is.na(diversity_metrics$simpson)))

richness_mean <- mean(diversity_metrics$richness, na.rm = TRUE)
richness_sd   <- sd(diversity_metrics$richness, na.rm = TRUE)
richness_se   <- richness_sd / sqrt(sum(!is.na(diversity_metrics$richness)))

pielou_mean <- mean(diversity_metrics$pielou, na.rm = TRUE)
pielou_sd   <- sd(diversity_metrics$pielou, na.rm = TRUE)
pielou_se   <- pielou_sd / sqrt(sum(!is.na(diversity_metrics$pielou)))

overall_diversity <- list(
  mean_shannon = shannon_mean, shannon_sd = shannon_sd, shannon_se = shannon_se,
  mean_simpson = simpson_mean, simpson_sd = simpson_sd, simpson_se = simpson_se,
  mean_richness = richness_mean, richness_sd = richness_sd, richness_se = richness_se,
  mean_pielou = pielou_mean, pielou_sd = pielou_sd, pielou_se = pielou_se,
  total_species = total_species_richness,
  total_sites = nrow(diversity_metrics)
)

# ── Species accumulation curves ─────────────────────────────
accum_matrix <- survey_processed %>%
  group_by(site_id, species_code) %>%
  summarise(present = any(observation_count > 0), .groups = "drop") %>%
  pivot_wider(id_cols = site_id, names_from = species_code,
              values_from = present, values_fill = FALSE)

accum_matrix_numeric <- as.matrix(select(accum_matrix, -site_id)) * 1

accum_curve       <- specaccum(accum_matrix_numeric, method = "random")
accum_curve_ordered <- specaccum(accum_matrix_numeric, method = "rarefaction")

# ── eBird data prep ─────────────────────────────────────────
ebird_renamed <- ebird %>%
  rename(lat = lat.x, lon = lon.x) %>%
  select(common_name, scientific_name, lat, lon, year, month, day, hour, minute, observation_count, everything())

ebird_species <- unique(na.omit(ebird_renamed[, c("common_name", "scientific_name")]))
all_spp_species <- unique(na.omit(all_spp_survey[, c("common_name", "scientific_name")]))

# ── Publication tables ──────────────────────────────────────
diversity_stats_table <- data.frame(
  Metric = c("Shannon Diversity", "Simpson Diversity", "Species Richness",
             "Pielou's Evenness", "Total Species", "Total Sites"),
  Value = c(round(shannon_mean, 3), round(simpson_mean, 3), round(richness_mean, 1),
            round(pielou_mean, 3), overall_diversity$total_species, overall_diversity$total_sites),
  SD = c(round(shannon_sd, 3), round(simpson_sd, 3), round(richness_sd, 1),
         round(pielou_sd, 3), NA, NA),
  SE = c(round(shannon_se, 3), round(simpson_se, 3), round(richness_se, 2),
         round(pielou_se, 3), NA, NA)
)

species_abundance_table <- species_metrics %>%
  select(species_code, common_name, scientific_name, total_abundance, mean_abundance,
         max_abundance, encounter_frequency, n_sites) %>%
  mutate(
    total_abundance = round(total_abundance, 1),
    mean_abundance = round(mean_abundance, 2),
    max_abundance = round(max_abundance, 1),
    encounter_frequency = round(encounter_frequency, 3),
    n_sites = as.integer(n_sites)
  )

site_richness_table <- site_richness %>%
  arrange(desc(species_richness)) %>%
  mutate(species_richness = as.integer(species_richness),
         total_abundance = round(total_abundance, 1))

# ── Plots ───────────────────────────────────────────────────
species_plot <- ggplot(species_abundance_table,
                       aes(x = reorder(species_code, -total_abundance), y = total_abundance)) +
  geom_bar(stat = "identity", fill = "#80b1d3") +
  labs(title = "(A) Total Abundance", x = "Species", y = "Total Abundance") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
        axis.text.y = element_text(size = 16),
        plot.title = element_text(face = "bold", size = 18),
        axis.title = element_text(size = 16))

status_counts <- all_spp_survey %>%
  count(status, name = "count") %>%
  mutate(status = factor(status, levels = unique(status)))

status_plot <- ggplot(status_counts, aes(x = status, y = count, fill = status)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = count), vjust = 1.5, size = 10) +
  labs(title = "(B) Status", x = NULL, y = "Species Count") +
  scale_fill_brewer(palette = "Set3") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 18),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title = element_text(size = 16))

accumulation_plot <- ggplot() +
  geom_line(data = data.frame(sites = accum_curve$sites, richness = accum_curve$richness),
            aes(x = sites, y = richness, color = "Random"), linewidth = 2) +
  geom_line(data = data.frame(sites = accum_curve_ordered$sites, richness = accum_curve_ordered$richness),
            aes(x = sites, y = richness, color = "Rarefaction"), linewidth = 2, linetype = "dashed") +
  labs(title = "(C) Species Accumulation Curve", x = "Number of Sites",
       y = "Species Richness", color = "Curve Type") +
  scale_color_manual(values = c("Random" = "#80b1d3", "Rarefaction" = "#fb8072")) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 18),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title = element_text(size = 16),
    legend.position = c(0.87, 0.1),
    legend.background = element_rect(fill = "white"),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16)
  )


## eBird observations by month (count of checklist records)
ebird_monthly <- ebird_renamed %>%
  group_by(month) %>%
  summarise(n_records = n(), .groups = "drop") %>%
  arrange(month) %>%
  mutate(month_label = month.abb[month])

ebird_monthly_plot <- ggplot(ebird_monthly, aes(x = factor(month_label, levels = month_label[order(month)]), y = n_records)) +
  geom_col(fill = "#bc80bd", width = 0.6) +
  labs(title = "(D) eBird Records by Month", x = "Month", y = "Checklist Records") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 18),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
        axis.text.y = element_text(size = 16),
        axis.title = element_text(size = 16))


## Combine into 2×2 panel
panel_plot <- wrap_plots(
  plots = list(species_plot, status_plot, accumulation_plot, ebird_monthly_plot),
  ncol = 2, nrow = 2
) +
  plot_layout(widths = c(1, 1), heights = c(1, 1))

ggsave("seabird_community_panel_4panel.png", plot = panel_plot, width = 14, height = 12, dpi = 300)

# ── Save results ────────────────────────────────────────────
write.csv(diversity_stats_table, "diversity_statistics_summary.csv", row.names = FALSE)
write.csv(species_abundance_table, "species_abundance_publication.csv", row.names = FALSE)
write.csv(site_richness_table, "site_richness_results.csv", row.names = FALSE)
write.csv(diversity_metrics, "site_diversity_metrics.csv", row.names = FALSE)
