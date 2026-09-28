# ── Packages ────────────────────────────────────────────────
library(ggplot2)
library(patchwork)
library(sf)
library(ggspatial)
library(rnaturalearth)
library(rnaturalearthdata)
library(tidyverse)
library(metR)
library(raster)
library(marmap)

# ── Data loading & prep ─────────────────────────────────────

# Bathymetry (cached after first download)
mar_bathy <- getNOAA.bathy(lon1 = 177, lon2 = 180, lat1 = -20, lat2 = -15)
mar_bathy_rast <- marmap::as.raster(mar_bathy)   # always works, regardless of load order

mar_bathy_rast[mar_bathy_rast >= 0] <- NA
bathy_df <- as.data.frame(mar_bathy_rast, xy = TRUE)
colnames(bathy_df)[3] <- "depth"
bathy_df <- bathy_df[!is.na(bathy_df$depth), ]

contour_levels <- c(-20, -50, -100, -200, -500, -1000, -2000, -3000)
label_levels   <- c(-20, -50, -100, -500, -1000)

# eBird & chumming data
ebird_fiji_final <- read_csv("ebird_fiji_final.csv")
chum             <- read_csv("chumming.csv")

chum_sf <- chum %>%
  filter(Group == "Position 1") %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)

plot_data_sf <- st_as_sf(ebird_fiji_final, coords = c("lon.x", "lat.x"), crs = 4326)

# Site-level summary for geographic diversity panel
site_summary <- survey_processed %>%
  group_by(site_id, lat, lon) %>%
  summarise(species_richness = n_distinct(species_code),
            n_obs            = n(),
            .groups          = "drop") %>%
  filter(!is.na(species_richness)) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

# Fiji coastline
fiji_map <- ne_countries(country = "Fiji", scale = 10, returnclass = "sf")
if (st_bbox(fiji_map)$xmin < 0) fiji_map <- st_shift_longitude(fiji_map)

# ── Shared theme ────────────────────────────────────────────
theme_pub <- theme_minimal() +
  theme(
    axis.text      = element_text(size = 11),
    axis.title     = element_blank(),
    legend.position= "none",
    panel.grid.major = element_line(colour = grey(0.8), linetype = "dashed", linewidth = 0.4),
    panel.background = element_rect(fill = "#e8f4fd"),
    plot.margin    = margin(4, 4, 4, 4)
  )

# ── Panel A: Geographic diversity (Gau inset) ───────────────
geo_diversity_plot <- ggplot() +
  geom_sf(data = fiji_map, fill = "grey85", colour = "grey40", linewidth = 0.3) +
  stat_contour(data = bathy_df, aes(x = x, y = y, z = depth),
               breaks = contour_levels, colour = "grey50", linewidth = 0.4, alpha = 0.6) +
  geom_text_contour(data = bathy_df, aes(x = x, y = y, z = depth),
                    breaks = label_levels,
                    colour = "grey40", size = 2.8, rotate = TRUE) +
  geom_sf(data = site_summary, aes(size = n_obs), colour = "#8dd3c7", alpha = 0.8) +
  geom_sf(data = chum_sf, shape = 17, colour = "#fdb462", size = 4, alpha = 1) +
  coord_sf(crs = 4326, xlim = c(179, 179.6), ylim = c(-18.65, -17.8), expand = FALSE) +
  scale_size_continuous(range = c(4, 10), name = NULL) +
  annotate("text", x = 179.4, y = -17.96, label = "Gau", size = 5.5) +
  annotate("text", x = 179.03, y = -17.83, label = "A", size = 7, fontface = "bold") +
  annotation_scale(location = "bl", style = "ticks", width_hint = 0.4, text_cex = 1.2, line_width = 1.5) +
  theme_pub +
  theme(panel.border = element_rect(colour = "#fb8072", linewidth = 1.2, fill = NA),
        axis.text.y = element_text(angle = 90))

# ── Panel B: eBird observations ─────────────────────────────
ebird_plot <- ggplot() +
  geom_sf(data = fiji_map, fill = "grey85", colour = "grey40", linewidth = 0.3) +
  stat_contour(data = bathy_df, aes(x = x, y = y, z = depth),
               breaks = contour_levels, colour = "grey50", linewidth = 0.4, alpha = 0.6) +
  geom_text_contour(data = bathy_df, aes(x = x, y = y, z = depth),
                    breaks = label_levels,
                    colour = "grey40", size = 2.8, rotate = TRUE) +
  geom_sf(data = plot_data_sf, colour = "#bc80bd", alpha = 0.8, size = 2.5) +
  coord_sf(xlim = c(178.1, 180), ylim = c(-18.9, -17), expand = FALSE) +
  annotate("text", x = 179.5, y = -17.96, label = "Gau", size = 5.5) +
  annotate("text", x = 178.2, y = -17.08, label = "B", size = 7, fontface = "bold") +
  annotate("rect", xmin = 179, xmax = 179.6, ymin = -18.65, ymax = -17.8,
           fill = NA, colour = "#fb8072", linewidth = 0.8) +
  annotation_scale(location = "bl", style = "ticks", width_hint = 0.4, text_cex = 1.2, line_width = 1.5) +
  theme_pub +
  theme(panel.border = element_rect(colour = "#bc80bd", linewidth = 1.2, fill = NA),
        axis.text.y = element_text(angle = 90))

# ── Panel C: Fiji overview ──────────────────────────────────
fiji_plot <- ggplot(fiji_map) +
  geom_sf(fill = "grey90", colour = "grey30") +
  geom_sf(data = st_graticule(fiji_map, lon = seq(177, 181, by = 1),
                              lat = seq(-19, -15, by = 1)),
          colour = grey(0.8), linewidth = 0.3, linetype = "dashed") +
  coord_sf(xlim = c(176.8, 181), ylim = c(-19.4, -15.9),
           datum = st_crs(fiji_map), expand = FALSE) +
  scale_x_continuous(breaks = seq(177, 181, by = 1),
                     labels = function(x) paste0(x, "\u00B0E")) +
  scale_y_continuous(breaks = seq(-19, -15, by = 1),
                     labels = function(y) paste0(abs(y), "\u00B0S")) +
  annotate("text", x = 178.5, y = -16.2, label = "Fiji", colour = "grey22", size = 7) +
  annotate("text", x = 178.4, y = -18,    label = "Suva",   size = 4.5) +
  annotate("text", x = 177.2, y = -17.75, label = "Nadi",   size = 4.5) +
  annotate("text", x = 179.2, y = -16.55, label = "Levuka", size = 4.5) +
  annotate("text", x = 179.5, y = -17.96, label = "Gau",    size = 4.5) +
  annotate("rect", xmin = 179, xmax = 179.6, ymin = -18.65, ymax = -17.8,
           fill = NA, colour = "#fb8072", linewidth = 0.8) +
  annotate("rect", xmin = 178, xmax = 180, ymin = -19, ymax = -17,
           fill = NA, colour = "#bc80bd", linewidth = 0.8) +
  annotation_scale(location = "bl", style = "ticks", width_hint = 0.4, text_cex = 1.2, line_width = 1.5) +
  theme_pub +
  theme(panel.background = element_rect(fill = "aliceblue"),
        axis.text.y = element_text(angle = 90))

# ── Panel D: Pacific context ────────────────────────────────
pacific_plot <- ggplot() +
  annotation_borders(data = ne_countries(scale = 50, returnclass = "sf"),
                     fill = "grey90", colour = "grey30", linewidth = 0.1) +
  annotate("text", x = 135, y = -25,  label = "Australia", colour = "grey30", size = 5) +
  annotate("text", x = 194, y = -40,  label = "Aotearoa",  colour = "grey30", size = 5) +
  annotate("text", x = 194, y = -16.5,label = "Fiji",      colour = "grey30", size = 6.5, fontface = "bold") +
  coord_cartesian(xlim = c(110, 210), ylim = c(-45, -10)) +
  annotate("rect", xmin = 176, xmax = 181, ymin = -20, ymax = -15,
           fill = NA, colour = "#fb8072", linewidth = 0.8) +
  theme_pub +
  theme(panel.background = element_rect(fill = "aliceblue"),
        axis.text.y = element_text(angle = 90))

# ── Combine & save ──────────────────────────────────────────
final_fig <-
  fiji_plot +
  pacific_plot +
  geo_diversity_plot +
  ebird_plot +
  plot_layout(ncol = 2, nrow = 2,
              widths  = c(1.2, 1.6),
              heights = c(.5, 1))

ggsave("seabird_panel_figure.png", plot = final_fig,
       width = 10, height = 9, dpi = 300, units = "in")
