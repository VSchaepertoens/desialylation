library(tidyverse)
library(janitor)



# Load data ---------------------------------------------------------------

mofi_results_sial <-
  read_csv(
    "data/041022_Hits MoFi N-glycan + acetylated 1%Cutoff_cysteinyl.csv",
    skip = 102
  ) %>%
  clean_names() %>%
  separate(
    id,
    into = c("peak_id", "hit_id", "perm_id"),
    sep = "-",
    remove = FALSE
  )

mofi_results_desial <-
  read_csv(
    "data/111022_Hits MoFi N-glycan + acetylated 1%Cutoff_cysteinyl_Sialidase.csv",
    skip = 72
  ) %>%
  clean_names() %>%
  separate(
    id,
    into = c("peak_id", "hit_id", "perm_id"),
    sep = "-",
    remove = FALSE
  )



# Desialylate -------------------------------------------------------------

# finds the mean of an interval given as string
# e.g. "(101,137]" -> 118
find_interval_mean <- function(interval) {
  span <- str_match(interval, "\\((.*),(.*)\\]")
  (as.numeric(span[,2]) + as.numeric(span[,3])) / 2
}

# calculates desialylated masses, compares and matches them to experimentally
# desialylated peaks, bins them with a given width (mass_tolerance), and sums
# intensity in each bin
desialylate <- function(peaks_sial,
                        peaks_desial_exp,
                        mass_desial_peak,
                        filter_peaks = TRUE,
                        filter_max_peaks = TRUE,
                        mass_tolerance = 5) {
  MASS_NEU5AC <- 291.256  # mass of a sialic acid
  MASS_AC <- 42.0106  # acetylation

  # experimental desialylation
  peaks_desial_comp <-
    peaks_sial %>%
    mutate(mass_desial = exp_mass - MASS_NEU5AC * neu5ac - MASS_AC * acetyl)

  # if desired, only keep MoFi hits where the computationally desialylated mass
  # corresponds to a peak in the experimentally desialylated spectrum (within
  # the given mass_tolerance)
  if (filter_peaks) {
    peaks_desial_comp <-
      map_dfr(
        unique(peaks_desial_exp$exp_mass),
        function(mass) {
          peaks_desial_comp %>%
            filter(abs(mass_desial - mass) < mass_tolerance)
        }
      ) %>%
      distinct(id, .keep_all = TRUE) %>%
      group_by(peak_id) %>%
      mutate(hit_score = hit_score / sum(hit_score) * 100) %>%
      ungroup()
  }

  n_bins <- round(
    (max(peaks_desial_comp$mass_desial) - min(peaks_desial_comp$mass_desial))
    / mass_tolerance
  )

  # if desired, only keep MoFi hits where the computationally desialylated mass
  # corresponds to a mass of a given peak (within the given mass_tolerance)
  if (filter_max_peaks){
    peaks_desial_comp %>%
      mutate(mass = mass_desial %>% cut(n_bins) %>% find_interval_mean()) %>%
      filter(near(mass, mass_desial_peak, tol = mass_tolerance))
  }

  # peaks_desial_comp %>%
  #     mutate(mass = mass_desial %>% cut(n_bins) %>% find_interval_mean()) %>%
  #     group_by(mass) %>%
  #     summarise(intensity = sum(percent)) %>%
  #     mutate(intensity = intensity / max(intensity) * 100)

}

###BUG In order to run this, uncomment lines 93-97
df_desial <- desialylate(mofi_results_sial, mofi_results_desial, filter_max_peaks = FALSE)


## To find the masses of the desired peaks
# sort(df_desial$intensity,decreasing=TRUE,index.return= TRUE)
#
# df_desial[c(108,159,167,83),]#display top 4 peaks and their masses
#
# max_peak2 <- df_desial[c(159),'mass']
# max_peak3 <- df_desial[c(167),'mass']

mass_peak2 = 110633.5
mass_peak3 = 110958.5
df_desial_peak2 <- desialylate(mofi_results_sial, mofi_results_desial,mass_desial_peak = mass_peak2)
df_desial_peak3 <- desialylate(mofi_results_sial, mofi_results_desial,mass_desial_peak = mass_peak3)

# write.table(df_desial_peak2, file = "df_desial_peak2.txt", sep = "\t",
#             row.names = TRUE, col.names = NA)
# write.table(df_desial_peak3, file = "df_desial_peak3.txt", sep = "\t",
#             row.names = TRUE, col.names = NA)


# Plot data ---------------------------------------------------------------

# bind_rows(
#   experimental =
#     mofi_results_desial %>%
#     group_by(peak_id) %>%
#     summarise(across(c(exp_mass, percent), first)) %>%
#     select(mass = exp_mass, intensity = percent),
#   computational =
#     df_desial %>%
#     mutate(intensity = intensity * -1),
#   .id = "desialylation"
# ) %>%
#   mutate(desialylation = fct_rev(desialylation)) %>%
#   ggplot(aes(mass, 0, xend = mass, yend = intensity)) +
#   geom_segment(aes(color = desialylation)) +
#   geom_hline(yintercept = 0) +
#   xlab("mass (Da)") +
#   ylab("relative intensity (%)") +
#   theme_bw() +
#   theme(panel.grid = element_blank())
#ggsave("desialylated_spectra_experimental_vs_computational.pdf")



# df_plot <-
#   bind_rows(
#     experimental =
#       mofi_results_sial %>%
#       group_by(peak_id) %>%
#       summarise(across(c(exp_mass, percent), first)) %>%
#       select(mass = exp_mass, intensity = percent),
#     desialylated =
#       df_desial %>%
#       mutate(intensity = intensity * -1),
#     .id = "spectrum"
#   ) %>%
#   mutate(spectrum = fct_rev(spectrum))
#
# ggplot(df_plot, aes(mass, 0, xend = mass, yend = intensity)) +
#   geom_segment(aes(color = spectrum)) +
#   geom_hline(yintercept = 0) +
#   xlab("mass (Da)") +
#   ylab("relative intensity (%)") +
#   theme_bw() +
#   theme(panel.grid = element_blank())
# ggsave("spectra_original_desialylated_filtered.pdf")
