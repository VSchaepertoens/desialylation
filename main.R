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


# calculates desialylated masses, (optional) compares and matches them to experimentally
# desialylated peaks, (optional) filters peaks of certain mass,(optional) filters hit_scores, 
# (optional) filters fractional abundances, bins them with a given width (mass_tolerance), 
# and sums intensity in each bin
desialylate <- function(peaks_sial,
                        peaks_desial_exp,
                        mass_desial_peak,
                        filter_peaks = TRUE,
                        filter_hit_score = TRUE,
                        filter_fract_abund = TRUE,
                        hit_score_cutoff = 0.01,
                        mass_tolerance = 5) {
  MASS_NEU5AC <- 291.256  # mass of a sialic acid
  MASS_AC <- 42.0106  # acetylation
  
  # experimental desialylation
  peaks_desial_comp <-
    peaks_sial %>%
    mutate(mass_desial = exp_mass - MASS_NEU5AC * neu5ac - MASS_AC * acetyl)
  
  #return(peaks_desial_comp)
  
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
  #return(peaks_desial_comp)
  
  # if desired, only keep MoFi hits where hit_score > hit_score_cutoff
  if (filter_hit_score) {
    peaks_desial_comp <-
      peaks_desial_comp %>%
      filter(hit_score > hit_score_cutoff) %>%
      group_by(peak_id) %>%
      mutate(hit_score = hit_score / sum(hit_score) * 100) %>%
      ungroup()
  }
  
  
  # if desired, calculate fractional abundance and filter 
  if (filter_fract_abund) {
    peaks_desial_comp <-
      peaks_desial_comp %>%
      group_by(peak_id) %>%
      mutate(fractional_abundance = (hit_score * percent) / sum(percent)) %>%
      filter(fractional_abundance > 0.01)
  }
  
  #return(peaks_desial_comp)
  
  n_bins <- round(
    (max(peaks_desial_comp$mass_desial) - min(peaks_desial_comp$mass_desial))
    / mass_tolerance
  )
  
  peaks_desial_comp %>%
    mutate(mass = mass_desial %>% cut(n_bins) %>% find_interval_mean()) %>%
    group_by(mass) %>%
    summarise(intensity = sum(percent)) %>%
    mutate(intensity = intensity / sum(intensity) * 100)
  
}

# 1.The intact glycoform annotations were computationally desialylated to obtain
# a desialylated in silico spectrum of Myozyme
df_desial <- desialylate(
  mofi_results_sial,
  mofi_results_desial,
  filter_hit_score = FALSE,
  filter_peaks = FALSE,
  filter_fract_abund = FALSE
)

#write_csv(df_desial_df, "analysis/desialylated_spectra_experimental_vs_computational_unfiltered.csv")

# 2. The _in silico_ desialylated masses were filtered based on their
# correspondence with the experimentally desialylated masses and the fractional
# abundances of the glycoform annotations were normalized to 100%
df_desial_filtered <- desialylate(
  mofi_results_sial,
  mofi_results_desial,
  filter_hit_score = FALSE,
  filter_peaks = TRUE,
  filter_fract_abund = FALSE
)

#write_csv(df_desial_filtered_df, "analysis/desialylated_spectra_experimental_vs_computational_filtered.csv")


# 3. The _in silico_ desialylated hit score were filtered with a
# cut-off (range from 0.01 to 1) and afterwards the fractional abundances of the
# glycoform annotations were normalized to 100%

# df_desial_filtered_cutoff <- desialylate(
#   mofi_results_sial,
#   mofi_results_desial,
#   filter_hit_score = TRUE,
#   filter_peaks = TRUE,
#   hit_score_cutoff = 0.01
# )


# 4. The fractional abundances were calculated for each peak
# and then filtered with a cutoff & normalized to 100%
df_desial_filtered_fract_abund <- desialylate(
  mofi_results_sial,
  mofi_results_desial,
  filter_hit_score = FALSE,
  filter_peaks = TRUE,
  filter_fract_abund = TRUE
)

#write_csv(df_desial_filtered_fract_abund_df, "analysis/desialylated_spectra_experimental_vs_computational_filtered_frac_abund_cutoff_0.01.csv")


# Plot data ---------------------------------------------------------------

plot_spectrum <- function(computational_data) {
  bind_rows(
    experimental =
      mofi_results_desial %>%
      group_by(peak_id) %>%
      summarise(across(c(exp_mass, percent), first)) %>%
      select(mass = exp_mass, intensity = percent) %>%
      mutate(intensity = intensity / sum(intensity) * 100),
    computational =
      computational_data %>%
      mutate(intensity = intensity * -1),
    .id = "desialylation"
  ) %>%
    mutate(desialylation = fct_rev(desialylation)) %>%
    ggplot(aes(mass, 0, xend = mass, yend = intensity)) +
    geom_segment(aes(color = desialylation)) +
    geom_hline(yintercept = 0) +
    scale_y_continuous(n.breaks = 10) +
    xlab("mass (Da)") +
    ylab("fractional abundance (%)") +
    theme_bw() +
    theme(panel.grid = element_blank())

}


plot_spectrum(df_desial)
ggsave("plots/revisions/desialylated_spectra_experimental_vs_computational_unfiltered.wmf")

plot_spectrum(df_desial_filtered)
ggsave("plots/revisions/desialylated_spectra_experimental_vs_computational_filtered.pdf")

#plot_spectrum(df_desial_filtered_cutoff)
#ggsave("plots/desialylated_spectra_experimental_vs_computational_filtered_cutoff0.01.wmf")

plot_spectrum(df_desial_filtered_fract_abund)
ggsave("plots/revisions/desialylated_spectra_experimental_vs_computational_filtered_frac_abund_cutoff_0.01.pdf")

