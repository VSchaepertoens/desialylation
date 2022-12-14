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
# desialylated peaks, (optional) filters hit_scores, (optional) filters peaks of certain mass,
# bins them with a given width (mass_tolerance), and sums intensity in each bin
desialylate <- function(peaks_sial,
                        peaks_desial_exp,
                        mass_desial_peak,
                        filter_peaks = TRUE,
                        filter_hit_score = TRUE,
                        hit_score_cutoff = 0.01,
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
  
  # if desired, only keep MoFi hits where hit_score > hit_score_cutoff
  if (filter_hit_score) {
    peaks_desial_comp <-
      peaks_desial_comp %>%
      filter(hit_score > hit_score_cutoff) %>%
      group_by(peak_id) %>%
      mutate(hit_score = hit_score / sum(hit_score) * 100) %>%
      ungroup()
  }
  
  n_bins <- round(
    (max(peaks_desial_comp$mass_desial) - min(peaks_desial_comp$mass_desial))
    / mass_tolerance
  )
  
  peaks_desial_comp %>%
    mutate(mass = mass_desial %>% cut(n_bins) %>% find_interval_mean())%>%
    group_by(mass) %>%
    summarise(intensity = sum(percent)) %>%
    mutate(intensity = intensity / max(intensity) * 100)
  
}

# 1.The intact glycoform annotations were computationally desialylated to obtain a desialylated in silico spectrum of Myozyme
df_desial <- desialylate(mofi_results_sial, mofi_results_desial, filter_hit_score=FALSE, filter_peaks = FALSE)

# 2. The _in silico_ desialylated masses were filtered based on their correspondence with the experimentally desialylated masses and the fractional abundances of the glycoform annotations were normalized to 100%
df_desial_filtered <- desialylate(mofi_results_sial, mofi_results_desial, filter_hit_score=FALSE, filter_peaks = TRUE)


# Plot data ---------------------------------------------------------------

bind_rows(
  experimental =
    mofi_results_desial %>%
    group_by(peak_id) %>%
    summarise(across(c(exp_mass, percent), first)) %>%
    select(mass = exp_mass, intensity = percent),
  computational =
    df_desial %>%
    mutate(intensity = intensity * -1),
  .id = "desialylation"
) %>%
  mutate(desialylation = fct_rev(desialylation)) %>%
  ggplot(aes(mass, 0, xend = mass, yend = intensity)) +
  geom_segment(aes(color = desialylation)) +
  geom_hline(yintercept = 0) +
  xlab("mass (Da)") +
  ylab("relative intensity (%)") +
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave("plots/desialylated_spectra_experimental_vs_computational_unfiltered.pdf")

bind_rows(
  experimental =
    mofi_results_desial %>%
    group_by(peak_id) %>%
    summarise(across(c(exp_mass, percent), first)) %>%
    select(mass = exp_mass, intensity = percent),
  computational =
    df_desial_filtered %>%
    mutate(intensity = intensity * -1),
  .id = "desialylation"
) %>%
  mutate(desialylation = fct_rev(desialylation)) %>%
  ggplot(aes(mass, 0, xend = mass, yend = intensity)) +
  geom_segment(aes(color = desialylation)) +
  geom_hline(yintercept = 0) +
  xlab("mass (Da)") +
  ylab("relative intensity (%)") +
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave("plots/desialylated_spectra_experimental_vs_computational_filtered.wmf")

