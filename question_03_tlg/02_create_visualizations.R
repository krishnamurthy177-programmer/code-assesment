#krishna Murthy
#creation of figures : AE severity distribution by treatment (bar chart or heatmap), 
#                      Top 10 most frequent AEs (with 95% CI for incidence rates).
#program name: 01_create_visualizations.R


dir.create("logs", showWarnings = FALSE)

log_file <- "logs/01_create_visualizations.log"

sink(log_file)

cat(" 01_create_visualizations program started\n")

#install packages
install.packages(c("admiral", "sdtm.oak", "gt", "ggplot2"))
install.packages("pharmaversesdtm")
install.packages("pharmaverseadam")
install.packages("dplyr")
install.packages("lubridate")
install.packages("stringr")
install.packages("tidyverse")
install.packages("gtsummary")

library(admiral)
library(dplyr, warn.conflicts = FALSE)
library(pharmaversesdtm)
library(pharmaverseadam)
library(lubridate)
library(stringr)
library(tidyverse)
library(gt)
library(dplyr)
library(gtsummary)

adsl <- pharmaverseadam::adsl
adae <- pharmaverseadam::adae
#########################
## Plot 1: AE severity distribution by treatment.
# Pre-processing --------------------------------------------
adae <- adae |>
  filter(
    # safety population
    SAFFL == "Y",
    # serious adverse events
    TRTEMFL == "Y"
  )

#Count AESEV
ae_data <- adae |> 
  group_by(TRT01A) |> 
  count(AESEV,name="Count") |> 
  mutate(Severity=AESEV,
         Treatment=TRT01A)

ae_data$Treatment <- factor(
  ae_data$Treatment,
  levels = c(
    "Placebo",
    "Xanomeline High Dose",
    "Xanomeline Low Dose"
  )
)

ae_data$Severity <- factor(
  ae_data$Severity,
  levels = c("SEVERE", "MODERATE", "MILD")
)

#  Plot with ggplot2

p <- ggplot(
  ae_data,
  aes(x = Treatment, y = Count, fill = Severity)
) +
  geom_bar(
    stat = "identity",
    width = 0.75
  ) +
  scale_fill_manual(
    values = c(
      "MILD" = "#F8766D",
      "MODERATE" = "#00BA38",
      "SEVERE" = "#619CFF"
    )
  ) +
  labs(
    title = "AE severity distribution by treatment",
    x = "Treatment Arm",
    y = "Count of AEs",
    fill = "Severity/Intensity"
  ) +
  theme_gray(base_size = 12) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5
    )
  )

p

#save png file
ggsave(
  filename = "AE_severity_distribution.png",
  plot = p,
  width = 10,
  height = 6,
  units = "in",
  dpi = 300
)

#########################
## Plot 2: Top 10 most frequent AEs (with 95% CI for incidence rates).
# Count frequency of AEs
ae_counts <- adae %>%
  group_by(AETERM) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(desc(n))

# Select top 10 most frequent AEs
top10_ae <- ae_counts %>%
  slice_head(n = 10)

# Calculate incidence rates and 95% CI
# Assuming denominator = total number of patients in ADAE
n_total <- length(unique(adae$USUBJID))

top10_ae <- top10_ae %>%
  mutate(incidence = n / n_total,
         se = sqrt((incidence * (1 - incidence)) / n_total),
         lower_ci = incidence - 1.96 * se,
         upper_ci = incidence + 1.96 * se)

#  Plot with ggplot2
p2 <-
  ggplot(top10_ae, aes(x = reorder(AETERM, incidence), y = incidence)) +
  geom_point(size = 3, color = "black") +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.2) +
  coord_flip() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(title = "Top 10 Most Frequent Adverse Events",
       subtitle = paste("N =", n_total, "subjects; 95% Clopper-Pearson CIs"),
       x = "Adverse Event Term",
       y = "Percentage of Patients (95% CI)") +
  theme_minimal(base_size = 14)


#save png file
ggsave("Top10_AEs.png", plot = p2, width = 8, height = 6, dpi = 300)


cat(" 01_create_visualizations program completed\n")

sink()


