#krishna Murthy
#creation of Table : summary table of treatment-emergent adverse events (TEAEs)
#program name: 01_create_ae_summary_table.R


dir.create("logs", showWarnings = FALSE)

log_file <- "logs/01_create_ae_summary_table_program.log"

sink(log_file)

cat(" AE summary table  program started\n")

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



# FDA TABLE 10 ---------------------------------------------------------------
# Self-contained template. Run top-to-bottom to build the objects below.
# Used by tests via run_template() and published in the catalog.

library(dplyr)
library(gtsummary)

adsl <- pharmaverseadam::adsl
adae <- pharmaverseadam::adae

# Pre-processing --------------------------------------------
adae1 <- adae |>
  filter(
    # safety population
    SAFFL == "Y",
    # serious adverse events
    TRTEMFL == "Y"
  )

adae2 <- adae |>
  filter(
    # safety population
    SAFFL == "Y",
    # serious adverse events
    TRTEMFL == "Y"
  ) |> 
  mutate(TRT01A="Total",
         TRT01AN=999)

adae <- bind_rows(adae1,adae2)

adsl1 <- adsl
adsl2 <- adsl |> 
  mutate(TRT01A="Total",
         TRT01AN=999) 

adsl <- bind_rows(adsl1,adsl2)

tbl <- adae |>
  tbl_hierarchical(
    variables = c(AESOC, AEDECOD),
    by = TRT01A,
    id = USUBJID,
    denominator = adsl,
    overall_row = TRUE,
    label = "..ard_hierarchical_overall.." ~ "Treatment Emergent AEs"
  )

tbl

ard <- gtsummary::gather_ard(tbl)

withr::local_options(width = 9999, tibble.print_max = Inf)
print(ard, columns = "all")


gtsummary::as_gt(tbl) %>%
  gt::gtsave("AE_Table_Summary.html")


cat("AE summary Table program completed\n")

sink()


