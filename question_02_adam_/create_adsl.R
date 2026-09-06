#krishna Murthy
#conversion various raw dataset to ADAM.ADSL
#program name: create_adsl.R



dir.create("logs", showWarnings = FALSE)

log_file <- "logs/ADSL_program.log"

sink(log_file)

cat("ADSL program started\n")

#install packages
install.packages(c("admiral", "sdtm.oak", "gt", "ggplot2"))
install.packages("pharmaversesdtm")
install.packages("dplyr")
install.packages("lubridate")
install.packages("stringr")
install.packages("tidyverse")
install.packages("arrow")

library(admiral)
library(dplyr, warn.conflicts = FALSE)
library(pharmaversesdtm)
library(lubridate)
library(stringr)
library(tidyverse)
library(arrow)



dm <- pharmaversesdtm::dm
ds <- pharmaversesdtm::ds
ex <- pharmaversesdtm::ex
ae <- pharmaversesdtm::ae
lb <- pharmaversesdtm::lb
vs <- pharmaversesdtm::vs

dm <- convert_blanks_to_na(dm)
ds <- convert_blanks_to_na(ds)
ex <- convert_blanks_to_na(ex)
ae <- convert_blanks_to_na(ae)
lb <- convert_blanks_to_na(lb)
vs <- convert_blanks_to_na(vs)

adsl <- dm %>%
  select(-DOMAIN)

adsl <- dm %>%
  mutate(TRT01P = ARM, TRT01A = ACTARM)

# Impute start and end time of exposure to first and last respectively,
# Do not impute date
ex_ext <- ex %>%
  derive_vars_dtm(
    dtc = EXSTDTC,
    new_vars_prefix = "EXST"
  ) %>%
  derive_vars_dtm(
    dtc = EXENDTC,
    new_vars_prefix = "EXEN",
    time_imputation = "last"
  )

adsl <- adsl %>%
  derive_vars_merged(
    dataset_add = ex_ext,
    filter_add = (EXDOSE > 0 |
                    (EXDOSE == 0 &
                       str_detect(EXTRT, "PLACEBO"))) & !is.na(EXSTDTM),
    new_vars = exprs(TRTSDTM = EXSTDTM, TRTSTMF = EXSTTMF),
    order = exprs(EXSTDTM, EXSEQ),
    mode = "first",
    by_vars = exprs(STUDYID, USUBJID)
  ) %>%
  derive_vars_merged(
    dataset_add = ex_ext,
    filter_add = (EXDOSE > 0 |
                    (EXDOSE == 0 &
                       str_detect(EXTRT, "PLACEBO"))) & !is.na(EXENDTM),
    new_vars = exprs(TRTEDTM = EXENDTM, TRTETMF = EXENTMF),
    order = exprs(EXENDTM, EXSEQ),
    mode = "last",
    by_vars = exprs(STUDYID, USUBJID)
  )

#This call returns the original data frame with the column TRTSDTM, TRTSTMF, TRTEDTM, and TRTETMF added.
# adsl <- adsl %>%
#   derive_vars_dtm_to_dt(source_vars = exprs(TRTSDTM, TRTEDTM))


adsl<- adsl %>%  derive_vars_dtm_to_dt(
  source_vars = exprs(TRTSDTM, TRTEDTM)
)





#Treatment duration (TRTDURD).
adsl <- adsl %>%
  derive_var_trtdurd()

adsl <- adsl %>%
  mutate(
    TRTSDTM = if_else(
      TRTSTMF == "H",
      paste0(TRTSDTM, " 00:00:00"),
      NA_character_
    )
  )

#To add the End of Study date (EOSDT) 
# Convert character date to numeric date without imputation
ds_ext <- derive_vars_dt(
  ds,
  dtc = DSSTDTC,
  new_vars_prefix = "DSST"
)

adsl <- adsl %>%
  derive_vars_merged(
    dataset_add = ds_ext,
    by_vars = exprs(STUDYID, USUBJID),
    new_vars = exprs(EOSDT = DSSTDT),
    filter_add = DSCAT == "DISPOSITION EVENT" & DSDECOD != "SCREEN FAILURE"
  )

format_eosstt <- function(x) {
  case_when(
    x %in% c("COMPLETED") ~ "COMPLETED",
    x %in% c("SCREEN FAILURE") ~ NA_character_,
    TRUE ~ "DISCONTINUED"
  )
}

adsl <- adsl %>%
  derive_vars_merged(
    dataset_add = ds,
    by_vars = exprs(STUDYID, USUBJID),
    filter_add = DSCAT == "DISPOSITION EVENT",
    new_vars = exprs(EOSSTT = format_eosstt(DSDECOD)),
    missing_values = exprs(EOSSTT = "ONGOING")
  )

#To derive the End of Study reason(s)

adsl <- adsl %>%
  derive_vars_merged(
    dataset_add = ds,
    by_vars = exprs(USUBJID),
    new_vars = exprs(DCSREAS = DSDECOD, DCSREASP = DSTERM),
    filter_add = DSCAT == "DISPOSITION EVENT" &
      !(DSDECOD %in% c("SCREEN FAILURE", "COMPLETED", NA))
  )
# adsl <- adsl %>%
#   derive_vars_merged(
#     dataset_add = ds,
#     by_vars = exprs(USUBJID),
#     new_vars = exprs(DCSREAS = DSDECOD),
#     filter_add = DSCAT == "DISPOSITION EVENT" &
#       DSDECOD %notin% c("SCREEN FAILURE", "COMPLETED", NA)
#   ) %>%
#   derive_vars_merged(
#     dataset_add = ds,
#     by_vars = exprs(USUBJID),
#     new_vars = exprs(DCSREASP = DSTERM),
#     filter_add = DSCAT == "DISPOSITION EVENT" & DSDECOD %in% "OTHER"
#   )
# Randomization Date (RANDDT)

adsl <- adsl %>%
  derive_vars_merged(
    dataset_add = ds_ext,
    filter_add = DSDECOD == "RANDOMIZED",
    by_vars = exprs(STUDYID, USUBJID),
    new_vars = exprs(RANDDT = DSSTDT)
  )
# Derive Birth Date and Analysis Age (BRTHDT, AAGE, AAGEU)

# Derive birth date from BRTHDTC
adsl <- adsl %>%
  derive_vars_dt(
    new_vars_prefix = "BRTH",
    dtc = BRTHDTC
  )

adsl <- adsl %>%
  derive_vars_aage(
    start_date = BRTHDT,
    end_date = RANDDT
  )

# Derive Death Variables
adsl <- adsl %>%
  derive_vars_dt(
    new_vars_prefix = "DTH",
    dtc = DTHDTC
  )
# Death Reason
adsl <- adsl %>%
  derive_vars_extreme_event(
    by_vars = exprs(STUDYID, USUBJID),
    events = list(
      event(
        dataset_name = "ae",
        condition = AEOUT == "FATAL",
        set_values_to = exprs(DTHCAUS = AEDECOD),
      ),
      event(
        dataset_name = "ds",
        condition = DSDECOD == "DEATH" & grepl("DEATH DUE TO", DSTERM),
        set_values_to = exprs(DTHCAUS = DSTERM),
      )
    ),
    source_datasets = list(ae = ae, ds = ds),
    tmp_event_nr_var = event_nr,
    order = exprs(event_nr),
    mode = "first",
    new_vars = exprs(DTHCAUS)
  )
# some traceability variables (e.g. DTHDOM would store the domain where the date of death is collected, and DTHSEQ would store the xxSEQ value of that domain).
adsl <- adsl %>%
  select(-DTHCAUS) %>% # Remove it before deriving it again
  derive_vars_extreme_event(
    by_vars = exprs(STUDYID, USUBJID),
    events = list(
      event(
        dataset_name = "ae",
        condition = AEOUT == "FATAL",
        set_values_to = exprs(DTHCAUS = AEDECOD, DTHDOM = "AE", DTHSEQ = AESEQ),
      ),
      event(
        dataset_name = "ds",
        condition = DSDECOD == "DEATH" & grepl("DEATH DUE TO", DSTERM),
        set_values_to = exprs(DTHCAUS = DSTERM, DTHDOM = "DS", DTHSEQ = DSSEQ),
      )
    ),
    source_datasets = list(ae = ae, ds = ds),
    tmp_event_nr_var = event_nr,
    order = exprs(event_nr),
    mode = "first",
    new_vars = exprs(DTHCAUS, DTHDOM, DTHSEQ)
  )

# grouping of Death
adsl <- adsl %>%
  mutate(DTHCGR1 = case_when(
    is.na(DTHDOM) ~ NA_character_,
    DTHDOM == "AE" ~ "ADVERSE EVENT",
    str_detect(DTHCAUS, "(PROGRESSIVE DISEASE|DISEASE RELAPSE)") ~ "PROGRESSIVE DISEASE",
    TRUE ~ "OTHER"
  ))

# Duration Relative to Death
adsl <- adsl %>%
  derive_vars_duration(
    new_var = DTHADY,
    start_date = TRTSDT,
    end_date = DTHDT
  )

# Elapsed Days from Last Dose to Death

adsl <- adsl %>%
  derive_vars_duration(
    new_var = LDDTHELD,
    start_date = TRTEDT,
    end_date = DTHDT,
    add_one = FALSE
  )
#Derive the Last Date Known Alive (LSTALVDT)

adsl <- adsl %>%
  derive_vars_extreme_event(
    by_vars = exprs(STUDYID, USUBJID),
    events = list(
      event(
        dataset_name = "ae",
        order = exprs(AESTDTC, AESEQ),
        condition = !is.na(AESTDTC),
        set_values_to = exprs(
          LSTALVDT = convert_dtc_to_dt(AESTDTC, highest_imputation = "M"),
          seq = AESEQ
        ),
      ),
      event(
        dataset_name = "ae",
        order = exprs(AEENDTC, AESEQ),
        condition = !is.na(AEENDTC),
        set_values_to = exprs(
          LSTALVDT = convert_dtc_to_dt(AEENDTC, highest_imputation = "M"),
          seq = AESEQ
        ),
      ),
      event(
        dataset_name = "vs",
        order = exprs(VSDTC, VSSEQ),
        condition = !is.na(VSDTC),
        set_values_to = exprs(
          LSTALVDT = convert_dtc_to_dt(VSDTC, highest_imputation = "M"),
          seq = VSSEQ
        ),
      ),
      event(
        dataset_name = "lb",
        order = exprs(LBDTC, LBSEQ),
        condition = !is.na(LBDTC),
        set_values_to = exprs(
          LSTALVDT = convert_dtc_to_dt(LBDTC, highest_imputation = "M"),
          seq = LBSEQ
        ),
      ),
      event(
        dataset_name = "ds",
        order = exprs(DSSTDTC, DSSEQ),
        condition = !is.na(DSSTDTC),
        set_values_to = exprs(
          LSTALVDT = convert_dtc_to_dt(DSSTDTC, highest_imputation = "M"),
          seq = DSSEQ
        ),
      ),
      event(
        dataset_name = "adsl",
        condition = !is.na(TRTEDT),
        set_values_to = exprs(LSTALVDT = TRTEDT, seq = 0),
      )
    ),
    source_datasets = list(ae = ae,vs=vs,lb = lb,ds=ds,adsl = adsl),
    tmp_event_nr_var = event_nr,
    order = exprs(LSTALVDT, seq, event_nr),
    mode = "last",
    new_vars = exprs(LSTALVDT)
  )

#Tracebility FOR LSTALVDT

adsl <- adsl %>%
  select(-LSTALVDT) %>% # Created in the previous call
  derive_vars_extreme_event(
    by_vars = exprs(STUDYID, USUBJID),
    events = list(
      event(
        dataset_name = "ae",
        order = exprs(AESTDTC, AESEQ),
        condition = !is.na(AESTDTC),
        set_values_to = exprs(
          LSTALVDT = convert_dtc_to_dt(AESTDTC, highest_imputation = "M"),
          LALVSEQ = AESEQ,
          LALVDOM = "AE",
          LALVVAR = "AESTDTC"
        ),
      ),
      event(
        dataset_name = "ae",
        order = exprs(AEENDTC, AESEQ),
        condition = !is.na(AEENDTC),
        set_values_to = exprs(
          LSTALVDT = convert_dtc_to_dt(AEENDTC, highest_imputation = "M"),
          LALVSEQ = AESEQ,
          LALVDOM = "AE",
          LALVVAR = "AEENDTC"
        ),
      ),
      event(
        dataset_name = "vs",
        order = exprs(VSDTC, VSSEQ),
        condition = !is.na(VSDTC),
        set_values_to = exprs(
          LSTALVDT = convert_dtc_to_dt(VSDTC, highest_imputation = "M"),
          LALVSEQ = VSSEQ,
          LALVDOM = "VS",
          LALVVAR = "VSDTC"
        ),
      ),
      event(
        dataset_name = "lb",
        order = exprs(LBDTC, LBSEQ),
        condition = !is.na(LBDTC),
        set_values_to = exprs(
          LSTALVDT = convert_dtc_to_dt(LBDTC, highest_imputation = "M"),
          LALVSEQ = LBSEQ,
          LALVDOM = "LB",
          LALVVAR = "LBDTC"
        ),
      ),
      event(
        dataset_name = "ds",
        order = exprs(DSSTDTC, DSSEQ),
        condition = !is.na(DSSTDTC),
        set_values_to = exprs(
          LSTALVDT = convert_dtc_to_dt(DSSTDTC, highest_imputation = "M"),
          LALVSEQ = DSSEQ,
          LALVDOM = "DS",
          LALVVAR = "DSSTDTC"
        ),
      ),
      event(
        dataset_name = "adsl",
        condition = !is.na(TRTEDT),
        set_values_to = exprs(LSTALVDT = TRTEDT, LALVSEQ = NA_integer_, LALVDOM = "ADSL", LALVVAR = "TRTEDTM"),
      )
    ),
    source_datasets = list(ae = ae,vs=vs, lb = lb,ds=ds, adsl = adsl),
    tmp_event_nr_var = event_nr,
    order = exprs(LSTALVDT, LALVSEQ, event_nr),
    mode = "last",
    new_vars = exprs(LSTALVDT, LALVSEQ, LALVDOM, LALVVAR)
  )

# Grouping (e.g. AGEGR1 or REGION1)

# Create lookup tables
agegr1_lookup <- exprs(
  ~condition,           ~AGEGR1,
  AGE < 18,               "<18",
  between(AGE, 18, 64), "18-64",
  AGE > 64,               ">64",
  is.na(AGE),         "Missing"
)

region1_lookup <- exprs(
  ~condition,                          ~REGION1,
  COUNTRY %in% c("CAN", "USA"), "North America",
  !is.na(COUNTRY),          "Rest of the World",
  is.na(COUNTRY),                     "Missing"
)

format_agegr1 <- function(var_input) {
  case_when(
    var_input < 18 ~ "<18",
    between(var_input, 18, 64) ~ "18-64",
    var_input > 64 ~ ">64",
    TRUE ~ "Missing"
  )
}
format_region1 <- function(var_input) {
  case_when(
    var_input %in% c("CAN", "USA") ~ "North America",
    !is.na(var_input) ~ "Rest of the World",
    TRUE ~ "Missing"
  )
}

adsl <- adsl %>%
  mutate(
    AGEGR1 = format_agegr1(AAGE),
    REGION1 = format_region1(COUNTRY)
  )

format_agegr9 <- function(var_input) {
  case_when(
    var_input < 18 ~ "<18",
    between(var_input, 18, 50) ~ "18-50",
    var_input > 50 ~ ">50",
    TRUE ~ "Missing"
  )
}
format_agegr91 <- function(var_input) {
  case_when(
    var_input < 18 ~ 1,
    between(var_input, 18, 50) ~ 2,
    var_input > 50 ~ 3,
    TRUE ~ 99
  )
}

adsl <- adsl %>%
  mutate(
    AGEGR9 = format_agegr9(AAGE),
    AGEGR9N = format_agegr91(AAGE),
  )

#POPULATION FLAGS
#SAFFL
adsl <- adsl %>%
  derive_var_merged_exist_flag(
    dataset_add = ex,
    by_vars = exprs(STUDYID, USUBJID),
    new_var = SAFFL,
    false_value = "N",
    missing_value = "N",
    condition = (EXDOSE > 0 | (EXDOSE == 0 & str_detect(EXTRT, "PLACEBO")))
  )
#ITTFL
adsl <- adsl %>%
  derive_var_merged_exist_flag(
    dataset_add = dm,
    by_vars = exprs(STUDYID, USUBJID),
    new_var = ITTFL,
    false_value = "N",
    missing_value = "N",
    condition = (ARM != "")
  )

attr(adsl$STUDYID, "label") <- "Study Identifier"
attr(adsl$USUBJID, "label") <- "Unique Subject Identifier"
attr(adsl$AGE, "label") <- "Age"
attr(adsl$SEX, "label") <- "Sex"
attr(adsl$RACE, "label") <- "Race"

attr(adsl$TRT01P, "label") <- "Planned Treatment for Period 01"
attr(adsl$TRT01A, "label") <- "Actual Treatment for Period 01"

attr(adsl$TRTSDT, "label") <- "Date of First Exposure to Treatment"
attr(adsl$TRTEDT, "label") <- "Date of Last Exposure to Treatment"

attr(adsl$TRTSTMF, "label") <- "Date of First Exposure to Treatment flag"
attr(adsl$TRTETMF, "label") <- "Date of Last Exposure to Treatment flag"


attr(adsl$TRTSDTM, "label") <- "Date/Time of First Exposure to Treatment"
attr(adsl$TRTEDTM, "label") <- "Date/Time of Last Exposure to Treatment"

attr(adsl$RANDDT, "label") <- "Randomized Date"
attr(adsl$BRTHDT, "label") <- "Birth Date"
attr(adsl$TRTDURD, "label") <- "Treatment Duration in Days"

attr(adsl$EOSDT, "label") <- "End of Study Date"
attr(adsl$EOSSTT, "label") <- "End of Study Status"

attr(adsl$SAFFL, "label") <- "Safety Population Flag"
attr(adsl$ITTFL, "label") <- "Intent-to-Treat Population Flag"

attr(adsl$AAGE,     "label") <- "Analysis Age"
attr(adsl$AAGEU,    "label") <- "Analysis Age Units"
attr(adsl$DTHDT,    "label") <- "Date of Death"
attr(adsl$DTHCAUS,  "label") <- "Cause of Death"
attr(adsl$DTHDOM,   "label") <- "Death Domain"
attr(adsl$DTHSEQ,   "label") <- "Death Sequence Number"
attr(adsl$DTHCGR1,  "label") <- "Death Cause Category 1"
attr(adsl$DTHADY,   "label") <- "Analysis Day of Death"
attr(adsl$LDDTHELD, "label") <- "Last Death Data Held"
attr(adsl$LSTALVDT, "label") <- "Last Known Alive Date"
attr(adsl$LALVSEQ,  "label") <- "Last Known Alive Sequence Number"
attr(adsl$LALVDOM,  "label") <- "Last Known Alive Domain"
attr(adsl$LALVVAR,  "label") <- "Last Known Alive Variable"
attr(adsl$AGEGR1,   "label") <- "Age Group 1"
attr(adsl$REGION1,  "label") <- "Region 1"
attr(adsl$AGEGR9,   "label") <- "Age Group 9"
attr(adsl$AGEGR9N,  "label") <- "Age Group 9 Numeric"

write_parquet(adsl, "adsl.parquet")
cat("ADSL program completed\n")

sink()

