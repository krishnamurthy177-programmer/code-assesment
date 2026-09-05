
#krishna Murthy
#conversion raw dataset to SDTM.DS
#program name: 01_create_ds_domain.R



dir.create("logs", showWarnings = FALSE)

log_file <- "logs/DS_program.log"

sink(log_file)

cat("DS program started\n")

#install packages
install.packages(c("admiral", "sdtm.oak", "gt", "ggplot2"))

#input dataset
install.packages("pharmaverseraw")

pharmaverseraw::ds_raw

ds_raw <- pharmaverseraw::ds_raw

install.packages("dplyr")
install.packages("tidyverse")
library(dplyr)
library(tidyverse)

study_ct <-
  data.frame(
    stringsAsFactors = FALSE,
    codelist_code = c("C66727","C66727",
                      "C66727","C66727","C66727","C66727","C66727","C66727",
                      "C66727","C66727"),
    term_code = c("C41331","C25250",
                  "C28554","C48226","C48227","C48250","C142185","C49628",
                  "C49632","C49634"),
    term_value = c("ADVERSE EVENT",
                   "COMPLETED","DEATH","LACK OF EFFICACY","LOST TO FOLLOW-UP",
                   "PHYSICIAN DECISION","PROTOCOL VIOLATION",
                   "SCREEN FAILURE","STUDY TERMINATED BY SPONSOR",
                   "WITHDRAWAL BY SUBJECT"),
    collected_value = c("Adverse Event",
                        "Complete","Dead","Lack of Efficacy","Lost To Follow-Up",
                        "Physician Decision","Protocol Violation",
                        "Trial Screen Failure","Study Terminated By Sponsor",
                        "Withdrawal by Subject"),
    term_preferred_term = c("AE","Completed","Died",
                            NA,NA,NA,"Violation",
                            "Failure to Meet Inclusion/Exclusion Criteria",NA,"Dropout"),
    term_synonyms = c("ADVERSE EVENT",
                      "COMPLETE","Death",NA,NA,NA,NA,NA,NA,
                      "Discontinued Participation")
  ) 



install.packages("sdtm.oak")
library(sdtm.oak)

ds_raw <- ds_raw %>%
  mutate(
    IT.DSTERM=if_else(is.na(OTHERSP),IT.DSTERM,OTHERSP)
  )

ds_raw <- ds_raw %>%
  generate_oak_id_vars(
    pat_var = "PATNUM",
    raw_src = "ds_raw"
  )

ds <-
  # Map topic variable
  assign_no_ct(
    raw_dat = ds_raw,
    raw_var = "IT.DSTERM",
    tgt_var = "DSTERM"
  )

ds <- ds %>%
  # Map qualifier DSDECOD
  assign_ct(
    raw_dat = ds_raw,
    raw_var = "IT.DSDECOD",
    tgt_var = "DSDECOD",
    ct_spec = study_ct,
    ct_clst = "C66727",
    id_vars = oak_id_vars()
  )
ds <- ds %>%
  # Map qualifier DSDECOD
  assign_ct(
    raw_dat = ds_raw,
    raw_var = "IT.DSDECOD",
    tgt_var = "DSDECOD",
    ct_spec = study_ct,
    ct_clst = "C66727",
    id_vars = oak_id_vars()
  )
# Map qualifier DSCAT
ds <- ds %>% 
  mutate(
    DSCAT=case_when(is.na(ds_raw$IT.DSDECOD) ~ "OTHER EVENT",
                    ds_raw$IT.DSDECOD=="Randomized" ~ "PROTOCOL MILESTONE",
                    TRUE ~"DISPOSITION EVENT" ),)

ds <- ds %>%
  # Map DSSTDTC. This function calls create_iso8601
  assign_datetime(
    raw_dat = ds_raw,
    raw_var = c("IT.DSSTDAT"),
    tgt_var = "DSSTDTC",
    raw_fmt = c(list(c("d-m-y", "dd mmm yyyy"))),
    raw_unk = c("UN", "UNK"),
    id_vars = oak_id_vars()
  )

ds <- ds %>%
  # Map DSDTC. This function calls create_iso8601
  assign_datetime(
    raw_dat = ds_raw,
    raw_var = c("DSDTCOL", "DSTMCOL"),
    tgt_var = "DSDTC",
    raw_fmt = c(list(c("d-m-y", "dd mmm yyyy")), "H:M"),
    raw_unk = c("UN", "UNK"),
    id_vars = oak_id_vars()
  )

# Map qualifier DSSTDY
install.packages("pharmaversesdtm")
pharmaversesdtm::dm

dm<-pharmaversesdtm::dm

ds <- ds %>% mutate(id=str_sub(ds_raw$STUDY, -2, -1))%>%
  mutate(
    STUDYID = ds_raw$STUDY,
    DOMAIN = "DS",
    DSCAT = DSCAT,
    USUBJID = paste0(id, "-",ds_raw$PATNUM)
  ) %>%
  derive_study_day(
    sdtm_in = .,
    dm_domain = dm,
    tgdt = "DSSTDTC",
    refdt = "RFXSTDTC",
    study_day_var = "DSSTDY"
  )%>%
  # derive sequence number
  derive_seq(tgt_var = "DSSEQ",
             rec_vars= c("USUBJID", "DSTERM"))


#ATTRIBUTES
install.packages("haven")
library(haven)

attr(ds$STUDYID, "label") <- "Study Identifier"
attr(ds$DOMAIN,  "label") <- "Domain Abbreviation"
attr(ds$USUBJID, "label") <- "Unique Subject Identifier"
attr(ds$DSSEQ,   "label") <- "Sequence Number"
attr(ds$DSTERM,  "label") <- "Reported Term for Disposition Event"
attr(ds$DSDECOD, "label") <- "Standardized Disposition Term"
attr(ds$DSCAT,   "label") <- "Category"
attr(ds$DSDTC,   "label") <- "Date/Time of Collection"
attr(ds$DSSTDTC, "label") <- "Start Date/Time of Disposition Event"
attr(ds$DSSTDY,  "label") <- "Study Day of Start of Disposition Event"


DS <- ds %>%
  select(
    STUDYID, DOMAIN, USUBJID, DSSEQ, DSTERM, DSDECOD, DSCAT, DSDTC,
    DSSTDTC, DSSTDY
  )

# as per SDTMIG there will be no VISIT,VISTNUM 

cat("DS program completed\n")

sink()








