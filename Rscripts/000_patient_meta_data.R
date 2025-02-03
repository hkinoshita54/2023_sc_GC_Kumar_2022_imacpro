library(tidyverse)
library(readxl)

pts_sample <- read_excel("meta_data/pts_sample.xlsx") %>% as.data.frame
names(pts_sample) <- c("patient_ID", "pri_T", "pri_NT", "peri_T", "peri_NT")
pts_sample <- pivot_longer(pts_sample, cols = c("pri_T", "pri_NT", "peri_T", "peri_NT"), names_to = "tissue_type", values_to = "sample_ID")
pts_sample <- filter(pts_sample, sample_ID != "-")
pts_sample$sample_ID <- sub(pattern = ".csv", replacement = "", pts_sample$sample_ID)

GSM_sample <- read_excel("meta_data/GSM_sample.xlsx", col_names = c("GSM", "sample_ID", "tissue")) %>% as.data.frame

GSM_sample_pts <- inner_join(GSM_sample, pts_sample, by = "sample_ID")

supp_table1 <- read_excel("meta_data/supp_table1.xlsx") %>% as.data.frame
supp_table1 <- supp_table1[,2:11]
names(supp_table1)[1] <- "patient_ID"

meta_data <- inner_join(GSM_sample_pts, supp_table1, by = "patient_ID")
meta_data$sample_ID <- sub(pattern = "sample", replacement = "", meta_data$sample_ID) %>% as.integer()

openxlsx2::write_xlsx(meta_data, "meta_data/patients_data.xlsx")
