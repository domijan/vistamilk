library(haven)
library(tidyverse)
mixed.data <- read_sas("dataset_for_mixed_model.sas7bdat")
mixed.data%>% glimpse()
y <- mixed.data$delta_bcs
x <- mixed.data %>% select(PM5:PM752)
glimpse(x)
matplot(t(x), type = "l")
