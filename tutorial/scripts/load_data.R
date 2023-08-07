dt_raw <- readRDS(file="data/osf_mini.Rda")
var_meta <- readxl::read_excel("data/osf_var_meta.xlsx")

qs_in <- qs::qread("data/osf_feature.qs")
list2env(qs_in, envir = .GlobalEnv)
rm(qs_in)


