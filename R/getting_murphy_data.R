library("tabulapdf")
# set Java memory limit to 600 MB (optional)
options(java.parameters = "-Xmx600m")
f <- system.file("data", "sms14075-sup-0003-tables2.pdf")
# extract table from first page of example PDF
tab1 <- extract_tables("data/sms14075-sup-0003-tables2.pdf", pages = 1)[[1]] |>
  slice_tail(n = 31) |>
  select(c(1,13)) |>
  rename(
    study = "...1",
    effect = "...13"
  ) |>
  separate(effect, into = c("outcome", "effect"), sep = ":") |>
  separate(effect, into = c("yi", "ci"), sep = "\\(") |>
  separate(ci, into = c("ci.lb", "ci.ub"), sep = ", ", extra = "drop") %>%
  mutate(ci.lb = str_remove_all(ci.lb, "\\("),
         ci.ub = str_remove_all(ci.ub, "\\)"))

write.csv(tab1, "data/tab1.csv")

tab2 <- extract_tables("data/sms14075-sup-0003-tables2.pdf", pages = 2)[[1]]  |>
  slice_tail(n = 36) |>
  select(c(1,13)) |>
  rename(
    study = "...1",
    effect = "...13"
  ) |>
  separate(effect, into = c("outcome", "effect"), sep = ":") |>
  separate(effect, into = c("yi", "ci"), sep = "\\(") |>
  separate(ci, into = c("ci.lb", "ci.ub"), sep = ", ", extra = "drop") %>%
  mutate(ci.lb = str_remove_all(ci.lb, "\\("),
         ci.ub = str_remove_all(ci.ub, "\\)"))

write.csv(tab2, "data/tab2.csv")

tab3 <- extract_tables("data/sms14075-sup-0003-tables2.pdf", pages = 3)[[1]]  |>
  slice_tail(n = 38) |>
  select(c(1,13)) |>
  rename(
    study = "...1",
    effect = "...13"
  ) |>
  separate(effect, into = c("outcome", "effect"), sep = ":") |>
  separate(effect, into = c("yi", "ci"), sep = "\\(") |>
  separate(ci, into = c("ci.lb", "ci.ub"), sep = ", ", extra = "drop") %>%
  mutate(ci.lb = str_remove_all(ci.lb, "\\("),
         ci.ub = str_remove_all(ci.ub, "\\)"))

write.csv(tab3, "data/tab3.csv")

tab4 <- extract_tables("data/sms14075-sup-0003-tables2.pdf", pages = 4)[[1]]  |>
  slice_tail(n = 12) |>
  select(c(1,13)) |>
  rename(
    study = "...1",
    effect = "...13"
  ) |>
  separate(effect, into = c("outcome", "effect"), sep = ":") |>
  separate(effect, into = c("yi", "ci"), sep = "\\(") |>
  separate(ci, into = c("ci.lb", "ci.ub"), sep = ", ", extra = "drop") %>%
  mutate(ci.lb = str_remove_all(ci.lb, "\\("),
         ci.ub = str_remove_all(ci.ub, "\\)"))

write.csv(tab4, "data/tab4.csv")

# Some manual editing for the study names

data_murphy_koehler <- bind_rows(
  read_csv("data/tab1.csv"),
  read_csv("data/tab2.csv"),
  read_csv("data/tab3.csv"),
  read_csv("data/tab4.csv")
)

write_csv(data_murphy_koehler, "data/murphy_koehler_study_data.csv")
  
