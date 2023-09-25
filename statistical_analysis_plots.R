## Header ---------------------------
#
# Script name: sano_stats.R
#
# Purpose of script: Compute statistics for Dr. Sano's tumor scoring project
#
# Author: Ethan N. Okoshi
#
# Date Created: 2023-06-22
#
# Copyright (c) Ethan Okoshi, 2023
# Email: ethanokoshi@gmail.com
#
# Notes:
#
# set working directory
#
# setwd()
#
# set options
#
options(scipen = 6, digits = 4)
Sys.setenv(LANG = "EN")
extrafont::loadfonts(device = "win")
#
# load functions from another file
#
# source()
#
## load packages -------------
#
require(tidyverse)
library(magrittr)
library(pROC)
library(gtsummary)
library(gt)
library(coin)
library(irr)
library(tidymodels)
#
## load data -------------
#
input_data_path <- "~/ethancoding/Sano/sanoproject_ENO.xlsx"
library(readxl)
sano <- read_excel(input_data_path, sheet = "sano")
tachibana <- read_excel(input_data_path, sheet = "tachibana")
brcic <- read_excel(input_data_path, sheet = "brcic")
fukuoka <- read_excel(input_data_path, sheet = "fukuoka")
consensus3 <- read_excel(input_data_path, sheet = "consensus3")
consensus4 <- read_excel(input_data_path, sheet = "consensus4")


no3sano <-
  sano[sano$target != 3, ] %>%
  mutate(across(intbronch:lymphoidagg, na_if, "NA"),
    across(intbronch:op, as.ordered),
    across(smokehx, na_if, "NA"),
    maligvsbenign = if_else(target == 1 | target == 2, 1, 0),
    .before = Location,
    maligvsbenign = as.factor(maligvsbenign),
    bronch_epithelialcells = as.factor(bronch_epithelialcells),
    across(target, factor,
      levels = c(5, 4, 2, 1),
      labels = c(
        "Benign",
        "Probably Benign",
        "Probably Malignant",
        "Malignant"
      ),
      ordered = TRUE
    ),
    Sex = as.factor(Sex)
  )

no3tachibana <-
  tachibana[tachibana$target != 3, ] %>%
  mutate(across(intbronch:lymphoidagg, na_if, "NA"),
    across(intbronch:op, factor,
      levels = c(0, 1, 2, 3),
      ordered = TRUE
    ),
    maligvsbenign = if_else(target == 1 | target == 2, 1, 0),
    maligvsbenign = as.factor(maligvsbenign)
  ) %>%
  select(!target)

no3brcic <- 
  brcic[brcic$target != 3, ] %>% 
  mutate(across(intbronch:lymphoidagg, na_if, "NA"),
         across(intbronch:op, factor,
                levels = c(0, 1, 2, 3),
                ordered = TRUE
         ),
         maligvsbenign = if_else(target == 1 | target == 2, 1, 0),
         maligvsbenign = as.factor(maligvsbenign)
  ) %>%
  select(!target)

no3fukuoka <- 
  fukuoka[fukuoka$target != 3, ] %>% 
  mutate(across(intbronch:lymphoidagg, na_if, "NA"),
         across(intbronch:op, factor,
                ordered = TRUE
         ),
         maligvsbenign = if_else(target == 1 | target == 2, 1, 0),
         maligvsbenign = as.factor(maligvsbenign)
  ) %>%
  select(!target)

no3consensus3 <- 
  consensus3[consensus3$target != 3, ] %>% 
  mutate(across(intbronch:lymphoidagg, na_if, "NA"),
         across(intbronch:op, factor,
                levels = c(0, 1, 2, 3),
                ordered = TRUE
         ),
         maligvsbenign = if_else(target == 1 | target == 2, 1, 0),
         maligvsbenign = as.factor(maligvsbenign)
  ) %>%
  select(!target)

no3consensus4 <- 
  consensus4[consensus4$target != 3, ] %>% 
  mutate(across(intbronch:lymphoidagg, na_if, "NA"),
         across(intbronch:op, factor,
                levels = c(0, 1, 2, 3),
                ordered = TRUE
         ),
         maligvsbenign = if_else(target == 1 | target == 2, 1, 0),
         maligvsbenign = as.factor(maligvsbenign)
  ) %>%
  select(!target)

get_expected_values <- function(data, variable, by) {
  chisq.test(data[[variable]], data[[by]])$expected
}

## End Header --------------


## Bronchial Epithelial Cells USELESS ----------------------------------------------

CPCOLS <- c("#FF5784", "#05EBC8")

ggplot(no3sano, aes(bronch_epithelialcells, fill = bronch_epithelialcells)) +
  geom_bar(stat = "count", show.legend = F) +
  scale_x_discrete(labels = c("Absent", "Present")) +
  scale_fill_manual(values = CPCOLS, ) +
  labs(x = "Bronchial Epithelium", y = "Count") +
  geom_text(
    aes(label = after_stat(count)),
    stat = "count",
    position = "identity",
    vjust = 1.3,
    size = 5,
    color = "white",
    fontface = "bold"
  ) +
  ggpubr::theme_classic2()
ggsave("sano_bronchEpcount.png", dpi = 320, width = 12, height = 12, units = "cm")

## Age -------------------------------------------
ggplot(data = no3sano) +
  geom_dotplot(
    binaxis = "y",
    stackdir = "center",
    binwidth = 1,
    aes(
      x = target,
      y = Age,
      group = target,
      fill = target
    )
  ) +
  scale_fill_brewer(
    "Diagnosis",
    palette = "Set1"
  ) +
  labs(x = "Diagnosis") +
  theme(axis.text.x = element_text(
    angle = 45,
    vjust = 1,
    hjust = 1
  ))
ggsave("sano_ageVSdiag.png", dpi = 320, width = 15, height = 15, units = "cm")

# roc curve
png(file = "ageauc.png")
plot.roc(
  maligvsbenign ~ Age,
  data = no3sano,
  legacy.axes = FALSE,
  print.auc = TRUE,
  print.auc.x = 0.3,
  print.auc.y = 0.1,
  direction = "<"
)
dev.off()

## Sex ----------------------------------
ggplot(data = no3sano, aes(target, fill = Sex)) +
  geom_bar(position = "dodge") +
  scale_fill_brewer("Sex",
    palette = "Accent",
    labels = c("Female", "Male")
  ) +
  labs(x = "Diagnosis", y = "Count") +
  theme(
    axis.text.x = element_text(
      angle = 45,
      vjust = 1,
      hjust = 1
    ),
    legend.background = element_rect(colour = "black", linewidth = 0.3)
  )
ggsave("sano_sexVSdiag.png", dpi = 320, width = 15, height = 15, units = "cm")

png(file = "sexauc.png")
plot.roc(
  maligvsbenign ~ as.ordered(Sex),
  data = no3sano,
  legacy.axes = FALSE,
  print.auc = TRUE,
  print.auc.x = 0.3,
  print.auc.y = 0.1
)
dev.off()

## intbronch ------------------

# linear-by-linear test
no3sano %>% filter(!is.na(intbronch)) %$%
  lbl_test(target ~ intbronch,
    distribution = approximate(nresample = 9999L)
  )

# fisher 2x2 test
no3sano %>%
  select(maligvsbenign, intbronch) %>%
  mutate(intbronch = if_else(intbronch == 0, "Negative", "Positive")) %>%
  table() %>%
  fisher.test()

# plot relative frequencies
ggplot(data = filter(no3sano, !is.na(intbronch)), aes(x = target, fill = intbronch)) +
  geom_bar(position = "fill") +
  geom_text(
    aes(label = after_stat(count)),
    stat = "count",
    position = "fill",
    vjust = 1.1,
    size = 3,
    color = "black",
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = (wesanderson::wes_palette("Darjeeling1")),
    name = "Int. Bronch Grade"
  ) +
  labs(title = "Rate of Interface Bronchitis by Diagnosis Group", 
       x = "Diagnosis", 
       y = "Int. Bronch Grade, % of Group") +
  theme_minimal() +
  theme(
    text = element_text(family = "Consolas"),
    plot.margin = margin(1, 0, 1, 0.2, "cm"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.title.x.bottom = element_text(margin = margin(10, 0, 0, 0, "pt")),
    axis.title.y.left = element_text(margin = margin(0, 10, 0, 0, "pt")),
    axis.text = element_text(size = 8),
    axis.text.x.bottom = element_text(margin = margin(-5, 0, 0, 0)),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.y.left = element_text(vjust = 0.5),
    legend.title = element_text(angle = 90, vjust = 1)
  )
ggsave("sano_daigVSintbronchFRACTION.png", dpi = 320, height = 15, width = 15, units = "cm")

# plot roc curve
png(file = "intbronchauc.png")
plot.roc(
  maligvsbenign ~ intbronch,
  data = no3sano,
  legacy.axes = TRUE,
  print.auc = TRUE,
  print.auc.x = 0.3,
  print.auc.y = 0.1,
  direction = ">"
)
dev.off()

## plasmacellinfil ------------------

# linear-by-linear test
no3sano %>% filter(!is.na(plasmacellinfil)) %$%
  lbl_test(target ~ plasmacellinfil,
    distribution = approximate(nresample = 9999L)
  )

# fisher 2x2 test
no3sano %>%
  select(maligvsbenign, plasmacellinfil) %>%
  mutate(plasmacellinfil = if_else(plasmacellinfil == 0, "Negative", "Positive")) %>%
  table() %>%
  fisher.test()

# plot relative frequencies
ggplot(
  data = filter(no3sano, !is.na(plasmacellinfil)),
  aes(x = target, fill = plasmacellinfil)
) +
  geom_bar(position = "fill") +
  geom_text(
    aes(label = after_stat(count)),
    stat = "count",
    position = "fill",
    vjust = 1.1,
    size = 3,
    color = "white",
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = (wesanderson::wes_palette("Darjeeling2")),
    name = "Plasma Cell Infil. Grade"
  ) +
  labs(
    title = "Rate of Plasma Cell Infiltration by Diagnosis Group",
    x = "Diagnosis",
    y = "Plasma Cell Infil. Grade, % of Group"
  ) +
  theme_minimal() +
  theme(
    text = element_text(family = "Consolas"),
    plot.margin = margin(1, 0, 1, 0.2, "cm"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.title.x.bottom = element_text(margin = margin(10, 0, 0, 0, "pt")),
    axis.title.y.left = element_text(margin = margin(0, 10, 0, 0, "pt")),
    axis.text = element_text(size = 8),
    axis.text.x.bottom = element_text(margin = margin(-5, 0, 0, 0)),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.y.left = element_text(vjust = 0.5),
    legend.title = element_text(angle = 90, vjust = 1)
  )
ggsave("sano_daigVSplasmacellinfilFRACTION.png", 
       dpi = 320, 
       height = 15, 
       width = 15, 
       units = "cm")

# plot roc curve
png(file = "plasmacellinfilauc.png")
plot.roc(
  maligvsbenign ~ plasmacellinfil,
  data = no3sano,
  legacy.axes = TRUE,
  print.auc = TRUE,
  print.auc.x = 0.3,
  print.auc.y = 0.1,
  direction = ">"
)
dev.off()

## eosinoinfil ------------------

# linear-by-linear test
no3sano %>% filter(!is.na(eosinoinfil)) %$%
  lbl_test(target ~ eosinoinfil,
    distribution = approximate(nresample = 9999L)
  )

# fisher 2x2 test
no3sano %>%
  select(maligvsbenign, eosinoinfil) %>%
  mutate(eosinoinfil = if_else(eosinoinfil == 0, "Negative", "Positive")) %>%
  table() %>%
  fisher.test()

# plot relative frequencies
ggplot(
  data = filter(no3sano, !is.na(eosinoinfil)),
  aes(x = target, fill = eosinoinfil)
) +
  geom_bar(position = "fill") +
  geom_text(
    aes(label = after_stat(count)),
    stat = "count",
    position = "fill",
    vjust = 1.1,
    size = 3,
    color = "white",
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = (wesanderson::wes_palette("BottleRocket1")[2:5]),
    name = "Eosinophil Infil. Grade"
  ) +
  labs(
    title = "Rate of Eosinophil Infiltration by Diagnosis Group",
    x = "Diagnosis",
    y = "Eosinophil Infil. Grade, % of Group"
  ) +
  theme_minimal() +
  theme(
    text = element_text(family = "Consolas"),
    plot.margin = margin(1, 0, 1, 0.2, "cm"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.title.x.bottom = element_text(margin = margin(10, 0, 0, 0, "pt")),
    axis.title.y.left = element_text(margin = margin(0, 10, 0, 0, "pt")),
    axis.text = element_text(size = 8),
    axis.text.x.bottom = element_text(margin = margin(-5, 0, 0, 0)),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.y.left = element_text(vjust = 0.5),
    legend.title = element_text(angle = 90, vjust = 1)
  )
ggsave("sano_daigVSeosinoinfilFRACTION.png", 
       dpi = 320, 
       height = 15, 
       width = 15, 
       units = "cm")

# plot roc curve
png(file = "eosinoinfilauc.png")
plot.roc(
  maligvsbenign ~ eosinoinfil,
  data = no3sano,
  legacy.axes = TRUE,
  print.auc = TRUE,
  print.auc.x = 0.3,
  print.auc.y = 0.1,
  direction = ">"
)
dev.off()

## lymphoidagg ------------------

# linear-by-linear test
no3sano %>% filter(!is.na(lymphoidagg)) %$%
  lbl_test(target ~ lymphoidagg,
    distribution = approximate(nresample = 9999L)
  )

# fisher 2x2 test
no3sano %>%
  select(maligvsbenign, lymphoidagg) %>%
  mutate(lymphoidagg = if_else(lymphoidagg == 0, "Negative", "Positive")) %>%
  table() %>%
  fisher.test()

# plot relative frequencies
ggplot(
  data = filter(no3sano, !is.na(lymphoidagg)),
  aes(x = target, fill = lymphoidagg)
) +
  geom_bar(position = "fill") +
  geom_text(
    aes(label = after_stat(count)),
    stat = "count",
    position = "fill",
    vjust = 1.1,
    size = 3,
    color = "white",
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = (wesanderson::wes_palette("GrandBudapest2")[2:5]),
    name = "Lymphoid Aggregation Grade"
  ) +
  labs(title = "Lymphoid Aggregation by Diagnosis Group", 
       x = "Diagnosis", 
       y = "Lymphoid Aggregation Grade, % of Group") +
  theme_minimal() +
  theme(
    text = element_text(family = "Consolas"),
    plot.margin = margin(1, 0, 1, 0.2, "cm"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.title.x.bottom = element_text(margin = margin(10, 0, 0, 0, "pt")),
    axis.title.y.left = element_text(margin = margin(0, 10, 0, 0, "pt")),
    axis.text = element_text(size = 8),
    axis.text.x.bottom = element_text(margin = margin(-5, 0, 0, 0)),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.y.left = element_text(vjust = 0.5),
    legend.title = element_text(angle = 90, vjust = 1)
  )
ggsave("sano_daigVSlymphoidaggFRACTION.png", dpi = 320, height = 15, width = 15, units = "cm")

png(file = "lymphoidaggauc.png")
plot.roc(
  maligvsbenign ~ lymphoidagg,
  data = no3sano,
  legacy.axes = TRUE,
  print.auc = TRUE,
  print.auc.x = 0.3,
  print.auc.y = 0.1,
  direction = ">"
)
dev.off()

## fibroelastosis ----------------------------------------------------------

# linear-by-linear association test (chi sq for ordered categorical variables)
lbl_test(
  maligvsbenign ~ fibroelastosis,
  data = no3sano,
  distribution = approximate(nresample = 9999L)
)

# fisher test for 2x2 contingency table comparing malig vs benign and fibro pos vs neg
no3sano %>%
  select(maligvsbenign, fibroelastosis) %>%
  mutate(fibroelastosis = if_else(fibroelastosis == 0, "Negative", "Positive")) %>%
  table() %>%
  fisher.test()

# plot relative frequencies
ggplot(no3sano, aes(target, fill = fibroelastosis)) +
  geom_bar(position = "fill") +
  scale_fill_brewer("Fibroelastosis", palette = "Accent") +
  labs(x = "Diagnosis", 
       y = "Fraction of Diagnosis Group", 
       title = "Fibroelastosis by Diagnosis Group") +
  theme(
    text = element_text(family = "Arial"),
    axis.text.x = element_text(
      angle = 45,
      vjust = 1,
      hjust = 1
    ),
    legend.background = element_rect(colour = "black", linewidth = 0.3)
  )
ggsave("sano_diagVSfibroFRACTION.png", 
       dpi = 320, 
       height = 12, 
       width = 12, 
       units = "cm")

# plot roc curve
png(file = "fibroauc.png")
plot.roc(
  maligvsbenign ~ fibroelastosis,
  data = no3sano,
  legacy.axes = TRUE,
  print.auc = TRUE,
  print.auc.x = 0.3,
  print.auc.y = 0.1,
  direction = "<"
)
dev.off()

## OP ------------------

# linear-by-linear association test (chi sq for ordered categorical variables)
lbl_test(
  maligvsbenign ~ op,
  data = no3sano,
  distribution = approximate(nresample = 9999L)
)

# fisher test for 2x2 contingency table comparing malig vs benign and fibro pos vs neg
no3sano %>%
  select(maligvsbenign, op) %>%
  mutate(op = if_else(op == 0, "Negative", "Positive")) %>%
  table() %>%
  fisher.test()

# plot relative frequencies
ggplot(no3sano, aes(target, fill = op)) +
  geom_bar(position = "fill") +
  geom_text(
    aes(label = after_stat(count)),
    stat = "count",
    position = "fill",
    vjust = 1.2,
    size = 3,
    color = "white",
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = (wesanderson::wes_palette("Zissou1"))[2:5],
    name = "OP Grade"
  ) +
  labs(title = "Rate of OP by Diagnosis Group", 
       x = "Diagnosis", 
       y = "OP Grade % of Group") +
  theme_minimal() +
  theme(
    text = element_text(family = "Consolas"),
    plot.margin = margin(1, 0, 1, 0.2, "cm"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.title.x.bottom = element_text(margin = margin(10, 0, 0, 0, "pt")),
    axis.title.y.left = element_text(margin = margin(0, 10, 0, 0, "pt")),
    axis.text = element_text(size = 8),
    axis.text.x.bottom = element_text(margin = margin(-5, 0, 0, 0)),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.y.left = element_text(vjust = 0.5)
  )
ggsave("sano_daigVSopFRACTION.png", 
       dpi = 320, 
       height = 15, 
       width = 15, 
       units = "cm")

# plot roc curve
png(file = "opauc.png")
plot.roc(
  maligvsbenign ~ op,
  data = no3sano,
  legacy.axes = TRUE,
  print.auc = TRUE,
  print.auc.x = 0.3,
  print.auc.y = 0.1,
  direction = ">"
)
dev.off()

## gtsummary clinical --------------

no3sano %>%
  select(Age, Sex, smokehx, maligvsbenign) %>%
  mutate(maligvsbenign = if_else(maligvsbenign == 1, "Malignant", "Benign")) %>%
  tbl_summary(
    by = maligvsbenign,
    missing = "ifany",
    missing_text = "Missing",
    label = list(smokehx ~ "Smoking History"),
    digits = list(Age ~ 2),
    statistic = list(Age ~ "{mean} ({sd})"),
  ) %>%
  add_p(
    test = all_continuous() ~ "kruskal.test",
    test.args = all_tests("fisher.test") ~ list(simulate.p.value = TRUE, B = 10000)
  ) %>%
  bold_p() %>%
  add_overall() %>%
  separate_p_footnotes() %>%
  modify_header(label = "**Factor**") %>%
  as_gt() %>%
  gt::gtsave("clinicalstats_table.png", zoom = 10)

## gtsummary sano -------------


no3sano %>%
  select(intbronch:op, Sex) %>%
  colnames() %>%
  map(get_expected_values, data = no3sano, by = "maligvsbenign")


no3sano %>%
  select(intbronch:op, maligvsbenign) %>%
  mutate(maligvsbenign = if_else(maligvsbenign == 1, "Malignant", "Benign")) %>%
  tbl_summary(
    by = maligvsbenign,
    missing = "ifany",
    missing_text = "Missing",
    label = list(
      intbronch ~ "Interface Bronchitis",
      plasmacellinfil ~ "Plasma Cell Infiltration",
      eosinoinfil ~ "Eosinophil Infiltration",
      lymphoidagg ~ "Lymphoid Aggregation",
      fibroelastosis ~ "Fibroelastosis",
      op ~ "Organized Pneumonia"
    ),
  ) %>%
  add_p(test.args = all_tests("fisher.test") ~ list(simulate.p.value = TRUE, B = 10000)) %>%
  bold_p() %>%
  add_overall() %>%
  separate_p_footnotes() %>%
  modify_header(label = "**Factor**") %>%
  as_gt() %>%
  gt::gtsave("sanostats_table.png", zoom = 10)



## gtsummary tachibana ---------------

no3tachibana %>%
  select(intbronch:op) %>%
  colnames() %>%
  map(get_expected_values, data = no3tachibana, by = "maligvsbenign")

no3tachibana %>%
  select(intbronch:op, maligvsbenign) %>%
  mutate(maligvsbenign = if_else(maligvsbenign == 1, "Malignant", "Benign")) %>%
  tbl_summary(
    by = maligvsbenign,
    missing = "ifany",
    missing_text = "Missing",
    label = list(
      intbronch ~ "Interface Bronchitis",
      plasmacellinfil ~ "Plasma Cell Infiltration",
      eosinoinfil ~ "Eosinophil Infiltration",
      lymphoidagg ~ "Lymphoid Aggregation",
      fibroelastosis ~ "Fibroelastosis",
      op ~ "Organized Pneumonia"
    ),
  ) %>%
  add_p(test.args = all_tests("fisher.test") ~ list(simulate.p.value = TRUE, B = 10000)) %>%
  bold_p() %>%
  add_overall() %>%
  separate_p_footnotes() %>%
  modify_header(label = "**Factor**") %>%
  as_gt() %>%
  gt::gtsave("tachibanastats_table.png", zoom = 10)


## gtsummary brcic --------------------

no3brcic %>%
  select(intbronch:op) %>%
  colnames() %>%
  map(get_expected_values, data = no3brcic, by = "maligvsbenign")

no3brcic %>%
  select(intbronch:op, maligvsbenign) %>%
  mutate(maligvsbenign = if_else(maligvsbenign == 1, "Malignant", "Benign")) %>%
  tbl_summary(
    by = maligvsbenign,
    missing = "ifany",
    missing_text = "Missing",
    label = list(
      intbronch ~ "Interface Bronchitis",
      plasmacellinfil ~ "Plasma Cell Infiltration",
      eosinoinfil ~ "Eosinophil Infiltration",
      lymphoidagg ~ "Lymphoid Aggregation",
      fibroelastosis ~ "Fibroelastosis",
      op ~ "Organized Pneumonia"
    ),
  ) %>%
  add_p(test.args = all_tests("fisher.test") ~ list(simulate.p.value = TRUE, B = 10000)) %>%
  bold_p() %>%
  add_overall() %>%
  separate_p_footnotes() %>%
  modify_header(label = "**Factor**") %>%
  as_gt() %>%
  gt::gtsave("brcicstats_table.png", zoom = 10)

## gtsummary fukuoka ----------------------

no3fukuoka %>%
  select(intbronch:op) %>%
  colnames() %>%
  map(get_expected_values, data = no3fukuoka, by = "maligvsbenign")

no3fukuoka %>%
  select(intbronch:op, maligvsbenign) %>%
  mutate(maligvsbenign = if_else(maligvsbenign == 1, "Malignant", "Benign")) %>%
  tbl_summary(
    by = maligvsbenign,
    missing = "ifany",
    missing_text = "Missing",
    label = list(
      intbronch ~ "Interface Bronchitis",
      plasmacellinfil ~ "Plasma Cell Infiltration",
      eosinoinfil ~ "Eosinophil Infiltration",
      lymphoidagg ~ "Lymphoid Aggregation",
      fibroelastosis ~ "Fibroelastosis",
      op ~ "Organized Pneumonia"
    ),
  ) %>%
  add_p(test.args = all_tests("fisher.test") ~ list(simulate.p.value = TRUE, B = 10000)) %>%
  bold_p() %>%
  add_overall() %>%
  separate_p_footnotes() %>%
  modify_header(label = "**Factor**") %>%
  as_gt() %>%
  gt::gtsave("fukuokastats_table.png", zoom = 10)


## gtsummary merged -------------

# vertically aligned table
tbl_merge(
  tbls = list(
    no3sano %>%
      select(intbronch:op, maligvsbenign) %>%
      mutate(maligvsbenign = if_else(maligvsbenign == 1, "Malignant", "Benign")) %>%
      tbl_summary(
        by = maligvsbenign,
        missing = "ifany",
        missing_text = "Missing",
        label = list(
          intbronch ~ "Interface Bronchitis",
          plasmacellinfil ~ "Plasma Cell Infiltration",
          eosinoinfil ~ "Eosinophil Infiltration",
          lymphoidagg ~ "Lymphoid Aggregation",
          fibroelastosis ~ "Fibroelastosis",
          op ~ "Organized Pneumonia"
        ),
      ) %>%
      add_p(test.args = all_tests("fisher.test") ~ list(simulate.p.value = TRUE, B = 10000)) %>%
      bold_p() %>%
      add_overall() %>%
      separate_p_footnotes() %>%
      modify_header(label = "**Factor**"),
    no3tachibana %>%
      select(intbronch:op, maligvsbenign) %>%
      mutate(maligvsbenign = if_else(maligvsbenign == 1, "Malignant", "Benign")) %>%
      tbl_summary(
        by = maligvsbenign,
        missing = "ifany",
        missing_text = "Missing",
        label = list(
          intbronch ~ "Interface Bronchitis",
          plasmacellinfil ~ "Plasma Cell Infiltration",
          eosinoinfil ~ "Eosinophil Infiltration",
          lymphoidagg ~ "Lymphoid Aggregation",
          fibroelastosis ~ "Fibroelastosis",
          op ~ "Organized Pneumonia"
        ),
      ) %>%
      add_p(test.args = all_tests("fisher.test") ~ list(simulate.p.value = TRUE, B = 10000)) %>%
      bold_p() %>%
      add_overall() %>%
      separate_p_footnotes() %>%
      modify_header(label = "**Factor**"),
    no3brcic %>%
      select(intbronch:op, maligvsbenign) %>%
      mutate(maligvsbenign = if_else(maligvsbenign == 1, "Malignant", "Benign")) %>%
      tbl_summary(
        by = maligvsbenign,
        missing = "ifany",
        missing_text = "Missing",
        label = list(
          intbronch ~ "Interface Bronchitis",
          plasmacellinfil ~ "Plasma Cell Infiltration",
          eosinoinfil ~ "Eosinophil Infiltration",
          lymphoidagg ~ "Lymphoid Aggregation",
          fibroelastosis ~ "Fibroelastosis",
          op ~ "Organized Pneumonia"
        ),
      ) %>%
      add_p(test.args = all_tests("fisher.test") ~ list(simulate.p.value = TRUE, B = 10000)) %>%
      bold_p() %>%
      add_overall() %>%
      separate_p_footnotes() %>%
      modify_header(label = "**Factor**"),
    no3fukuoka %>%
      select(intbronch:op, maligvsbenign) %>%
      mutate(maligvsbenign = if_else(maligvsbenign == 1, "Malignant", "Benign")) %>%
      tbl_summary(
        by = maligvsbenign,
        missing = "ifany",
        missing_text = "Missing",
        label = list(
          intbronch ~ "Interface Bronchitis",
          plasmacellinfil ~ "Plasma Cell Infiltration",
          eosinoinfil ~ "Eosinophil Infiltration",
          lymphoidagg ~ "Lymphoid Aggregation",
          fibroelastosis ~ "Fibroelastosis",
          op ~ "Organized Pneumonia"
        ),
      ) %>%
      add_p(test.args = all_tests("fisher.test") ~ list(simulate.p.value = TRUE, B = 10000)) %>%
      bold_p() %>%
      add_overall() %>%
      separate_p_footnotes() %>%
      modify_header(label = "**Factor**"),
    no3consensus3 %>%
      select(intbronch:op, maligvsbenign) %>%
      mutate(maligvsbenign = if_else(maligvsbenign == 1, "Malignant", "Benign")) %>%
      tbl_summary(
        by = maligvsbenign,
        missing = "ifany",
        missing_text = "Missing",
        label = list(
          intbronch ~ "Interface Bronchitis",
          plasmacellinfil ~ "Plasma Cell Infiltration",
          eosinoinfil ~ "Eosinophil Infiltration",
          lymphoidagg ~ "Lymphoid Aggregation",
          fibroelastosis ~ "Fibroelastosis",
          op ~ "Organized Pneumonia"
        ),
      ) %>%
      add_p(test.args = all_tests("fisher.test") ~ list(simulate.p.value = TRUE, B = 10000)) %>%
      bold_p() %>%
      add_overall() %>%
      separate_p_footnotes() %>%
      modify_header(label = "**Factor**"),
    no3consensus4 %>%
      select(intbronch:op, maligvsbenign) %>%
      mutate(maligvsbenign = if_else(maligvsbenign == 1, "Malignant", "Benign")) %>%
      tbl_summary(
        by = maligvsbenign,
        missing = "ifany",
        missing_text = "Missing",
        label = list(
          intbronch ~ "Interface Bronchitis",
          plasmacellinfil ~ "Plasma Cell Infiltration",
          eosinoinfil ~ "Eosinophil Infiltration",
          lymphoidagg ~ "Lymphoid Aggregation",
          fibroelastosis ~ "Fibroelastosis",
          op ~ "Organized Pneumonia"
        ),
      ) %>%
      add_p(test.args = all_tests("fisher.test") ~ list(simulate.p.value = TRUE, B = 10000)) %>%
      bold_p() %>%
      add_overall() %>%
      separate_p_footnotes() %>%
      modify_header(label = "**Factor**")
  ),
  tab_spanner = c("**Dr. Sano**", 
                  "**Dr. Tachibana**", 
                  "**Dr. Brcic**", 
                  "**Dr.Fukuoka**", 
                  "**3-Way Consensus**", 
                  "**4-Way Consensus**")
) %>%
  # as_flex_table() %>% # use this section if you want to save as .docx
  # flextable::save_as_docx(path = "test.docx",
  #                         pr_section = officer::prop_section(
  #                           page_size = officer::page_size(width = 20,
  #                                                          height = 20,
  #                                                          orient = "portrait")
  #                         ))
  as_gt() %>% # uncomment this section if you want to save as .png 
  tab_header(title = md("**Histological Scoring by Pathologist**")) %>% 
  tab_options(
    data_row.padding = 3,
    data_row.padding.horizontal = 0,
    heading.title.font.size = 36,
    heading.padding = 10) %>% 
  gtsave("ts_stats_merged_table_tall.png",
    zoom = 5,
    vheight = 5000,
    vwidth = 4500,
    delay = 0.5,
    expand = 10
  ) # sometimes the table clips off at the bottom when saving
# i added the above delay and expand and it doesnt happen anymore

# horizontally aligned -- just playing around with different formats
t1 <- no3sano %>%
  pivot_longer(cols = intbronch:op, names_to = "Finding", values_to = "Grade") %>%
  select(Finding, Grade, maligvsbenign) %>%
  mutate(maligvsbenign = if_else(maligvsbenign == 1, "Malignant", "Benign")) %>%
  tbl_strata(
    strata = Finding,
    .tbl_fun =
      ~ .x %>%
        tbl_summary(
          by = maligvsbenign,
          missing = "ifany",
          missing_text = "Missing"
        ) %>%
        add_p(test.args = all_tests("fisher.test") ~ list(
          simulate.p.value = TRUE,
          B = 10000
        )) %>%
        bold_p()
  )

t2 <- no3tachibana %>%
  pivot_longer(cols = intbronch:op, names_to = "Finding", values_to = "Grade") %>%
  select(Finding, Grade, maligvsbenign) %>%
  mutate(maligvsbenign = if_else(maligvsbenign == 1, "Malignant", "Benign")) %>%
  tbl_strata(
    strata = Finding,
    .tbl_fun =
      ~ .x %>%
        tbl_summary(
          by = maligvsbenign,
          missing = "ifany",
          missing_text = "Missing"
        ) %>%
        add_p(test.args = all_tests("fisher.test") ~ list(
          simulate.p.value = TRUE,
          B = 10000
        )) %>%
        bold_p()
  )

t3 <- tbl_stack(
  tbls = list(t1, t2),
  group_header = c("Dr. Sano", "Dr. Tachibana")
) %>%
  as_gt()

updated_spanners <- list(
  gt::md("**Eosinophil Infiltration**"),
  gt::md("**Fibroelastosis**"),
  gt::md("**Interface Bronchitis**"),
  gt::md("**Lymphoid Aggregation**"),
  gt::md("**Organized Pneumonia**"),
  gt::md("**Plasma Cell Infiltration**")
)
t3[["_spanners"]][["spanner_label"]] <- updated_spanners

gt::gtsave(t3, "ts_stats_merged_table_wide.png",
  vwidth = 4500,
  vheight = 4500,
  zoom = 5,
  delay = 0.5,
  expand = 50
)

## weighted kappa ------------------

# this organizes all the scoring data neatly
kappamatrix <-
  left_join(no3tachibana, no3sano, by = "ID", suffix = c(".t", ".s"), keep = FALSE) %>%
  left_join(no3brcic, by = "ID", keep = FALSE) %>%
  rename_with(~ paste0(.x,".b"), intbronch:op) %>% 
  left_join(no3fukuoka, by = "ID", keep = FALSE) %>%
  rename_with(~ paste0(.x,".f"), intbronch:op) %>% 
  select(where(is.ordered)) %>%
  mutate(
    target = NULL,
    target.t = NULL,
    target.s = NULL,
    target.b = NULL,
    target.f = NULL,
    bronch_epithelialcells = NULL
  )

# list of pathologist pairs to iterate through
pairlist <- list(
  c(".s", ".t"),
  c(".s", ".b"),
  c(".s", ".f"),
  c(".t", ".b"),
  c(".t", ".f"),
  c(".b", ".f")
)

# creates the tibble of results from the kappa calculation of each pair for each finding
cohenkapparesults <- 
  tribble(~finding,
          "intbronch",
          "plasmacellinfil",
          "eosinoinfil",
          "lymphoidagg",
          "fibroelastosis",
          "op") %>%
  mutate(kappa =
           {map(.$finding, \(finding)
                kappamatrix %>%
                  select(starts_with(finding)) %>% 
                  {map(pairlist, \(pair) select(., ends_with(pair)))} %>% 
                  {map(., kappa2, "squared")}
           )}
  ) %>% 
  unnest_longer(col = kappa) %>% 
  mutate(value =
           {
             map(.$kappa, \(row)
                 row[[5]])
           }) %>% 
  mutate(`p-value` =
           {
             map(.$kappa, \(row)
                 row[[8]])
           }) %>% 
  rowwise() %>% 
  mutate(value = round(value, digits = 3),
         `p-value` = case_when(`p-value`!=0 ~ format(`p-value`, scientific = TRUE, digits = 3))
  ) %>%
  select(!kappa) %>%
  bind_cols(rep(c("S-T",
                  "S-B",
                  "S-F",
                  "T-B",
                  "T-F",
                  "B-F"), 6)) %>%
  rename(pair = `...4`)

# create a gt
cohenkapparesults %>%
  # pivot_wider(names_from = finding, # this commented out block was from when the findings were the col spanners
  #             values_from = c(value, p),
  #             names_glue = "{finding}_{.value}",
  #             ) %>% 
  # relocate(intbronch_p, .after = intbronch_value) %>% 
  # relocate(plasmacellinfil_p, .after = plasmacellinfil_value) %>% 
  # relocate(eosinoinfil_p, .after = eosinoinfil_value) %>% 
  # relocate(lymphoidagg_p, .after = lymphoidagg_value) %>% 
  # relocate(fibroelastosis_p, .after = fibroelastosis_value) %>% 
  # relocate(op_p, .after = op_value) %>% 
  pivot_wider(names_from = pair,
              values_from = c(value, `p-value`),
              names_glue = "{pair}_{.value}") %>% 
  select(order(colnames(.), decreasing = TRUE)) %>% 
  group_by(finding) %>% 
  gt(row_group_as_column = TRUE) %>% 
  tab_spanner_delim(delim = "_") %>% 
  text_replace(locations = cells_row_groups(),
               pattern = "intbronch",
               replacement = "Interface Bronchitis") %>% 
  text_replace(locations = cells_row_groups(),
               pattern = "plasmacellinfil",
               replacement = "Plasma Cell Infiltration") %>% 
  text_replace(locations = cells_row_groups(),
               pattern = "eosinoinfil",
               replacement = "Eosinofil Infiltration") %>% 
  text_replace(locations = cells_row_groups(),
               pattern = "lymphoidagg",
               replacement = "Lymphoid Aggregation") %>% 
  text_replace(locations = cells_row_groups(),
               pattern = "fibroelastosis",
               replacement = "Fibroelastosis") %>% 
  text_replace(locations = cells_row_groups(),
               pattern = "op",
               replacement = "OP") %>% 
  text_replace(locations = cells_body(),
               pattern = "0e+00",
               replacement = "0") %>% 
  sub_missing(missing_text = "<1e-20") %>% 
  tab_style(cell_text(align = "center"),
            locations = cells_column_labels()) %>% 
  tab_options(data_row.padding.horizontal = 10) %>% 
  tab_style(style = cell_text(weight = "bold"),
            locations = cells_column_spanners()) %>% 
  tab_style(style = cell_text(weight = "bold"),
            locations = cells_column_labels()) %>% 
  tab_style(style =cell_text(align = "right"),
            locations = cells_body()) %>% 
  tab_header(title = "Weighted Kappa Values Between Each Pathologist Pair") %>% 
  tab_options(heading.padding = 15,
              source_notes.padding = 5) %>% 
  tab_source_note(source_note = "<0, no agreement; 0-0.20, 
                  slight agreement; 0.21-0.40, 
                  fair agreement; 0.41-0.60, 
                  moderate agreement; 0.61-0.80, 
                  substantial agreement; 0.81-1.0, 
                  nearly perfect agreement") %>% 
  gtsave("kappatable.png",
         zoom = 10,
         delay = 0.5,
         vwidth = 4000,
         expand = 10)

# this is experimenting with the vcd::Kappa function as opposed to the irr:kappa2 function
# cohenkapparesults2 <- 
#   tribble(~finding,
#           "intbronch",
#           "plasmacellinfil",
#           "eosinoinfil",
#           "lymphoidagg",
#           "fibroelastosis",
#           "op") %>%
#   mutate(kappa =
#            {map(.$finding, \(finding)
#                 kappamatrix %>%
#                   select(starts_with(finding)) %>% 
#                   {map(pairlist, \(pair) select(., ends_with(pair)))} %>% 
#                   {map(., table)} %>%
#                   {map(., vcd::Kappa)} %>%
#                   {map(., print)}
#            )}
#   )         

## fukuoka ROC analysis ----------------------

AUC_analysis <- function(dat) {
  roc_list <- list(
    "Intbronch" = roc(dat$maligvsbenign, dat$intbronch,
                      direction = ">"),
    "Plasma Cell Infil" = roc(dat$maligvsbenign, dat$plasmacellinfil,
                              direction = ">"),
    "Eosinophil Infil" = roc(dat$maligvsbenign, dat$eosinoinfil,
                             direction = ">"),
    "Lymphoid Agg" = roc(dat$maligvsbenign, dat$lymphoidagg,
                         direction = ">"),
    "Fibroelastosis" = roc(dat$maligvsbenign, dat$fibroelastosis,
                           direction = "<"),
    "OP" = roc(dat$maligvsbenign, dat$op,
               direction = ">")
  )
  
  auc_list <- 
    map(roc_list, \(x) auc(x)) %>% 
    tibble("name" = .) %>% 
    mutate(auc = paste0(round(as.numeric(name), digits = 3))) %>% 
    mutate(name = c("Intbronch",
                    "Plasma Cell Infil",
                    "Eosinophil Infil",
                    "Lymphoid Agg",
                    "Fibroelastosis",
                    "OP"))
  
  ggroc(roc_list,
        show.legend = FALSE) +
    facet_wrap(~name) +
    theme_minimal() +
    geom_abline(slope = 1,
                intercept = 1) +
    geom_text(data = auc_list,
              aes(
                  label = auc,
                  0,0,
                  vjust = 0,
                  hjust = 1
              ),
              show.legend = FALSE)
}
AUC_analysis(no3fukuoka)
ggsave("fukuoka_roc.png",
         dpi = 320,
         width = 8,
         height = 6)

## sano ROC analysis ----------------
AUC_analysis(no3sano)
ggsave("sano_auc.png",
       dpi = 320,
       height = 6,
       width = 8)
## consensus analysis -------------

no3consensus3 %>%
  select(intbronch:op) %>%
  colnames() %>%
  map(get_expected_values, data = no3consensus3, by = "maligvsbenign")

no3consensus3 %>%
  select(intbronch:op, maligvsbenign) %>%
  mutate(maligvsbenign = if_else(maligvsbenign == 1, "Malignant", "Benign")) %>%
  tbl_summary(
    by = maligvsbenign,
    missing = "ifany",
    missing_text = "Missing",
    label = list(
      intbronch ~ "Interface Bronchitis",
      plasmacellinfil ~ "Plasma Cell Infiltration",
      eosinoinfil ~ "Eosinophil Infiltration",
      lymphoidagg ~ "Lymphoid Aggregation",
      fibroelastosis ~ "Fibroelastosis",
      op ~ "Organized Pneumonia"
    ),
  ) %>%
  add_p(test.args = all_tests("fisher.test") ~ list(simulate.p.value = TRUE, B = 10000)) %>%
  bold_p() %>%
  add_overall() %>%
  separate_p_footnotes() %>%
  modify_header(label = "**Factor**") %>%
  as_gt() %>%
  gt::gtsave("consensus3stats_table.png", zoom = 10)

no3brcic %>%
  select(intbronch:op) %>%
  colnames() %>%
  map(get_expected_values, data = no3brcic, by = "maligvsbenign")

no3consensus4 %>%
  select(intbronch:op, maligvsbenign) %>%
  mutate(maligvsbenign = if_else(maligvsbenign == 1, "Malignant", "Benign")) %>%
  tbl_summary(
    by = maligvsbenign,
    missing = "ifany",
    missing_text = "Missing",
    label = list(
      intbronch ~ "Interface Bronchitis",
      plasmacellinfil ~ "Plasma Cell Infiltration",
      eosinoinfil ~ "Eosinophil Infiltration",
      lymphoidagg ~ "Lymphoid Aggregation",
      fibroelastosis ~ "Fibroelastosis",
      op ~ "Organized Pneumonia"
    ),
  ) %>%
  add_p(test.args = all_tests("fisher.test") ~ list(simulate.p.value = TRUE, B = 10000)) %>%
  bold_p() %>%
  add_overall() %>%
  separate_p_footnotes() %>%
  modify_header(label = "**Factor**") %>%
  as_gt() %>%
  gt::gtsave("consensus4stats_table.png", zoom = 10)

AUC_analysis(no3consensus3)
ggsave("consensus3_auc.png",
       dpi = 320,
       height = 6,
       width = 8)
AUC_analysis(no3consensus4)
ggsave("consensus4_auc.png",
       dpi = 320,
       height = 6,
       width = 8)