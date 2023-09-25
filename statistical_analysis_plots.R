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
# Email: ethanokoshi+gh <at> gmail <dot> com
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
library(corrplot)
#
## load data -------------
#
input_data_path <- "~/ethancoding/Sano/sanoproject_ENO.xlsx"
library(readxl)
data <- read_excel(input_data_path, sheet = "fukuoka")

no3data <-
  data[data$target != 3, ] %>%
  mutate(across(intbronch:lymphoidagg, na_if, "NA"),
         across(intbronch:op, factor,
                ordered = TRUE
         ),
         maligvsbenign = if_else(target == 1 | target == 2, 1, 0),
         maligvsbenign = as.factor(maligvsbenign)
  ) %>%
  select(!target)

## End Header --------------


## gtsummary clinical --------------

no3data %>%
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
  as_gt()
  gt::gtsave("clinicalstats_table.png", zoom = 10)
## fig. 4 AUC curves ------------------

AUC_analysis <- function(dat) {
  roc_list <- list(
    "IB/B" = roc(dat$maligvsbenign, dat$intbronch,
                 direction = ">"),
    "PLC" = roc(dat$maligvsbenign, dat$plasmacellinfil,
                direction = ">"),
    "Eo" = roc(dat$maligvsbenign, dat$eosinoinfil,
               direction = ">"),
    "Ly" = roc(dat$maligvsbenign, dat$lymphoidagg,
               direction = ">"),
    "FE" = roc(dat$maligvsbenign, dat$fibroelastosis,
               direction = "<"),
    "OP" = roc(dat$maligvsbenign, dat$op,
               direction = ">")
  )
  
  auc_list <- 
    map(roc_list, \(x) auc(x)) %>% 
    tibble("name" = .) %>% 
    mutate(auc = paste0(round(as.numeric(name), digits = 3))) %>% 
    mutate(name = c("IB/B",
                    "PLC",
                    "Eo",
                    "Ly",
                    "FE",
                    "OP"))
  
  ggroc(roc_list,
        show.legend = FALSE,
        # color = "grey20",
        linewidth = 1
  ) +
    facet_wrap(~name) +
    geom_abline(slope = 1,
                intercept = 1,
    ) +
    geom_text(data = auc_list,
              aes(
                label = paste0("AUC = ",auc),
                0,0,
                vjust = 0,
                hjust = 1
              ),
              show.legend = FALSE,
              color = "grey20",
              family = "Barlow Medium") +
    # labs(caption = "Figure 4. Receiver Operating Curve of the Scored Histological Features") +
    xlab("Specificity") +
    ylab("Sensitivity") +
    scale_x_continuous(trans = "reverse",
                       breaks = c(1,0.5,0)) +
    scale_y_continuous(breaks = c(1,0.5,0)) +
    theme_bw() +
    theme(plot.caption = element_text(size = 12, hjust = 0),
          text = element_text(family = "Barlow Medium"),
          axis.title.y.left = element_text(margin = margin(0,8,0,0, unit = "pt"))
    )
}
AUC_analysis(no3data)

## fig. 4 bar plots

fig4barplotdata <- 
  select(no3data, intbronch:op, maligvsbenign) %>% 
  mutate(across(intbronch:op, as.character)) %>%
  rename("IB/B" = intbronch,
         "Eo" = eosinoinfil,
         "Ly" = lymphoidagg,
         "PLC" = plasmacellinfil,
         "FE" = fibroelastosis,
         "OP" = op) %>% 
  pivot_longer(cols = !maligvsbenign, names_to = "feature", values_to = "value") %>% 
  mutate(across(value, as.factor)) %>% 
  drop_na()

barplots <-
  ggplot(data = fig4barplotdata, aes(x = maligvsbenign, fill = value)) +
  facet_wrap("feature", scales = "free_x") +
  geom_bar(position = "fill") +
  scale_x_discrete(labels = c("Benign", "Malignant")) +
  scale_fill_brewer("Score", palette = "Accent") +
  scale_y_continuous(limits = c(0,1.15), breaks = c(0,0.25,0.5,0.75,1)) +
  guides(fill = guide_legend(reverse = TRUE)) +
  labs(x = NULL,
       y = "Score Frequency") +
  theme_gray() +
  theme(text = element_text(family = "Barlow Medium"),
        axis.title.y.left = element_text(margin = margin(0,8,0,0, unit = "pt"))
  )

rocplots <- 
  AUC_analysis(no3data)

fig4 <- gridExtra::grid.arrange(barplots, rocplots, ncol=1)

ggsave("fig4.png",
       plot = fig4,
       dpi = 320,
       width = 5,
       height = 8)

## corrplot -----------------

test <- no3data %>% 
  select(intbronch,plasmacellinfil,eosinoinfil,lymphoidagg,fibroelastosis,op) %>% 
  mutate(across(everything(), as.numeric)) %>% 
  filter(!is.na(intbronch)) %>% 
  cor.mtest(method = "kendall")

no3data %>% 
  select(intbronch,plasmacellinfil,eosinoinfil,lymphoidagg,fibroelastosis,op) %>% 
  mutate(across(everything(), as.numeric)) %>% 
  filter(!is.na(intbronch)) %>% 
  cor(method = "kendall") %>% 
  corrplot(p.mat = test$p, insig = "p-value", sig.level = -1)

## Supplemental fig 1 corr plot -----------------


test <- no3data %>% 
  select(intbronch,plasmacellinfil,eosinoinfil,lymphoidagg,fibroelastosis,op) %>% 
  mutate(across(everything(), as.numeric)) %>% 
  filter(!is.na(intbronch)) %>% 
  rename("IB/B" = intbronch,
         "Eo" = eosinoinfil,
         "Ly" = lymphoidagg,
         "PLC" = plasmacellinfil,
         "FE" = fibroelastosis,
         "OP" = op) %>% 
  cor.mtest(method = "kendall")

no3data %>% 
  select(intbronch,plasmacellinfil,eosinoinfil,lymphoidagg,fibroelastosis,op) %>% 
  mutate(across(everything(), as.numeric)) %>% 
  filter(!is.na(intbronch)) %>% 
  rename("IB/B" = intbronch,
         "Eo" = eosinoinfil,
         "Ly" = lymphoidagg,
         "PLC" = plasmacellinfil,
         "FE" = fibroelastosis,
         "OP" = op) %>% 
  cor(method = "kendall") %>% 
  ggcorrplot::ggcorrplot(method = "square", 
                         p.mat = test$p,
                         insig = "blank",
                         hc.order = TRUE,
                         outline.color = "white",
                         lab = T,
                         lab_size = 4) +
  theme(legend.direction = "horizontal",
        legend.position = "top")
ggsave("suppfig1.jpeg", dpi = 1000, height = 12, width = 12, units = "cm")

## Supplemental table 1 classifier comparison -------------
classifier_compare_metrics <- read_csv("~/Sano_TBLB/classifier_compare_metrics.csv")
classifier_compare_metrics %>% 
  rename(
    "Algorithm" = `0`,
    "Accuracy" = `1`,
    "AUC" = `2`,
    "F1" = `3`,
    "Precision" = `4`,
    "Recall" = `5`,
    "Time" = `6`,
  ) %>% 
  group_by(Algorithm) %>% 
  summarise(across(everything(), mean)) %>% 
  gt() %>% 
  fmt_number(decimals = 3) %>% 
  gtsave("supptable1.docx")
