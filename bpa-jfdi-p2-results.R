#' Main functions for generating plots for the BLOODPAC JFDI
#' publication.
require(tidyverse)
require(ggpubr)
require(combinat)

source("utils/jfdi-utils.R")

# Commandline options and parsing
args <- commandArgs(trailingOnly = TRUE)

# validate
if(length(args) != 1) {
  stop("Usage: Rscript bpa-jfdi-p2-results.R <input-data.csv>\nYou must provide the path to the input CSV.")
}

input.csv <- args[1]
if(!file.exists(input.csv)) {
  stop(paste("Unable to find input CSV:", input.csv, ". Make sure file exists."))
}

## Load dataset
message(paste("Loading data set", input.csv, "..."))
jfdi.detect <- read.csv(input.csv, stringsAsFactors = FALSE, na.strings=c("NA", "")) %>%
  mutate(expected_af=forcats::fct_relevel(expected_af, c("WT", "0.1%", "0.5%", "1.0%", "5.0%")))

## Setup output directories
dir.create("figures/", showWarnings = FALSE, recursive = TRUE)
dir.create("csvs/", showWarnings = FALSE, recursive = TRUE)

## figures
BaseSize <- 6
PointSize <- 0.07
ColorPalette <- c("#969696","#000000")

#' do.fig02
#' Generates the figure showing boxplots of observed AF across the experiment.
#'
#' @return figure 02 from publication
do.fig02 <- function(){
  jfdi.fig02.a <- jfdi.detect %>%
    filter(!is.na(Observed) & sequencing_platform != "dPCR" & contrived_source=="Horizon") %>%
    rename(Participant=participant) %>%
    ggboxplot(x="expected_af", y="observed_af_pct", color="Participant",
              palette="jco",
              outlier.size=0.2,
              size=0.2,
              facet.by = c("dna_source", "contrived_source"),
              ggtheme = theme_pubclean(base_size=5),
              ylab="Observed VAF (%)",
              xlab="Expected VAF") +
    theme(strip.text = element_text(face = "bold"),
          strip.placement = "outside",
          strip.background=element_rect(color=NULL, fill="white"),
          legend.key.width=unit(.25, "cm"),
          legend.key.height=unit(.25, "cm"))

  jfdi.fig02.b <- jfdi.detect %>%
    filter(!is.na(Observed) & sequencing_platform != "dPCR" & contrived_source=="SeraCare") %>%
    rename(Participant=participant) %>%
    ggboxplot(x="expected_af", y="observed_af_pct", color="Participant",
              palette="jco",
              outlier.size=0.2,
              size=0.2,
              facet.by = c("dna_source", "contrived_source"),
              ggtheme = theme_pubclean(base_size=5),
              ylab="Observed VAF (%)",
              xlab="Expected VAF") +
    theme(strip.text = element_text(face = "bold"),
          strip.placement = "outside",
          strip.background=element_rect(color=NULL, fill="white"),
          legend.key.width=unit(.25, "cm"),
          legend.key.height=unit(.25, "cm"))
  fig <- ggarrange(
    jfdi.fig02.a,
    jfdi.fig02.b,
    ncol=2,
    common.legend = TRUE
  )
  return(fig)
}

#' do.fig03
#' Subsets the dataframe to variants covered by both plasma and dna formats within a
#' participant and is covered by at least 2 participants.
#'
#' @return figure 03
do.fig03 <- function(){
  sel <- jfdi.detect %>%
    filter(!(is.na(Observed)) & sequencing_platform != "dPCR") %>%
    distinct(contrived_source, dna_source, expected_af, participant) %>%
    count(contrived_source, expected_af, participant) %>%
    filter(n==2) %>%
    select(-n)

  sel.plasma <- jfdi.detect %>%
    filter(!(is.na(Observed)) & sequencing_platform != "dPCR" & dna_source == "Plasma") %>%
    inner_join(sel) %>%
    distinct(contrived_source, expected_af, variant, participant) %>%
    count(contrived_source, expected_af, variant)

  sel.dna_source <- jfdi.detect %>%
    filter(!(is.na(Observed)) & sequencing_platform != "dPCR" & dna_source == "DNA") %>%
    inner_join(sel) %>%
    distinct(contrived_source, expected_af, variant, participant) %>%
    count(contrived_source, expected_af, variant) %>%
    inner_join(sel.plasma) %>%
    ungroup() %>%
    filter(n>2)

  jfdi.fig03.a <- jfdi.detect %>%
    filter(!(is.na(Observed)) & sequencing_platform != "dPCR" & contrived_source == "Horizon") %>%
    inner_join(sel) %>%
    inner_join(sel.dna_source) %>%
    rename(`DNA Format`=dna_source) %>%
    ggbarplot(x="variant", y="observed_af_pct", fill="DNA Format",
              palette=ColorPalette, facet.by=c("expected_af", "contrived_source"),
              scales="free",
              add="mean_se",
              add.params=list(size=.2),
              size=.2,
              position=position_dodge(0.8),
              ggtheme=theme_pubclean(base_size=BaseSize),
              ylab=FALSE,
              xlab=FALSE) +
    rotate_x_text(angle=45) +
    theme(strip.text = element_text(face = "bold"),
          strip.placement = "outside",
          strip.background=element_rect(color=NULL, fill="white"),
          legend.key.width=unit(.25, "cm"),
          legend.key.height=unit(.25, "cm"))

  jfdi.fig03.b <- jfdi.detect %>%
    filter(!(is.na(Observed)) & sequencing_platform != "dPCR" & contrived_source == "SeraCare") %>%
    inner_join(sel) %>%
    inner_join(sel.dna_source) %>%
    rename(`DNA Format`=dna_source) %>%
    ggbarplot(x="variant", y="observed_af_pct", fill="DNA Format",
              palette=ColorPalette, facet.by=c("expected_af", "contrived_source"),
              scales="free",
              add="mean_se",
              add.params=list(size=.2),
              size=.2,
              position=position_dodge(0.8),
              ggtheme=theme_pubclean(base_size=BaseSize),
              ylab=FALSE,
              xlab=FALSE) +
    rotate_x_text(angle=45) +
    theme(strip.text = element_text(face = "bold"),
          strip.placement = "outside",
          strip.background=element_rect(color=NULL, fill="white"),
          legend.key.width=unit(.25, "cm"),
          legend.key.height=unit(.25, "cm"))

  ggarrange(
    jfdi.fig03.a,
    jfdi.fig03.b,
    ncol=2,
    common.legend = TRUE,
    widths = c(0.4, 0.6)
  ) %>% annotate_figure(left=text_grob("Observed VAF (%)", rot=90, size = 8))
}

#' do.fig04
#' Loads the sensitivity and specificity data frames; Fills in some missing data as NAs for
#' participant c; generates plots
#'
#' @return Figure 04
do.fig04 <- function(){
  jfdi.sensitivity <- get.sensitivity.dataframe(jfdi.detect %>% filter(sequencing_platform != "dPCR"))
  jfdi.specificity <- get.specificity.dataframe(jfdi.detect %>% filter(sequencing_platform != "dPCR"))

  ## Group by sens and spec, seracare plasma/seracare dna etc.
  c.meta <- jfdi.sensitivity %>%
    filter(participant == "c") %>%
    distinct(participant, sequencing_platform, library_method, contrived_source, dna_source)
  c.meta$expected_af <- "WT"
  c.meta$sample_replicate <- NA
  c.meta$total <- NA
  c.meta$TN <- NA
  c.meta$FP <- NA
  c.meta$FPR <- NA
  c.meta$TNR <- NA

  jfdi.specificity.2 <- rbind(jfdi.specificity, c.meta)

  jfdi.fig04.a <- jfdi.sensitivity %>%
    filter(contrived_source == "Horizon") %>%
    rename(`DNA Format`=dna_source) %>%
    mutate(sample_replicate=as.factor(sample_replicate)) %>%
    ggline(x="expected_af", y="TPR", linetype="DNA Format",
           ggtheme=theme_pubclean(base_size=BaseSize),
           facet.by = "participant",
           point.size=PointSize,
           shape="DNA Format",
           size=.3,
           add="mean_se",
           ylab="Sensitivity", xlab=FALSE,
           nrow=1) +
    theme(strip.text = element_text(face = "bold"),
          strip.placement = "outside",
          strip.background=element_rect(color=NULL, fill="white"),
          legend.key.width=unit(.25, "cm"),
          legend.key.height=unit(.25, "cm"))

  jfdi.fig04.b <- jfdi.specificity.2 %>%
    filter(contrived_source == "Horizon") %>%
    rename(`DNA Format`=dna_source) %>%
    ggbarplot(x='DNA Format', y='TNR', fill='DNA Format',
              palette=ColorPalette,
              ggtheme=theme_pubclean(base_size=BaseSize),
              size=.2,
              facet.by = "participant",
              add="mean_se", nrow=1, ylab="Specificity", xlab=FALSE) +
    theme(strip.text = element_text(face = "bold"),
          strip.placement = "outside",
          strip.background=element_rect(color=NULL, fill="white"),
          legend.key.width=unit(.25, "cm"),
          legend.key.height=unit(.25, "cm"))

  jfdi.fig04.c <- jfdi.sensitivity %>%
    filter(contrived_source == "SeraCare") %>%
    rename(`DNA Format`=dna_source) %>%
    mutate(sample_replicate=as.factor(sample_replicate)) %>%
    ggline(x="expected_af", y="TPR", linetype="DNA Format",
           palette=ColorPalette,
           ggtheme=theme_pubclean(base_size=BaseSize),
           size=.3,
           facet.by = "participant",
           point.size=PointSize,
           shape="DNA Format",
           add="mean_se",
           ylab="Sensitivity", xlab=FALSE,
           scales="free_x",
           nrow=1) +
    theme(strip.text = element_text(face = "bold"),
          strip.placement = "outside",
          strip.background=element_rect(color=NULL, fill="white"),
          legend.key.width=unit(.25, "cm"),
          legend.key.height=unit(.25, "cm"))

  jfdi.fig04.d <- jfdi.specificity.2 %>%
    filter(contrived_source == "SeraCare") %>%
    rename(`DNA Format`=dna_source) %>%
    ggbarplot(x='DNA Format', y='TNR', fill='DNA Format',
              palette=ColorPalette,
              ggtheme=theme_pubclean(base_size=BaseSize),
              size=.2,
              facet.by = "participant",
              add="mean_se", nrow=1, ylab="Specificity", xlab=FALSE) +
    theme(strip.text = element_text(face = "bold"),
          strip.placement = "outside",
          strip.background=element_rect(color=NULL, fill="white"),
          legend.key.width=unit(.25, "cm"),
          legend.key.height=unit(.25, "cm"))

  ggarrange(ggarrange(jfdi.fig04.a,
                      jfdi.fig04.b,
                      nrow = 2,
                      common.legend=TRUE),
            ggarrange(jfdi.fig04.c, jfdi.fig04.d, nrow=2, common.legend=TRUE),
            ncol=1,
            labels = c("Horizon", "SeraCare"),
            font.label = list(size=8),
            common.legend = TRUE)

}

#' do.fig05
#' Creates metadata df; loads w/i and b/w concordance; merges into a single df; generates plot
#'
#' @return figure 05
do.fig05 <- function(){
  # Metadata to add back
  lib.meta <- jfdi.detect %>% distinct(participant, library_method, sequencing_platform)
  jfdi.wi.concordance <- get.within.concordance.df(jfdi.detect %>%
                                                     filter(expected_af != "WT" & sequencing_platform != "dPCR")) %>%
    inner_join(lib.meta)
  jfdi.bw.concordance <- get.between.concordance.df(jfdi.detect %>%
                                                      filter(expected_af != "WT" & sequencing_platform != "dPCR"))


  jfdi.bw.concordance.fmt <- jfdi.bw.concordance %>%
    mutate(comparison="cross-participant") %>%
    inner_join(lib.meta, by=c("participant.a"="participant")) %>%
    select(contrived_source, dna_source, sequencing_platform, library_method, participant.a, participant.b, comparison, expected_af, rc)

  jfdi.wi.concordance.fmt <- jfdi.wi.concordance %>%
    mutate(participant.a=participant, participant.b=participant, comparison="intra-participant") %>%
    select(contrived_source, dna_source, sequencing_platform, library_method, participant.a, participant.b, comparison, expected_af, rc)

  jfdi.concord.combined <- rbind(jfdi.bw.concordance.fmt, jfdi.wi.concordance.fmt) %>%
    mutate(comparison=factor(comparison, levels=c("intra-participant", "cross-participant")))

  jfdi.fig.05.a <- jfdi.concord.combined %>%
    filter(contrived_source == "Horizon") %>%
    rename(Grouping=comparison) %>%
    ggline(x="expected_af", y="rc", color="black",
           linetype="Grouping",
           repel=TRUE,
           facet.by=c("dna_source", "participant.a"),
           ggtheme=theme_pubclean(base_size=BaseSize),
           point.size=PointSize,
           shape="Grouping",
           size=0.3,
           add="mean_se",
           ylab="Concordance",
           xlab=FALSE) +
    theme(strip.text = element_text(face = "bold"),
          strip.placement = "outside",
          strip.background=element_rect(color=NULL, fill="white"),
          legend.key.width=unit(.25, "cm"),
          legend.key.height=unit(.25, "cm"))

  jfdi.fig.05.b <- jfdi.concord.combined %>%
    filter(contrived_source == "SeraCare") %>%
    rename(Grouping=comparison) %>%
    ggline(x="expected_af", y="rc", color="black",
           linetype="Grouping",
           facet.by=c("dna_source", "participant.a"),
           ggtheme=theme_pubclean(base_size=BaseSize),
           point.size=PointSize,
           shape="Grouping",
           repel=TRUE,
           size=0.3,
           add="mean_se",
           ylab="Concordance",
           xlab=FALSE) +
    theme(strip.text = element_text(face = "bold"),
          strip.placement = "outside",
          strip.background=element_rect(color=NULL, fill="white"),
          legend.key.width=unit(.25, "cm"),
          legend.key.height=unit(.25, "cm"))

  ggarrange(ggarrange(jfdi.fig.05.a),
            ggarrange(jfdi.fig.05.b),
            nrow = 2,
            common.legend=TRUE,
            labels = c("Horizon", "SeraCare"),
            font.label=list(size=8))
}

# Stop pesky Rplots.pdf from being generated
if(!interactive()) pdf(NULL)
message("Creating figures/jfdi.fig02.jpeg ...")
jfdi.fig02 <- do.fig02()
ggsave("jfdi.fig02.jpeg",
       plot=jfdi.fig02,
       path="figures",
       width=8,
       height=8,
       units="cm",
       dpi=300)

message("Creating figures/jfdi.fig03.jpeg ...")
jfdi.fig03 <- do.fig03()
ggsave("jfdi.fig03.jpeg",
       plot=jfdi.fig03,
       path="figures",
       width=17,
       height=10,
       units="cm",
       dpi=300)

message("Creating figures/jfdi.fig04.jpeg ...")
jfdi.fig04 <- do.fig04()
ggsave("jfdi.fig04.jpeg",
       plot=jfdi.fig04,
       path="figures",
       width=17,
       height=10,
       units="cm",
       dpi=300)

message("Creating figures/jfdi.fig05.jpeg ...")
jfdi.fig05 <- do.fig05()
ggsave("jfdi.fig05.jpeg",
       plot=jfdi.fig05,
       path="figures",
       width=17,
       height=10,
       units="cm",
       dpi=300)
message("All figures generated.")

message("Creating csvs/jfdi.table_S3.csv ...")
jfdi.detect %>%
  dplyr::filter(!is.na(Observed) & participant != "i") %>%
  dplyr::distinct(contrived_source, dna_source, expected_af, variant, participant) %>%
  dplyr::count(contrived_source, dna_source, expected_af, participant) %>%
  tidyr::pivot_wider(names_from=expected_af, values_from = n, values_fill = NA) %>%
  dplyr::arrange(contrived_source, dna_source, participant) %>%
  write.csv(file="csvs/jfdi.table_S3.csv", row.names=FALSE, na="NA")

message("Creating csvs/jfdi.table_S4.csv ...")
jfdi.detect %>%
  dplyr::filter(!is.na(Observed) & participant != "i") %>%
  dplyr::distinct(contrived_source, dna_source, expected_af, participant) %>%
  dplyr::count(contrived_source, expected_af, dna_source) %>%
  tidyr::pivot_wider(names_from=expected_af, values_from = n, values_fill = NA) %>%
  dplyr::arrange(contrived_source, dna_source) %>%
  write.csv(file="csvs/jfdi.table_S4.csv", row.names=FALSE, na="NA")

message("Creating csvs/jfdi.table_S5.csv ...")
jfdi.detect %>%
  dplyr::filter(!is.na(Observed) & participant != "i") %>%
  dplyr::distinct(contrived_source, dna_source, expected_af, variant, participant) %>%
  dplyr::count(contrived_source, dna_source, expected_af, variant) %>%
  tidyr::pivot_wider(names_from=expected_af, values_from = n, values_fill = NA) %>%
  dplyr::arrange(contrived_source, dna_source, variant) %>%
  write.csv(file="csvs/jfdi.table_S5.csv", row.names=FALSE, na="NA")

message("Done. All supplemental csvs generated.")

print(sessionInfo())
