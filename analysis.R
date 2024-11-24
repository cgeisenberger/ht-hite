library(tidyverse)
library(minfi)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
library(conumee2.0)
library(random)
library(Rsamtools)
library(DescTools)
library(vcfR)
library(MutationalPatterns)
library(BSgenome)

ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)


source("./functions.R")


# color scales
#cols_fig1 <- RColorBrewer::brewer.pal(9, "Oranges")[4:9]
#cols_fig2 <- RColorBrewer::brewer.pal(9, "Greens")[4:9]
#cols_fig3 <- RColorBrewer::brewer.pal(9, "Purples")

colorscheme <- c("#c8553d","#588b8b","#0a2463","#e2e2e2","#f28f3b", "#798071")



# Comparison DNA for manual HiTE vs. reference ---------------------------------

manual <- readxl::read_xlsx(path = "./input/20240220_comparison_manual.xlsx")

# t-test: yield vs. incubation time
t.test(qubit ~ incubation, data = manual) # p-value = 2.028e-05

# t-test: yield vs. method for 24 hrs
manual %>% 
  filter(incubation == "24h") %>% 
  with(., anova(lm(qubit ~ method)))

# avg DNA yield 
manual %>% 
  filter(incubation == "24h") %>% 
  group_by(method) %>% 
  summarise(avg_yield = mean(qubit))

# t-test: purity vs. method for 24 hrs
manual %>% 
  filter(incubation == "24h") %>% 
  with(., anova(lm(a260_a280 ~ method)))

manual %>% 
  ggplot(aes(method, qubit, col = as.factor(method))) +
  geom_boxplot(lwd = 1.5) +
  ylim(0, 25) +
  theme_bw(base_size = 20) +
  theme(legend.position = "none", legend.title = element_blank(), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  facet_wrap(facets = vars(incubation)) +
  labs(x = NULL, y = "DNA (ng/ul)") +
  scale_color_manual(values = colorscheme)
ggsave(filename = "./plots/fig_1a.pdf", width = 7, height = 5.5)


manual %>% 
  ggplot(aes(method, a260_a280, col = as.factor(method))) +
  geom_boxplot(lwd = 1) +
  theme_bw(base_size = 20) +
  geom_hline(yintercept = 1.8, col = "grey", lty = 2) +
  theme(legend.position = "none", legend.title = element_blank(), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  facet_wrap(facets = vars(incubation)) +
  labs(x = NULL, y = "DNA Purity (A260/A280)") +
  scale_color_manual(values = colorscheme)
ggsave(filename = "./plots/fig_1b.pdf", width = 7, height = 5.5)


manual %>% 
  ggplot(aes(qubit, nanodrop,  col = cleanup, shape = cleanup)) +
  geom_point(size = 5) +
  theme_bw(base_size = 20) +
  theme(legend.position = "top", legend.title = element_blank()) +
  geom_abline(slope = 1, intercept = 0, lty = 2, col = "grey") +
  labs(x = "Qubit ng/ul", y = "Nanodrop (ng/ul)") +
  xlim(0, 25) +
  ylim(0, 80) +
  scale_color_manual(values = rep(colorscheme[1:3], each = 2))
ggsave(filename = "./plots/fig_1c.pdf", width = 7, height = 5.5)

manual %>% 
  with(., cor(qubit, nanodrop))

manual %>% 
  with(., cor(qubit, nanodrop, method = "spearman"))

manual %>% 
  mutate(fold_x = nanodrop/qubit) %>% 
  filter(cleanup != "Reference_1h") %>% 
  filter(qubit > 10) %>% 
  summarise(mean_fold = mean(fold_x), 
            median_fold = median(fold_x))

manual %>% 
  mutate(fold_x = nanodrop/qubit) %>% 
  filter(cleanup != "Reference_1h") %>%
  ggplot(aes(x = as.factor(1), y = fold_x, col = method)) +
  geom_boxplot(lwd = 2, col = "darkgrey", alpha = 0.7) +
  geom_jitter(size = 4, width = 0.2) +
  theme_bw(base_size = 24) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  scale_color_manual(values = colorscheme) +
  theme(legend.position = "none") +
  labs(x = NULL, y = "Overestimation Nanodrop (fold)")
ggsave(filename = "./plots/fig_1d.pdf", width = 3.8, height = 5.5)

manual %>% 
  mutate(fold_x = nanodrop/qubit) %>% 
  filter(cleanup != "Reference_1h") %>%
  ggplot(aes(qubit, fold_x, col = method)) +
  geom_smooth(col = "steelblue", fill = "darkgrey") +
  geom_point(size = 4) +
  theme_bw(base_size = 24) +
  theme(legend.position = "none") +
  scale_color_manual(values = colorscheme) +
  labs(x = "Qubit (ng/ul)", y = "Overestimation Nanodrop (fold)")
ggsave(filename = "./plots/fig_1e.pdf", width = 4.4, height = 5.5)



# bioanalyzer traces -----------------------------------------------------------

library(bioanalyzeR)
eph_object <- read.bioanalyzer(xml.file = "./input/bioanalyzer/Extraction Test 2.xml")

eph_anno <- eph_object$samples %>% 
  as_tibble() %>% 
  dplyr::select(-c(batch, sample.comment, ladder.well)) %>% 
  dplyr::rename("sample.index" = well.number, 
                "method" = sample.observations) %>% 
  dplyr::mutate(sample.index = as.integer(sample.index))

eph_anno <- eph_anno %>% 
  mutate(sample.id = c(rep(paste0("Sample ", c(2, 3, 1)), each = 3), "Ladder"),
         method = c(rep(c("Reference", "Beads", "Column"), 3), NA)) %>% 
  mutate(sample.name = paste0(sample.id, " ", method))

eph_data <- eph_object$data %>% as_tibble

eph_data <- left_join(eph_data, eph_anno)

  
eph_data %>% 
  filter(sample.index != 10) %>% 
  group_by(sample.index) %>% 
  ggplot(aes(time, fluorescence, col = as.factor(method))) +
  geom_line() +
  facet_grid(cols = vars(sample.id), 
             rows = vars(method)) +
  theme_bw(base_size = 20) +
  theme(legend.position = "none") +
  labs(x = "Time (s)", y = "Fluorescence (AU)") +
  scale_color_manual(values = colorscheme[1:3])
ggsave(filename = "./plots/fig_2_bioanalyzer.pdf", width = 12, height = 9)



# plate-based processing -------------------------------------------------------

plates <- readxl::read_xlsx(path = "./input/20240522_plate_workflow.xlsx")

plates <- plates %>% 
  mutate(overestimation = nanodrop/qubit)

plates %>% 
  filter(is.na(evaporation_beads)) %>% 
  dplyr::select(-c("well", "sample", "replicate", "evaporation_beads")) %>% 
  GGally::ggpairs()

# overall yield (qubit)
plates %>%
  filter(is.na(evaporation_beads)) %>% 
  filter(control == "no") %>% 
  ggplot(aes(qubit)) + 
  geom_histogram(alpha = 0.7, bins = 40, fill = colorscheme[1]) +
  theme_bw(base_size = 20) +
  labs(x = "Qubit (ng/ul)", y = "Observations (n)")
ggsave(filename = "./plots/fig_3b.pdf", width = 8.6, height = 2.8)

# overall yield (nanodrop)
plates %>%
  filter(is.na(evaporation_beads)) %>% 
  filter(control == "no") %>% 
  ggplot(aes(nanodrop)) + 
  geom_histogram(aes(y = after_stat(density)),
                 alpha = 0.7, bins = 40, fill = colorscheme[2]) +
  geom_density(lwd = 1, color = colorscheme[2]) +
  theme_bw(base_size = 20) +
  labs(x = "Nanodrop (ng/ul)", y = "Density")
ggsave(filename = "./plots/fig_s1a.pdf", width = 8.6, height = 2.8)

# mean and median yield
plates %>%
  filter(is.na(evaporation_beads)) %>% 
  filter(control == "no") %>% 
  summarise(mean_qubit = mean(qubit, na.rm = TRUE), 
            median_qubit = median(qubit, na.rm = TRUE),
            mean_nanodrop = mean(nanodrop, na.rm = TRUE), 
            median_nanodrop = median(nanodrop, na.rm = TRUE))
# mean qubit = 86.6
# median = 72

# difference in concentration
plates %>% 
  filter(is.na(evaporation_beads)) %>% 
  filter(control == "no") %>% 
  pivot_longer(cols = c(5, 6), names_to = "measurement", values_to = "conc") %>% 
  mutate(measurement = str_to_title(measurement), 
         method = str_to_title(method)) %>% 
  ggplot(aes(method, conc, col = method)) +
  geom_boxplot(lwd = 1) +
  theme_bw(base_size = 24) +
  theme(legend.position = "none") +
  facet_wrap(facets = vars(measurement), scales = "free") +
  labs(x = NULL, y = "DNA yield (ng/ul)") +
  scale_color_manual(values = colorscheme[c(1, 2)])
ggsave(filename = "./plots/fig_3c.pdf", width = 8.6, height = 4.3)


# samples vs. controls (qubit)
plates %>%
  filter(is.na(evaporation_beads)) %>% 
  with(., t.test(qubit ~ control)$p.value)

plates %>%
  filter(is.na(evaporation_beads)) %>% 
  ggplot(aes(control, qubit, col = control)) +
  geom_boxplot(lwd = 1) +
  theme_bw(base_size = 24) +
  theme(legend.position = "none") +
  labs(x = "Control", y = "Qubit (ng/ul)")
ggsave(filename = "./plots/fig_s1b.pdf", width = 7, height = 6)

# beads vs. columns (qubit)
plates %>%
  filter(is.na(evaporation_beads)) %>% 
  filter(control == "no") %>% 
  with(., t.test(qubit ~ method))

# beads vs. columns (nanodrop)
plates %>%
  filter(is.na(evaporation_beads)) %>% 
  filter(control == "no") %>% 
  with(., t.test(nanodrop ~ method))
# p = 0.02342 for qubit
# p = 0.008597 for nanodrop (NEEDS TO BE UPDATED)

# difference in A260 & A280 ratio
plates %>% 
  filter(is.na(evaporation_beads)) %>% 
  filter(control == "no") %>% 
  filter(qubit > 2) %>% 
  mutate(method = str_to_title(method)) %>% 
  ggplot(aes(method, a260_a280, col = method)) +
  geom_boxplot(lwd = 1) +
  geom_hline(yintercept = 1.8, lwd = 1, lty = 2, col = "grey75") +
  theme_bw(base_size = 24) +
  theme(legend.position = "none") +
  labs(x = NULL, y = "Purity (A260/A280)") +
  scale_color_manual(values = colorscheme)
ggsave(filename = "./plots/fig_3d.pdf", width = 5.2, height = 5.4)

# mean a260/a280 ratio
plates %>% 
  filter(is.na(evaporation_beads)) %>% 
  filter(control == "no") %>% 
  group_by(method) %>% 
  summarise(mean_ratio = mean(a260_a280))

# qubit vs. nanodrop 
plates %>% 
  filter(is.na(evaporation_beads)) %>% 
  filter(control == "no") %>% 
  mutate(method = str_to_title(method)) %>% 
  ggplot(aes(qubit, nanodrop, col = method)) +
  geom_point(size = 2, alpha = 1) +
  theme_bw(base_size = 24) +
  theme(legend.position = "none") +
  labs(x = "Qubit (ng/ul)", y = "Nanodrop (ng/ul)") +
  scale_color_manual(values = colorscheme) +
  geom_abline(slope = 1, intercept = 0, lty = 2, col = "grey75")
ggsave(filename = "./plots/fig_3e.pdf", width = 5.4, height = 5.4)

plates %>% 
  filter(is.na(evaporation_beads)) %>% 
  filter(control == "no") %>% 
  with(., cor(qubit, nanodrop))

# nanodrop overestimation
plates %>% 
  filter(is.na(evaporation_beads)) %>%
  mutate(overestimation = nanodrop/qubit) %>% 
  ggplot(aes(qubit, overestimation)) +
  geom_point() +
  theme_bw(base_size = 24) +
  theme(legend.position = "none") +
  labs(x = "Qubit (ng/ul)", y = "Overstimation DNA content")
ggsave(filename = "./plots/fig_s1c.pdf", width = 5.4, height = 5.4)


# association sample age / overestimation
plates %>% 
  filter(sample != "empty") %>% 
  filter(is.na(evaporation_beads)) %>%
  mutate(year = str_extract(sample, pattern = "^(e|ek|a|ak)[0-9]{2}")) %>% 
  mutate(year = str_extract(year, pattern = "[0-9]{2}")) %>% 
  mutate(year = 2000 + as.numeric(year)) %>% 
  ggplot(aes(year, overestimation)) +
  geom_point() + 
  geom_smooth() +
  geom_hline(yintercept = 1, col = "darkgrey", lty = 2) +
  theme_bw(base_size = 20) +
  labs(x = "Sample age (Year)", y = "Overestimation (Nanodrop / Qubit)")
ggsave(filename = "./plots/fig_s1d.pdf", width = 7, height = 6)

plates %>% 
  filter(sample != "empty") %>% 
  filter(is.na(evaporation_beads)) %>%
  mutate(year = str_extract(sample, pattern = "^(e|ek|a|ak)[0-9]{2}")) %>% 
  mutate(year = str_extract(year, pattern = "[0-9]{2}")) %>% 
  mutate(year = 2000 + as.numeric(year)) %>% 
  ggplot(aes(method, overestimation)) +
  geom_boxplot(fill = "darkgrey") +
  theme_bw(base_size = 24) +
  labs(x = "Method", y = "Overestimation (Nanodrop / Qubit)")
ggsave(filename = "./plots/fig_s1e.pdf", width = 7, height = 6)


# p = 0.4169
plates %>% 
  filter(sample != "empty") %>% 
  filter(is.na(evaporation_beads)) %>%
  with(., t.test(overestimation ~ method))

# correction model for Nanodrop overestimation ---------------------------------

# determine best correction model

set.seed(23452341)
n_reps <- 200

for (i in 1:n_reps) {
  print(i)
  
  input <- split_data()
  
  oe <- input$train$nanodrop / input$train$qubit
  oe_mean <- mean(oe)
  oe_median <- median(oe)
  oe_gmean <- exp(mean(log(oe)))
  
  # fit linear model, predict test samples
  model_lm <- lm(qubit ~ nanodrop + a260_a280, data = input$train)
  fit_lm <- predict(model_lm, newdata = input$test)
  
  performance <- input$test %>% 
    dplyr::mutate(fit_lm = fit_lm, 
           fit_mean = nanodrop/oe_mean, 
           fit_median = nanodrop/oe_median, 
           fit_gmean = nanodrop/oe_gmean) %>% 
    pivot_longer(cols = starts_with("fit_")) %>% 
    dplyr::mutate(deviation_squared = (qubit - value)^2) %>% 
    group_by(name) %>% 
    dplyr::summarise(mse = mean(deviation_squared)) %>% 
    mutate(rep = i)
  
  if (i == 1) {
    model_comp <- performance
  } else {
    model_comp <- rbind(model_comp, performance)
  }
}

rm(input, oe_gmean, oe_mean, oe_median, oe, model_lm, fit_lm, performance, i)

model_comp

model_comp_summarized <- model_comp %>% 
  group_by(rep) %>% 
  dplyr::mutate(rank = rank(mse)) %>% 
  ungroup %>% 
  group_by(name, rank) %>%
  dplyr::summarise(n = n()/n_reps*100) %>%
  ungroup() %>% 
  dplyr::mutate(name = str_replace(name, pattern = "fit_", replacement = ""))

model_comp_summarized %>% 
  mutate(rank = paste0("Rank ", rank)) %>% 
  ggplot(aes(name, n, fill = name)) +
    geom_col() +
  facet_wrap(facets = vars(rank), nrow = 4) +
  labs(x = NULL, y = "Proportion (%)") +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
       legend.position = "none", legend.title = element_blank()) +
  scale_fill_manual(values = colorscheme)
ggsave(filename = "./plots/fig_3f.pdf", width = 2.7, height = 7.9)


# how many samples are needed to asses geometric mean?

bootstraps = 200
sample_size = c(5:15, seq(from = 20, to = 50, by = 5))
bootstrap_results <- matrix(data = NaN, nrow = bootstraps*length(sample_size), ncol = 5)

for (i in 1:bootstraps){
  bootstrapped_data <- split_data(p_train = 0.75)$train
  print(i)
  
  gm_real <- geom_mean(bootstrapped_data$nanodrop / bootstrapped_data$qubit)
  
  for (j in 1:length(sample_size)){
    indeces_rand <- sample(1:nrow(bootstrapped_data), size = sample_size[j], replace = FALSE)
    indeces_rank <- order(bootstrapped_data$nanodrop)[ceiling(seq(1, to = nrow(bootstrapped_data), 
                                                                  length.out = sample_size[j]))]
    oe_rand <- with(bootstrapped_data, geom_mean(nanodrop[indeces_rand] / qubit[indeces_rand]))
    oe_rank <- with(bootstrapped_data, geom_mean(nanodrop[indeces_rank] / qubit[indeces_rank]))
    bootstrap_results[(i-1)*length(sample_size) + j, ] <- c(i, sample_size[j], oe_rand, oe_rank, gm_real)
  }
}

colnames(bootstrap_results) <- c("bootstrap", "sample_size", "random", "ranked", "gm_real")
bootstrap_results <- as_tibble(bootstrap_results)

bootstrap_results %>% 
  pivot_longer(cols = c(random, ranked)) %>% 
  mutate(value_norm = value/gm_real) %>% 
  ggplot(aes(sample_size, value_norm)) +
  geom_point() +
  facet_wrap(facets = vars(name))
ggsave(filename = "./plots/fig_s2.pdf", width = 8.3, height = 7.9)

bootstrap_results %>% 
  pivot_longer(cols = c(random, ranked)) %>% 
  group_by(sample_size, name) %>% 
  dplyr::summarise(var = var(value)) %>% 
  ggplot(aes(sample_size, var, col = name)) +
  geom_smooth(lwd = 1.5, se = FALSE, formula = y ~ log(x), col = "steelblue") +
  geom_point(size = 3, col = "darkgrey") +
  theme_bw(base_size = 24) +
  facet_wrap(facets = vars(name)) +
  theme(legend.position = "none") +
  labs(x = "Sample Size", y = "Spread (variance)")
ggsave(filename = "./plots/fig_3g.pdf", width = 6.8, height = 5.5)


# calculate correction factor for actual plates

calculate_correction <- function(nd, qb, n){
  ind <- order(nd)[ceiling(seq(1, to = length(nd), length.out = n))]
  cf <- geom_mean(nd[ind]/qb[ind])
  return(cf)
}

cf <- plates %>% 
  filter(control == "no" & is.na(evaporation_beads) & replicate == "1") %>% 
  group_by(method) %>% 
  mutate(cf = calculate_correction(nanodrop, qubit, n = 20)) %>% 
  mutate(pred = nanodrop / cf) %>% 
  mutate(err_rel = pred / qubit) %>% 
  mutate(within_bounds = ifelse(err_rel > 2 | err_rel < 0.5, 0, 1))

cf %>% 
  dplyr::summarise(wb = sum(within_bounds)/n())

cf %>% 
  pivot_longer(cols = c(overestimation, err_rel)) %>% 
  ggplot(aes(name, value)) + 
  geom_boxplot()

cf %>% 
  pivot_longer(cols = c(overestimation, err_rel)) %>% 
  group_by(name) %>% 
  dplyr::summarise(median = median(value))


# geom_boxplot())


## DNA methylation data --------------------------------------------------------

array_anno <- readxl::read_xlsx(path = "./input/array_annotation.xlsx")

# different EPIC array version were used, load separately
idats_v1 <- array_anno %>% 
  filter(array == "epic_v1") %>% 
  pull(basename) %>% 
  paste0("./input/idat/lmu/", .)

idats_v2 <- array_anno %>% 
  filter(array == "epic_v2") %>% 
  pull(basename) %>% 
  paste0("./input/idat/lmu/", .)

raw_v1 <- minfi::read.metharray(basenames = idats_v1, force = TRUE)
raw_v2 <- minfi::read.metharray(basenames = idats_v2, force = TRUE)

# change v1 annotation to hg38 for compatibility with v2
raw_v1@annotation = c(array = "IlluminaHumanMethylationEPIC", annotation = "ilm10b5.hg38")

# normalize data
processed_v1 <- preprocessNoob(rgSet = raw_v1)

processed_v2 <- preprocessNoob(rgSet = raw_v2)
processed_v2 <- remove_epic_suffix(processed_v2)

# extract beta values
betas_v1 <- getBeta(processed_v1)
betas_v2 <- getBeta(processed_v2)

# shrink data to probes available on both platforms
common_probes <- intersect(rownames(betas_v1), rownames(betas_v2))
betas_v1 <- betas_v1[common_probes, ]
betas_v2 <- betas_v2[common_probes, ]

betas <- cbind(betas_v1, betas_v2)
betas <- betas[, array_anno$basename]
colnames(betas) <- array_anno$id
rm(betas_v1, betas_v2)


# read control data
controls_v1 <- list.files(path = "./input/idat/controls_v1/", pattern = "Grn.idat", full.names = TRUE)
controls_v1 <- minfi::read.metharray(basenames = controls_v1) %>% 
  preprocessNoob

controls_v2 <- list.files(path = "./input/idat/controls_v2/", pattern = "Grn.idat", full.names = TRUE)
controls_v2 <- minfi::read.metharray(basenames = controls_v2) %>% 
  preprocessNoob
controls_v2 <- remove_epic_suffix(controls_v2)

# failed pobes per purification scheme (based on detection p-value > 0.05)
detp_v1 <- detectionP(raw_v1)
detp_v1 <- apply(detp_v1, 2, function(x) sum(x >= 0.05)/length(detp_v1))*100

detp_v2 <- detectionP(raw_v2)
detp_v2 <- apply(detp_v2, 2, function(x) sum(x >= 0.05)/length(detp_v2))*100

detp <- c(detp_v1, detp_v2)
detp <- detp[array_anno$basename]

array_anno <- array_anno %>% 
  mutate(det_p = detp)

array_anno <- array_anno %>% 
  mutate(purification = ifelse(purification == "gold_standard", "reference", purification)) %>% 
  mutate(purification = str_to_title(purification))

array_anno %>% 
  mutate(id = str_extract(id, pattern = "s[0-9]{1}")) %>% 
  arrange(id) %>% 
  ggplot(aes(x = id, y = det_p, col = purification)) +
  geom_jitter(size = 4, width = 0.3, alpha = 1) +
  scale_color_manual(values = colorscheme) +
  theme_bw(base_size = 20) +
  labs(x = "Sample", y = "Failed probes (%)") +
  theme(legend.position = "top", legend.justification.top = 1, legend.title = element_blank())
ggsave(filename = "./plots/fig_4a.pdf", width = 4.2, height = 6.3)


# density plot / beta value distributions
betas %>% 
  as_tibble() %>% 
  pivot_longer(cols = everything(), names_to = "sample", values_to = "beta") %>% 
  mutate(purification = str_extract(sample, pattern = "[^_]*$"), 
         sample = str_extract(sample, pattern = "[^_]*")) %>% 
  mutate(purification = ifelse(purification == "col", "Column", ifelse(purification == "std", "Reference", "Beads"))) %>% 
  mutate(purification = str_to_title(purification)) %>% 
  ggplot(aes(sample, beta, fill = purification)) +
  geom_violin(alpha = .4) +
  theme_bw(base_size = 18) +
  scale_fill_manual(values = colorscheme) +
  theme(legend.position = "top", legend.title = element_blank()) +
  labs(x = "Sample", y = "Methylation (beta)")
ggsave(filename = "./plots/fig_4b.pdf", width = 7.5, height = 5.1)

# example plot of beta values beads vs gold standard
betas %>% 
  as_tibble %>% 
  slice_sample(n = 100000) %>% 
  ggplot(aes(x = s1_std, y = s1_beads)) +
  geom_point(alpha = 0.05, col = "darkgrey") +
  geom_abline(slope = 1, intercept = 0, colour = "black", lty = 2) +
  theme_bw(base_size = 20) +
  labs(x = "Beta (reference)", y ="Beta (beads)")
ggsave(filename = "./plots/fig_4d.tiff", width = 6.5, height = 6.5)


# heatmap of beta values
betas_var <- apply(betas, 1, var)
betas_var_top <- order(betas_var, decreasing = TRUE)[1:5000]

betas_subset <- betas[betas_var_top, ]
colnames(betas_subset) <- str_replace(colnames(betas_subset), "std", "ref")

# include dendrograms
superheat::superheat(X = betas_subset,
                     col.dendrogram = TRUE, 
                     row.dendrogram = TRUE, 
                     pretty.order.rows = TRUE, 
                     pretty.order.cols = TRUE, 
                     heat.pal = rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")),
                     left.label.col = "white", 
                     bottom.label.col = "white", 
                     bottom.label.text.angle = 90)

# all sample correlations
betas_cor <- cor(betas, method = "spearman")
colnames(betas_cor) <- str_replace(colnames(betas_cor), "std", "ref")
rownames(betas_cor) <- str_replace(rownames(betas_cor), "std", "ref")
  
superheat::superheat(X = as.matrix(betas_cor),
                     X.text = round(betas_cor, 2),
                     col.dendrogram = TRUE, 
                     row.dendrogram = TRUE, 
                     pretty.order.rows = TRUE, 
                     pretty.order.cols = TRUE, 
                     heat.lim = c(0.85, 1), 
                     heat.pal = rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")),
                     left.label.col = "white", 
                     bottom.label.col = "white", 
                     bottom.label.text.angle = 90)

# CNV plots

conumee_anno_v1 <- CNV.create_anno(array_type = "EPIC", chrXY = FALSE)
conumee_anno_v2 <- CNV.create_anno(array_type = "EPICv2", chrXY = FALSE)

conumee_data_v1 <- CNV.load(input = processed_v1)
conumee_data_v2 <- CNV.load(input = processed_v2)

conumee_ctrl_v1 <- CNV.load(input = controls_v1)
conumee_ctrl_v2 <- CNV.load(input = controls_v2)

conumee_v1 <- CNV.fit(query = conumee_data_v1, anno = conumee_anno_v1, ref = conumee_ctrl_v1)
conumee_v1 <- CNV.bin(conumee_v1)
conumee_v1 <- CNV.segment(conumee_v1)

conumee_v2 <- CNV.fit(query = conumee_data_v2, anno = conumee_anno_v2, ref = conumee_ctrl_v2)
conumee_v2 <- CNV.bin(conumee_v2)
conumee_v2 <- CNV.segment(conumee_v2)

dir_cnv_plots <- "./cnv_plots/"
CNV.genomeplot(object = conumee_v1, directory = dir_cnv_plots, output = "png", width = 12, height = 4, 
               cols = rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))
CNV.genomeplot(object = conumee_v2, directory =  dir_cnv_plots, output = "png", width = 12, height = 4, 
               cols = rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))


# rename the resulting files
cnv_files <- list.files(path = dir_cnv_plots, pattern = ".png", full.names = TRUE)

for (i in 1:nrow(array_anno)){
  cb <- array_anno$basename[i]
  crep <- array_anno$id[i]
  f_in <- grep(cb, cnv_files, value = TRUE)
  f_out <- str_replace(string = f_in, pattern = cb, crep)
  file.move(from = f_in, to = f_out)
}



# NGS: Exome sequencing --------------------------------------------------------


## Basics stats (reads, base quality, duplication rate) -------------------------

# read raw sequencing data quality stats

base_qualities <- read_tsv("./input/ngs/quality_scores/fastq_quality_scores.tsv", 
                        col_names = c("id", "reads_n", "sampled_n", "q20", "q30"))

base_qualities <- base_qualities %>% 
  mutate(read = ifelse(str_detect(id, "R1"), "R1", "R2"), 
         method = str_extract(id, pattern = "(Beads|Columns|Promega)")) %>% 
  mutate(method = ifelse(method == "Promega", "Reference", method))  %>% 
  mutate(sample = str_extract(id, pattern = "E[0-9]{2}_[0-9]*")) %>% 
  mutate(q20 = q20 * 100, 
         q30 = q30 * 100) %>% 
  dplyr::select(id, sample, read, method, reads_n, sampled_n, q20, q30)

base_qualities <- base_qualities %>% 
  group_by(sample, method) %>% 
  dplyr::summarise(q20 = round(mean(q20), 2),
                   q30 = round(mean(q30), 2), 
                   reads_n = min(reads_n), 
                   sampled_n = min(sampled_n)) %>% 
  ungroup


# deduplication rate

deduplication_rate <- read_csv(file = "./input/ngs/deduplication_stats/summarised.csv")

deduplication_rate <- deduplication_rate %>% 
  mutate(sample = str_extract(library, pattern = "E[0-9]{2}_[0-9]*"), 
         method = str_extract(library, pattern = "(Beads|Columns|Promega)")) %>% 
  mutate(method = ifelse(method == "Promega", "Reference", method)) %>% 
  mutate(percent_duplication = percent_duplication*100) %>% 
  drop_na() %>% 
  dplyr::select(sample, method, percent_duplication)


# samtools flagstat summaries

flagstat_files <- list.files(path = "./input/ngs/flagstat/", 
                             pattern = "flagstat.tsv", full.names = TRUE)

flagstat <- flagstat_files %>% 
  as.list %>% 
  map(.f = read_flagstat) %>% 
  Reduce(f = bind_rows, x = .)

flagstat <- flagstat %>% 
  pivot_wider(names_from = measure, values_from = data)
rm(flagstat_files)

# combined into one dataframe
ngs_stats <- Reduce(x = list(flagstat, deduplication_rate, base_qualities), 
                    f = function(x, y) left_join(x, y, by = c("sample", "method")))
rm(flagstat, base_qualities, deduplication_rate)



## Read insert size distribution ------------------------------------------------

bam_files <- list.files(path = "./input/ngs/subsampled_bams/", full.names = TRUE, 
                        pattern = "bam")

# Extract insert sizes for each BAM file and store in a list
insert_size_data <- lapply(bam_files, function(file) {
  sizes <- extract_insert_sizes(file)
  data.frame(file = basename(file), insert_size = sizes)
})

# Combine results into a single data frame
insert_size_df <- do.call(rbind, insert_size_data)
insert_size_df <- as_tibble(insert_size_df)
rm(insert_size_data)

insert_size_df <- insert_size_df %>% 
  mutate(method = str_extract(file, pattern = "(Beads|Columns|Promega|Control)"), 
         sample = str_extract(file, pattern = "E[0-9]{2}_[0-9]+"))

insert_size_df <- insert_size_df %>% 
  mutate(method = ifelse(method == "Promega", "Reference", method))

insert_size_df <- insert_size_df %>% 
  filter(insert_size != 0) # remove unpaired reads

insert_size_df %>% 
  mutate(insert_size = abs(insert_size)) %>% 
  group_by(sample, method) %>% 
  dplyr::summarise(q10 = quantile(insert_size, probs = 0.1), 
            q50 = median(insert_size),
            q90 = quantile(insert_size, probs = 0.9))

insert_size_df %>% 
  slice_sample(by = method, n = 20000) %>% 
  mutate(insert_size = abs(insert_size)) %>% 
  filter(insert_size < 500) %>% 
  filter(insert_size != 0) %>% 
  ggplot(aes(insert_size, fill = method)) +
  geom_histogram() +
  facet_grid(cols = vars(sample), rows = vars(fct_relevel(method,'Beads','Columns','Reference','Control'))) +
  theme_grey(base_size = 16) +
  theme(legend.position = "none") +
  scale_fill_manual(values = c(colorscheme[1:2], "grey50", colorscheme[3])) +
  labs(x = "Insert size (bp)", y = "Count")
ggsave(filename = "./plots/fig_5a.pdf", width = 10, height = 8)



# NGS: on-target rate and coverage uniformity ----------------------------------

bed_path <- "./input/ngs/coverage/"
bed_data <- list.files(path = bed_path, full.names = TRUE)
bed_data <- bed_data %>% 
  as.list %>% 
  map(.f = read_bed)

bed_df <- Reduce(f = bind_rows, x = bed_data)

bed_df <- bed_df %>% 
  mutate(method = str_extract(id, pattern = "(Beads|Columns|Promega|Control)"), 
         sample = str_extract(id, pattern = "E[0-9]{2}_[0-9]+")) %>% 
  mutate(method = ifelse(method == "Promega", "Reference", method))

# calculate between-sample coverage
coverage_cor <- bed_df %>% 
  select(id, chr, start, end, coverage) %>% 
  pivot_wider(id_cols = c(chr, start, end), names_from = id, values_from = coverage) %>% 
  select(4:11) %>% 
  as.matrix() %>% 
  cor(use = "pairwise")

mean(coverage_cor[upper.tri(coverage_cor)]) # 0.9703945
coverage_cor[upper.tri(coverage_cor)] %>% 
  as.data.frame() %>% 
  ggplot(aes(.)) +
  geom_boxplot(fill = "steelblue") + 
  coord_flip() +
  theme_bw(base_size = 20) +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())

# coverage distribution
bed_df %>% 
  filter(id == "E14_53319_Beads") %>% 
  slice_sample(n = 50000) %>% 
  ggplot(aes(coverage)) +
  geom_histogram(fill = "steelblue", alpha = 0.7, bins = 100) +
  scale_x_log10() +
  labs(x = "Coverage", y = "Count") +
  theme_bw(base_size = 28)
ggsave(filename = "./plots/fig_5b.pdf", width = 8, height = 7)

# relationship between feature size and coverage
bed_df %>% 
  filter(id == "E14_53319_Beads") %>% 
  slice_sample(n = 50000) %>% 
  ggplot(aes(length, coverage)) +
  #geom_bin_2d(bins = 80) +
  geom_point(col = "steelblue", alpha = 0.2) +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "Feature size (bp)", y = "Coverage (no. of reads)") +
  theme_bw(base_size = 28)
ggsave(filename = "./plots/fig_5c.pdf", width = 8, height = 7)

# plot Lorentz curve
bed_df %>% 
  group_by(sample, method) %>% 
  filter(method != "Control") %>% 
  slice_sample(n = 1000) %>% 
  arrange(coverage) %>% 
  mutate(cumulative_coverage = cumsum(coverage)) %>% 
  mutate(rank = rank(coverage)) %>%
  mutate(cumululative_coverage_norm = cumulative_coverage/max(cumulative_coverage), 
         rank_norm = rank/max(rank)) %>% 
  ungroup %>% 
  ggplot(aes(rank_norm, cumululative_coverage_norm, col = method)) +
  geom_line(lwd = 1) +
  geom_abline(slope = 1, intercept = 0, col = "grey75", lty = 2) +
  facet_grid(cols = vars(sample), rows = vars(method)) +
  labs(x = "Target regions (ranked by coverage)", y = "Cumulative Coverage") +
  scale_colour_manual(values = colorscheme) +
  theme_bw(base_size = 16) +
  theme(legend.position = "none")

# overall on-target rate
on_target <- bed_df %>% 
  group_by(sample, method) %>% 
  dplyr::summarise(n_target = sum(coverage))

ngs_stats <- left_join(ngs_stats, on_target) %>% 
  mutate(perc_target = n_target/mapped*100)

# calculate uniformity measures 
uniformity <- bed_df %>% 
  group_by(id) %>% 
  dplyr::summarise(cv = mean(coverage)/sd(coverage), 
                   gini = Gini(coverage),
                   p90_p10 = quantile(coverage, probs = 0.9)/quantile(coverage, probs = 0.1),  
                   p75_p25 = quantile(coverage, probs = 0.75)/quantile(coverage, probs = 0.25))

bed_df %>% 
  group_by(id) %>% 
  dplyr::summarise(p90 = quantile(coverage, probs = 0.9), 
                   p10 = quantile(coverage, probs = 0.1)) %>% 
  ungroup %>% 
  filter(!str_detect(id, pattern = "Control")) %>% 
  summarise(mean_p90 = mean(p90), 
            mean_p10 = mean(p10))


uniformity <- uniformity %>% 
  mutate(sample = str_extract(id, pattern = "E[0-9]{2}_[0-9]+"), 
         method = str_extract(id, "(Beads|Columns|Promega|Control)")) %>% 
  mutate(method = ifelse(method == "Promega", "Reference", method))
mean(uniformity$p90_p10) # 5.199427

uniformity %>% pivot_longer(cols = c(cv, gini, p90_p10, p75_p25)) %>% 
  filter(method != "Control") %>% 
  filter(name %in% c("cv", "gini")) %>% 
  ggplot(aes(factor(sample), value, fill = method)) +
  geom_col(position = "dodge") +
  facet_wrap(nrow = 4, facets = vars(name), scales = "free") +
  theme_bw(base_size = 18) +
  labs(x = NULL, y = NULL) +
  scale_fill_manual(values = colorscheme) +
  theme(legend.title = element_blank(), 
        legend.position = "top")
ggsave(filename = "./plots/fig_5d.pdf", width = 5, height = 8)




# NGS: Variant calling ---------------------------------------------------------

vcf_files = list.files(path = "./input/ngs/vcf_files/", 
                       pattern = "snvs.vcf.gz$", 
                       full.names = TRUE, 
                       recursive = TRUE)
vcf_id <- str_extract(string = vcf_files, pattern = "E[0-9]{2}_[0-9]{5}_(Beads|Columns|Promega)")
vcf_id <- str_replace(vcf_id, pattern = "Promega", replacement = "Reference")
vcf_samples <- str_extract(string = vcf_id, pattern = "E[0-9]{2}_[0-9]{5}")
vcf_method <- str_extract(string = vcf_id, pattern = "(Beads|Columns|Reference)")

# read data
vcf_data <- vcf_files %>% 
  as.list %>% 
  map(.f = read_vcf_file)

# add pseudo id for each variant
vcf_data <- vcf_data %>% 
  map(.f = function(x){
    x <- x %>% 
      mutate(var_id = paste(chr, pos, ref, alt, sep = "_"))
  })

vcf_data_pass <- vcf_data %>% 
  map(.f = function(x){
    x <- x %>% 
      filter(filter == "PASS")
  })
  

# calculate overlap between variants
min_depth = 100
mind_depth_probs <- seq(0.0, 0.9, by = 0.1)
min_depth_quantiles <- quantile(vcf_data[[i]]$depth, probs = mind_depth_probs)

variant_accordance <- as.list(min_depth_quantiles) %>% 
  map(calculate_variants_overlap, vcf_list = vcf_data, sample_names = vcf_id)

variant_accordance_thresholds <- variant_accordance %>% 
  map(.f = function(x) return(c(mean(x[1,2], x[1,3], x[2,3]), mean(x[5,4], x[6,4], x[6,5])))) %>% 
  Reduce(f = rbind)

colnames(variant_accordance_thresholds) <- c("E14_53319", "E15_46004")
variant_accordance_thresholds <- variant_accordance_thresholds %>% 
  as_tibble() %>% 
  mutate(cutoff_q = mind_depth_probs, 
         cutoff_n = min_depth_quantiles)

variant_accordance_thresholds %>% 
  pivot_longer(cols = 1:2) %>% 
  ggplot(aes(cutoff_n, value, col = name)) +
  geom_line(lwd = 2) +
  theme_bw(base_size = 24) +
  labs(x = "Cutoff (n reads)", y = "Average Jaccard index") +
  theme(legend.position = "top", legend.title = element_blank())

# calculate variance accordance at read depth 100
variance_accordance_100 <- calculate_variants_overlap(vcf_list = vcf_data, 
                                                      min_depth = 100, 
                                                      sample_names = vcf_id)

superheat::superheat(X = variance_accordance_100,
                     X.text = round(variance_accordance_100, 2),
                     col.dendrogram = TRUE, 
                     row.dendrogram = TRUE, 
                     pretty.order.rows = TRUE, 
                     pretty.order.cols = TRUE, 
                     heat.pal = rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")),
                     left.label.col = "white", 
                     bottom.label.col = "white", 
                     bottom.label.text.angle = 90)




## NGS: mutational signatures --------------------------------------------------

sn <- str_extract(string = vcf_files, pattern = "E[0-9]{2}_[0-9]{5}_(Beads|Columns|Promega)")
sn <- str_replace(sn, pattern = "Promega", replacement = "Reference")

grl <- readRDS(file = "./input/vcf_granges.rds")
#grl <- read_vcfs_as_granges(vcf_files, sample_names = sn , ref_genome)
#saveRDS(object = grl, file = "./input/vcf_granges.rds")

# extract only mutations which passed Strelka filters
grl2 <- grl %>% 
  map(.f = function(x) return(x[x$FILTER == "PASS", ]))

# plot mutation occurences between extraction methods
type_occurrences <- mut_type_occurrences(grl2, ref_genome)

plot_spectrum(type_occurrences, CT = TRUE, by = vcf_method)
ggsave(filename = './plots/fig_5f.pdf', width = 10, height = 4)


# Chi-Square test for mutation occurences
type_occurrences %>% 
  as_tibble(rownames = "sample") %>% 
  add_column(method = str_extract(.$sample, pattern = "(Beads|Columns|Reference)"), .before = 1) %>% 
  select(-"sample") %>% 
  pivot_longer(cols = 2:9) %>% 
  group_by(method, name) %>% 
  summarise(value = mean(value)) %>% 
  ungroup %>% 
  pivot_wider(names_from = name, values_from = value) %>% 
  select(-method) %>% 
  as.matrix() %>% 
  chisq.test() 
# p-value = 0.9964


