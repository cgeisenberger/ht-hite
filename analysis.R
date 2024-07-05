library(tidyverse)
library(minfi)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
library(conumee2.0)
library(random)

source("./functions.R")


# color scales
cols_fig1 <- RColorBrewer::brewer.pal(9, "Oranges")[4:9]
cols_fig2 <- RColorBrewer::brewer.pal(9, "Greens")[4:9]
cols_fig3 <- RColorBrewer::brewer.pal(9, "Purples")



# validation of manual HiTE ----------------------------------------------------

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
  scale_color_manual(values = cols_fig1)
#ggsave(filename = "./plots/fig_1a.pdf", width = 7, height = 5.5)


manual %>% 
  ggplot(aes(method, a260_a280, col = as.factor(method))) +
  geom_boxplot(lwd = 1) +
  theme_bw(base_size = 20) +
  geom_hline(yintercept = 1.8, col = "grey", lty = 2) +
  theme(legend.position = "none", legend.title = element_blank(), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  facet_wrap(facets = vars(incubation)) +
  labs(x = NULL, y = "DNA Purity (A260/A280)") +
  scale_color_manual(values = cols_fig1)
#ggsave(filename = "./plots/fig_1b.pdf", width = 7, height = 5.5)


manual %>% 
  ggplot(aes(qubit, nanodrop,  col = cleanup, shape = cleanup)) +
  geom_point(size = 5) +
  theme_bw(base_size = 20) +
  theme(legend.position = "top", legend.title = element_blank()) +
  geom_abline(slope = 1, intercept = 0, lty = 2, col = "grey") +
  labs(x = "Qubit ng/ul", y = "Nanodrop (ng/ul)") +
  xlim(0, 25) +
  ylim(0, 80) +
  scale_color_manual(values = cols_fig1)
#ggsave(filename = "./plots/fig_1c.pdf", width = 7, height = 5.5)

manual %>% 
  #filter(cleanup != "Promega_1h") %>% 
  with(., cor(qubit, nanodrop))

manual %>% 
  #filter(cleanup != "Promega_1h") %>% 
  with(., cor(qubit, nanodrop, method = "spearman"))

manual %>% 
  mutate(fold_x = nanodrop/qubit) %>% 
  filter(cleanup != "Promega_1h") %>% 
  filter(qubit > 10) %>% 
  summarise(mean_fold = mean(fold_x), 
            median_fold = median(fold_x))

manual %>% 
  mutate(fold_x = nanodrop/qubit) %>% 
  filter(cleanup != "Promega_1h") %>%
  ggplot(aes(x = as.factor(1), y = fold_x)) +
  geom_boxplot(lwd = 2, col = cols_fig1[1], alpha = 0.7) +
  geom_jitter(size = 4, width = 0.2, col = cols_fig1[5]) +
  theme_bw(base_size = 24) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  labs(x = NULL, y = "Overestimation Nanodrop (fold)")
#ggsave(filename = "./plots/fig_1d.pdf", width = 3.8, height = 5.5)

manual %>% 
  mutate(fold_x = nanodrop/qubit) %>% 
  filter(cleanup != "Promega_1h") %>%
  ggplot(aes(qubit, fold_x)) +
  geom_smooth(col = cols_fig1[3], fill = cols_fig1[1]) +
  geom_point(size = 4, col = cols_fig1[5]) +
  theme_bw(base_size = 24) +
  labs(x = "Qubit (ng/ul)", y = "Overestimation Nanodrop (fold)")
#ggsave(filename = "./plots/fig_1e.pdf", width = 4.4, height = 5.5)



# plate-based processing -------------------------------------------------------

plates <- readxl::read_xlsx(path = "./input/20240522_plate_workflow.xlsx")

plates <- plates %>% 
  mutate(overestimation = nanodrop/qubit)

plates %>% 
  filter(is.na(evaporation_beads)) %>% 
  select(-c("well", "sample", "replicate", "evaporation_beads")) %>% 
  GGally::ggpairs()

# overall yield (qubit)
plates %>%
  filter(is.na(evaporation_beads)) %>% 
  filter(control == "no") %>% 
  ggplot(aes(qubit)) + 
  geom_histogram(aes(y = after_stat(density)),
                 alpha = 0.7, bins = 40, fill = cols_fig2[6]) +
  geom_density(lwd = 1, color = cols_fig2[6]) +
  theme_bw(base_size = 20)
ggsave(filename = "./plots/fig_2b.pdf", width = 8.6, height = 2.8)

# overall yield (nanodrop)
plates %>%
  filter(is.na(evaporation_beads)) %>% 
  filter(control == "no") %>% 
  ggplot(aes(nanodrop)) + 
  geom_histogram(aes(y = after_stat(density)),
                 alpha = 0.7, bins = 40, fill = cols_fig2[6]) +
  geom_density(lwd = 1, color = cols_fig2[6]) +
  theme_bw(base_size = 20)
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
  ggplot(aes(method, conc, col = method)) +
  geom_boxplot(lwd = 1) +
  theme_bw(base_size = 24) +
  theme(legend.position = "none") +
  facet_wrap(facets = vars(measurement), scales = "free") +
  labs(x = NULL, y = "DNA ng/ul") +
  scale_color_manual(values = cols_fig2[c(1, 4)])
ggsave(filename = "./plots/fig_2c.pdf", width = 8.6, height = 4.3)


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
  scale_color_manual(values = cols_fig2[c(1, 4)])
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
  ggplot(aes(method, a260_a280, col = method)) +
  geom_boxplot(lwd = 1) +
  geom_hline(yintercept = 1.8, lwd = 1, lty = 2, col = "grey75") +
  theme_bw(base_size = 24) +
  theme(legend.position = "none") +
  labs(x = NULL, y = "A260/A280 ratio") +
  scale_color_manual(values = cols_fig2[c(1, 4)])
ggsave(filename = "./plots/fig_2e.pdf", width = 5.2, height = 5.4)

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
  ggplot(aes(qubit, nanodrop, col = method)) +
  geom_point(size = 2, alpha = 1) +
  theme_bw(base_size = 24) +
  theme(legend.position = "none") +
  labs(x = "Qubit (ng/ul)", y = "Nanodrop (ng/ul)") +
  scale_color_manual(values = cols_fig2[c(1, 4)]) +
  geom_abline(slope = 1, intercept = 0, lty = 2, col = "grey75")
ggsave(filename = "./plots/fig_2d.pdf", width = 5.4, height = 5.4)

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
    mutate(#fit = fit$predicted,
           fit_lm = fit_lm, 
           fit_mean = nanodrop/oe_mean, 
           fit_median = nanodrop/oe_median, 
           fit_gmean = nanodrop/oe_gmean) %>% 
    pivot_longer(cols = starts_with("fit")) %>% 
    mutate(se = (qubit - value)^2) %>% 
    group_by(name) %>% 
    summarise(mse = mean(se)) %>% 
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
  mutate(rank = rank(mse)) %>% 
  ungroup %>% 
  group_by(name, rank) %>%
  summarise(n = n()/n_reps*100) %>%
  ungroup() %>% 
  mutate(name = str_replace(name, pattern = "fit_", replacement = ""))

model_comp_summarized %>% 
  ggplot(aes(name, n, fill = name)) +
    geom_col() +
  facet_wrap(facets = vars(rank), nrow = 4) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
       legend.position = "none", legend.title = element_blank()) +
  scale_fill_manual(values = cols_fig2)
ggsave(filename = "./plots/fig_s2f.pdf", width = 2.7, height = 7.9)


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
  summarise(var = var(value)) %>% 
  ggplot(aes(sample_size, var, col = name)) +
  geom_smooth(lwd = 1.5, se = FALSE, formula = y ~ log(x), col = cols_fig2[1]) +
  geom_point(size = 3, col = cols_fig2[4]) +
  theme_bw(base_size = 24) +
  facet_wrap(facets = vars(name)) +
  theme(legend.position = "none") +
  labs(x = NULL, y = "Spread (variance)")
ggsave(filename = "./plots/fig_2g.pdf", width = 6.8, height = 5.5)


# calculate correction factor for actual plates

calculate_correction <- function(nd, qb, n){
  ind <- order(nd)[ceiling(seq(1, to = length(nd), length.out = n))]
  cf <- geom_mean(nd[ind]/qb[ind])
  return(cf)
}

cf <- plates_filtered %>% 
  group_by(method) %>% 
  mutate(cf = calculate_correction(nanodrop, qubit, n = 20)) %>% 
  mutate(pred = nanodrop / cf) %>% 
  mutate(err_rel = pred / qubit) %>% 
  mutate(within_bounds = ifelse(err_rel > 2 | err_rel < 0.5, 0, 1))

cf %>% 
  summarise(wb = sum(within_bounds)/n())



# comparison of methylation data -----------------------------------------------

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

array_anno %>% 
  arrange(purification) %>%
  mutate(purification = ifelse(purification == "gold_standard", "std", purification)) %>% 
  ggplot(aes(x = purification, y = det_p, col = purification)) +
  geom_jitter(size = 4, width = 0.3, alpha = 0.7) +
  scale_color_manual(values = cols_fig3[c(5,7,9)]) +
  theme_bw(base_size = 20) +
  labs(x = NULL, y = "Failed probes (%)") +
  theme(legend.position = "top", legend.title = element_blank()) +
  scale_x_discrete(labels = str_replace(arrange(array_anno, purification)$id, pattern = "_.*", replacement = ""))
ggsave(filename = "./plots/fig_3a.pdf", width = 4.2, height = 6.3)


# density plot / beta value distributions
betas %>% 
  as_tibble() %>% 
  pivot_longer(cols = everything(), names_to = "sample", values_to = "beta") %>% 
  mutate(purification = str_extract(sample, pattern = "[^_]*$"), 
         sample = str_extract(sample, pattern = "[^_]*")) %>% 
  ggplot(aes(sample, beta, fill = purification)) +
  geom_violin(alpha = .4) +
  theme_bw(base_size = 20) +
  scale_fill_manual(values = cols_fig3[c(3,6,9)]) +
  theme(legend.position = "top")
ggsave(filename = "./plots/fig_3b.pdf", width = 7.5, height = 5.1)

# example plot of beta values column bs gold standard
betas %>% 
  as_tibble %>% 
  slice_sample(n = 100000) %>% 
  ggplot(aes(x = s1_std, y = s1_col)) +
  geom_point(alpha = 0.05, col = "#542788") +
  geom_abline(slope = 1, intercept = 0, colour = "grey75", lty = 2) +
  theme_bw(base_size = 20) +
  labs(x = "Beta (gold standard)", y ="Beta (column)")
ggsave(filename = "./plots/fig_3d.pdf", width = 6.5, height = 6.5)


# heatmap of beta values
betas_var <- apply(betas, 1, var)
betas_var_top <- order(betas_var, decreasing = TRUE)[1:5000]

superheat::superheat(X = as.matrix(betas[betas_var_top, ]),
                     pretty.order.rows = TRUE, 
                     pretty.order.cols = TRUE, 
                     heat.pal = cols_fig3, 
                     left.label.col = "white", 
                     bottom.label.col = "white", 
                     bottom.label.text.angle = 90)

superheat::superheat(X = as.matrix(betas[betas_var_top, ]),col.dendrogram = TRUE, 
                     pretty.order.rows = TRUE, 
                     pretty.order.cols = TRUE, 
                     heat.pal = cols_fig3, 
                     left.label.col = "white", 
                     bottom.label.col = "white", 
                     bottom.label.text.angle = 90)

# all sample correlations
betas_cor <- cor(betas)
betas_cor

superheat::superheat(X = as.matrix(betas_cor),
                     pretty.order.rows = TRUE, 
                     pretty.order.cols = TRUE, 
                     heat.lim = c(0.8, 1), 
                     heat.pal = cols_fig3, 
                     left.label.col = "white", 
                     bottom.label.col = "white", 
                     bottom.label.text.angle = 90)

superheat::superheat(X = as.matrix(betas_cor),
                     pretty.order.rows = TRUE, 
                     pretty.order.cols = TRUE, 
                     heat.lim = c(0.8, 1), 
                     heat.pal = cols_fig3, 
                     left.label.col = "white", 
                     bottom.label.col = "white", 
                     bottom.label.text.angle = 90, 
                     X.text = round(as.matrix(betas_cor), 2), 
                     X.text.col = "white")

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

CNV.genomeplot(object = conumee_v1, directory = "./cnv_plots2/", output = "png", width = 12, height = 4, cols = cols_fig3)
CNV.genomeplot(object = conumee_v2, directory = "./cnv_plots2/", output = "png", width = 12, height = 4, cols = cols_fig3)

