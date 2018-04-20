library(tidyverse)
library(broom)
library(psyphy)
library(cowplot)
library(quickpsy)

set.seed(999) 

### Load functions and parameters ####
list.files("R", full.names = TRUE) %>% walk(source)
source("parameters.R")

### Read data ####
dat_resp_ipad <- quickreadfiles(path = "data", 
                                participant = str_pad(1:13, 2, pad = "0"),
                                platform = c("iPad"), 
                                session = c("01", "02")) %>%
  dplyr::select(session, platform, participant, Duration, 
                Size, Direction, Correct, Contrast) %>% 
  rename(duration = Duration, size = Size, contrast = Contrast, 
         direction = Direction, correct = Correct) %>% 
  mutate(platform = as.character(platform)) 

dat_resp_crt <- quickreadfiles(path = "data", 
                               participant =str_pad(1:13, 2, pad = "0"),
                               platform = c("CRT"), 
                               session = c("01", "02")) %>%
  dplyr::select(session, platform, participant, duration, 
                size, direction, correct, contrast) %>% 
  mutate(size = if_else(size == 90, 4, 1),
         platform = as.character(platform)) 

dat_resp <- dat_resp_crt %>% 
  bind_rows(dat_resp_ipad) %>% 
  mutate(size = if_else(size == 1, "Small", "Large"),
         duration = round(duration, 5)) # to match iPad and CRT

# Eliminating the shortest duration for the iPad (see paper)
dat_resp <- dat_resp %>% filter(!(duration == 0.0100 & platform == "iPad"))

### Calculate proportions ####
prob <- calculate_proportions(dat_resp, correct, duration, 
                              platform, size, participant) %>% 
  mutate(r = n - k, log10_duration = log10(duration)) %>% 
  ungroup()

### Models #### 
models <- prob %>% 
  group_by(participant, platform) %>% 
  nest() %>% 
  mutate(
    dif_slope = map(data, 
                    ~glm(cbind(k, r) ~ size / log10_duration - 1, 
                         data = ., family = binomial(mafc.logit(2)))),
    same_slope = map(data, 
                     ~glm(cbind(k, r) ~ size + log10_duration - 1, 
                          data = ., family = binomial(mafc.logit(2))))
  )

model_comparisons <- models %>% 
  group_by(participant, platform) %>% 
  mutate(anov2 = map2(dif_slope, same_slope, anova, test = "Chisq"),
         p.value = map_dbl(anov2, ~.$`Pr(>Chi)`[2]),
         significant = p.value < alpha) 
  
model_comparisons %>% filter(p.value < alpha)

# In 3 models, the p.value is NA. In these models glm does not fit the data well 
# as the deviance of them is larger for the model that includes  
# more parameters (different slope models). We fitted these models by direct 
# maximization of the likelihood and found that the model with identical slopes 
# was better (not shown) that the models with the same slope. These are the 3
# models: 
models %>% 
  mutate(dev_dif_slope = map_dbl(dif_slope, "deviance"),
          dev_same_slope = map_dbl(same_slope, "deviance"),
          poor_fits = dev_dif_slope > dev_same_slope) %>% filter(poor_fits)

model_same_slope <- models %>% dplyr::select(-dif_slope)

#### Duration sequences ####
log10_duration_seq <- prob %>% 
  distinct(size) %>% 
  crossing(tibble(log10_duration = seq(log_duration_min, 
                                       log_duration_max, length.out = 100)))
### Psychometric functions ####
predictions <- model_same_slope %>% 
  calculate_predictions(prob %>% distinct(size, log10_duration))

psychometric_functions <- model_same_slope %>% 
  calculate_predictions(log10_duration_seq)

### Thresholds ####
thresholds <- model_same_slope %>% 
  calculate_thresholds(participant, platform)

### Parametric bootstrap ####
predictions_n <- prob %>% 
  group_by(participant, platform) %>% 
  summarise(n = first(n)) %>% 
  left_join(predictions) 

prob_samples <- tibble(sample = 1:B, prob = list(predictions_n)) %>% 
  unnest() %>% 
  group_by(participant, platform, sample) %>% 
  mutate(k = rbinom(n(), size = n, prob = .fitted), 
         r = n - k, 
         prob = k /n)

model_same_slope_boot <- prob_samples %>% 
  group_by(participant, platform, sample) %>% 
  nest() %>% 
  mutate(
    same_slope = map(data, 
                     ~glm(cbind(k, r) ~ size + log10_duration - 1, 
                          data = ., family = binomial(mafc.logit(2))))
  )

### Thresholds bootstrap #### 
thresholds_boot <- model_same_slope_boot %>% 
  calculate_thresholds(participant, platform, sample)

# Notice that some bootstrap samples produce Inf threshold 
# for participants 3 and 5
thresholds_boot %>% 
  group_by(participant, platform, size) %>% 
  filter(is.infinite(threshold)) %>% 
  summarise(n = n())

### Confidence intervals #### 
conf_int <- thresholds_boot %>% 
  group_by(participant, platform, size) %>% 
  summarise(threshold_min = quantile(threshold, alpha /2),
            threshold_max = quantile(threshold, 1 - alpha /2)) 

# includding a jitter for the plot
conf_int_size <- conf_int %>% 
  mutate(prob = if_else(size == "Large", .75, .73))

conf_int_platform <- conf_int %>% 
  mutate(prob = if_else(platform == "CRT", .75, .73))

differences_size <- thresholds_boot %>%    
  dplyr::select(-log_threshold) %>% 
  spread(size, threshold) %>% 
  mutate(dif = Large - Small) %>% 
  group_by(participant, platform) %>% 
  summarise(dif_min = quantile(dif, alpha /2), 
            dif_max = quantile(dif, 1 - alpha /2), 
            significant = if_else(dif_min * dif_max > 0, "*",""))

differences_platform <- thresholds_boot %>% 
  dplyr::select(-log_threshold) %>% 
  spread(platform, threshold) %>% 
  mutate(dif = CRT - iPad) %>% 
  group_by(participant, size) %>% 
  summarise(dif_min = quantile(dif, alpha /2), 
            dif_max = quantile(dif, 1 - alpha /2), 
            significant = if_else(dif_min * dif_max > 0, "*",""))

### Plot psychometric functions ####
p_size <- ggplot(prob) +
  facet_grid(platform~ participant, scales = "free") +
  geom_line(data = psychometric_functions, size = .5 * size_line, 
            aes(x = duration,y = .fitted,  color = size)) +
  
  geom_point(aes(x = duration, y = prob, 
                 color = size, shape = size), size = size_point) +
  geom_segment(data = conf_int_size, 
               aes(x = threshold_min, xend = threshold_max, 
                   y = prob, yend = prob, color = size),
               size = size_line) +
  geom_text(data = differences_size, aes(label = significant, 
                                         x = .01, y = .9)) +
  scale_color_brewer(labels = name_size, palette = "Set1") +
  scale_fill_brewer(labels = name_size, palette = "Set1") +
  scale_shape_discrete(labels = name_size) +
  scale_x_log10(breaks = c(.01, .04, .16), labels = c(".01", ".04", ".16")) +
  scale_y_continuous(breaks = seq(0, 1,.5),
                     limits = c(.1, 1)) +
  coord_cartesian(xlim = c(0.005, .3)) +
  labs(x = label_duration, y = label_proportion,
       color = label_size, shape = label_size, fill = label_size) +
  theme(legend.position = "top",
        legend.text = element_text(size = 9))

p_platform <- ggplot(prob) +
  facet_grid(size ~ participant, scales = "free") +
  geom_line(data = psychometric_functions, size = .5 * size_line, 
            aes(x = duration,y = .fitted,  color = platform)) +
  geom_point(aes(x = duration, y = prob, 
                 color = platform, shape = platform), size = size_point) +
  geom_segment(data = conf_int_platform, 
               aes(x = threshold_min, xend = threshold_max, 
                   y = prob, yend = prob, color = platform),
               size = size_line) +
  geom_text(data = differences_platform, aes(label = significant, 
                                             x = .01, y = .9)) +
  scale_color_brewer(labels = name_size, palette = "Dark2") +
  scale_fill_brewer(labels = name_size, palette = "Dark2") +
  scale_shape_discrete(labels = name_size) +
  scale_x_log10(breaks = c(.01, .04, .16), labels = c(".01", ".04", ".16")) +
  scale_y_continuous(breaks = seq(0, 1,.5),
                     limits = c(.1, 1)) +
  coord_cartesian( xlim = c(0.005, .3)) +
  labs(x = label_duration, y = label_proportion,
       color = label_platform, shape = label_platform, fill = label_platform) +
  theme(legend.position = "top",
        legend.text = element_text(size = 9))

p_psycho <- plot_grid(p_size, p_platform, ncol = 1, labels = "AUTO")

ggsave("figures/figure2.pdf", p_psycho, width = two_columns_width, height = 5) 


### Goodness of fit ####
likelihoods <- prob %>% 
  left_join(predictions %>% 
              filter(!(log10_duration == -2 & platform == "iPad")), 
            by = c("log10_duration", "platform", "size", "participant")) %>% 
  mutate(log_like_saturated_for_each_dur = dbinom(k, n, prob, log = TRUE),
         log_like_model_for_each_dur = dbinom(k, n, .fitted, log = TRUE)) %>% 
  group_by(participant, platform) %>% 
  summarise(log_like_saturated = sum(log_like_saturated_for_each_dur), 
            log_like_model = sum(log_like_model_for_each_dur)) %>% 
  group_by(participant, platform) %>% 
  transmute(log_lrt = log_like_model - log_like_saturated)

likelihoods_boot <- prob_samples %>% 
  mutate(log_like_saturated_for_each_dur = dbinom(k, n, prob, log = TRUE),
         log_like_model_for_each_dur = dbinom(k, n, .fitted, log = TRUE)) %>% 
  group_by(sample, participant, platform) %>% 
  summarise(log_like_saturated = sum(log_like_saturated_for_each_dur), 
            log_like_model = sum(log_like_model_for_each_dur)) %>% 
  group_by(sample, participant, platform) %>% 
  transmute(log_lrt_boot = log_like_model - log_like_saturated) 

all_likelihoods <- likelihoods_boot %>% 
  left_join(likelihoods)

all_likelihoods %>% 
  group_by(participant, platform) %>% 
  summarise(p_value = mean(log_lrt_boot < log_lrt)) %>% 
  filter(p_value < alpha)


### Correlations across thresholds ####
thresholds_long <- thresholds %>% 
  dplyr::select(-log_threshold) %>% 
  spread(platform, threshold)  

thresholds_mean_ci <- thresholds %>% 
  group_by(size, platform) %>% 
  nest() %>% 
  mutate(t = map(data, ~t.test(.$log_threshold, conf.level = 1 - alpha) 
                 %>% tidy())) %>% 
  unnest(t) %>% 
  select(size, platform, estimate, conf.low, conf.high) %>% 
  group_by(size, platform) %>% 
  mutate_all(~10^.)

thresholds_mean_ci_crt <- thresholds_mean_ci %>% 
  ungroup() %>% 
  filter(platform == "CRT") %>% 
  dplyr::select(-platform) %>% 
  rename(CRT = estimate, CRTmin = conf.low, CRTmax = conf.high) 

thresholds_mean_ci_ipad <- thresholds_mean_ci %>% 
  ungroup() %>% 
  filter(platform == "iPad") %>% 
  dplyr::select(-platform) %>% 
  rename(iPad = estimate, iPadmin = conf.low, iPadmax = conf.high) 

thresholds_mean_ci_long <- thresholds_mean_ci_crt %>% 
  left_join(thresholds_mean_ci_ipad)

log_thresholds_long <- thresholds %>% 
  dplyr::select(-threshold) %>% 
  spread(platform, log_threshold)  

linear_model <- log_thresholds_long %>% 
  group_by(size) %>% 
  nest() %>% 
  mutate(cor = map(data, ~cor.test(.$CRT, .$iPad, conf.level = 1  - alpha)), 
         model = map(data, ~lm(iPad ~ CRT, data = . )), 
         ci = map(model, confint, level = 1  - alpha)) 

p_cor_size <- ggplot(thresholds_long) +
  geom_abline(color = "grey", size = size_line) +
  geom_smooth(aes(x = CRT, y = iPad, color = size), 
              method = "lm", se = FALSE, size = size_line) +
  geom_point(aes(x = CRT, y = iPad, color = size, 
                 shape = size), size = size_point_cor) +
  geom_point(data = thresholds_mean_ci_long,
             aes(x = CRT, y = iPad, shape = size), color = "black",
             size = size_point_cor, show.legend = FALSE) +
  geom_errorbarh(data = thresholds_mean_ci_long , height = .1, size = size_line, 
                 aes(x = CRT, xmin = CRTmin, xmax = CRTmax, 
                     y = iPad, group = size)) +
  geom_errorbar(data = thresholds_mean_ci_long , width = .1, size = size_line,
                aes( ymin = iPadmin, ymax = iPadmax, 
                     x = CRT, group = size)) +
  geom_point(data = thresholds_mean_ci_long,
             aes(x = CRT, y = iPad, shape = size, color = size), 
             size = .5 *size_point_cor) +
  scale_color_brewer(labels = name_size, palette = "Set1") +
  scale_shape_discrete(labels = name_size) +
  coord_equal() +
  scale_x_log10(breaks = c(.01, .04, .16),
                labels = c(".01", ".04", ".16"),
                limits = c(.008,.3)) +
  scale_y_log10(breaks = c(.01, .04, .16),
                labels = c(".01", ".04", ".16"),
                limits = c(.008,.3)) +
  labs(x = label_crt, y = label_ipad,
       color = label_size, shape = label_size) +
  theme(legend.text = element_text(size = 9))

### Anova #####
aov(log_threshold ~ platform * size  + 
      Error(participant / (platform * size)), data = thresholds) %>% 
  summary()

### Correlations across suppression ####
ss <- thresholds %>% 
  dplyr::select(-threshold) %>% 
  spread(size, log_threshold) %>% 
  mutate(ss = Large - Small) %>% 
  dplyr::select(-Large, -Small) 

ss_long <- ss %>% 
  spread(platform, ss)

ss_mean_ci <- ss %>% 
  group_by(platform) %>% 
  nest() %>% 
  mutate(t = map(data, ~t.test(.$ss, conf.level = 1 - alpha) %>% tidy())) %>% 
  unnest(t) %>% 
  select(platform, estimate, conf.low, conf.high)

ss_mean_ci_crt <- ss_mean_ci %>% 
  filter(platform == "CRT") %>% 
  dplyr::select(-platform) %>% 
  rename(CRT = estimate, CRTmin = conf.low, CRTmax = conf.high)

ss_mean_ci_ipad <- ss_mean_ci %>% 
  filter(platform == "iPad") %>% 
  dplyr::select(-platform) %>% 
  rename(iPad = estimate, iPadmin = conf.low, iPadmax = conf.high)

ss_mean_ci_long <- ss_mean_ci_crt %>% 
  bind_cols(ss_mean_ci_ipad)

cor.test(ss_long$CRT, ss_long$iPad)

p_ss <- ggplot(ss_long, aes(x = CRT, y = iPad)) +
  geom_abline(color = "grey", size = size_line) +
  geom_point(size = size_point_cor, color = "#2ca25f") +
  geom_smooth(method = "lm", se = FALSE, size = size_line, 
              color = "#2ca25f") +
  geom_point(data = ss_mean_ci_long,
             aes(x = CRT, y = iPad), color = "black",
             size = size_point_cor, show.legend = FALSE) +
  geom_errorbarh(data = ss_mean_ci_long , height = .1, size = size_line, 
                 aes(x = CRT, xmin = CRTmin, xmax = CRTmax, 
                     y = iPad)) +
  geom_errorbar(data = ss_mean_ci_long , width = .1, size = size_line,
                aes( ymin = iPadmin, ymax = iPadmax, 
                     x = CRT)) +
  geom_point(data = ss_mean_ci_long,
             aes(x = CRT, y = iPad), color = "#2ca25f",
             size = .5 * size_point_cor, show.legend = FALSE) +
  geom_line(aes(lty = "")) +
  scale_linetype_manual(values = 0) +
  coord_equal() +
  scale_color_manual(values = "black")  +
  scale_x_continuous(breaks = seq(0, 1, .25), 
                     limits = c(0, 1)) +
  scale_y_continuous(breaks = seq(0, 1, .25),
                     limits = c(0, 1)) +
  labs(x = label_crt_ss, y = label_ipad_ss) +
  theme(legend.title = element_blank())

### t-test ####
t.test(ss_long$CRT, ss_long$iPad, paired = TRUE)

p_cor <- plot_grid(p_cor_size, 
                   p_ss, 
                   ncol = 1, 
                   labels = "AUTO",
                   align = "v")

ggsave("figures/figure3.pdf", p_cor, width = single_column_width, height = 3.5) 



