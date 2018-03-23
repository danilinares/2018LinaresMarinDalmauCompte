library(tidyverse)
library(broom)
library(psyphy)
library(cowplot)
library(quickpsy)
#library(rlang)
#library(stringr)

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
  mutate(platform = as.character(platform)) #%>% 
 # select(session, platform, participant, duration, size, direction, correct, contrast) 


dat_resp_crt <- quickreadfiles(path = "data", 
                               participant =str_pad(1:13, 2, pad = "0"),
                               platform = c("CRT"), 
                               session = c("01", "02")) %>%
  dplyr::select(session, platform, participant, duration, 
                size, direction, correct, contrast) %>% 
  mutate(size = if_else(size == 90, 4, 1),
         platform = as.character(platform)) #%>% 
  #select(session, platform, participant,duration, size, direction, correct, contrast)


dat_resp <- dat_resp_crt %>% 
  bind_rows(dat_resp_ipad) %>% 
  mutate(size = if_else(size == 1, "Small", "Large"),
         duration = round(duration, 5)) # to match iPad and CRT

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
thresholds <- calculate_thresholds(model_same_slope)

### Parametric bootstrap ####
predictions_n <- prob %>% 
  group_by(participant, platform) %>% 
  summarise(n = first(n)) %>% 
  left_join(predictions) 

# prob_samples <- tibble(sample = 1:B, prob = list(predictions_n)) %>% 
#   unnest() %>% 
#   group_by(participant, platform, sample) %>% 
#   nest() %>% 
#   mutate(data = map(data, . %>% 
#                       rowwise() %>% 
#                       mutate(k = rbinom(1, size = n, prob = .fitted), 
#                              r = n -k, 
#                              prob = k /n)))

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

### Thresholds bootstrap #### ESO VA LENTO QUE FLIPAS!!
thresholds_boot <- model_same_slope_boot %>% 
  calculate_thresholds()

# checking bootstrap thresholds
thresholds_boot %>% 
  group_by(participant, platform, size) %>% 
  summarise(max_threshold = max(threshold))

# hay dos infinitos para large que van acompañados de .2 para el pequeño
thresholds_boot %>% 
  group_by(participant, platform, size) %>% 
  filter(threshold > .32) %>% 
  summarise(n = n())

#quitamos samples de 03 y 05 (no muchas )
thresholds_boot_ok <- thresholds_boot %>% 
  filter(threshold < .32) 

ggplot(data = thresholds_boot_ok, 
       aes(x = threshold, fill = size)) +
  facet_grid(platform ~ participant) +
  geom_histogram()


