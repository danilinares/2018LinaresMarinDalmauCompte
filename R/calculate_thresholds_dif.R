calculate_thresholds_dif <- function(df, ...) {
  group_by <- quos(...)
  
  df %>% 
    mutate(dif_slope_tidy = map(dif_slope, tidy)) %>% 
    select(-data, -dif_slope) %>% 
    unnest() %>% 
    select(-std.error, -statistic, -p.value) %>% 
    group_by(!!!group_by) %>% 
    spread(term, estimate) %>% 
    mutate(Small = - sizeSmall / log10_duration, 
           Large = - sizeLarge / log10_duration) %>% 
    dplyr::select(-log10_duration, -sizeLarge, -sizeSmall) %>% 
    gather(size, log_threshold, Small, Large) %>% 
    mutate(threshold = 10^log_threshold) 
}