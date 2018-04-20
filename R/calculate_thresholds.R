calculate_thresholds <- function(df, ...) {
  group_by <- quos(...)
  
  df %>% 
    mutate(same_slope_tidy = map(same_slope, tidy)) %>% 
    select(-data, -same_slope) %>% 
    unnest() %>% 
    select(-std.error, -statistic, -p.value) %>% 
    group_by(!!!group_by) %>% 
    spread(term, estimate) %>% 
    mutate(Small = - sizeSmall / log10_duration, 
           Large = - sizeLarge / log10_duration) %>% 
    dplyr::select(-log10_duration, -sizeLarge, -sizeSmall) %>% 
    gather(size, log_threshold, Small, Large) %>% 
    mutate(threshold = 10^log_threshold)
  # df %>% 
  #   mutate(thresholds = map(same_slope, 
  #                           . %>% tidy() %>% 
  #                             dplyr::select(term, estimate) %>% 
  #                             spread(term, estimate) %>% 
  #                             mutate(Small = - sizeSmall / log10_duration, 
  #                                    Large = - sizeLarge / log10_duration) %>% 
  #                             dplyr::select(Small, Large) %>% 
  #                             gather(size, log_threshold) %>% 
  #                             mutate(threshold = 10^log_threshold)
  #   )) %>% 
  #   dplyr::select(-data, -same_slope) %>% 
  #   unnest()  
}