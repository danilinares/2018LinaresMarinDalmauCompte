calculate_thresholds <- function(df) {
  df %>% 
    mutate(thresholds = map(same_slope, 
                            . %>% tidy() %>% 
                              dplyr::select(term, estimate) %>% 
                              spread(term, estimate) %>% 
                              mutate(Small = - sizeSmall / log10_duration, 
                                     Large = - sizeLarge / log10_duration) %>% 
                              dplyr::select(Small, Large) %>% 
                              gather(size, log_threshold) %>% 
                              mutate(threshold = 10^log_threshold)
    )) %>% 
    dplyr::select(-data, -same_slope) %>% 
    unnest()  
}