calculate_thresholds_r <- function(df, ...) {
  group_by <- quos(...)
  
  df %>% 
    mutate(same_slope_tidy = map(same_slope, tidy)) %>% 
    select(-data, -same_slope) %>% 
    unnest() %>% 
    select(-std.error, -statistic, -p.value) %>% 
    group_by(!!!group_by) %>% 
    spread(term, estimate) %>% 
    mutate(s1 = - session01 / log10_duration, 
           s2 = - session02 / log10_duration) %>% 
    dplyr::select(-log10_duration, -session02, -session01) %>% 
    gather(session, log_threshold, s1, s2) %>% 
    mutate(threshold = 10^log_threshold) 
}