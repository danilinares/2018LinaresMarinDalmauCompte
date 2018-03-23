calculate_predictions <- function (df, x_seq){
  df %>% 
    mutate(prob = map(same_slope, augment, 
                      newdata = x_seq, 
                      type.predict = "response")) %>% 
    dplyr::select(-data, -same_slope) %>% 
    unnest() %>% 
    mutate(duration = 10^log10_duration)
}
