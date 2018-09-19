calculate_predictions_dif <- function (df, x_seq){
  df %>% 
    mutate(prob = map(dif_slope,
                      ~tibble(.fitted = predict(., 
                                                newdata = x_seq, 
                                                type = "response")) %>% 
                        bind_cols(x_seq)
    )
    ) %>% 
    dplyr::select(-data, -dif_slope) %>% 
    unnest() %>% 
    mutate(duration = 10^log10_duration)
  
  # df %>% 
  #   mutate(prob = map(same_slope, augment, 
  #                     newdata = x_seq, 
  #                     type.predict = "response")) %>% 
  #   dplyr::select(-data, -same_slope) %>% 
  #   unnest() %>% 
  #   mutate(duration = 10^log10_duration)
}
