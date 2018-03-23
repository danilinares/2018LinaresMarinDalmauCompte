calculate_proportions <- function(df, resp, ...) {
  resp <- enquo(resp)
  group_by <- quos(...)
  df %>%
    group_by(!!!group_by) %>%
    summarise(n = n(), k = sum(!!resp),
              prob = mean(!!resp),
              pvalue = binom.test(k, n)$p.value,
              signif = if_else(pvalue < .05, TRUE, FALSE),
              ymin = binom.test(k, n)$conf.int[1],
              ymax = binom.test(k, n)$conf.int[2]
              )
}
