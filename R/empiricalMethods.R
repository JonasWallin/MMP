##
#'r.emperical
#'
#' Emperical estimate of consensus parameter
#' D: 2020-03-18
#' @param y     (n x 1) the data
#' @param group (n x 1) which group the data belongs to
#' @param time  (n x 1) time point which the data is record on
#'
##
r.emperical <- function(y, group, time){
  data <- data.frame(y=y,
                     group = group,
                     time  = time )
  data_gt <-data %>% dplyr::group_by(group, time) %>%
                     dplyr::summarise(v = var(y)) %>%
                     dplyr::ungroup()
  r      <-   data_gt %>%
              dplyr::group_by(time) %>%
              dplyr::summarise(r = mean(v)) %>%
              dplyr::ungroup()
  r <- data.frame(r)
  r$r <- sqrt(r$r/r$r[1])
  return(r)
}
