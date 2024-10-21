#' Johnson et al. (2015) student cohesion data
#'
#' A dataset on team cohesion from a study on student teamwork, collected in 2007.
#' The students are third-year bachelor students in a French business school
#' following a mandatory five month course in venture creation with the objective 
#' to, as a team, write and present a business plan. 
#' The students were surveyed about their teamwork on four occasions (week 1, 9, 17, and 21). 
#' The measure of team cohesion uses a 7-point Likert scale ranging from 1 "Never" to 7 "Always".
#'
#'
#' @format A data frame with 881 rows and 6 variables:
#' \describe{
#'
#'\tabular{llll}{
#' [,1] \tab person  \tab numeric  \tab Participant ID within a group \cr
#'[,2] \tab group   \tab numeric  \tab Group identifier \cr
#'[,3] \tab time    \tab numeric  \tab Measurment occasion, 0-3 corresponding to week 1, 9, 17 and 21, respectively.\cr
#'[,4] \tab qs76       \tab numeric  \tab Item: "There is a strong team spirit" \cr
#'[,5] \tab qs77   \tab numeric  \tab Item: "We stick together in challenging situations"\cr
#'[,6] \tab qs78   \tab numeric  \tab Item: "We are loyal to each other" \cr
#'}
#'}
#'
#' @references
#' Johnson, A. R., Van de Schoot, R., Delmar, F., & Crano, W. D. (2015). Social influence interpretation of interpersonal processes and team performance over time using bayesian model selection. Journal of Management, 41 (2), 574–606
#'
"cohesiondat"

