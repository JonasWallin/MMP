#' Sherif (1935) autokinetic effect data
#'
#' A dataset containing estimates of movement length (in inches) of a
#' light in a completely dark room. Eight groups of three individuals
#' provided four estimates for a total of 96 observations. Participants
#' in groups 1-4 first made estimates individually prior to providing
#' estimates as a group (time = -1). Participants in groups 5-8 started as groups
#' and provided indvidually made estimates after group estimates (time = 3).
#'
#'
#' @format A data frame with 96 rows and 7 variables:
#' \describe{
#'
#'\tabular{llll}{
#' [,1] \tab person  \tab numeric  \tab Participant ID within a group\cr
#'[,2] \tab time    \tab numeric  \tab Measurment occasion\cr
#'[,3] \tab group   \tab numeric  \tab Group identifier \cr
#'[,4] \tab y       \tab numeric  \tab Estimate of movement length in inches \cr
#'[,5] \tab condition   \tab numeric  \tab Experimental condition for either starting individually (1) or as a group (0)\cr
#'[,6] \tab g.mean   \tab numeric  \tab Group mean at each measurment occasion \cr
#'[,7] \tab y.centered   \tab numeric  \tab Within-group mean centered values of y for each measurment occasion \cr
#'}
#'}
#'
#' @references
#' Sherif, M. (1935). A study of some social factors in perception: Chapter 3. Archives of Psychology, 27, 23- 46.
#'
"sherifdat"

