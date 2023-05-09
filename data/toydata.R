#'
#' Small, artificially-generated toy data that contains a count variable and a grouping variable
#'
#' @docType data
#'
#' @usage data(toydata)
#'
#' @format An object of class \code{"data.frame"}
#' \describe{
#'  \item{grp}{sample grouping variable)}
#'  \item{nb_count}{A count variable (integer format) that follows a negative binomial distribution
#'  \item{norm_count}{A variable that follows a normal distribution}}
#' }
#' @references This data set was artificially created for the whichdist package.
#' @keywords datasets
#' @examples
#'
#' data(toydata)
#' head(toydata)
#'
"toydata"

grp_1 = data.frame(rep(1,500),
                   rnbinom(n = 500, mu = 6*1.03, size = .2),
                   rnorm(n = 500, mean = 6*1.03, sd = 1))
names(grp_1) = c("grp", "nb_count", "norm_count")

grp_0 = data.frame(rep(0,500),
                   rnbinom(n = 500, mu = 6, size = .2),
                   rnorm(n = 500, mean = 6, sd = 1))
names(grp_0) = c("grp", "nb_count", "norm_count")

toydata = rbind(grp_1, grp_0)

save(toydata, file = "~/Documents/whichdist/data/toydata.RData")
