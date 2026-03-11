#' Example Dataset for Cronbach's Alpha
#'
#' A simulated dataset containing responses for a simple person-by-item design,
#' suitable for demonstrating Cronbach's alpha and simple G-studies.
#'
#' @format A data frame with 2000 rows and 4 variables:
#' \describe{
#'   \item{Employee_ID}{Integer identifier for the employee}
#'   \item{i}{Item identifier (e.g., "Item_1")}
#'   \item{y}{Response value (numeric)}
#'   \item{p}{Person identifier (e.g., "p 1")}
#' }
#' @source Simulated data.
"alpha_example"

#' Brennan's (2001) Generalizability Theory Example Data
#'
#' A classic dataset from Brennan (2001) illustrating a Person crossed with Task
#' and Rater design (p x t x r). Commonly used to demonstrate complex 
#' generalizability theory analyses.
#'
#' @format A data frame with 120 rows and 5 variables:
#' \describe{
#'   \item{Task}{Factor with 3 levels indicating the task}
#'   \item{Person}{Factor with 10 levels indicating the person/examinee}
#'   \item{Rater}{Factor with 12 levels indicating the rater}
#'   \item{Score}{Integer score}
#'   \item{y}{Integer score (identical to Score)}
#' }
#' @source Brennan, R. L. (2001). \emph{Generalizability theory}. Springer-Verlag.
"brennan"

#' Rajaratnam, Cronbach, & Gleser (1965) Example Data
#'
#' A classic dataset illustrating a 
#' Person crossed with Subtest and Item design.
#'
#' @format A data frame with 64 rows and 4 variables:
#' \describe{
#'   \item{Person}{Factor with 8 levels indicating the person}
#'   \item{Subtest}{Factor with 3 levels indicating the subtest}
#'   \item{Item}{Factor with 4 levels indicating the item within subtest}
#'   \item{Score}{Integer score}
#' }
#' @source Rajaratnam, N., Cronbach, L. J., & Gleser, G. C. (1965). 
#' \emph{Generalizability of stratified-parallel tests.} Psychometrika, 30(1), 39-56.
"rajaratnam"
