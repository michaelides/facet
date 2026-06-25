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
#' A classic dataset from Brennan (2001) used to illustrate a Person crossed
#' with Task and Rater-nested-within-Task design. The design equation is
#' \eqn{p \times t \times (r:t)}: Person is crossed with Task, and each Task
#' has its own panel of 4 Raters. Raters 1--4 only see Task 1, raters 5--8
#' only see Task 2, and raters 9--12 only see Task 3. The 120 rows
#' correspond to 10 persons x 3 tasks x 4 raters-per-task (one observation
#' per Person x Rater:Task cell).
#'
#' Because the \code{Rater} factor has 12 levels and these labels are unique
#' within each Task, \code{(1 | Rater)} and \code{(1 | Rater:Task)} are
#' mathematically equivalent random-effect specifications in this dataset;
#' the canonical G-study formula uses the un-nested form \code{(1 | Rater)}
#' because it makes the "all main effects + every estimable interaction"
#' structure of the model more transparent. The canonical formula is
#' \code{Score ~ (1 | Person) + (1 | Task) + (1 | Rater) + (1 | Person:Task)};
#' the Person x Rater:Task interaction is confounded with the residual (one
#' observation per cell) and so is absorbed into the error term rather than
#' being estimated as a separate component.
#'
#' @format A data frame with 120 rows and 5 variables:
#' \describe{
#'   \item{Task}{Factor with 3 levels indicating the task}
#'   \item{Person}{Factor with 10 levels indicating the person/examinee}
#'   \item{Rater}{Factor with 12 levels indicating the rater (nested within Task)}
#'   \item{Score}{Integer score}
#'   \item{y}{Integer score (identical to Score)}
#' }
#' @source Brennan, R. L. (2001). \emph{Generalizability theory}. Springer-Verlag.
"brennan"

#' Rajaratnam, Cronbach, & Gleser (1965) Example Data
#'
#' A classic dataset illustrating a Person crossed with Subtest and Item
#' design, where items are nested within subtests. The design equation is
#' \eqn{p \times s \times (i:s)}: Person is crossed with Subtest, and each
#' Subtest has its own panel of items. The 64 rows correspond to 8
#' persons x 3 subtests with one observation per Person x Item:Subtest
#' cell (Subtest 1 has 2 unique items, Subtest 2 has 4 unique items,
#' Subtest 3 has 2 unique items).
#'
#' Two item-related columns are provided. \code{Item} is the original
#' 4-level factor with reused labels (e.g., "Item 1" in Subtest 1 is a
#' *different* physical item from "Item 1" in Subtest 2, and Items 3, 4
#' exist only in Subtest 2). \code{ItemId} is an 8-level factor that
#' disambiguates items by combining the Subtest and Item labels (e.g.,
#' "1.1", "1.2", "2.1", "2.2", "2.3", "2.4", "3.1", "3.2"). The
#' \code{ItemId} column makes the nesting structure explicit and is the
#' recommended column for the canonical G-study formula:
#' \code{Score ~ (1 | Person) + (1 | Subtest) + (1 | ItemId) + (1 | Person:Subtest)}.
#' Because \code{ItemId} is unique within \code{Subtest}, the
#' specifications \code{(1 | ItemId)} and \code{(1 | Item:Subtest)} are
#' mathematically equivalent. The Person x Item:Subtest interaction is
#' confounded with the residual (one observation per cell) and is
#' absorbed into the error term. The original \code{Item} column is
#' retained for backward compatibility with code that expects the
#' 4-level factor.
#'
#' @format A data frame with 64 rows and 5 variables:
#' \describe{
#'   \item{Person}{Factor with 8 levels indicating the person}
#'   \item{Subtest}{Factor with 3 levels indicating the subtest}
#'   \item{Item}{Factor with 4 levels; the original Rajaratnam et al. (1965)
#'     item labels, which are reused across subtests}
#'   \item{ItemId}{Factor with 8 levels; the disambiguated item identifier
#'     formed as \code{paste(Subtest, Item, sep = ".")}. Items 1, 2 in
#'     Subtest 1 are "1.1", "1.2"; Items 1, 2, 3, 4 in Subtest 2 are
#'     "2.1", "2.2", "2.3", "2.4"; Items 1, 2 in Subtest 3 are "3.1", "3.2".}
#'   \item{Score}{Integer score}
#' }
#' @source Rajaratnam, N., Cronbach, L. J., & Gleser, G. C. (1965).
#' \emph{Generalizability of stratified-parallel tests.} Psychometrika, 30(1), 39-56.
"rajaratnam"
