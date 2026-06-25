## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(facet)
data(brennan)
head(brennan)

## ----gstudy-------------------------------------------------------------------
g_obj <- gstudy(Score ~ (1 | Person) + (1 | Task) + (1 | Rater) +
  (1 | Person:Task),
  data = brennan)

# View the variance components
print(g_obj)

## ----plot---------------------------------------------------------------------
plot(g_obj)

## ----dstudy-------------------------------------------------------------------
# Project a design with 3 Tasks and 4 Raters
d_obj <- dstudy(g_obj, n = list(Task = 3, Rater = 4))

print(d_obj)

