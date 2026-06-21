# Brennan's (2001) Generalizability Theory Example Data

A classic dataset from Brennan (2001) used to illustrate a Person
crossed with Task and Rater-nested-within-Task design. The design
equation is \\p \times t \times (r:t)\\: Person is crossed with Task,
and each Task has its own panel of 4 Raters. Raters 1–4 only see Task 1,
raters 5–8 only see Task 2, and raters 9–12 only see Task 3. Because
Rater is nested within Task (not crossed), the random-effects formula
must use `(1 | Rater:Task)` rather than `(1 | Rater)`. The 120 rows
correspond to 10 persons x 3 tasks x 4 raters-per-task.

## Usage

``` r
brennan
```

## Format

A data frame with 120 rows and 5 variables:

- Task:

  Factor with 3 levels indicating the task

- Person:

  Factor with 10 levels indicating the person/examinee

- Rater:

  Factor with 12 levels indicating the rater (nested within Task)

- Score:

  Integer score

- y:

  Integer score (identical to Score)

## Source

Brennan, R. L. (2001). *Generalizability theory*. Springer-Verlag.
