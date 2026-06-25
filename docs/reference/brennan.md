# Brennan's (2001) Generalizability Theory Example Data

A classic dataset from Brennan (2001) used to illustrate a Person
crossed with Task and Rater-nested-within-Task design. The design
equation is \\p \times t \times (r:t)\\: Person is crossed with Task,
and each Task has its own panel of 4 Raters. Raters 1–4 only see Task 1,
raters 5–8 only see Task 2, and raters 9–12 only see Task 3. The 120
rows correspond to 10 persons x 3 tasks x 4 raters-per-task (one
observation per Person x Rater:Task cell).

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

## Details

Because the `Rater` factor has 12 levels and these labels are unique
within each Task, `(1 | Rater)` and `(1 | Rater:Task)` are
mathematically equivalent random-effect specifications in this dataset;
the canonical G-study formula uses the un-nested form `(1 | Rater)`
because it makes the "all main effects + every estimable interaction"
structure of the model more transparent. The canonical formula is
`Score ~ (1 | Person) + (1 | Task) + (1 | Rater) + (1 | Person:Task)`;
the Person x Rater:Task interaction is confounded with the residual (one
observation per cell) and so is absorbed into the error term rather than
being estimated as a separate component.
