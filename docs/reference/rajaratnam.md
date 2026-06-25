# Rajaratnam, Cronbach, & Gleser (1965) Example Data

A classic dataset illustrating a Person crossed with Subtest and Item
design, where items are nested within subtests. The design equation is
\\p \times s \times (i:s)\\: Person is crossed with Subtest, and each
Subtest has its own panel of items. The 64 rows correspond to 8 persons
x 3 subtests with one observation per Person x Item:Subtest cell
(Subtest 1 has 2 unique items, Subtest 2 has 4 unique items, Subtest 3
has 2 unique items).

## Usage

``` r
rajaratnam
```

## Format

A data frame with 64 rows and 5 variables:

- Person:

  Factor with 8 levels indicating the person

- Subtest:

  Factor with 3 levels indicating the subtest

- Item:

  Factor with 4 levels; the original Rajaratnam et al. (1965) item
  labels, which are reused across subtests

- ItemId:

  Factor with 8 levels; the disambiguated item identifier formed as
  `paste(Subtest, Item, sep = ".")`. Items 1, 2 in Subtest 1 are "1.1",
  "1.2"; Items 1, 2, 3, 4 in Subtest 2 are "2.1", "2.2", "2.3", "2.4";
  Items 1, 2 in Subtest 3 are "3.1", "3.2".

- Score:

  Integer score

## Source

Rajaratnam, N., Cronbach, L. J., & Gleser, G. C. (1965).
*Generalizability of stratified-parallel tests.* Psychometrika, 30(1),
39-56.

## Details

Two item-related columns are provided. `Item` is the original 4-level
factor with reused labels (e.g., "Item 1" in Subtest 1 is a *different*
physical item from "Item 1" in Subtest 2, and Items 3, 4 exist only in
Subtest 2). `ItemId` is an 8-level factor that disambiguates items by
combining the Subtest and Item labels (e.g., "1.1", "1.2", "2.1", "2.2",
"2.3", "2.4", "3.1", "3.2"). The `ItemId` column makes the nesting
structure explicit and is the recommended column for the canonical
G-study formula:
`Score ~ (1 | Person) + (1 | Subtest) + (1 | ItemId) + (1 | Person:Subtest)`.
Because `ItemId` is unique within `Subtest`, the specifications
`(1 | ItemId)` and `(1 | Item:Subtest)` are mathematically equivalent.
The Person x Item:Subtest interaction is confounded with the residual
(one observation per cell) and is absorbed into the error term. The
original `Item` column is retained for backward compatibility with code
that expects the 4-level factor.
