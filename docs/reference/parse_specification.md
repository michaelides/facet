# Parse Specification for Object or Error Components

Parses a specification for object of measurement or error components.
Accepts:

- A single character string: "p"

- A character vector: c("p", "p:d")

- A formula: obj ~ p + p:d (LHS is ignored for error spec, or can be
  used for object name)

- A one-sided formula: ~ p + p:d (extracts variables from RHS)

## Usage

``` r
parse_specification(x)
```

## Arguments

- x:

  Specification (character, character vector, or formula).

## Value

A character vector of component names.
