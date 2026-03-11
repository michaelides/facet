# Design: dstudy Variance Components Display Enhancement

## Overview

This document describes the design for modifying the variance components table in dstudy objects to:
1. Remove diagnostic columns (`2.5% CI`, `97.5% CI`, `SE`, `Rhat`, `Bulk_ESS`, `Tail_ESS`)
2. Add scaled variance estimates (`var_scaled`, `pct_scaled`)
3. Keep unscaled estimates (`var_unscaled`, `pct_unscaled`)

This change applies to all backends (brms, lme4, mom).

---

## Current Implementation

### Variance Components Table Structure

#### brms backend
Currently displays:
- `component` - Name of variance component
- `dim` - Dimension/response variable
- `estimate` - Variance estimate
- `pct` - Percentage of total variance
- `2.5% CI`, `97.5% CI` - Credible intervals
- `SE` - Standard error
- `Rhat`, `Bulk_ESS`, `Tail_ESS` - MCMC diagnostics

#### lme4 backend
Currently displays:
- `component`, `dim`, `estimate`, `pct`
- `2.5% CI`, `97.5% CI` - Confidence intervals (when CIs requested)

#### mom backend
Currently displays:
- `component`, `dim`, `estimate`, `pct`
- `2.5% CI`, `97.5% CI` - Confidence intervals
- `SE` - Standard error

### Current Scaling Behavior

The `calculate_dstudy_variance()` function in `R/coefficients.R` already handles scaled variance calculations when:
1. `n_provided = TRUE` (user specifies sample sizes)
2. `aggregation` is specified

Currently, when `n` is not provided by the user, the variance components are NOT scaled.

---

## Proposed Changes

### 1. New Column Structure for dstudy

**For all backends:**

| Column | Description | Type |
|--------|-------------|------|
| `component` | Name of variance component | character |
| `var_unscaled` | Unscaled variance estimate (from G-study) | numeric |
| `pct_unscaled` | Percentage of total unscaled variance | numeric |
| `var_scaled` | Scaled variance (divided by D-study sample sizes) | numeric |
| `pct_scaled` | Percentage of total scaled variance | numeric |
| `dim` | Dimension/response variable | character |

**Removed columns:**
- `2.5% CI`, `97.5% CI` (lower/upper)
- `SE` / `error`
- `Rhat`, `Bulk_ESS`, `Tail_ESS` (brms-specific)

### 2. Scaling Logic

Based on existing `calculate_dstudy_variance()` logic:

**For each variance component:**

1. **Object of measurement** (e.g., `person`):
   - `var_scaled = var_unscaled` (not scaled)

2. **Main effects of non-object facets** (e.g., `item`, `rater`):
   - `var_scaled = var_unscaled / n_facet`

3. **Interactions with object** (e.g., `person:item`, `person:rater`):
   - `var_scaled = var_unscaled / n_facet` (where facet is the non-object component)

4. **Interactions without object** (e.g., `item:rater`):
   - `var_scaled = var_unscaled / (n_facet1 * n_facet2)`

5. **Residual**:
   - `var_scaled = var_unscaled / (product of all non-object facet sample sizes)`

**Percentage calculations:**
- `pct_unscaled = (var_unscaled / sum(var_unscaled)) * 100`
- `pct_scaled = (var_scaled / sum(var_scaled)) * 100`

### 3. Example

With `object = "person"`, `n = list(item = 10, rater = 3)`:

| component | var_unscaled | pct_unscaled | var_scaled | pct_scaled |
|-----------|--------------|--------------|------------|------------|
| person | 2.500 | 41.7% | 2.500 | 71.4% |
| item | 1.800 | 30.0% | 0.180 | 5.1% |
| person:item | 0.500 | 8.3% | 0.050 | 1.4% |
| rater | 0.600 | 10.0% | 0.200 | 5.7% |
| Residual | 0.600 | 10.0% | 0.020 | 0.6% |

---

## Implementation Strategy

### File Modifications

| File | Function | Changes |
|------|----------|---------|
| `R/dstudy.R` | `dstudy()` | Always compute both `var_unscaled` and `var_scaled` in variance_components |
| `R/variance-components.R` | `summarize_vc()` | Add method to detect and handle dstudy variance components |
| `R/methods.R` | `print.dstudy()`, `summary.dstudy()` | Use modified summary logic for display |
| `R/coefficients.R` | `calculate_dstudy_variance()` | Ensure both scaled and unscaled are returned |

### Key Changes

#### 1. Modify `dstudy()` function

**Current behavior:**
- When `n_provided = FALSE`: uses unscaled variance components
- When `n_provided = TRUE`: computes scaled variance components
- When `aggregation` specified: computes both unscaled and scaled

**New behavior:**
- Always compute both `var_unscaled` and `var_scaled`
- Always compute both `pct_unscaled` and `pct_scaled`
- Store in `result$variance_components`

**Implementation approach:**
```r
# In dstudy() function, after getting vc from G-study:
# Always compute scaled variance using D-study sample sizes
d_vc <- calculate_dstudy_variance(
  vc = vc,
  n = n,
  object = object,
  aggregation = aggregation,
  n_provided = TRUE,  # Force scaling calculation
  residual_is = residual_is_effective,
  facet_n = gstudy_obj$facet_n
)

# Add unscaled columns
d_vc$var_unscaled <- vc$var
d_vc$pct_unscaled <- vc$pct

# Rename scaled columns for clarity
d_vc$var_scaled <- d_vc$var
d_vc$pct_scaled <- d_vc$pct

# Remove temporary columns
d_vc$var <- NULL
d_vc$pct <- NULL
```

#### 2. Create `summarize_vc_dstudy()` function

New function specifically for dstudy variance components:

```r
summarize_vc_dstudy <- function(vc, digits = 3) {
  # Check if this is a dstudy variance components tibble
  required_cols <- c("component", "var_unscaled", "pct_unscaled",
                     "var_scaled", "pct_scaled", "dim")
  missing_cols <- setdiff(required_cols, names(vc))
  if (length(missing_cols) > 0) {
    stop("dstudy variance components tibble is missing columns: ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  # Select and format columns
  summary_df <- vc[, required_cols]

  # Round numeric columns
  summary_df$var_unscaled <- round(summary_df$var_unscaled, digits)
  summary_df$pct_unscaled <- round(summary_df$pct_unscaled, digits)
  summary_df$var_scaled <- round(summary_df$var_scaled, digits)
  summary_df$pct_scaled <- round(summary_df$pct_scaled, digits)

  summary_df
}
```

#### 3. Modify `summarize_vc()` to detect dstudy objects

Add detection logic:

```r
summarize_vc <- function(vc, digits = 3, scale = c("variance", "sd")) {
  # Check if this is a dstudy variance components tibble
  if (all(c("var_unscaled", "var_scaled") %in% names(vc))) {
    return(summarize_vc_dstudy(vc, digits = digits))
  }

  # Original implementation for gstudy variance components
  # ... (existing code)
}
```

#### 4. Update print and summary methods

**In `print.dstudy()`:**

Replace the variance components display section:

```r
# Variance Components:
cat("Variance Components:\n")
vc_summary <- summarize_vc_dstudy(x$variance_components, digits = digits)
print(vc_summary, row.names = FALSE, ...)
```

**In `summary.dstudy()`:**

Similar replacement in the summary output.

---

## Edge Cases

### 1. When n is not provided

Current behavior: uses G-study sample sizes as default
New behavior: same, but always shows both unscaled and scaled estimates

### 2. Multivariate models

Ensure `dim` column is preserved and displayed correctly for each dimension.

### 3. Sweeping (multiple n values)

When `is_sweep = TRUE`, the variance_components table is the base (unscaled) version.
Need to handle this case appropriately - may not show scaled estimates for sweep.

### 4. Zero variance components

Handle gracefully if some variance components are zero or near-zero.

---

## Testing Strategy

### Unit Tests

1. **Test for brms backend:**
   - Verify `var_unscaled`, `var_scaled`, `pct_unscaled`, `pct_scaled` columns present
   - Verify scaling is correct
   - Verify no diagnostic columns (CI, SE, Rhat, ESS) in display

2. **Test for lme4 backend:**
   - Same verification as brms
   - Test with and without CIs

3. **Test for mom backend:**
   - Same verification

4. **Test scaling logic:**
   - Object component: not scaled
   - Non-object main effects: scaled by single n
   - Interactions: scaled by product of facet n's
   - Residual: scaled by all non-object facet n's

5. **Test percentages:**
   - Verify `pct_unscaled` sums to 100
   - Verify `pct_scaled` sums to 100

### Integration Tests

1. Test that `print.dstudy()` displays correctly
2. Test that `summary.dstudy()` displays correctly
3. Test with multivariate models
4. Test with sweep (multiple sample sizes)

---

## Backward Compatibility

### Potential Breaking Changes

1. **Column names change:** Users accessing `variance_components$var` will need to use `variance_components$var_unscaled` or `variance_components$var_scaled`

2. **Diagnostic columns removed:** Users expecting CI or diagnostic columns in variance_components will not find them

### Mitigation

1. **Document changes clearly** in NEWS.md
2. **Provide migration guide** for common use cases
3. **Consider deprecation period:**
   - Keep `var` and `pct` columns as aliases to `var_unscaled` and `pct_unscaled` initially
   - Add deprecation warning if user accesses them
   - Remove in future version

---

## Documentation Updates

### Files to Update

1. `man/dstudy.Rd` - Update return value documentation
2. `man/summarize_vc.Rd` - Add dstudy variance components handling
3. `vignettes/introduction.Rmd` - Update examples to show new column structure
4. `NEWS.md` - Document breaking changes
5. `README.md` - Update examples if needed

---

## Implementation Checklist

- [ ] Modify `calculate_dstudy_variance()` to always return both scaled and unscaled
- [ ] Update `dstudy()` to use modified variance calculation
- [ ] Create `summarize_vc_dstudy()` function
- [ ] Modify `summarize_vc()` to detect dstudy objects
- [ ] Update `print.dstudy()` method
- [ ] Update `summary.dstudy()` method
- [ ] Add unit tests for all backends
- [ ] Add integration tests
- [ ] Update documentation
- [ ] Add NEWS entry
- [ ] Test with real-world examples
