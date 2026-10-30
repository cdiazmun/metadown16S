# Rarefaction Analysis for 16S rRNA Amplicon Data: A Practical Guide

## Contents

1. What rarefaction is
2. The two objectives of our analysis
3. Step-by-step guide with code
4. Reporting checklist
5. References and further reading

---

## 1. What rarefaction is

Amplicon sequencing does not sample every organism in a sample equally. The number of reads obtained from a sample (its **sequencing depth**) depends on library preparation, pooling, and sequencing run, and it is rarely the same across samples. Diversity estimates depend on depth: a sample sequenced to 300,000 reads will usually show more ASVs than an identical community sequenced to 5,000 reads, simply because more rare variants are detected.

**Rarefaction** addresses this by asking what diversity would look like if every sample had been sequenced to the same depth. It works in two ways:

- **Rarefaction curves.** For each sample, the number of observed ASVs (or genera, or species) is computed after subsampling 1, 2, 3, ... up to *N* reads. Plotting observed richness against reads gives a curve. A curve that flattens (plateaus) indicates that additional sequencing is unlikely to reveal many new taxa, so the sample is said to be **saturated**. A curve still rising at the sample's maximum depth means that the sample is **under-sampled**.
- **Rarefying the count table.** Each sample is randomly subsampled, without replacement, to a fixed depth (the cutoff). Samples with fewer reads than the cutoff are discarded. The result is a count table in which all samples have the same total, which allows fair comparison of alpha and beta diversity.

Rarefaction has been debated. McMurdie and Holmes (2014) argued that discarding reads wastes information and that rarefying is not appropriate for every analysis, in particular differential abundance testing. The current practice, which we follow, is:

- Use **rarefied** data for **alpha diversity** and **beta diversity**.
- Use **unrarefied** counts (with a normalization or model suited to count data) for **differential abundance** testing.

Note that rarefaction is a statistical equalization of sampling effort, not a way to recover information that was not sequenced. Richness estimates for under-sampled samples are underestimates, and this should be stated when reporting results.

### Key concepts

| Term | Meaning |
|---|---|
| ASV | Amplicon sequence variant, produced by DADA2 |
| Sequencing depth | Total number of reads in a sample |
| Rarefaction curve | Observed richness versus number of reads subsampled |
| Saturation depth | The depth at which a sample's rarefaction curve has flattened |
| Slope (rareslope) | Rate at which new ASVs appear per additional read at a given depth |
| Cutoff | The common depth to which all retained samples are rarefied |

### Sources

- Gotelli, N. J., & Colwell, R. K. (2001). Quantifying biodiversity: procedures and pitfalls in the measurement and comparison of species richness. *Ecology Letters*, 4, 379–391. https://doi.org/10.1046/j.1461-0248.2001.00230.x
- Hurlbert, S. H. (1971). The nonconcept of species diversity: a critique and alternative parameters. *Ecology*, 52, 577–586. https://doi.org/10.2307/1934145
- McMurdie, P. J., & Holmes, S. (2014). Waste not, want not: why rarefying microbiome data is inadmissible. *PLoS Computational Biology*, 10, e1003531. https://doi.org/10.1371/journal.pcbi.1003531
- Chiu, C.-H., & Chao, A. (2014). Individual-based representation and high-throughput sequencing: a unified approach to species diversity. *Ecology Letters* (see also the iNEXT package for sample-size-based rarefaction and extrapolation). https://chao.shinyapps.io/iNEXTOnline/
- vegan package documentation (`rarefy`, `rrarefy`, `rarecurve`, `rareslope`, `specnumber`): https://vegandevs.github.io/vegan/reference/
- Callahan, B. J. et al. (2016). DADA2: high-resolution sample inference from Illumina amplicon data. *Nature Methods*, 13, 581–583. https://doi.org/10.1038/nmeth.3869
- Weiss, S. et al. (2017). Normalization and microbial differential abundance strategies depend upon data characteristics. *Microbiome*, 5, 27. https://doi.org/10.1186/s40168-017-0237-y

---

## 2. The two objectives

The analysis has two separate objectives, and they use different criteria.

### Objective 1: Per-sample adequacy

*Question:* Has the sequencing depth of each sample been enough to represent the diversity present in that sample?

*Criterion:* A sample is adequate if it is saturated within its own reads, that is, its saturation depth is no greater than its total depth. On a plot of sequencing depth (x) against saturation depth (y), adequate samples fall **below** the 1:1 line.

*Action:* Drop samples on or above the line. This decision does not depend on any common cutoff.

### Objective 2: A common cutoff

*Question:* At what depth can every retained sample be rarefied so that alpha and beta diversity are computed on comparable, adequately sampled data?

*Criterion:* For a candidate cutoff *c*, a sample is usable if (a) its total depth is at least *c*, so that it can be subsampled to *c*, and (b) its saturation depth is at most *c*, so that the subsample captures its richness.

*Action:* Count retained samples across candidate cutoffs, choose a cutoff where the count stops improving or where the loss of samples is no longer worth the gain in depth, then rarefy.

Objective 1 is applied first, and objective 2 is applied to the samples that pass it.

---

## 3. Step-by-step guide

The examples use R with the `vegan` package. The input is a DADA2 sequence table with **samples as columns and ASVs as rows**, which is the layout produced by `dada2::mergeSequenceTables()` or by iterative joins across sequencing runs. Adapt the transposition (`t()`) if your layout differs.

### Step 0: Load packages and prepare the table

```r
library(vegan)

# seqtab: ASVs x samples (rows = ASVs, columns = samples)
# Full joins across runs introduce NAs for ASVs absent from a run.
# Absent ASVs must be zero counts, not missing values.
seqtab[is.na(seqtab)] <- 0
seqtab <- as.matrix(seqtab)
mode(seqtab) <- "numeric"

stopifnot(sum(is.na(seqtab)) == 0, all(seqtab >= 0))
```

Check for duplicated sequences. If the same sequence received different ASV identifiers in different runs, merge them by summing counts:

```r
length(unique(rownames(seqtab)))   # identifiers
# If you have the sequences as a vector 'seqs' (one per row):
length(unique(seqs))               # should match
counts <- rowsum(seqtab, group = seqs)   # sums rows with identical sequences
sum(seqtab) == sum(counts)               # total reads must be unchanged
```

Remove mock communities and negative controls before computing depths, since they distort the cutoff.

```r
# Example: drop control samples by name
controls <- c("MOCK", "NEG1")
seqtab <- seqtab[, !colnames(seqtab) %in% controls]
```

### Step 1: Compute per-sample depth

Depth is the column sum, because samples are columns.

```r
depth <- colSums(seqtab)
summary(depth)
```

Inspect the spread. Very low-depth samples are often failed libraries, and they should be identified before setting a cutoff.

### Step 2: Rarefaction curves (visual overview)

`vegan::rarecurve()` expects samples as rows.

```r
X <- t(seqtab)
rare <- rarecurve(X, step = 1000, sample = min(depth), label = FALSE)
```

Each curve shows observed ASVs against reads. Curves that flatten indicate saturation. This plot is useful for an overview, but the saturation depth is calculated in Step 3.

### Step 3: Slope-based saturation depth

`rareslope()` returns the slope of the rarefaction curve at a given depth. A slope near zero means the sample is saturated at that depth.

Because the function is evaluated for each sample at many depths, use a coarse grid (for example, every 1,000 reads) to keep runtime manageable.

```r
depth_grid <- seq(1000, 50000, by = 1000)

# Slope matrix: samples as rows, depth grid as columns
slopes <- t(sapply(seq_len(nrow(X)), function(i) {
  sapply(depth_grid, function(n) {
    if (n >= depth[i]) return(NA)          # depth beyond the sample's total
    rareslope(X[i, ], sample = n)
  })
}))
rownames(slopes) <- rownames(X)
colnames(slopes) <- depth_grid

# Numerical noise to zero
slopes[abs(slopes) < 1e-8] <- 0
```

Choose a slope threshold from the data, not in advance. A common starting point is 0.01 (about 10 new ASVs per 1,000 reads), but check the distribution first:

```r
hist(slopes[, which.min(abs(depth_grid - 10000))], breaks = 50,
     main = "Slope at 10,000 reads", xlab = "Slope")

# A few samples, including low- and high-richness ones
matplot(depth_grid, t(slopes[c(1, 100, 300), ]), type = "l", log = "x",
        xlab = "Reads", ylab = "Slope")
```

Compute the saturation depth for each sample: the smallest grid depth at which the slope falls below the threshold.

```r
threshold <- 0.01
sat_depth <- apply(slopes, 1, function(s) {
  i <- which(s < threshold)[1]
  if (is.na(i)) NA else depth_grid[i]
})
names(sat_depth) <- rownames(X)

sum(is.na(sat_depth))              # samples that never reached the threshold
summary(sat_depth)
```

Watch for two issues:

- **Non-monotone curves.** Some curves flatten, rise again, then flatten. Check a few samples visually before trusting `which()`.
- **Grid ceiling.** Samples with `NA` saturation may need a larger grid. Check whether any values equal the maximum grid depth, since these may not have saturated.

### Step 4: Objective 1 — per-sample adequacy

Compare each sample's depth with its saturation depth. A sample is adequate if its depth is at least its saturation depth.

```r
keep1 <- !is.na(sat_depth) & depth >= sat_depth
sum(!keep1)                         # samples dropped
names(depth)[!keep1]

plot(depth, sat_depth, log = "xy",
     xlab = "Sample depth (reads)", ylab = "Saturation depth (reads)")
abline(0, 1)                        # points below this line are adequate
```

Samples on or above the line are dropped. Samples with `NA` saturation are excluded by this rule, so resolve them first if you want to retain them.

### Step 5: Objective 2 — choosing a common cutoff

For each candidate cutoff, count samples that are both deep enough and saturated at that depth.

```r
cands <- seq(5000, 50000, by = 5000)
retained <- sapply(cands, function(cu) {
  sum(keep1 & depth >= cu & sat_depth <= cu)
})
data.frame(cutoff = cands, retained = retained)
```

Plot the cutoff with the depth-saturation data to see which samples are excluded:

```r
plot(depth, sat_depth, log = "xy",
     xlab = "Sample depth (reads)", ylab = "Saturation depth (reads)")
abline(0, 1)
abline(h = cutoff, col = "red")       # replace 'cutoff' with your chosen value
```

Choose the cutoff where the retained count levels off, or where the next step in cutoff would cost more samples than it gains in depth. Check the sensitivity of this choice to the slope threshold by repeating Steps 3–5 with thresholds such as 0.005 and 0.02. If the chosen cutoff changes substantially, report that the result depends on the threshold.

Then retain the samples that pass both objectives:

```r
cutoff <- 30000                        # example value; choose from Step 5
keep <- keep1 & depth >= cutoff & sat_depth <= cutoff
X_final <- X[keep, ]
```

### Step 6: Rarefy and compute diversity

Rarefying is random, so repeat it many times and average the results. `rrarefy()` returns a rarefied count matrix, with one random subsample per sample and each row summing to the cutoff.

```r
set.seed(123)
n_iter <- 100

# Alpha diversity: observed richness, averaged over iterations
richness <- replicate(n_iter, specnumber(rrarefy(X_final, sample = cutoff)))
mean_richness <- rowMeans(richness)

# Shannon diversity, averaged over iterations
shannon <- replicate(n_iter, diversity(rrarefy(X_final, sample = cutoff), index = "shannon"))
mean_shannon <- rowMeans(shannon)
```

For beta diversity, compute the distance for each iteration and average the distance matrices:

```r
dist_mats <- lapply(seq_len(n_iter), function(i) {
  as.matrix(vegdist(rrarefy(X_final, sample = cutoff), method = "bray"))
})
mean_dist <- as.dist(Reduce(`+`, dist_mats) / n_iter)
```

Use `mean_dist` for ordination (for example, `cmdscale()` or `capscale()`) and for permutation tests (for example, `adonis2()`), and `betadisper()` to check dispersion homogeneity.

For **differential abundance** analyses (Wilcoxon tests, MaAsLin3, ANCOM), do not use rarefied counts. Use the unrarefied table for the retained samples, with a method suited to count data.

### Step 7 (optional): Genus-level analysis

If you analyse diversity at the genus level, the sequencing depth is still a number of reads. What changes is the richness: genera saturate at lower depths than ASVs. Build the curves and saturation depths on the genus table, then choose the cutoff from those.

```r
# taxa: assignTaxonomy output, one row per ASV, with a Genus column
genus <- taxa[rownames(seqtab), "Genus"]
keep_asv <- !is.na(genus)              # decide how to treat unassigned ASVs

genus_counts <- rowsum(seqtab[keep_asv, ], group = genus[keep_asv])
Xg <- t(genus_counts)
depth_g <- rowSums(Xg)
rarecurve(Xg, step = 100, sample = min(depth_g), label = FALSE)
```

Repeat Steps 3–6 on `Xg`. Decide whether unassigned ASVs are dropped or kept as an "unassigned" category, and apply that decision consistently. Dropping them changes the depths, so check the cutoff afterwards.

---

## 4. Reporting checklist

When reporting a rarefaction-based analysis, state:

- The sequencing depth range before filtering (minimum, median, maximum).
- Which samples were removed, and why (controls, low-depth, failed Objective 1).
- The slope threshold and the grid used for saturation depth.
- The cutoff chosen and the number of samples retained.
- The number of rarefaction iterations and the random seed.
- That richness of under-sampled samples is an underestimate, if any retained samples fall close to the saturation boundary.
- Which analyses used rarefied data (alpha and beta diversity) and which used unrarefied data (differential abundance).
- Whether unassigned taxa were dropped, at which taxonomic level diversity was computed.

---

## 5. References and further reading

- Gotelli & Colwell (2001), Hurlbert (1971): the theory of rarefaction for species richness.
- McMurdie & Holmes (2014): the argument against rarefying for differential abundance, and when it is acceptable.
- Callahan et al. (2016): the DADA2 pipeline that produces the sequence table.
- Weiss et al. (2017): normalization choices for microbiome count data.
- vegan documentation: https://vegandevs.github.io/vegan/reference/
- iNEXT (sample-size-based rarefaction and extrapolation with confidence intervals): https://chao.shinyapps.io/iNEXTOnline/
- DADA2 documentation, including `mergeSequenceTables()`: https://benjjneb.github.io/dada2/
