# Rarefaction analysis

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
- Use **unrarefied** counts for **differential abundance** testing.

Note that rarefaction is a statistical equalization of sampling effort, not a way to recover information that was not sequenced. Richness estimates for under-sampled samples are underestimates, and this should be stated when reporting results.

### Key concepts

| Term | Meaning |
|---|---|
| ASV | Amplicon sequence variant, produced by DADA2 |
| Sequencing depth | Total number of reads in a sample |
| Rarefaction curve | Observed richness versus number of reads subsampled |
| Saturation depth | The depth at which a sample's rarefaction curve has flattened |
| Rarefaction slope | Rate at which new ASVs appear per additional read at a given depth |
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

*Criterion:* A sample is adequate if it is saturated within its own reads, that is, its saturation depth is lower than its total depth. On a plot of sequencing depth (x) against saturation depth (y), adequate samples fall below the 1:1 line.

*Action:* Drop samples on the line or too close. This decision does not depend on any common cutoff.

*Applies to:* Any downstream analysis. 

### Objective 2: A common cutoff

*Question:* At what depth can every retained sample be rarefied so that alpha and beta diversity are computed on comparable, adequately sampled data?

*Criterion:* For a candidate cutoff *c*, a sample is usable if (a) its total depth is at least *c*, so that it can be subsampled to *c*, and (b) its saturation depth is at most *c*, so that the subsample captures its richness.

*Action:* Count retained samples across candidate cutoffs, choose a cutoff where the count stops improving or where the loss of samples is no longer worth the gain in depth, then rarefy.

*Applies to:* Diversity analyses. 

Objective 1 is applied first, and objective 2 is applied to the samples that pass it.

---

## 3. Step-by-step guide

Rarefaction is typically performed using functions embedded in the `vegan` package. The input is the DADA2 sequence table with **samples as columns and ASVs as rows**.

### Step 0: Load packages and prepare the table

```r
library(vegan)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
setwd("/path/to/your/working/directory/")

# Import the metadata files 
meta <- read_delim("metadata/metadata_file.txt")
# Clean the metadata file (if needed) to get only true samples
meta_clean <- meta |> 
  filter_out(SampleID %in% c("NTC", "MOCK", "PrepCo"))
# Import the sequence table
filereads <- list.files(path = "dada2/output/", pattern = "dada2_seqtab*", full.names = T)
run_id <- sapply(strsplit(basename(filereads), split = "_"),`[`, 3)
run_id <- sapply(strsplit(run_id, split = "[.]"),`[`, 1)

seqtab <- NULL

for (i in 1:length(run_id)){
  df_r <- read_delim(filereads[i])
  if (is.null(seqtab)){
    seqtab<-df_r
  } else {
    seqtab <- seqtab |> 
      full_join(df_r)
  }
}
```
Full joins across runs introduce NAs for ASVs absent from a run. Absent ASVs must be zero counts, not missing values.
``` r
# Transform NAs into 0s
seqtab[is.na(seqtab)] <- 0
```
We need to work with a numerical matrix, so let's remove the first column (ASVs)
```r
seqtab_rows <- seqtab[-1]
row.names(seqtab_rows) <- seqtab$Sequence
```
Remove mock communities and other controls before computing depths, since they distort the cutoff.
```r
# Example
no_samples <- meta$`NGS-Code`[meta$SampleID %in% c("MOCK", "NTC", "PrepCo")]

seqtab_final <- seqtab_rows |> 
  select(!(all_of(no_samples)))
```

### Step 1: Compute per-sample depth

Get the per-sample read depth, which is basically the column sum, because samples are columns.

```r
depth <- colSums(seqtab)
summary(depth)
sort(depth)
```

Inspect the spread. Very low-depth samples are often failed libraries, and they should be identified before setting a cutoff.

### Step 2: Rarefaction curves (visual overview)

`vegan::rarecurve()` expects samples as rows, so we need to transpose the matrix/data frame.

```r
tseqtab <- as.data.frame(t(seqtab_final))
colnames(tseqtab) <- rownames(seqtab_final)
# Plot the rarefaction curves
rare <- rarecurve(tseqtab, step = 1000, sample = min(depth), label = FALSE,
                  xlab = "Number of reads", ylab = "ASV Richness", col = "blue")
```

Each curve shows observed ASVs against reads. Curves that flatten indicate saturation. This plot is useful for an overview, but the saturation depth is calculated in Step 3.

### Step 3: Slope-based saturation depth

`rareslope()` returns the slope of the rarefaction curve at a given depth. A slope near zero means the sample is saturated at that depth.

Because the function is evaluated for each sample at many depths, you can use a coarse grid (for example, every 1,000 reads) to keep runtime manageable. I didn't because I'm loco, so:

```r
slope <- rareslope(tseqtab, sample = depth - 1)

# Clean the slope matrix
depths <- as.numeric(sub("^N", "", colnames(slope)))
ord <- order(depths)
slope_ordered <- slope[, ord]
depths_ordered <- depths[ord]
```
If a sample has a depth of 50,000, it doesn't make much sense to calculate the slope for N = 100,000. In our matrix it's 0, but it should be `NA`. We can use this to test if everything went well:
```r
# Depth beyond the sample's total: undefined
slope_ordered[slope_ordered == 0] <- NA

# Check the NA pattern against the sample depths
last_valid <- apply(slope_ordered, 1, function(s) max(depths_ordered[!is.na(s)]))
head(cbind(last_valid, depth))
```

Choose a slope threshold from the data, not in advance. A common starting point is 0.01 (about 10 new ASVs per 1,000 reads), but check the distribution first:

```r
# Change the number of reads to see (5000, 10000, 20000, ...)
col_idx <- which.min(abs(depths_ordered - 10000))
hist(slope_ordered[, col_idx], breaks = 50, main = "Slope at 10,000 reads", xlab = "Slope")
```

Check how many samples reach saturation (slope < 0.01) with that read cutoff. Now we can compute the saturation depth for each sample: the smallest grid depth at which the slope falls below the threshold.

```r
threshold <- 0.01
sat_depth <- apply(slope_ordered, 1, function(s) {
  i <- which(s < threshold)[1]
  if (is.na(i)) NA else depths_ordered[i]
})
summary(sat_depth)

# Inspect the sat_depth by depth
head(cbind(sat_depth, depth))
# Put it on a data frame
sat_df <- as.data.frame(cbind(sat_depth, depth)) |> 
  mutate(sampleID = rownames(sat_df),
         enough_depth = ifelse(sat_depth <= depth, TRUE, FALSE))
```

### Step 4: Objective 1 — per-sample adequacy

Compare each sample's depth with its saturation depth. A sample is adequate if its depth is greater than its saturation depth.

```r
plot(sat_df$depth, sat_df$sat_depth, log = "xy",
     xlab = "Sequencing depth (ASVs)", ylab = "Saturation depth (ASVs)")
abline(0,1)                        # points below this line are adequate
```

Samples on the line or very close to it should be dropped.

### Step 5: Objective 2 — choosing a common cutoff

For each candidate cutoff, count samples that are both deep enough and saturated at that depth. If your cutoff is *e.g.,* 10,000 reads, then you'll only retain samples that fulfill two criteria: sequencing depth &ge; 10,000 AND saturation depth &le; 10,000. In our case, we are selecting 20,000 reads, which is a very good compromise to retain most of the samples. It's actually 19,693, which is the closest depth to 20,000. 

```r
sat_df <- sat_df |> 
  mutate(sat_20k = ifelse(sat_depth <= 19693 & depth >= 19693, "pass", "fail"))
```

Plot the cutoff with the depth-saturation data to see which samples are excluded:

```r
ggplot(sat_df, aes(x = depth/1000, y = sat_depth/1000, color = sat_20k)) + 
  geom_point() + 
  geom_hline(yintercept = 19.693, linetype = "dashed", color = "darkred") +
  geom_vline(xintercept = 19.693, linetype = "dashed", color = "darkred") +
  geom_abline(slope = 1) +
  scale_y_log10() +
  scale_x_log10() + 
  scale_color_manual(values = c("darkgray", "darkred")) +
  theme_bw() +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.1, 0.9),
    panel.grid = element_blank(),
    axis.text = element_text(color = "black")
  ) + 
  labs(x="Sequencing depth (ASVs, in thousands)",
       y="Saturation depth (ASVs, in thousands)")

ggsave("rarefaction/seqdepth_vs_satdepth_threshold.png", width = 8.5, height = 6, units = "in", dpi = 300)
```

Choose the cutoff where the retained count levels off, or where the next step in cutoff would cost more samples than it gains in depth. Check the sensitivity of this choice to the slope threshold by repeating Steps 3–5 with thresholds such as 0.005 and 0.02. If the chosen cutoff changes substantially, report that the result depends on the threshold.

Then retain the samples that pass both objectives:

```r
cutoff <- 20000
tseqtab_kept <- tseqtab[depth >= cutoff, ]
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
