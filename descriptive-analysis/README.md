# Descriptive analysis
The contents of this part of the analysis can be found in the [ASVsPlot](https://github.com/cdiazmun/ASVsPlot/tree/main) repository. I explained in detail every step needed to go from the two files that DADA2 outputs to various plots of relative abundance at genus and species level, and an extra ASV-level analysis in case we would be interested in following specific "strains". The data there consists of a much smaller number of samples, from a complex microbial community during a cocoa fermentation process. That's a very appropriate workflow if your data consists of temporal evolutions of a few samples. If your data consists of hundreds of independent samples, belonging to different factors, then a descriptive analysis as the one proposed there is maybe not that informative, or even feasible. An example of the latter scenario could be 16S bacterial community obtained from stool samples of patients and healthy controls. Nevertheless, the first steps will be needed for the subsequent downstream analyses, so I will detail down below an adapted version of the [ASVsPlot](https://github.com/cdiazmun/ASVsPlot/tree/main) initial steps to get relative abundances at different taxonomic levels. 

## Relative abundances at genus level

This part of the script will execute some commands to adapt the two tables described above (and that you should obtain after [fastq-processing](https://github.com/cdiazmun/metadown16S/tree/main/fastq-processing)) to a format that we can work with, transform the data, group the ASVs at genus level, and then perform any downstream analysis.

### Set up the basics of the script

First of all, let's set the working directory and the load the libraries needed for this script.

```r
setwd("C:/Path/to/your/working/directory")

library(reshape2) # Flexibly Reshape Data. To manipulate data. 
library(tidyverse) # A package of packages. To manipulate data, transform it to tidy-tables, and accomodate it to ggplot2.
library(ggthemes) # Extra Themes, Scales and Geoms for 'ggplot2'.
```

Next, we are going to define and load into the working space the three functions that will be used afterwards. The first one is to assign missing taxonomic assignments to the closest, upper rank that could be assigned for that ASV. The second one is to calculate the per-column percentage of relative abundance of every row and the third is to group very low occuring rows through the whole dataset under "Minorities". Why do I speak about rows and not ASVs? Because the functions do not care about what is in there, it will take whatever is in the rows, and this can be ASvs themselves, species or genera. And this is very important for the interpretation of the data, but you will see it more clearly when we actually use the functions.

``` r
# Function to fill taxa ranks that are unclassified
fill_rank <- function(x, parent) {
  coalesce(x, if_else(str_starts(parent, "Unclassified_"),
                      parent,
                      str_c("Unclassified_", parent)))
}

# Function to calculate the relative abundances
percentage <- function(d){
  d_per <- apply(d[2:ncol(d)], 2, function(x) x/sum(x)*100)
  d_per <- cbind(d[1],d_per)
  return(d_per)
}

# Function to put the genus/species that appear in less than "threshold" in others/minorities
minoritize_loose <- function(d,threshold){
  d_min <- d[!(apply(d[2:ncol(d)] < threshold, 1, function(x) all(x))),]
  minorities <- c("Minorities", apply(d_min[2:ncol(d_min)], 2, function(x) 100-sum(x)))
  d_min <- rbind(d_min, minorities)
  d_min[,2:ncol(d_min)] <- apply(d_min[,2:ncol(d_min)], 2, as.numeric)
  return(d_min)
}
```

### Import the input files and get the relative abundance table
Read the taxonomy table (`taxa`) and the sequence table (`reads`). In this case, since I ran DADA2 per batches (*i.e.,* per sequencing run), I have to incorporate them and join them iteratively:

```r
filereads <- list.files(pattern = "dada2_seqtab*")
filetaxa <- list.files(pattern = "dada2_taxa*")
run_id <- sapply(strsplit(basename(filereads), split = "_"),`[`, 3)
run_id <- sapply(strsplit(run_id, split = "[.]"),`[`, 1)

reads <- NULL
taxa <- NULL

for (i in 1:length(run_id)){
  df_r <- read_delim(filereads[i])
  df_t <- read_delim(filetaxa[i])
  if (is.null(reads)){
    reads<-df_r
    taxa<-df_t
  } else {
    reads <- reads |> 
      full_join(df_r)
    taxa <- taxa |> 
      full_join(df_t)
  }
}
```
Merge the two data frames using the first column (`Sequence`) as the common vector for the merging. 

```r
bacteria <- merge(taxa, reads, by="Sequence", all=T)
```

Next we are going to fill unclassified taxonomic ranks to their respective lowest classified rank. Notice that we also remove ASVs belonging to mithocondrial DNA (from human samples, which we are not interested, of course).

``` r
bacteria_filled <- bacteria |>
  filter_out(Family == "Mitochondria") |> # Remove the 16S belonging to mitochondria
  filter_out(is.na(Kingdom)) |> # Remove things we don't even know the Kingdom
  mutate(
    Phylum = fill_rank(Phylum,  Kingdom),
    Class  = fill_rank(Class,  Phylum),
    Order  = fill_rank(Order,  Class),
    Family = fill_rank(Family, Order),
    Genus  = fill_rank(Genus,  Family)
  )
```

Since we are only interested in the microbial composition at genus-level, we are going to select only the `Genus` column and discard the `Species` column. It is in general a **good practice** to create a new object every time we make a change to the data frame. This way if there is an error in the next steps or if we something weird in the data frame, we can go back step by step and see where the mistake was. Also, I **recommend** viewing the data frames as we transform them by clicking on them at the right-hand side, under Global Environment (in RStudio).

```r
bacteria_genus <- bacteria_filled[,c(7, 9:ncol(bacteria_filled))]
```

In the next steps we are going to perform some calculations, so we better transform any `NA` value to 0 to avoid problems. We do this step now, because there are not any NA anymore at the Genus column. 

```r
bacteria_genus[is.na(bacteria_genus)] <- 0
```

Transform the reads in percentage of relative abundance per sample and merge the ASV that correspond to the same genus.

```r
# Calculate percentage of relative abundance
bacteria_perc <- percentage(bacteria_genus)
# Aggregate all ASVs under the same genus by summing their respective ASV counts 
bacteria_aggr <- aggregate(bacteria_perc[,2:ncol(bacteria_perc)],
                            by=list(bacteria_perc[,"Genus"]), FUN=sum)
colnames(bacteria_aggr)[1] <- "Taxon"
```

As a sanity check, we make sure each genus label map to exactly one lineage:

```r
bacteria_aggr |> count(Taxon) |> filter(n > 1)
```
If everything's fine, then we can save the dataframe as such or transpose it. Depending on the type of downstream analysis we will need it in one way or the other.

```r
# Save the data frame
write.table(bacteria_aggr, "bacteria_genus_rel_abundance_filled.txt", 
            row.names = F, col.names = T, sep = "\t", quote = F)

# Transform to matrix and transpose it
bacteria_matrix <- as.matrix(bacteria_aggr[-1])
row.names(bacteria_matrix) <- bacteria_aggr$Taxon
tbac <- t(bacteria_matrix)

# Save the data frame
write.table(tbac, "bacteria_genus_rel_abundance_filled_transposed.txt", 
            row.names = T, col.names = T, sep = "\t", quote = F)
```

We can then use these data frames to perform a:
- Study the microbial community dynamics --> [ASVsPlot](https://github.com/cdiazmun/ASVsPlot/tree/main)
- Study the within-sample richness and evenness --> Alpha-diversity
- Study the differences in microbial composition across samples --> Beta-diversity 
- Study the differences in predicted metabolic potential of the microbiota across samples --> Functional analysis

### Calculate minorities

A further step we can perform is to group the genera that are in less than xx % of relative abundance in any sample. In other words, if there is only one sample in which a genus has more than xx %, it will be included as such, otherwise grouped under minorities. It is better to apply minoritize after Genus melting, so we can represent more accurately the Genera that are present in the samples. The threshold to set depends on your sample diversity and which story you want to tell. If you are interested in watching the big picture and the main microorganisms in your samples, you can increase that threshold. If you want to know all the small microbial groups that are present in, for example, some of your samples, then you can decrease the threshold to include more genera. 

```r
bacteria_min <- minoritize_loose(bacteria_aggr, 0.1)
# Inspect the resulting data frame
head(bacteria_min[1:10,1:10])
# Save the data frame
write_excel_csv2(bacteria_min, "bacteria_genus_relative_filled_minoritized.csv")
```
