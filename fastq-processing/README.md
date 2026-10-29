# Processing of fastq reads
## A bit of background
Probably a colleague, a collaborator, or yourself performed the microbial DNA extraction from a matrix (*e.g.* plant material, stool samples, soil), which was subsequently purified and subjected to PCR amplification of a barcoding region of interest. This region will typically be the v3v4 region of the 16S rRNA gene for bacteria or the ITS1-ITS2 region of the rRNA locus for fungi (commonly yeasts are studied within this kingdom). These regions are hypervariable, meaning that phylogenetic differences among different taxonomic groups are clearly visible here, which make them ideal to characterize complex microbial communities. Amplicon-based sequencing of these regions has been employed for over a decade already, making use of short-read Illumina sequencing technologies (*e.g.*, [MiSeq](https://www.illumina.com/systems/sequencing-platforms/miseq.html), [NovaSeq](https://www.illumina.com/systems/sequencing-platforms/novaseq.html)). However, long-read sequencing technologies from Oxford Nanopore Technologies (*e.g.,* [MinION](https://nanoporetech.com/products/sequence/minion)) and PacBio (*e.g.,* [HiFi sequencing](https://www.pacb.com/technology/hifi-sequencing/)) have been applied to the metabarcoding field to increase the taxonomic resolution (typically limited to genus-level) and accuracy. This pipeline will tackle the analysis of 16S amplicon sequencing using short-read Illumina sequencing, but the steps and principles of the analysis remain the same if you're studying yeasts, or if you're using a different sequencing technology. 

## Quality check
The initial part of the pipeline involves the quality check and trimming of the fastq reads. **THE** software to analyze the overall quality of your reads is [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). It is GUI-based, so it is very easy to manually inspect each sample independently. However, if you have hundreds of samples, this becomes impractical. You can download it and install it in your machine and run it from the command-line to process many samples at once. The output typically consists of an `html` that you can inspect, again, sample-by-sample. So, this step is mostly focused on detecting if there is something wrong in the sequencing run, which is rarely the case. And if this happens, the sequencing facility (if you're using one) would have noticed before sending the fastq files to you. But still, it's recommended to have a look to a few forward (R1) and reverse (R2) reads as a sanity-check. Actually, by looking at the fastQC files, you'll probably be able to see if the adaptors have been removed or they are still there. 

## PCR primers removal
The most common scenario is that the Illumina barcodes used to multiplex during sequencing are already removed from the fastq samples you received but the PCR primers are still there. So, next step is to remove the PCR primers used during the DNA amplification. We can use [Trimmomatic](https://github.com/usadellab/trimmomatic) or [Cutadapt](https://cutadapt.readthedocs.io/en/stable/) for that. There will probably be more options, but these are the two widely used tools for adapter removal and trimming (if needed) of sequencing reads. For now, we will only focus on primer removal, since quality filtering and trimming is done in a subsequent step of the pipeline.

So, let's install `Cutadapt`:
``` bash
conda config --add channels bioconda
conda create -n cutadapt cutadapt
conda activate cutadapt
```
This tool has a good manual so just do `--help` or click on the link I left above to know more about the options. I created a short `bash-script` named `cutadapt.sh` to process several samples at once. I ran it per batches, being each sequencing run a batch of samples to process iteratively. Please note that the `FWD` and `REV` primers can change, depending on the PCR conditions used, so you will need to adapt the nucleotide sequence. 

``` bash
#!/usr/bin bash

# 341F
FWD=CCTACGGGNGGCWGCAG

# 785R
REV=GACTACHVGGGTATCTAATCC

# Set the PATH
INPATH=/path/to/your/working/directory/fastq_reads
OUTPATH=/path/to/your/working/directory/cutadapt

cd ${INPATH}

for r1 in *SeqRunID*_R1_001.fastq.gz; do
    r2=${r1/_R1/_R2}

    sample=$(basename "$r1")
    sample=${sample%%_R1*}

    echo "Processing ${sample}"

    cutadapt \
        -g "${FWD}" \
        -G "${REV}" \
        -a "${REV}" \
        -A "${FWD}" \
        --discard-untrimmed \
        --cores 0 \
        --report=minimal \
        -o "${OUTPATH}/${sample}_R1_001.trimmed.fastq.gz" \
        -p "${OUTPATH}/${sample}_R2_001.trimmed.fastq.gz" \
        "$r1" "$r2"
done
```
Then, I used the standard output generated to save it in a way I can easily track the percentage of reads kept after trimming:

``` bash
# Edit the script for the Run ID of interest
nano cutadapt.sh
# Execute the script
bash cutadapt.sh
# Copy the standard output from the terminal and paste inside the newly created log file
nano SeqRunID_A.log
# Take the lines containing the relevant info and paste them on a new file
grep OK SeqRunID_A.log > SeqRunID_A_processed.log
# Do the same in case there are lines with a WARNING instead than an OK
grep WARN SeqRunID_A.log >> SeqRunID_A_processed.log
```

## Amplicon sequence variant (ASV) inference
Long-story short, one would process amplicon sequencing reads to ensure a minimum quality and length, then merge forward and reverse reads and cluster the final, merged reads into operational taxonomic units (OTUs). OTUs would cluster very similar sequences into a single taxon (at a level of family, genus, species depending on the phylogeny), based on nucleotide identity percentage (threshold of 97%). This was later replaced by ASVs, which icorporates a model that learns from the (Illumina) sequencing errors on your samples to differentiate them from true genetic variation. Thisw way, exact sequences can be infered, referred to as amplicon sequence variants (ASVs), resolving differences as little as 1 nucleotide. This pipeline was developed by [Benjamin Callahan and colleagues](https://www.nature.com/articles/nmeth.3869) and implemented in the widely used [DADA2](https://benjjneb.github.io/dada2/) pipeline. That's the gold standard protocol for amplicon sequencing processing until the current day. The website has very nice tutorials for [bacterial 16S](https://benjjneb.github.io/dada2/tutorial.html), [fungal ITS](https://benjjneb.github.io/dada2/ITS_workflow.html), and [long-read PacBio full-length 16S](https://github.com/benjjneb/LRASManuscript) sequencing data. Since the tutorials have all the information needed, I won't detail any step, since they will basically repeat what it is already there. For any further question, error, or concern I recommend visiting their [GitHub Issues](https://github.com/benjjneb/dada2/issues), there will be someone that had that problem before you, for sure. 

Just To have something to compare with, these are the parameters I used and the tracking of reads throughout the pipeline:
- **Filter and Trim**. truncLen: 270; maxN:0; maxEE:2; truncQ:2.
- **Learn Error Rates**. *122,641,020 total bases in 454,226 reads from 2 samples will be used for learning the error rates*.
- **Sample Inference**. *1,990 sequence variants were inferred from 24,275 input unique sequences*.
- **Merge paired reads**. *378,981 paired-reads (in 11,427 unique pairings) successfully merged out of 409,055 (in 25,532 pairings) input*. 
- **Construct sequence table**. *93,082 ASVs in 33 samples*.
- **Remove chimeras**. *Identified 81688 bimeras out of 93082 input sequences. 11,394 ASVs in 33 samples*.

You should always aim to have more than 50% of your total input reads retained by the end of the pipeline. If you removed too many during **Filter and Trim**, relax the parameters; if you loose a lot during **Merge paired reads** check that you did not trim too much (ITS are tricky here and you may even need to proceed just with FW reads); if you loose a lot during **Remove chimeras** then check tha primer removal or...maybe the PCR was not so good on the first place. These were my numbers for >700 samples (which included MOCK communites, PrepCo, NTC) from around 10 different sequencing batches. Notice the huge standard deviations... this means that there is a lot of variability between sequencing batches and between samples, so it's key to keep the same conditions for all the samples that we will later on process together, so that they are comparable. 

|Step    |Percentage of reads lost|
|--------|------------------------|
|Filter and trim    |13.44 &plusmn; 10.1|
|Merge paired reads    |15.61 &plusmn; 6.8|
|Remove chimeras    |26.25 &plusmn; 9.6|

## Taxonomy assignment
Still within the DADA2 pipeline we will assign the taxonomy to the ASVs we just inferred. To do that, DADA2 uses a bootstrapping approach (not pure similarity as OTUs) to assign a specific taxonomic level (ideally spesies, most often genus, and sometimes just family level or above). To do that, we make use of an already formatted database of bacteria, namely the Silva database, which is available at their [website](https://benjjneb.github.io/dada2/training.html) for download and use. For yeasts, I have always used the [UNITE database](https://unite.ut.ee/repository.php), although other exists. Another possibility is to use the Genome Taxonomy Database (GTDB), also available already formatted for DADA2 at their website. As DADA2 developers mention themselves: *The dada2 package also implements a method to make species level assignments based on exact matching between ASVs and sequenced reference strains. Recent analysis suggests that exact matching (or 100% identity) is the only appropriate way to assign species to 16S gene fragments.* Still, species-level assignmente using short reads must be interpreted with caution. If you really want to go to specis level safely, the use long-read sequencing of the full-length 16S or do metagenomics instead. 

So, by the end of the pipeline you'll be able to export two files:
- **The sequence table**: First column is the ASV sequence, rest of the columns are the counts of that ASV in each sample. Each row is an ASV
- **The taxonomy table**: First column is the ASV sequence, rest of the columns are the taxonomic classification from Kingdom to Genus/Species. Each row is an ASV


