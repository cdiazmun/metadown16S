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
