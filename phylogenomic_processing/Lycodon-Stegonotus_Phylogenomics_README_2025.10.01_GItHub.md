# Lycodon and Stegonotus Phylogenomics README

This project will involve a phylogenomic assesment of the snake genera *Lycodon* and *Stegonotus* using Ultraconserved Elements (UCEs).

## Processing the Target Capture Data

First we will need to clear and trim the raw read data to remove barcode and adapter sequences using [illumiprocessor](https://github.com/faircloth-lab/illumiprocessor/). We will use the [Phyluce](https://phyluce.readthedocs.io/en/latest/installation.html) pipeline by Brant Faircloth to clean and processing the raw sequence data. Please refer to the Phyluce documentation for installation (this README will involve Phyluce being run on a computer cluster).

### Illumiprocessor: Adapter and Barcode Removal

The first part of illumiprocessor will be setting up a configuration (config) file (`illumiprocessor_lyco.conf`). This file is shown below for just two samples:

```bash
#illumiprocessor.conf for January 2025 - First round of Lycdon SqCL
# this is the section where you list the adapters you used.  the asterisk
# will be replaced with the appropriate index for the sample.
[adapters]
i7:AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC*ATCTCGTATGCCGTCTTCTGCTT
i5:GTGTAGATCTCGGTGGTCGCCGTATCATT*AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

# this is the list of indexes we used; you MUST include i7 or i5 at the beginning of the tags.
[tag sequences]
i7-RS04A_L0035-1:GTCAGTAC
i7-RS04A_L0036-1:AGATGTAC
i5-RS04A_L0035-2:TTCCATTG
i5-RS04A_L0036-2:GTCCTGTT

# this is how each index maps to each set of reads
[tag map]
RS04A_L0035:i7-RS04A_L0035-1,i5-RS04A_L0035-2
RS04A_L0036:i7-RS04A_L0036-1,i5-RS04A_L0036-2

# we want to rename our read files something a bit more nice - so we will
# rename Alligator_mississippiensis_GGAGCTATGG to alligator_mississippiensis
[names]
RS04A_L0035:Lyco_capucinus_BM8141_SG
RS04A_L0036:Lyco_aulicus_BM8156_SG

```

Note that you will need to know the adapters and barcode sequences for your samples to successfully remove these. The default regex for how illumiprocessor looks for your sequence data files is:

```
“{}_(?:.*)_(R1|READ1|Read1|read1)_\\d+.fastq(?:.gz)*"
“{}_(?:.*)_(R1|READ2|Read2|read2)_\\d+.fastq(?:.gz)*"
```

However, our sequence data is formatted as the following (how files are named varies by sequencing company):

```
RS04A_L0035_R1.fastq.gz
RS04A_L0035_R2.fastq.gz
RS04A_L0036_R1.fastq.gz
RS04A_L0036_R2.fastq.gz
```

Thus, we need to edit the regex formula that illumiprocessor uses. You can do this by using the `--r1-pattern` and `--r2-pattern` flags as follows:

```bash
illumiprocessor \
    --input /home/jbernstein/nas5/data_nas5/lyconotus/Lycodon_rawreads \
    --output /home/jbernstein/nas5/data_nas5/lyconotus/lycodon_phyluce/illumiprocessor/clean-fastq \
    --config /home/jbernstein/nas5/data_nas5/lyconotus/lycodon_phyluce/illumiprocessor/illumiprocessor_lyco2.conf \
    --cores 12 \
    --r1-pattern '{}_R1.fastq(?:.gz)*' \
    --r2-pattern '{}_R2.fastq(?:.gz)*'
```

### Check Quality of Sequences

We can check how our sequecnes did after illumiprocessor. Keep in mind, we need to use the **full** path

```bash
# move into illumiprocessor directory
cd /home/jbernstein/nas5/data_nas5/lyconotus/lycodon_phyluce/illumiprocessor/clean-fastq

# run this script against all directories of reads
for i in *;
do
    phyluce_assembly_get_fastq_lengths --input /home/jbernstein/nas5/data_nas5/lyconotus/lycodon_phyluce/illumiprocessor/clean-fastq/$i/split-adapter-quality-trimmed/ --csv; 
done
```

## Data Assembly: Spades

To assemble our contigs, we will use Spades in Phyluce. This process will be the most memory intensive and longest running part of the phyluce pipeline.  For this, we have to make another configuration file (only three samples shown):

```
[samples]
Lyco_capucinus_BM8141_SG:/home/jbernstein/nas5/data_nas5/lyconotus/lycodon_phyluce/illumiprocessor/clean-fastq/Lyco_capucinus/split-adapter-quality-trimmed/
Lyco_aulicus_BM8156_SG:/home/jbernstein/nas5/data_nas5/lyconotus/lycodon_phyluce/illumiprocessor/clean-fastq/Lyco_aulicus/split-adapter-quality-trimmed/
Lyco_zawi_CAS215570_MM:/home/jbernstein/nas5/data_nas5/lyconotus/lycodon_phyluce/illumiprocessor/clean-fastq/Lyco_zawi/split-adapter-quality-trimmed/
```

The paths need to be exact, as does the naming scheme you used. An important thing to note is that here, Spades' default maximum memory capacity is 8 GB. So even if, in the cluster or on a personal computer, you are giving more memory to your analysis, you need to utilize the `--memory` flag and allow Spades to ues the full amount of memory you want (I am giving is 244 GB since it is a 256 GB cluster node). Even if you specify 244 GB in your Slurm or PBS job submission script, Spades will still cap out at 8 GB of memory without the `--memory` flag:

To run Spades through Phyluce, we can run:

```bash
phyluce_assembly_assemblo_spades \
    --conf /home/jbernstein/nas5/data_nas5/lyconotus/lycodon_phyluce/spades/lycodon_spades1.conf \
    --output /home/jbernstein/nas5/data_nas5/lyconotus/lycodon_phyluce/spades/spades-assemblies1 \
    --cores 32 \
    --memory 244
```

### Assembly Quality Check

Similar to before, we can do another QC using the following:

```bash
for i in spades-assemblies/contigs/*.fasta;
do
    phyluce_assembly_get_fasta_lengths --input $i --csv;
done
```

**Note from the phyluce documentation:** The process of read assembly often differs by operating system and sometimes by OS version, and some of these differences are due to libraries that underlie many of the assembly programs. Expect to see differences. You should not expect for them to be huge.

**In addition:** If you see max-contig sizes around 16KB (for vertebrates), that is commonly the entire or almost-entire mtDNA genome. You do not tend to see entire mtDNA assemblies when the input DNA was extracted from a source having few mitochondria (e.g. blood).

## Finding UCE Loci

Next, we will need to actually find the targeted loci (UCEs). This involves separating contigs that are assocaited with UCEs from those that are not. We will need the probe fasta file for whichever prove set was used for your capture (whether that be UCEs 5k, UCEs 2.5k, SqCL, etc.):

```bash
wget https://raw.githubusercontent.com/faircloth-lab/uce-probe-sets/master/uce-5k-probe-set/uce-5k-probes.fasta
```

With this file in the same directory, run the command to match contigs to probes:

```bash
phyluce_assembly_match_contigs_to_probes \
    --contigs spades-assemblies/contigs \
    --probes uce-5k-probes.fasta \
    --output uce-search-results
```

Note: Make sure to run this on the post-illumiprocessor contigs. I accidentally used the files for the PNG_2020 dataset that were already run through *phyluce_assembly_match_contigs_to_probes* and the headers of the loci in the contig files had '`uce-XXXX`' in them, rather than '`NODE_`' (post-illumiprocessor format as the 'NODEs' have not been assigned to specific UCEs yet). If you made this mistake, and do not feel like going back and running illumiprocessor and spades on the old dataset again, you can add a 'uce:' line to the [headers] part of the phyluce config (found in whatever your version of `/home/jbernstein/nas4/anaconda3/envs/phyluce-1.7.3/phyluce/config/phyluce.conf` is):

```bash
[headers]
trinity:comp\d+_c\d+_seq\d+|c\d+_g\d+_i\d+|TR\d+\|c\d+_g\d+_i\d+|TRINITY_DN\d+_c\d+_g\d+_i\d+
velvet:node_\d+
abyss:node_\d+
idba:contig-\d+_\d+
spades:NODE_\d+_length_\d+_cov_\d+.\d+
uce:uce-\d+_\d+_length_\d+_cov_\d+\.\d+
```

###### A Note on Combining Data from Separate Sequencing Efforts

At this point, if you have done target capture/sequencing efforts before (other UCEs, SqCL, etc.), you can add the `*.contigs.fasta` files from those efforts into the `contigs` directory that `phyluce_assembly_match_contigs_to_probes` is searching in to perform the probe matching. As long as the sequencing efforts were targeting the same loci using the same probe set (whether that is partially or completely), then this should work fine. This was done for this project, but for the sake of keeping things simple, I will just be showing the *Lycodon* sequencing plate.

## Extracting UCE Loci

When extracting the UCE loci, we can (and need to) specify which samples we want to do this for, using the taxon names we used in prior parts of the workflow. We can create a another config file for this (subset of data shown):

```bash
[all]
Stegonotus__PNG26
Stegonotus_PNG27
Stegonotus_PNG28
Stegonotus_PNG29
```

We also need to make a directory for the outputs:

```bash
mkdir -p taxon-sets/all
```

Then we can run *phyluce_assembly_get_match_counts* to make the list of loci for each taxon*:*

```bash
# create the data matrix configuration file
phyluce_assembly_get_match_counts \
    --locus-db /home/jbernstein/nas5/data_nas5/lyconotus/lycodon_phyluce/spades/uce-search-results/probe.matches.sqlite \
    --taxon-list-config /home/jbernstein/nas5/data_nas5/lyconotus/lycodon_phyluce/spades/taxon-sets.conf \
    --taxon-group 'all' \
    --incomplete-matrix \
    --output /home/jbernstein/nas5/data_nas5/lyconotus/lycodon_phyluce/spades/taxon-sets/all/all-taxa-incomplete.conf
```

Because I used those files with different headers (`uce-` and `NODE_`), I was worried that there would be twice as many 'loci' if something got screwed up during the analysis (e.g., `uce-` and `NODE_` did end up with the same naming scheme and we essentially end up with double the loci).

You can run the following on the all-taxa-incomplete.conf file:

```bash
# print out all unique loci (phyluce prints each on a separate line) and then count the number of lines
sort all-taxa-incomplete.conf | uniq | wc -l

# Cool, we get 5197. Let's just confirm by specifically printing the duplicates
sort all-taxa-incomplete.conf | uniq -d 

# Cool again. We get no results. I.e., no duplicates!
```

Important note: This part of the phyluce pipeline seems to run into issues when there are hyphens/dashes in the names instead of underscores. You will get errors saying that {sample_name} was not found. If this happens, you should copy your contigs into a different directory I called this `contigs-combined-underscores` and rename them. For example, if we had the below samples, we would change them as follows:

```bash
mv Boig_dendrophila-BRN Boig_dendrophila_BRN
mv Boig_drapiezii-BRN Boig_drapiezii_BRN
```

Finally, we can get our FASTAs from our targeted UCEs!

```bash
# go into the /taxon-sets/all directory
cd taxon-sets/all

# make a log directory for the log file we are about to make
mkdir log

# get FASTA data for taxa in our taxon set
phyluce_assembly_get_fastas_from_match_counts \
    --contigs ../../spades-assemblies/contigs-combined-underscores \
    --locus-db ../../uce-search-results/probe.matches.sqlite \
    --match-count-output all-taxa-incomplete.conf \
    --output all-taxa-incomplete.fasta \
    --incomplete-matrix all-taxa-incomplete.incomplete \
    --log-path log
```

## Exploding the monolithic FASTA file

Next, we can separate our data by taxon/samples or by locus. This is called 'exploding' the monolithic FASTA file. Below, I explode it by taxon, but I actually do this by both in case I want alignment files for each locus. This step takes a few seconds to do, so its no waste of time and very beneficial to downstread analysis.

```bash
# explode the monolithic FASTA by taxon (you can also do by locus by commenting out '--by-taxon')
phyluce_assembly_explode_get_fastas_file \
    --input all-taxa-incomplete.fasta \
    --output exploded-fastas-taxon 
    --by-taxon                    # If you comment this line out, the monolith will be exploded by locus

# # get summary stats on the FASTAS
for i in exploded-fastas-taxon/*.fasta;
do
    phyluce_assembly_get_fasta_lengths --input $i --csv;
done
```

## Aligning the UCE loci

Now that we have our UCEs, we can align them. We will need to clean these up by trimming them, in which there are multiple options. You can align the UCEs and use the alignments with no trimming, you can do edge-trimming following an algorithm, or you can end+internally trim alignments following an algorithm. The phyluce tutorial suggests that when taxa are “closely” related (< 30-50 MYA, perhaps), edge-trimming alignments is reasonable. When the taxa span a wider range of divergence times (> 50 MYA), you may want to think about internal trimming. We will go with edge-trimming for now, but I will add the internal trimming steps below, as I checked these to compare the phylognies.

### Edge Trimming

```bash
# make sure we are in the correct directory
cd uce-tutorial/taxon-sets/all

# align the data
phyluce_align_seqcap_align \
    --input all-taxa-incomplete.fasta \
    --output mafft-nexus-edge-trimmed \
    --taxa 4 \
    --aligner mafft \
    --cores 32 \
    --incomplete-matrix \
    --log-path log
```

We can get some stats on our trimming as well:

```bash
# get edge trim stats
phyluce_align_get_align_summary_data \
    --alignments mafft-nexus-edge-trimmed \
    --cores 12 \
    --log-path log
```

### Internal Trimming

Just to check, we will also do some internal trimming as well. Start by aligning the data but do not trim them.

```bash
phyluce_align_seqcap_align \
    --input all-taxa-incomplete.fasta \
    --output mafft-nexus-internal-trimmed \
    --taxa 185 \
    --aligner mafft \
    --cores 32 \
    --incomplete-matrix \
    --output-format fasta \
    --no-trim \
    --log-path log
```

Now we perform internal trimming using [Gblocks](https://www.biologiaevolutiva.org/jcastresana/Gblocks.html#:~:text=Gblocks%20is%20a%20computer%20program,of%20DNA%20or%20protein%20sequences.):

```bash
# run gblocks trimming on the alignments
phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed \
--alignments mafft-nexus-internal-trimmed \
--output mafft-nexus-internal-trimmed-gblocks \
--cores 32 \
--log log
```

Get stats:

```bash
# get internal trim stats
phyluce_align_get_align_summary_data \
    --alignments mafft-nexus-internal-trimmed \
    --cores 12 \
    --log-path log
```

## Alignment Cleaning

One last thing we need to do before we run analyses on the alignments is to 'clean' them. This is beause the headers of our sequences contain both the locus and taxon names (we only want the taxon names). So let's clean those up:

```bash
# align the data - turn off trimming and output FASTA
phyluce_align_remove_locus_name_from_files \
    --alignments mafft-nexus-edge-trimmed \
    --output mafft-nexus-edge-trimmed-clean \
    --cores 12 \
    --log-path log
```

## Final Data Matrices

With our alignments now cleaned, we can make our final data matrices. This is done by 'completeness' in which, say, a 75% completeness matrix would denote that 75% percent of your taxa are in all alignments. We will make a 75% and 95% completeness matrix. It is important to note that these matrices are ambiguous to which taxa are in the alignments, so in a 75% completeness matrix, this does not mean that a given taxon will have data in all 75 alignments (if, say, out of 100 alignments).

```bash
# the integer following --taxa is the number of TOTAL taxa
# and I use "75p" to denote the 75% complete matrix. For a 95% completeness matrix, just change percent to 0.95 and change the output directory name
phyluce_align_get_only_loci_with_min_taxa \
    --alignments mafft-nexus-edge-trimmed-clean \
    --taxa 185 \
    --percent 0.75 \
    --output mafft-nexus-edge-trimmed-clean-75p \
    --cores 12 \
    --log-path log
```

## Preparing Data for Downstream Analysis

With our data matrices, we can now prepare these for an analysis like IQ-TREE2. Let's concatenate our alignments (shown below for the 75% completeness matrices):

```bash
# build the concatenated data matrix
phyluce_align_concatenate_alignments \
    --alignments mafft-nexus-edge-trimmed-clean-75p \
    --output mafft-nexus-edge-trimmed-clean-75p-raxml \
    --phylip \
    --log-path log
```

## Extracting UCEs from a Published Genome

For our divergence dating, we will want to get the split between Colubroidea and Psammodynastidae. We will download the genome sequence of *Psammodynastes pulverulentus.* We will use NCBI and use `wget` with the FTP link for the genome sequence.

```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/025/802/295/GCA_025802295.1_ASM2580229v1/GCA_025802295.1_ASM2580229v1_genomic.fna.gz
```

If we unzip the resulting file (`gunzip GCA_025802295.1_ASM2580229v1_genomic.fna.gz`) and look at the header of the sequence file, it is very text-filled, and has spaces:

```bash
gunzip GCA_025802295.1_ASM2580229v1_genomic.fna.gz 
head -n1 GCA_025802295.1_ASM2580229v1_genomic.fna 
```

Output:

```
>JAOYMU010000001.1 Psammodynastes pulverulentus voucher FMNH:273629 1, whole genome shotgun sequence
```

Let's convert this file to 2bit format using `faToTwoBit` from the [Kent Source Archive](http://hgdownload.soe.ucsc.edu/admin/exe/), which will remove everything following the first space in the header line:

```bash
/home/jbernstein/nas5/data_nas5/lyconotus/lycodon_phyluce/uce-genomes/bin/faToTwoBit GCA_025802295.1_ASM2580229v1_genomic.fna psamDyna_FMNH273629.2bit

/home/jbernstein/nas5/data_nas5/lyconotus/lycodon_phyluce/uce-genomes/bin/twoBitInfo psamDyna_FMNH273629.2bit sizes.tab
```

```bash
head -n 5 sizes.tab
```

Output:

```
JAOYMU010000001.1       159523122
JAOYMU010000002.1       59101030
JAOYMU010000003.1       21441
JAOYMU010000004.1       17570
JAOYMU010000005.1       18686
```

### Finding UCE Loci

First, get the probe file:

```bash
wget https://raw.githubusercontent.com/faircloth-lab/uce-probe-sets/master/uce-5k-probe-set/uce-5k-probes.fasta
```

Next, make a name for the database to create (here `tutorial3.sqlite`), the name of an output to store the lastz search results, the path to the genome sequences, the name of the genome sequences, the path to the probe file, and a number of compute cores to use:

```bash
# run the search
phyluce_probe_run_multiple_lastzs_sqlite \
    --db tutorial3.sqlite \
    --output psammodynastes-genome-lastz \
    --scaffoldlist psamDyna_FMNH273629 \
    --genome-base-path ./ \
    --probefile uce-5k-probes.fasta \
    --cores 12
```

### Extracting FASTA Sequence Matching UCE Loci from Genome Sequences

Once thes search is done, we can extract the identified loci from the genomes. First, we need to make a config file that states where each genome is found.

```
[scaffolds]
psamDyna_FMNH273629:/home/jbernstein/nas5/data_nas5/lyconotus/lycodon_phyluce/uce-genomes/psamDyna_FMNH273629/psamDyna_FMNH273629.2bit
```

Now, we can extract FASTA data from each genome for each UCE locus (so in this case, just the one sample). To do this, I input the path to the lastz files from above, the path to the conf file the `genome_lastz.conf` file, the amount of flanking sequence (to each side) that we want to slice, a name pattern, matching the lastz files that we would like to use, and the name of the output directory we want to create:

```bash
phyluce_probe_slice_sequence_from_genomes \
    --lastz psammodynastes-genome-lastz \
    --conf  genome_lastz.conf \
    --flank 500 \
    --name-pattern "uce-5k-probes.fasta_v_{}.lastz.clean" \
    --output psammodynastes-genome-fasta
```

One last note from the Phyluce tutorial: "The easiest way for you to use the extracted sequences is to basically pretend like they are “newly assembled contigs” and place the fasta files (`psamdyna_fmnh273629.fasta` as in the above) into a contigs directory. Alternatively, you can symlink them into a new or existing contigs folder (that resulted from a [PHYLUCE](https://github.com/faircloth-lab/phyluce) assembly process) and then proceed with the [Finding UCE loci](https://phyluce.readthedocs.io/en/latest/tutorials/tutorial-1.html#finding-uce-loci) procedure." Although we have already extracted the UCE loci from each genome sequence and even though it seems redundant to go back through the Finding UCE loci process, this is the best path forward.

## Phylogenetic Analysis (Gene Trees and Concatenated)

For the fresh samples, I removed some taxa by running the follow code: (note, this will remove taxa, but not change the number of taxa stated in the nexus file, so it is best to do this on FASTA alignments):

Taxa were removed due to data missingness causing long branches. I ran all gene trees in the computer cluster, as well as the concatenated analyses, with [IQ-TREE2 v2.2.20](http://www.iqtree.org/doc/):

```bash
iqtree2 -s {ALIGNMENT} -bb 1000 -nt AUTO
```

The above code uses ModelFinder2 to find the best substitution model and partitioning scheme. For the conatenated analyses, I used a GTR+G model:

```bash
iqtree2 -s {ALIGNMENT} -bb 1000 -m GTR+G -nt AUTO
```

## Species Tree Analysis

Using [ASTRAL-III](https://github.com/smirarab/ASTRAL), we can use the individual gene trees to make a species tree.

We need a file that has all of our gene trees listed, one on each line. We will call this `75p_fresh_iqtrees.tre`. The easiest way to obtain this is to run `cat *.treefile >> 75p_fresh_iqtrees.tre` in a directory with all of your finished IQ-TREE2 analyses. The file should look like this:

You will also need to down the astral `.jar` file. The `-o` flag will be your outfile name, and because we want to collapse individuals into single species tips, we will use a mapping file with the `-a` flag, Though, newer versions of ASTRAL can handle this, but for now we will use a mapping file called `mapping_file.txt`.

Finally, we can run ASTRAL!

```bash
java -jar astral.5.7.8.jar -i 75p_fresh_iqtrees.tre -a mapping_file.txt -o Homalopsidae_75p_ALLLoci_ATRAL_FreshOnly.tre
```

## Divergence Dating with treePL

To perform a divergence date estimation on our species tree (ASTRAL), we will use a phylogenetic penalized likelihood approach with [treePL](https://github.com/blackrim/treePL). The analysis is rather simple to execute and requires a Newick format of your species tree (you can convert it to Newick format in FigTree by exporting it) and a configuration (config) file. This config file will look like this:

```
treefile = E75p_ASTRAL.newick.tre
smooth = 100
numsites = 3249042 
mrca = COLUBROIDEA P_pu T_d
min = COLUBROIDEA 35.2 
max = COLUBROIDEA 46.76
mrca = COLU-NATR T_de Ly_f
min = COLU-NATR 39.9 
max = COLU-NATR 43.0
mrca = DAV-BAT Ly_da St_b
min = DAV-BAT 17.8 
max = DAV-BAT 17.8
outfile = stegonotus-lycodon_fresh_treePL.out.tre 
thorough 
# prime
opt = 2
moredetail
optad = 2
moredetailad
optcvad = 2
moredetailcvad
# cv 
# randomcv 
# cvoutfile = file 
# cvstart = number 
# cvstop = number 
# nthreads = integer 
#cviter = number 
#pliter = number 
#cvsimaniter = number 
#plsimaniter = number 
#log_pen 
#seed = number
```

You will need to provide the number of sites (numsites) that you made the ASTRAL tree based off of. Note that in this tree we are using the stem-Colubroidea fossil calibration from [Smith (2013)](https://doi.org/10.1016/j.jcz.2012.05.006) to set a minimum date of 35.2 million years ago, and used [Burbrink et al.&#39;s (2020)](https://doi.org/10.1093/sysbio/syz062) upper end of the 95% highest posterior density interval for the upper bound. I also set a secondary calibration using the 95% HPD from Burbrink et al. (2020) for the MRCA of Colubridae and Natricidae, and another secondary calibration from [Zaher et al.&#39;s (2019)](https://doi.org/10.1371/journal.pone.0216148) study to set a strict date for the *L. davisonii*-*S.batajanensis* node. This script was ran with the *prime* and *thorough* commands first, then I added parameters suggested by *prime*, and then I commented out *prime* and uncommonted the *randomcv* option. The output of that run will provide a `cv.out` file, in which the lowest error will give you the recommended smooth value (1000 in this case). Finally, after entering the smoothing value, make sure *cvout* and *prime* are commented out (leave the *thorough* command uncommented), and run the analysis to get your dated tree.

You can easily visualize this with FigTree, or make it look a bit fancier with some R code:

```R
# set working directory
setwd("C:/Users/Justin/Documents/Publications/Lycodon_Stegonotus/treePL")

# load libraries
library(ape)
library(strap)

# read time tree in
t <- read.tree("stegonotus-lycodon_fresh_treePL.out.tre")

# set the root time
t$root.time <- 44.6664

# plot the chronogram
geoscalePhylo(tree=ladderize(t,right=TRUE), units=c("Period", "Epoch"), boxes="Epoch", 
              cex.tip =0.8, cex.age=0.7, cex.ts=0.7, label.offset=0, x.lim=c(-15,45), lwd=3, width=2) 

# save pdf
pdf("LycoSteg_Chronogram.pdf")
geoscalePhylo(tree=ladderize(t,right=TRUE), units=c("Period", "Epoch"), boxes="Epoch", 
              cex.tip =0.8, cex.age=0.7, cex.ts=0.7, label.offset=0, x.lim=c(-15,45), lwd=3, width=2) 
dev.off()
```


