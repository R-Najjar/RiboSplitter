# RiboSplitter
RNA alternative splicing protein translation prediction

In order to gain meaningful insight into the biological consequences of alternative splicing events, we need to know the relative differences between the two isoforms of the alternative splicing event at the protein level, which exons are involved relative to the full transcript, and what protein domains are encoded by the involved exons. RiboSplitter was created to address these issues.

RiboSplitter will:
- Apply filtering to alternative splicing events and run differential testing
- Predict relative protein differences between the two isoforms of each event
- Create 3 figures: a zoomed-in view to show details and protein changes of the event, a zoomed-out view to show the location of the event relative to the full transcript, and a figure with protein domains aligned to exons. 

## Requirements
- R (4.2.1) libraries: tidyverse (2.0.0), Biostrings (2.66.0), biomaRt (2.54.1), patchwork (1.1.3), rhdf5 (2.42.1), aod (1.3.2)
- bedtools (2.30.0)


## Use
First run SplAdder (3.0.3)

Clone the repository, or load functions from raw url
```
source ([raw url ribosplitter.R])
```
Two datasets are required to start: 
- A dataset of details of alternative splicing events including gene name, event ID, genomic locations of involved exons. This will be called 'details'. These are SplAdder confirmed events. Provide the dir path to SplAdder output
```
details= read_details (dir)
```


- A dataset of read counts supporting isoforms 1 and 2 for each sample per event plus percent spliced-in (PSI). This will be called 'isoforms'. PSI is set to missing when isoform counts are <10. Provide the dir path to SplAdder output and prefix of the sample names
```
isoforms= read_isoforms (dir, 'sample_prefix')
```
Note 1: SplAdder's naming convention is to name the longer isoform as isoform 2. This works for most event types but I think it creates confusion for mutually exclusive events where either of the middle exons can be longer, or they can have the same length. I solve this by renaming them so that isoform 1 is always the first one in genomic location and by adjusting the PSI accordingly. 

Note 2: If using a different tool than SplAdder, you can still run RiboSplitter by creating the two datasets above from the tool's output. Include the following variables in the details dataset: gene_name, event_id, chrm, strand, type; start:stop genomic positions for involved exons (e1, e2, e3, e4); and for multiple exon skip events only, e2_starts contains colon separated start positions for all skipped exons, and e2_ends contains the end positions for multiple skipped exons. Format the isoforms dataset to include: event_id, sample, iso1, iso2, iso_total, and psi.

Next, merge your metadata with the isoforms dataset (adding a variable called "group" that classifies the samples into "disease" and "control")


## Running pipeline
Use ribosplitter function to run full pipeline: filtering, statistical testing, multiple comparison adjustment, reading frame prediction, protein translation, amino acid differences between isoforms, and generate event genomic names.

Input parameters
- "isoforms" and "details" dataframes. "isoforms" must contain a "group" variables with 2 values: "disease" and "control"
- min_disease: alternative splicing events with this minimum number of disease samples will be included
- min_control: alternative splicing events with this minimum number of control samples will be included
- sd_cutoff: The lowest PSI SD for an event to be included
- dir: path to save output files
- ref_fasta: path to reference genome fasta (needed to extract exon DNA sequences)
- q_cutoff: alpha value for q values (p values adjusted for multiple comparisons)

After filtering (events that are non-coding, have low variability, or have many missing values), RiboSplitter uses a beta binomial model with overdispersion to allow for variability by group, with p values adjusted for multiple comparisons using the BH methodology.

Not all events will be present in all samples. This could reflect real biology (e.g. disease heterogeneity), or incomplete/low sequencing coverage in some samples for specific splicing events. Therefore, we must decide how many samples per group are required for the event to be included. Additionally, low-variability events are excluded, based on PSI SD, I recommend using a low cutoff (e.g. sd_cutoff=0.05)


Example
```
de_events= ribosplitter (isoforms_df=isoforms, details_df=details,
                min_disease=10, min_control=10, sd_cutoff=0.05,
                dir='~/rna/sjogrenB',
                ref_fasta='~/rna/nw/GRCh38.primary_assembly.genome.fa',
                q_cutoff=0.05)

```



### Output
The output dataset will have the following elements:

| Variable(s) | Description |
| --- | --- |
| gene_name, gene_id, event_id | names for genes and splicing events |
| genomic_name | splicing event name based on chromosomal positions |
| novel | which isoforms are annotated vs not |
| n_control, n_disease | number of samples included for this event |
| avg_control, avg_disease | average percent spliced-in (PSI) by group |
| sd_control, sd_disease | standard deviation by group |
| pval, adj_pval | pvalues: raw and adjusted for multiple comparisons |
| frame, validframe | reading frames and their description (frame=9 means ambiguous reading frame) |
| iso1_dna, iso2_dna | DNA code for isoforms 1 and 2 |
| iso1_aa, iso2_aa | translated amino acid sequences for isoforms 1 and 2 |
| iso1_stop, iso2_stop | distance in amino acid letters from start of sequence to first stop codon |
| frameshift | relative frameshift introduced by the middle exons. 0 means no frameshift, 1 or 2 means a frameshift by 1 or 2 nucleotides |
| same_start | amino acid sequences that are the same in the two isoforms at the beginning |
| same_end | amino acid sequences that are the same in the two isoforms at the end |
| diff_iso1 | amino acid sequences that are different in isoform 1 (compared to 2) |
| diff_iso2 | amino acid sequences that are different in isoform 2 (compared to 1) |

A dataset called "positions" will be saved in the provided "dir". This contains event genomic positions and is needed to generate figures.

## Visualization

This function will create a zoomed-in view of alternative splicing events. Fig2 can be 'jitter' or 'dot'
```
f= splicing_figure(final_df, positions, isoforms, fig2='jitter')

# You can either save an individual figure or output multiple in one large figure as below
ggsave(paste0(out_dir,'/exons_figure.png'), plot=wrap_plots(f[1:30], ncol=3), 
       dpi=300,units='cm', height=45, width=60)
```
This figure will show the gene name, splicing event ID, chromosome, strand, and adjusted p value (q value). The x axis is genomic position. The upper left corner shows the relative frameshift between the two isoforms. The counts above and below the isoforms are reads supporting either isoform. Light blue exons mean they have the same amino acid sequences in both isoforms, while light green exons indicate relative differences in amino acid sequences. The adjacent figure shows PSI for all samples by group. Note that a PSI of 1 means there is 100% isoform 2 expression in that sample.
![fig1](https://github.com/R-Najjar/RiboSplitter/assets/119631106/daa371d2-e20f-43b8-83d1-9a21b903bb1d)




A 2nd figure can be generated as below that will zoom out and show the event parallel to its best fit transcript, so we can see the location of the splicing event. 
```
f= splice_figure_ref_TX (final_df , positions3)
ggsave(paste0(out_dir,'/exons_ref_transcript.png'), plot=wrap_plots(f[1:30], ncol=3),
           dpi=300,units='cm', height=40, width=60)
```
![fig2](https://github.com/R-Najjar/RiboSplitter/assets/119631106/65c74baa-ddbc-46bd-a219-e694918116f5)



A 3rd figure aligns exons (grey) to protein domains
```
domain_fig('ENSG00000064012.24')

# you can also create domain figures for all genes of interest as below:
domains_figs= domain_fig (final_df$gene_name)
```
![fig3](https://github.com/R-Najjar/RiboSplitter/assets/119631106/1dc529f0-2222-4af0-8664-f4106d659101)



