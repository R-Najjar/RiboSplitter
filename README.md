# RiboSplitter
RNA alternative splicing protein translation prediction

In order to gain meaningful insight into the biological consequences of alternative splicing events, we need to know the relative differences between the two isoforms of the alternative splicing event at the protein level, which exons are involved relative to the full transcript, and what protein domains are encoded by the involved exons. RiboSplitter was created to address these issues.

RiboSplitter will:
- Apply filtering to alternative splicing events and run differential testing
- Predict relative protein differences between the two isoforms of each event
- Create 3 figures: a zoomed-in view to show details and protein changes of the event, a zoomed-out view to show the location of the event relative to the full transcript, and a figure with protein domains aligned to exons. 

## Requirements
- R libraries: tidyverse, Biostrings, biomaRt, patchwork, rhdf5, aod
- bedtools


## Use
First run SplAdder

Load functions
```
source ([raw url ribosplitter.R])
```
Two datasets are required to start: 
- A dataset of details of alternative splicing events including gene name, event ID, genomic locations of involved exons. This will be called 'details'. These are SplAdder confirmed events. Provide the dir path to SplAdder output
```
details= read_details (dir)
```


- A dataset of read counts supporting isoforms 1 and 2 for each sample per event plus percent spliced-in (PSI). This will be called 'isoforms'. PSI is set to missing when isoform counts are <10. Provide the dir path to SplAdder output 
```
isoforms= read_isoforms (dir)
```
Note 1: SplAdder's naming convention is to name the longer isoform as isoform 2. This works for most event types but I think it creates confusion for mutually exclusive events where either of the middle exons can be longer, or they can have the same length. I solve this by renaming them so that isoform 1 is always the first one in genomic location and by adjusting the PSI accordingly. 

Note 2: If using a different tool than SplAdder, you can still run RiboSplitter by creating the two datasets above from the tool's output. Include the following variables in the details dataset: gene_name, event_id, chrm, strand, type; start:stop genomic positions for involved exons (e1, e2, e3, e4); and for multiple exon skip events only, e2_starts contains colon separated start positions for all skipped exons, and e2_ends contains the end positions for multiple skipped exons. Format the isoforms dataset to include: event_id, sample, iso1, iso2, iso_total, and psi.

Next, merge your metadata with the isoforms dataset (adding a variable called 'group' that classifies the samples into disease and control)


### Filtering and statistical testing
After filtering (events that are non-coding, have low variability, or have many missing values), RiboSplitter uses a beta binomial model with overdispersion to allow for variability by group, with p values adjusted for multiple comparisons using the BH methodology.

Not all events will be present in all samples. This could reflect real biology (e.g. disease heterogeneity), or incomplete/low sequencing coverage in some samples for specific splicing events. Therefore, we must decide how many samples per group are required for the event to be included (n_dis and n_hc below). Additionally, a value for sd is needed, which is the lowest PSI SD for an event to be included, I recommend using a low value (e.g. 0.05)

```
events= isoforms %>%
    group_by(event_id) %>%
    summarise(sd= sd(psi), n_dis= sum (group=='disease'), n_hc= sum(group=='control') ) %>%
    filter (n_dis >=[number] & n_nt >=[number]) %>%
    filter (sd >=[number]) 
```

Exclude events in non-coding genes
```
g= gene_info (details$gene_name)
details2= details %>%
  filter (event_id %in% events) %>%
  left_join(g, by='gene_name') %>%
  mutate (type= str_extract(event_id,'[^.]+')) %>%
  filter (gene_biotype== 'protein_coding')
```

Create a new dataset 'positions' with each row being a single exon with its genomic positions. Then create a bed file to extract DNA code of each exon
```
positions=exon_positions (details2)

bed= positions %>%
  mutate (start=as.character (as.numeric(start)-1) ,
          score=300, id= paste0(chrm,':',start,"-",end,'(',strand,')')) %>%
  select (chrm, start, end, id, score, strand) %>%
  distinct ()

write_delim (bed, [path diff_events.bed'], col_names =F, delim='\t')
```

Use bedtools
```
bedtools getfasta -fi GRCh38.primary_assembly.genome.fa -bed ./diff_events.bed  \
  -fo ./diff_seq -s -tab
```

Read file output from bedtools
```
s= read.table(paste0(dir,'/diff_seq'))
colnames(s)= (c('bed_id','seq'))

# format and merge back to positions
positions2= positions %>%
  mutate (start=as.character (as.numeric(start)-1), 
          bed_id= paste0(chrm,':', start, "-",end,'(',strand,')'),
          start=as.character (as.numeric(start)+1)) %>%
  left_join(s,by='bed_id') 
```

### Reading frame prediction, part I
```
# Identify events with a single open reading frame in first 5 prime exon
exon1frames= first_exon (positions2)

# This excludes noncoding transcripts
possibleframe= filter (exon1frames, validframe !='no protein')$event_id

```

### Statistical testing
```
samples= filter (isoforms, event_id %in% possibleframe ) %>%
  mutate (grp= as.factor (group))

# This will run beta binomial model, adjust for multiple hypotheses testing, and create event-level dataset (each row is a splicing event)
diff_events= event_level (samples, details2)
```

### Reading frame prediction, part II
Continue reading frame prediction, this time by comparing translated peptides to known gene products
```
t= filter (exon1frames, event_id %in% allsle2$event_id) 
frames= peptide_match (t, details2)

table (frames[,3:4])
table (frames[,c(2,5)])
```

Next, run function to put exons on one row per event to create isoform 1 and isoform 2 in DNA sequence 
then translate them to amino acid sequences 
```
t= filter (positions2, event_id %in% allsle2$event_id)
iso_dna= stitch_exons (t, frames)
```

### Protein changes prediction
This function will extract amino acid sequences that are the same at the beginning and end of the peptides from the two isoforms, and sequences that are different between the two isoforms
```
protein_diff= protein_changes (iso_dna)
```

Create names for splicing events based on genomic coordinates (chromosome, strand, start and stop positions of exons)
```
positions3= filter (positions2, event_id %in% diff_events$event_id)
names= splice_name(positions3)
```


Finally, put all together
```
final_df= left_join(diff_events, names,by='event_id') %>%
  relocate (genomic_name, .after=event_id) %>%
  left_join(protein_diff, by=c('gene_name', 'event_id'))
```
Note: save the positions and isoforms datasets as these will be needed for visualizations

### Output
The output dataset 'final_df' will have the following elements:

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


## Visualization

This function will create a zoomed-in view of alternative splicing events. Fig2 can be 'jitter' or 'dot'
```
f= splicing_figure(final_df, positions3, isoforms, fig2='jitter')

# You can either save an individual figure or output multiple in one large figure as below
ggsave(paste0(out_dir,'/exons_figure.png'), plot=wrap_plots(f[1:30], ncol=3), 
       dpi=300,units='cm', height=45, width=60)
```
This figure will show the gene name, splicing event ID, chromosome, strand, and adjusted p value (q value). The x axis is genomic position. The upper left corner shows the relative frameshift between the two isoforms. The counts above and below the isoforms are reads supporting either isoform. Light blue exons mean they have the same amino acid sequences in both isoforms, while light green exons indicate relative differences in amino acid sequences. The adjacent figure shows PSI for all samples by group. Note that a PSI of 1 means there is 100% isoform 2 expression in that sample.
![fig1](https://github.com/R-Najjar/RiboSplitter/assets/119631106/9880ba30-9748-4fb9-b877-c598e9e603e0)



A 2nd figure can be generated as below that will zoom out and show the event with its best fit transcript, so we can see the location of the splicing event. 
```
f= splice_figure_ref_TX (final_df , positions3)
ggsave(paste0(out_dir,'/exons_ref_transcript.png'), plot=wrap_plots(f[1:30], ncol=3),
           dpi=300,units='cm', height=40, width=60)
```
![fig2](https://github.com/R-Najjar/RiboSplitter/assets/119631106/405637aa-fa1e-46ef-9643-ff28ff59e4a4)


A 3rd figure aligns exons (grey) to protein domains
```
domain_fig('ENSG00000064012.24')

# you can also create domain figures for all genes of interest as below:
domains_figs= domain_fig (final_df$gene_name)
```
![fig3](https://github.com/R-Najjar/RiboSplitter/assets/119631106/370cf247-695b-4490-992c-4d87e1af6f87)

