
# RiboSplitter
# Written by Rayan Najjar

suppressPackageStartupMessages({
  library (tidyverse)
  library (patchwork)
  library (rhdf5)
  library (Biostrings)
  library (aod)
})

# read files. spladder confirmed events with exon start and stop positions
readspladder= function (event_name) {
  df= read.table (paste0('merge_graphs_',event_name,'.confirmed.txt.gz'), header=T) %>%
    select (gene_name,event_id,chrm,strand,is_annotated, starts_with('e'))
  return (df)
}
read_details= function (dir) {
  setwd(dir)
  t1= readspladder('exon_skip_C3')
  t2= readspladder('alt_3prime_C3')
  t3= readspladder('alt_5prime_C3')
  t4= readspladder('intron_retention_C3')
  t5= readspladder('mutex_exons_C3')
  t6= readspladder('mult_exon_skip_C3')
  det= bind_rows (t1,t2,t3,t4,t5,t6) 
  details= det %>%
    mutate (type= str_extract(event_id,'[^.]+')) %>%
    ungroup() %>%
    mutate (e1= paste(e1_start,e1_end,sep=':'), e2= paste(e2_start,e2_end,sep=':'),
            e3= paste(e3_start,e3_end,sep=':'), e4= paste(e4_start,e4_end,sep=':')) %>%
    mutate (condition= ifelse (type=='mutex_exons' & e3_start<e2_start,1,0) ,
            e2= ifelse (condition==1, paste(e3_start,e3_end,sep=':'), e2), 
            e3= ifelse (condition==1, paste(e2_start,e2_end,sep=':'), e3) ) %>%
    dplyr::rename (flipped_mxe=condition)
}

# isoforms 
read_isoforms= function (dir, prefix) {
  setwd(dir)
  hf1= H5Fopen ('merge_graphs_exon_skip_C3.counts.hdf5')
  hf2= H5Fopen ('merge_graphs_alt_3prime_C3.counts.hdf5')
  hf3= H5Fopen ('merge_graphs_alt_5prime_C3.counts.hdf5')
  hf4= H5Fopen ('merge_graphs_intron_retention_C3.counts.hdf5')
  hf5= H5Fopen ('merge_graphs_mult_exon_skip_C3.counts.hdf5')
  hf6= H5Fopen ('merge_graphs_mutex_exons_C3.counts.hdf5')
  iso= function (hf,event_type) {
    conf=hf$confirmed
    confirmed= as.data.frame (conf) %>%
      mutate (event_id=paste0(event_type,'.',row_number())) %>%
      filter (conf >0)
    iso1= hf$iso1 
    colnames(iso1)= hf$samples  
    iso1= as.data.frame (iso1) %>%
      mutate (event_id=paste0(event_type,'.',row_number())) %>%
      pivot_longer(cols=starts_with(prefix),names_to ='sample',values_to ='iso1')
    iso2= hf$iso2 
    colnames(iso2)= hf$samples  
    iso2= as.data.frame (iso2) %>%
      mutate (event_id=paste0(event_type,'.',row_number())) %>%
      pivot_longer(cols=starts_with(prefix),names_to ='sample',values_to ='iso2')
    psi= hf$psi 
    colnames(psi)= hf$samples  
    psi= as.data.frame (psi) %>%
      mutate (event_id=paste0(event_type,'.',row_number())) %>%
      pivot_longer(cols=starts_with(prefix),names_to ='sample',values_to ='psi')
    fin= full_join(iso1,iso2,by=c('event_id','sample')) %>%
      full_join(psi,by=c('event_id','sample')) %>%
      filter (event_id %in% confirmed$event_id )
    return (fin)
  }
  t1=iso (hf1,'exon_skip')
  t2=iso (hf2,'alt_3prime')
  t3=iso (hf3,'alt_5prime')
  t4=iso (hf4,'intron_retention')
  t5=iso (hf5,'mult_exon_skip')
  t6=iso (hf6,'mutex_exons')
  iso= rbind (t1,t2,t3,t4,t5,t6)
  h5closeAll()
  t= select (details, event_id,flipped_mxe)
  isoforms= left_join(iso, t, by='event_id') %>%
    mutate (old1=iso1, old2=iso2, oldpsi=psi, flipped_mxe=replace_na(flipped_mxe,0),
            iso1=ifelse (flipped_mxe==1,old2,iso1),
            iso2=ifelse (flipped_mxe==1,old1,iso2),
            psi= ifelse (flipped_mxe==1,1-oldpsi,psi)) %>%
    select(-c(old1,old2,oldpsi))
  return (isoforms)
}

# wrapper function for pipeline
# isoforms must contain a "group" variable with 2 values only: "disease" and "control"
ribosplitter= function (isoforms_df, details_df, min_disease, min_control, sd_cutoff=0.05,
                        dir, ref_fasta, q_cutoff, ensembl_version) {
  suppressWarnings({
    events= isoforms_df %>%
      group_by(event_id) %>%
      summarise(sd= sd(psi), n_dis= sum (group=='disease'), n_hc= sum(group=='control') ) %>%
      filter (n_dis >=min_disease & n_hc >=min_control) %>%
      filter (sd >=sd_cutoff) 
    
    g= gene_info (details_df$gene_name, ensembl_version)
    details2= details_df %>%
      filter (event_id %in% events$event_id) %>%
      left_join(g, by='gene_name') %>%
      mutate (type= str_extract(event_id,'[^.]+')) %>%
      filter (gene_biotype== 'protein_coding')
    
    positions=exon_positions (details2)
    
    bed= positions %>%
      mutate (start=as.character (as.numeric(start)-1) ,
              score=300, id= paste0(chrm,':',start,"-",end,'(',strand,')')) %>%
      select (chrm, start, end, id, score, strand) %>%
      distinct ()
    
    write_delim (bed, paste0(dir,'/exons.bed'), col_names =F, delim='\t')
    
    system2 (command='bedtools', 
             args=c('getfasta', '-fi', ref_fasta, 
                    '-bed', paste0(dir,'/exons.bed'),
                    '-fo', paste0(dir,'/exons_seq'), 
                    '-s', '-tab'))
    
    s= read.table(paste0(dir,'/exons_seq'))
    colnames(s)= (c('bed_id','seq'))
    
    positions2= positions %>%
      mutate (start=as.character (as.numeric(start)-1), 
              bed_id= paste0(chrm,':', start, "-",end,'(',strand,')'),
              start=as.character (as.numeric(start)+1)) %>%
      left_join(s,by='bed_id') 
    
    exon1frames= first_exon (positions2, ensembl_version)
    possibleframe= filter (exon1frames, !validframe %in% c('no protein','stops+nearbyTSS'))$event_id
    
    diff_samples= filter (isoforms, event_id %in% possibleframe) %>%
      mutate (grp= as.factor (group))
    
    all_pvals= event_level (diff_samples, details2)
    
    diff2= filter (all_pvals, adj_pval<q_cutoff)
    
    events2= diff2$event_id
    t= filter (exon1frames, event_id %in% events2) 
    frames= peptide_match (t, details2, ensembl_version)
    
    t= filter (positions2, event_id %in% events2)
    iso_dna= stitch_exons (t, frames)
    protein_diff= protein_changes (iso_dna)
    
    positions3= filter (positions2, event_id %in% events2)
    samples3= filter (isoforms, event_id %in% events2)
    
    event_names= splice_name(positions3)
    
    ## final data 
    diff3= left_join(diff2, event_names,by='event_id') %>%
      relocate (genomic_name, .after=event_id) %>%
      left_join(protein_diff, by=c('gene_name', 'event_id')) %>%
      mutate (frameshift= ifelse (type=='mutex_exons' & e2shift != e3shift, abs(e2shift-e3shift),
                                  ifelse (type !='mutex_exons' & e2shift>0,e2shift,0)),
              delta_avg= abs(avg_disease-avg_control)) 
    
  })
  saveRDS(positions3, paste0(dir,'/positions.rds'))
  saveRDS(diff3, paste0(dir,'/differential_events.rds'))
  file.remove(paste0(dir,'/exons.bed'))
  file.remove(paste0(dir,'/exons_seq'))
  
  return (diff3)
}

# input gene names in ensembl gene ID version format, 
# function will classify them into protein coding vs not
gene_info= function (gene_ensembl_id_version, ensembl_version) {
  library (biomaRt)
  # human genes
  ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version=ensembl_version)
  gg= getBM (mart=ensembl, filters='ensembl_gene_id_version', unique (gene_ensembl_id_version), 
             attributes =c('ensembl_gene_id_version' ,'ensembl_gene_id','external_gene_name', 'gene_biotype'))
  gg= dplyr::rename (gg, gene_name=ensembl_gene_id_version, ensembl=ensembl_gene_id, 
                     gene_id=external_gene_name) %>%
    mutate (gene_id= ifelse (gene_id=='',ensembl,gene_id))
  remove(ensembl)
  detach("package:biomaRt", unload = TRUE)
  return(gg)
}

# function that takes details file and puts one exon per row with genomic positions
exon_positions= function (details_df) {
  pos1= details_df %>%
    ungroup() %>%
    select (-c(type, ensembl,flipped_mxe,e1_start, e1_end, e2_start, e2_end, e3_start, e3_end, e4_start, e4_end, 
               e2_starts, e2_ends,gene_biotype,is_annotated,gene_id)) %>%
    pivot_longer(cols=c(e1,e2,e3,e4), names_to ='exons', values_to ='position') %>%
    filter (position !='NA:NA') 
  mult.skip= details_df %>%
    filter (str_starts(event_id,'mult_exon_skip')==T) %>%
    ungroup() %>%
    select (event_id, e2_starts, e2_ends) 
  begin= str_split_fixed (mult.skip$e2_starts,':',Inf)
  finish= str_split_fixed (mult.skip$e2_ends,':',Inf)
  df=data.frame (event_id=character(), position=character())
  for (r in 1:nrow(begin)) {
    id= mult.skip$event_id[r]
    for (c in 1: ncol(begin)) {
      exon= paste(begin[r,c], finish[r,c],sep=':')
      if (exon==':') {break} 
      df= add_row(df, event_id=id, position=exon)
    }
  }
  pos2= left_join(df,details_df,by='event_id') %>%
    mutate (exons='e2') %>%
    select(gene_name, event_id, chrm, strand, exons, position )
  positions= rbind (pos1,pos2) %>%
    mutate (start=str_extract(position,'[^:]+'), end=str_extract(position,'(?<=:)[:digit:]+')) %>%
    arrange (gene_name, event_id, exons, start)
  return (positions)
}

## input data with id variables and dna seq. function will translate in 3 frames and count stop codons
readframes= function (dnaseq, id) {
  df= data.frame(event_id=character(),dna=character(), r1=character(),r2=character(),r3=character(),
                 s1=numeric(),s2=numeric(),s3=numeric(),frame=numeric())
  for (i in 1:length(dnaseq)) {
    event_id=id[i]
    dna= dnaseq[i]
    t= DNAString(dna)
    r1= toString (translate (t), no.init.codon=T)
    r2= toString (translate (subseq (t,2)), no.init.codon=T)
    r3= toString (translate (subseq (t,3)), no.init.codon=T)
    r1= ifelse(r1=='','*',r1)
    r2= ifelse(r2=='','*',r2)
    r3= ifelse(r3=='','*',r3)
    s1= str_count (r1,"\\*")
    s2= str_count (r2,"\\*")
    s3= str_count (r3,"\\*")
    frame=ifelse(s1==0 & min(s2,s3)>0, 1,
                 ifelse(s2==0 & min(s1,s3)>0, 2,
                        ifelse(s3==0 & min(s1,s2)>0, 3, 9)))
    df= add_row(df,event_id,dna,r1,r2,r3,s1,s2,s3,frame)
  }
  return (df)
}

# identify the correct reading frame based on first exon (in 5' to 3' direction) with a single open reading frame
first_exon = function (exon_positions, ensembl_version) {
  library (Biostrings)
  neg= exon_positions %>%
    filter (strand=='-') %>%
    group_by(event_id) %>%
    arrange (event_id, exons, start) %>%
    slice_tail(n=1)
  pos= exon_positions %>%
    filter (strand=='+') %>%
    group_by(event_id) %>%
    arrange (event_id, exons, start) %>%
    slice_head(n=1)
  fiveprime_e1= rbind (pos,neg) %>%
    mutate (len= str_length(seq)) %>%
    filter (len >2)
  frames= readframes(fiveprime_e1$seq, fiveprime_e1$event_id)
  exonframes= frames %>%
    mutate (validframe=ifelse(frame<9,'singleframe', 
                              ifelse (s1>0 & s2>0 & s3>0,'no protein',
                                      ifelse (s1+s2+s3==0,'all3',
                                              ifelse(s1+s2==0,'1and2',
                                                     ifelse(s1+s3==0,'1and3',
                                                            ifelse(s2+s3==0,'2and3','other')))))))
  fiveprime= dplyr::select(fiveprime_e1, gene_name, event_id, chrm, strand, start, end)
  noprotein= filter (exonframes, validframe=='no protein') %>%
    left_join(fiveprime, by='event_id') %>%
    dplyr::select (gene_name, event_id,chrm, strand, start, end)
  library (biomaRt)
  ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version=ensembl_version)
  tx= getBM (mart=ensembl, filters='ensembl_gene_id_version', noprotein$gene_name, 
             attributes =c('ensembl_gene_id_version' , 'ensembl_transcript_id', 
                           'transcription_start_site', 'transcript_biotype') ) %>%
    filter (transcript_biotype=='protein_coding')  %>%
    dplyr::select (ensembl_gene_id_version, transcription_start_site) %>%
    distinct () %>%
    dplyr::rename (gene_name=ensembl_gene_id_version, tss=transcription_start_site)
  tss= inner_join(noprotein, tx, by='gene_name') %>%
    mutate (start= as.numeric(start), end=as.numeric(end), 
            
            tss_pos_distance= ifelse (strand=='+', 
                                      ifelse (tss < start, paste0('upstream',start-tss),
                                              ifelse (tss > end,  paste0('downstream',tss-end),
                                                      'within0' )) , 
                                      ifelse (tss > end, paste0('upstream',tss-end),
                                              ifelse (tss < start, paste0('downstream',start-tss), 
                                                      'within0' ))        # neg strand 
            ),
            tss_pos= str_extract(tss_pos_distance, '[:alpha:]+'),
            tss_distance= as.numeric (str_extract(tss_pos_distance, '[:digit:]+')) ,
            nearby_upstream= ifelse (tss_pos=='upstream' & tss_distance <=100,1,0), 
            inside= ifelse ( tss >= start & tss <= end, 1,0 ) ) %>%
    group_by(gene_name,event_id) %>%
    summarise(upstream= sum(tss_pos=='upstream'), downstream=sum(tss_pos=='downstream'), 
              inside=sum(tss_pos=='within'), nearby_up=sum(nearby_upstream==1), tss_dist=min(tss_distance) ) %>%
    mutate (close2tss= ifelse (inside+downstream+nearby_up >0, 1, 0)) %>%
    ungroup () %>%
    filter (close2tss==1) %>%
    dplyr::select (event_id, close2tss)
  exonframes2= left_join(exonframes,tss,by='event_id') %>%
    mutate (close2tss=replace_na(close2tss,0),
            validframe= ifelse (close2tss==1,'stops+nearbyTSS',validframe)) %>%
    dplyr::select(-close2tss)
  remove(ensembl)
  detach("package:biomaRt", unload = TRUE)
  return (exonframes2)
}

# beta binomial model with adjusted p value. input samples with isoform counts for 2 groups
betabinomial= function (samples_df) {
  t= droplevels(samples_df)
  p= data.frame (event_id=character(), beta=numeric() )
  for (i in unique (t$event_id)) {
    df=filter (t,event_id==i) %>%
      drop_na(psi) 
    beta= betabin (cbind(iso2, iso1) ~ grp, ~ grp, data=df, link='logit')
    p= add_row(p,event_id=i, beta= summary(beta)@Coef[2,4])
  }
  pval= p %>%
    mutate (adj.beta= p.adjust(beta, method='BH'),
            sig= ifelse (adj.beta <0.05,1,0),
            sig= replace_na(sig,0)) %>%
    dplyr::rename (pval=beta, adj_pval=adj.beta) 
  return (pval)
}

# create event level file. avg psi per group. add event details
event_level= function (samples_df, details_df) {
  pvalue= betabinomial(samples_df)
  avg= filter (samples_df, event_id %in% pvalue$event_id) %>%
    group_by(event_id, grp) %>%
    summarise(avg= mean (psi), sd=sd(psi), n=n()) %>%
    ungroup() %>%
    pivot_wider(names_from =grp, values_from = c(avg,sd,n))
  df= avg %>%
    inner_join(pvalue, by='event_id') %>%
    left_join(details_df, by='event_id') %>%
    mutate (novel= ifelse(is_annotated==0,'both isoforms are novel',
                          ifelse(is_annotated==1,'isoform 2 is novel',
                                 ifelse(is_annotated==2,'isoform 1 is novel',
                                        ifelse(is_annotated==3,'both isoforms are not novel',' '))))) %>%
    select (gene_name,gene_id,event_id,type, novel, starts_with('n_'), starts_with('avg'), starts_with('sd_'), 
            pval, adj_pval, chrm,strand, is_annotated) %>%
    arrange (adj_pval )
  return (df)
}

# function to do exact matching of first exon amino acid sequence with known proteins of the gene
# input is data from first_exon function above and details file
peptide_match = function (e1frames, details_df, ensembl_version) {
  library (biomaRt)
  ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version=ensembl_version)
  frames_aa= e1frames %>%
    dplyr::select (event_id,r1,r2,r3) %>%
    pivot_longer(cols=c(r1,r2,r3), names_to = 'frame', values_to ='aa') %>%
    mutate (stop_codon=str_count (aa,"\\*"), len=str_length(aa)) %>%
    filter (stop_codon==0 & len >1)
  t= dplyr::select(details_df, gene_name, event_id)
  frames_aa= left_join(frames_aa,t, by='event_id')
  genes=unique (frames_aa$gene_name)
  pep= getSequence(mart=ensembl, type='ensembl_gene_id_version', id=genes, 
                   seqType ='peptide')
  pep= filter(pep, peptide !='Sequence unavailable') %>%
    dplyr::rename (gene_name=ensembl_gene_id_version) %>%
    arrange (gene_name)
  mdf= inner_join(frames_aa,pep,by='gene_name',relationship = "many-to-many")
  df= data.frame(event_id=character(), frame=character(), matched=numeric())
  for (i in 1:length(mdf$aa)) {
    event_id=mdf$event_id[i]
    frame= mdf$frame[i]
    a= mdf$aa[i]
    matched= countPattern(pattern=AAString(a), subject=AAString(mdf$peptide[i]), min.mismatch =0,
                          max.mismatch =0, with.indels =F, algorithm ='auto')
    df= add_row(df, event_id,frame, matched)
  }
  pepframe= df %>%
    group_by(event_id, frame) %>%
    summarise(n=sum(matched==1)) %>%
    pivot_wider(names_from = 'frame', values_from ='n', values_fill =0) %>%
    mutate (peptideframe= ifelse (r1>0 & max(r2,r3)==0,1,
                                  ifelse (r2>0 & max(r1,r3)==0,2,
                                          ifelse(r3>0 & max(r1,r2)==0,3,9)))) 
  pepframe= dplyr::select (pepframe, event_id, peptideframe)
  frames_df= left_join(e1frames,pepframe,by='event_id') %>%
    dplyr::rename (exonframe=frame) %>%
    mutate (peptideframe= replace_na(peptideframe,9)) %>%
    dplyr::select (c(event_id,validframe, exonframe, peptideframe)) %>%
    mutate (frame= ifelse(exonframe==9 & peptideframe !=9, peptideframe, exonframe),
            validframe= ifelse (frame<9,'singleframe',validframe))
  remove(ensembl)
  detach("package:biomaRt", unload = TRUE)
  return (frames_df)
}

# function to stitch exons together to form isoforms
stitch_exons = function (positions_df, frames_df) {
  skip= filter (positions_df, str_starts (event_id,'mult_exon_skip')==T & exons=='e2') %>%
    arrange (gene_name, event_id, start) %>%
    dplyr::select (gene_name, event_id,exons, start,seq)
  df= data.frame (event_id=character(), e2_mskip=character())
  for (i in unique (skip$event_id)) {
    t= filter (skip, event_id==i)
    exon=character()
    for (k in t$seq) {
      exon= paste0(exon,k)
    }
    df= add_row(df, event_id=i, e2_mskip=exon)
  }
  t= positions_df %>%
    anti_join(skip, by=c('event_id','exons')) %>%
    dplyr::select (gene_name,strand, event_id, exons, seq) %>%
    group_by (gene_name, event_id) %>%
    pivot_wider(names_from =exons, values_from = seq) %>%
    left_join(frames_df,by='event_id')
  iso_dna= left_join(t,df,by='event_id') %>%
    mutate (type= str_extract(event_id,'[^.]+'),
            e2len=ifelse (type=='mult_exon_skip', nchar(e2_mskip,type='chars'), nchar(e2,type='chars')),
            e3len=ifelse (type=='mutex_exons', nchar(e3,type='chars'), NA ),
            e2shift= e2len%%3, e3shift= e3len%%3,
            
            iso1= ifelse (strand=='-',  
                          ifelse (type=='mutex_exons' , paste0(e4,e3,e1), paste0(e3,e1)), # for - strand
                          ifelse (type=='mutex_exons' , paste0(e1,e2,e4), paste0(e1,e3)) ), # for + strand
            
            iso2= ifelse (strand=='-', 
                          ifelse (type=='mutex_exons', paste0(e4,e2,e1), 
                                  ifelse (type=='mult_exon_skip', paste0(e3,e2_mskip,e1), 
                                          paste0(e3,e2,e1)) ),
                          ifelse (type=='mutex_exons', paste0(e1,e3,e4),                   # for + strand
                                  ifelse (type=='mult_exon_skip', paste0(e1,e2_mskip,e3), 
                                          paste0(e1,e2,e3)) )) ) %>%
    dplyr::select (gene_name,event_id, type, frame, validframe, iso1, iso2, e2len, e3len, e2shift, e3shift)
  t= filter (iso_dna, frame==9 & validframe !='stops+nearbyTSS') %>%
    ungroup () %>%
    dplyr::select (event_id, iso1, iso2) %>%
    pivot_longer(!event_id, names_to ='iso', values_to ='seq') %>%
    mutate (id= ifelse (iso=='iso1',paste0(event_id,'_1'),paste0(event_id,'_2')))
  updatedframe= readframes(t$seq, t$id) %>%
    filter (frame !=9) %>%
    select (event_id,frame) %>%
    mutate (event_id=str_sub(event_id,1,-3)) %>%
    group_by(event_id) %>%
    summarise(mn=min(frame),mx=max(frame)) %>%
    filter (mn==mx) %>%
    dplyr::rename (updatedframe=mn) %>%
    select (event_id, updatedframe)
  iso_dna2= left_join(iso_dna,updatedframe,by='event_id') %>%
    ungroup () %>%
    mutate (updatedframe=replace_na(updatedframe,9) , frame= ifelse (frame==9 & updatedframe<9,updatedframe,frame),
            validframe= ifelse (frame<9,'singleframe',validframe))
  t= filter (iso_dna2, frame !=9) 
  df= data.frame(event_id=character(), iso1_aa=character(),iso2_aa=character(),
                 iso1_stop=numeric(),iso2_stop=numeric())
  for (i in 1:length(t$iso1)) {
    event_id= t$event_id[i]
    f= t$frame[i]
    iso1_aa= toString (translate (subseq (DNAString(t$iso1[i]), f ),if.fuzzy.codon='solve'), no.init.codon=T)
    iso2_aa= toString (translate (subseq (DNAString(t$iso2[i]), f ),if.fuzzy.codon='solve'), no.init.codon=T)
    iso1_stop= str_locate (iso1_aa,"\\*")[1]
    iso2_stop= str_locate (iso2_aa,"\\*")[1]
    df= add_row(df, event_id,iso1_aa,iso2_aa,iso1_stop,iso2_stop)
  }
  iso_dna3= left_join(iso_dna2, df, by='event_id') %>%
    dplyr::rename (iso1_dna=iso1, iso2_dna=iso2, full_isoform_frame=updatedframe) 
  return (iso_dna3)
}

# predict protein changes between 2 isoforms
# input data with iso1_aa and iso2_aa that have amino acid sequences
protein_changes= function (translated_isoforms) {
  tt= filter (translated_isoforms, frame !=9) %>%
    ungroup () %>%
    dplyr::select (event_id, iso1_aa, iso2_aa) 
  df= data.frame (event_id=character(), first_mismatched=numeric(),last_mismatched=numeric(), 
                  same_start=character(), same_end=character(), diff_iso1=character(), diff_iso2=character())
  for (k in 1:length (tt$event_id)) {
    event_id= tt$event_id[k]
    str1= tt$iso1_aa[k]
    str2= tt$iso2_aa[k]
    a= as.character (str_split_fixed (str1, '', Inf))
    b= as.character (str_split_fixed (str2, '', Inf))
    a_len= length(a)
    b_len= length(b)
    mn= min (a_len, b_len)
    for (i in 1:mn) {
      if (a[i] != b[i] ) {
        first_mismatched=i
        break
      } else if (mn==i) {
        first_mismatched=i+1
        break
      }
    }
    for (i in 0:mn) {
      if (mn==i) {
        last_mismatched=i
        break
      } else if (a[a_len-i] != b[b_len-i] ) {
        last_mismatched=i
        break
      } 
    }
    same_start= str_sub (str1,1,first_mismatched-1)
    same_end= ifelse (last_mismatched==0,' ', str_sub (str1,-last_mismatched)) # zero is when last AA is different
    diff_iso1 = str_sub (str1,first_mismatched,a_len-last_mismatched)
    diff_iso2 = str_sub (str2,first_mismatched,b_len-last_mismatched)
    df= add_row (df, event_id, first_mismatched, last_mismatched, same_start, same_end,
                 diff_iso1, diff_iso2) 
    remove (first_mismatched, last_mismatched)
  }
  df2= left_join(translated_isoforms, df, by='event_id') %>%
    select (-type)
  return (df2)
}

## function to create event names based on genomic positions. input positions
splice_name= function (positions_df) {
  pos= positions_df %>%
    group_by (event_id) %>%
    arrange (event_id, exons, start) %>%
    select (event_id, chrm,strand,position)
  df= data.frame (event_id=character(), exons=character())
  for (i in unique (pos$event_id)) {
    t= filter (pos, event_id==i)
    exon=character()
    for (k in t$position) {
      exon= paste(exon,k, sep=',')
    }
    df= add_row(df, event_id=i, exons=exon)
  }
  df=mutate (df,exons=str_sub (exons,2)) 
  pos2= pos %>%
    select (-position) %>%
    distinct () %>%
    left_join(df, by='event_id') %>%
    mutate (type= str_extract(event_id,'[^.]+'), genomic_name= paste0(chrm,strand,exons)) %>%
    select (event_id, genomic_name)
  return (pos2)
}

# find genomic position
find_genomic = function (df) {
  new_df=data.frame (id=character(), genomic_position=numeric())
  neg= filter (df, strand=='-') %>%
    arrange (id, desc(exons), desc(end))
  for (i in unique(neg$id)) {
    genomic_position=numeric()
    t= filter (neg,id==i) 
    id=i
    location= t$distance[1]
    for (k in 1:length (t$id)) {
      len=length(t$id)
      exon_len= abs(t$end[k]-t$start[k])+1
      if  (exon_len-location >=0) {
        genomic_position= t$end[k]+1-location
        break
      } else { location=location-exon_len }
    }
    new_df= add_row (new_df, id,  genomic_position)
  }
  pos= filter (df, strand=='+') %>%
    arrange (id, exons,start)
  for (i in unique(pos$id)) {
    genomic_position=numeric()
    t= filter (pos,id==i) 
    id=i
    location= t$distance[1]
    for (k in 1:length (t$id)) {
      len=length(t$id)
      exon_len= t$end[k]-t$start[k]+1 # genomic positions include diff+1 
      if (exon_len-location >=0) {
        genomic_position=t$start[k]-1+location
        break
      } else { location=location-exon_len }
    }
    new_df= add_row (new_df, id,  genomic_position)
  }
  return (new_df)
}

find_last_mismatch = function (df) {
  new_df=data.frame (id=character(), genomic_position=numeric())
  neg= filter (df, strand=='-') %>%
    arrange (id, exons, end)
  for (i in unique(neg$id)) {
    genomic_position=numeric()
    t= filter (neg,id==i) 
    id=i
    location= t$distance[1]
    for (k in 1:length (t$id)) {
      exon_len= abs(t$end[k]-t$start[k])+1
      if (location==0) {
        genomic_position=t$start[k]
        break
      } else if (exon_len-location >=0) {
        if (k==1 & t$frameshift[k]==0) {genomic_position=t$end[k]+1
        } else {genomic_position= t$start[k]+location}
        break
      } else { location=location-exon_len }
    }
    new_df= add_row (new_df, id,  genomic_position)
  }
  pos= filter (df, strand=='+') %>%
    arrange (id, desc(exons), desc(start))
  for (i in unique(pos$id)) {
    genomic_position=numeric()
    t= filter (pos,id==i) 
    id=i
    location= t$distance[1]
    for (k in 1:length (t$id)) {
      exon_len= t$end[k]-t$start[k]+1 
      if (location==0) {
        genomic_position=t$end[k]
        break
      } else if (exon_len-location >=0) {
        if (k==1 & t$frameshift[k]==0) {genomic_position=t$start[k]-1
        } else {genomic_position=t$end[k]-location}
        break
      } else { location=location-exon_len }
    }
    new_df= add_row (new_df, id,  genomic_position)
  }
  return (new_df)
}

# function to create zoomed-in illustration of alternative splicing event with annotated protein changes
# input positions data, final events data and samples df
# output list with figures with first stop codon (red line), differences in proteins between isoforms (light green)
# whether there was a frameshift relative between isoforms
# choose whether 2nd figure is a jitter or dot plot
splicing_figure= function (events_df, positions_df, samples_df, fig2) {
  t= positions_df %>%
    filter (event_id %in% events_df$event_id) %>%
    arrange (gene_name, event_id, exons, start) %>%
    mutate (type= str_extract(event_id,'[^.]+')) %>%
    dplyr::select (gene_name, event_id, type, chrm, strand, exons, start, end) 
  iso2= t %>%
    mutate (incl=ifelse (type=='mutex_exons' & ((strand=='+' & exons !='e2') | (strand=='-' & exons !='e3')),1,0),
            iso='Isoform 2') %>%
    filter (type !='mutex_exons' | incl==1) %>%
    select (-incl)
  iso1= t %>%
    mutate (incl=ifelse (type=='mutex_exons' & ((strand=='+' & exons !='e3') | (strand=='-' & exons !='e2')),1,0), 
            iso='Isoform 1') %>%
    filter ((type !='mutex_exons' & exons !='e2') | incl==1) %>%
    select (-incl)
  t= events_df %>%
    ungroup () %>%
    dplyr::select(gene_name, gene_id) %>%
    distinct ()
  both= rbind(iso1, iso2) %>%
    mutate (start=as.numeric(start), end=as.numeric(end)) %>%
    left_join(t, by='gene_name') %>%
    group_by(gene_id) %>%
    mutate (y= as.numeric (factor(iso)), event_id=as.factor(event_id)) %>%
    arrange (gene_id, event_id, iso, start)
  t= select (events_df, event_id, iso1_stop,iso2_stop,frame) %>%
    filter (frame !=9) %>%
    mutate (stop1= frame-1+(iso1_stop*3), stop2=frame-1+(iso2_stop*3) ) %>%
    select (event_id, stop1,stop2) %>%
    pivot_longer(!event_id, names_to ='iso',values_to ='distance') %>%
    mutate (iso= ifelse (iso=='stop1','Isoform 1','Isoform 2'))
  stops= both %>%
    left_join(t,by=c('event_id','iso')) %>%
    filter (is.na(distance)==F) %>%
    mutate (id= paste(event_id,y,sep='_'))
  stops2= find_genomic(stops) %>%
    dplyr::rename (genomic_stop=genomic_position) %>%
    full_join(stops,by='id') %>%
    select (event_id,y,genomic_stop) %>%
    distinct ()
  t= select (events_df, event_id ,first_mismatched,last_mismatched,frame, frameshift, iso1_aa, iso2_aa) %>%
    ungroup () %>%
    filter (frame !=9) %>%
    mutate (last_mis= last_mismatched*3) %>%        
    select (event_id, last_mis, frameshift)
  mismatches= both %>%
    left_join(t,by='event_id') %>%
    filter (is.na(last_mis)==F) %>%
    mutate (id= paste(event_id,y,sep='_'), distance=last_mis)
  last_mis= find_last_mismatch(mismatches) %>%
    dplyr::rename (genomic_last_mismatch=genomic_position) %>%
    full_join(mismatches,by='id') %>%
    select (event_id,y,strand,exons,start,end,distance, genomic_last_mismatch) %>%
    mutate (dif= end-start,inside=ifelse (genomic_last_mismatch>=start & genomic_last_mismatch<=end,1,0)) %>%
    select (event_id, y,genomic_last_mismatch) %>%
    distinct ()
  t= select (events_df, event_id ,first_mismatched,frame) %>%
    filter (frame !=9) %>%
    mutate (first_mis= frame-1+ (first_mismatched*3)) %>%         # convert AA count to DNA count
    select (event_id,first_mis) %>%
    dplyr::rename (distance=first_mis)
  mismatches= both %>%
    left_join(t,by='event_id') %>%
    filter (is.na(distance)==F) %>%
    mutate (id= paste(event_id,y,sep='_'))
  first_mis= find_genomic(mismatches) %>%
    dplyr::rename (genomic_first_mismatch=genomic_position) %>%
    full_join(mismatches,by='id') %>%
    select (event_id,y,strand,start,end,distance, genomic_first_mismatch) %>%
    mutate (dif= end-start,inside=ifelse (genomic_first_mismatch>=start & genomic_first_mismatch<=end,1,0)) %>%
    select (event_id, y, genomic_first_mismatch) %>%
    distinct ()
  iso_label= full_join(first_mis, last_mis, by=c('event_id','y')) %>%
    left_join(both, by=c('event_id','y')) %>%
    select (event_id, strand, y, exons, start, end, genomic_first_mismatch, genomic_last_mismatch) 
  pos= filter (iso_label, strand=='+') %>%
    filter (start < genomic_last_mismatch & end > genomic_first_mismatch) %>%
    mutate (newstart= pmax(start, genomic_first_mismatch), newend= pmin (end, genomic_last_mismatch) )
  neg= filter (iso_label, strand=='-') %>%
    filter (start < genomic_first_mismatch & end > genomic_last_mismatch) %>%
    mutate (newstart= pmax (start, genomic_last_mismatch) ,newend= pmin(end, genomic_first_mismatch)  )
  iso_label2= rbind(pos, neg)
  start_stop= both %>%
    group_by(event_id) %>%
    summarise(mn=min(start), mx=max(end)) %>%
    mutate (mid= mn+(mx-mn)/2, fs=mn+(mid-mn)/10,  amb=mx-(mx-mid)/5 )
  reads= samples_df %>%
    filter (event_id %in% events_df$event_id) %>%
    group_by(event_id) %>%
    summarise (iso1=sum(iso1), iso2=sum(iso2)) %>%
    mutate (iso1=paste0('n=',format (iso1, big.mark=',', trim=T)), 
            iso2=paste0('n=',format (iso2, big.mark=',', trim=T)) )
  yes_p= any(str_detect (names(events_df),'adj_pval'))
  if (yes_p==T) {events_df$p='yes'} else {events_df$p='no'}
  event_label= events_df  %>%
    mutate (p_labl= ifelse (p=='no','', format(adj_pval,scientific=T,digits=3)) ) %>%
    mutate (label= paste0(gene_id,':',event_id,':',chrm,'(',strand,')', p_labl), 
            fshift=paste0('Frameshift=',frameshift), 
            ambig=ifelse (frame==9, 'Ambiguous frame','')) %>%
    select (gene_id, event_id,type,chrm,strand, frameshift, label, fshift,ambig) %>%
    arrange (gene_id, event_id)
  f=list()
  for (i in 1:length (event_label$event_id)) {
    df= filter (both, event_id==event_label$event_id[i])
    ss= filter (stops2, event_id==event_label$event_id[i])
    iso= filter (iso_label2, event_id==event_label$event_id[i])
    sts= filter (start_stop, event_id==event_label$event_id[i])
    rd= filter (reads, event_id==event_label$event_id[i])
    f1=ggplot(df, aes(xmin=start, xmax=end, ymin=y, ymax=y+0.7)) +
      geom_rect(fill='lightblue', color='black') + 
      scale_y_continuous(limits=c(0.6,3.1),breaks=unique(df$y)+0.35,labels=unique(df$iso)) + 
      ggtitle(event_label$label[i]) +
      annotate("rect", xmin=ss$genomic_stop-2, xmax=ss$genomic_stop+2, ymin =ss$y, ymax =ss$y+0.7,
               color='red', fill="red") +
      annotate("rect", xmin=iso$newstart, xmax=iso$newend, 
               ymin =iso$y, ymax =iso$y+0.7, fill="yellow", alpha=0.2) +
      annotate ('text', x=sts$mid, y=0.7, label=rd$iso1,size=4) +
      annotate ('text', x=sts$mid, y=3, label=rd$iso2,size=4) +
      annotate ('text', x=sts$amb, y=0.7, label=event_label$ambig[i],size=4) +
      annotate ('text', x=sts$fs, y=3, label=event_label$fshift[i],size=4) + xlab(NULL)+ylab(NULL) + 
      theme(plot.title=element_text(size=16), axis.text.x=element_text(size=10, color='black'),
            axis.text.y=element_text(size=11, color='black'))
    df=  filter(samples_df, event_id==event_label$event_id[i])
    if (fig2=='jitter') {
      f2=ggplot (df, aes(y=psi, x=group, color=group))+ geom_jitter(width=0.2, size=2, height=0) +
        scale_y_continuous(limits =c(0,1)) + xlab (NULL) + ylab ('PSI') + 
        scale_color_manual (values=c('darkblue', 'darkred')) + 
        theme(legend.position = "none", axis.title.y=element_text(size=11,color='black'), 
              axis.text=element_text(size=12, color='black')) 
    } else if (fig2=='dot') {
      f2=ggplot (df, aes(x=psi))+ geom_dotplot()+ scale_x_continuous(limits=c(0,1))+
        theme(legend.position = "none") + ylab (NULL) + xlab ('PSI') 
    } 
    f[[i]]= wrap_plots(f1,f2, widths=c(4,1)) 
  }
  return (f)
}

# function to get transcript as reference for splicing event. it will be either canonical or the one that
# matches to more exons of isoform 1 (shared exons)
get_tx= function (events_df, iso1_df, ensembl_version ) {
  library (biomaRt)
  ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version=ensembl_version)
  df= dplyr::select(events_df, gene_name, gene_id, event_id)
  genes= unique (df$gene_name)
  tx= getBM (mart=ensembl, filters='ensembl_gene_id_version', genes, 
             attributes =c('ensembl_gene_id_version' , 'ensembl_transcript_id','transcript_start','transcript_end', 
                           'transcription_start_site', 'transcript_is_canonical', 'transcript_length',
                           'transcript_count', 'transcript_biotype') )
  tx_canon= filter (tx, transcript_is_canonical==1) %>%
    dplyr::select(ensembl_gene_id_version,ensembl_transcript_id, transcript_is_canonical) %>%
    dplyr::rename (gene_name=ensembl_gene_id_version, tx_id=ensembl_transcript_id, canonical=transcript_is_canonical)
  txlist= unique (tx$ensembl_transcript_id)
  exons= getBM(mart =ensembl, filters='ensembl_transcript_id', txlist,
               attributes =c('ensembl_transcript_id','ensembl_gene_id_version','chromosome_name', 'strand',
                             'ensembl_exon_id','exon_chrom_start','exon_chrom_end', 'is_constitutive')) %>%
    dplyr::select(ensembl_gene_id_version,ensembl_transcript_id,exon_chrom_start,exon_chrom_end,
                  chromosome_name,strand) %>%
    distinct () %>%
    dplyr::rename (gene_name=ensembl_gene_id_version, tx_id=ensembl_transcript_id, 
                   start=exon_chrom_start, end=exon_chrom_end ) %>%
    mutate (start=as.numeric(start), end=as.numeric(end))
  ex= iso1_df %>%
    ungroup () %>%
    dplyr::select (gene_name,event_id, start,end, exons) %>%
    filter (exons !='e2') %>%
    mutate (start=as.numeric(start), end=as.numeric(end), exons=ifelse(exons=='e4','e3',exons)) %>%
    distinct() %>%
    inner_join(exons, by=c('gene_name','start','end'),relationship="many-to-many") %>%
    dplyr::select (-c(chromosome_name,strand))
  tx2= ex %>%
    group_by(gene_name,event_id,tx_id,exons) %>%
    summarise (n=n()) %>%
    pivot_wider(names_from ='exons', values_from ='n') %>%
    mutate (n=sum(e1,e3,na.rm =T)) %>%
    left_join(tx_canon, by=c('gene_name','tx_id')) %>%
    ungroup () %>%
    arrange (gene_name, event_id, desc(n),desc(canonical)) %>%
    group_by(gene_name,event_id) %>%
    slice_head (n=1) %>%
    dplyr::select(-c(e1,e3))
  notx= dplyr::select(iso1_df, gene_name, event_id) %>%
    distinct () %>%
    anti_join(tx2,by='event_id') %>%
    left_join(tx_canon,by='gene_name') %>%
    mutate (n=0)
  tx3= rbind (tx2,notx)
  ref= exons %>%
    inner_join(tx3,by=c('gene_name','tx_id'),relationship="many-to-many") %>%
    mutate (canonical=replace_na(canonical,0), 
            strand= ifelse(strand==1,'+','-'), chrm=paste0('chr',chromosome_name),
            id=ifelse (canonical==1,'Canonical',tx_id)) %>%
    dplyr::select (gene_name,id, event_id, chrm, start, end, strand)
  t= events_df %>%
    ungroup () %>%
    dplyr::select(gene_name, gene_id) %>%
    distinct ()
  ref2= ref %>%
    left_join(t, by='gene_name') %>%
    arrange (gene_id, start)
  remove (ensembl)
  detach("package:biomaRt", unload = TRUE)
  return (ref2)
}

# zoomed-out figure with transcript 
splice_figure_ref_TX= function (events_df , positions_df, ensembl_version) {
  pos= positions_df %>%
    filter (event_id %in% events_df$event_id) %>%
    arrange (gene_name, event_id, exons, start) %>%
    mutate (type= str_extract(event_id,'[^.]+')) %>%
    dplyr::select (gene_name, event_id, type, chrm, strand, exons, start, end) 
  iso2= pos %>%
    mutate (incl= ifelse (type=='mutex_exons' & ((strand=='+' & exons !='e2') | (strand=='-' & exons !='e3')),1,0),
            iso='Isoform 2') %>%
    filter (type !='mutex_exons' | incl==1) %>%
    dplyr::select (-incl)
  iso1= pos %>%
    mutate (incl=ifelse (type=='mutex_exons' & ((strand=='+' & exons !='e3') | (strand=='-' & exons !='e2')),1,0), 
            iso='Isoform 1') %>%
    filter ((type !='mutex_exons' & exons !='e2') | incl==1) %>%
    dplyr::select (-incl)
  ours= rbind(iso1, iso2) %>%
    mutate (id= ifelse (iso=='Isoform 1', paste0(event_id,'_iso1'), paste0(event_id,'_iso2')),
            start=as.numeric(start), end=as.numeric(end)) %>%
    dplyr::select(gene_name, id, event_id, chrm, start, end, strand, iso) %>%
    dplyr::arrange (gene_name, id, start) %>%
    group_by(event_id) %>%
    mutate (y= as.numeric (factor(iso)), event_id=as.factor(event_id)) 
  ref2=get_tx (events_df, iso1, ensembl_version)
  event_label= events_df  %>%
    mutate (label= paste0(gene_id,':',event_id,':',chrm,'(',strand,')')) %>%
    dplyr::select (gene_id, event_id,type,chrm,strand, label) %>%
    arrange (gene_id, event_id)
  f=list()
  for (i in 1:length (event_label$event_id)) {
    df= filter (ours, event_id==event_label$event_id[i])
    r= filter (ref2, event_id==event_label$event_id[i])
    ev= filter (event_label, event_id==event_label$event_id[i])
    f[[i]]=ggplot(df, aes(xmin=start, xmax=end, ymin=y+1, ymax=y+1.8)) +
      geom_rect(fill='darkgoldenrod1', color='black') + 
      scale_y_continuous(breaks=c(1.4,unique(df$y)+1.4),labels=c(unique(r$id),unique(df$iso))) + 
      ggtitle(ev$label) +
      annotate("rect", xmin=r$start, xmax=r$end, ymin =1, ymax =1.8,
               color='black', fill="darkgrey")+
      theme (plot.title =element_text(size=16), axis.text.x =element_text(size=10,color='black'),
             axis.text.y=element_text(size=11, color='black'))
  }
  return (f)
}

# Protein domains. input gene names in ensembl gene ID version
get_domains= function (gene_names, ensembl_version) {
  library (biomaRt)
  # human genes 
  ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version=ensembl_version)
  tx= getBM (mart=ensembl, filters='ensembl_gene_id_version', unique (gene_names), 
             attributes =c('ensembl_gene_id_version', 'external_gene_name' , 'ensembl_transcript_id', 'transcript_is_canonical') ) %>%
    filter (transcript_is_canonical==1)
  pp= getBM (mart=ensembl, filters='ensembl_transcript_id', unique (tx$ensembl_transcript_id),
             attributes =c('ensembl_transcript_id','external_gene_name' ,
                           'interpro' , 'interpro_short_description',
                           'interpro_description', 'interpro_start', 'interpro_end'))
  domains= filter (pp, interpro !='') %>%
    dplyr::rename (tx_id=ensembl_transcript_id) %>%
    arrange (tx_id, interpro, interpro_start,interpro_end) %>%
    group_by(tx_id, interpro) %>%
    mutate (dif=interpro_start-lag(interpro_end )) %>%
    filter (is.na(dif)==T | dif >=0) %>%
    dplyr::select (-dif) %>%
    ungroup () %>%
    mutate (id=paste0(tx_id,interpro,interpro_start, interpro_end)) %>%
    arrange (tx_id, interpro_start, interpro_end)
  df= data.frame (id=character(), y=numeric() )
  for (j in unique (domains$tx_id)) {
    t=filter (domains, tx_id==j) %>%
      ungroup () %>%
      arrange (tx_id, interpro_start, interpro_end)
    for (i in 1:nrow(t)) {
      id= t$id[i]
      if (i==1) {
        y=2
        end2=t$interpro_end[i]
        end3=0
        end4=0
        end5=0
        end6=0
        end7=0
        end8=0
      } else {
        if (t$interpro_start[i]>=end2) {
          y=2
          end2=t$interpro_end[i]
        } else if (t$interpro_start[i]>=end3) {
          y=3
          end3=t$interpro_end[i]
        } else if (t$interpro_start[i]>=end4) {
          y=4
          end4=t$interpro_end[i]
        } else if (t$interpro_start[i]>=end5) {
          y=5
          end5=t$interpro_end[i]
        } else if (t$interpro_start[i]>=end6) {
          y=6
          end6=t$interpro_end[i]
        } else if (t$interpro_start[i]>=end7) {
          y=7
          end7=t$interpro_end[i]
        } else {
          y=8
          end8=t$interpro_end[i]
        }
      }
      df= add_row(df, id=id, y=y )
    }
  }
  domains2= left_join(domains, df, by='id') %>%
    dplyr::select (-id) 
  coding=getSequence(mart=ensembl, type='ensembl_transcript_id', id=unique(domains2$tx_id), seqType ='coding')
  cdna=getSequence(mart=ensembl, type='ensembl_transcript_id', id=unique(domains2$tx_id), seqType ='cdna')
  t= inner_join(cdna,coding,by='ensembl_transcript_id')
  tt= as.data.frame(str_locate (t$cdna,t$coding))
  tss= bind_cols(t,tt) %>%
    mutate (tss_aa=round (start/3,2)) %>%
    dplyr::select (ensembl_transcript_id, tss_aa)
  domains3= left_join(domains2,tss,by=c('tx_id'='ensembl_transcript_id')) %>%
    mutate (interpro_start=interpro_start+tss_aa, interpro_end=interpro_end+tss_aa)
  remove(ensembl)
  detach("package:biomaRt", unload = TRUE)
  return (domains3)
}

# Process exons
get_exons= function (gene_names, ensembl_version) {
  library (biomaRt)
  ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version=ensembl_version)
  tx= getBM (mart=ensembl, filters='ensembl_gene_id_version', unique (gene_names), 
             attributes =c('ensembl_gene_id_version', 'external_gene_name' , 'ensembl_transcript_id', 
                           'transcript_is_canonical','chromosome_name', 'strand') ) %>%
    filter (transcript_is_canonical==1)
  exons= getBM(mart =ensembl, filters='ensembl_transcript_id', unique(tx$ensembl_transcript_id),
               attributes =c('ensembl_transcript_id','ensembl_gene_id_version','external_gene_name','chromosome_name', 'strand',
                             'ensembl_exon_id','exon_chrom_start','exon_chrom_end', 'is_constitutive')) %>%
    dplyr::select(ensembl_gene_id_version,external_gene_name,ensembl_transcript_id,exon_chrom_start,exon_chrom_end,
                  chromosome_name,strand) %>%
    distinct () %>%
    dplyr::rename (gene_name=ensembl_gene_id_version, tx_id=ensembl_transcript_id, 
                   start=exon_chrom_start, end=exon_chrom_end ) %>%
    mutate (start=as.numeric(start), end=as.numeric(end), 
            strand= ifelse(strand==1,'+','-'), chrm=paste0('chr',chromosome_name)) %>%
    dplyr::select (-chromosome_name) %>%
    mutate (len= end-start)
  minus= filter(exons, strand=='-')
  minus_df= data.frame (tx_id=character(), new_start=numeric(),new_end=numeric() )
  for (t in unique (minus$tx_id)) {
    df=filter (minus, tx_id==t) %>%
      ungroup %>%
      arrange (desc(start))
    for (i in 1:length (df$len)) {
      if (i==1) {
        newstart=0
        newend=df$len[i]
      } else {
        newstart=newend
        newend=newend+df$len[i]
      }
      minus_df= add_row(minus_df, tx_id=t, new_start=newstart, new_end=newend )
    }
  }
  plus= filter(exons, strand=='+')
  plus_df= data.frame (tx_id=character(), new_start=numeric(),new_end=numeric() )
  for (t in unique (plus$tx_id)) {
    df=filter (plus, tx_id==t) %>%
      ungroup %>%
      arrange (start)
    for (i in 1:length (df$len)) {
      if (i==1) {
        newstart=0
        newend=df$len[i]
      } else {
        newstart=newend
        newend=newend+df$len[i]
      }
      plus_df= add_row(plus_df, tx_id=t, new_start=newstart, new_end=newend )
    }
  }
  out_df= rbind (plus_df, minus_df)
  ex= out_df %>%
    group_by(tx_id) %>%
    mutate (y= as.numeric (factor(tx_id))) %>%
    left_join(tx,by=c('tx_id'='ensembl_transcript_id')) %>%
    dplyr::rename (gene_name=ensembl_gene_id_version, gene_id=external_gene_name) %>%
    mutate (strand= ifelse(strand==1,'+','-'), chrm=paste0('chr',chromosome_name)) %>%
    dplyr::select (-c(transcript_is_canonical, chromosome_name))
  remove(ensembl)
  detach("package:biomaRt", unload = TRUE)
  return (ex)
}

# protein domain figure
domain_fig= function (ensembl_genes, ensembl_version ) {
  domains_df= get_domains(ensembl_genes, ensembl_version)
  exons_df= get_exons (ensembl_genes, ensembl_version)
  exons= exons_df %>%
    mutate (aa_start=new_start/3, aa_end=new_end/3)
  domains= domains_df %>%
    mutate (interpro_description=paste0(interpro,':',interpro_description))
  gene_lbl= as.data.frame (exons) %>%
    ungroup() %>%
    dplyr::select (gene_id, tx_id, chrm, strand) %>%
    distinct () %>%
    mutate(lbl= paste0(gene_id,':',chrm,'(',strand,')')) %>%
    arrange (gene_id)
  dd_lbl= domains %>%
    ungroup () %>%
    dplyr::select (tx_id, interpro, interpro_description) %>%
    distinct ()
  f= list()
  for (i in 1:length(gene_lbl$tx_id)) {
    ex= filter (exons, tx_id==gene_lbl$tx_id[i]) 
    dd= filter (domains, tx_id==gene_lbl$tx_id[i]) 
    gg= filter (gene_lbl, tx_id==gene_lbl$tx_id[i])
    dd_lbls= filter (dd_lbl, tx_id==gene_lbl$tx_id[i])
    f[[i]]= ggplot(dd, aes(xmin=interpro_start, xmax=interpro_end, ymin=y, ymax =y+0.8, fill=interpro)) +
      geom_rect (alpha=1) + 
      scale_fill_discrete (name=NULL, limits=dd_lbls$interpro, labels=dd_lbls$interpro_description) + 
      theme(legend.position = "bottom", legend.direction="vertical",legend.justification='left') +
      ggtitle (paste0(gg$gene_id,':Exons+Protein Domains')) +
      annotate("rect", xmin=ex$aa_start, xmax=ex$aa_end, ymin =ex$y, ymax =ex$y+0.8, 
               color='black',fill='darkgrey')+
      scale_y_continuous(breaks = NULL)+ scale_x_continuous(breaks = NULL)+
      theme (plot.title =element_text(size=16), legend.text =element_text(size=14))
  }
  return (f)
}
