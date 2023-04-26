setwd("~/Documents/GitHub/GUTS_MGS_Psych/")

####Load libraries####

library(tidyverse)
library(vegan)
library(ape)
library(MetBrewer) #Lets plot with some art :)
library(ggrepel)
ggplot <- function(...) ggplot2::ggplot(...) + scale_color_manual(values=met.brewer("Manet",5)[c(1,5, 3)] ) +
  scale_fill_manual(values=met.brewer("Manet", 5)[c(1,5)] )

###############
#### DATA ####
##############
Clean_path = function(C, Keep){
  'Clean pathway IDs'
  colnames(C)[1] = "ID"
  sapply(colnames(C), function(x){ str_split(x, "_knead")[[1]][1] } ) -> New_names
  colnames(C) = New_names
  #C %>% select(c("ID", Keep$ID)) -> C
  return(C)
}
Format_path = function(Path){
  'Format Pathways from HM3 to format for script'
  Path %>% filter(! grepl("#", ID)) -> Path
  Path %>% filter(! grepl("\\|",ID )) -> Path
  t(Path) %>% as.data.frame() %>% rownames_to_column("ID") %>% as_tibble() %>% filter(! ID == "ID") -> Path_2
  colnames(Path_2) = c("ID",  Path$ID)
  return(Path_2)
}
Format_MP3_table = function( Path_input, Path_output ){
  'Format MP3 table to script input'
  Path_input = read_tsv(Path_input, skip = 1) 
  Path_input %>% select(-clade_taxid) %>% as.data.frame() %>% column_to_rownames("clade_name") %>% t() %>% as_tibble() %>% mutate(ID = colnames(Path_input)[3:dim(Path_input)[[2]]], .before=1 ) -> Path_input
  colnames(Path_input) %>% sapply( function(x){ y = str_split(x, "\\|")[[1]] ;if (length(y) < 2 ){ return(y) } else{ return(y[length(y)]) }  } ) -> New_names
  colnames(Path_input) = as.vector(New_names)
  Path_input$ID =  str_replace(Path_input$ID, "_metaphlan", "")
  write_tsv(Path_input, Path_output)
}

#motus DB 2.6.0, motus 3.0.3 on python 3.7.15

#0. Samples info
##0.1 Metadata samples
###Cases
Covariates_GUTS = read_tsv("Data/MetaData/Metadata_GUTS_final2.tsv") #%>% dplyr::rename(Sample_ID = ID)
####Keep only baseline samples
Covariates_GUTS %>% dplyr::filter(Treatment_time == "V2") -> Covariates_GUTS
####Check repeated samples  
Covariates_GUTS %>%  group_by(Sample_ID, Batch) %>% summarise(N = n() ) %>% arrange(desc(N)) -> Repeated #19 repeated samples, 18 in batch 2, 1 in batch 1
###Controls
Covariates_control = read_tsv("Data/MetaData/Covariates_population.tsv") %>% drop_na() #dplyr::rename(Sequencing_ID = ID, Read_number = META.DNA.postclean.reads, Batch=META.BATCH, Sex=Gender) %>% drop_na()
##0.2 Sequencing info
Reads_info = read_csv("Data/Data2/All_kneaddata_reads_stats2_F.csv") %>% dplyr::rename(ID = Sample)
Reads_info %>% dplyr::select( c("ID",colnames(Reads_info)[grepl("p1", colnames(Reads_info))]) ) %>% gather(Source, Reads, 2:6) -> Reads_info2
##0.3 Taxonomy abundance
###Previous work used MP3 profiles, preprocissing is included inside a function. Not used for publication
Previous_BioBakery = function(){
  #1. Taxonomy abundance: Metaphlan 3
  ####MP3 format has been previously pre-processed. Format for the script is: ID column, participant per row, column per taxa. Unknown estiamtion included. Only final part of the taxonmy reported
  Case1 = read_tsv("Data/MP3/Cases_batch1.tsv") #58 samples, 508 taxonomies (in at least 1 sample)
  Case1 %>% select(-ID) %>% apply(2, function(x){ length(x[!x==0]) } )  %>% as.data.frame() %>% rownames_to_column("Taxonomy") %>% as_tibble() %>% dplyr::rename(N = ".") %>% arrange(N) 
  Case2 = read_tsv("Data/MP3/Cases_batch2.tsv") #150 samples. 570 taxonomies
  Case2 %>% select(-ID) %>% apply(2, function(x){ length(x[!x==0]) } )  %>% as.data.frame() %>% rownames_to_column("Taxonomy") %>% as_tibble() %>% dplyr::rename(N = ".") %>% arrange(N) 
  Control =  read_tsv("Data/MP3/Controls_mp3.tsv") %>% drop_na()
  colnames(Control) = sapply(colnames(Control), function(x) {  y = str_split(x, "\\|" )[[1]] ; y[length(y)] })
  # Merging the sequecing data and metadata
  Case1$ID %>% sapply(. , function(x){ substr(str_split(x, "_")[[1]][2],2 ,5) } ) %>%
    as.data.frame() %>% rownames_to_column("Sequencing_ID") %>% as_tibble() %>% dplyr::rename(Sample_ID = ".") %>% mutate(Batch = 1) -> ID_info1
  Case2$ID %>% sapply(. , function(x){ if (grepl("V2", x) ){ substr(str_split(x, "_")[[1]][2],2 ,5) } else{ substr(str_split(x, "_")[[1]][1], 3,6) }  } ) %>%
    as.data.frame() %>% rownames_to_column("Sequencing_ID") %>% as_tibble() %>% dplyr::rename(Sample_ID = ".") %>% mutate(Batch = 2) -> ID_info2
  rbind(ID_info1, ID_info2) %>% mutate(Treatment_time = ifelse(grepl("V5",Sequencing_ID ), "V5", "V2" ) ) -> ID_info
  
  plyr::rbind.fill(Case1, Case2) %>%as_tibble() %>% filter(ID %in% ID_info$Sequencing_ID) -> Case
  ID_info %>% left_join(.,Covariates_GUTS) -> Covariates_GUTS
  write_tsv(Covariates_GUTS, "Data/MetaData/Metadata_GUTS_final2.tsv")
  #Include sequencing depth
  read_tsv("Data/MetaData/readdepth_batch1.tsv") %>% rename(Sequencing_ID = ID) -> Batch1_reads
  read_csv("Data/MetaData/results_kneaddata_readdepth_batch2.csv") %>% select(Sample, `Final (p1)`) %>% rename(Sequencing_ID = Sample, Read_number= `Final (p1)`) -> Batch2_reads
  left_join( Covariates_GUTS, rbind(Batch1_reads, Batch2_reads) ) -> Covariates_GUTS
  Covariates_GUTS %>% ggplot(aes(x=Read_number, fill=as.factor(Batch))) + geom_density(alpha=0.4) + theme_bw()
  
  #2. Pathway abundance
  Case_ptw1 = Format_path(Clean_path(read_tsv("Data/HM3/Pathways_cases_batch1.txt") )) %>% filter(ID %in% Covariates_GUTS$Sequencing_ID )
  Case_ptw2 = Format_path(Clean_path(read_tsv("Data/HM3/Pathways_cases_batch2.txt") )) %>% filter(ID %in% Covariates_GUTS$Sequencing_ID )
  plyr::rbind.fill(Case_ptw1, Case_ptw2) %>%as_tibble() -> Case_pathway
  Case_pathway %>%select(-ID) %>% apply(2, as.numeric) %>% as_tibble() %>% mutate(ID = Case_pathway$ID, .before=1) -> Case_pathway
  
  #Data cases
  Case_pathway[is.na(Case_pathway)] = 0
  Case[is.na(Case)] = 0
  
}
###Info about which sample is matching with each batch, in case analysis is done separately per batch
read_tsv("Data/MetaData/Match_with_info_batchTime.tsv") %>% select(ID, Match)  -> BatchInfoMatch
####mOTUS 3.0.3 taxonomic profiles (min marker number 4)
read_tsv("Data/Data2/Merged_taxonomy2.tsv") -> Taxonomic_profile
Taxonomic_profile %>% t()  %>% as.data.frame()  %>% rownames_to_column("ID") %>% as_tibble() %>% filter(! grepl("#", ID)) %>% `colnames<-`(c("ID",Taxonomic_profile$`#mOTUs2_clade`)) -> Taxonomic_profile2
Taxonomic_profile2 %>% select(-ID) %>% apply(2, as.numeric) %>% as_tibble() %>% mutate(ID = Taxonomic_profile2$ID, .before=1) -> Taxonomic_profile2

####OTUS 3.0.3 taxonomic profiles (min marker number 4) --> For replication
read_tsv("Data/Data2/Replication/taxonomy_replication.tsv") -> Taxonomic_replication
Taxonomic_replication %>% t()  %>% as.data.frame()  %>% rownames_to_column("ID") %>% as_tibble() %>% filter(! grepl("#", ID)) %>% `colnames<-`(c("ID",Taxonomic_replication$`#mOTUs2_clade`)) -> Taxonomic_replication2
Taxonomic_replication2 %>% select(-ID) %>% apply(2, as.numeric) %>% as_tibble() %>% mutate(ID = Taxonomic_replication2$ID, .before=1) -> Taxonomic_replication2
readxl::read_excel("Data/Data2/Replication/41467_2020_15457_MOESM4_ESM.xlsx") -> Cov_rep
colnames(Cov_rep)[1] = "ID_sample"
Names = read_tsv("Data/Data2/Replication/Change_names.tsv", col_names = F, ) -> Names_rep ; colnames(Names_rep) = c("ID", "ID_sample") 
left_join(Cov_rep, Names_rep) -> Cov_rep
Cov_rep %>% select(ID, Age,  `Gender\r\n(1:male, 2:female)`, `BMI...8`, Group) %>% mutate(Status =  ifelse(Group=="HC", 0, 1) ) %>% rename(Sex=`Gender\r\n(1:male, 2:female)`, BMI=`BMI...8` ) -> Cov_rep
read_csv("Data/Data2/Replication/Replication_kneaddata_reads_stats.csv") %>% rename(ID = Sample) -> reads_replication
left_join(Cov_rep , reads_replication %>%  dplyr::rename(Reads_raw =`Raw reads (p1)`, Reads_clean=`Final (p1)`, Reads_human=`Contaminants (p1)`) %>% select(ID, Reads_clean, Reads_raw, Reads_human) ) -> Cov_rep


#### replication using DMP participants with SZ or BPD
read_tsv("Data/Data2/Replication2/Merged_replication2.tsv") -> Taxonomic_replication_DMP
Taxonomic_replication_DMP %>% t()  %>% as.data.frame()  %>% rownames_to_column("ID") %>% as_tibble() %>% filter(! grepl("#", ID)) %>% `colnames<-`(c("ID",Taxonomic_replication_DMP$`#mOTUs2_clade`)) -> Taxonomic_replication_DMP2
Taxonomic_replication_DMP2 %>% select(-ID) %>% apply(2, as.numeric) %>% as_tibble() %>% mutate(ID = Taxonomic_replication_DMP2$ID, .before=1) -> Taxonomic_replication_DMP2
read_csv("Data/Data2/Replication2/Stats_knead.csv") %>% rename(ID = Sample) -> reads_replication2
tibble(ID = Taxonomic_replication_DMP2$ID,  Age = 0, Sex=0, BMI=0, Status=1 ) %>% left_join(. , reads_replication2 %>%  dplyr::rename(Reads_raw =`Raw reads (p1)`, Reads_clean=`Final (p1)`, Reads_human=`Contaminants (p1)`) %>% select(ID, Reads_clean, Reads_raw, Reads_human)   ) -> Cov_rep2


##0.4 Merge taxonomy table and metadata
Taxonomic_profile2 %>% filter(ID %in% c(Covariates_GUTS$Sequencing_ID, Covariates_control$Sequencing_ID) ) -> Taxonomic_profile3
##0.5 Selection of replicates: take the one with more alpha diversity 
Taxonomic_profile3 %>% filter(ID %in% Covariates_GUTS$Sequencing_ID) %>% select(colnames(Taxonomic_profile3)[grepl("s__", colnames(Taxonomic_profile3))])  %>% diversity(index = "shannon") -> Alpha_diversity
left_join(Covariates_GUTS, tibble(Sequencing_ID = filter(Taxonomic_profile3, ID %in% Covariates_GUTS$Sequencing_ID)$ID, Shannon = Alpha_diversity)) -> Covariates_GUTS
Repeated %>% filter(N==2) -> To_select
No_Choices = c()
for (I in unique(To_select$Sample_ID)){
  Covariates_GUTS %>% filter(Sample_ID == I) %>% select(Sequencing_ID,Shannon) %>% arrange(desc(Shannon)) %>% tail(1) -> No_Choice
  No_Choices = c(No_Choices, No_Choice$Sequencing_ID)
}
Covariates_GUTS %>% filter(! Sequencing_ID %in% No_Choices) -> Covariates_GUTS
Taxonomic_profile3 %>% filter(! ID %in% No_Choices) -> Taxonomic_profile4
##0.6 Case/Control analysis: Compile covariates
Covariates_GUTS %>% select(Sequencing_ID, Batch, Age, Sex, BMI) %>% mutate(Batch = paste("GUTS_",Batch, sep=""), Status=1 ) -> Covariates_GUTS
Covariates_control %>% select(Sequencing_ID, Batch, Age, Sex, BMI) %>% mutate(Sex=ifelse(Sex=="F", 2, 1), Status=0 ) -> Covariates_control
Covariates_all = rbind(Covariates_GUTS, Covariates_control)
left_join(Covariates_all, Reads_info %>% dplyr::rename(Sequencing_ID = ID, Reads_raw =`Raw reads (p1)`, Reads_clean=`Final (p1)`, Reads_human=`Contaminants (p1)`) %>% select(Sequencing_ID, Reads_clean, Reads_raw, Reads_human) ) -> Covariates_all
###Check that matching worked (no significant differences in covariates)
Covariates_all %>% glm( as.factor(Status) ~ Age + Sex + BMI, . , family=binomial()  ) %>% summary() #Data is matched
###Check differences in read number
####Check outliers
Covariates_all %>% arrange(Reads_clean) %>% select(Sequencing_ID, Reads_clean, Reads_human, Reads_raw) %>% filter(Reads_clean < mean(Covariates_all$Reads_clean) - 3*sd(Covariates_all$Reads_clean) ) -> Outliers
####Plots
Covariates_all %>% mutate(Batch = ifelse(grepl("dag3", Batch), "DMP", Batch ) ) %>% filter(!Sequencing_ID %in% Outliers$Sequencing_ID) %>% ggplot(aes(x=Reads_clean, fill=Batch)) + geom_density(alpha=0.4) + theme_bw() + scale_fill_manual(values=met.brewer("Manet", 5)[c(1,5, 3)] )
Covariates_all %>% mutate(Batch = ifelse(grepl("dag3", Batch), "DMP", Batch ) ) %>% filter(!Sequencing_ID %in% Outliers$Sequencing_ID) %>% ggplot(aes(x=Reads_raw, fill=Batch)) + geom_density(alpha=0.4) + theme_bw() + scale_fill_manual(values=met.brewer("Manet", 5)[c(1,5, 3)] )
####Statistical testing
Covariates_all %>% mutate(Batch = ifelse(grepl("dag3", Batch), "DMP", Batch ) ) %>% aov(log10(Reads_clean) ~ Batch, data = .) -> Read_differences
Read_differences %>% summary() ; TukeyHSD(Read_differences)

#1. Data exploration
##1.1 Compute prevalences
Comp_Prevalence = function(DF){
  #Compute prevalence of different taxa
  apply(DF,2, function(x){ sum(!as.numeric(x)==0)/length(x) } ) -> Prevalence_vector
  apply(DF,2, function(x){ sum(!as.numeric(x)==0) } ) -> N
  Prevalence_df = tibble(Bug = names(Prevalence_vector), Prevalence = Prevalence_vector, N = N)
  return(Prevalence_df)
}
Comp_Prevalence( select( filter(Taxonomic_profile4, ID %in% Covariates_GUTS$Sequencing_ID )  , -ID) ) -> Prevalence_cases
Comp_Prevalence( select(filter(Taxonomic_profile4, ID %in% Covariates_control$Sequencing_ID ), -ID) ) -> Prevalence_controls
###Not the same taxonomy detected, so we keep only taxa seen in both
full_join(Prevalence_cases, Prevalence_controls, by="Bug", suffix = c("_case", "_control") ) %>% drop_na() -> Prevalence
Prevalence %>% ggplot(aes(x=Prevalence_case, y=Prevalence_control)) + geom_point() +theme_bw()+ geom_abline()
Prevalence %>% filter(! (Prevalence_case == 0 | Prevalence_control == 0) ) -> Prevalence
Taxonomic_profile4 %>% select(c(ID, Prevalence$Bug)) -> Data
###Add annotation about which cohort does it belong to
Data %>% left_join(. , dplyr::rename(Covariates_all, ID = Sequencing_ID  )) %>% drop_na()  -> Data_c
Data %>% filter(ID %in% Data_c$ID) -> Data

colnames(Data) %>% sapply(function(x){ str_split(x, "\\|")[[1]] -> y ; y[length(y)] } ) -> New_names
##1.2 Transform data
###Transformation of the data. Microbiome data is compositional, and thus, log-ratios might be used to get all microbial features in the same reference framework
Pseudocount = function(Taxa){
  return(min(Taxa[!Taxa==0])/2)
}
Geom_mean = function(x){
  exp(mean(log(x)))
}
CLR = function(D){
  log2(D / Geom_mean(D) )
  
}
rCLR = function(D){
  log(D / Geom_mean(D[!D==0]) )
  
}
Compute_CLR_taxonomy = function(Data, Taxonomy = "s", Keep_column = "UNKNOWN", Do_pseudo=F){
  if (!Taxonomy == "Path"){
    paste(c(Taxonomy, "__"), collapse="" ) -> Ta
    colnames(select(Data, - one_of(c("ID", Keep_column )))) -> Taxa
    #Taxa[grepl(Ta, Taxa)] -> Taxa
    Select_taxonomy_level(Taxa, Ta) -> Taxa
    Data %>% select( one_of(c("ID", Keep_column, Taxa) )) -> Data_taxonomy
  } else { Data_taxonomy = Data}
  
  if (Do_pseudo == T){
  #computation of a different pseudocoount per taxa.
  apply(select(Data_taxonomy, -ID), 2, function(taxa){
    if( sum(taxa==0) >= 1){ pscount =  Pseudocount(taxa) 
    } else { pscount = 0 }
    pscount } ) -> Pseudocounts
  apply(select(Data_taxonomy, -ID), 1, function(Participant){ 
    P = CLR(Participant + Pseudocounts)
    P
  }) %>% t() %>% as_tibble() %>% mutate(ID = Data_taxonomy$ID, .before=1) -> Data_taxonomy
  } else { 
    apply(select(Data_taxonomy, -ID), 1, function(Participant){ 
      P = rCLR(Participant) } ) %>% t() %>% as_tibble() %>% mutate(ID = Data_taxonomy$ID, .before=1) -> Data_taxonomy
    }  
  Data_taxonomy %>% select(-one_of(Keep_column)) -> Data_taxonomy
  return(Data_taxonomy)
}  
Select_taxonomy_level = function(Data, Pattern){
  Data %>% sapply(function(x){ str_split(x, "\\|")[[1]] -> y ; y[length(y)] } ) -> New_names
  Data[grepl(Pattern, New_names)] -> Data
  return(Data)
}
Do_CLR_all_taxonomy = function(Data){
  Compute_CLR_taxonomy(Data, Taxonomy = "s", Do_pseudo = T, Keep_column = "") -> Species_data
  All_taxonomy = Species_data
  for (i in c("p", "c", "o", "f", "g")){
    Compute_CLR_taxonomy(Data, Taxonomy = i, Do_pseudo = T, Keep_column = "") -> Transformed_data
    left_join(All_taxonomy, Transformed_data, by = "ID") -> All_taxonomy
  }    
  All_taxonomy <- Map(function(x) replace(x, is.infinite(x), NA), All_taxonomy) %>% as_tibble()
  Species_data <- Map(function(x) replace(x, is.infinite(x), NA), Species_data) %>% as_tibble()
  return(list(All_taxonomy, Species_data))
}

#CLR on data
Do_CLR_all_taxonomy(Data) -> Data_CLR
All_taxonomy = Data_CLR[[1]] ; Species_data = Data_CLR[[2]]
#CLR on validation
Do_CLR_all_taxonomy(Taxonomic_replication2) -> Data_CLR_replication1
Do_CLR_all_taxonomy(Taxonomic_replication_DMP2) -> Data_CLR_replication2



##1.3 Beta diversity
###Compute beta diversity at the species level
Beta_taxonomy = function(Data_t, Data_c, Data, Meta = "Reads_clean", Distance="robust_ait"){
  "s__" -> Ta
  colnames(select(Data, -c("ID"))) -> Taxa
  Taxa[grepl(Ta, Taxa)] -> Taxa
  Data %>% select( c("ID", Taxa) ) -> Data_taxonomy
  
  
  #Compute diversity
  if (Distance == "robust_ait"){
    vegdist(select(Data_taxonomy, -c(ID)), method = "robust.aitchison") -> Beta_diversity
  } else {
    vegdist(select(Data_taxonomy, -c(ID)), method = "euclidean") -> Beta_diversity
  }
  
  Compare_BetaDiv = F
  #Comparison of beta divesities using different metrics
  if (Compare_BetaDiv == T){
    #Euclidean distance in CLR-transformed data corresponds to Aitchison distance.
    vegdist(select(Data_t, -ID), method = "euclidean") -> Beta_diversity1
    #New vegan includes robust aitchison computation. rCLR data is used with USV^T for decomposition.
    vegdist(select(Data_taxonomy, -c(ID)), method = "robust.aitchison") -> Beta_diversity2
    #Check Bray-Curtis
    vegdist(select(Data_taxonomy, -c(ID)), method = "bray") -> Beta_diversity3
    #Check Aitchison with imputation
    Data_im = zCompositions::cmultRepl(as.data.frame(select(Data_taxonomy, -ID))  %>% dplyr::select(one_of((Prevalence %>% filter(Prevalence_case> 0.1))$Bug)))
    vegdist( Data_im , method = "aitchison") -> Beta_diversity4
  
    mantel(Beta_diversity1, Beta_diversity2) -> Corr
    mantel(Beta_diversity1, Beta_diversity3) -> Corr2
    mantel(Beta_diversity2, Beta_diversity3) -> Corr3
    mantel(Beta_diversity2, Beta_diversity4) -> Corr4
    mantel(Beta_diversity1, Beta_diversity4) -> Corr5
    mantel(Beta_diversity3, Beta_diversity4) -> Corr6
  
  
    tibble(Comparison = c("Ait-rAit","Ait-BC", "rAit-BC", "rAit-imAit", "Ait-imAit", "BC-imAit"), Rho = c(Corr$statistic,Corr2$statistic,Corr3$statistic,Corr4$statistic,Corr5$statistic,Corr6$statistic ),
         P=c(Corr$signif,Corr2$signif,Corr3$signif, Corr4$signif,Corr5$signif,Corr6$signif ) ) -> Comparison_betas
  }
  
  
  #A PCOA
  print("Computing PCoA with both cases and controls")
  pcoa(Beta_diversity) -> summary_pcoa
  PCs = as_tibble(summary_pcoa$vectors)
  Variab = as_tibble(summary_pcoa$values %>% rownames_to_column("PC") )
  PCs %>% mutate(ID = Data_c$ID, Status=Data_c$Status) -> PCs
  PCs <- merge(PCs,aggregate(cbind(mean.x=Axis.1, mean.y=Axis.2)~Status,PCs,mean),by="Status")
  PCs %>% ggplot(aes(x=Axis.1, y=Axis.2, col=as.factor(Status))) + geom_point() +
  xlab(paste( c("PC1 (", as.character(100*round(Variab$Relative_eig[1],3)), "%)" ), collapse="" ) ) +
  ylab(paste( c("PC2 (", as.character(100*round(Variab$Relative_eig[2],3)), "%)" ), collapse="" ) ) + 
    theme_bw() + stat_ellipse() +
    geom_point(aes(x=mean.x,y=mean.y),size=5)+
    geom_segment(aes(x=mean.x, y=mean.y, xend=Axis.1, yend=Axis.2), alpha=0.5) -> Fig
  print(Fig)
  
  #PCOA cases only
  print("Computing PCoA with only cases")
  Data %>% filter(ID %in%  (Data_c %>% filter(Status==1))$ID ) ->Data_cases
  if (Distance == "robust_ait"){
  vegdist(select(Data_cases, -ID), method = "robust.aitchison") -> Beta_diversity2
  } else {
    vegdist(select(Data_cases, -ID), method = "euclidean") -> Beta_diversity2
  }
  pcoa(Beta_diversity2) -> summary_pcoa2
  PCs2 = as_tibble(summary_pcoa2$vectors)
  Variab2 = as_tibble(summary_pcoa2$values %>% rownames_to_column("PC") )
  PCs2%>% mutate(ID = (Data_c %>% filter(Status==1))$ID, Batch =(Data_c %>% filter(Status==1))$Batch) -> PCs2
  PCs2 <- merge(PCs2,aggregate(cbind(mean.x=Axis.1, mean.y=Axis.2)~Batch,PCs2,mean),by="Batch")
  PCs2 %>% ggplot(aes(x=Axis.1, y=Axis.2, col=as.factor(Batch))) + geom_point() +
    xlab(paste( c("PC1 (", as.character(100*round(Variab2$Relative_eig[1],2)), "%)" ), collapse="" ) ) +
    ylab(paste( c("PC2 (", as.character(100*round(Variab2$Relative_eig[2],2)), "%)" ), collapse="" ) ) + 
    theme_bw() + stat_ellipse() +
    geom_point(aes(x=mean.x,y=mean.y),size=5)+
    geom_segment(aes(x=mean.x, y=mean.y, xend=Axis.1, yend=Axis.2), alpha=0.5) -> Fig2
  print(Fig2)
  
  
  
  
  #B PERMANOVA
  print("Computing PERMANOVA with cases and controls")
  Fo = paste0("Beta_diversity ~", paste(Meta, "Status", sep="+" )) 
  adonis2( as.formula(Fo) ,Data_c, permutations = 5000) %>% print()
  
  #B PERMANOVA batch
  print("Computing PERMANOVA only with cases")
  (Data_c %>% filter(Status==1))-> Data_m
  #Fo = paste0("Beta_diversity2 ~", paste(Meta, "Batch", sep="+" )) 
  #adonis2(Fo, Data_m , permutations = 5000) %>% print()
  adonis2(Beta_diversity2 ~  Data_m$Reads_clean + Data_m$Batch  , permutations = 5000) %>% print()
  
}
Do_beta_validation = function(Data, Data_r1=Taxonomic_replication2, Data_r2=Taxonomic_replication_DMP2, Cov=Covariates_all, Cov_r1=Cov_rep, Cov_r2=Cov_rep2){
  ###Beta diversity replication
  Data %>% select(-ID) %>% Comp_Prevalence() -> Prevalence_guts
  Data_r1 %>% select(-ID) %>% Comp_Prevalence() -> Prevalence_replication
  Data_r2 %>% select(-ID) %>% Comp_Prevalence() -> Prevalence_replication2
  
  ###Not the same taxonomy detected, so we keep only taxa seen in both
  full_join(Prevalence_replication, Prevalence_guts, by="Bug", suffix = c("_guts", "_validation") ) %>% drop_na() -> Prevalence_replication
  Prevalence_replication %>% ggplot(aes(x=Prevalence_guts, y=Prevalence_validation)) + geom_point() +theme_bw()+ geom_abline()
  Prevalence_replication %>% filter(! (Prevalence_guts == 0 | Prevalence_validation == 0) ) -> Prevalence_replication
  
  rbind(select(Data, c("ID", Prevalence_replication$Bug)) , select(Data_r1, c("ID", Prevalence_replication$Bug)) ) -> Data_integrated
  
  #Include second validation
  Prevalence_replication$Bug[Prevalence_replication$Bug %in%  colnames(Taxonomic_replication_DMP2) ] -> Bugs_rep2
  select(Data_integrated, c("ID", Bugs_rep2)) %>% rbind(. , select(Taxonomic_replication_DMP2, c("ID", Bugs_rep2)) ) -> Data_integrated2
  
  #Beta diversity using only the first validation
  rbind(Cov_rep %>%select(-Group) %>% mutate(Batch="Zhu"), Cov %>% rename(ID = Sequencing_ID)  %>% mutate(Batch="GUTS+DMP") ) -> Cov_integrated
  Beta_taxonomy(NA, Data_c = left_join(Data_integrated,Cov_integrated) , Data_integrated, Meta="Batch")
  #Beta diversity using both validations
  rbind( Cov_rep %>%select(-Group) %>% mutate(Batch="Zhu") ,
         Cov %>% rename(ID = Sequencing_ID)  %>% mutate(Batch = ifelse( Status == 1,"GUTS", "DMP")) ) %>% 
    rbind( . , mutate(Cov_r2, Batch="DMP") ) %>% filter(ID %in% Data_integrated2$ID ) -> Cov_integrated2
  Beta_taxonomy(NA, Data_c = left_join(Data_integrated2,Cov_integrated2) , Data_integrated2, Meta="Batch")
}
Beta_taxonomy(Species_data, Data_c, Data, Meta="Reads_clean")
Do_beta_validation(Data)

##1.4 Alpha diversity
###Compute alpha diversity at the species level ; need to do this in non-transformed data
Alpha_taxonomy = function(Data, Data_c, Taxonomy = "s", Meta = Covariates_all, Remove_names = c("ID", "UNKNOWN"), Cov_name=c("Reads_clean"), Check_batch=T ){
  paste( c(Taxonomy, "__"), collapse="" ) -> Ta
  colnames(select(Data, -one_of(Remove_names))) -> Taxa
  Select_taxonomy_level(Taxa, Ta) -> Taxa
  Data %>% select(Taxa)  -> Data_taxonomy
  diversity(Data_taxonomy, "shannon" ) -> Diversity
  Data_c %>% select( one_of(c("ID", "Status", "Batch", Cov_name))) %>% mutate(Diversity = Diversity) -> Data_div

  Formula = as.formula(paste0("Diversity ~ Status +", paste(Cov_name, collapse=" + ") ))
  summary(lm(Formula,Data_div)) %>% print()
  
  ggplot(Data_div, aes(x=as.factor(Status), y=Diversity, col=as.factor(Status)))  + geom_violin() +
  theme_bw() + coord_flip() + geom_jitter() +
    stat_summary(fun = "median",geom = "crossbar", aes(color = as.factor(Status))) +
    stat_summary(fun = "mean", geom = "point", color= "black") + ggtitle("Disease differences") -> Plot1
  
  if (Check_batch ==T){
  print("Running for batches, if no batches, this will fail")
  summary(lm( as.formula(paste0("Diversity ~ Batch +",paste(Cov_name, collapse=" + ") ))  ,filter(Data_div, Status == 1 ) )) %>% print()
  ggplot( filter(Data_div, Status==1) , aes(x=as.factor(Batch), y=Diversity, col=as.factor(Batch)))  + geom_violin() +
    theme_bw() + coord_flip() + geom_jitter() +
    stat_summary(fun = "median",geom = "crossbar", aes(color = as.factor(Batch))) +
    stat_summary(fun = "mean", geom = "point", color= "black")  + ggtitle("Batch differences (GUTS)")-> Plot2
  
  return(list(Plot1, Plot2))
  }
  else ( return(Plot1) )
  
}
Alpha_taxonomy(Data, Data_c, "s") ; Alpha_taxonomy(Data, Data_c, "g")
#Validation
Alpha_taxonomy(Taxonomic_replication2, left_join(Taxonomic_replication2, Cov_rep) , "g", Cov_name=c("Age", "Sex", "BMI"), Check_batch = F) #This is the analysis ran in the original paper, not replicated if genera and no covariates



#2. Biomarker discovery

##2.1 Use coda core to predict disease
library(codacore)
library(tensorflow)
library(pROC)

predict.codacore = function(object, newx, asLogits=TRUE, numLogRatios=NA, ...) {
  # Throw an error if zeros are present
  if (any(newx == 0)) {
    if (object$logRatioType == 'A') {
      warning("The data contain zeros. An epsilon is used to prevent divide-by-zero errors.")
    } else if (object$logRatioType == 'B') {
      stop("The data contain zeros. Balances cannot be used in this case.")
    }
  }
  
  x = .prepx(newx)
  yHat = rep(0, nrow(x))
  
  if (is.na(numLogRatios)) {
    numLogRatios = length(object$ensemble)
  }
  
  for (i in 1:numLogRatios) {
    cdbl = object$ensemble[[i]]
    yHat = yHat + object$shrinkage * predict.CoDaBaseLearner(cdbl, x)
  }
  
  if (object$objective == 'binary classification') {
    if (asLogits) {
      return(yHat)
    } else {
      return(1 / (1 + exp(-yHat)))
    }
  } else if (object$objective == 'regression') {
    return(yHat * object$yScale + object$yMean)
  }
}
.prepx = function(x) {
  if (class(x)[1] == 'tbl_df') {x = as.data.frame(x)}
  if (class(x)[1] == 'data.frame') {x = as.matrix(x)}
  if (is.integer(x)) {x = x * 1.0}
  
  # If the data is un-normalized (e.g. raw counts),
  # we normalize it to ensure our learning rate is well calibrated
  x = x / rowSums(x)
  return(x)
}
predict.CoDaBaseLearner = function(cdbl, x, asLogits=TRUE) {
  logRatio = computeLogRatio.CoDaBaseLearner(cdbl, x)
  eta = cdbl$slope * logRatio + cdbl$intercept
  if (asLogits) {
    return(eta)
  } else {
    if (cdbl$objective == 'regression') {
      stop("Logits argument should only be used for classification, not regression.")
    }
    return(1 / (1 + exp(-eta)))
  }
}
computeLogRatio.CoDaBaseLearner = function(cdbl, x) {
  
  if (!any(cdbl$hard$numerator) | !any(cdbl$hard$denominator)) {
    logRatio = rowSums(x * 0)
  } else { # we have a bona fide log-ratio
    if (cdbl$logRatioType == 'A') {
      epsilon = cdbl$optParams$epsilonA
      pvePart = rowSums(x[, cdbl$hard$numerator, drop=FALSE]) # drop=FALSE to keep as matrix
      nvePart = rowSums(x[, cdbl$hard$denominator, drop=FALSE])
      logRatio = log(pvePart + epsilon) - log(nvePart + epsilon)
    } else if (cdbl$logRatioType == 'B') {
      pvePart = rowMeans(log(x[, cdbl$hard$numerator, drop=FALSE])) # drop=FALSE to keep as matrix
      nvePart = rowMeans(log(x[, cdbl$hard$denominator, drop=FALSE]))
      logRatio = pvePart - nvePart
    }
  }
  
  return(logRatio)
}

Run_balance_analysis = function(Balance_input, Info = Data_c, lambda=1, add_pseudo = T, Strategy="Split"){
  set.seed(50)
  Labels = Data_c$Status
  Balance_input %>% select(-ID) -> Balance_input2
  if (add_pseudo == T){
    Pseudo = apply(Balance_input2, 2, Pseudocount)
    Balance_input2 = as.data.frame(Balance_input2) + Pseudo
  } else if (add_pseudo == "Imput"){
    Balance_input2 = zCompositions::cmultRepl(as.data.frame(Balance_input2) )
  }
  tf$random$set_seed(0)
  if (Strategy == "Split"){
    #Split dataset
    trainIndex <- sample(1:nrow( Balance_input2), 0.8 * nrow(Balance_input2))
    #Train set
    xTrain <- Balance_input2[trainIndex,]
    yTrain <- as.factor(Labels[trainIndex])
    #Test set
    xTest <- Balance_input2[-trainIndex,]
    yTest <-as.factor(Labels[-trainIndex])
  } else if (Strategy  == "Batch"){
    xTrain1 = Balance_input2 %>% mutate(ID = Balance_input$ID) %>% filter(ID %in% filter(Data_c, Batch=="GUTS_2")$ID ) %>% select(-one_of("UNKNOWN","ID"))
    xTest1 = Balance_input2 %>% mutate(ID = Balance_input$ID) %>% filter(ID %in% filter(Data_c, Batch=="GUTS_1")$ID ) %>% select(-one_of("UNKNOWN","ID"))
    
    IDs_controls =  filter(Data_c, Status == 0)$ID
    trainIndex <- sample(size = dim(xTest1)[1], x = 1:length(IDs_controls) , replace = F )
    IDTrain2 = IDs_controls[-trainIndex]
    IDTest2 = IDs_controls[trainIndex]
    
    xTrain2 =  select(filter( mutate(Balance_input2, ID=Balance_input$ID), ID %in% IDTrain2 ), -c(ID, UNKNOWN))
    xTest2 = select(filter( mutate(Balance_input2, ID=Balance_input$ID),  ID %in% IDTest2 ), -c(ID, UNKNOWN))
    
    rbind(xTrain1, xTrain2) -> xTrain
    rbind(xTest1, xTest2) -> xTest
    
    c(rep("Case", dim(xTrain1)[1]) , rep("Control", dim(xTrain2)[1]) ) -> yTrain
    c(rep("Case", dim(xTest1)[1]) , rep("Control", dim(xTest2)[1]) ) -> yTest
  }
  yTest = as.factor(yTest)
  yTrain = as.factor(yTrain)
  model=codacore( xTrain ,  yTrain , logRatioType = 'balances', lambda = lambda) #offset in logit space

  print("Train results")
  codacoreAUC = model$ensemble[[1]]$AUC
  cat("Train set AUC (model 1) =", codacoreAUC, "\n")
  if (length(model$ensemble) > 1){ 
    codacoreAUC = model$ensemble[[2]]$AUC
    cat("Train set AUC (model 2) =", codacoreAUC, "\n")
  }
  ##Test
  yHat <- predict.codacore(model, newx = xTest,asLogits = F)
  yHat2 <- predict.codacore(model, newx = xTest, numLogRatios = 1, asLogits=F)
  
  print("Test results ratio 2")
  testAUC <- pROC::auc(pROC::roc(yTest, yHat, quiet=T))
  cat("Test set AUC =", testAUC, "\n")
  print("Test results ratio 1")
  testAUC <- pROC::auc(pROC::roc(yTest, yHat2, quiet=T))
  cat("Test set AUC =", testAUC, "\n")
  
  
  failure <- yHat < 0.5 ; success <- yHat >= 0.5
  yHat[failure] <- levels(yTest)[1] ; yHat[success] <- levels(yTest)[2]
  cat("Classification accuracy on test set, ratio 2 =", round(mean(yHat == yTest), 2))
  
  failure <- yHat2 < 0.5 ; success <- yHat2 >= 0.5
  yHat2[failure] <- levels(yTest)[1] ; yHat2[success] <- levels(yTest)[2]
  cat("Classification accuracy on test set, ratio1 =", round(mean(yHat2 == yTest), 2))
  
  
  
  plot(model)
  plotROC(model)
  Numerators = colnames(xTrain)[getNumeratorParts(model, 1)] 
  Denominators = colnames(xTrain)[getDenominatorParts(model, 1)]
 
  if (length(model$ensemble) > 1){ 
    Numerators2 = colnames(xTrain)[getNumeratorParts(model, 2)] 
    Denominators2 = colnames(xTrain)[getDenominatorParts(model, 2)]
    return(list(Numerators, Denominators, Numerators2, Denominators2, Balance_input[trainIndex,]$ID, model))
    
  } else { return(list(Numerators, Denominators,  Balance_input2[trainIndex,]$ID ), model ) }
}   

#Get only prevalent bacteria for the model
Prevalence %>% filter(N_case > 20 & N_control > 20) -> Keep
Run_balance_analysis( select(Data, c("ID", Keep$Bug)) ) -> Balances
#Balance using batches as different train/test sets
#Run_balance_analysis( select(Data, c("ID", Keep$Bug)), add_pseudo = "Imput" , Strategy="Batch"  ) -> Balances2
Balances[[4]] %>% sapply(function(x){ str_split(x, "\\|")[[1]] -> y ; y[length(y)]  } ) %>% as_vector() %>% as.vector()

Data %>% select(Balances[[1]]) %>% apply(1, function(x){ y = x + 1; Geom_mean(y) } ) -> Numerator ; Data %>% select(Balances[[3]]) %>% apply(1, function(x){ y = x + 1; Geom_mean(y) } ) -> Numerator2
Data %>% select(Balances[[2]]) %>% apply(1, function(x){ y = x + 1; Geom_mean(y) } ) -> Denominator ; Data %>% select(Balances[[4]]) %>% apply(1, function(x){ y = x + 1; Geom_mean(y) } ) -> Denominator2
Data %>% mutate( Balance = log10(Numerator/Denominator), Balance2=log10(Numerator2/Denominator2), Status = Data_c$Status, Read_number=Data_c$Reads_clean ) ->Data_b
Data_b %>% ggplot(aes(x=as.factor(Status), y= Balance, col=as.factor(Status)))  + geom_violin() + theme_bw() + coord_flip() + geom_jitter() +
  stat_summary(fun = "median",geom = "crossbar", aes(color = as.factor(Status))) + stat_summary(fun = "mean", geom = "point", color= "black")
summary(lm(Balance ~ Status + Read_number , filter( Data_b, ! abs(Balance) == Inf )))
Data_b %>% ggplot(aes(x=as.factor(Status), y= Balance2, col=as.factor(Status)))  + geom_violin() + theme_bw() + coord_flip() + geom_jitter() +
  stat_summary(fun = "median",geom = "crossbar", aes(color = as.factor(Status))) + stat_summary(fun = "mean", geom = "point", color= "black")
summary(lm(Balance2 ~ Status + Read_number , filter( Data_b, ! abs(Balance2) == Inf )))

Data_b %>% select(ID, Balance, Balance2, Read_number) %>% write_tsv("Results/Balances.tsv")

#Replication
        
Comp_Prevalence( select( Taxonomic_replication2  , -ID) ) -> Prevalence_replication
Prevalence_replication %>% filter(Bug %in% unlist(Balances[1:4]) ) %>% group_by(Prevalence>0.1) %>% summarise(n()) #19/19 found

Taxonomic_replication2 %>% select( one_of(Balances[[1]])) %>% apply(1, function(x){ y = x + 1; Geom_mean(y) } ) -> Numerator ; Taxonomic_replication2 %>% select(one_of(Balances[[3]])) %>% apply(1, function(x){ y = x + 1; Geom_mean(y) } ) -> Numerator2
Taxonomic_replication2 %>% select(one_of(Balances[[2]])) %>% apply(1, function(x){ y = x + 1; Geom_mean(y) } ) -> Denominator ; Taxonomic_replication2 %>% select(one_of(Balances[[4]])) %>% apply(1, function(x){ y = x + 1; Geom_mean(y) } ) -> Denominator2
Taxonomic_replication2 %>% select(ID) %>% mutate(Balance = log10(Numerator/Denominator), Balance2=log10(Numerator2/Denominator2)) %>% left_join(. , Cov_rep) -> Data_b_replication
Data_b_replication %>% ggplot(aes(x=as.factor(Status), y= Balance, col=as.factor(Status)))  + geom_violin() + theme_bw() + coord_flip() + geom_jitter() +
  stat_summary(fun = "median",geom = "crossbar", aes(color = as.factor(Status))) + stat_summary(fun = "mean", geom = "point", color= "black")
summary(lm(Balance ~ Status + Age + Sex , filter( Data_b_replication, ! abs(Balance) == Inf )))
Data_b_replication %>% ggplot(aes(x=as.factor(Status), y= Balance2, col=as.factor(Status)))  + geom_violin() + theme_bw() + coord_flip() + geom_jitter() +
  stat_summary(fun = "median",geom = "crossbar", aes(color = as.factor(Status))) + stat_summary(fun = "mean", geom = "point", color= "black")
summary(lm(Balance2 ~ Status + Age + Sex  , filter( Data_b_replication, ! abs(Balance2) == Inf )))

#Replication with DMP
Comp_Prevalence( select( Taxonomic_replication_DMP2  , -ID) ) -> Prevalence_replication2
Prevalence_replication2 %>% filter(Bug %in% unlist(Balances[1:4]) ) %>% group_by(Prevalence>0.1) %>% summarise(n()) #19/19 found

#Compute balances in new data
Taxonomic_replication_DMP2 %>% select( one_of(Balances[[1]])) %>% apply(1, function(x){ y = x + 1; Geom_mean(y) } ) -> Numerator ; Taxonomic_replication_DMP2 %>% select(one_of(Balances[[3]])) %>% apply(1, function(x){ y = x + 1; Geom_mean(y) } ) -> Numerator2
Taxonomic_replication_DMP2 %>% select(one_of(Balances[[2]])) %>% apply(1, function(x){ y = x + 1; Geom_mean(y) } ) -> Denominator ; Taxonomic_replication_DMP2 %>% select(one_of(Balances[[4]])) %>% apply(1, function(x){ y = x + 1; Geom_mean(y) } ) -> Denominator2
Taxonomic_replication_DMP2 %>% select(ID) %>% mutate(Balance = log10(Numerator/Denominator), Balance2=log10(Numerator2/Denominator2)) %>% mutate(Status=2) -> Data_b_replicationDMP

#Merge previous data and new data
Data_b %>% select(ID, Balance, Balance2, Status) %>% rbind(. , Data_b_replicationDMP) -> Data_b_replicationDMP

#lm(Balance2 ~ as.factor(Status), filter( Data_b_replicationDMP, ! abs(Balance2) == Inf ),contrasts = matrix(c(1, -1/2, -1/2, 0, .5, -.5), ncol = 2)
# ) %>% summary()
#Using both train and test data
lm(Balance ~ as.factor(Status), filter( Data_b_replicationDMP, ! abs(Balance2) == Inf ) ) %>% summary()



Data_b_replicationDMP %>% ggplot(aes(x=as.factor(Status), y= Balance2, col=as.factor(Status)))  + geom_violin() + theme_bw() + coord_flip() + geom_jitter() +
  stat_summary(fun = "median",geom = "crossbar", aes(color = as.factor(Status))) + stat_summary(fun = "mean", geom = "point", color= "black") + scale_color_manual(values=met.brewer("Manet", 5)[c(1,5, 3)] )

Data_b_replicationDMP %>% mutate( Cohort = ifelse(Status == 0 | Status== 2, "DMP", "GUTS" ), Status = ifelse(Status==2, 1, Status)  ) %>% rbind(. , select(Data_b_replication, c(ID, Balance, Balance2, Status)) %>% mutate(Cohort = "Zhu")) -> Data_b_all

Data_b_all %>% lmerTest::lmer(Balance2 ~   Status + (1|Cohort) , . ) %>% summary()
Data_b_all %>% lmerTest::lmer(Balance ~  Status + (1|Cohort) , . ) %>% summary()

Data_b_all%>% filter(! is.infinite(Balance2) | is.na(Balance2) ) %>% mutate(Status= as.factor(Status)) %>% ggplot(aes(y=Balance2, x= Status, col=Cohort ))  + geom_violin() + theme_bw()  + ggforce::geom_sina() +
  stat_summary(fun = "mean", geom = "point", color= "black") + scale_color_manual(values=met.brewer("Manet", 5)[c(1,5, 3)] ) + 
  stat_summary(aes(group = factor(paste0(Status, Cohort)), color=Cohort) ,fun = "median",geom = "crossbar")
#Removing training (saved in Balances[[5]])
Data_b_all %>% filter(! ID %in% Balances[[5]]) %>% lmerTest::lmer(Balance2 ~  Cohort + Status + (1|Cohort) , . ) %>% summary()
Data_b_all%>% filter(! ID %in% Balances[[5]]) %>% filter(! is.infinite(Balance2) | is.na(Balance2) ) %>% mutate(Status= as.factor(Status)) %>% ggplot(aes(y=Balance2, x= Status, col=Cohort ))  + geom_violin() + theme_bw()  + ggforce::geom_sina() +
  stat_summary(fun = "mean", geom = "point", color= "black") + scale_color_manual(values=met.brewer("Manet", 5)[c(1,5, 3)] ) + 
  stat_summary(aes(group = factor(paste0(Status, Cohort)), color=Cohort) ,fun = "median",geom = "crossbar")
Data_b_all%>% filter(! ID %in% Balances[[5]]) %>% filter(Cohort == "DMP") %>% lm(Balance2 ~ as.factor(Status), . ) %>% summary()
Data_b_all%>% filter(! ID %in% Balances[[5]]) %>% filter(Cohort == "DMP") %>% lm(Balance ~ as.factor(Status), . ) %>% summary()

#Conclusions:  Without training DMP alone replicates. Mixed-model with all three cohorts is also significant.



#2.2 Association analysis: linear model on CLR-transformed data 
##Inverse-rank normal transformation, also known as INT --> enforce normality of the data. Applied after sample-specific normalization factors applied in CLR
Inverse_rank_normal = function(Measurement){
  qnorm((rank(Measurement,na.last="keep")-0.5)/sum(!is.na(Measurement)))
}

Association_analysis = function(Prevalence, DF, Data_c, FILTER=20, Meta=c("Read_number"), X = "Status") {
  Total_results = tibble()
  Prevalence %>% filter(N_case > FILTER & N_control > FILTER) %>% filter(! Bug == "UNKNOWN") -> To_Test
  for (Bug in To_Test$Bug){
    if (! Bug %in% colnames(DF) ){ next }
    DF %>% select( one_of(Bug) ) %>% as_vector() -> vector_Bug
    #IF we want to normalize it, apply function here
    normalized_Bug = Inverse_rank_normal(vector_Bug)
    #invers rank normal transf?
    Data_c %>% select(one_of(c("ID", X, "Batch", Meta))) %>% mutate(B =  vector_Bug, B_n = normalized_Bug) -> Model_input
    Formula =   paste0(c("B ~ ", paste0(c(X, Meta) , collapse = "+")), collapse="")
    Formula2 = paste0(c("B_n ~ ", paste0(c(X, Meta) , collapse = "+")), collapse="")
    
    lm(Formula, Model_input  ) -> Model_out
    lm(Formula2, Model_input  ) -> Model_out_n
    #If not normalized
    Normality = shapiro.test(Model_out$residuals)$p.value
    as.data.frame(summary(Model_out)$coefficients)[X,] %>% as_tibble() %>%
      mutate(Bug = Bug, Shapiro_p = Normality, .before=1) -> results
    #Normalized
    Normality = shapiro.test(Model_out_n$residuals)$p.value
    as.data.frame(summary(Model_out_n)$coefficients)[X,] %>% as_tibble() %>%
      mutate(Bug = Bug, Shapiro_p = Normality, .before=1) -> Normalized_results
    
    left_join(results, Normalized_results, by= "Bug", suffix = c("","_norm")) -> results
    rbind(Total_results, results) -> Total_results
}
  return(Total_results)
}
Association_analysis_all = function(Prevalence, DF, Data_c, FILTER=20, Meta=c("Read_number")) {
  Total_results = tibble()
  Prevalence %>% filter(N_case > FILTER & N_control > FILTER ) %>% filter(! Bug == "UNKNOWN") -> To_Test
  for (Bug in To_Test$Bug){
    if (! Bug %in% colnames(DF) ){ next }
    DF %>% select( one_of(Bug) ) %>% as_vector() -> vector_Bug
    #IF we want to normalize it, apply function here
    normalized_Bug = Inverse_rank_normal(vector_Bug)
    #invers rank normal transf?
    Data_c %>% select(one_of(c("ID", "Status", "Batch", Meta))) %>% mutate(B =  vector_Bug, B_n = normalized_Bug) -> Model_input
    Formula2 = paste(c("B_n ~ Status", Meta, "(1|Batch)"), collapse="+")
    
    #lm(Formula, Model_input  ) -> Model_out
    lmerTest::lmer(Formula2, Model_input  ) -> Model_out_n
    #If not normalized
    #Normality = shapiro.test(Model_out$residuals)$p.value
    #as.data.frame(summary(Model_out)$coefficients)["Status",] %>% as_tibble() %>%
    #  mutate(Bug = Bug, Shapiro_p = Normality, .before=1) -> results
    #Normalized
    #Normality = shapiro.test(Model_out_n$residuals)$p.value
    as.data.frame(summary(Model_out_n)$coefficients)["Status",] %>% as_tibble() %>%
      mutate(Bug = Bug, .before=1) %>% rename(`Pr(>|t|)_norm`  = `Pr(>|t|)`) -> Normalized_results
    
    Normalized_results -> results
    rbind(Total_results, results) -> Total_results
  }
  return(Total_results)
}  

Association_analysis_readNumber = function(Prevalence, DF, Data_c, FILTER=20, Meta = "Reads_clean" ){
  Total_results = tibble()
  Prevalence %>% filter(N_case > FILTER & N_control > FILTER) %>% filter(! Bug == "UNKNOWN") -> To_Test
  for (Bug in To_Test$Bug){
    if (! Bug %in% colnames(DF) ){ next }
    DF %>% select( one_of(Bug) ) %>% as_vector() -> vector_Bug
    #IF we want to normalize it, apply function here
    normalized_Bug = Inverse_rank_normal(vector_Bug)
    #invers rank normal transf?
    Data_c %>% select(one_of(c("ID", "Status", "Batch", Meta))) %>% mutate(B =  vector_Bug, B_n = normalized_Bug) -> Model_input
    Formula = paste( c("B ~ ", Meta), collapse = "+")
    Formula2= paste( c("B_n ~ ",Meta),collapse = "+")
    
    lm(Formula, Model_input  ) -> Model_out
    lm(Formula2, Model_input  ) -> Model_out_n
    #If not normalized
    Normality = shapiro.test(Model_out$residuals)$p.value
    as.data.frame(summary(Model_out)$coefficients)[Meta,] %>% as_tibble() %>%
      mutate(Bug = Bug, Shapiro_p = Normality, .before=1) -> results
    #Normalized
    Normality = shapiro.test(Model_out_n$residuals)$p.value
    as.data.frame(summary(Model_out_n)$coefficients)[Meta,] %>% as_tibble() %>%
      mutate(Bug = Bug, Shapiro_p = Normality, .before=1) -> Normalized_results
    
    left_join(results, Normalized_results, by= "Bug", suffix = c("","_norm")) -> results
    rbind(Total_results, results) -> Total_results
  }
  return(Total_results)
}  
Run_Associations = function( Prevalence, DF, Data_c, FILTER_n=20, Meta_n=c("Reads_clean"), Function =  Association_analysis, Run_permutations = F, Permutation_number=100, Independent="Status", Permute_p="Pr(>|t|)_norm" ){
  #Per taxonomic level, run Function and perform 100 permutations to obtain taxonomic-level-specific FDR
  Summary_statistics = tibble()
  for (Taxonomic_level in c("p", "c", "o", "f", "g", "s")){
    print( paste0("Subsetting features of taxonomic level: ", Taxonomic_level) )
    paste(c(Taxonomic_level, "__"), collapse="" ) -> Ta
  
    colnames(select(DF, - one_of(c("ID" )))) -> Taxa
    Select_taxonomy_level(Taxa, Ta) -> Taxa
    
    DF %>% select( one_of(c("ID", Taxa) )) -> DF_taxonomyX

    print("Calling association function" )
    Summary_statistics_taxonomyX = Function( Prevalence, DF_taxonomyX, Data_c, FILTER=FILTER_n, Meta=Meta_n, X=Independent )
    Summary_statistics_taxonomyX %>% mutate(FDR_bh =p.adjust(`Pr(>|t|)`, "fdr") ,FDR_bh_norm = p.adjust(`Pr(>|t|)_norm`, "fdr") ) -> Summary_statistics_taxonomyX
    print("Association succesful" )
    
    if (Run_permutations == T){
        print("Running permutations" )
        Null_distribution = c()
        for (Permutation_n in seq(Permutation_number) ){
              Data_perm = Data_c
              Data_perm[Independent] = sample( as_vector(Data_c[Independent]), replace = F,size = dim(Data_c)[1] )
              Association_analysis(Prevalence, DF_taxonomyX, Data_perm, FILTER= FILTER_n, Meta=Meta_n, X=Independent) -> Association_results_p
              Null_distribution = c(Null_distribution, as_vector(Association_results_p[Permute_p]) )
        }

      #sapply( Summary_statistics_taxonomyX$`Pr(>|t|)_norm`, function(x){ sum( Null_distribution <= x  )/Permutation_number  } ) -> FDR_perm
      Summary_statistics_taxonomyX["P"] = as_vector(Summary_statistics_taxonomyX[Permute_p])
      FDR_perm = Compute_FDR_null(Summary_statistics_taxonomyX , Null_distribution, fdr_threshold=0.05  ) 
      FDR_perm[FDR_perm > 1] = 1
      Summary_statistics_taxonomyX %>% mutate(FDR_permutation = FDR_perm ) -> Summary_statistics_taxonomyX
      print("Permutations succesful" )
    }
    Summary_statistics = rbind ( Summary_statistics, mutate(Summary_statistics_taxonomyX, Taxonomic_level = Taxonomic_level )   )
  }
  return(Summary_statistics)
}
Compute_FDR_null =  function(DF, Permutation_results,fdr_threshold = 0.05 ){
  # calculate FDR for each observed P-value
  fdr_values <- sapply(DF$P, function(p) {
    # count number of false positives (i.e., permuted P < observed P)
    false_positives <- sum(Permutation_results <= p)
    # calculate FDR
    fdr <- false_positives / sum(Permutation_results <= fdr_threshold)
    # ensure FDR is between 0 and 1
    fdr <- ifelse(fdr < 0, 0, ifelse(fdr > 1, 1, fdr))
    return(fdr)
  })
  return(fdr_values)
}  

##2.2.1 Check the effect of read number  per batch
Association_results_rn = tibble()
for (batch in unique(Data_c$Batch )){
  print(paste0("Running in batch:", batch))
  Data_c %>% filter(Batch == batch) -> Subset_c
  All_taxonomy %>% filter(ID %in% Subset_c$ID) -> Subset
  Run_Associations(Prevalence, Subset, Subset_c, FILTER=10, Function=Association_analysis_readNumber) %>% mutate(Batch=batch) %>% rbind(Association_results_rn, . ) -> Association_results_rn
  #Association_analysis_readNumber(Prevalence, Subset, Subset_c, FILTER=10) %>% mutate(Batch=batch) %>% rbind(Association_results_rn, . ) -> Association_results_rn
}
Association_results_rn_meta = tibble()
#Batch meta analysis
for ( Bug_name in unique(Association_results_rn$Bug) ){
  Association_results_rn %>% filter(Bug == Bug_name) -> For_meta
  meta::metagen(TE=Estimate, seTE=`Std. Error`, data=For_meta, comb.fixed = T, comb.random= T ) -> meta_value
  sub_result = tibble(Bug = Bug_name, MetaP=meta_value$pval.fixed, Treatment_effect=meta_value$TE.fixed,
                      SE_meta=meta_value$seTE.fixed,Heterogeneity_Q=meta_value$Q,Heterogeneity_Pvalue=meta_value$pval.Q,  
                      Meta_random_P=meta_value$pval.random, Treatment_effect_random=meta_value$TE.random, SE_meta_random=meta_value$seTE.random, Taxonomic_level=unique(For_meta$Taxonomic_level)) 
  rbind(Association_results_rn_meta, sub_result) -> Association_results_rn_meta
}
Association_results_rn_meta %>% mutate(FDR = p.adjust(Meta_random_P, "fdr") ) %>% ggplot(aes(x=Treatment_effect_random, y = -log10(Meta_random_P), col=FDR<0.05 )) + 
  geom_point() + theme_bw() + geom_hline(yintercept = -log10(0.05)) + 
  ggtitle("Meta analysis of read count effect\non taxa abundance per batch") + facet_wrap(~Taxonomic_level)
#At the species level, no FDR significant hits

##2.2.2 Run association analysis between disease and gut taxa

Run_Associations(Prevalence, All_taxonomy, Data_c, FILTER=10, Function=Association_analysis, Run_permutations = F) -> Association_results
Association_results %>% filter(Taxonomic_level == "s") %>% ggplot(aes(x=-log10(FDR_bh), y = -log10(FDR_bh_norm) )) + geom_point() + theme_bw() + geom_abline() + geom_vline(xintercept = -log10(0.05) ) + geom_hline(yintercept = -log10(0.05) )


###2.2.3 Add info from balances
Balances_results = tibble(Bug = unique(c(Balances[[1]], Balances[[3]]) ), Direction_balance = "Positive")
Balances_results = rbind(Balances_results, tibble(Bug = unique(c(Balances[[2]], Balances[[4]] ) ), Direction_balance = "Negative") )
left_join(Association_results,Balances_results) -> Association_results
Association_results %>% filter(Taxonomic_level == "s") %>% ggplot(aes(x=-log10(FDR_bh), y = -log10(FDR_bh_norm), col=Direction_balance )) + geom_point() + theme_bw() + geom_abline() + geom_vline(xintercept = -log10(0.05) ) + geom_hline(yintercept = -log10(0.05) )


###2.2.4 Plot
Association_results %>% ggplot(aes(x=Estimate_norm, y= -log10(`Pr(>|t|)_norm`), col=FDR_permutation<0.05, shape= is.na(Direction_balance) )) +
  geom_hline( yintercept = -log10(0.05), color="red" ) + geom_point() + theme_bw() + facet_wrap(~Taxonomic_level) +
  geom_label_repel(data = Association_results %>% filter(FDR_permutation<0.05) %>% unique(),  aes(label = Bug), size=1.5, color = "black")

#Run associaton with all data
Comp_Prevalence( select( filter(Taxonomic_replication_DMP2 )  , -ID) ) -> Prevalence_replication2
Comp_Prevalence( select( filter(Taxonomic_replication2 )  , -ID) ) -> Prevalence_replication
filter(Prevalence, Bug %in% Prevalence_replication2$Bug & Bug %in% Prevalence_replication$Bug)   -> To_keep

rbind(rbind( dplyr::select(All_taxonomy, one_of(c("ID", To_keep$Bug))) ,  dplyr::select(Data_CLR_replication1[[1]], one_of(c("ID", To_keep$Bug) ) )), dplyr::select(Data_CLR_replication2[[1]], one_of(c("ID", To_keep$Bug) ) )) -> taxa_test
Comp_Prevalence( select( taxa_test, -ID) ) -> Prevalence_all
rbind(rbind(Covariates_all %>% rename(ID = Sequencing_ID) %>% mutate(Batch = ifelse(Status==1, "GUTS", "DMP"  )) , Cov_rep %>% select(-Group) %>% mutate(Batch="Zhu") ), Cov_rep2 %>% mutate(Batch="DMP")) %>% arrange(ID) -> Covariates_test
Association_results_all = Run_Associations(To_keep, taxa_test, left_join(taxa_test, Covariates_test), FILTER_n=20, Meta_n=c("Reads_clean"), Function =  Association_analysis_all, Run_permutations = T, Permutation_number=100) 
Association_results_all %>% ggplot(aes(x=Estimate, y= -log10(`Pr(>|t|)_norm`), col=FDR_permutation<0.05 )) +
  geom_hline( yintercept = -log10(0.05), color="red" ) + geom_point() + theme_bw() + facet_wrap(~Taxonomic_level) +
  geom_label_repel(data = Association_results_all %>% filter(FDR_permutation<0.05) %>% unique(),  aes(label = Bug), size=1.5, color = "black")
Association_results_all %>% filter(FDR_permutation < 0.05) %>% View()



#
##2.2.3 Run association analysis between disease and gut taxa using Aldex2
Aldex_results = tibble()
for (Taxonomic_level in c("p", "c", "o", "f", "g", "s")){
  
  Prevalence %>% filter(N_case > 10 & N_control > 10) -> To_run 
  
  print( paste0("Subsetting features of taxonomic level: ", Taxonomic_level) )
  paste(c(Taxonomic_level, "__"), collapse="" ) -> Ta
  
  colnames(dplyr::select(Data, - one_of(c("ID" )))) -> Taxa
  Select_taxonomy_level(Taxa, Ta) -> Taxa
  
  Data %>% dplyr::select( one_of(c("ID", Taxa) )) %>% select_if( colnames(.) %in% c("ID", To_run$Bug)) -> Data_taxonomyX
  Data_taxonomyX %>% as.data.frame() %>% column_to_rownames("ID") %>% t() -> Input_aldex
  print("Calling association function" )
  x <- ALDEx2::aldex.clr(Input_aldex, as.character(as.factor(Data_c$Status)), mc.samples=160, denom="all", verbose=FALSE)
  x.tt <- ALDEx2::aldex.ttest(x, paired.test=FALSE, verbose=FALSE)
  x.effect <- ALDEx2::aldex.effect(clr = x, CI=T, verbose=FALSE )

  x.all <- data.frame(x.tt,x.effect)
  x.all %>% rownames_to_column("Bug") %>% as_tibble() %>% mutate(Taxonomic_level = Taxonomic_level) -> Result_taxa
  print("Association succesful" )
  rbind(Aldex_results, Result_taxa) -> Aldex_results
}

Aldex_results %>% left_join(. ,Balances_results) %>% ggplot(. , aes(x=effect, y= -log10(`wi.ep`), col=`wi.eBH`<0.05, shape= is.na(Direction_balance) )) +
  geom_hline( yintercept = -log10(0.05), color="red" ) + geom_point() + theme_bw() + facet_wrap(~Taxonomic_level) +
  geom_label_repel(data = . %>% filter(wi.eBH<0.05) %>% unique(),  aes(label = Bug), size=1.5, color = "black")

left_join(Association_results, Aldex_results, by = c("Bug", "Taxonomic_level")) %>% left_join(., Prevalence) -> Association_results_all
Association_results_all %>%
    write_tsv( ., "Results/Summary_statistics_TaxonomyAnalysis.tsv")  

Association_results_all %>% filter(Taxonomic_level == "s") %>% ggplot(aes(x=-log10(FDR_bh), y = -log10(we.eBH ), col=Direction_balance )) + geom_point() + theme_bw() + geom_abline() + geom_vline(xintercept = -log10(0.05) ) + geom_hline(yintercept = -log10(0.05) ) -> All_comparisons1
Association_results_all %>% filter(Taxonomic_level == "s") %>% ggplot(aes(x=-log10(we.eBH), y = -log10(FDR_bh_norm), col=Direction_balance )) + geom_point() + theme_bw() + geom_abline() + geom_vline(xintercept = -log10(0.05) ) + geom_hline(yintercept = -log10(0.05) )-> All_comparisons2

Association_results_all %>% filter(Taxonomic_level == "s") %>% ggplot(aes(x=-log10(FDR_bh), y = -log10(wi.eBH ), col=Direction_balance )) + geom_point() + theme_bw() + geom_abline() + geom_vline(xintercept = -log10(0.05) ) + geom_hline(yintercept = -log10(0.05) )-> All_comparisons3
Association_results_all %>% filter(Taxonomic_level == "s") %>% ggplot(aes(x=-log10(wi.eBH ), y = -log10(FDR_bh_norm), col=Direction_balance )) + geom_point() + theme_bw() + geom_abline() + geom_vline(xintercept = -log10(0.05) ) + geom_hline(yintercept = -log10(0.05) )-> All_comparisons4

( All_comparisons1 | All_comparisons2 )  / (All_comparisons3 | All_comparisons4 ) 


###Replication
Taxonomic_replication2 %>% select(one_of(c("ID", Aldex_results$Bug))) %>% arrange(ID) -> replication_for_aldex
#Taxonomic_replication2 %>%arrange(ID) -> replication_for_aldex
Cov_rep %>% arrange(ID) %>% filter(ID %in% Taxonomic_replication2$ID) -> Cov_rep
Aldex_replication = tibble()
for (Taxonomic_level in c("p", "c", "o", "f", "g", "s")){
  print( paste0("Subsetting features of taxonomic level: ", Taxonomic_level) )
  paste(c(Taxonomic_level, "__"), collapse="" ) -> Ta
  
  colnames(dplyr::select(replication_for_aldex, - one_of(c("ID" )))) -> Taxa
  Select_taxonomy_level(Taxa, Ta) -> Taxa
  
  replication_for_aldex %>% dplyr::select( one_of(c("ID", Taxa) ))  -> Data_taxonomyX
  Data_taxonomyX %>% as.data.frame() %>% column_to_rownames("ID") %>% t() -> Input_aldex
  print("Calling association function" )
  x <- ALDEx2::aldex.clr(Input_aldex, as.character(as.factor(Cov_rep$Status)), mc.samples=160, denom="all", verbose=FALSE)
  x.tt <- ALDEx2::aldex.ttest(x, paired.test=FALSE, verbose=FALSE)
  x.effect <- ALDEx2::aldex.effect(clr = x, CI=T, verbose=FALSE )
  
  x.all <- data.frame(x.tt,x.effect)
  x.all %>% rownames_to_column("Bug") %>% as_tibble() %>% mutate(Taxonomic_level = Taxonomic_level) -> Result_taxa
  print("Association succesful" )
  rbind(Aldex_replication, Result_taxa) -> Aldex_replication
}

#17/17 taxa are available for replication
#1/15 reach P-value < 0.05 , Ruthenibacterium   same direction , almost significant: Actinobacteria same direction
Aldex_replication %>% filter(Bug %in%  filter(Aldex_results, wi.eBH < 0.05)$Bug ) %>% ggplot(. , aes(x=effect, y= -log10(`wi.ep`), col=`wi.eBH`<0.05 )) +
  geom_hline( yintercept = -log10(0.05), color="red" ) + geom_point() + theme_bw() + facet_wrap(~Taxonomic_level) +
  geom_label_repel(data = . %>% filter(wi.ep<0.05) %>% unique(),  aes(label = Bug), size=1.5, color = "black")
Aldex_replication %>% ggplot(. , aes(x=effect, y= -log10(`wi.ep`), col=`wi.eBH`<0.05 )) +
  geom_hline( yintercept = -log10(0.05), color="red" ) + geom_point() + theme_bw() + facet_wrap(~Taxonomic_level) +
  geom_label_repel(data = . %>% filter(`wi.eBH`<0.05) %>% unique(),  aes(label = Bug), size=1.5, color = "black")
Aldex_replication %>% arrange(wi.ep) %>% filter(Taxonomic_level == "s") %>% select(Bug, wi.ep, effect)

Aldex_replication_cov = tibble()
for (Taxonomic_level in c("p", "c", "o", "f", "g", "s")){
  print( paste0("Subsetting features of taxonomic level: ", Taxonomic_level) )
  paste(c(Taxonomic_level, "__"), collapse="" ) -> Ta
  
  colnames(dplyr::select(replication_for_aldex, - one_of(c("ID" )))) -> Taxa
  Select_taxonomy_level(Taxa, Ta) -> Taxa
  
  replication_for_aldex %>% dplyr::select( one_of(c("ID", Taxa) ))  -> Data_taxonomyX
  Data_taxonomyX %>% as.data.frame() %>% column_to_rownames("ID") %>% t() -> Input_aldex
  print("Calling association function" )
  mm <- model.matrix(~ as.character(Status) + Sex + Age + BMI, Cov_rep)
  x <- ALDEx2::aldex.clr(Input_aldex,mm, mc.samples=160, denom="all", verbose=FALSE)
  glm.test <- ALDEx2::aldex.glm(x, mm)
  glm.effect <- ALDEx2::aldex.glm.effect(x)
  
  x.all <- data.frame(glm.test,glm.effect)
  x.all %>% rownames_to_column("Bug") %>% as_tibble() %>% mutate(Taxonomic_level = Taxonomic_level) -> Result_taxa
  print("Association succesful" )
  rbind(Aldex_replication_cov, Result_taxa) -> Aldex_replication_cov
}
Aldex_replication_cov %>% ggplot(. , aes(x=`model.as.character.Status.1.Estimate`, y= -log10(`model.as.character.Status.1.Pr...t..`), col=`model.as.character.Status.1.Pr...t...BH`<0.05 )) +
  geom_hline( yintercept = -log10(0.05), color="red" ) + geom_point() + theme_bw() + facet_wrap(~Taxonomic_level) +
  geom_label_repel(data = . %>% filter( `model.as.character.Status.1.Pr...t...BH`<0.05) %>% unique(),  aes(label = Bug), size=1.5, color = "black")

#nothing significant when accounting for covariates
Aldex_replication_cov %>% filter(Bug %in%  filter(Aldex_results, wi.eBH < 0.05)$Bug ) %>% 
ggplot(. , aes(x=`model.as.character.Status.1.Estimate`, y= -log10(`model.as.character.Status.1.Pr...t..`), col=`model.as.character.Status.1.Pr...t...BH`<0.05 )) +
  geom_hline( yintercept = -log10(0.05), color="red" ) + geom_point() + theme_bw() + facet_wrap(~Taxonomic_level) +
  geom_label_repel(data = . %>% filter( `model.as.character.Status.1.Pr...t...BH`<0.05) %>% unique(),  aes(label = Bug), size=1.5, color = "black")

#Compare effects
left_join( Aldex_results %>% select(Bug, effect), Aldex_replication%>% select(Bug, effect), by="Bug", suffix=c("_Study", "_Zhu") ) %>% ggplot(aes(x=effect_Study, y=effect_Zhu)) + geom_point() + theme_bw() + geom_abline()+ geom_smooth(method = "lm")


##check in DMP2 
Taxonomic_replication_DMP2 %>% select(one_of(c("ID", Aldex_results$Bug))) %>% arrange(ID) -> replication_for_aldex
#Taxonomic_replication2 %>%arrange(ID) -> replication_for_aldex
rbind(Covariates_control, tibble(Sequencing_ID = replication_for_aldex$ID, Batch = 0, Age=0, Sex=0, BMI=0, Status=1))  -> Cov_rep
rbind( Data %>% filter(ID %in% Cov_rep$Sequencing_ID) %>% select(one_of(colnames(replication_for_aldex))) ,  replication_for_aldex) -> replication_for_aldex

Aldex_replication2 = tibble()
for (Taxonomic_level in c("p", "c", "o", "f", "g", "s")){
  print( paste0("Subsetting features of taxonomic level: ", Taxonomic_level) )
  paste(c(Taxonomic_level, "__"), collapse="" ) -> Ta
  
  colnames(dplyr::select(replication_for_aldex, - one_of(c("ID" )))) -> Taxa
  Select_taxonomy_level(Taxa, Ta) -> Taxa
  
  replication_for_aldex %>% dplyr::select( one_of(c("ID", Taxa) ))  -> Data_taxonomyX
  Data_taxonomyX %>% as.data.frame() %>% column_to_rownames("ID") %>% t() -> Input_aldex
  print("Calling association function" )
  x <- ALDEx2::aldex.clr(Input_aldex, as.character(as.factor(Cov_rep$Status)), mc.samples=160, denom="all", verbose=FALSE)
  x.tt <- ALDEx2::aldex.ttest(x, paired.test=FALSE, verbose=FALSE)
  x.effect <- ALDEx2::aldex.effect(clr = x, CI=T, verbose=FALSE )
  
  x.all <- data.frame(x.tt,x.effect)
  x.all %>% rownames_to_column("Bug") %>% as_tibble() %>% mutate(Taxonomic_level = Taxonomic_level) -> Result_taxa
  print("Association succesful" )
  rbind(Aldex_replication2, Result_taxa) -> Aldex_replication2
}

# 3 / 17 are replicated ( g__Odoribacter, g__Clostridiales, s__Odoribacter, f__Odoribacteraceae   )
Aldex_replication2 %>% filter(Bug %in%  filter(Aldex_results, wi.eBH < 0.05)$Bug ) %>% ggplot(. , aes(x=effect, y= -log10(`wi.ep`), col=`wi.eBH`<0.05 )) +
  geom_hline( yintercept = -log10(0.05), color="red" ) + geom_point() + theme_bw() + facet_wrap(~Taxonomic_level) +
  geom_label_repel(data = . %>% filter(wi.ep<0.05) %>% unique(),  aes(label = Bug), size=1.5, color = "black")


left_join( Aldex_results %>% select(Bug, effect), Aldex_replication2%>% select(Bug, effect), by="Bug", suffix=c("_Study", "_validation") ) %>% ggplot(aes(x=effect_Study, y=effect_validation)) + geom_point() + theme_bw() + geom_abline()+ geom_smooth(method = "lm")






##Association of Biomarkers
readxl::read_excel("Data/Biomarkers.xlsx") -> Biomarkers
Data_c$ID %>% sapply( function(x){ if( grepl("V2", x)  ){ substring(str_split(x, "_")[[1]][2], 2, 5)   }else( return(x) )  }  ) -> IDs

Data_c %>% mutate(ID = as.vector(IDs) ) %>% filter(!grepl("LL", ID) ) -> Data_c2
Data_c2 %>% dplyr::left_join(., Biomarkers, by="ID") %>% drop_na() %>% arrange(ID) -> Data_c2

#Run_Associations(Prevalence, All_taxonomy %>% mutate(ID = IDs) %>% filter(ID %in% Data_c2$ID) %>% arrange()  , Data_c2, FILTER=10, Function=Association_analysis, Run_permutations = F, Independent = "Zonulin", Permute_p = "Pr(>|t|)" ) -> Association_results_zn2

Run_Associations(Prevalence, All_taxonomy %>% mutate(ID = IDs) %>% filter(ID %in% Data_c2$ID) %>% arrange()  , Data_c2, FILTER=10, Function=Association_analysis, Run_permutations = F, Independent = "Zonulin") -> Association_results_zn
Run_Associations(Prevalence, All_taxonomy %>% mutate(ID = IDs) %>% filter(ID %in% Data_c2$ID) %>% arrange()  , Data_c2, FILTER=10, Function=Association_analysis, Run_permutations = F, Independent = "Alfa1antitrypsine") -> Association_results_alfa
Run_Associations(Prevalence, All_taxonomy %>% mutate(ID = IDs) %>% filter(ID %in% Data_c2$ID) %>% arrange()  , Data_c2, FILTER=10, Function=Association_analysis, Run_permutations = F, Independent = "Calprotectin") -> Association_results_cal


rbind(rbind(Association_results_zn %>% mutate(Biomarker = "Zonulin"), Association_results_alfa%>% mutate(Biomarker = "Alfa1antitrypsine")), Association_results_cal %>%  mutate(Biomarker = "Calprotectin") ) -> Results_biomarkers
#Control for FDR between all test per taxa level
Results_biomarkers2 =tibble()
for (TL in unique(Results_biomarkers$Taxonomic_level) ){
  Results_biomarkers %>% filter(Taxonomic_level == TL ) %>% mutate(FDR_all_bh = p.adjust(`Pr(>|t|)`, "fdr"), FDR_all_bh_norm = p.adjust(`Pr(>|t|)_norm`, "fdr")   ) %>% rbind(Results_biomarkers2 , . ) -> Results_biomarkers2
}

Results_biomarkers2 %>% write_tsv(. , "Results/Biomarkers_summarystats.tsv")

Results_biomarkers %>% left_join(Balances_results) %>% filter(Taxonomic_level == "s") %>% ggplot(aes(x=-log10(FDR_bh), y = -log10(FDR_bh_norm ), col=Direction_balance )) + geom_point() + theme_bw() + geom_abline() + geom_vline(xintercept = -log10(0.05) ) + geom_hline(yintercept = -log10(0.05) ) + facet_wrap(~Biomarker)
Results_biomarkers2 %>% left_join(Balances_results) %>% filter(Taxonomic_level == "s") %>% ggplot(aes(x=-log10(FDR_all_bh), y = -log10(FDR_all_bh_norm ), col=Direction_balance )) + geom_point() + theme_bw() + geom_abline() + geom_vline(xintercept = -log10(0.05) ) + geom_hline(yintercept = -log10(0.05) ) + facet_wrap(~Biomarker)




filter(Results_biomarkers, FDR_permutation < 0.05)$Bug %>% sapply(function(x){ if ( grepl("s__", x) ) { str_split(x, "\\|s__")[[1]][2] } }) -> short_names_keep
colnames(All_taxonomy) %>% sapply(function(x){ if ( grepl("s__", x) ) { str_split(x, "\\|s__")[[1]][2] } else(return(x)) }) -> short_names
All_taxonomy_figure = All_taxonomy
colnames(All_taxonomy_figure) = short_names

All_taxonomy_figure %>% mutate(ID = IDs) %>% filter(ID %in% Data_c2$ID) %>% arrange() %>% dplyr::left_join(., Biomarkers, by="ID") %>% drop_na() -> All_taxonomy_figure

All_taxonomy_figure %>%
 select("ID", as.vector(short_names_keep), "Zonulin", "Alfa1antitrypsine", "Calprotectin" ) %>% gather(Bimarker, concentration, Zonulin:Calprotectin) %>% gather("Bacterium", "CLR_abundance" ,2:8) %>% ggplot(aes(x=concentration, y=CLR_abundance)) + geom_point() + 
  theme_bw() + facet_grid(Bacterium~Bimarker, scales = "free") + geom_smooth(se=F, method = lm) + theme(strip.text.y = element_text(size = 4))


All_taxonomy %>% ggplot(aes(x=`k__Bacteria|p__Actinobacteria|c__Coriobacteriia|o__Eggerthellales|f__Eggerthellaceae|g__Eggerthellaceae gen. incertae sedis|s__Eggerthellaceae species incertae sedis [ext_mOTU_v3_16295]`)) + geom_histogram() + theme_bw()
Data_c %>% ggplot(aes(x=`k__Bacteria|p__Actinobacteria|c__Coriobacteriia|o__Eggerthellales|f__Eggerthellaceae|g__Eggerthellaceae gen. incertae sedis|s__Eggerthellaceae species incertae sedis [ext_mOTU_v3_16295]`)) + geom_histogram() + theme_bw()
All_taxonomy %>% ggplot(aes(x=Inverse_rank_normal(`k__Bacteria|p__Actinobacteria|c__Coriobacteriia|o__Eggerthellales|f__Eggerthellaceae|g__Eggerthellaceae gen. incertae sedis|s__Eggerthellaceae species incertae sedis [ext_mOTU_v3_16295]`) )) + geom_histogram() + theme_bw()

Run_aldex = function(Data, Bm="Zonulin"){
  Aldex_results_biom = tibble()
  Data$ID %>% sapply( function(x){ if( grepl("V2", x)  ){ substring(str_split(x, "_")[[1]][2], 2, 5)   }else( return(x) )  }  ) -> IDs
  for (Taxonomic_level in c("p", "c", "o", "f", "g", "s")){
    Prevalence %>% filter(N_case > 20) -> To_run 
  
    print( paste0("Subsetting features of taxonomic level: ", Taxonomic_level) )
    paste(c(Taxonomic_level, "__"), collapse="" ) -> Ta
  
    colnames(dplyr::select(Data, - one_of(c("ID" )))) -> Taxa
    Select_taxonomy_level(Taxa, Ta) -> Taxa
  
    Data %>% mutate(ID = IDs) %>% dplyr::select( one_of(c("ID", Taxa) )) %>% select_if( colnames(.) %in% c("ID", To_run$Bug)) %>% filter(ID %in% Data_c2$ID ) %>% arrange(ID) -> Data_taxonomyX
    Data_taxonomyX %>% as.data.frame() %>% column_to_rownames("ID") %>% t() -> Input_aldex
    #print("Calling association function" )
    mm <- model.matrix(as.formula(paste0("~ Age + Sex + ", Bm)) , Data_c2)
    x <- ALDEx2::aldex.clr(Input_aldex, mm, mc.samples=160, denom="all", verbose=FALSE)
    x.tt <- Check(x, mm, verbose=FALSE) #this function might not work... I wrote a modified version in which the glm call does not have the "..." and works, find at the bottom
    #x.tt2 = ALDEx2::aldex.corr(x, Data_c2$Zonulin )
    #x.effect <- ALDEx2::aldex.glm.effect(x, mm, CI=T, verbose=FALSE )
    
    x.all <- data.frame(x.tt)
    x.all %>% rownames_to_column("Bug") %>% as_tibble() %>% mutate(Taxonomic_level = Taxonomic_level) -> Result_taxa
    colnames(Result_taxa)[14:17] = c("Estimate", "SE", "T_value", "P")
    colnames(Result_taxa)[21] = c("FDR")
    #print("Association succesful" )
    rbind(Aldex_results_biom, Result_taxa) -> Aldex_results_biom
  }
  return(Aldex_results_biom)
}  
Check = function (clr, verbose = FALSE, ...) {
  conditions <- clr@conds
  lr2glm <- function(lr, conditions, ...) {
    if (!is(conditions, "matrix") && !("assign" %in% names(attributes(conditions)))) {
      stop("Please define the aldex.clr object for a model.matrix 'conditions'.")
    }
    if (nrow(lr) != nrow(conditions)) {
      stop("Input data and 'model.matrix' should have same number of rows.")
    }
    model. <- conditions
    #print("Applying model")
    #print(as_tibble(lr))
    #print(model.)
    #print("now...")
    glms <- apply(lr, 2, function(x) {
      glm(x ~ model.)
    })
    extract <- function(model) {
      x <- coef(summary(model))
      coefs <- lapply(1:nrow(x), function(i) {
        y <- x[i, , drop = FALSE]
        colnames(y) <- paste(rownames(y), colnames(y))
        y
      })
      do.call("cbind", coefs)
    }
    #print("Extracting coefficients")
    extracts <- lapply(glms, extract)
    #print("Binding")
    df <- do.call("rbind", extracts)
    rownames(df) <- colnames(lr)
    df <- as.data.frame(df)
    pvals <- colnames(df)[grepl("Pr\\(>", colnames(df))]
    df.bh <- df[, pvals]
    colnames(df.bh) <- paste0(colnames(df.bh), ".BH")
    for (j in 1:ncol(df.bh)) {
      df.bh[, j] <- p.adjust(df.bh[, j])
    }
    cbind(df, df.bh)
  }
  if (verbose) 
    message("running tests for each MC instance:")
  mc <- ALDEx2::getMonteCarloInstances(clr)
  k <- ALDEx2::numMCInstances(clr)
  r <- 0
  for (i in 1:k) {
    if (verbose == TRUE) 
      numTicks <- progress(i, k, numTicks)
    mci_lr <- t(sapply(mc, function(x) x[, i]))
    r <- r + lr2glm(mci_lr, conditions, ...)
  }
  r/k
}


Alex_results_biomarkers = tibble()  
for (Biomarker in colnames(Biomarkers) ){
  if (Biomarker %in% c("ID", "Biovis_Batch") ){ next }
  Run_aldex(Data, Biomarker) %>% mutate(Biomarker = Biomarker) %>% rbind(. , Alex_results_biomarkers) -> Alex_results_biomarkers

}
Interest= c("Bug",  "Estimate", "SE", "T_value", "P", "FDR" )

left_join(Alex_results_biomarkers, Results_biomarkers, by=c("Bug", "Biomarker", "Taxonomic_level")) %>% left_join(Balances_results) %>% filter(Taxonomic_level == "s") %>% 
  ggplot(aes(x=-log10(FDR_bh), y = -log10(FDR), col=Direction_balance )) + geom_point() + theme_bw() + geom_abline() + geom_vline(xintercept = -log10(0.05) ) + geom_hline(yintercept = -log10(0.05) ) + facet_wrap(~Biomarker)
left_join(Alex_results_biomarkers, Results_biomarkers, by=c("Bug", "Biomarker", "Taxonomic_level")) %>% left_join(Balances_results) %>% filter(Taxonomic_level == "s") %>% 
  ggplot(aes(x=-log10(FDR_bh_norm), y = -log10(FDR), col=Direction_balance )) + geom_point() + theme_bw() + geom_abline() + geom_vline(xintercept = -log10(0.05) ) + geom_hline(yintercept = -log10(0.05) ) + facet_wrap(~Biomarker)  
  
  




left_join(Results_biomarkers, Alex_results_biomarkers, by = c("Bug","Biomarker", "Taxonomic_level")) %>% left_join(., Prevalence) -> Merged_biomarkers
Merged_biomarkers %>%  write_tsv( ., "Results/Summary_statistics_TaxonomyAnalysis_biomarkers.tsv")  


Merged_biomarkers %>% filter(FDR_permutation < 0.05) %>% select(N_case)
Merged_biomarkers %>% filter(FDR_permutation < 0.05) %>% select(Bug, Biomarker)
#Checking significant assocations.
Check_plots_assocation = function(Long_name, Short_name, Marker){
  
  x1 = as_vector(All_taxonomy_figure[Short_name])
  x2 = Inverse_rank_normal(as_vector(All_taxonomy_figure[Short_name]))
  x3 = as_vector(Data_c2[Long_name]) > 0
  y1 = as_vector(All_taxonomy_figure[Marker])
  
  New_v = tibble(Bug1 = x1, Bug2 = x2 , Bug3 =  x3,  Marker = y1)
  Plot1 = New_v %>% ggplot(aes(x=Bug1, y = Marker)) + geom_point() + theme_bw() + geom_smooth(se=F, method = lm) + xlab(Short_name) + ylab(Marker) + ggtitle("CLR")
  Plot2 = New_v %>% ggplot(aes(x=Bug2, y = Marker)) + geom_point() + theme_bw() + geom_smooth(se=F, method = lm) + xlab(Short_name) + ylab(Marker) + ggtitle("InvRank_CLR")
  Plot3 = New_v %>% ggplot(aes(x=Bug3, y = Marker)) + geom_boxplot() + geom_jitter() + theme_bw() + geom_smooth(se=F, method = lm) + xlab(Short_name) + ylab(Marker) + ggtitle("Precence/Absence")
  Plot_bug = Plot1 | Plot2 | Plot3
  print(Plot_bug)
  
}

#1. Eggerthela and Zonulin.
Check_plots_assocation("k__Bacteria|p__Actinobacteria|c__Coriobacteriia|o__Eggerthellales|f__Eggerthellaceae|g__Eggerthellaceae gen. incertae sedis|s__Eggerthellaceae species incertae sedis [ext_mOTU_v3_16295]","Eggerthellaceae species incertae sedis [ext_mOTU_v3_16295]", "Zonulin")
Check_plots_assocation("k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Bacteroidaceae|g__Bacteroides|s__Bacteroides intestinalis [ref_mOTU_v3_02809]","Bacteroides intestinalis [ref_mOTU_v3_02809]", "Zonulin") #Prevalence <0.2
Check_plots_assocation("k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Prevotellaceae|g__Prevotella|s__Prevotella species incertae sedis [meta_mOTU_v3_12510]","Prevotella species incertae sedis [meta_mOTU_v3_12510]", "Zonulin") #Prevalence < 0.2
#Not clear differences in presence/absence, not sure 
Check_plots_assocation("k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Clostridiales fam. incertae sedis|g__Clostridiales gen. incertae sedis|s__Clostridiales species incertae sedis [ext_mOTU_v3_26621]","Clostridiales species incertae sedis [ext_mOTU_v3_26621]", "Zonulin") #Prevalence < 0.2
Check_plots_assocation("k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Lachnospiraceae|g__Blautia|s__[Ruminococcus] torques [ref_mOTU_v3_03703]","[Ruminococcus] torques [ref_mOTU_v3_03703]", "Alfa1antitrypsine")
Check_plots_assocation("k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Ruminococcaceae|g__Ruminococcaceae gen. incertae sedis|s__Ruminococcaceae species incertae sedis [meta_mOTU_v3_13455]","Ruminococcaceae species incertae sedis [meta_mOTU_v3_13455]", "Alfa1antitrypsine")
#Strange association in CLR
Check_plots_assocation("k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Lachnospiraceae|g__Blautia|s__[Ruminococcus] torques [ref_mOTU_v3_03703]","Ruminococcaceae species incertae sedis [meta_mOTU_v3_13455]", "Calprotectin")
Check_plots_assocation("k__Bacteria|p__Firmicutes|c__Negativicutes|o__Acidaminococcales|f__Acidaminococcaceae|g__Acidaminococcus|s__Acidaminococcus intestini [ref_mOTU_v3_01949]","Acidaminococcus intestini [ref_mOTU_v3_01949]", "Calprotectin") #Prevalence < 0.2

#Significant in non-inv rank
Check_plots_assocation("k__Bacteria|p__Actinobacteria|c__Coriobacteriia|o__Eggerthellales|f__Eggerthellaceae|g__Eggerthella|s__Eggerthella lenta [ref_mOTU_v3_00719]","Eggerthella lenta [ref_mOTU_v3_00719]", "Zonulin")
Check_plots_assocation("k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Oscillospiraceae|g__Oscillibacter|s__Oscillibacter species incertae sedis [ext_mOTU_v3_17772]","Oscillibacter species incertae sedis [ext_mOTU_v3_17772]", "Zonulin")
Check_plots_assocation("k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Clostridiales fam. incertae sedis|g__Clostridiales gen. incertae sedis|s__Clostridiales species incertae sedis [meta_mOTU_v3_13390]","Clostridiales species incertae sedis [meta_mOTU_v3_13390]", "Zonulin")



#Test presence/absence
Binary_association = function(D,Items, Prevalence_t){
  Taxonomic_level_result = tibble()
  for (Bug in dplyr::filter( dplyr::filter(Prevalence, Bug %in% Items) , N_case>Prevalence_t & N_case < 103-Prevalence_t  )$Bug ){
    Formula = paste0(Biomarker, "~ `" , Bug, "`")
    Data_c3 = D
    Data_c3[Bug] = D[Bug] > 0
    lm( Formula ,  Data_c3) %>% summary() -> Model
    as.data.frame(Model$coefficients)[paste0("`", Bug, "`"),] %>% rownames_to_column("Bug") %>% mutate(Biomarker= Biomarker, Taxonomic_level=str_replace(Taxa_level, "__", "") ) %>% rbind(Taxonomic_level_result, .) -> Taxonomic_level_result
  }
  Taxonomic_level_result %>% mutate(FDR_bh_level = p.adjust(`Pr(>|t|)`, "fdr") ) -> Taxonomic_level_result
  return(Taxonomic_level_result)
}
Results_binary = tibble()
Permutation = F
for (Biomarker in colnames(Biomarkers)){
  if (Biomarker %in% c("ID", "Biovis_Batch") ){ next }
  print(Biomarker)
  for ( Taxa_level in c("p__", "c__", "o__", "f__", "g__", "s__") ){
    Select_taxonomy_level(Prevalence$Bug, Taxa_level) -> Items
    Taxonomic_level_result = tibble()
    Binary_association(Data_c2, Items, 10) %>% as_tibble() %>% drop_na() -> Summary_statistics_taxonomyX 
    
  
    if (Permutation == T){
      print("Running permutations" )
      Null_distribution = c()
    for (Permutation_n in seq(100) ){
        print(paste0("running permutation: ", Permutation_n))
        Data_perm = Data_c2
        Data_perm[Biomarker] = sample( as_vector(Data_c2[Biomarker]), replace = F,size = dim(Data_c2)[1] )
        Binary_association(Data_perm, Items, 20) %>% as_tibble()  %>% drop_na()  -> Association_results_p
        Null_distribution = c(Null_distribution, Association_results_p$`Pr(>|t|)` ) 
      }
    Summary_statistics_taxonomyX["P"] = as_vector(Summary_statistics_taxonomyX["Pr(>|t|)"])
    FDR_perm = Compute_FDR_null(Summary_statistics_taxonomyX , Null_distribution, fdr_threshold=0.05  ) 
    FDR_perm[FDR_perm > 1] = 1
    Summary_statistics_taxonomyX %>% mutate(FDR_permutation = FDR_perm ) -> Summary_statistics_taxonomyX
    print("Permutations succesful" )
    }
    
    
    Results_binary = rbind ( Results_binary, Summary_statistics_taxonomyX )
      
    }
  }    

#Control for FDR between all test per taxa level
Results_binary2 =tibble()
for (TL in unique(Results_binary$Taxonomic_level) ){
  Results_binary %>% filter(Taxonomic_level == TL ) %>% mutate(FDR_all_bh = p.adjust(`Pr(>|t|)`, "fdr")  ) %>% rbind(Results_binary2 , . ) -> Results_binary2
}
Results_binary2$Bug =  str_replace(Results_binary2$Bug, "TRUE", "") %>% str_replace_all(. , "`", "") 
left_join(Results_binary2, Results_biomarkers2, by=c("Bug", "Biomarker", "Taxonomic_level"), suffix=c("_Binary", "_Continuous") ) %>% left_join(Balances_results) %>% filter(Taxonomic_level == "s") %>% 
  ggplot(aes(x=-log10(FDR_all_bh_Binary), y = -log10(FDR_all_bh_Continuous), col=Direction_balance )) + geom_point() + theme_bw() + geom_abline() + geom_vline(xintercept = -log10(0.05) ) + geom_hline(yintercept = -log10(0.05) ) + facet_wrap(~Biomarker) -> Plot1
left_join(Results_binary2, Results_biomarkers2, by=c("Bug", "Biomarker", "Taxonomic_level"), suffix=c("_Binary", "_Continuous") ) %>% left_join(Balances_results) %>% filter(Taxonomic_level == "s") %>% 
  ggplot(aes(x=-log10(FDR_all_bh_Binary), y = -log10(FDR_all_bh_norm), col=Direction_balance )) + geom_point() + theme_bw() + geom_abline() + geom_vline(xintercept = -log10(0.05) ) + geom_hline(yintercept = -log10(0.05) ) + facet_wrap(~Biomarker) -> Plot2
Plot1 | Plot2 
Results_binary2 %>%  write_tsv( ., "Results/Summary_statistics_biomarkers_binary.tsv")  



colnames(Data_c2)[grepl("s__", colnames(Data_c2))] -> For_beta
Data_c2 %>% select(For_beta) %>% vegdist(method = "jaccard") -> JD
Data_c2 %>% select(For_beta) %>% vegdist(method = "robust.aitchison") -> AD
#Binary matrix
adonis2(JD ~ Zonulin , Data_c2) #0.017 / R2 0.012
adonis2(JD ~ Alfa1antitrypsine , Data_c2) #0.106 / R2 0.011
adonis2(JD ~ Calprotectin , Data_c2) #0.148 / R2 0.11
#Aitchison distance
adonis2(AD ~ Zonulin , Data_c2) #0.033 / R2 0.012
adonis2(AD ~ Alfa1antitrypsine , Data_c2) #0.226 / R2 0.10
adonis2(AD ~ Calprotectin , Data_c2) #0.335 / R2 0.010

All_taxonomy_figure$ID = IDs
All_taxonomy_figure %>% left_join(Biomarkers) %>% drop_na() -> All_taxonomy_figure
Check_plots_assocation("k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Ruminococcaceae|g__Ruminococcaceae gen. incertae sedis|s__Ruminococcaceae species incertae sedis [meta_mOTU_v3_12259]","Ruminococcaceae species incertae sedis [meta_mOTU_v3_12259]", "Zonulin")





#X. Compare Gut-brain modules
read_tsv("Data/Data2/GBM_quant.tsv") -> GBM
GBM[is.na(GBM)] = 0
GBM %>% select(-ID)  %>% apply(1, function(x) { CLR(x+1)  }) %>% t() %>% as_tibble() %>% mutate(ID = GBM$ID) -> GBM

GBM_associaton = function(GBM, Covariates_all){
  GBM %>% left_join(dplyr::rename(Covariates_all, ID=Sequencing_ID), by="ID") -> GBM_test
  Results_gbm = tibble()
  for (Module in colnames(GBM)){
    if (Module == "ID"){ next }
    Formula = paste0( Module, " ~ Status + Reads_clean"  )
    lm(Formula, GBM_test) %>% summary() -> result
    as.data.frame(result$coefficients)["Status",] %>% as_tibble() %>% mutate(Module = Module) %>% rbind(Results_gbm, .) -> Results_gbm
  }
  return(Results_gbm)
}
GBM_associaton(GBM, Covariates_all) -> Results_gbm
Null_gbm = tibble()
N_perm = 100
#Permutations
for (i in seq(1:N_perm)){
  sample(Covariates_all$Sequencing_ID, length(Covariates_all$Sequencing_ID), replace = F) -> NID
  Covariates_perm = Covariates_all
  Covariates_perm$Sequencing_ID = NID
  GBM_associaton(GBM, Covariates_perm) -> Permuted
  rbind(Null_gbm, Permuted) -> Null_gbm
}
FDR_p = c()
for (P in Results_gbm$`Pr(>|t|)`){
  FDR_p = c(FDR_p, sum(Null_gbm$`Pr(>|t|)` <= P )/N_perm)
}
FDR_p[FDR_p > 1] = 1
Results_gbm %>% mutate(FDR_p = FDR_p) %>% arrange(`Pr(>|t|)`) -> Results_gbm
Results_gbm %>% ggplot(aes(x=Estimate, y=-log10(`Pr(>|t|)`), col=FDR_p<0.05 )) + geom_point() + geom_hline(yintercept = -log10(0.05)) + theme_bw() +
  geom_label_repel(data = . %>% filter( FDR_p<0.05) %>% unique(),  aes(label = Module), size=1.5, color = "black")
#MGB052	Butyrate synthesis I --> significantly reduced 
#MGB048	Propionate synthesis I --> increased
#MGB022	GABA synthesis III --> decreased
#MGB054	Propionate synthesis II --> increased
read_tsv("Data/Data2/Annotation_GBM.tsv", col_names = F) -> Annotation
colnames(Annotation) = c("Module","Annotation")
left_join(Results_gbm, Annotation) -> Results_gbm
write_tsv(Results_gbm, "Results/Sumamry_stats_GBM.tsv")



#3. Network analysis
##3.1 Generate input for sparseCC
###1. Filter species level
Data %>% dplyr::select( ID, colnames(Data)[grepl("s__",colnames(Data)) ]  ) -> Data_for_network
###2. Filter prevalence
Prevalence %>% filter(N_case > 10 & N_control > 10) -> To_keep
Data_for_network %>% dplyr::select(one_of(c("ID", To_keep$Bug))) -> Data_for_network
###2. Format
Data_for_network %>% as.data.frame() %>% column_to_rownames("ID") %>% t() %>% as.data.frame() %>% rownames_to_column("Taxa") %>% as_tibble() -> Data_for_network
###3. Split Case/Control and save
dplyr::select(Data_for_network, c("Taxa", Covariates_control$Sequencing_ID) ) %>% write_tsv(. , "Data/SparseCC_input/Controls.tsv")
dplyr::select(Data_for_network, one_of(c("Taxa", Covariates_GUTS$Sequencing_ID)) ) %>% write_tsv(. , "Data/SparseCC_input/Cases.tsv")
##3.2 Running SparseCC
###To run sparseCC you need to execute: >bash Run_SparseCC.sh  On the terminal
###Correlation between two taxa can be inferred from its ALR, since the var(ALR) among samples is informative of the correlation between two taxa. Instead of fractions or raw counts, probabilities from a Dirichlet are used 
##3.3 Reading correlation matrices
Cor_cases = read.table("Results/SparseCC/Cases/Corr_matrix_cases.txt", header = TRUE, check.names = F, sep = "\t", row.names = 1)
Cor_controls = read.table("Results/SparseCC/Controls/Corr_matrix_controls.txt", header = TRUE, check.names = F, sep = "\t", row.names = 1)
colnames(Cor_cases) %>% sapply(function(x){ str_split(x, "\\|")[[1]] -> y ;  y[length(y)] -> z ;  str_split(z, " ")[[1]] -> W ; return(W[length(W)])    }) -> New_names
colnames(Cor_cases) = New_names ; rownames(Cor_cases) = colnames(Cor_cases)
colnames(Cor_controls) = New_names ; rownames(Cor_controls) = colnames(Cor_controls)

###Get Pvalues, filter correlations to be only significant
Ps_cases = read.table("Results/SparseCC/Cases/P_matrix.txt", header = TRUE, check.names = F, sep = "\t", row.names = 1)
Ps_controls = read.table("Results/SparseCC/Controls/P_matrix.txt", header = TRUE, check.names = F, sep = "\t", row.names = 1)
###Filter out non-significant correlations
Cor_cases[!Ps_cases<0.005] = 0
Cor_controls[!Ps_controls<0.005] = 0
###filter out weak correlations
Cor_cases[Cor_cases<0.25] = 0
Cor_controls[Cor_controls<0.25] = 0


##3.5 Build network
igraph::graph_from_adjacency_matrix( as.matrix(Cor_controls), mode='undirected', weighted = 'correlation',add.rownames = T) -> EdgeObject_controls
igraph::graph_from_adjacency_matrix( as.matrix(Cor_cases), mode='undirected', weighted = 'correlation', add.rownames = T) -> EdgeObject_cases
igraph::write_graph(EdgeObject_cases, "Results/SparseCC/Network_cases.txt", format="ncol")
igraph::write_graph(EdgeObject_controls, "Results/SparseCC/Network_controls.txt", format="ncol") #E(EdgeObject_controls) %>% length() should be < 1000 to work in webserver
##3.5 Run Netshift
###Go to https://web.rniapps.net/netshift/ and upload the two networks
###Output:
###Jaccard similarity between cases and controls: 0.314
read_tsv("Results/SparseCC/Netshift_node_stats.tsv") -> Node_stats
Node_stats %>% left_join(. , Results) %>% mutate(Log_p = -log10(wi.ep)) %>% select(`NESH-score`, DelBet, Log_p) %>% cor(method = "spearman") %>% corrplot::corrplot(diag=F)
#There's a certain positive association between:
#NESH-score: Similarity in number of connections between a node in two networks, weigthed by the exclusive connections in cases. The more the more different
#Delta Betweenes: betweenes --> centrality of a certain node ; this is the difference of centraility in the node in cases in comparison with controls
sapply(Node_stats$S_ID, function(x){ str_replace(x, "\\]", "") %>% str_replace("\\[", "") -> y ; colnames(Data)[grepl( y, colnames(Data) )] }) -> Names
Node_stats %>% mutate(Bug = Names) -> Node_stats
summary(Node_stats$DelBet) -> Between_stats 
summary(Node_stats$`NESH-score`) -> NESH_stats
Node_stats %>% left_join(. , Results) %>% ggplot(aes(x=`NESH-score`, y=DelBet, col=-log10(wi.ep))) + geom_point() + theme_bw() + geom_hline(yintercept = Between_stats[5]) + geom_vline(xintercept = NESH_stats[5]) + scale_colour_gradient()
Node_stats %>% filter(DelBet>Between_stats[5] & `NESH-score`>NESH_stats[5]) %>% select(Bug) 


#4. Subspecies level analysis

Subspecies_analysis = function(File, Distance="mann", FILTER=10, Covariates_all = Covariates_all, P = "Data/Data2/metasnv_distances-m2-d5-b80-c5-p0.9/"){
  #if (! Distance == "mann" ){ Distance = "allele" }
  
  File_path = paste0(P, File)
  mOTU = str_replace(File, paste0(".filtered.",Distance,".dist"), "")
  read.table(File_path, header = TRUE, check.names = F, sep = "\t", row.names = 1) -> Distance_matrix
  rownames(Distance_matrix) %>% sapply(., function(x){str_split(x,".bam")[[1]][1] } ) %>%as.vector() -> Sample_names 
  rownames(Distance_matrix) = Sample_names ; colnames(Distance_matrix) = Sample_names
  Distance_matrix %>% rownames_to_column("ID") %>% as_tibble() %>% filter(ID %in% Covariates_all$Sequencing_ID) %>% select_if( colnames(.) %in% c("ID", Covariates_all$Sequencing_ID) ) %>% column_to_rownames("ID") %>% as.data.frame() -> Distance_matrix
  
  Covariates_all %>% filter(Sequencing_ID %in% Sample_names) %>% group_by(Status) %>% summarise(N = n()) -> Counts
  if ( Counts$N[1] < FILTER  |   Counts$N[2] < FILTER ){ 
    return( list(tibble(mOTU=mOTU, Stat=NA, N_control= Counts$N[1],  N_cases = Counts$N[2], Distance=Distance ), NA ) )
  }
  
  Covariates_all %>% filter(Sequencing_ID %in% Sample_names) -> Covariates_sort
  Covariates_sort[match(rownames(Distance_matrix), Covariates_sort$Sequencing_ID),] -> Covariates_sort
  
  as.matrix(Distance_matrix) %>% as.dist() -> Distance_matrix
  
  capscale(Distance_matrix ~ Status, Covariates_sort   ) -> RDA
  summary(RDA)$sites %>% as_tibble() %>% cbind(. , Covariates_sort) %>% ggplot(aes(x=CAP1, y=MDS1, col=as.factor(Status))) + geom_point() + theme_bw() + ggtitle(mOTU) -> Plot
  anova(RDA, permutations = 999) -> Stat
  tibble(mOTU = mOTU, Stat =  Stat$`Pr(>F)`[1], N_control= Counts$N[1],  N_cases = Counts$N[2],  Distance=Distance ) -> Result
  if( Stat$`Pr(>F)`[1]<0.05) { print(Plot) }
  return(list(Result, Plot))
  
}


read_delim("Data/Data2/metasnv_distances-m2-d5-b80-c5-p0.9/Number_samples_distance.txt", delim=" ", col_names = F) -> Samples_discovery
read_delim("Data/Data2/distances-m2-d5-b80-c5-p0.9_withRep/Number_samples_distance.txt", delim=" ", col_names = F) %>% mutate(X1 = as.numeric(X1)) -> Samples_discovery


Samples_discovery %>% mutate(X1 = as.numeric(X1)) %>% filter(X1 > 10)  -> Manhattan_to_check

Permanova_stats= tibble()
for (File in Manhattan_to_check$X2){
  if (grepl("allele", File)){ D = "allele"
  } else {
    D = "mann"
  }
  print(File)
  Subspecies_analysis(File, Distance=D, FILTER=5, Covariates_all = Cov_integrated %>% rename(Sequencing_ID  = ID), P="Data/Data2/distances-m2-d5-b80-c5-p0.9_withRep/" ) -> Res
  Res[[1]] %>% rbind(Permanova_stats, . ) -> Permanova_stats
  
}
sapply(Permanova_stats$mOTU, function(x){ y = colnames(Data)[grepl( x, colnames(Data) )] ; if (length(y)==0){ return(NA) }else{ return(y) }   }) -> Names
Permanova_stats %>% mutate(Taxa = unlist(Names)) %>% arrange(Stat) -> Permanova_stats
Permanova_stats %>% filter(Distance == "allele") %>% mutate(FDR = p.adjust(Stat, "fdr"))
#Nominal significance evidence of subspecies differences:
# s__Sutterella species incertae sedis [meta_mOTU_v3_13005]
# Clostridiales species incertae sedis [meta_mOTU_v3_12974]
# Not in differential abundance or network  






