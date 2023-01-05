setwd("~/Documents/GitHub/GUTS_MGS_Psych/")

####Load libraries####

library(tidyverse)
library(vegan)
library(ape)
library(MetBrewer) #Lets plot with some art :)
library(ggrepel)
ggplot <- function(...) ggplot2::ggplot(...) + scale_color_manual(values=met.brewer("Manet",5)[c(1,5)] ) +
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

#0. Samples info
##0.1 Metadata samples
###Cases
Covariates_GUTS = read_tsv("Data/MetaData/Metadata_GUTS_final2.tsv") #%>% dplyr::rename(Sample_ID = ID)
####Keep only baseline samples
Covariates_GUTS %>% dplyr::filter(Treatment_time == "V2") -> Covariates_GUTS
####Check repeated samples  
Covariates_GUTS %>%  group_by(Sample_ID, Batch) %>% summarise(N = n() ) %>% arrange(desc(N)) -> Repeated #19 repeated samples, 18 in batch 2, 1 in batch 1
###Controls
Covariates_control = read_tsv("Data/MetaData/Covariates_population.tsv") %>% dplyr::rename(Sequencing_ID = ID, Read_number = META.DNA.postclean.reads, Batch=META.BATCH, Sex=Gender) %>% drop_na()
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

Compute_CLR_taxonomy(Data, Taxonomy = "s", Do_pseudo = T, Keep_column = "") -> Species_data
All_taxonomy = Species_data
for (i in c("p", "c", "o", "f", "g")){
  Compute_CLR_taxonomy(Data, Taxonomy = i, Do_pseudo = T, Keep_column = "") -> Transformed_data
  left_join(All_taxonomy, Transformed_data, by = "ID") -> All_taxonomy
}    
All_taxonomy <- Map(function(x) replace(x, is.infinite(x), NA), All_taxonomy) %>% as_tibble()
Species_data <- Map(function(x) replace(x, is.infinite(x), NA), Species_data) %>% as_tibble()

##1.3 Beta diversity
###Compute beta diversity at the species level
Beta_taxonomy = function(Data_t, Data_c, Data, Meta = "Reads_clean"){
  "s__" -> Ta
  colnames(select(Data, -c("ID"))) -> Taxa
  Taxa[grepl(Ta, Taxa)] -> Taxa
  Data %>% select( c("ID", Taxa) ) -> Data_taxonomy
  
  
  #Compute diversity
  vegdist(select(Data_taxonomy, -c(ID)), method = "robust.aitchison") -> Beta_diversity
  
  
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
  vegdist(select(Data_cases, -ID), method = "robust.aitchison") -> Beta_diversity2
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
Beta_taxonomy(Species_data, Data_c, Data, Meta="Reads_clean")

##1.4 Alpha diversity
###Compute alpha diversity at the species level ; need to do this in non-transformed data
Alpha_taxonomy = function(Data, Data_c, Taxonomy = "s", Meta = Covariates_all, Remove_names = c("ID", "UNKNOWN"), Cov_name=c("Reads_clean") ){
  paste( c(Taxonomy, "__"), collapse="" ) -> Ta
  colnames(select(Data, -one_of(Remove_names))) -> Taxa
  Select_taxonomy_level(Taxa, Ta) -> Taxa
  Data %>% select(Taxa)  -> Data_taxonomy
  diversity(Data_taxonomy, "shannon" ) -> Diversity
  Data_c %>% select(c("ID", "Status", "Batch", Cov_name)) %>% mutate(Diversity = Diversity) -> Data_div

  Formula = as.formula(paste0("Diversity ~ Status +", paste(Cov_name, collapse=" + ") ))
  summary(lm(Formula,Data_div)) %>% print()
  summary(lm( as.formula(paste0("Diversity ~ Batch +",paste(Cov_name, collapse=" + ") ))  ,filter(Data_div, Status == 1 ) )) %>% print()
  
  
  ggplot(Data_div, aes(x=as.factor(Status), y=Diversity, col=as.factor(Status)))  + geom_violin() +
  theme_bw() + coord_flip() + geom_jitter() +
    stat_summary(fun = "median",geom = "crossbar", aes(color = as.factor(Status))) +
    stat_summary(fun = "mean", geom = "point", color= "black") + ggtitle("Disease differences") -> Plot1
  
  ggplot( filter(Data_div, Status==1) , aes(x=as.factor(Batch), y=Diversity, col=as.factor(Batch)))  + geom_violin() +
    theme_bw() + coord_flip() + geom_jitter() +
    stat_summary(fun = "median",geom = "crossbar", aes(color = as.factor(Batch))) +
    stat_summary(fun = "mean", geom = "point", color= "black")  + ggtitle("Batch differences (GUTS)")-> Plot2
  
  return(list(Plot1, Plot2))
  
}
Alpha_taxonomy(Data, Data_c, "s")


#2. Biomarker discovery

##2.1 Use coda core to predict disease
library(codacore)
library(tensorflow)
library(pROC)
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
  yHat <- predict(model, newx = xTest)
  print("Test results")
  testAUC <- pROC::auc(pROC::roc(yTest, yHat, quiet=T))
  cat("Test set AUC =", testAUC, "\n")
  
  failure <- yHat < 0.5 ; success <- yHat >= 0.5
  yHat[failure] <- levels(yTest)[1] ; yHat[success] <- levels(yTest)[2]
  cat("Classification accuracy on test set =", round(mean(yHat == yTest), 2))
  
  plot(model)
  plotROC(model)
  Numerators = colnames(xTrain)[getNumeratorParts(model, 1)] 
  Denominators = colnames(xTrain)[getDenominatorParts(model, 1)]
 
  if (length(model$ensemble) > 1){ 
    Numerators2 = colnames(xTrain)[getNumeratorParts(model, 2)] 
    Denominators2 = colnames(xTrain)[getDenominatorParts(model, 2)]
    return(list(Numerators, Denominators, Numerators2, Denominators2))
    
  } else { return(list(Numerators, Denominators)) }
}   

#Get only prevalent bacteria for the model
Prevalence %>% filter(N_case > 20 & N_control > 20) -> Keep
Run_balance_analysis( select(Data, c("ID", Keep$Bug)) ) -> Balances
#Balance using batches as different train/test sets
#Run_balance_analysis( select(Data, c("ID", Keep$Bug)), add_pseudo = "Imput" , Strategy="Batch"  ) -> Balances2
Balances[[4]] %>% sapply(function(x){ str_split(x, "\\|")[[1]] -> y ; y[length(y)]  } ) %>% as_vector() %>% as.vector()
Data %>% select(Balances[[1]]) %>% apply(1, sum) -> Numerator ; Data %>% select(Balances[[3]]) %>% apply(1, sum) -> Numerator2
Data %>% select(Balances[[2]]) %>% apply(1, sum) -> Denominator ; Data %>% select(Balances[[4]]) %>% apply(1, sum) -> Denominator2
Data %>% mutate( Balance = log10(Numerator/Denominator), Balance2=log10(Numerator2/Denominator2), Status = Data_c$Status, Read_number=Data_c$Reads_clean ) ->Data_b
Data_b %>% ggplot(aes(x=as.factor(Status), y= Balance, col=as.factor(Status)))  + geom_violin() + theme_bw() + coord_flip() + geom_jitter() +
  stat_summary(fun = "median",geom = "crossbar", aes(color = as.factor(Status))) + stat_summary(fun = "mean", geom = "point", color= "black")
summary(lm(Balance ~ Status + Read_number , filter( Data_b, ! abs(Balance) == Inf )))
Data_b %>% ggplot(aes(x=as.factor(Status), y= Balance2, col=as.factor(Status)))  + geom_violin() + theme_bw() + coord_flip() + geom_jitter() +
  stat_summary(fun = "median",geom = "crossbar", aes(color = as.factor(Status))) + stat_summary(fun = "mean", geom = "point", color= "black")
summary(lm(Balance2 ~ Status + Read_number , filter( Data_b, ! abs(Balance2) == Inf )))

#2.2 Association analysis: linear model on CLR-transformed data 
##Inverse-rank normal transformation, also known as INT --> enforce normality of the data. Applied after sample-specific normalization factors applied in CLR
Inverse_rank_normal = function(Measurement){
  qnorm((rank(Measurement,na.last="keep")-0.5)/sum(!is.na(Measurement)))
}

Association_analysis = function(Prevalence, DF, Data_c, FILTER=20, Meta=c("Read_number")) {
  Total_results = tibble()
  Prevalence %>% filter(N_case > FILTER & N_control > FILTER) %>% filter(! Bug == "UNKNOWN") -> To_Test
  for (Bug in To_Test$Bug){
    if (! Bug %in% colnames(DF) ){ next }
    DF %>% select( one_of(Bug) ) %>% as_vector() -> vector_Bug
    #IF we want to normalize it, apply function here
    normalized_Bug = Inverse_rank_normal(vector_Bug)
    #invers rank normal transf?
    Data_c %>% select(one_of(c("ID", "Status", "Batch", Meta))) %>% mutate(B =  vector_Bug, B_n = normalized_Bug) -> Model_input
    Formula = paste( c("B ~ Status", Meta), collapse = "+")
    Formula2 = paste(c("B_n ~ Status", Meta), collapse="+")
    
    lm(Formula, Model_input  ) -> Model_out
    lm(Formula2, Model_input  ) -> Model_out_n
    #If not normalized
    Normality = shapiro.test(Model_out$residuals)$p.value
    as.data.frame(summary(Model_out)$coefficients)["Status",] %>% as_tibble() %>%
      mutate(Bug = Bug, Shapiro_p = Normality, .before=1) -> results
    #Normalized
    Normality = shapiro.test(Model_out_n$residuals)$p.value
    as.data.frame(summary(Model_out_n)$coefficients)["Status",] %>% as_tibble() %>%
      mutate(Bug = Bug, Shapiro_p = Normality, .before=1) -> Normalized_results
    
    left_join(results, Normalized_results, by= "Bug", suffix = c("","_norm")) -> results
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
Run_Associations = function( Prevalence, DF, Data_c, FILTER_n=20, Meta_n=c("Reads_clean"), Function =  Association_analysis, Run_permutations = F, Permutation_number=100 ){
  #Per taxonomic level, run Function and perform 100 permutations to obtain taxonomic-level-specific FDR
  Summary_statistics = tibble()
  for (Taxonomic_level in c("p", "c", "o", "f", "g", "s")){
    print( paste0("Subsetting features of taxonomic level: ", Taxonomic_level) )
    paste(c(Taxonomic_level, "__"), collapse="" ) -> Ta
  
    colnames(select(DF, - one_of(c("ID" )))) -> Taxa
    Select_taxonomy_level(Taxa, Ta) -> Taxa
    
    DF %>% select( one_of(c("ID", Taxa) )) -> DF_taxonomyX

    print("Calling association function" )
    Summary_statistics_taxonomyX = Function( Prevalence, DF_taxonomyX, Data_c, FILTER=FILTER_n, Meta=Meta_n )
    Summary_statistics_taxonomyX %>% mutate(FDR_bh = p.adjust(`Pr(>|t|)_norm`, "fdr") ) -> Summary_statistics_taxonomyX
    print("Association succesful" )
    
    if (Run_permutations == T){
        print("Running permutations" )
        Null_distribution = c()
        for (Permutation_n in seq(Permutation_number) ){
              Data_perm = Data_c
              Data_perm$Status = sample(Data_c$Status, replace = F,size = dim(Data_c)[1] )
              Association_analysis(Prevalence, DF_taxonomyX, Data_perm, FILTER= FILTER_n, Meta=Meta_n) -> Association_results_p
              Null_distribution = c(Null_distribution, Association_results_p$`Pr(>|t|)_norm`)
        }

      sapply( Summary_statistics_taxonomyX$`Pr(>|t|)_norm`, function(x){ sum( Null_distribution <= x  )/Permutation_number  } ) -> FDR_perm
      FDR_perm[FDR_perm > 1] = 1
      Summary_statistics_taxonomyX %>% mutate(FDR_permutation = FDR_perm ) -> Summary_statistics_taxonomyX
      print("Permutations succesful" )
    }
    Summary_statistics = rbind ( Summary_statistics, mutate(Summary_statistics_taxonomyX, Taxonomic_level = Taxonomic_level )   )
  }
  return(Summary_statistics)
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

Run_Associations(Prevalence, All_taxonomy, Data_c, FILTER=10, Function=Association_analysis, Run_permutations = T) -> Association_results

###2.2.3 Add info from balances
Balances_results = tibble(Bug = unique(c(Balances[[1]], Balances[[3]]) ), Direction_balance = "Positive")
Balances_results = rbind(Balances_results, tibble(Bug = unique(c(Balances[[2]], Balances[[4]] ) ), Direction_balance = "Negative") )
left_join(Association_results,Balances_results) -> Association_results

###2.2.4 Plot
Association_results %>% ggplot(aes(x=Estimate_norm, y= -log10(`Pr(>|t|)_norm`), col=FDR_permutation<0.05, shape= is.na(Direction_balance) )) +
  geom_hline( yintercept = -log10(0.05), color="red" ) + geom_point() + theme_bw() + facet_wrap(~Taxonomic_level) +
  geom_label_repel(data = Top_results %>% filter(FDR_permutation<0.05) %>% unique(),  aes(label = Bug), size=1.5, color = "black")


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

left_join(Association_results, Aldex_results, by = c("Bug", "Taxonomic_level")) %>% left_join(., Prevalence) %>%
    write_tsv( ., "Results/Summary_statistics_TaxonomyAnalysis.csv")  



Old_pathway = function(){

  ######################
  ## Pathway analysis###
  ######################
  
  #Get prevalences
  Comp_Prevalence( select(Case_ptw, -ID) ) -> Prevalence_cases_ptw
  Comp_Prevalence( select(Control_ptw, -ID) ) -> Prevalence_controls_ptw
  full_join(Prevalence_cases_ptw, Prevalence_controls_ptw, by="Bug", suffix = c("_case", "_control") ) %>% drop_na() -> Prevalence_ptw
  
  #Prepare data
  rbind( select(Control_ptw, c(ID, Prevalence_ptw$Bug) ), select(Case_ptw, c(ID, Prevalence_ptw$Bug) ))  -> Data_ptw
  Data_ptw %>% select(-ID) %>% apply(2, as.numeric) %>% as.data.frame() %>% mutate(ID = Data_ptw$ID) -> Data_ptw
  
  Data_ptw %>% mutate(Cohort = ifelse(grepl("LL", ID), "Control", "Case" )) -> Data_c_ptw
  Data_c_ptw$Cohort = factor(Data_c_ptw$Cohort, levels= c("Control", "Case"))
  
  
  #Prediction## Using only fully prevalent taxa
  Prevalence_ptw %>% filter(Prevalence_case == 1 & Prevalence_control == 1) %>% filter(! Bug %in% c("UNMAPPED", "UNINTEGRATED")) -> Keep_ptw
  Run_balance_analysis( Balance_input = select(Data_ptw, c("ID", Keep_ptw$Bug)) , Labels = Data_c_ptw$Cohort, lambda= 1, add_pseudo = F  ) -> Balances_ptw
  Data_ptw %>% select(Balances_ptw[[1]]) %>% apply(1, sum) -> Numerator_ptw ; Data_ptw %>% select(Balances_ptw[[3]]) %>% apply(1, sum) -> Numerator2_ptw ; Data_ptw %>% select(Balances_ptw[[2]]) %>% apply(1, sum) -> Denominator_ptw ; Data_ptw %>% select(Balances_ptw[[4]]) %>% apply(1, sum) -> Denominator2_ptw
  Data_ptw %>% mutate( Balance = log10(Numerator_ptw/Denominator_ptw), Balance2 = log10(Numerator2_ptw/Denominator2_ptw), Cohort = Data_c_ptw$Cohort ) -> Data_b_ptw
  Data_b_ptw %>% ggplot(aes(x=Cohort, y= Balance, col=Cohort))  + geom_violin() + theme_bw() + coord_flip() + geom_jitter() + stat_summary(fun = "median",geom = "crossbar", aes(color = Cohort)) + stat_summary(fun = "mean", geom = "point", color= "black") 
  Data_b_ptw %>% ggplot(aes(x=Cohort, y= Balance2, col=Cohort))  + geom_violin() + theme_bw() + coord_flip() + geom_jitter() + stat_summary(fun = "median",geom = "crossbar", aes(color = Cohort)) + stat_summary(fun = "mean", geom = "point", color= "black")
  
  #Linear model
  Compute_CLR_taxonomy(Data_ptw, Keep_column = c("UNMAPPED", "UNINTEGRATED")) -> Data_ptw2
  Association_analysis( Prevalence_ptw, Data_ptw2, Data_c_ptw, Meta=Read_number ,FILTER=30) -> Association_results_ptw
  
  #Check if association results are reproduced in the Balance analysis 
  Balances_results_ptw = tibble(Bug = unique(c(Balances_ptw[[1]], Balances_ptw[[3]])), Direction_balance = "Positive")
  Balances_results_ptw = rbind(Balances_results_ptw, tibble(Bug = unique(c(Balances_ptw[[2]], Balances_ptw[[4]])), Direction_balance = "Negative") )
  #2/4 bugs in balances are nominally significant
  left_join(Association_results_ptw,Balances_results_ptw) -> Association_results_ptw
  
  #Permutations
  Null_distribution_ptw= c()
  for (Permutation_n in seq(Permutation_number) ){
    Data_perm = Data_c_ptw
    Data_perm$Cohort = sample(Data_c_ptw$Cohort, replace = F,size = dim(Data_c_ptw)[1] )
    Association_analysis(Prevalence_ptw, Data_ptw2, Data_perm, FILTER=30) -> Association_results_p
    Null_distribution = c(Null_distribution_ptw, Association_results_p$`Pr(>|t|)_norm`)
  }
  
  sapply( Association_results_ptw$`Pr(>|t|)_norm`, function(x){ sum( Null_distribution <= x  )/Permutation_number  } ) -> FDR_perm_ptw
  FDR_perm_ptw[FDR_perm_ptw > 1] = 1
  Association_results_ptw %>% mutate(FDR_permutation = FDR_perm_ptw ) -> Association_results_ptw
  
  
  Association_results_ptw %>% ggplot(aes(x=Estimate_norm, y= -log10(`Pr(>|t|)_norm`), shape= is.na(Direction_balance) )) +
    geom_hline( yintercept = -log10(0.05), color="red" ) + geom_point() + theme_bw() +
    geom_label_repel(data = Top_results %>% filter(FDR_permutation<0.05) %>% unique(), 
                     aes(label = Bug),
                     color = "black")
  write_tsv(Association_results_ptw, "Summary_statistics_PathwayAnalysis.csv")  
}



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


################################
#Bayesian regression analysis##
###############################

packageurl <- "https://cran.r-project.org/src/contrib/Archive/fido/fido_1.0.0.tar.gz"
install.packages(packageurl, repos=NULL, type="source")

library(fido)
library(phyloseq)

Fit_pibble_bayesian =  function(Data, Read_number, Prevalence,Plot_n, Prev_thresho=20, Taxonomy_level = "s"){
  if (! Taxonomy_level == "Path"){
    paste( c(Taxonomy_level, "__"), collapse="" ) -> Ta
    colnames(dplyr::select(Data, -c(ID, UNKNOWN ))) -> Taxa
    Taxa[grepl(Ta, Taxa)] -> Taxa
    Data %>% dplyr::select(Taxa)  -> Data_taxonomy
  } else { Data_taxonomy = dplyr::select(Data, -c(ID,UNMAPPED, UNINTEGRATED)) }
  #Make data in proportions
  Y = as.data.frame(Data_taxonomy)
  rownames(Y) = Data$ID
  apply(Y, 1, function(x){ x/sum(x) } ) %>% t() -> Y
  Prevalence %>% filter(N_case > Prev_thresho & N_control > Prev_thresho)  -> To_Test
  #Seems that this one cannot work in low prevalence bugs... (it needs samples as rows)
  Y <- zCompositions::cmultRepl(as.data.frame(Y)  %>% dplyr::select(one_of(To_Test$Bug))) # multiplicative 0 replacement 
  Y = t(Y)
  #Remove Nor prevalent taxa
  #Prevalence %>% filter(N_case > Prev_thresho & N_control > Prev_thresho)  -> To_Test #I had to do this so that cmultRepl would work
  #Y[rownames(Y) %in% To_Test$Bug, ] -> Y

  
  #Covariates
  left_join(dplyr::select(Data, ID),  Read_number) -> Covariates
  Covariates %>% mutate(Cohort = factor(Cohort, levels = c("Control", "Case")))  -> Covariates
  
  # your linear model
  f <- as.formula( " ~ Cohort + Read_number ")
  # construct a design/model matrix
  X <- t(model.matrix(f, data=Covariates))


  ##-----------------------------------------------------------------------------------------------------------------------
  # Because pi_j = alrInv(eta_j) it also holds that eta_j = alr(pi_j) (see Pibble documentation)
  # pi_j denote multinomial category j which exist in the simplex, and eta_j denote the transformed category j which exists in real space  
  ##-----------------------------------------------------------------------------------------------------------------------
  N <- ncol(Y) # number of samples
  D <- nrow(Y) # number of multinomial categories (i.e. features)

 #Eta are the transformed data abundances. Give einitial values
  eta_init <- t(driver::alr(t(Y))) #What does alr divide for by default?
  eta_array <- array(eta_init, dim=c(nrow(eta_init), ncol(eta_init), 2000)) # needs to be an array of dimension (D-1) x N x iter 
 #Priors
  upsilon <- D+3
  Omega <- diag(D)
  G <- cbind(diag(D-1), -1)
  Xi <- (upsilon-D)*G%*%Omega%*%t(G)
  Theta <- matrix(0, D-1, nrow(X))
  Gamma <- diag(nrow(X))
  #Prior model
  priors <- pibble(NULL, X, upsilon, Theta, Gamma, Xi)
  #priors2 <- to_clr(priors) 
  #Make some checks on the priors
  #fido::summary(priors2, pars="Lambda")
  #names_covariates(priors2) <- rownames(X)
  #fido::plot(priors, par="Lambda") + ggplot2::xlim(c(-10, 10))
  
  #Posterior model in logistic normal (no multinomial)
  posterior <- uncollapsePibble(eta_array, priors$X, priors$Theta, priors$Gamma, priors$Xi, priors$upsilon, seed=4302) # # uncollapsePibble can be fitted on proportions 

  # Attach dimnames
  dimnames(posterior$Lambda)[[2]] <- rownames(X)
  dimnames(posterior$Lambda)[[1]] <- rownames(Y)[-length(rownames(Y))]
  dimnames(posterior$Sigma)[[1]] <- dimnames(posterior$Sigma)[[2]] <- rownames(Y)[-length(rownames(Y))]

  posterior <- pibblefit(D=D,
                       N=N,
                       Q=nrow(X),
                       coord_system="alr",
                       iter=2000L,
                       alr_base=D,
                       Eta=eta_array,
                       Lambda=posterior$Lambda,
                       Sigma=posterior$Sigma,
                       Y=Y,
                       X=X,
                       names_categories=rownames(Y),
                       names_samples=colnames(Y),
                       names_covariates=rownames(X))
  #QC, not sure if it works
  #ppc_summary(posterior)
  
  fido::summary(posterior, pars="Lambda")$Lambda  %>% filter(covariate == "CohortCase") %>% mutate( Significant = ifelse(  sign(p2.5) * sign(p97.5) == 1, "yes", "no" )   ) %>% arrange(desc(abs(mean)))  -> Result_bay
  
  focus <- Result_bay[sign(Result_bay$p2.5) == sign(Result_bay$p97.5),]
  if(! length(focus) == 0){
    focus <- unique(focus$coord)
    fido::plot(posterior, par="Lambda", focus.coord = focus, focus.cov = rownames(X)[2])
  }
  # Change to CLR transform
  posterior_clr <- to_clr(posterior)
  # Attached dimnames
  dimnames(posterior_clr$Lambda)[[2]] <- rownames(X)
  dimnames(posterior_clr$Lambda)[[1]] <- rownames(Y)
  dimnames(posterior_clr$Sigma)[[1]] <- dimnames(posterior_clr$Sigma)[[2]] <- rownames(Y)

  fido::summary(posterior_clr, pars="Lambda")$Lambda %>% filter(covariate == "CohortCase") %>% mutate( Significant = ifelse(  sign(p2.5) * sign(p97.5) == 1, "yes", "no" )   ) -> Result_bay_clr

  focus <- Result_bay_clr[sign(Result_bay_clr$p2.5) == sign(Result_bay_clr$p97.5),]
  if (! length(focus) == 0){
    focus <- unique(focus$coord)
    fido::plot(posterior_clr, par="Lambda", focus.coord = focus, focus.cov = rownames(X)[2]) + geom_vline(xintercept = 0) -> To_plot
    To_plot %>% print()
    ggsave(Plot_n, To_plot)
  }
  
  return(Result_bay_clr)
}


Bayesian_results = tibble()
for (i in c("p", "c", "o", "f", "g", "s")){
  Plot_name = paste(c("Plots/Bayes/Taxonomy_", i, ".pdf" ), collapse="")
  Fit_pibble_bayesian(Data, Read_number, Prevalence, Prev_thresho = 20, Taxonomy_level = i, Plot_n = Plot_name ) -> Results_b
  Bayesian_results = rbind(Bayesian_results, Results_b)
  
}  
write_tsv(Bayesian_results,"Summary_statistics_TaxonomyAnalysisBayesian.csv")

Plot_name = "Plots/Bayes/Taxonomy_pathway.pdf" 
Fit_pibble_bayesian(Data_ptw, Read_number, Prevalence_ptw, Prev_thresho = 20, Plot_n = Plot_name, Taxonomy_level = "Path")



Covariates_GUTS 
Covariates_control

Covariates_all %>% mutate(Age = scale(Age), BMI = scale(BMI) ) %>% select(Age, Sex, BMI) %>% dist(method = "euclidean") -> D 
All_matches = tibble()
for (i in  1:104){
  as.matrix(D)[i,]  -> z
  which.min(z[104:length(z)]) + 104 -> Answer
  Matches = tibble(N1 =Covariates_all$Sequencing_ID[i] , N2 =Covariates_all$Sequencing_ID[Answer], Distance=min(z[104:length(z)])  )
  rbind(All_matches, Matches) -> All_matches
}



read_csv("/Users/sergio/Documents/GitHub/GUTS_MGS_Psych/Data/MetaData/Stats_knead.csv") -> KN
KN %>% filter(Sample == "LL10_B11_856" )

colnames(KN)[1] = "Sequencing_ID"
left_join(KN, Covariates_all) -> KN

KN %>% ggplot(aes(x=`Final (p1)`, y=Read_number)) + 
  geom_point() + facet_wrap(~Batch)
KN  %>% mutate(Batch = ifelse(grepl("dag3",Batch), "DMP", Batch ) ) %>% filter(!is.na(Batch)) %>%
  ggplot(aes(x=log10(`Final (p1)`), fill=Batch )) + geom_density(alpha=0.4) + scale_fill_manual(values=met.brewer("Manet", 5)[c(1,5, 3)] )
