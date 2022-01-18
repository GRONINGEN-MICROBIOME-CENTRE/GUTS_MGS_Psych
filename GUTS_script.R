setwd("~/Documents/PhD/Psychiatry/DATA/Data/")
library(tidyverse)
library(vegan)
library(ape)
library(MetBrewer) #Lets plot with some art :)
library(ggrepel)
ggplot <- function(...) ggplot2::ggplot(...) + scale_color_manual(values=met.brewer("Manet",5)[c(1,5)] ) +
  scale_fill_manual(values=met.brewer("Manet", 5)[c(1,5)] )

Clean_path = function(C, Keep){
  colnames(C)[1] = "ID"
  sapply(colnames(C), function(x){ str_split(x, "_knead")[[1]][1] } ) -> New_names
  colnames(C) = New_names
  C %>% select(c("ID", Keep$ID)) -> C
  return(C)
}
Format_path = function(Path){
  Path %>% filter(! grepl("#", ID)) -> Path
  Path %>% filter(! grepl("\\|",ID )) -> Path
  t(Path) %>% as.data.frame() %>% rownames_to_column("ID") %>% as_tibble() %>% filter(! ID == "ID") -> Path_2
  colnames(Path_2) = c("ID",  Path$ID)
  return(Path_2)
}

Case = read_tsv("Cases.tsv")
Control = read_tsv("Matched_controls.tsv")

#Pathway abundances
Case_ptw = Format_path(Clean_path(read_tsv("merged_Pathway_table_cases.txt"), Keep = Case ))
Control_ptw = Format_path(Clean_path(read_tsv("merged_Pathway_table_controls.txt"), Keep = Control)) 

#Covariate Read number
Case_readN = read_tsv("Read_number.tsv") %>%  mutate(Cohort = "Case") 
Control_readN = read_tsv("DAG3_readnumber.tsv") %>% mutate(Cohort = "Control") %>% filter(ID %in% Control$ID)
rbind(Case_readN, Control_readN) -> Read_number
Read_number %>% ggplot(aes(x=Read_number, fill=Cohort )) + geom_histogram(alpha=0.5, position="identity") + theme_bw()
#1. Compute prevalences
Comp_Prevalence = function(DF){
  #Compute prevalence of different taxa
  apply(DF,2, function(x){ sum(!as.numeric(x)==0)/length(x) } ) -> Prevalence_vector
  apply(DF,2, function(x){ sum(!as.numeric(x)==0) } ) -> N
  Prevalence_df = tibble(Bug = names(Prevalence_vector), Prevalence = Prevalence_vector, N = N)
  return(Prevalence_df)
}
Comp_Prevalence( select(Case, -ID) ) -> Prevalence_cases
Comp_Prevalence( select(Control, -ID) ) -> Prevalence_controls
#Not the same taxonomy detected, so we keep only taxa seen in both
full_join(Prevalence_cases, Prevalence_controls, by="Bug", suffix = c("_case", "_control") ) %>% drop_na() -> Prevalence
Prevalence %>% ggplot(aes(x=Prevalence_case, y=Prevalence_control)) + geom_point() + theme_bw()+ geom_abline()
rbind( select(Control, c(ID, Prevalence$Bug) ), select(Case, c(ID, Prevalence$Bug) ))  -> Data
#Add annotation about which cohort does it belong to
Data %>% mutate(Cohort = ifelse(grepl("LL", ID), "Control", "Case" )) -> Data_c
Data_c$Cohort = factor(Data_c$Cohort, levels= c("Control", "Case"))

#Check if proportion of unknown is different between cohorts
summary(lm(UNKNOWN ~ Cohort, Data_c)) #On average, there seems to be a bit more of unknown in controls.

#Transformation of the data. Microbiome data is compositional, and thus, log-ratios might be used to get all microbial features in the same reference framework
Pseudocount = function(Taxa){
  return(min(Taxa[!Taxa==0])/2)
}
Geom_mean = function(x){
  exp(mean(log(x)))
}
CLR = function(D){
  log(D / Geom_mean(D) )
  
}
Compute_CLR_taxonomy = function(Data, Taxonomy = "s", Keep_column = "UNKNOWN"){
  if ("Taxonomy" == "Path"){
    paste( c(Taxonomy, "__"), collapse="" ) -> Ta
    colnames(select(Data, -c("ID", Keep_column ))) -> Taxa
    Taxa[grepl(Ta, Taxa)] -> Taxa
    Data %>% select( c("ID", Keep_column, Taxa) ) -> Data_taxonomy
  } else { Data_taxonomy = Data}
  #computation of a different pseudocoount per taxa.
  apply(select(Data_taxonomy, -ID), 2, function(taxa){
    if( sum(taxa==0) >= 1){ pscount =  Pseudocount(taxa) 
    } else { pscount = 0 }
    pscount } ) -> Pseudocounts
  apply(select(Data_taxonomy, -ID), 1, function(Participant){ 
    P = CLR(Participant + Pseudocounts)
    P
  }) %>% t() %>% as_tibble() %>% mutate(ID = Data_taxonomy$ID, .before=1) -> Data_taxonomy
  Data_taxonomy %>% select(-Keep_column) -> Data_taxonomy
  return(Data_taxonomy)
}  

Compute_CLR_taxonomy(Data, Taxonomy = "s") -> Species_data
All_taxonomy = Species_data
for (i in c("p", "c", "o", "f", "g")){
  Compute_CLR_taxonomy(Data, Taxonomy = i) -> Transformed_data
  left_join(All_taxonomy, Transformed_data, by = "ID") -> All_taxonomy
}    


#Compute beta diversity at the species level
Beta_taxonomy = function(Data_t, Data_c, Meta = Read_number){
  #Compute diversity
  #Euclidean distance in CLR-transformed data corresponds to Aitchison distance.
  vegdist(select(Data_t, -ID), method = "euclidean") -> Beta_diversity
  
  #A PCOA
  pcoa(Beta_diversity) -> summary_pcoa
  PCs = as_tibble(summary_pcoa$vectors)
  Variab = as_tibble(summary_pcoa$values %>% rownames_to_column("PC") )
  PCs %>% mutate(ID = Data_c$ID, Cohort=Data_c$Cohort) -> PCs
  PCs <- merge(PCs,aggregate(cbind(mean.x=Axis.1, mean.y=Axis.2)~Cohort,PCs,mean),by="Cohort")
  PCs %>% ggplot(aes(x=Axis.1, y=Axis.2, col=Cohort)) + geom_point() +
  xlab(paste( c("PC1 (", as.character(100*round(Variab$Relative_eig[1],2)), "%)" ), collapse="" ) ) +
  ylab(paste( c("PC2 (", as.character(100*round(Variab$Relative_eig[2],2)), "%)" ), collapse="" ) ) + 
    theme_bw() + stat_ellipse() +
    geom_point(aes(x=mean.x,y=mean.y),size=5)+
    geom_segment(aes(x=mean.x, y=mean.y, xend=Axis.1, yend=Axis.2), alpha=0.5) -> Fig
  print(Fig)
  
  #B PERMANOVA
  left_join(Data_c, Meta) -> Data_m
  adonis2(Beta_diversity ~ Data_m$Read_number  + Data_m$Cohort, permutations = 5000) %>% print()
  
}
Beta_taxonomy(Species_data, Data_c)

#Compute alpha diversity at the species level ; need to do this in non-transformed data
Alpha_taxonomy = function(Data, Data_c, Taxonomy = "s", Meta = Read_number){
  paste( c(Taxonomy, "__"), collapse="" ) -> Ta
  colnames(select(Data, -c(ID, UNKNOWN ))) -> Taxa
  Taxa[grepl(Ta, Taxa)] -> Taxa
  Data %>% select(Taxa)  -> Data_taxonomy
  diversity(Data_taxonomy, "shannon" ) -> Diversity
  Data_c %>% select(ID, Cohort) %>% mutate(Diversity = Diversity) -> Data_div
  left_join(Data_div, Meta) -> Data_div
  summary(lm(Diversity ~ Cohort + Read_number,Data_div)) %>% print()
  ggplot(Data_div, aes(x=Cohort, y=Diversity, col=Cohort))  + geom_violin() +
  theme_bw() + coord_flip() + geom_jitter() +
    stat_summary(fun = "median",geom = "crossbar", aes(color = Cohort)) +
    stat_summary(fun = "mean", geom = "point", color= "black")
  
}
Alpha_taxonomy(Data, Data_c, "s")


#Use coda core to predict Cohort
library(codacore)
library(tensorflow)
library(pROC)
Run_balance_analysis = function(Balance_input, Labels = Data_c$Cohort, lambda=1, add_pseudo = T){
  set.seed(50)
  Balance_input %>% select(-ID) -> Balance_input2
  if (add_pseudo == T){
    Pseudo = apply(Balance_input2, 2, Pseudocount)
    Balance_input2 = as.data.frame(Balance_input2) + Pseudo
  }
  #Split dataset
  tf$random$set_seed(0)
  trainIndex <- sample(1:nrow( Balance_input2), 0.8 * nrow(Balance_input2))
  #Train set
  xTrain <- Balance_input2[trainIndex,]
  yTrain <- as.factor(Labels[trainIndex])
  #Test set
  xTest <- Balance_input2[-trainIndex,]
  yTest <-as.factor(Labels[-trainIndex])
  
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
    
  }else { return(list(Numerators, Denominators)) }
}   

#Get only prevalent bacteria for the model
Prevalence %>% filter(N_case > 20 & N_control > 20) -> Keep
Run_balance_analysis( select(Data, c("ID", Keep$Bug)) ) -> Balances
Data %>% select(Balances[[1]]) %>% apply(1, sum) -> Numerator ; Data %>% select(Balances[[3]]) %>% apply(1, sum) -> Numerator2
Data %>% select(Balances[[2]]) %>% apply(1, sum) -> Denominator ; Data %>% select(Balances[[4]]) %>% apply(1, sum) -> Denominator2
Data %>% mutate( Balance = log10(Numerator/Denominator), Balance2 = log10(Numerator2/Denominator2), Cohort = Data_c$Cohort ) ->Data_b
Data_b %>%
  ggplot(aes(x=Cohort, y= Balance, col=Cohort))  + geom_violin() + theme_bw() + coord_flip() + geom_jitter() +
  stat_summary(fun = "median",geom = "crossbar", aes(color = Cohort)) + stat_summary(fun = "mean", geom = "point", color= "black")
summary(lm(Balance ~ Cohort, filter(Data_b, ! abs(Balance) == Inf )))
Data_b %>%
  ggplot(aes(x=Cohort, y= Balance2, col=Cohort))  + geom_violin() + theme_bw() + coord_flip() + geom_jitter() +
  stat_summary(fun = "median",geom = "crossbar", aes(color = Cohort)) + stat_summary(fun = "mean", geom = "point", color= "black")
summary(lm(Balance2 ~ Cohort, filter(Data_b, ! abs(Balance2) == Inf )))

#Association analysis
Inverse_rank_normal = function(Measurement){
  qnorm((rank(Measurement,na.last="keep")-0.5)/sum(!is.na(Measurement)))
}
Association_analysis = function(Prevalence, DF, Data_c, FILTER=20, Meta=Read_number){
  Total_results = tibble()
  Prevalence %>% filter(N_case > FILTER & N_control > FILTER) %>% filter(! Bug == "UNKNOWN") -> To_Test
  for (Bug in To_Test$Bug){
    if (! Bug %in% colnames(DF) ){ next }
    DF %>% select( one_of(Bug) ) %>% as_vector() -> vector_Bug
    #IF we want to normalize it, apply function here
    normalized_Bug = Inverse_rank_normal(vector_Bug)
    #invers rank normal transf?
    Data_c %>% select(ID, Cohort) %>% mutate(B =  vector_Bug, B_n = normalized_Bug) -> Model_input
    left_join(Model_input, Meta) -> Model_input
    Model_input %>% mutate(Cohort = factor(Cohort,  levels = c("Control", "Case") )) -> Model_input
    lm(B ~ Cohort+ Read_number, Model_input  ) -> Model_out
    lm(B_n ~ Cohort+ Read_number, Model_input  ) -> Model_out_n
    #If not normalized
    Normality = shapiro.test(Model_out$residuals)$p.value
    as.data.frame(summary(Model_out)$coefficients)["CohortCase",] %>% as_tibble() %>%
      mutate(Bug = Bug, Shapiro_p = Normality, .before=1) -> results
    #Normalized
    Normality = shapiro.test(Model_out_n$residuals)$p.value
    as.data.frame(summary(Model_out_n)$coefficients)["CohortCase",] %>% as_tibble() %>%
      mutate(Bug = Bug, Shapiro_p = Normality, .before=1) -> Normalized_results
    
    left_join(results, Normalized_results, by= "Bug", suffix = c("","_norm")) -> results
    rbind(Total_results, results) -> Total_results
}
  return(Total_results)
}  

#Just on species
Association_analysis(Prevalence, Species_data, Data_c, FILTER=10) -> Association_results
#In all taxonomy
Association_analysis(Prevalence, All_taxonomy, Data_c, FILTER=10) -> Association_results2

#Compare the distribution of Ps, in both cases it seems alright
ggplot() + geom_histogram( aes(x=Association_results$`Pr(>|t|)_norm`), alpha = 0.3, fill="blue") + 
  geom_histogram( aes(x=Association_results2$`Pr(>|t|)_norm`), alpha = 0.3, fill="green") +  theme_bw()
  
#Compare results normalized vs no normalized
Association_results2 %>% ggplot(aes(x = -log10(`Pr(>|t|)`), y = -log10(`Pr(>|t|)_norm`) )) + 
  theme_bw() + geom_point() + geom_abline() + geom_vline(xintercept = -log10(0.05), color = "red") + geom_hline(yintercept = -log10(0.05), color = "red")
Association_results2 %>% arrange(`Pr(>|t|)_norm`) %>% select(-c(`Std. Error`, `t value`, `Std. Error_norm`, `t value_norm` )) -> Top_results
#Check if association results are reproduced in the Balance analysis 
Balances_results = tibble(Bug = unique(c(Balances[[1]], Balances[[3]])), Direction_balance = "Positive")
Balances_results = rbind(Balances_results, tibble(Bug = unique(c(Balances[[2]], Balances[[4]])), Direction_balance = "Negative") )
# 42 nominally significant results, 8/16 bugs used in balances
left_join(Top_results,Balances_results) -> Top_results
Top_results %>% select(Bug, Estimate_norm, `Pr(>|t|)_norm`, Direction_balance ) -> Associated


#Compute FDR using label permutations
Null_distribution = c()
Permutation_number = 100
for (Permutation_n in seq(Permutation_number) ){
  Data_perm = Data_c
  Data_perm$Cohort = sample(Data_c$Cohort, replace = F,size = dim(Data_c)[1] )
  Association_analysis(Prevalence, All_taxonomy, Data_perm, FILTER=10) -> Association_results_p
  Null_distribution = c(Null_distribution, Association_results_p$`Pr(>|t|)_norm`)
}

sapply( Top_results$`Pr(>|t|)_norm`, function(x){ sum( Null_distribution <= x  )/Permutation_number  } ) -> FDR_perm
FDR_perm[FDR_perm > 1] = 1
Top_results %>% mutate(FDR_permutation = FDR_perm ) -> Top_results
Top_results %>% arrange(`Pr(>|t|)_norm`) %>% select(`Pr(>|t|)_norm`, FDR_permutation)

Top_results %>% ggplot(aes(x=Estimate_norm, y= -log10(`Pr(>|t|)_norm`), col=FDR_permutation<0.05, shape= is.na(Direction_balance) )) +
  geom_hline( yintercept = -log10(0.05), color="red" ) + geom_point() + theme_bw() +
  geom_label_repel(data = Top_results %>% filter(FDR_permutation<0.05) %>% unique(), 
                   aes(label = Bug),
                   color = "black")
write_tsv(Top_results, "Summary_statistics_TaxonomyAnalysis.csv")  


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
