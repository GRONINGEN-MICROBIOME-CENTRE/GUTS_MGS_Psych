setwd("~/Documents/PhD/Psychiatry/DATA/Data/")
library(tidyverse)
library(vegan)
library(ape)

Case = read_tsv("Cases.tsv")
Control = read_tsv("Matched_controls.tsv")

#1. Compute prevalences
Comp_Prevalence = function(DF){
  #Compute prevalence of different taxa
  apply(DF,2, function(x){ sum(!x==0)/length(x) } ) -> Prevalence_vector
  apply(DF,2, function(x){ sum(!x==0) } ) -> N
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
Data %>% mutate(Cohort = ifelse(grepl("LL", ID), "Case", "Control" )) -> Data_c

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
Compute_CLR_taxonomy = function(Data, Taxonomy = "s"){
  paste( c(Taxonomy, "__"), collapse="" ) -> Ta
  colnames(select(Data, -c(ID, UNKNOWN ))) -> Taxa
  Taxa[grepl(Ta, Taxa)] -> Taxa
  Data %>% select( c("ID","UNKNOWN", Taxa) ) -> Data_taxonomy
  
  #computation of a different pseudocoount per taxa.
  apply(select(Data_taxonomy, -ID), 2, function(taxa){
    if( sum(taxa==0) >= 1){ pscount =  Pseudocount(taxa) 
    } else { pscount = 0 }
    pscount } ) -> Pseudocounts
  apply(select(Data_taxonomy, -ID), 1, function(Participant){ 
    P = CLR(Participant + Pseudocounts)
    P
  }) %>% t() %>% as_tibble() %>% mutate(ID = Data_taxonomy$ID, .before=1) -> Data_taxonomy
  Data_taxonomy %>% select(-UNKNOWN) -> Data_taxonomy
  return(Data_taxonomy)
}  

Compute_CLR_taxonomy(Data, Taxonomy = "s") -> Species_data
All_taxonomy = Species_data
for (i in c("p", "c", "o", "f", "g")){
  Compute_CLR_taxonomy(Data, Taxonomy = i) -> Transformed_data
  left_join(All_taxonomy, Transformed_data, by = "ID") -> All_taxonomy
}    


#Compute beta diversity at the species level
Beta_taxonomy = function(Data_t, Data_c){
  #Compute diversity
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
  
  
  #B PERMANOVA
  adonis2(Beta_diversity ~ Data_c$Cohort, permutations = 5000) %>% print()
  
}
Beta_taxonomy(Species_data, Data_c)

#Compute alpha diversity at the species level ; need to do this in non-transformed data
Alpha_taxonomy = function(Data, Data_c, Taxonomy = "s"){
  paste( c(Taxonomy, "__"), collapse="" ) -> Ta
  colnames(select(Data, -c(ID, UNKNOWN ))) -> Taxa
  Taxa[grepl(Ta, Taxa)] -> Taxa
  Data %>% select(Taxa)  -> Data_taxonomy
  diversity(Data_taxonomy, "shannon" ) -> Diversity
  Data_c %>% select(ID, Cohort) %>% mutate(Diversity = Diversity) -> Data_div
  summary(lm(Diversity ~ Cohort,Data_div)) %>% print()
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
Run_balance_analysis = function(Balance_input, Labels = Data_c$Cohort, lambda=1){
  set.seed(50)
  Balance_input %>% select(-ID) -> Balance_input2 
  Pseudo = apply(Balance_input2, 2, Pseudocount)
  Balance_input2 = as.data.frame(Balance_input2) + Pseudo
  
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
  cat("Train set AUC =", codacoreAUC, "\n")

  ##Test
  yHat <- predict(model, newx = xTest)
  print("Test results")
  testAUC <- pROC::auc(pROC::roc(yTest, yHat, quiet=T))
  cat("Test set AUC =", testAUC, "\n")
  
  plot(model)
  plotROC(model)
  Numerators = colnames(xTrain)[getNumeratorParts(model, 1)]
  Denominators = colnames(xTrain)[getDenominatorParts(model, 1)]
  
  return(list(Numerators, Denominators))
}   
Run_balance_analysis(Data) -> Balances
Data %>% select(Balances[[1]]) %>% apply(1, sum) -> Numerator
Data %>% select(Balances[[2]]) %>% apply(1, sum) -> Denominator
Data %>% mutate( Balance = log10(Numerator/Denominator), Cohort = Data_c$Cohort ) ->Data_b
Data_b %>%
  ggplot(aes(x=Cohort, y= Balance, col=Cohort))  + geom_violin() + theme_bw() + coord_flip() + geom_jitter() +
  stat_summary(fun = "median",geom = "crossbar", aes(color = Cohort)) + stat_summary(fun = "mean", geom = "point", color= "black")
summary(lm(Balance ~ Cohort, filter(Data_b, ! abs(Balance) == Inf )))

#Association analysis
Inverse_rank_normal = function(Measurement){
  qnorm((rank(Measurement,na.last="keep")-0.5)/sum(!is.na(Measurement)))
}
Association_analysis = function(Prevalence, DF, Data_c, FILTER=20){
  Total_results = tibble()
  Prevalence %>% filter(N_case > FILTER & N_control > FILTER) %>% filter(! Bug == "UNKNOWN") -> To_Test
  for (Bug in To_Test$Bug){
    if (! Bug %in% colnames(DF) ){ next }
    DF %>% select( one_of(Bug) ) %>% as_vector() -> vector_Bug
    #IF we want to normalize it, apply function here
    normalized_Bug = Inverse_rank_normal(vector_Bug)
    #invers rank normal transf?
    Data_c %>% select(ID, Cohort) %>% mutate(B =  vector_Bug, B_n = normalized_Bug) -> Model_input
    lm(B ~ Cohort, Model_input ) -> Model_out
    lm(B_n ~ Cohort, Model_input ) -> Model_out_n
    #If not normalized
    Normality = shapiro.test(Model_out$residuals)$p.value
    as.data.frame(summary(Model_out)$coefficients)["CohortControl",] %>% as_tibble() %>%
      mutate(Bug = Bug, Shapiro_p = Normality, .before=1) -> results
    #Normalized
    Normality = shapiro.test(Model_out_n$residuals)$p.value
    as.data.frame(summary(Model_out_n)$coefficients)["CohortControl",] %>% as_tibble() %>%
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
Balances_results = tibble(Bug = Balances[[1]], Direction = "Positive")
Balances_results = rbind(Balances_results, tibble(Bug = Balances[[2]], Direction = "Negative") )
Top_results %>% select(Bug, Estimate_norm, `Pr(>|t|)_norm` ) -> Associated
#seems that most of the bacteria used for the ratios are low prevalence (maybe batch?)
left_join(Associated,Balances_results) -> Assoc_Balance_info


