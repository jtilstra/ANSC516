#set working directory
list.files()


library(tidyverse)
library(qiime2R)
library(ggpubr)

##Load Data
meta<-read_q2metadata("Proj_metadata.txt")
str(meta)

#create variables for each of the alpha diversity metrics
evenness = read_qza("core-metrics-results/evenness_vector.qza")
evenness<-evenness$data %>% rownames_to_column("SampleID") # this moves the sample names to a new column that matches the metadata and allows them to be merged

observed_features = read_qza("core-metrics-results/observed_features_vector.qza")
observed_features<-observed_features$data %>% rownames_to_column("SampleID") # this moves the sample names to a new column that matches the metadata and allows them to be merged

shannon = read_qza("core-metrics-results/shannon_vector.qza")
shannon<-shannon$data %>% rownames_to_column("SampleID") # this moves the sample names to a new column that matches the metadata and allows them to be merged

faith_pd = read_qza("core-metrics-results/faith_pd_vector.qza")
faith_pd<-faith_pd$data %>% rownames_to_column("SampleID") # this moves the sample names to a new column that matches the metadata and allows them to be merged\

## Clean up the data

str(meta)
#observed_features$observed_features_num <- lapply(observed_features$observed_features, as.numeric)
observed_features$observed_features <- as.numeric(observed_features$observed_features)
str(observed_features)
shannon$shannon <-as.numeric(shannon$shannon)
str(shannon)


###Alpha Diversity tables
# These tables will be merged for convenience and added to the 
# metadata table as the original tutorial was organized.

alpha_diversity = merge(x=faith_pd, y=evenness, by.x = "SampleID", by.y = "SampleID")
alpha_diversity = merge(alpha_diversity, observed_features, by.x = "SampleID", by.y = "SampleID")
alpha_diversity = merge(alpha_diversity, shannon, by.x = "SampleID", by.y = "SampleID")
meta = merge(meta, alpha_diversity, by.x = "SampleID", by.y = "SampleID")
row.names(meta) <- meta$SampleID
#meta = meta[,-1]
str(meta)


#Run the ANOVA and save it as an object
aov.shannon.sample_time = aov(shannon_entropy ~ Sample_time, data=meta)
aov.evenness.Sample_time = aov(pielou_evenness ~ Sample_time, data=meta)
aov.Observed_features.Sample_time = aov(observed_features ~ Sample_time, data=meta)
aov.Faith.Sample_time = aov (faith_pd ~ Sample_time, data = meta)
#Call for the summary of that ANOVA, which will include P-values
summary(aov.shannon.sample_time)
summary(aov.evenness.Sample_time)
summary(aov.Observed_features.Sample_time)
summary(aov.Faith.Sample_time)

#do Tukey test to make pairwise comparisons
TukeyHSD(aov.Observed_features.Sample_time)
TukeyHSD(aov.evenness.Sample_time)
TukeyHSD(aov.Faith.Sample_time)


levels(meta$Sample_time)
#Re-order the groups 
meta$Sample_time.ord = factor(meta$Sample_time, c("pre-season", "during-season"))
levels(meta$Sample_time.ord)


evenness_summary <- meta %>% # the names of the new data frame and the data frame to be summarised
  group_by(Sample_time) %>%   # the grouping variable
  summarise(mean_evenness = mean(pielou_evenness),  # calculates the mean of each group
            sd_evenness = sd(pielou_evenness), # calculates the standard deviation of each group
            n_evenness = n(),  # calculates the sample size per group
            se_evenness = sd(pielou_evenness)/sqrt(n())) # calculates the standard error of each group

# We can now make a bar plot of means vs sampling time, with standard 
# deviations or standard errors as the error bar. The following code 
# uses the standard deviations.

evenness_se <- ggplot(evenness_summary, aes(Sample_time, mean_evenness, fill = Sample_time)) + 
  geom_col() + 
  geom_errorbar(aes(ymin = mean_evenness - se_evenness, ymax = mean_evenness + se_evenness), width=0.2) + 
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.title = element_blank()) +
  labs(y="Pielou's evenness  ± s.e.", x = "") 

ggsave("output/evenness_se.png", evenness_se, height = 2.5, width = 3)


#observed features
ObFeat_summary <- meta %>% # the names of the new data frame and the data frame to be summarised
  group_by(Sample_time) %>%   # the grouping variable
  summarise(mean_ObFeat = mean(observed_features),  # calculates the mean of each group
            sd_ObF = sd(observed_features), # calculates the standard deviation of each group
            n_ObF = n(),  # calculates the sample size per group
            se_ObF = sd(observed_features)/sqrt(n())) # calculates the standard error of each group

ObFeat_se <- ggplot(ObFeat_summary, aes(Sample_time, mean_ObFeat, fill = Sample_time)) + 
  geom_col() + 
  geom_errorbar(aes(ymin = mean_ObFeat - se_ObF, ymax = mean_ObFeat + se_ObF), width=0.2) + 
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.title = element_blank()) +
  labs(y="Observed Features  ± s.e.", x = "") 

ggsave("output/Observed_Features_se.png", ObFeat_se, height = 2.5, width = 3)


#Shannon
Shannon_summary <- meta %>% # the names of the new data frame and the data frame to be summarised
  group_by(Sample_time) %>%   # the grouping variable
  summarise(mean_Shannon = mean(shannon_entropy),  # calculates the mean of each group
            sd_Shan = sd(shannon_entropy), # calculates the standard deviation of each group
            n_Shan = n(),  # calculates the sample size per group
            se_Shan = sd(shannon_entropy)/sqrt(n())) # calculates the standard error of each group

Shannon_se <- ggplot(Shannon_summary, aes(Sample_time, mean_Shannon, fill = Sample_time)) + 
  geom_col() + 
  geom_errorbar(aes(ymin = mean_Shannon - se_Shan, ymax = mean_Shannon + se_Shan), width=0.2) + 
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.title = element_blank()) +
  labs(y="Shannon Entropy  ± s.e.", x = "") 

ggsave("output/Shannon_Entroyp_se.png", Shannon_se, height = 2.5, width = 3)

#faith PD
Faith_summary <- meta %>% # the names of the new data frame and the data frame to be summarised
  group_by(Sample_time) %>%   # the grouping variable
  summarise(mean_Faith = mean(faith_pd),  # calculates the mean of each group
            sd_Faith = sd(faith_pd), # calculates the standard deviation of each group
            n_Faith = n(),  # calculates the sample size per group
            se_Faith = sd(faith_pd)/sqrt(n())) # calculates the standard error of each group

Faith_se <- ggplot(Faith_summary, aes(Sample_time, mean_Faith, fill = Sample_time)) + 
  geom_col() + 
  geom_errorbar(aes(ymin = mean_Faith - se_Faith, ymax = mean_Faith + se_Faith), width=0.2) + 
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.title = element_blank()) +
  labs(y="Faith Phylogenetic Distance  ± s.e.", x = "") 

ggsave("output/Faith_PD_se.png", Faith_se, height = 2.5, width = 3)






