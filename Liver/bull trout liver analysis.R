#CBest 2023
#Bull Trout Liver Analysis - Transcriptional Changes due to Acclimation Temperature
#Includes qPCR analysis and PCA

#MCMC.qpcr tutorial here (it is very good) https://matzlab.weebly.com/uploads/7/6/2/2/76229469/mcmc.qpcr.tutorial.v1.2.4.pdf
#MCMC.qpcr RDocumentation here https://www.rdocumentation.org/packages/MCMC.qpcr/versions/1.2.4

library(tidyverse)
library(readxl)
library(writexl)
library(FactoMineR)
library(factoextra)
library(MCMC.qpcr)
library(cowplot)
library(psych)
library(paletteer)

####IMPORT DATA####

#Import cleaned data with technical replicates retained
#The MCMC.qpcr tutorial says to retain your technical replicates for analysis
liver_data_clean_wells <- read_excel("liver_data_qc_omit.xlsx")

#Need to uniquely identify the duplicates - can use the array position
#i.e. if duplicate 1 was in B1 and 2 was in D7, this should be the same for every gene
#Take the last 2 characters (subarray position) off of well, leaving the array position and add that to sample
liver_data_clean_wells <- liver_data_clean_wells %>%
  mutate(arraypos = str_sub(Well, end=-3)) %>%
  mutate(sample.unique = paste0(Sample.Name,"_",arraypos))

#Format for correct input to mcmc.qpcr package
liver_data_clean_wells_wide <- liver_data_clean_wells %>% 
  select(Gene.Name, sample.unique, Sample.Name, Biological.Group.Name, Crt) %>% #select a subset of wells
  arrange(Biological.Group.Name) %>% #arrange by group (makes life easier later)
  pivot_wider(names_from = "Gene.Name", values_from = "Crt") %>% #make 1 gene per column of Crt
  rename(sample = Sample.Name) %>% #rename column "sample" required for input
  rename(treatment = Biological.Group.Name) %>% #rename column "treatment" required for input
  mutate(treatment = as.factor(treatment)) %>% #convert treatment to a factor
  mutate(sample.unique=NULL) #remove the "unique identifier" column

#Convert to a data frame
liver_for_analysis <- as.data.frame(liver_data_clean_wells_wide)  

#Import efficiency list as a data frame
liver_efficiency_list <- read.csv("liver_linreg_efficiencies.csv")

#If you get an error message for the upcoming line, be sure to use DATA FRAMES

#Convert Crts to counts using poisson distribution
qs_liver=cq2counts(data=liver_for_analysis,effic=liver_efficiency_list,
                   genecols=c(3:45),condcols=c(1:2),Cq1=37)

#Set the reference treatment to 9 degC (used for plotting model in 'relative' mode)
qs_liver$treatment=relevel(qs_liver$treatment,ref="9deg")


####NAIVE MODEL####

###Run the naive model (no reference genes)
  #Should run this first to check for potential global effects
  #Global effects: all genes being up/downregulated under some conditions
  #This will likely be most obvious in the reference genes

naive_liver=mcmc.qpcr(data=qs_liver, fixed="treatment")
options(max.print=2000)
summary(naive_liver)

#Use summary to check model convergence
  #When the model is complex and/or the data sparse, the chain might not “mix well”,
  #meaning that parameter samples will be not independent (autocorrelated) and/or show a trend
  #Look at eff.samp column - expected (nitt-burnin)/thin = 1000
    #nitt - number of iterations, default 13000
    #burnin - number of initial iterations to discard, default 3000
    #thin - thinning interval, number of iterations beween when parameters are sampled
  #If eff.samp is much lower than 1000, then try nitt=110000, thin=100, burnin=10000
#OK

#You can also check model convergence using these plots
  #Should be random noise around a mean value, no waves/trends etc.
  #Plot function below will generate many plots!
  plot(naive_liver)

####NAIVE/ABSOLUTE####
#Create HPDsummary object of the model, produces 3 things:
  #$summary: summary Table containing calculated abundances, their SD and 95% credible limits
  #$ggPlot: Summary Plot
  #$geneWise: P-Values & fold-changes in a matrix

#This produces the HPDsummary object (this is the absolute abundances)
naive_liver_abs <- HPDsummary(model=naive_liver,data=qs_liver)

#$summary: summary Table containing calculated abundances, their SD and 95% credible limits
naive_liver_abs_table <- naive_liver_abs$summary

#Summary Plot
summaryPlot(naive_liver_abs_table, "treatment", facet = "gene",type="line", 
            x.order=c("6deg", "9deg", "12deg","15deg", "18deg","21deg"))+
  theme(legend.position="none") +
  theme(axis.text.x = element_text(angle=90, vjust=0)) +
  facet_wrap(.~gene, scales="free_y") #can omit this to see everything on same scale

#$geneWise: P-Values & fold-changes in a matrix
#upper triangle - pairwise differences between factor combinations
#lower triangles - corresponding p-values
naive_liver_abs$geneWise
cat(capture.output(print(naive_liver_abs$geneWise),file="naive_liver_abs_geneWise_test.txt"))

#CORRECTION FOR MULTIPLE COMPARISONS
#Benjamini-Hotchberg method
naive_liver_abs_adj = padj.hpdsummary(naive_liver_abs, controls = NULL, method="BH")

#table/plot will be same as without adjustment - only changes p-values
naive_liver_abs_adj_table <- naive_liver_abs_adj$summary

summaryPlot(naive_liver_abs_adj_table, "treatment", facet = "gene", type="line",
            x.order=c("6deg", "9deg", "12deg","15deg", "18deg","21deg"))+
  theme(legend.position="none") +
  theme(axis.text.x = element_text(angle=90, vjust=0)) +
  facet_wrap(.~gene, scales="free_y") #can omit this to see everything on same scale

naive_liver_abs_adj$geneWise
cat(capture.output(print(naive_liver_abs_adj$geneWise),file="naive_liver_abs_adj_geneWise_test.txt"))

#Check genewise p-values - OK for ref genes, do not use ef1a for this or any tissue
#Check plots - no global trends

#There don't appear to by any obvious global effects
#Can use informed model to sharpen credible intervals

####INFORMED MODEL####

#Specify reference genes to sharpen credible intervals
informed_liver_model=mcmc.qpcr(data=qs_liver, fixed="treatment", controls = c("SALRPL13A","SALRPL7","SALRPS9"))
options(max.print=2000)
summary(informed_liver_model)

#Use summary to check model convergence
  #When the model is complex and/or the data sparse, the chain might not “mix well”,
  #meaning that parameter samples will be not independent (autocorrelated) and/or show a trend
  #Look at eff.samp column - expected (nitt-burnin)/thin = 1000
    #nitt - number of iterations, default 13000
    #burnin - number of initial iterations to discard, default 3000
    #thin - thinning interval, number of iterations beween when parameters are sampled
  #If eff.samp is much lower than 1000, then try nitt=110000, thin=100, burnin=10000
#OK

#You can also check model convergence using these plots
  #Should be random noise around a mean value, no waves/trends etc.
  #Plot function below will generate many plots!
  plot(informed_liver_model)

####INFORMED/ABSOLUTE####
informed_liver_abs <- HPDsummary(model=informed_liver_model, data=qs_liver)

informed_liver_abs_table <- informed_liver_abs$summary

summaryPlot(informed_liver_abs_table,xgroup="treatment", facet="gene",type="line", 
            x.order=c("6deg","9deg","12deg","15deg","18deg","21deg"))+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle=90, vjust=0))+
  facet_wrap(.~gene, scales="free_y") #can omit this to see everything on same scale

informed_liver_abs$geneWise
cat(capture.output(print(informed_liver_abs$geneWise),file="informed_liver_abs_geneWise_test.txt"))

#CORRECTION FOR MULTIPLE COMPARISONS
#Benjamini-Hotchberg method
informed_liver_abs_adj = padj.hpdsummary(informed_liver_abs, controls = NULL, method="BH")

#table/plot will be same without adjustment - only changes p-values
informed_liver_abs_adj_table <- informed_liver_abs_adj$summary

summaryPlot(informed_liver_abs_adj_table, "treatment", facet = "gene", type="line",
            x.order=c("6deg", "9deg", "12deg","15deg", "18deg","21deg"))+
  theme(legend.position="none") +
  theme(axis.text.x = element_text(angle=90, vjust=0)) +
  facet_wrap(.~gene, scales="free_y") #can omit this to see everything on same scale

informed_liver_abs_adj$geneWise
cat(capture.output(print(informed_liver_abs_adj$geneWise),file="informed_liver_abs_adj_geneWise_test.txt"))


####MODEL DIAGNOSTICS####

#Ensure modelling approach is valid (meets assumptions of linear model)
#Rerun model with pr/pl true and use diagnostic plots to check assumptions
DIAG_informed_liver_model=mcmc.qpcr(data=qs_liver, fixed="treatment", controls = c("SALRPL13A","SALRPL7","SALRPS9"), pr=T, pl=T)
#pr: should the posterior distribution of random effects be saved? (F by default)
#pl: should the posterior distribution of latent variables be saved? (F by default)
diagnostic.mcmc(
  model=DIAG_informed_liver_model,
  col="grey50",
  cex=0.8
)
#1 - the residuals of the model should not show a trend depending on the predicted value (data are linear)
#2 - the size of the residuals should not change depending on the predicted value (data have homoscedastisity)
#3 - residuals should be approximately normally distributed
#OK


####PLOTS####

#Creating plots from table

#Use the summary table from the HPDSummary object to plot as a ggplot object
#Summary table containing mean (log 2 fold-change or abundance), their SD and 95% credible limits
informed_liver_abs_table

#Convert to tibble and format to facilitate plotting
informed_liver_abs_tibble <- as_tibble(informed_liver_abs_table) %>%
  mutate(temp = str_sub(treatment, end=-4)) %>%
  mutate(temp = as.numeric(temp)) %>%
  arrange(gene, temp) %>%
  mutate(temp = as_factor(temp))

#Swap out gene names (from SALGENE to gene)
gene_function_key <- read_excel("gene_function_key_revised.xlsx")
gene_function_key <- gene_function_key %>% rename(gene = Gene)
gene_function_key
informed_liver_abs_tibble_SYMBOLS <- left_join(informed_liver_abs_tibble,gene_function_key,by="gene")

#LINE GRAPHS:

#Line graph (one gene - don't forget to change TITLE, and y-axis limits if needed)
ggplot(subset(informed_liver_abs_tibble_SYMBOLS, Symbol %in% "cpt1a"), aes(x=temp, y=mean, colour=temp)) +
  scale_colour_manual(values = c("#053061","#4393C3","#F4A582","#D6604D","#B2182B","#67001F"))+
  geom_errorbar(aes(ymin=lower, ymax=upper), width = 0.5)+
  geom_point(position = position_dodge(), stat="identity",size=3) +
  theme_classic()+
  theme(legend.position="none")+
  labs(y=expression(log[2](abundance)), x="Acclimation Temp. (°C)", title="LIVER-CPT1A")+
  ylim(8.5,14)+
  theme(text = element_text(size=10),
        axis.title.y = element_text(margin=margin(0,5,0,0)),
        axis.title.x = element_text(margin=margin(5,0,0,0)),
        axis.text = element_text(colour="black"))
  #copy paste as 200 x 200 pixel for panel

#Line graph (all target genes)
ggplot(subset(informed_liver_abs_tibble_SYMBOLS, !(Symbol %in% c("ef1a","rpl13a","rpl7","rps9"))), aes(x=temp, y=mean, colour=temp)) +
  scale_colour_manual(values = c("#053061","#4393C3","#F4A582","#D6604D","#B2182B","#67001F"))+
  geom_errorbar(aes(ymin=lower, ymax=upper), width = 0.5)+
  geom_point(position = position_dodge(), stat="identity",size=3) +
  theme_classic()+
  theme(legend.position="none")+
  labs(y=expression(log[2](abundance)), x="Acclimation Temp. (°C)", title = "Liver")+
  theme(text = element_text(size=10),
        axis.title.y = element_text(margin=margin(0,5,0,0)),
        axis.title.x = element_text(margin=margin(5,0,0,0)),
        axis.text = element_text(colour="black"))+
  facet_wrap(.~Symbol, scales="free_y")+
  theme(strip.text = element_text(face="italic",size=12),
        strip.background = element_rect(colour="white",fill="grey95"))
  #Save as 1024x768 tiff file or a letter size PDF

####PCA - IMPORT AND PROCESS DATA####

#Use liver_data_clean (just has means of replicates)
liver_data_clean <- read_excel("liver_data_clean.xlsx")

#Ct of target gene minus geomean of reference Cts
#Select just the reference genes, pivot wide, and calculate geometric mean
liver_geomean_refs <- liver_data_clean %>% select(Gene.Name, Sample.Name, Mean.Crt) %>%
  filter(Gene.Name=="SALRPL13A"|Gene.Name=="SALRPL7"|Gene.Name=="SALRPS9") %>%
  pivot_wider(names_from="Gene.Name",values_from="Mean.Crt") %>%
  rowwise() %>%
  mutate(gm = geometric.mean(c(SALRPL13A,SALRPL7,SALRPS9),na.rm=FALSE)) #SPECIFICALLY not na.rm=TRUE because will mess up normalization, just omit sample

#Calculate dCrt
#Add the geometric mean for each sample to liver_data_clean (like vlookup)
liver_dCrt <- left_join(liver_data_clean,liver_geomean_refs,by="Sample.Name") %>%
  rowwise() %>%
  mutate(dCrt = Mean.Crt - gm) %>%
  select(Gene.Name,Biological.Group.Name,Sample.Name,dCrt)

#make wide table, drop ref genes since they are used in the normalization
liver_dCrt_wide_for_pca <- liver_dCrt %>% 
  pivot_wider(names_from = "Gene.Name", values_from = "dCrt") %>%
  mutate(SALRPL13A=NULL, SALRPL7=NULL, SALRPS9=NULL, SALEF1A=NULL) #remove ref genes

#Remove the sample ids column from the file (works better - code knows each row is an individual)
#so col 1 is treatment, and other columns are just genes with counts
liver_dCrt_wide_for_pca
liver_dCrt_wide_for_pca_nosamp <- mutate(liver_dCrt_wide_for_pca, Sample.Name=NULL)
liver_dCrt_wide_for_pca_nosamp

#To help put the temps in the correct order in the plot
#Convert the treatment column to a number, then sort, then change it back to a factor
liver_dCrt_wide_for_pca_nosamp_ORD <- liver_dCrt_wide_for_pca_nosamp %>%
  rename(treatment = Biological.Group.Name) %>%
  mutate(treatment = str_sub(treatment, end=-4)) %>%
  mutate(treatment = as.numeric(treatment)) %>%
  arrange(treatment) %>%
  mutate(treatment = as_factor(treatment))

liver_dCrt_wide_for_pca_nosamp_ORD


####PCA - RUN####

#This is the PCA function from the FactoMineR package
liver.dCrt.pca <- PCA(liver_dCrt_wide_for_pca_nosamp_ORD, scale.unit = TRUE, ncp = 5, 
                      quali.sup = c(1), graph = TRUE, axes = c(1,2))

#Looking at the results of the PCA (contributing to PC1/Dim1, PC2/Dim2)
summary.PCA(liver.dCrt.pca)

#Look at eigenvalues, and the percent of variance accounted for by each component
liver.dCrt.pca$eig

#Look at the contributions of each gene to the principal components
liver.dCrt.pca$var$contrib


####PCA - PLOT####

#First check the plot function from the FactoMiner package (hard to modify appearance)
#Compare with subsequent PCAs to confirm was done correctly
plot(liver.dCrt.pca, graph.type="ggplot",label="none", habillage = 1,
     col.hab = c("#053061","#4393C3","#F4A582","#D6604D","#B2182B","#67001F"),
     col.quali = c("#053061","#4393C3","#F4A582","#D6604D","#B2182B","#67001F"))+
  theme_classic()

#Plots the PCA with confidence ellipses
plotellipses(liver.dCrt.pca, label="none", level=0.95)

#Making nicer (and more editable) PCA using ggplot

#Extract the following to plot the data
individuals.dCrt <- rownames(liver.dCrt.pca$ind$coord)
pc1.dCrt <- liver.dCrt.pca$ind$coord[,1]
pc2.dCrt <- liver.dCrt.pca$ind$coord[,2]
treatment.dCrt <- liver_dCrt_wide_for_pca_nosamp_ORD$treatment
#This last one must be the data set that went into the PCA function.
#If this is done correctly, your PCA should look the same as the one generated
#by the FactoMiner plot functions above!

#Collect these vectors into a tibble
#60 individuals for liver
liver.dCrt_pca_plot_data <- tibble(Ind = individuals.dCrt, Dim1 = pc1.dCrt, Dim2 = pc2.dCrt, Treatment = treatment.dCrt)

#Check the relative contribution to variance in the data (percentage) of the first 2 principal components
#Round the %s so they can be automatically entered in the axis label - just in case!
liver.dCrt.pca$eig
pcent_var_pc1.dCrt <- liver.dCrt.pca$eig[1,2] %>% round(digits=2)
pcent_var_pc2.dCrt <- liver.dCrt.pca$eig[2,2] %>% round(digits=2)

#Add column Temperature with degrees symbol
liver.dCrt_pca_plot_data <- liver.dCrt_pca_plot_data %>%
  mutate(Temperature = paste0(Treatment,"°C")) %>%
  mutate(Temperature = as_factor(Temperature))

# Plot the PCA with ggplot
liver.dCrt_PCA_graph <- ggplot(data = liver.dCrt_pca_plot_data, aes(x = Dim1, y = Dim2, colour = Temperature)) +
  geom_point(alpha=0.8,size = 2) +
  geom_hline(yintercept=0, lty=2) +
  geom_vline(xintercept=0, lty=2) +
  theme_minimal() +
  theme(text = element_text(size = 16))+
  xlab(paste0("PC1 (",pcent_var_pc1.dCrt,"%)"))+
  ylab(paste0("PC2 (",pcent_var_pc2.dCrt,"%)"))+
  stat_ellipse(geom="polygon",aes(fill=Temperature),alpha=0.1,level=0.95)+
  scale_colour_manual(values = c("#053061","#4393C3","#F4A582","#D6604D","#B2182B","#67001F"))+
  scale_fill_manual(values = c("#053061","#4393C3","#F4A582","#D6604D","#B2182B","#67001F"))+
  ggtitle("liver")

liver.dCrt_PCA_graph


####PCA - LOADINGS####

gene_function_key <- read_excel("gene_function_key_revised.xlsx")
gene_function_key

#Loadings - see FactorMiner's FAQ page http://factominer.free.fr/question/FAQ.html
#Loadings (i.e. standard coordinates) are not given by FactoMineR's methods. They return principal coordinates.
#You can calculate them by dividing variables' coordinates on a dimension by this dimension's eigenvalue's square root.
liver.dCrt.loadings <- sweep(liver.dCrt.pca$var$coord,2,sqrt(liver.dCrt.pca$eig[1:ncol(liver.dCrt.pca$var$coord),1]),FUN="/")

liver.dCrt.loadings.tib <- as_tibble(liver.dCrt.loadings, rownames="Gene")
liver.dCrt.loadings.tib

liver.dCrt.loadings.tib <- left_join(liver.dCrt.loadings.tib,gene_function_key,by="Gene")

#Plot PC1 loadings just to get legend for panel
forlegend.dCrt<-ggplot(data=liver.dCrt.loadings.tib, aes(y=reorder(Symbol,Dim.1), x=Dim.1, fill=Function))+
  scale_fill_paletteer_d("ggthemes::Superfishel_Stone")+
  geom_bar(position = position_dodge(), stat="identity")+
  theme_classic()+
  theme(axis.text.y = element_text(face = "italic"))+
  labs(x="Dim.1 Loading",y="Gene")
legend<-get_legend(forlegend.dCrt)

#PC1 Loadings Plot
a.dCrt<-ggplot(data=liver.dCrt.loadings.tib, aes(y=reorder(Symbol,Dim.1), x=Dim.1, fill=Function))+
  scale_fill_paletteer_d("ggthemes::Superfishel_Stone")+
  geom_bar(position = position_dodge(), stat="identity")+
  theme_classic()+
  theme(axis.text.y = element_text(face = "italic"))+
  theme(legend.position="none")+
  labs(x="PC1 Loading",y="Gene")

#PC2 Loadings Plot
b.dCrt<-ggplot(data=liver.dCrt.loadings.tib, aes(y=reorder(Symbol,Dim.2), x=Dim.2, fill=Function))+
  scale_fill_paletteer_d("ggthemes::Superfishel_Stone")+
  geom_bar(position = position_dodge(), stat="identity")+
  theme_classic()+
  theme(axis.text.y = element_text(face = "italic"))+
  theme(legend.position="none")+
  labs(x="PC2 Loading",y="")

#Plot as panel - PC1, PC2, and one legend
plot_grid(a.dCrt,b.dCrt,legend, nrow=1)


####FIN####