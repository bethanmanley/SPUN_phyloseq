library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(decontam); packageVersion("decontam")
library(DT)
library(tidyverse)
library(readxl)
library(dplyr)

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

setwd("/home/bethan/seqdata/alerce/lotus2_SSU_ASVs")

### Load phyloseq object (physeq<-readRDS("physeq_decontam.Rdata"))
physeq <- loadRData("phyloseq.Rdata")

## First add Read_depth as a variable
sample_data(physeq)$Read_depth <- sample_sums(physeq)

#List read depth
get_variable(physeq, "Read_depth")

## Remove OTU counts based on a percentage read count per sample - IF NO CONTROL
## source("/home/bethan/R_scripts/R_functions/filter_OTU_per_sample.R")
## physeq <- filter_OTU_per_sample(physeq,0.01)

### Inspect library sizes

df <- as.data.frame(sample_data(physeq)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(physeq)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()


## Identify contaminants - prevalence - https://bioconductor.org/packages/devel/bioc/vignettes/decontam/inst/doc/decontam_intro.html#identifying-contaminants-in-marker-gene-and-metagenomics-data
sample_data(physeq)$is.neg <- sample_data(physeq)$Sample_or_Control == "Control"
contamdf.prev.1 <- isContaminant(physeq, method="prevalence", neg="is.neg", threshold=0.1)
table(contamdf.prev.1$contaminant)

## Import required files to assess taxonomy of removed OTUs - for the whole process, we need the 'hiera_BLAST.txt' file and the phyloseq object
classification <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
otu_taxonomy <- read_delim("hiera_BLAST.txt")
otu_taxonomy_samp <- as.data.frame(otu_taxonomy)
otu_taxonomy <- otu_taxonomy_samp[,-1]
rownames(otu_taxonomy) <- otu_taxonomy_samp[,1]
colnames(otu_taxonomy) <- classification


# Make phyloseq object of presence-absence in negative controls and true samples
physeq.pa <- transform_sample_counts(physeq, function(abund) 1*(abund>0))
physeq.pa.neg <- prune_samples(sample_data(physeq.pa)$Sample_or_Control == "Control", physeq.pa)
physeq.pa.pos <- prune_samples(sample_data(physeq.pa)$Sample_or_Control == "True sample", physeq.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(physeq.pa.pos), pa.neg=taxa_sums(physeq.pa.neg),
                    contaminant=contamdf.prev.1$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")



## Extract the taxonomic classifications of the identified contaminants

row_indices <- which(contamdf.prev.1$contaminant) #grab the row indices that correspond with identified contaminants to locate taxonomic information in the corresponding OTU file

taxonomy_table <- tibble()

for (i in row_indices){
  loc <-  contamdf.prev.1[i, 0]
  tax_key <- row.names(loc)
  tax_value <- otu_taxonomy[tax_key, ]
  taxonomy_table <- rbind(taxonomy_table, tax_value)
}

names(taxonomy_table) <- classification
datatable(taxonomy_table)


## Prune contaminant taxa from the phyloseq tax_table
physeq_decontam <- prune_taxa(!contamdf.prev.1$contaminant, physeq)

## Save file. Note: this must be opened in R using: phyoseq_decontam <- readRDS("physeq_decontam.Rdata")
saveRDS(physeq_decontam, file="physeq_decontam.Rdata")

## Move the pre-decontaminated phyloseq object into a folder labelled 'pre_decontam'
dir.create("pre_decontam")
file.rename(from="phyloseq.Rdata",to="pre_decontam/phyloseq.Rdata")


## Subset filtered phloseq object to include only the three classes of Mucoromycota that are AMF: "Glomeromycetes", "Archaeosporomycetes" and "Paraglomeromycetes" 

amf_physeq <- physeq%>% subset_taxa(Class =="Glomeromycetes" | Class ==  "Archaeosporomycetes" | Class ==  "Paraglomeromycetes" )

## Save file. Note: this must be opened in R using: amf_physeq <- readRDS("physeq_decontam.Rdata")
saveRDS(amf_physeq, file="amf_physeq.Rdata")

## Open file. amf_physeq<-readRDS("amf_physeq.Rdata")

### OPTIONAL: Test and explore
plot_bar(amf_physeq, fill="Genus")




# Adding all metadata to the phyloseq object

# Read the Excel file
metadata <- read.csv("/home/bethan/seqdata/lesotho/SL-Lesotho23_metadata.csv", row.names = 1)

# Extract the sample data from the phyloseq object
sample_data_existing <- sample_data(amf_physeq)

# Convert sample_data_existing to a data frame for merging
sample_data_existing_df <- as.data.frame(sample_data_existing)

# Ensure the metadata has row names matching the sample names in phyloseq object
sample_names <- rownames(sample_data_existing_df)
metadata <- metadata[match(sample_names, rownames(metadata)), ]

# Combine the existing sample data with the new metadata
combined_metadata <- cbind(sample_data_existing_df, metadata)

# Ensure the combined metadata has rownames matching the sample names
rownames(combined_metadata) <- sample_names

# Convert back to a sample_data object
new_sample_data <- sample_data(combined_metadata)

# Update the phyloseq object with the new sample data
sample_data(amf_physeq) <- new_sample_data

# Verify the update
sample_data(amf_physeq)

saveRDS(amf_physeq, file="amf_physeq.Rdata")


