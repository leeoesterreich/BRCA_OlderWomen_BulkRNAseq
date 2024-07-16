---
title: "AgeStudy_Human_BRCA_BulkRNAseq"
author: "SanghoonLee"
format: qmd
editor: R
---

# Human breast cancer bulk RNA-seq data analysis
Bulk RNA-seq data analysis to run GSVA and PROGENy

![image](https://github.com/user-attachments/assets/73f104e5-8df6-47b3-b495-f38b427d8809)

## Gene Set Variation Analysis (GSVA) and Pathway RespOnsive GENes for activity inference (PROGENy)

-   GSVA, it calculates gene set or pathway scores on a per-sample basis (Hänzelmann et al. 2013a). GSVA transforms a gene by sample gene expression matrix into a gene set by sample pathway enrichment matrix (Hänzelmann et al. 2013. BMC bioinformatics)
-   PROGENy, it is resource that leverages a large compendium of publicly available signaling perturbation experiments to yield a common core of pathway responsive genes for human and mouse. These, coupled with any statistical method, can be used to infer pathway activities from bulk or single-cell transcriptomics (Schubert et al. 2018. Nature Communication)

### Necessary input data files. (See Step2 below)

-   bulk RNA-seq expression count data
-   bulk RNA-seq expression TPMLog2 data
-   bulk RNA-seq sample annotation data
-   Human gene annotation data
-   Gene list up-regulated by E1 hormone

## Preparation 1. Install and load the necessary R packages.

```{r}
library(tidyr)
library(DESeq2)
library(dplyr)
library(qusage)
library(data.table)
library(stringr)
library(ggplot2)

library(GSVA)
library(progeny)
```

## Preparation 2. Set working directory and direct your input data files

```{r}
# Your working directory that this code file
dir <- dirname(rstudioapi::getSourceEditorContext()$path); 
setwd(dir); print(dir)

## You can download this file from NCBI GEO GSEOOOOOOOOO
CountFile <- "HumanERpAge_39404g168s_FeatureCount.txt"

## You can download this file from NCBI GEO GSEOOOOOOOOO
TPMLog2File <- "HumanERpAge_39404g168s_TPMlog2.txt"

## You can download this file from this Github page or NCBI GEO GSEOOOOOOOOO
SmpAnnotFile <- "HumanERpAge_BulkRNAseq_SampleInformation.txt"  

## You can download this file from NCBI FTP, https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/ 
GeneAnnotFile <- "Homo_sapiens.gene_info.txt"  
 
# HallMarkReactome <- "Hallmark"  # or "Reactome"

## You can get this file from "Qureshi et al. 2020., Cell Metabolism"
GeneUpregByE1File <- "SuppleTable1_UpregulatedByE1_NotE2_407g.txt"

InterestingGene <- c("PAK4", "HSD17B7","GREB1","PGR","ESR1","TFF1","CYP19A1","HSD17B2","SAA1","RAB19")
```

# Section I. Read expression count data and run GSVA

## Step1. Read Gene Annotation file

```{r}
GeneAnnot <- fread(GeneAnnotFile, header=TRUE, stringsAsFactors=FALSE); dim(GeneAnnot); GeneAnnot[1:2,]
GeneAnnot_GeneSymbType <- GeneAnnot %>% dplyr::select(Symbol,type_of_gene); head(GeneAnnot_GeneSymbType)
# Filter in only 'protein-coding' and filter out "LOCOOO"
GeneAnnot_ProtCoding <- GeneAnnot_GeneSymbType %>% dplyr::filter(type_of_gene=="protein-coding", !grepl("LOC\\d",Symbol)); dim(GeneAnnot_ProtCoding) # 19492     2
```

## Step2. Read RNA-seq count data and clinical data file.

```{r}
ExpCountData <- fread(CountFile, header=TRUE, stringsAsFactors=FALSE); colnames(ExpCountData)[1] <- "GeneSymb"; dim(ExpCountData); ExpCountData[1:2,1:3] #   39404   169
colnames(ExpCountData) <- gsub("_LEE(.*)","", colnames(ExpCountData)); 
ExpCountData_ProtCoding <- ExpCountData %>% dplyr::filter(ExpCountData$GeneSymb %in% GeneAnnot_ProtCoding$Symbol); dim(ExpCountData_ProtCoding) # 17685 169
ExpCountData_ProtCoding[1,1:3]
#         GeneSymb AL_001 AL_002
#           <char>  <int>  <int>
# 1:        OR4F5      0      0


ExpCountData_ProtCoding_Tp <- ExpCountData_ProtCoding %>% tibble::column_to_rownames("GeneSymb") %>% t %>% data.frame %>% tibble::rownames_to_column("SampleName")
dim(ExpCountData_ProtCoding_Tp); ExpCountData_ProtCoding_Tp[1,1:3]; # 168 18118
#     SampleName OR4F5 LOC112268260
# 1     AL_001     0            0

### Annotation data. 
AnnotData <- fread(SmpAnnotFile, stringsAsFactors=FALSE, header=TRUE); 
colnames(AnnotData)[2] <- "SampleName"; 
AnnotData <- AnnotData %>% dplyr::mutate(AgeRange_Group=paste0(AgeRange,"_",Group), SampleNameGroup=paste0(SampleName,"_",AgeRange,Group)); dim(AnnotData); AnnotData[1,] # 168  8 
# Elderly  Middle   Young 
# 43      85      40 
#       SampleID SampleName AgeRange   Age  TPNumber    Group AgeRange_Group      SampleNameGroup
#       <char>     <char>   <char> <int>    <char>   <char>         <char>               <char>
# 1:      1_1     AL_001    Young    37 TP06_1074    Tumor    Young_Tumor    AL_001_YoungTumor

table(AnnotData$Group)
# Tumor TumorAdj 
# 83       85

## inner_join SampleAnnot and Count_Tp data 
AnnotCount_AllAgeGroup <- dplyr::inner_join(AnnotData[, c("SampleName","SampleNameGroup")], ExpCountData_ProtCoding_Tp) %>% dplyr::select(-SampleName)
dim(AnnotCount_AllAgeGroup); AnnotCount_AllAgeGroup[1,1:5] # 168  17686  
#         SampleNameGroup OR4F5 LOC112268260 OR4F29 LOC105378947
#                   <char>  <int>        <int>  <int>        <int>
# 1:    AL_001_YoungTumor     0            0      0            2

```

## Step3. Make DESeq input data and run DESeq2 to normalize expression count data

```{r}
DESeqMetaData_AgeGroup <- data.frame(AnnotCount_AllAgeGroup$SampleNameGroup); colnames(DESeqMetaData_AgeGroup)[1] <- "AgeGroup";
DESeqMetaData_AgeGroup[1:2,]; dim(DESeqMetaData_AgeGroup) # 168     1
AnnotCount_AllAgeGroup_17685g168s <- AnnotCount_AllAgeGroup %>% tibble::column_to_rownames("SampleNameGroup") %>% t %>% data.frame; 
dim(AnnotCount_AllAgeGroup_17685g168s); AnnotCount_AllAgeGroup_17685g168s[1:2,1:5];  # 17685  168

# Create a `DESeqDataSet` object 
dds_AllAgeGroup <- DESeqDataSetFromMatrix(
  countData = AnnotCount_AllAgeGroup_17685g168s, # Our prepped data frame with counts
  colData = DESeqMetaData_AgeGroup, # Data frame with annotation for our samples.   Just first columnn should have sample IDs. That is enough
  design = ~1 # Here we are not specifying a model
)

dds_norm_AllAgeGroup <- varianceStabilizingTransformation(dds_AllAgeGroup)   # This takes several sectons
saveRDS(dds_norm_AllAgeGroup , file="dds_norm_byIndividualSmp_AllAgeGroup.rds")   ## Save .rds file 

# Retrieve the normalized data from the `DESeqDataSet`   
vst_df_AllAgeGroup <- assay(dds_norm_AllAgeGroup) %>% as.data.frame() %>% tibble::rownames_to_column("GeneSymb") # Make Gene IDs into their own column
vst_df_AllAgeGroup$GeneSymb <- gsub("\\.(.*)","",vst_df_AllAgeGroup$GeneSymb)

## Remove rows of duplicated genes.
vst_df_NoDupGene <- vst_df_AllAgeGroup %>% filter(!duplicated(vst_df_AllAgeGroup$GeneSymb)); dim(vst_df_NoDupGene)  # 17617g 169s

```

## Step4. Import Gene Sets - Gene upregulated by E1, Hallmark Reactome, Biocarta, PID, WIKIPATHWAYS, and GO BP.

```{r}
########## ========= 407 genes upregulated by E1, not E2  ========= ##################
E1UpRegGene <- fread(GeneUpregByE1File, header=TRUE, stringsAsFactors=FALSE); dim(E1UpRegGene); head(E1UpRegGene) # 407 7
colnames(E1UpRegGene) <- gsub("\\ ", "", colnames(E1UpRegGene)); 

## Process gene symbols
E1UpRegGene$Gene <- gsub("-(.*)","", E1UpRegGene$Gene)
E1UpRegGene_NoDup <- E1UpRegGene %>% data.frame %>% dplyr::filter(!duplicated(E1UpRegGene$Gene))
dim(E1UpRegGene_NoDup); head(E1UpRegGene_NoDup) # 366  7
#         Gene E1log2FoldChange pvalue qvalue E2log2FoldChange pvalue qvalue
# 1       RP11            6.773  0.000  0.000            5.112    0.175    0.236
# 2        CTD            6.730  0.000  0.000            5.102    0.047    0.073

E1UpRegGene_NoDup$Gene[which(!E1UpRegGene_NoDup$Gene %in% vst_df_NoDupGene$GeneSymb)]
E1UpRegGeneList <- list(E1UpRegGene_NoDup$Gene); names(E1UpRegGeneList)<-"E1UpRegGene"

table(E1UpRegGene_NoDup$Gene %in% vst_df_NoDupGene$GeneSymb)  # FALSE 88  TRUE 278

########## ========= Hallmark GeneSets  ========= ##################
hallmark_gene_sets <- msigdbr::msigdbr(species = "Homo sapiens", # Can change this to what species you need    
  category = "H" # Only hallmark gene sets
)
dim(hallmark_gene_sets)  # 8209  15    This has just 8209 genes. 
table(hallmark_gene_sets$gene_symbol %in% vst_df_AllAgeGroup$GeneSymb)  # FALSE 505  TRUE 7704

### Convert tibble to list 
hallmarks_list <- split(hallmark_gene_sets$gene_symbol, # The genes we want split into pathways    # I will use gene_symbol
  hallmark_gene_sets$gs_name # The pathways made as the higher levels of the list
)
HallmarkList_Subset <- hallmarks_list[ names(hallmarks_list) %in% c("HALLMARK_ESTROGEN_RESPONSE_EARLY", "HALLMARK_ESTROGEN_RESPONSE_LATE")  ]

## Concatenate Hallmark subset and E1UpRegGeneList
hallmarks_E1UpReg <- c(HallmarkList_Subset,E1UpRegGeneList ); length(hallmarks_E1UpReg) # 3

# ########## ========= BIOCARTA GeneSets  ========= ##################
BIOCARTA_gene_sets <- msigdbr::msigdbr(species = "Homo sapiens", # Can change this to what species you need
                                       category = "C2",subcategory="BIOCARTA")
table(BIOCARTA_gene_sets$gene_symbol %in% vst_df_AllAgeGroup$GeneSymb)  # FALSE 3948  TRUE 12335
vst_df_AllAgeGroup$GeneSymb[which(!vst_df_AllAgeGroup$GeneSymb %in% BIOCARTA_gene_sets$gene_symbol)[1:5]]

### Convert tibble to list
BIOCARTA_list <- split(BIOCARTA_gene_sets$gene_symbol, # The genes we want split into pathways    # I will use gene_symbol
                       BIOCARTA_gene_sets$gs_name # The pathways made as the higher levels of the list
)


# ########## ========= KEGG GeneSets  ========= ##################
KEGG_gene_sets <- msigdbr::msigdbr(species = "Homo sapiens", # Can change this to what species you need
                                   category = "C2",subcategory="KEGG")
table(KEGG_gene_sets$gene_symbol %in% vst_df_AllAgeGroup$GeneSymb)  # FALSE 3948  TRUE 12335
vst_df_AllAgeGroup$GeneSymb[which(!vst_df_AllAgeGroup$GeneSymb %in% KEGG_gene_sets$gene_symbol)[1:5]]
grep("FAM", KEGG_gene_sets$gene_symbol, value=TRUE)
grep("^RP1", KEGG_gene_sets$gene_symbol, value=TRUE)

# ### Convert tibble to list 
KEGG_list <- split(KEGG_gene_sets$gene_symbol, # The genes we want split into pathways    # I will use gene_symbol
                   KEGG_gene_sets$gs_name # The pathways made as the higher levels of the list
)

########## ========= Reactome GeneSets  ========= ##################
Reactome_gene_sets <- msigdbr::msigdbr(species = "Homo sapiens", # Can change this to what species you need
  category = "C2",subcategory="REACTOME")
table(Reactome_gene_sets$gene_symbol %in% vst_df_AllAgeGroup$GeneSymb)  # FALSE 9284  TRUE 90567
vst_df_AllAgeGroup$GeneSymb[which(!vst_df_AllAgeGroup$GeneSymb %in% Reactome_gene_sets$gene_symbol)[1:5]]
grep("FAM", Reactome_gene_sets$gene_symbol, value=TRUE)
grep("^RP1", Reactome_gene_sets$gene_symbol, value=TRUE)

### Convert tibble to list 
Reactome_list <- split(Reactome_gene_sets$gene_symbol, # The genes we want split into pathways    # I will use gene_symbol
  Reactome_gene_sets$gs_name # The pathways made as the higher levels of the list
)

REACTOME_Subset <- Reactome_list[names(Reactome_list) %in% c("REACTOME_ESTROGEN_DEPENDENT_GENE_EXPRESSION")]; length(REACTOME_Subset) # 1 

########## ========= PID GeneSets  ========= ##################
PID_gene_sets <- msigdbr::msigdbr(species = "Homo sapiens", # Can change this to what species you need
                                   category = "C2",subcategory="PID")
table(PID_gene_sets$gene_symbol %in% vst_df_AllAgeGroup$GeneSymb)  # FALSE 400  TRUE 8322
vst_df_AllAgeGroup$GeneSymb[which(!vst_df_AllAgeGroup$GeneSymb %in% PID_gene_sets$gene_symbol)[1:5]]

### Convert tibble to list
PID_list <- split(PID_gene_sets$gene_symbol, # The genes we want split into pathways    # I will use gene_symbol
                  PID_gene_sets$gs_name # The pathways made as the higher levels of the list
)
# 
# # PID_FromBulk <- PID_list[names(PID_list) %in% c("PID_IL2_PI3K_PATHWAY","PID_IL8_CXCR1_PATHWAY", "PID_IL6_7_PATHWAY", "PID_IL3_PATHWAY",
# #                             "PID_IL2_1PATHWAY", "PID_IL4_2PATHWAY", "PID_CXCR4_PATHWAY","PID_IL23_PATHWAY", "PID_IL2_STAT5_PATHWAY", "PID_IL27_PATHWAY") ] # 10
# # length(PID_FromBulk)

######### =========== WIKIPATHWAYS Pathway ========= ##################
WIKIPATHWAYS_gene_sets <- msigdbr::msigdbr(species = "Homo sapiens", # Can change this to what species you need
                                   category = "C2",subcategory="WIKIPATHWAYS")
table(WIKIPATHWAYS_gene_sets$gene_symbol %in% vst_df_AllAgeGroup$GeneSymb)  # FALSE $101  TRUE 29748
vst_df_AllAgeGroup$GeneSymb[which(!vst_df_AllAgeGroup$GeneSymb %in% WIKIPATHWAYS_gene_sets$gene_symbol)[1:5]]

### Convert tibble to list 
WIKIPATHWAYS_list <- split(WIKIPATHWAYS_gene_sets$gene_symbol, # The genes we want split into pathways    # I will use gene_symbol
                           WIKIPATHWAYS_gene_sets$gs_name # The pathways made as the higher levels of the list
)

WIKIPATHWAYS_Subset <- WIKIPATHWAYS_list[names(WIKIPATHWAYS_list) %in% c("WP_ESTROGEN_SIGNALING_PATHWAY")]; length(WIKIPATHWAYS_Subset) # 1 

######### =========== GO Biological Pathway ========= ##################
GOBP_gene_sets <- msigdbr::msigdbr(species = "Homo sapiens", # Can change this to what species you need
                                       category = "C5",subcategory="BP")
table(GOBP_gene_sets$gene_symbol %in% vst_df_AllAgeGroup$GeneSymb)  # FALSE 130479  TRUE 590900
vst_df_AllAgeGroup$GeneSymb[which(!vst_df_AllAgeGroup$GeneSymb %in% GOBP_gene_sets$gene_symbol)[1:5]]

### Convert tibble to list 
GOBP_list <- split(GOBP_gene_sets$gene_symbol, # The genes we want split into pathways    # I will use gene_symbol
                   GOBP_gene_sets$gs_name # The pathways made as the higher levels of the list
)

GOBP_Subset <- GOBP_list[names(GOBP_list) %in% c("GOBP_INTRACELLULAR_ESTROGEN_RECEPTOR_SIGNALING_PATHWAY", "GOBP_CELLULAR_RESPONSE_TO_ESTROGEN_STIMULUS")]; length(GOBP_Subset) # 2
```

## Step5. Make Estrogen-related pathway list and run GSVA

```{r}
## Make Estrogen-related pathway list
EstrogenPathway_List <- c(hallmarks_E1UpReg,REACTOME_Subset,WIKIPATHWAYS_Subset, GOBP_Subset); length(EstrogenPathway_List); names(EstrogenPathway_List)  # 7

## Convert data.frame to matrix
vst_df_RowGeneSymb <- vst_df_NoDupGene %>% tibble::column_to_rownames("GeneSymb") %>% as.matrix # GeneSymb is rownames. 
saveRDS(vst_df_RowGeneSymb, file="vst_df_RowGeneSymb_IndividualSample.rds")

####### Run GSVA with estrogen pathway 
GSVAResult <- GSVA::gsva(gsvaParam(vst_df_RowGeneSymb, EstrogenPathway_List,
                              # Appropriate for our vst transformed data
                              kcdf = "Gaussian",  ## when 'filltered_mapped_matrix is continuous value. vst_df has continous values.  #kcdf = "Poisson" on integer counts
                              maxDiff = TRUE))
dim(GSVAResult); GSVAResult[1:5,1:5]  ## 7 168
saveRDS(GSVAResult, "GSVAOut_BulkRNAseq_EstrogenPathway_AllAgeGroup.rds")
```

# Section II. Read expression TPMLog2 data and run PROGENy

## Step6. Read gene expression TPMLog2 data

```{r}
TPMLog2Data <- fread(TPMLog2File, stringsAsFactors=FALSE, header=TRUE); TPMLog2Data[1:3,1:3]
colnames(TPMLog2Data)[1] <- "GeneSymb"; dim(TPMLog2Data); TPMLog2Data[1:3,1:3] # 40826  67
### filter by only protein-coding genes 
TPMLog2Data_ProtCoding <- TPMLog2Data %>% dplyr::filter(TPMLog2Data$GeneSymb %in% GeneAnnot_ProtCoding$Symbol); dim(TPMLog2Data_ProtCoding) # 18117 169
colnames(TPMLog2Data_ProtCoding) <- gsub("_LEE(.*)","", colnames(TPMLog2Data_ProtCoding)); TPMLog2Data_ProtCoding[1,1:3]
# GeneSymb AL_001 AL_002
# 1:        OR4F5      0      0

```

## Step7. Make Progeny input data from TPMLog2 data - All YoungElderlyTumor

```{r}
TPMLog2Data_ProtCoding_Tp <- TPMLog2Data_ProtCoding %>% tibble::column_to_rownames("GeneSymb") %>% t %>% data.frame %>% tibble::rownames_to_column("SampleName")
dim(TPMLog2Data_ProtCoding_Tp); TPMLog2Data_ProtCoding_Tp[1:3,1:3]; # 168 17686
#   SampleName OR4F5 LOC112268260
# 1   AL_001     0   0.00000000
# 2   AL_002     0   0.00000000

## inner_join SampleAnnot and TPMLOg2 data 
AnnotTPMLog2_ElderlyYoungTumor <- dplyr::inner_join(AnnotData[, c("SampleName","SampleNameGroup")], TPMLog2Data_ProtCoding_Tp) %>% dplyr::select(-SampleName)
dim(AnnotTPMLog2_ElderlyYoungTumor); AnnotTPMLog2_ElderlyYoungTumor[1:3,1:5] # 168  17686 
#       SampleNameGroup OR4F5 LOC112268260 OR4F29 LOC105378947
# 1: AL_001_YoungTumor     0   0.00000000      0   0.05604033
# 2: AL_003_YoungTumor     0   0.04086972      0   0.00000000

AnnotTPMLog2_ElderlyYoungTumor_Tp <- AnnotTPMLog2_ElderlyYoungTumor %>% tibble::column_to_rownames("SampleNameGroup") %>% t %>% data.frame; 
AnnotTPMLog2_ElderlyYoungTumor_Tp[1,1:3]; dim(AnnotTPMLog2_ElderlyYoungTumor_Tp) # 17686 168
#               AL_001_YoungTumor AL_003_YoungTumor AL_005_YoungTumor
# OR4F5                        0        0.00000000        0.00000000

```

## Step8. Run PROGENy YoungElderlyTumor

```{r}
# Run PROGENy and calculate pathway activities
PathwayActivity <- progeny(as.matrix(AnnotTPMLog2_ElderlyYoungTumor_Tp), scale=FALSE, organism="Human", top=100, perm=1000); PathwayActivity[1,]; dim(PathwayActivity) # 168 14
saveRDS(PathwayActivity , file="PROGENy_PathwayActivity_WholeYoungMidElderly168s_Top100Perm1000.rds")
```

# Section III. inner_join, GSVA result, PROGENy result, and gene expression TPMLog2 data.

## Step9. Process GSVA result

```{r}
# GSVAResult <- readRDS("GSVAOut_BulkRNAseq_EstrogenPathway_AllAgeGroup.rds"); dim(GSVAResult); GSVAResult[1,1:3] # 7 168 
#                                 AL_001_YoungTumor AL_002_YoungTumorAdj AL_003_YoungTumor
# HALLMARK_ESTROGEN_RESPONSE_EARLY       -0.05506254          -0.09394596        -0.4510458

GSVAResult_Tp <- GSVAResult %>% t %>% data.frame %>% tibble::rownames_to_column("SampleNameGroup");dim(GSVAResult_Tp); GSVAResult_Tp[1:2,1:3] # 168 8
```

## Step10. inner_join TPMlog2 and AnnotData, and then with Progeny and GSVA result.

```{r}
## inner_join SampleAnnot Age and TPMLog2 data, and select only interesting 9 genes
AnnotTPMLog2_Age <- dplyr::inner_join(AnnotData[, c("SampleName","Age", "SampleNameGroup")], TPMLog2Data_ProtCoding_Tp) 
AnnotTPMLog2_Age_InterestGene <- AnnotTPMLog2_Age %>% dplyr::select_if(colnames(AnnotTPMLog2_Age) %in% c("Age", "SampleNameGroup", InterestingGene)); 
dim(AnnotTPMLog2_Age_InterestGene); AnnotTPMLog2_Age_InterestGene[1,]  # 168  11
#       Age   SampleNameGroup  HSD17B7    GREB1     ESR1  SAA1      PGR CYP19A1 HSD17B2     PAK4     TFF1
#       <int>            <char>    <num>    <num>    <num> <num>    <num>   <num>   <num>    <num>    <num>
#   1:    37 AL_001_YoungTumor 1.715383 2.944336 1.861138     0 3.517439       0 1.46526 2.463698 3.892977

## ============= inner_join with GeneExp, Pathway activity, GSVA result data. - only Tumor samples With MidAge ================
PathwayActivity_YesMidAgeTumor_Estrogen <- PathwayActivity %>% data.frame %>% dplyr::filter(grepl("YoungTumor$|MiddleTumor$|ElderlyTumor$",rownames(PathwayActivity))) %>% 
                            dplyr::select("Estrogen") %>% tibble::rownames_to_column("SampleNameGroup") #  %>% dplyr::mutate(SampleName=gsub("_YoungTumor|_ElderlyTumor","",SampleName));
dim(PathwayActivity_YesMidAgeTumor_Estrogen);PathwayActivity_YesMidAgeTumor_Estrogen[1:2,] # 83 2

TPMLog2AgeGeneExp_ProgenyGSVA_WithMidAge <-  dplyr::inner_join(AnnotTPMLog2_Age_InterestGene, PathwayActivity_YesMidAgeTumor_Estrogen) %>% dplyr::inner_join(GSVAResult_Tp) %>% data.frame 
dim(TPMLog2AgeGeneExp_ProgenyGSVA_WithMidAge); TPMLog2AgeGeneExp_ProgenyGSVA_WithMidAge[1:2,] # 83 18
#  Age     SampleNameGroup  HSD17B7     GREB1      ESR1       SAA1       PGR    CYP19A1   HSD17B2      PAK4     TFF1 Estrogen HALLMARK_ESTROGEN_RESPONSE_EARLY HALLMARK_ESTROGEN_RESPONSE_LATE E1UpRegGene REACTOME_ESTROGEN_DEPENDENT_GENE_EXPRESSION
# 1 37       AL_001_YoungTumor 1.715383 2.9443357 1.8611384 0.00000000 3.5174386 0.00000000 1.4652595 2.4636982 3.892977    0.862                      -0.05506254                      0.01541189 -0.06346937                                   0.2015138

## ============ inner_join with GeneExp, Pathway activity, GSVA result data.  - Tumor samples Without MidAge ================
PathwayActivity_NoMidAgeTumor_Estrogen <- PathwayActivity %>% data.frame %>% dplyr::filter(grepl("YoungTumor$|ElderlyTumor$",rownames(PathwayActivity))) %>% 
  dplyr::select("Estrogen") %>% tibble::rownames_to_column("SampleNameGroup") # %>% dplyr::mutate(SampleName=gsub("_YoungTumor|_ElderlyTumor","",SampleName));
dim(PathwayActivity_NoMidAgeTumor_Estrogen) # 41 18
TPMLog2AgeGeneExp_ProgenyGSVA_NoMidAge <-  dplyr::inner_join(AnnotTPMLog2_Age_InterestGene, PathwayActivity_NoMidAgeTumor_Estrogen) %>% dplyr::inner_join(GSVAResult_Tp) %>% data.frame 
dim(TPMLog2AgeGeneExp_ProgenyGSVA_NoMidAge); TPMLog2AgeGeneExp_ProgenyGSVA_NoMidAge[1:2,] # 83 18
```

## Step11. Calculate correlation between PROGENy activity and gene expression - including or excluding Middle Age

```{r}
MyGene <- c("Age", InterestingGene ); GeneLength <- length(MyGene); print(GeneLength); print(MyGene) # 11 # [1] "Age"     "PAK4"    "HSD17B7" "GREB1"   "PGR"     "ESR1"    "TFF1"    "CYP19A1" "HSD17B2" "SAA1"  "RAB19"  
MyPathwayName <- colnames(TPMLog2AgeGeneExp_ProgenyGSVA_NoMidAge)[(GeneLength+2):ncol(TPMLog2AgeGeneExp_ProgenyGSVA_NoMidAge)]; print(MyPathwayName)
# [1] "Estrogen"                                               "HALLMARK_ESTROGEN_RESPONSE_EARLY"                       "HALLMARK_ESTROGEN_RESPONSE_LATE"                        "E1UpRegGene"                                           
# [5] "REACTOME_ESTROGEN_DEPENDENT_GENE_EXPRESSION"            "WP_ESTROGEN_SIGNALING_PATHWAY"                          "GOBP_CELLULAR_RESPONSE_TO_ESTROGEN_STIMULUS"            "GOBP_INTRACELLULAR_ESTROGEN_RECEPTOR_SIGNALING_PATHWAY"    

Func_CalculateCorrel <- function(BRCAAge_GeneCellType) {
    All_CorrTestSummary_Gene <- data.frame()    ; GeneCount=0;   #  data.frame(matrix(ncol=2, nrow=length(MyCellType)))
    #rownames(All_CorrTestSummary_Gene) <- MyCellType
    for(EachGene in MyGene) {
        # EachGene <- MyGene[1]; print(paste0("EachGene: ", EachGene))
        GeneCount=GeneCount+1;
        print(paste0("current GeneCount: ", GeneCount))
        
        CorrX <- BRCAAge_GeneCellType[ , EachGene]; length(CorrX)  # 24
        
        All_CorrTestSummary <- data.frame(); PathwayCount<-0
        for(EachPathway in MyPathwayName) {
            # EachPathway <- MyPathwayName[1]
            PathwayCount=PathwayCount+1;
            print(paste0("current PathwayCount: ", PathwayCount))
            CorrY <- BRCAAge_GeneCellType[, EachPathway];
            CorrTest <- cor.test(CorrX, CorrY, method=c( "spearman"), exact=FALSE)
            # plot(CorrX, CorrY)
            
            if(is.na(CorrTest$p.value)) {
              CorrTest$estimate <- "NotAvail"
              CorrTest$p.value<-"NotAvail"
            }
            
            CorrTestSummary_Cyto <- data.frame(t(c(CorrTest$estimate, CorrTest$p.value, EachGene) ))
            colnames(CorrTestSummary_Cyto) <- c(paste0("Spearman_Rho"),paste0("Spearman_pval"), "GeneSymb")
            rownames(CorrTestSummary_Cyto) <- EachPathway
            
            All_CorrTestSummary <- rbind(All_CorrTestSummary, CorrTestSummary_Cyto)
        }
        
        All_CorrTestSummary <- All_CorrTestSummary %>% tibble::rownames_to_column("PathwayName") %>% dplyr::filter(Spearman_Rho!="NotAvail"); head(All_CorrTestSummary)
        
        
        #if(GeneCount==1) {
        #All_CorrTestSummary_Gene <- All_CorrTestSummary
        #} else {
        All_CorrTestSummary_Gene <- rbind(All_CorrTestSummary_Gene, All_CorrTestSummary)
        #}
    }
    
    #colnames(BRCAAge_GeneCellType)[colnames(BRCAAge_GeneCellType) == "CurrentGene"] <- EachGene
    return(All_CorrTestSummary_Gene)     
}

## This is Tumor only with or without MidAge group
All_CorrTestSummary_ProgenyGSVA_WithMidAge <- Func_CalculateCorrel(BRCAAge_GeneCellType=TPMLog2AgeGeneExp_ProgenyGSVA_WithMidAge); dim(All_CorrTestSummary_ProgenyGSVA_WithMidAge) # 88 4
All_CorrTestSummary_ProgenyGSVA_NoMidAge <- Func_CalculateCorrel(BRCAAge_GeneCellType=TPMLog2AgeGeneExp_ProgenyGSVA_NoMidAge); dim(All_CorrTestSummary_ProgenyGSVA_NoMidAge) # 88 4

```

## Step12. Make bubble plot input - tidyr format, and draw the plot.

```{r}
### ======= $$$$$$$$$$  ====== With MidAge group, in only Tumor ======== $$$$$$$$$$  ========= ####<<== This is better
All_CorrTestSummary_ProgenyGSVA_WithMidAge$Spearman_Rho <- as.numeric(All_CorrTestSummary_ProgenyGSVA_WithMidAge$Spearman_Rho)
All_CorrTestSummary_ProgenyGSVA_WithMidAge$Spearman_pval <- as.numeric(All_CorrTestSummary_ProgenyGSVA_WithMidAge$Spearman_pval)

# ## If you want to remove Age and E1UpRegGene signature - Tumor
# All_CorrTestSummary_ProgenyGSVA_WithMidAge_NoAgeNoE1UpReg <- All_CorrTestSummary_ProgenyGSVA_WithMidAge %>% dplyr::filter(PathwayName!="E1UpRegGene",GeneSymb!="Age"); 
# All_CorrTestSummary_ProgenyGSVA_WithMidAge_NoAgeNoE1UpReg$GeneSymb <- factor(All_CorrTestSummary_ProgenyGSVA_WithMidAge_NoAgeNoE1UpReg$GeneSymb, 
#                                                                       levels=c("HSD17B2","HSD17B7","CYP19A1","RAB19","SAA1","PAK4","TFF1","PGR","GREB1","ESR1"))
# dim(All_CorrTestSummary_ProgenyGSVA_WithMidAge_NoAgeNoE1UpReg); All_CorrTestSummary_ProgenyGSVA_WithMidAge_NoAgeNoE1UpReg[1:2,] # 63  4

### ==========  Without MidAge group
All_CorrTestSummary_ProgenyGSVA_NoMidAge$Spearman_Rho <- as.numeric(All_CorrTestSummary_ProgenyGSVA_NoMidAge$Spearman_Rho)
All_CorrTestSummary_ProgenyGSVA_NoMidAge$Spearman_pval <- as.numeric(All_CorrTestSummary_ProgenyGSVA_NoMidAge$Spearman_pval)

## Remove "RAB19" gene and  Adjust the order of genes in Y-axis 
All_CorrTestSummary_ProgenyGSVA_NoMidAge <- All_CorrTestSummary_ProgenyGSVA_NoMidAge %>% dplyr::filter(!grepl("RAB19", GeneSymb))
All_CorrTestSummary_ProgenyGSVA_NoMidAge$GeneSymb <- factor(All_CorrTestSummary_ProgenyGSVA_NoMidAge$GeneSymb, 
                                                            levels=c("TFF1",  "SAA1","PGR","PAK4","HSD17B7","HSD17B2","GREB1","ESR1", "CYP19A1", "Age"))



my_palette <- colorRampPalette(c( "blue", "dodgerblue","yellow",  "orange", "red"), alpha=TRUE)(n=399) # "#FFFF99"

MyDotlot_Front <- ggplot(All_CorrTestSummary_ProgenyGSVA_NoMidAge, aes(x=GeneSymb,y=PathwayName)) +
  geom_point(aes(size=-log10(Spearman_pval),color=Spearman_Rho)) +
  scale_color_gradientn('Spearman Rho', colors=my_palette) +
  theme_bw() +   # coord_flip() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text=element_text(size=14, colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(size=12, colour = "black"),
        axis.title=element_blank(),
        panel.border = element_rect(linewidth = 0.7, linetype = "solid", colour = "black")) + coord_flip()
MyDotlot_Front

## save your bubble plot into a pdf file
ggsave(MyDotlot_Front, filename=paste0("OutBubbleplot_Correlation_NoMidAgeGeneExp_ProgenyGSVAActivity.pdf"), width=12, height=7, limitsize=F)  

```
You need to polish this figure with your own skills or methods.

![image](https://github.com/user-attachments/assets/a2207bac-9511-418e-b3df-be9f53073a31)



## The end. 
