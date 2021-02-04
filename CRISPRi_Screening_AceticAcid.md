CRISPRi PHENOMICS DATA ANALYSIS
================
Vaskar Mukherjee
2/3/2021

  - [IMPORT SCAN-O-MATIC RAW DATA](#import-scan-o-matic-raw-data)
      - [PURPOSE](#purpose)
      - [IMPORTING THE METADATA FILE](#importing-the-metadata-file)
      - [GENERATE BASAL **ABSOLUTE**
        DATASET](#generate-basal-absolute-dataset)
      - [GENERATE ACETIC ACID **ABSOLUTE**
        DATASET](#generate-acetic-acid-absolute-dataset)
      - [GENERATE BASAL **NORMALIZED**
        DATASET](#generate-basal-normalized-dataset)
      - [GENERATE ACETIC ACID **NORMALIZED**
        DATASET](#generate-acetic-acid-normalized-dataset)
      - [COMBINE THE DATASETS TO OBTAIN FINAL
        DATAFRAME](#combine-the-datasets-to-obtain-final-dataframe)
      - [IMPORT RESULTS FROM ROUND1](#import-results-from-round1)
      - [COMBINE THE DATASETS of ROUND 1 AND
        2](#combine-the-datasets-of-round-1-and-2)
  - [SCAN-O-MATIC PHENOMICS ANALYSIS](#scan-o-matic-phenomics-analysis)
      - [ESTIMATE THE LOG PHENOTYPIC INDEX (LPI)
        VALUES](#estimate-the-log-phenotypic-index-lpi-values)
      - [PERFORM PLATE-WISE BATCH
        CORRECTION](#perform-plate-wise-batch-correction)
      - [ESTIMATE THE BATCH CORRECTED LOG PHENOTYPIC INDEX (LPI)
        VALUES](#estimate-the-batch-corrected-log-phenotypic-index-lpi-values)
      - [SETTING THE NAMES OF THE NEW
        COLUMNS](#setting-the-names-of-the-new-columns)
      - [EXTRACT ONLY THE BATCH CORRECTED
        COLUMNS](#extract-only-the-batch-corrected-columns)
      - [CONSTRUCT A NEW DATA
        STRUCTURE](#construct-a-new-data-structure)
          - [REMOVE ROWS](#remove-rows)

# IMPORT SCAN-O-MATIC RAW DATA

The phenotypic data generated in scan-o-matic screening in .csv format.
We extract both the absolute and the normalized phenotypes.

The CRISPRi strains in the library were arrayed in 24 plates in 384
format. Each CRISPRi plate was subjected to two different condition
(Basal and 150 mM of Acetic acid). Therefore, for each plate four
different files are generated. All files generated in a single
independent experimental round are stored in a single folder.

  - **SOM\_SCR\_R001** : Raw data for round1

  - **SOM\_SCR\_R002** : Raw data for round2

**ABSOLUTE DATA**

The **Absolute** dataset gives the extracted phenotypes without any
spatial normalization

**NORMALIZED DATA**

The **Normalized** dataset is generated after removal of any spatial
bias. This is in log2 scale and referred as Log Strain Coefficient (LSC)
values

**FILE NAMING**

Each file is named with the plate identifier in such a way so that it
can be easily called programmatically

Eg. Plate 1 absolute data in basal (Ctrl) condition have the following
string  
*Ctrl1.phenotypes.Absolute*  
AND  
Plate 1 Normalized data in acetic acid (aa) stress have the string  
*aa1.phenotypes.Normalized*

## PURPOSE

At the end of this data import session, a single data.frame will be
generated with the data of 24 plates. The whole dataset will be labeled
with the strains attributes using the metadata key file (provided in the
COMPILED\_DATA folder). The data import below is shown for only Round2
dataset. Round1 can be generated modifying the folder location

**METADATA KEY FILE** : library\_keyfile1536.csv

## IMPORTING THE METADATA FILE

``` r
Metadata_CRISPRi <- read.csv("COMPILED_DATA/library_keyfile1536.csv", na.strings = "#N/A", stringsAsFactors = FALSE)
str(Metadata_CRISPRi)
```

    ## 'data.frame':    36864 obs. of  11 variables:
    ##  $ SL_No             : num  1 2 3 4 5 6 7 8 9 10 ...
    ##  $ gRNA_name         : chr  "FBP26-TRg-1" "FBP26-TRg-1" "HMI1-NRg-1" "HMI1-NRg-1" ...
    ##  $ Seq               : chr  "GCTTATCATACATTTACATC" "GCTTATCATACATTTACATC" "AAAAATTCTGACACATCACA" "AAAAATTCTGACACATCACA" ...
    ##  $ SOURCEPLATEID     : chr  "R2877.H.001" "R2877.H.001" "R2877.H.001" "R2877.H.001" ...
    ##  $ SOURCEDENSITY     : chr  "384A" "384A" "384A" "384A" ...
    ##  $ SOURCECOLONYCOLUMN: int  1 1 2 2 3 3 4 4 5 5 ...
    ##  $ SOURCECOLONYROW   : chr  "A" "A" "A" "A" ...
    ##  $ border            : logi  TRUE TRUE TRUE TRUE TRUE TRUE ...
    ##  $ GENE              : chr  "FBP26" "FBP26" "HMI1" "HMI1" ...
    ##  $ Control.gRNA      : int  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ Location_1536     : chr  "A1" "A2" "A3" "A4" ...

## GENERATE BASAL **ABSOLUTE** DATASET

``` r
m <- vector(mode = "character", length = 0)
file.names<-vector(mode = "character", length = 0)
temp_df<-data.frame()
data_Ctrl_Abs <- data.frame()
for(i in 1:24){
  m <- paste0("Ctrl", i, ".phenotypes.Absolute") 
  file.names[i] <- dir("RAW_DATA/SOM_SCR_R002/", pattern = m, full.names = TRUE)
  temp_df <- read.csv(file.names[i], na.strings = "NoGrowth")
  data_Ctrl_Abs <- rbind(data_Ctrl_Abs, temp_df)
}
str(data_Ctrl_Abs)
```

    ## 'data.frame':    36864 obs. of  18 variables:
    ##  $ Plate                                   : int  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ Row                                     : int  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ Column                                  : int  0 1 2 3 4 5 6 7 8 9 ...
    ##  $ Phenotypes.InitialValue                 : num  80116 83495 92467 97104 99622 ...
    ##  $ Phenotypes.ExperimentBaseLine           : num  81419 84504 94039 98641 101354 ...
    ##  $ Phenotypes.ExperimentEndAverage         : num  8615209 6710686 6037319 5554230 5234758 ...
    ##  $ Phenotypes.ColonySize48h                : num  6487316 5329111 4946774 4684166 4440201 ...
    ##  $ Phenotypes.ChapmanRichardsParam2        : num  14.5 12.8 42.7 17.6 18.3 ...
    ##  $ Phenotypes.ChapmanRichardsParam3        : num  -2.74 -2.65 -2.58 -2.52 -2.47 ...
    ##  $ Phenotypes.ChapmanRichardsParamXtra     : num  16 16 16.2 16.3 16.3 ...
    ##  $ Phenotypes.ChapmanRichardsParam1        : num  1.95 1.89 1.84 1.8 1.77 ...
    ##  $ Phenotypes.ChapmanRichardsParam4        : num  -31.33 -31.47 -3.93 -2.78 -2.48 ...
    ##  $ Phenotypes.GenerationTimeStErrOfEstimate: num  0.012879 0.002174 0.000633 0.001667 0.000759 ...
    ##  $ Phenotypes.ExperimentGrowthYield        : num  8533790 6626182 5943280 5455589 5133404 ...
    ##  $ Phenotypes.GenerationTime               : num  2.53 2.48 2.51 2.56 2.53 ...
    ##  $ Phenotypes.ExperimentPopulationDoublings: num  6.73 6.31 6 5.82 5.69 ...
    ##  $ Phenotypes.ChapmanRichardsFit           : num  0.999 0.999 0.999 0.999 0.999 ...
    ##  $ Phenotypes.GenerationTimeWhen           : num  6.14 4.1 4.1 3.76 3.76 ...

Several phenotypes are extracted. However, the most useful for this
study will be,

  - Column No: 14 i.e. **Phenotypes.ExperimentGrowthYield**
  - Column No: 15 i.e. **Phenotypes.GenerationTime**

Extract only this two column in the final data.frame  
Rename the column names to prevent any ambiguity

``` r
data_Ctrl_Abs_Trim <- data_Ctrl_Abs[, 14:15]
colnames(data_Ctrl_Abs_Trim) <- c("CTRL_Y_ABS", "CTRL_GT_ABS")
str(data_Ctrl_Abs_Trim)
```

    ## 'data.frame':    36864 obs. of  2 variables:
    ##  $ CTRL_Y_ABS : num  8533790 6626182 5943280 5455589 5133404 ...
    ##  $ CTRL_GT_ABS: num  2.53 2.48 2.51 2.56 2.53 ...

## GENERATE ACETIC ACID **ABSOLUTE** DATASET

Following the same strategy as above

``` r
m <- vector(mode = "character", length = 0)
file.names<-vector(mode = "character", length = 0)
temp_df<-data.frame()
data_AA_Abs <- data.frame()
for(i in 1:24){
  m <- paste0("aa", i, ".phenotypes.Absolute") 
  file.names[i] <- dir("RAW_DATA/SOM_SCR_R002/", pattern = m, full.names = TRUE)
  temp_df <- read.csv(file.names[i], na.strings = "NoGrowth")
  data_AA_Abs <- rbind(data_AA_Abs, temp_df)
}
data_AA_Abs_Trim <- data_AA_Abs[, 14:15]
colnames(data_AA_Abs_Trim) <- c("AA_Y_ABS", "AA_GT_ABS")
```

## GENERATE BASAL **NORMALIZED** DATASET

``` r
m <- vector(mode = "character", length = 0)
file.names<-vector(mode = "character", length = 0)
temp_df<-data.frame()
data_Ctrl_Norm <- data.frame()
for(i in 1:24){
  m <- paste0("Ctrl", i, ".phenotypes.Normalized") 
  file.names[i] <- dir("RAW_DATA/SOM_SCR_R002/", pattern = m, full.names = TRUE)
  temp_df <- read.csv(file.names[i], na.strings = "NoGrowth")
  data_Ctrl_Norm <- rbind(data_Ctrl_Norm, temp_df)
}
str(data_Ctrl_Norm)
```

    ## 'data.frame':    36864 obs. of  8 variables:
    ##  $ Plate                                   : int  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ Row                                     : int  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ Column                                  : int  0 1 2 3 4 5 6 7 8 9 ...
    ##  $ Phenotypes.ExperimentGrowthYield        : num  0.755 0.39 0.233 0.484 0.396 ...
    ##  $ Phenotypes.GenerationTime               : num  0.0505 0.0204 0.0389 -0.0173 -0.0327 ...
    ##  $ Phenotypes.ExperimentPopulationDoublings: num  0.1907 0.0991 0.0272 0.113 0.0818 ...
    ##  $ Phenotypes.ExperimentBaseLine           : num  -0.0884 -0.0347 0.1195 0.0367 0.0758 ...
    ##  $ Phenotypes.ColonySize48h                : num  0.584 0.301 0.193 0.397 0.32 ...

The most useful for this study will be,

  - Column No: 4 i.e. Phenotypes.ExperimentGrowthYield
  - Column No: 5 i.e. Phenotypes.GenerationTime

Extract only this two column

``` r
data_Ctrl_Norm_Trim <- data_Ctrl_Norm[, 4:5]
colnames(data_Ctrl_Norm_Trim) <- c("CTRL_Y_NORM", "CTRL_GT_NORM")
```

## GENERATE ACETIC ACID **NORMALIZED** DATASET

Same as above

``` r
m <- vector(mode = "character", length = 0)
file.names<-vector(mode = "character", length = 0)
temp_df<-data.frame()
data_AA_Norm <- data.frame()
for(i in 1:24){
  m <- paste0("aa", i, ".phenotypes.Normalized") 
  file.names[i] <- dir("RAW_DATA/SOM_SCR_R002/", pattern = m, full.names = TRUE)
  temp_df <- read.csv(file.names[i], na.strings = "NoGrowth")
  data_AA_Norm <- rbind(data_AA_Norm, temp_df)
}
data_AA_Norm_Trim <- data_AA_Norm[, 4:5]
colnames(data_AA_Norm_Trim) <- c("AA_Y_NORM", "AA_GT_NORM")
```

## COMBINE THE DATASETS TO OBTAIN FINAL DATAFRAME

Trimmed datasets are combined to obtain the final data.frame. The
combined data frame is labeled as data from ROUND2

``` r
R <- rep("2nd_round", 36864)
Round_ID <- data.frame(R, stringsAsFactors = FALSE)
whole_data_R2 <- cbind(Metadata_CRISPRi, 
                       Round_ID, 
                       data_Ctrl_Abs_Trim, 
                       data_AA_Abs_Trim, 
                       data_Ctrl_Norm_Trim, 
                       data_AA_Norm_Trim)
colnames(whole_data_R2)[12] <- "Round_ID"
str(whole_data_R2)
```

    ## 'data.frame':    36864 obs. of  20 variables:
    ##  $ SL_No             : num  1 2 3 4 5 6 7 8 9 10 ...
    ##  $ gRNA_name         : chr  "FBP26-TRg-1" "FBP26-TRg-1" "HMI1-NRg-1" "HMI1-NRg-1" ...
    ##  $ Seq               : chr  "GCTTATCATACATTTACATC" "GCTTATCATACATTTACATC" "AAAAATTCTGACACATCACA" "AAAAATTCTGACACATCACA" ...
    ##  $ SOURCEPLATEID     : chr  "R2877.H.001" "R2877.H.001" "R2877.H.001" "R2877.H.001" ...
    ##  $ SOURCEDENSITY     : chr  "384A" "384A" "384A" "384A" ...
    ##  $ SOURCECOLONYCOLUMN: int  1 1 2 2 3 3 4 4 5 5 ...
    ##  $ SOURCECOLONYROW   : chr  "A" "A" "A" "A" ...
    ##  $ border            : logi  TRUE TRUE TRUE TRUE TRUE TRUE ...
    ##  $ GENE              : chr  "FBP26" "FBP26" "HMI1" "HMI1" ...
    ##  $ Control.gRNA      : int  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ Location_1536     : chr  "A1" "A2" "A3" "A4" ...
    ##  $ Round_ID          : chr  "2nd_round" "2nd_round" "2nd_round" "2nd_round" ...
    ##  $ CTRL_Y_ABS        : num  8533790 6626182 5943280 5455589 5133404 ...
    ##  $ CTRL_GT_ABS       : num  2.53 2.48 2.51 2.56 2.53 ...
    ##  $ AA_Y_ABS          : num  2090439 2241914 1861277 1920070 1957912 ...
    ##  $ AA_GT_ABS         : num  8.92 8.43 9.19 8.91 8.87 ...
    ##  $ CTRL_Y_NORM       : num  0.755 0.39 0.233 0.484 0.396 ...
    ##  $ CTRL_GT_NORM      : num  0.0505 0.0204 0.0389 -0.0173 -0.0327 ...
    ##  $ AA_Y_NORM         : num  0.373 0.474 0.205 0.386 0.415 ...
    ##  $ AA_GT_NORM        : num  0.019 -0.0614 0.0621 -0.2135 -0.219 ...

## IMPORT RESULTS FROM ROUND1

The results from Round1 is already compiled to a .csv file in
COMPILED\_DATA folder **Results 1st Round** :
20190903\_CRISPRi\_Screen\_aa\_1st\_round.csv

Import the dataset and label as data from ROUND1

``` r
First_round <- read.csv("COMPILED_DATA/20190903_CRISPRi_Screen_aa_1st_round.csv", 
                        na.strings = c("#N/A", "NoGrowth"), 
                        stringsAsFactors = FALSE)
R <- rep("1st_round", 36864)
Round_ID <- data.frame(R, stringsAsFactors = FALSE)
whole_data_R1 <- cbind(Metadata_CRISPRi, Round_ID, First_round[, 12:19])
colnames(whole_data_R1)[12] <- "Round_ID"
str(whole_data_R1)
```

    ## 'data.frame':    36864 obs. of  20 variables:
    ##  $ SL_No             : num  1 2 3 4 5 6 7 8 9 10 ...
    ##  $ gRNA_name         : chr  "FBP26-TRg-1" "FBP26-TRg-1" "HMI1-NRg-1" "HMI1-NRg-1" ...
    ##  $ Seq               : chr  "GCTTATCATACATTTACATC" "GCTTATCATACATTTACATC" "AAAAATTCTGACACATCACA" "AAAAATTCTGACACATCACA" ...
    ##  $ SOURCEPLATEID     : chr  "R2877.H.001" "R2877.H.001" "R2877.H.001" "R2877.H.001" ...
    ##  $ SOURCEDENSITY     : chr  "384A" "384A" "384A" "384A" ...
    ##  $ SOURCECOLONYCOLUMN: int  1 1 2 2 3 3 4 4 5 5 ...
    ##  $ SOURCECOLONYROW   : chr  "A" "A" "A" "A" ...
    ##  $ border            : logi  TRUE TRUE TRUE TRUE TRUE TRUE ...
    ##  $ GENE              : chr  "FBP26" "FBP26" "HMI1" "HMI1" ...
    ##  $ Control.gRNA      : int  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ Location_1536     : chr  "A1" "A2" "A3" "A4" ...
    ##  $ Round_ID          : chr  "1st_round" "1st_round" "1st_round" "1st_round" ...
    ##  $ CTRL_Y_ABS        : num  7046465 5380541 4841666 4517281 4209174 ...
    ##  $ CTRL_GT_ABS       : num  3.15 3.64 3.18 3.04 3.13 ...
    ##  $ AA_Y_ABS          : num  833657 844666 708042 734725 717147 ...
    ##  $ AA_GT_ABS         : num  13.1 13.6 14.5 14.7 13.6 ...
    ##  $ CTRL_Y_NORM       : num  0.899 0.51 0.358 0.565 0.463 ...
    ##  $ CTRL_GT_NORM      : num  -0.10076 0.10855 -0.08699 0.00504 0.0462 ...
    ##  $ AA_Y_NORM         : num  -0.362 -0.343 -0.598 -0.642 -0.677 ...
    ##  $ AA_GT_NORM        : num  0.273 0.331 0.423 0.425 0.313 ...

## COMBINE THE DATASETS of ROUND 1 AND 2

``` r
whole_data_CRISPRi_aa <- rbind(whole_data_R1, whole_data_R2)
```

# SCAN-O-MATIC PHENOMICS ANALYSIS

In this study most of the downstream analysis was performed using the
phenotype Generation\_time(GT)

## ESTIMATE THE LOG PHENOTYPIC INDEX (LPI) VALUES

LPI of strain is the difference of its normalized Generation\_Time(GT) /
Yield(Y) (LSC, see [IMPORT SCAN-O-MATIC RAW
DATA](#import-scan-o-matic-raw-data)) on acetic acid stress plate to the
basal condition. It gives a relative estimate of how a strain performed
under acetic acid stress relative to the basal condition.

i.e. LPI\_GT = LSC\_GT\_Acetic\_Acid - LSC\_GT\_Basal

``` r
whole_data_CRISPRi_aa[, 21] <- whole_data_CRISPRi_aa[, 19]-whole_data_CRISPRi_aa[, 17]
whole_data_CRISPRi_aa[, 22] <- whole_data_CRISPRi_aa[, 20]-whole_data_CRISPRi_aa[, 18]
colnames(whole_data_CRISPRi_aa)[21] <- "LPI_Y"
colnames(whole_data_CRISPRi_aa)[22] <- "LPI_GT"
```

## PERFORM PLATE-WISE BATCH CORRECTION

Plate-wise batch correction was conducted by subtracting the median of
LSC GT values of all the individual colonies on a plate from the
individual LSC GT values of the colonies growing on that plate.

i.e. if strainX is growing in Basal condition on plate Z, the corrected
LSC\_GT value for strainX in the Basal condition is the following;

  - LSC\_GT\_Basal\_Corrected<sub>strainX</sub> =
    (LSC\_GT\_Basal<sub>strainX</sub>) - Median(LSC\_GT
    Basal<sub>PlateZ</sub>)

<!-- end list -->

``` r
plate_ID <- as.character(unique(whole_data_CRISPRi_aa$SOURCEPLATEID))
whole_data_CRISPRi_aa_corrected <- whole_data_CRISPRi_aa
med_LogLSCctrl_RND1_GT <- vector(mode = "integer", length = 0)
med_LogLSCaa_RND1_GT <- vector(mode = "integer", length = 0)
med_LogLSCctrl_RND2_GT <- vector(mode = "integer", length = 0)
med_LogLSCaa_RND2_GT <- vector(mode = "integer", length = 0)
med_LogLSCctrl_RND1_Y <- vector(mode = "integer", length = 0)
med_LogLSCaa_RND1_Y <- vector(mode = "integer", length = 0)
med_LogLSCctrl_RND2_Y <- vector(mode = "integer", length = 0)
med_LogLSCaa_RND2_Y <- vector(mode = "integer", length = 0)

for(i in 1:24){
med_LogLSCctrl_RND1_GT[i] <- median(whole_data_CRISPRi_aa_corrected$CTRL_GT_NORM[which(whole_data_CRISPRi_aa_corrected$SOURCEPLATEID==plate_ID[i]
                                                                                       & !is.na(whole_data_CRISPRi_aa_corrected$CTRL_GT_NORM) 
                                                                                         & whole_data_CRISPRi_aa_corrected$Round_ID=="1st_round")])

med_LogLSCaa_RND1_GT[i] <- median(whole_data_CRISPRi_aa_corrected$AA_GT_NORM[which(whole_data_CRISPRi_aa_corrected$SOURCEPLATEID==plate_ID[i]
                                                                                     & !is.na(whole_data_CRISPRi_aa_corrected$AA_GT_NORM) 
                                                                                     & whole_data_CRISPRi_aa_corrected$Round_ID=="1st_round")])

med_LogLSCctrl_RND2_GT[i] <- median(whole_data_CRISPRi_aa_corrected$CTRL_GT_NORM[which(whole_data_CRISPRi_aa_corrected$SOURCEPLATEID==plate_ID[i]
                                                                                         & !is.na(whole_data_CRISPRi_aa_corrected$CTRL_GT_NORM) 
                                                                                         & whole_data_CRISPRi_aa_corrected$Round_ID=="2nd_round")])

med_LogLSCaa_RND2_GT[i] <- median(whole_data_CRISPRi_aa_corrected$AA_GT_NORM[which(whole_data_CRISPRi_aa_corrected$SOURCEPLATEID==plate_ID[i]
                                                                                     & !is.na(whole_data_CRISPRi_aa_corrected$AA_GT_NORM) 
                                                                                     & whole_data_CRISPRi_aa_corrected$Round_ID=="2nd_round")])

med_LogLSCctrl_RND1_Y[i] <- median(whole_data_CRISPRi_aa_corrected$CTRL_Y_NORM[which(whole_data_CRISPRi_aa_corrected$SOURCEPLATEID==plate_ID[i]
                                                                                       & !is.na(whole_data_CRISPRi_aa_corrected$CTRL_Y_NORM) 
                                                                                       & whole_data_CRISPRi_aa_corrected$Round_ID=="1st_round")])

med_LogLSCaa_RND1_Y[i] <- median(whole_data_CRISPRi_aa_corrected$AA_Y_NORM[which(whole_data_CRISPRi_aa_corrected$SOURCEPLATEID==plate_ID[i]
                                                                                   & !is.na(whole_data_CRISPRi_aa_corrected$AA_Y_NORM) 
                                                                                   & whole_data_CRISPRi_aa_corrected$Round_ID=="1st_round")])

med_LogLSCctrl_RND2_Y[i] <- median(whole_data_CRISPRi_aa_corrected$CTRL_Y_NORM[which(whole_data_CRISPRi_aa_corrected$SOURCEPLATEID==plate_ID[i]
                                                                                       & !is.na(whole_data_CRISPRi_aa_corrected$CTRL_Y_NORM) 
                                                                                       & whole_data_CRISPRi_aa_corrected$Round_ID=="2nd_round")])

med_LogLSCaa_RND2_Y[i] <- median(whole_data_CRISPRi_aa_corrected$AA_Y_NORM[which(whole_data_CRISPRi_aa_corrected$SOURCEPLATEID==plate_ID[i]
                                                                                   & !is.na(whole_data_CRISPRi_aa_corrected$AA_Y_NORM) 
                                                                                   & whole_data_CRISPRi_aa_corrected$Round_ID=="2nd_round")])
  
whole_data_CRISPRi_aa_corrected[which(whole_data_CRISPRi_aa_corrected$SOURCEPLATEID==plate_ID[i] &
                                          whole_data_CRISPRi_aa_corrected$Round_ID=="1st_round") , 23] <- whole_data_CRISPRi_aa_corrected[which(whole_data_CRISPRi_aa_corrected$SOURCEPLATEID==plate_ID[i] &
                                                                                                                                                  whole_data_CRISPRi_aa_corrected$Round_ID=="1st_round"), 17] - med_LogLSCctrl_RND1_Y[i]
  whole_data_CRISPRi_aa_corrected[which(whole_data_CRISPRi_aa_corrected$SOURCEPLATEID==plate_ID[i] &
                                          whole_data_CRISPRi_aa_corrected$Round_ID=="1st_round") , 24] <- whole_data_CRISPRi_aa_corrected[which(whole_data_CRISPRi_aa_corrected$SOURCEPLATEID==plate_ID[i] &
                                                                                                                                                  whole_data_CRISPRi_aa_corrected$Round_ID=="1st_round"), 18] - med_LogLSCctrl_RND1_GT[i]
  whole_data_CRISPRi_aa_corrected[which(whole_data_CRISPRi_aa_corrected$SOURCEPLATEID==plate_ID[i] &
                                          whole_data_CRISPRi_aa_corrected$Round_ID=="2nd_round") , 23] <- whole_data_CRISPRi_aa_corrected[which(whole_data_CRISPRi_aa_corrected$SOURCEPLATEID==plate_ID[i] &
                                                                                                                                                  whole_data_CRISPRi_aa_corrected$Round_ID=="2nd_round"), 17] - med_LogLSCctrl_RND2_Y[i]
  whole_data_CRISPRi_aa_corrected[which(whole_data_CRISPRi_aa_corrected$SOURCEPLATEID==plate_ID[i] &
                                          whole_data_CRISPRi_aa_corrected$Round_ID=="2nd_round") , 24] <- whole_data_CRISPRi_aa_corrected[which(whole_data_CRISPRi_aa_corrected$SOURCEPLATEID==plate_ID[i] &
                                                                                                                                                  whole_data_CRISPRi_aa_corrected$Round_ID=="2nd_round"), 18] - med_LogLSCctrl_RND2_GT[i]
  whole_data_CRISPRi_aa_corrected[which(whole_data_CRISPRi_aa_corrected$SOURCEPLATEID==plate_ID[i] &
                                          whole_data_CRISPRi_aa_corrected$Round_ID=="1st_round") , 25] <- whole_data_CRISPRi_aa_corrected[which(whole_data_CRISPRi_aa_corrected$SOURCEPLATEID==plate_ID[i] &
                                                                                                                                                  whole_data_CRISPRi_aa_corrected$Round_ID=="1st_round"), 19] - med_LogLSCaa_RND1_Y[i]
  whole_data_CRISPRi_aa_corrected[which(whole_data_CRISPRi_aa_corrected$SOURCEPLATEID==plate_ID[i] &
                                          whole_data_CRISPRi_aa_corrected$Round_ID=="1st_round") , 26] <- whole_data_CRISPRi_aa_corrected[which(whole_data_CRISPRi_aa_corrected$SOURCEPLATEID==plate_ID[i] &
                                                                                                                                                  whole_data_CRISPRi_aa_corrected$Round_ID=="1st_round"), 20] - med_LogLSCaa_RND1_GT[i]
  whole_data_CRISPRi_aa_corrected[which(whole_data_CRISPRi_aa_corrected$SOURCEPLATEID==plate_ID[i] &
                                          whole_data_CRISPRi_aa_corrected$Round_ID=="2nd_round") , 25] <- whole_data_CRISPRi_aa_corrected[which(whole_data_CRISPRi_aa_corrected$SOURCEPLATEID==plate_ID[i] &
                                                                                                                                                  whole_data_CRISPRi_aa_corrected$Round_ID=="2nd_round"), 19] - med_LogLSCaa_RND2_Y[i]
  whole_data_CRISPRi_aa_corrected[which(whole_data_CRISPRi_aa_corrected$SOURCEPLATEID==plate_ID[i] &
                                          whole_data_CRISPRi_aa_corrected$Round_ID=="2nd_round") , 26] <- whole_data_CRISPRi_aa_corrected[which(whole_data_CRISPRi_aa_corrected$SOURCEPLATEID==plate_ID[i] &
                                                                                                                                                  whole_data_CRISPRi_aa_corrected$Round_ID=="2nd_round"), 20] - med_LogLSCaa_RND2_GT[i]
}
```

## ESTIMATE THE BATCH CORRECTED LOG PHENOTYPIC INDEX (LPI) VALUES

Estimate the corrected LPI values (see [ESTIMATE THE LOG PHENOTYPIC
INDEX (LPI) VALUES](#estimate-the-log-phenotypic-index-lpi-values))
based on the corrected LSC values

i.e. LPI\_GT<sub>corrected</sub> =
LSC\_GT\_Acetic\_Acid<sub>corrected</sub> -
LSC\_GT\_Basal<sub>corrected</sub>

**Estimate the corrected LPI\_Y**

``` r
whole_data_CRISPRi_aa_corrected[, 27] <- whole_data_CRISPRi_aa_corrected[, 25] - whole_data_CRISPRi_aa_corrected[, 23]
```

**Estimate the corrected LPI\_GT**

``` r
whole_data_CRISPRi_aa_corrected[, 28] <- whole_data_CRISPRi_aa_corrected[, 26] - whole_data_CRISPRi_aa_corrected[, 24] 
```

## SETTING THE NAMES OF THE NEW COLUMNS

``` r
colnm <- colnames(whole_data_CRISPRi_aa)[17:22]
colnm <- paste0(colnm, "_CR")
colnames(whole_data_CRISPRi_aa_corrected)[23:28] <- colnm
str(whole_data_CRISPRi_aa_corrected)
```

    ## 'data.frame':    73728 obs. of  28 variables:
    ##  $ SL_No             : num  1 2 3 4 5 6 7 8 9 10 ...
    ##  $ gRNA_name         : chr  "FBP26-TRg-1" "FBP26-TRg-1" "HMI1-NRg-1" "HMI1-NRg-1" ...
    ##  $ Seq               : chr  "GCTTATCATACATTTACATC" "GCTTATCATACATTTACATC" "AAAAATTCTGACACATCACA" "AAAAATTCTGACACATCACA" ...
    ##  $ SOURCEPLATEID     : chr  "R2877.H.001" "R2877.H.001" "R2877.H.001" "R2877.H.001" ...
    ##  $ SOURCEDENSITY     : chr  "384A" "384A" "384A" "384A" ...
    ##  $ SOURCECOLONYCOLUMN: int  1 1 2 2 3 3 4 4 5 5 ...
    ##  $ SOURCECOLONYROW   : chr  "A" "A" "A" "A" ...
    ##  $ border            : logi  TRUE TRUE TRUE TRUE TRUE TRUE ...
    ##  $ GENE              : chr  "FBP26" "FBP26" "HMI1" "HMI1" ...
    ##  $ Control.gRNA      : int  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ Location_1536     : chr  "A1" "A2" "A3" "A4" ...
    ##  $ Round_ID          : chr  "1st_round" "1st_round" "1st_round" "1st_round" ...
    ##  $ CTRL_Y_ABS        : num  7046465 5380541 4841666 4517281 4209174 ...
    ##  $ CTRL_GT_ABS       : num  3.15 3.64 3.18 3.04 3.13 ...
    ##  $ AA_Y_ABS          : num  833657 844666 708042 734725 717147 ...
    ##  $ AA_GT_ABS         : num  13.1 13.6 14.5 14.7 13.6 ...
    ##  $ CTRL_Y_NORM       : num  0.899 0.51 0.358 0.565 0.463 ...
    ##  $ CTRL_GT_NORM      : num  -0.10076 0.10855 -0.08699 0.00504 0.0462 ...
    ##  $ AA_Y_NORM         : num  -0.362 -0.343 -0.598 -0.642 -0.677 ...
    ##  $ AA_GT_NORM        : num  0.273 0.331 0.423 0.425 0.313 ...
    ##  $ LPI_Y             : num  -1.262 -0.854 -0.956 -1.207 -1.14 ...
    ##  $ LPI_GT            : num  0.374 0.223 0.51 0.42 0.266 ...
    ##  $ CTRL_Y_NORM_CR    : num  0.918 0.529 0.376 0.583 0.482 ...
    ##  $ CTRL_GT_NORM_CR   : num  -0.0877 0.1216 -0.074 0.0181 0.0592 ...
    ##  $ AA_Y_NORM_CR      : num  -0.0736 -0.0547 -0.3092 -0.3534 -0.3883 ...
    ##  $ AA_GT_NORM_CR     : num  0.19 0.248 0.34 0.342 0.229 ...
    ##  $ LPI_Y_CR          : num  -0.991 -0.583 -0.686 -0.937 -0.87 ...
    ##  $ LPI_GT_CR         : num  0.278 0.127 0.414 0.323 0.17 ...

## EXTRACT ONLY THE BATCH CORRECTED COLUMNS

``` r
whole_data_CRISPRi_aa_2 <- whole_data_CRISPRi_aa_corrected[, c(1:16, 23:28)]
colnames(whole_data_CRISPRi_aa_2)[17:22] <- colnames(whole_data_CRISPRi_aa)[17:22]
str(whole_data_CRISPRi_aa_2)
```

    ## 'data.frame':    73728 obs. of  22 variables:
    ##  $ SL_No             : num  1 2 3 4 5 6 7 8 9 10 ...
    ##  $ gRNA_name         : chr  "FBP26-TRg-1" "FBP26-TRg-1" "HMI1-NRg-1" "HMI1-NRg-1" ...
    ##  $ Seq               : chr  "GCTTATCATACATTTACATC" "GCTTATCATACATTTACATC" "AAAAATTCTGACACATCACA" "AAAAATTCTGACACATCACA" ...
    ##  $ SOURCEPLATEID     : chr  "R2877.H.001" "R2877.H.001" "R2877.H.001" "R2877.H.001" ...
    ##  $ SOURCEDENSITY     : chr  "384A" "384A" "384A" "384A" ...
    ##  $ SOURCECOLONYCOLUMN: int  1 1 2 2 3 3 4 4 5 5 ...
    ##  $ SOURCECOLONYROW   : chr  "A" "A" "A" "A" ...
    ##  $ border            : logi  TRUE TRUE TRUE TRUE TRUE TRUE ...
    ##  $ GENE              : chr  "FBP26" "FBP26" "HMI1" "HMI1" ...
    ##  $ Control.gRNA      : int  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ Location_1536     : chr  "A1" "A2" "A3" "A4" ...
    ##  $ Round_ID          : chr  "1st_round" "1st_round" "1st_round" "1st_round" ...
    ##  $ CTRL_Y_ABS        : num  7046465 5380541 4841666 4517281 4209174 ...
    ##  $ CTRL_GT_ABS       : num  3.15 3.64 3.18 3.04 3.13 ...
    ##  $ AA_Y_ABS          : num  833657 844666 708042 734725 717147 ...
    ##  $ AA_GT_ABS         : num  13.1 13.6 14.5 14.7 13.6 ...
    ##  $ CTRL_Y_NORM       : num  0.918 0.529 0.376 0.583 0.482 ...
    ##  $ CTRL_GT_NORM      : num  -0.0877 0.1216 -0.074 0.0181 0.0592 ...
    ##  $ AA_Y_NORM         : num  -0.0736 -0.0547 -0.3092 -0.3534 -0.3883 ...
    ##  $ AA_GT_NORM        : num  0.19 0.248 0.34 0.342 0.229 ...
    ##  $ LPI_Y             : num  -0.991 -0.583 -0.686 -0.937 -0.87 ...
    ##  $ LPI_GT            : num  0.278 0.127 0.414 0.323 0.17 ...

## CONSTRUCT A NEW DATA STRUCTURE

Construct a new data structure where data from each strain (have a
unique sgRNA) is in a separate row and the replicates from first and
second round are side by side. Add also the mean, median and standard
deviation statistics for each phenotype

### REMOVE ROWS
