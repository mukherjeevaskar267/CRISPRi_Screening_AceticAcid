#CRISPRi_Screening_AceticAcid_rscript

AA_titration_data <- read.csv("COMPILED_DATA/20210120_AA_titration_absolute_compiled.csv", na.strings = "NoGrowth")

library(reshape)
AA_titration_data_reshape <- reshape(data=AA_titration_data, idvar="gRNA_name",
                                     varying = colnames(AA_titration_data)[3:7],
                                     v.name=c("Generation_time"),
                                     new.row.names = 1:30000,
                                     direction="long",
                                     timevar = "Condition",
                                     times = colnames(AA_titration_data)[3:7])
library(ggplot2)
library(reshape)
library(ggridges)
plt0 <- ggplot(AA_titration_data_reshape, aes(x = Generation_time, y = Condition, height = stat(density))) + 
  geom_density_ridges2(stat = "binline", bins = 200, scale = 2, draw_baseline = FALSE)+
  theme_ridges()
suppressWarnings(print(plt0))

Metadata_CRISPRi <- read.csv("COMPILED_DATA/library_keyfile1536.csv", na.strings = "#N/A", stringsAsFactors = FALSE)
str(Metadata_CRISPRi)

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

data_Ctrl_Abs_Trim <- data_Ctrl_Abs[, 14:15]
colnames(data_Ctrl_Abs_Trim) <- c("CTRL_Y_ABS", "CTRL_GT_ABS")
str(data_Ctrl_Abs_Trim)

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

data_Ctrl_Norm_Trim <- data_Ctrl_Norm[, 4:5]
colnames(data_Ctrl_Norm_Trim) <- c("CTRL_Y_NORM", "CTRL_GT_NORM")

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

First_round <- read.csv("COMPILED_DATA/20190903_CRISPRi_Screen_aa_1st_round.csv", 
                        na.strings = c("#N/A", "NoGrowth"), 
                        stringsAsFactors = FALSE)
R <- rep("1st_round", 36864)
Round_ID <- data.frame(R, stringsAsFactors = FALSE)
whole_data_R1 <- cbind(Metadata_CRISPRi, Round_ID, First_round[, 12:19])
colnames(whole_data_R1)[12] <- "Round_ID"
str(whole_data_R1)

whole_data_CRISPRi_aa <- rbind(whole_data_R1, whole_data_R2)

whole_data_CRISPRi_aa[, 21] <- whole_data_CRISPRi_aa[, 19]-whole_data_CRISPRi_aa[, 17]
whole_data_CRISPRi_aa[, 22] <- whole_data_CRISPRi_aa[, 20]-whole_data_CRISPRi_aa[, 18]
colnames(whole_data_CRISPRi_aa)[21] <- "LPI_Y"
colnames(whole_data_CRISPRi_aa)[22] <- "LPI_GT"

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



whole_data_CRISPRi_aa_corrected[, 27] <- whole_data_CRISPRi_aa_corrected[, 25] - whole_data_CRISPRi_aa_corrected[, 23]

whole_data_CRISPRi_aa_corrected[, 28] <- whole_data_CRISPRi_aa_corrected[, 26] - whole_data_CRISPRi_aa_corrected[, 24]

colnm <- colnames(whole_data_CRISPRi_aa)[17:22]
colnm <- paste0(colnm, "_CR")
colnames(whole_data_CRISPRi_aa_corrected)[23:28] <- colnm
str(whole_data_CRISPRi_aa_corrected)

whole_data_CRISPRi_aa_2 <- whole_data_CRISPRi_aa_corrected[, c(1:16, 23:28)]
colnames(whole_data_CRISPRi_aa_2)[17:22] <- colnames(whole_data_CRISPRi_aa)[17:22]
str(whole_data_CRISPRi_aa_2)

Data_CRISPRi_aa <- subset(whole_data_CRISPRi_aa_2, whole_data_CRISPRi_aa_2$gRNA_name!="SP_Ctrl_CC23")

df_unique_sgRNA <- data.frame(table(Data_CRISPRi_aa$gRNA_name))

R1<-vector(mode = "integer", length = 0)
R2<-vector(mode = "integer", length = 0)
test2<-data.frame()
n<-nrow(df_unique_sgRNA)
for(i in 1:n){
  R1 <- which(Data_CRISPRi_aa$gRNA_name==df_unique_sgRNA$Var1[i] & Data_CRISPRi_aa$Round_ID=="1st_round")
  R2 <- which(Data_CRISPRi_aa$gRNA_name==df_unique_sgRNA$Var1[i] & Data_CRISPRi_aa$Round_ID=="2nd_round")
  test1 <- Data_CRISPRi_aa[c(R1, R2), ]
  test2[i, c(1:8)]<-test1[1, c(2:4, 6:7, 9:11)]
  test2[i, c(9:14)] <- test1$CTRL_GT_NORM
  test2[i, 15] <- mean(test1$CTRL_GT_NORM[1:3])
  test2[i, 16] <- mean(test1$CTRL_GT_NORM[4:6])
  test2[i, 17] <- sd(test1$CTRL_GT_NORM[1:3])
  test2[i, 18] <- sd(test1$CTRL_GT_NORM[4:6])
  test2[i, 19] <- mean(test1$CTRL_GT_NORM[1:6])
  test2[i, 20] <- median(test1$CTRL_GT_NORM[1:6])
  test2[i, 21] <- sd(test1$CTRL_GT_NORM[1:6])
  test2[i, c(22:27)] <- test1$AA_GT_NORM
  test2[i, 28] <- mean(test1$AA_GT_NORM[1:3])
  test2[i, 29] <- mean(test1$AA_GT_NORM[4:6])
  test2[i, 30] <- sd(test1$AA_GT_NORM[1:3])
  test2[i, 31] <- sd(test1$AA_GT_NORM[4:6])
  test2[i, 32] <- mean(test1$AA_GT_NORM[1:6])
  test2[i, 33] <- median(test1$AA_GT_NORM[1:6])
  test2[i, 34] <- sd(test1$AA_GT_NORM[1:6])
  test2[i, c(35:40)] <- test1$LPI_GT
  test2[i, 41] <- mean(test1$LPI_GT[1:3])
  test2[i, 42] <- mean(test1$LPI_GT[4:6])
  test2[i, 43] <- sd(test1$LPI_GT[1:3])
  test2[i, 44] <- sd(test1$LPI_GT[4:6])
  test2[i, 45] <- mean(test1$LPI_GT[1:6])
  test2[i, 46] <- median(test1$LPI_GT[1:6])
  test2[i, 47] <- sd(test1$LPI_GT[1:6])
  test2[i, c(48:53)] <- test1$CTRL_Y_NORM
  test2[i, 54] <- mean(test1$CTRL_Y_NORM[1:3])
  test2[i, 55] <- mean(test1$CTRL_Y_NORM[4:6])
  test2[i, 56] <- sd(test1$CTRL_Y_NORM[1:3])
  test2[i, 57] <- sd(test1$CTRL_Y_NORM[4:6])
  test2[i, 58] <- mean(test1$CTRL_Y_NORM[1:6])
  test2[i, 59] <- median(test1$CTRL_Y_NORM[1:6])
  test2[i, 60] <- sd(test1$CTRL_Y_NORM[1:6])
  test2[i, c(61:66)] <- test1$AA_Y_NORM
  test2[i, 67] <- mean(test1$AA_Y_NORM[1:3])
  test2[i, 68] <- mean(test1$AA_Y_NORM[4:6])
  test2[i, 69] <- sd(test1$AA_Y_NORM[1:3])
  test2[i, 70] <- sd(test1$AA_Y_NORM[4:6])
  test2[i, 71] <- mean(test1$AA_Y_NORM[1:6])
  test2[i, 72] <- median(test1$AA_Y_NORM[1:6])
  test2[i, 73] <- sd(test1$AA_Y_NORM[1:6])
  test2[i, c(74:79)] <- test1$LPI_Y
  test2[i, 80] <- mean(test1$LPI_Y[1:3])
  test2[i, 81] <- mean(test1$LPI_Y[4:6])
  test2[i, 82] <- sd(test1$LPI_Y[1:3])
  test2[i, 83] <- sd(test1$LPI_Y[4:6])
  test2[i, 84] <- mean(test1$LPI_Y[1:6])
  test2[i, 85] <- median(test1$LPI_Y[1:6])
  test2[i, 86] <- sd(test1$LPI_Y[1:6])
}

column_names <- read.table("COMPILED_DATA/Column_names.txt", header = FALSE, sep = "\t", as.is = TRUE)
colnames(test2) <- column_names$V1
Analysis_CRISPRi_aa_Complete <- test2
str(Analysis_CRISPRi_aa_Complete)

Analysis_Final <- Analysis_CRISPRi_aa_Complete

for(i in 1:nrow(Analysis_Final)){
  x1 <- as.numeric(Analysis_Final[i, 9:11][which(!is.na(Analysis_Final[i, 9:11]))])
  x2 <- as.numeric(Analysis_Final[i, 12:14][which(!is.na(Analysis_Final[i, 12:14]))])
  if(length(x1)==0){
    Analysis_Final$CTRL_GT_RND1_MEAN[i] <- NA
  } else{
    Analysis_Final$CTRL_GT_RND1_MEAN[i] <- as.numeric(mean(x1))
  }
  if(length(x2)==0){
    Analysis_Final$CTRL_GT_RND2_MEAN[i] <- NA
  } else{
    Analysis_Final$CTRL_GT_RND2_MEAN[i] <- as.numeric(mean(x2))
  }
  if(sum(is.na(c(Analysis_Final$CTRL_GT_RND1_MEAN[i], Analysis_Final$CTRL_GT_RND2_MEAN[i])))==0){
    Analysis_Final$CTRL_GT_RND1_2_MEAN[i] <- as.numeric(mean(c(Analysis_Final$CTRL_GT_RND1_MEAN[i], Analysis_Final$CTRL_GT_RND2_MEAN[i])))
    Analysis_Final[i, 87] <- as.numeric(sd(c(Analysis_Final$CTRL_GT_RND1_MEAN[i], Analysis_Final$CTRL_GT_RND2_MEAN[i])))
  } else{
    Analysis_Final$CTRL_GT_RND1_2_MEAN[i] <- NA
    Analysis_Final[i, 87] <- NA
  }
}
colnames(Analysis_Final)[87] <- "CTRL_GT_MEAN_RND1_2_SD"

for(i in 1:nrow(Analysis_Final)){
  x1 <- as.numeric(Analysis_Final[i, 22:24][which(!is.na(Analysis_Final[i, 22:24]))])
  x2 <- as.numeric(Analysis_Final[i, 25:27][which(!is.na(Analysis_Final[i, 25:27]))])
  if(length(x1)==0){
    Analysis_Final$AA_GT_RND1_MEAN[i] <- NA
  } else{
    Analysis_Final$AA_GT_RND1_MEAN[i] <- as.numeric(mean(x1))
  }
  if(length(x2)==0){
    Analysis_Final$AA_GT_RND2_MEAN[i] <- NA
  } else{
    Analysis_Final$AA_GT_RND2_MEAN[i] <- as.numeric(mean(x2))
  }
  if(sum(is.na(c(Analysis_Final$AA_GT_RND1_MEAN[i], Analysis_Final$AA_GT_RND2_MEAN[i])))==0){
    Analysis_Final$AA_GT_RND1_2_MEAN[i] <- as.numeric(mean(c(Analysis_Final$AA_GT_RND1_MEAN[i], Analysis_Final$AA_GT_RND2_MEAN[i])))
    Analysis_Final[i, 88] <- as.numeric(sd(c(Analysis_Final$AA_GT_RND1_MEAN[i], Analysis_Final$AA_GT_RND2_MEAN[i])))
  } else{
    Analysis_Final$AA_GT_RND1_2_MEAN[i] <- NA
    Analysis_Final[i, 88] <- NA
  }
}
colnames(Analysis_Final)[88] <- "AA_GT_MEAN_RND1_2_SD"

for(i in 1:nrow(Analysis_Final)){
  x1 <- as.numeric(Analysis_Final[i, 35:37][which(!is.na(Analysis_Final[i, 35:37]))])
  x2 <- as.numeric(Analysis_Final[i, 38:40][which(!is.na(Analysis_Final[i, 38:40]))])
  if(length(x1)==0){
    Analysis_Final$LPI_GT_RND1_MEAN[i] <- NA
  } else{
    Analysis_Final$LPI_GT_RND1_MEAN[i] <- as.numeric(mean(x1))
  }
  if(length(x2)==0){
    Analysis_Final$LPI_GT_RND2_MEAN[i] <- NA
  } else{
    Analysis_Final$LPI_GT_RND2_MEAN[i] <- as.numeric(mean(x2))
  }
  if(sum(is.na(c(Analysis_Final$LPI_GT_RND1_MEAN[i], Analysis_Final$LPI_GT_RND2_MEAN[i])))==0){
    Analysis_Final$LPI_GT_RND1_2_MEAN[i] <- as.numeric(mean(c(Analysis_Final$LPI_GT_RND1_MEAN[i], Analysis_Final$LPI_GT_RND2_MEAN[i])))
    Analysis_Final[i, 89] <- as.numeric(sd(c(Analysis_Final$LPI_GT_RND1_MEAN[i], Analysis_Final$LPI_GT_RND2_MEAN[i])))
  } else{
    Analysis_Final$LPI_GT_RND1_2_MEAN[i] <- NA
    Analysis_Final[i, 89] <- NA
  }
}
colnames(Analysis_Final)[89] <- "LPI_GT_MEAN_RND1_2_SD"

box_stat_LPI_GT_R1_2_mean <- boxplot(Analysis_Final$LPI_GT_RND1_2_MEAN, cex=0.3)
box_stat_LPI_GT_R1_2_mean$stats

Intermediate_50 <- Analysis_Final$LPI_GT_RND1_2_MEAN[which(Analysis_Final$LPI_GT_RND1_2_MEAN>=-0.02428792
                                                           &Analysis_Final$LPI_GT_RND1_2_MEAN<=0.07255828)]
summary(Intermediate_50)

for(i in 1:nrow(Analysis_Final)){
  if(sum(is.na(c(Analysis_Final$LPI_GT_RND1_MEAN[i], Analysis_Final$LPI_GT_RND2_MEAN[i])))==0){
    P_value <- t.test(Intermediate_50, c(Analysis_Final$LPI_GT_RND1_MEAN[i], Analysis_Final$LPI_GT_RND2_MEAN[i]))
    Analysis_Final[i, 90] <- P_value$p.value
  } else{
    Analysis_Final[i, 90] <- NA
  }
}
colnames(Analysis_Final)[90] <- "P_value_M1"

Analysis_Final[which(!is.na(Analysis_Final$P_value_M1)), 91] <- p.adjust(Analysis_Final$P_value_M1[which(!is.na(Analysis_Final$P_value_M1))], 
                                                                         method = "BH", 
                                                                         n = length(Analysis_Final$P_value_M1[which(!is.na(Analysis_Final$P_value_M1))]))
colnames(Analysis_Final)[91] <- "P.adjusted_M1"

length(Analysis_Final$P_value_M1[which(Analysis_Final$P_value_M1<=0.05)])
length(Analysis_Final$P.adjusted_M1[which(Analysis_Final$P.adjusted_M1<=0.05)])
length(Analysis_Final$P_value_M1[which(Analysis_Final$P_value_M1<=0.1)])
length(Analysis_Final$P.adjusted_M1[which(Analysis_Final$P.adjusted_M1<=0.1)])

par(mfrow=c(2,2))
hist(Analysis_Final$P_value_M1,
     breaks = 20,
     xlab = "P-value", 
     ylab = "Frequency", 
     main = "P-values of all strains", 
     col = "skyblue",
     xlim = c(0, 1),
     ylim = c(0, 1000))
hist(Analysis_Final$P_value_M1[which(Analysis_Final$Control.gRNA==1)],
     breaks = 10,
     xlab = "P-value", 
     ylab = "Frequency", 
     main = "P-values of control strains", 
     col = "skyblue",
     xlim = c(0, 1),
     ylim = c(0, 5))
hist(Analysis_Final$P.adjusted_M1,
     breaks = 20,
     xlab = "P.adj", 
     ylab = "Frequency", 
     main = "P.adjusted values of all strains", 
     col = "skyblue",
     xlim = c(0, 1),
     ylim = c(0, 1000))
hist(Analysis_Final$P.adjusted_M1[which(Analysis_Final$Control.gRNA==1)],
     breaks = 2,
     xlab = "P.adj", 
     ylab = "Frequency", 
     main = "P.adjusted values of control strains", 
     col = "skyblue",
     xlim = c(0, 1),
     ylim = c(0, 10))

Analysis_Final_2 <- Analysis_Final
str(Analysis_Final_2)

CRISPRi_Ctrl_Round1 <- whole_data_CRISPRi_aa_2$LPI_GT[which(whole_data_CRISPRi_aa_2$Control.gRNA == 1 
                                                            & whole_data_CRISPRi_aa_2$Round_ID=="1st_round")]
summary(CRISPRi_Ctrl_Round1)

CRISPRi_Ctrl_Round2 <- whole_data_CRISPRi_aa_2$LPI_GT[which(whole_data_CRISPRi_aa_2$Control.gRNA == 1 
                                                            & whole_data_CRISPRi_aa_2$Round_ID=="2nd_round")]
summary(CRISPRi_Ctrl_Round2)

for(i in 1:nrow(Analysis_Final_2)){
  test1 <- t(Analysis_Final_2[i, 35:37])
  test2 <- t(Analysis_Final_2[i, 38:40])
  if(sum(!is.na(test1[, 1]))>=2){
    P_value_RND1 <- t.test(CRISPRi_Ctrl_Round1, test1[which(!is.na(test1[, 1]))])
    Analysis_Final_2[i, 92] <- P_value_RND1$p.value
  } else {
    Analysis_Final_2[i, 92] <- NA
  }
  if(sum(!is.na(test2[, 1]))>=2){
    P_value_RND2 <- t.test(CRISPRi_Ctrl_Round2, test2[which(!is.na(test2[, 1]))])
    Analysis_Final_2[i, 93] <- P_value_RND2$p.value
  } else {
    Analysis_Final_2[i, 93] <- NA
  }
}
colnames(Analysis_Final_2)[92:93] <- c("P_value_RND1_M2", "P_value_RND2_M2")

Analysis_Final_2[which(!is.na(Analysis_Final_2$P_value_RND1_M2)), 94] <- p.adjust(Analysis_Final_2$P_value_RND1_M2[which(!is.na(Analysis_Final_2$P_value_RND1_M2))], 
                                                                                  method = "BH", 
                                                                                  n = length(Analysis_Final_2$P_value_RND1_M2[which(!is.na(Analysis_Final_2$P_value_RND1_M2))]))

Analysis_Final_2[which(!is.na(Analysis_Final_2$P_value_RND2_M2)), 95] <- p.adjust(Analysis_Final_2$P_value_RND2_M2[which(!is.na(Analysis_Final_2$P_value_RND2_M2))], 
                                                                                  method = "BH", 
                                                                                  n = length(Analysis_Final_2$P_value_RND2_M2[which(!is.na(Analysis_Final_2$P_value_RND2_M2))]))

colnames(Analysis_Final_2)[94:95] <- c("P.adjusted_RND1_M2", "P.adjusted_RND2_M2")

length(Analysis_Final_2$P_value_RND1_M2[which(Analysis_Final_2$P_value_RND1_M2<=0.05)])
length(Analysis_Final_2$P.adjusted_RND1_M2[which(Analysis_Final_2$P.adjusted_RND1_M2<=0.05)])
length(Analysis_Final_2$P_value_RND1_M2[which(Analysis_Final_2$P_value_RND1_M2<=0.1)])
length(Analysis_Final_2$P.adjusted_RND1_M2[which(Analysis_Final_2$P.adjusted_RND1_M2<=0.1)])

par(mfrow=c(2,2))
hist(Analysis_Final_2$P_value_RND1_M2,
     breaks = 20,
     xlab = "P-value", 
     ylab = "Frequency", 
     main = "P-values of all strains ROUND1", 
     col = "skyblue",
     ylim = c(0, 5000))
hist(Analysis_Final_2$P_value_RND1_M2[which(Analysis_Final_2$Control.gRNA==1)],
     breaks = 20,
     xlab = "P-value", 
     ylab = "Frequency", 
     main = "P-values of control strains ROUND1", 
     col = "skyblue",
     ylim = c(0, 10))
hist(Analysis_Final_2$P.adjusted_RND1_M2,
     breaks = 20,
     xlab = "P.adj", 
     ylab = "Frequency", 
     main = "P.adjusted values of all strains ROUND1", 
     col = "skyblue",
     ylim = c(0, 5000))
hist(Analysis_Final_2$P.adjusted_RND1_M2[which(Analysis_Final_2$Control.gRNA==1)],
     breaks = 20,
     xlab = "P.adj", 
     ylab = "Frequency", 
     main = "P.adjusted values of control strains ROUND1", 
     col = "skyblue",
     ylim = c(0, 10))

length(Analysis_Final_2$P_value_RND2_M2[which(Analysis_Final_2$P_value_RND2_M2<=0.05)])
length(Analysis_Final_2$P.adjusted_RND2_M2[which(Analysis_Final_2$P.adjusted_RND2_M2<=0.05)])
length(Analysis_Final_2$P_value_RND2_M2[which(Analysis_Final_2$P_value_RND2_M2<=0.1)])
length(Analysis_Final_2$P.adjusted_RND2_M2[which(Analysis_Final_2$P.adjusted_RND2_M2<=0.1)])

par(mfrow=c(2,2))
hist(Analysis_Final_2$P_value_RND2_M2,
     breaks = 20,
     xlab = "P-value", 
     ylab = "Frequency", 
     main = "P-values of all strains ROUND2", 
     col = "skyblue",
     ylim = c(0, 5000))
hist(Analysis_Final_2$P_value_RND2_M2[which(Analysis_Final_2$Control.gRNA==1)],
     breaks = 20,
     xlab = "P-value", 
     ylab = "Frequency", 
     main = "P-values of control strains ROUND2", 
     col = "skyblue",
     ylim = c(0, 10))
hist(Analysis_Final_2$P.adjusted_RND2_M2,
     breaks = 20,
     xlab = "P.adj", 
     ylab = "Frequency", 
     main = "P.adjusted values of all strains ROUND2", 
     col = "skyblue",
     ylim = c(0, 5000))
hist(Analysis_Final_2$P.adjusted_RND2_M2[which(Analysis_Final_2$Control.gRNA==1)],
     breaks = 20,
     xlab = "P.adj", 
     ylab = "Frequency", 
     main = "P.adjusted values of control strains ROUND2", 
     col = "skyblue",
     ylim = c(0, 10))

Analysis_Final_3 <- Analysis_Final_2

boxplot_stat_LPI_GT <- boxplot(Data_CRISPRi_aa$LPI_GT, cex = 0.3)
boxplot_stat_LPI_GT$stats 

Intermediate_50_M3 <- Data_CRISPRi_aa$LPI_GT[which(Data_CRISPRi_aa$LPI_GT >=-0.04373846
                                                   &Data_CRISPRi_aa$LPI_GT<=0.09804938)]
summary(Intermediate_50_M3)

Crispri_control_M4 <- Data_CRISPRi_aa$LPI_GT[which(Data_CRISPRi_aa$Control.gRNA==1)]
summary(Crispri_control_M4)

for(i in 1:nrow(Analysis_Final_3)){
  test1 <- t(Analysis_Final_3[i, 9:14])
  test2 <- t(Analysis_Final_3[i, 35:40])
  x1 <- sum(!is.na(test1[, 1]))
  x2 <- sum(!is.na(test2[, 1]))
  CTRL_GT_Mean_temp <- mean(test1[which(!is.na(test1[, 1]))])
  LPI_GT_Mean_temp <- mean(test2[which(!is.na(test2[, 1]))])
  Analysis_Final_3[i, 96] <- CTRL_GT_Mean_temp
  Analysis_Final_3[i, 97] <- x1
  Analysis_Final_3[i, 98] <- LPI_GT_Mean_temp
  Analysis_Final_3[i, 99] <- x2
}
colnames(Analysis_Final_3)[96:99] <- c("CTRL_GT_Mean_all", "n_CTRL", "LPI_GT_Mean_all", "n_LPI")

for(i in 1:nrow(Analysis_Final_3)){
  test <- t(Analysis_Final_3[i, 35:40])
  x <- sum(!is.na(test[, 1]))
  if(x>2){
    P.value_temp_M3 <- t.test(Intermediate_50_M3, test[which(!is.na(test[, 1]))])
    P.value_temp_M4 <- t.test(Crispri_control_M4, test[which(!is.na(test[, 1]))])
    Analysis_Final_3[i, 100] <- P.value_temp_M3$p.value
    Analysis_Final_3[i, 101] <- P.value_temp_M4$p.value
  } else {
    Analysis_Final_3[i, 100] <- NA
    Analysis_Final_3[i, 101] <- NA
  }
}
colnames(Analysis_Final_3)[100:101] <- c("P.value_M3", "P.value_M4")

Analysis_Final_3[which(!is.na(Analysis_Final_3$P.value_M3)), 102] <- p.adjust(Analysis_Final_3$P.value_M3[which(!is.na(Analysis_Final_3$P.value_M3))], 
                                                                              method = "BH", 
                                                                              n = length(Analysis_Final_3$P.value_M3[which(!is.na(Analysis_Final_3$P.value_M3))]))
Analysis_Final_3[which(!is.na(Analysis_Final_3$P.value_M4)), 103] <- p.adjust(Analysis_Final_3$P.value_M4[which(!is.na(Analysis_Final_3$P.value_M4))], 
                                                                              method = "BH", 
                                                                              n = length(Analysis_Final_3$P.value_M4[which(!is.na(Analysis_Final_3$P.value_M4))]))
colnames(Analysis_Final_3)[102:103] <- c("P.adjusted_M3", "P.adjusted_M4")

length(Analysis_Final_3$P.value_M3[which(Analysis_Final_3$P.value_M3<=0.05)])
length(Analysis_Final_3$P.adjusted_M3[which(Analysis_Final_3$P.adjusted_M3<=0.05)])
length(Analysis_Final_3$P.value_M3[which(Analysis_Final_3$P.value_M3<=0.1)])
length(Analysis_Final_3$P.adjusted_M3[which(Analysis_Final_3$P.adjusted_M3<=0.1)])

par(mfrow=c(2,2))
hist(Analysis_Final_3$P.value_M3,
     breaks = 20,
     xlab = "P-value", 
     ylab = "Frequency", 
     main = "P-value of all strains", 
     col = "skyblue",
     xlim = c(0, 1),
     ylim = c(0, 3000),
     cex.lab= 1.5)
hist(Analysis_Final_3$P.value_M3[which(Analysis_Final_3$Control.gRNA==1)],
     breaks = 20,
     xlab = "P-value", 
     ylab = "Frequency", 
     main = "P-value of control strains", 
     col = "skyblue",
     xlim = c(0, 1),
     ylim = c(0, 10),
     cex.lab= 1.5)
hist(Analysis_Final_3$P.adjusted_M3,
     breaks = 20,
     xlab = "P.value adjusted", 
     ylab = "Frequency", 
     main = "P.adjusted values of all strains", 
     col = "skyblue",
     xlim = c(0, 1),
     ylim = c(0, 3000),
     cex.lab= 1.5)
hist(Analysis_Final_3$P.adjusted_M3[which(Analysis_Final_3$Control.gRNA==1)],
     breaks = 20,
     xlab = "P.value adjusted", 
     ylab = "Frequency", 
     main = "P.adjusted values of control strains", 
     col = "skyblue",
     xlim = c(0, 1),
     ylim = c(0, 10),
     cex.lab= 1.5)

length(Analysis_Final_3$P.value_M4[which(Analysis_Final_3$P.value_M4<=0.05)])
length(Analysis_Final_3$P.adjusted_M4[which(Analysis_Final_3$P.adjusted_M4<=0.05)])
length(Analysis_Final_3$P.value_M4[which(Analysis_Final_3$P.value_M4<=0.1)])
length(Analysis_Final_3$P.adjusted_M4[which(Analysis_Final_3$P.adjusted_M4<=0.1)])

par(mfrow=c(2,2))
hist(Analysis_Final_3$P.value_M4,
     breaks = 20,
     xlab = "P-value", 
     ylab = "Frequency", 
     main = "P-value of all strains", 
     col = "skyblue",
     xlim = c(0, 1),
     ylim = c(0, 3000),
     cex.lab= 1.5)
hist(Analysis_Final_3$P.value_M4[which(Analysis_Final_3$Control.gRNA==1)],
     breaks = 20,
     xlab = "P-value", 
     ylab = "Frequency", 
     main = "P-value of control strains", 
     col = "skyblue",
     xlim = c(0, 1),
     ylim = c(0, 10),
     cex.lab= 1.5)
hist(Analysis_Final_3$P.adjusted_M4,
     breaks = 20,
     xlab = "P.value adjusted", 
     ylab = "Frequency", 
     main = "P.adjusted values of all strains", 
     col = "skyblue",
     xlim = c(0, 1),
     ylim = c(0, 3000),
     cex.lab= 1.5)
hist(Analysis_Final_3$P.adjusted_M4[which(Analysis_Final_3$Control.gRNA==1)],
     breaks = 20,
     xlab = "P.value adjusted", 
     ylab = "Frequency", 
     main = "P.adjusted values of control strains", 
     col = "skyblue",
     xlim = c(0, 1),
     ylim = c(0, 10),
     cex.lab= 1.5)

length(Analysis_Final_3$P.adjusted_M3[which(Analysis_Final_3$P.adjusted_M3 <= 0.1)])

max(Analysis_Final_3$LPI_GT_Mean_all[which(Analysis_Final_3$Control.gRNA==1)])

min(Analysis_Final_3$LPI_GT_Mean_all[which(Analysis_Final_3$Control.gRNA==1)])

candidate_padj_0.1_FIT_M3 <- which((Analysis_Final_3$LPI_GT_Mean_all < -0.03680838 & Analysis_Final_3$P.adjusted_M3<= 0.1))
length(candidate_padj_0.1_FIT_M3)

Fit_M3_complete <- Analysis_Final_3[candidate_padj_0.1_FIT_M3, ]
Fit_M3_complete <- Fit_M3_complete[order(Fit_M3_complete$LPI_GT_Mean_all, decreasing = FALSE), ]
str(Fit_M3_complete)

whole_Gene_list_Final <- read.csv("COMPILED_DATA/Gene_List_CRISPRi_lib.csv", na.strings = "", stringsAsFactors = FALSE)
rownames(whole_Gene_list_Final) <- whole_Gene_list_Final$LIB_ID

Fit_all_M3 <- data.frame(sort(table(Analysis_Final_3$GENE[candidate_padj_0.1_FIT_M3]), decreasing = TRUE))
y <- as.character(Fit_all_M3$Var1)
x <- whole_Gene_list_Final[y, ]
Fit_all_M3_description <- cbind(Fit_all_M3, x[, -1])
str(Fit_all_M3_description)
nrow(Fit_all_M3_description)

super_sen_M3 <- Analysis_Final_3[which(!is.na(Analysis_Final_3$CTRL_GT_Mean_all)
                                       &(Analysis_Final_3$n_LPI<3)
                                       &(
                                         is.na(Analysis_Final_3$LPI_GT_Mean_all)
                                         |(Analysis_Final_3$LPI_GT_Mean_all> 0.165662)
                                       )
), ]
nrow(super_sen_M3)

candidate_padj_0.1_SEN_M3 <- which((Analysis_Final_3$LPI_GT_Mean_all > 0.165662 & Analysis_Final_3$P.adjusted_M3<= 0.1))
length(candidate_padj_0.1_SEN_M3)

Sen_M3_complete <- rbind(super_sen_M3, Analysis_Final_3[candidate_padj_0.1_SEN_M3, ])
Sen_M3_complete <- Sen_M3_complete[order(Sen_M3_complete$LPI_GT_Mean_all, decreasing = TRUE), ]
nrow(Sen_M3_complete)

Sen_all_M3 <- data.frame(sort(table(c(Analysis_Final_3$GENE[candidate_padj_0.1_SEN_M3], super_sen_M3$GENE)), decreasing = TRUE))
y <- as.character(Sen_all_M3$Var1)
x <- whole_Gene_list_Final[y, ]
Sen_all_M3_description <- cbind(Sen_all_M3, x[, -1])
nrow(Sen_all_M3_description)

Fit_unique_M3 <- as.character(Fit_all_M3$Var1)
x <- whole_Gene_list_Final[Fit_unique_M3, ]
Fit_unique_M3_SGD_ID <- x$SGD_DB_ID
str(Fit_unique_M3_SGD_ID)

Sen_unique_M3 <- as.character(Sen_all_M3$Var1)
x <- whole_Gene_list_Final[Sen_unique_M3, ]
Sen_unique_M3_SGD_ID <- x$SGD_DB_ID
str(Sen_unique_M3_SGD_ID)

Growth_curve_data <- read.csv("COMPILED_DATA/Data_for_Representative_GC_SOM.csv", sep = "\t", header = TRUE)
str(Growth_curve_data)

Growth_curve_data_selected <- Growth_curve_data[, c("X", "X0_20_3", "X2_20_3", "X0_4_3", "X2_4_3", "X1_30_22", "X3_30_22")]
colnames(Growth_curve_data_selected) <- c("Time", "POL2-NRg-1_Basal", "POL2-NRg-1_Acetic", "RRP15-TRg-4_Basal", "RRP15-TRg-4_Acetic", "CC23_Basal", "CC23_Acetic")
str(Growth_curve_data_selected)

Growth_curve_data_selected[, 1] <- Growth_curve_data_selected[, 1]*20/60

library(reshape)
Growth_curve_data_selected_long <- reshape(data=Growth_curve_data_selected, idvar="Time",
                                           varying = colnames(Growth_curve_data_selected)[2:7],
                                           v.name=c("Population_size"),
                                           new.row.names = 1:30000,
                                           direction="long",
                                           timevar = "gRNA_condition",
                                           times = colnames(Growth_curve_data_selected)[2:7])
str(Growth_curve_data_selected_long)

library(ggplot2)
ggplot(Growth_curve_data_selected_long, aes(x=Time, y=Population_size, color=gRNA_condition)) +
  geom_smooth() +
  scale_y_continuous(trans='log10') + 
  theme_classic()

plot(Analysis_Final_3$LPI_GT_RND1_MEAN, Analysis_Final_3$LPI_GT_RND2_MEAN, 
     pch = 16, 
     cex = 0.5, 
     col = "black", 
     main = "Correlation between mean relative generation time (LPI GT) of Round 1 and 2", 
     xlab = "LPI GT Round1", 
     ylab = "LPI GT Round2", 
     xlim = c(-0.5, 2), 
     ylim = c(-0.5, 2),
     cex.lab=1.3,
     cex.axis=1.3)
points(Analysis_Final_3$LPI_GT_RND1_MEAN[candidate_padj_0.1_FIT_M3], 
       Analysis_Final_3$LPI_GT_RND2_MEAN[candidate_padj_0.1_FIT_M3],
       pch=16,
       cex = 0.7, 
       col = "blue")
points(Analysis_Final_3$LPI_GT_RND1_MEAN[candidate_padj_0.1_SEN_M3], 
       Analysis_Final_3$LPI_GT_RND2_MEAN[candidate_padj_0.1_SEN_M3], 
       pch=16,
       cex = 0.7, 
       col = "red")
points(Analysis_Final_3$LPI_GT_RND1_MEAN[which(Analysis_Final_3$Control.gRNA==1)], 
       Analysis_Final_3$LPI_GT_RND2_MEAN[which(Analysis_Final_3$Control.gRNA==1)], 
       pch=16,
       cex = 0.7, 
       col = "green")
stats_LPI_GT_Mean_RND1vsRND2_M3 <- lm(LPI_GT_RND2_MEAN ~ LPI_GT_RND1_MEAN, data = Analysis_Final_3)
stats_LPI_GT_Mean_RND1vsRND2_M3_selected <- lm(LPI_GT_RND2_MEAN[c(candidate_padj_0.1_SEN_M3, candidate_padj_0.1_FIT_M3)] ~ LPI_GT_RND1_MEAN[c(candidate_padj_0.1_SEN_M3, candidate_padj_0.1_FIT_M3)], data = Analysis_Final_3)

abline(stats_LPI_GT_Mean_RND1vsRND2_M3, lty=2, lwd=2)
abline(stats_LPI_GT_Mean_RND1vsRND2_M3_selected, col="red", lty=2, lwd=2)

summary(stats_LPI_GT_Mean_RND1vsRND2_M3)
cor(Analysis_Final_3$LPI_GT_RND1_MEAN, 
    Analysis_Final_3$LPI_GT_RND2_MEAN,  
    method = "pearson", 
    use = "complete.obs")

summary(stats_LPI_GT_Mean_RND1vsRND2_M3_selected)
cor(Analysis_Final_3$LPI_GT_RND1_MEAN[c(candidate_padj_0.1_SEN_M3, candidate_padj_0.1_FIT_M3)], 
    Analysis_Final_3$LPI_GT_RND2_MEAN[c(candidate_padj_0.1_SEN_M3, candidate_padj_0.1_FIT_M3)],  
    method = "pearson", 
    use = "complete.obs")

plot(Analysis_Final_3$CTRL_GT_Mean_all, Analysis_Final_3$LPI_GT_Mean_all, 
     pch = 16, 
     cex = 0.5, 
     col = "black", 
     main = "Selection of sensitive and tolerant strains", 
     xlab = "Normalized generation time (LSC GT) Basal.condition", 
     ylab = "Relative generation time (LPI GT) in 150mM Acetic acid", 
     xlim = c(-0.5, 2), 
     ylim = c(-0.5, 2),
     yaxt="n",
     xaxt="n",
     cex.lab=1.5)
points(Analysis_Final_3$CTRL_GT_Mean_all[candidate_padj_0.1_FIT_M3], 
       Analysis_Final_3$LPI_GT_Mean_all[candidate_padj_0.1_FIT_M3], 
       pch = 16, 
       cex = 0.5, 
       col = "blue")
points(Analysis_Final_3$CTRL_GT_Mean_all[candidate_padj_0.1_SEN_M3], 
       Analysis_Final_3$LPI_GT_Mean_all[candidate_padj_0.1_SEN_M3], 
       pch = 16, 
       cex = 0.5, 
       col = "red")
points(super_sen_M3$CTRL_GT_Mean_all, 
       super_sen_M3$LPI_GT_Mean_all, 
       pch = 16, 
       cex = 0.7, 
       col = "red")
points(Analysis_Final_3$CTRL_GT_Mean_all[which(Analysis_Final_3$Control.gRNA==1)], 
       Analysis_Final_3$LPI_GT_Mean_all[which(Analysis_Final_3$Control.gRNA==1)], 
       pch = 16, 
       cex = 0.6, 
       col = "green")
axis(side = 2, 
     at = c(-0.5, 0, 0.5, 1, 1.5, 2),
     cex.axis = 1.2,
     labels = c("-0.5", "0", "0.5", "1", "1.5", "2"), 
     tick = 0.05)
axis(side = 1, 
     at = c(-0.5, 0, 0.5, 1, 1.5, 2),
     cex.axis = 1.2,
     labels = c("-0.5", "0", "0.5", "1", "1.5", "2"), 
     tick = 0.05)
abline(h=c(-0.03680838, 0.165662), col="gray", lty=2, lwd=2)

Violin_LPI_Mean_M3 <- data.frame()
R <- length(which((Analysis_Final_3$Control.gRNA==0)
                  &(!is.na(Analysis_Final_3$LPI_GT_Mean_all))
))
Violin_LPI_Mean_M3[1:R, 1] <- Analysis_Final_3$LPI_GT_Mean_all[which((Analysis_Final_3$Control.gRNA==0)
                                                                     &(!is.na(Analysis_Final_3$LPI_GT_Mean_all)))]
Violin_LPI_Mean_M3[1:R, 2] <- "ALL"
R2 <- length(which(Analysis_Final_3$Control.gRNA==1))
Violin_LPI_Mean_M3[(R+1):(R+R2), 1] <- Analysis_Final_3$LPI_GT_Mean_all[which(Analysis_Final_3$Control.gRNA==1)]
Violin_LPI_Mean_M3[(R+1):(R+R2), 2] <- "CONTROL"
colnames(Violin_LPI_Mean_M3)[1:2] <- c("Mean", "Label")

library(ggplot2)
p_gg_M3 <- ggplot(Violin_LPI_Mean_M3, aes(x=Label, y=Mean, fill=Label)) + 
  geom_violin(trim=FALSE) + 
  geom_boxplot(width=0.1, fill="white") +
  labs(title="Violin plot",x="Data Type", y = "LPI_GT") +
  scale_fill_manual(values=c("white", "green"))
p_gg_M3 + theme_classic()

library("wordcloud")
wordcloud(words = Fit_all_M3$Var1, 
          freq = Fit_all_M3$Freq, 
          min.freq = 2,
          random.order=FALSE, rot.per=0.35, 
          colors=c("black", "red", "dark green", "blue"))

library("wordcloud")
wordcloud(words = Sen_all_M3$Var1, 
          freq = Sen_all_M3$Freq, 
          min.freq = 2,
          random.order=FALSE, rot.per=0.35, 
          colors=c("black", "red", "dark green", "blue"))

row.names(Analysis_Final_3) <- Analysis_Final_3$gRNA_name

Smith_Yepg_data <- read.csv("COMPILED_DATA/smith_YEPGdata.csv", na.strings = "")
str(Smith_Yepg_data)

for(i in 1:nrow(Analysis_Final_3)){
  x <- which(row.names(Analysis_Final_3)[i]==Smith_Yepg_data$guide_id)
  if(length(x)==0){
    Analysis_Final_3[i, 104:107] <- NA
  } else {
    Analysis_Final_3[i, 104:107] <- Smith_Yepg_data[x, 3:6]
  }
}
colnames(Analysis_Final_3)[104:107] <- colnames(Smith_Yepg_data)[3:6]
str((Analysis_Final_3)[104:107])

gRNA_Freq <- data.frame(sort(table(Analysis_Final_3$GENE), decreasing = TRUE))

par(mfrow=c(2,1))
hist(gRNA_Freq$Freq,
     breaks = 20,
     xlim = c(0, 20),
     ylim = c(0, 300),
     xlab = "Strains per gene",
     ylab = "Genes per bin",
     col = "skyblue",
     main = "Total number of strains per Gene")
hist(Analysis_Final_3$Midpoint_TSS_dist,
     breaks = 20,
     xlim=c(-300, 300),
     ylim=c(0, 1200),
     xlab = "Distance of gRNA relative to TSS",
     ylab = "Number of gRNA per bin",
     col = "skyblue",
     main= "Frequency of gRNA distance from TSS")

for(i in 1:nrow(Analysis_Final_3)){
  test1 <- t(Analysis_Final_3[i, 22:27])
  x1 <- sum(!is.na(test1[, 1]))
  AA_GT_Mean_temp <- mean(test1[which(!is.na(test1[, 1]))])
  Analysis_Final_3[i, 32] <- AA_GT_Mean_temp
}

par(mfrow=c(2,1))
hist(Analysis_Final_3$CTRL_GT_Mean_all,
     breaks = 25,
     xlim = c(-0.5, 2),
     ylim = c(0, 3000),
     xlab = "",
     ylab= "Strains per bin at Basal condition",
     col = "#A0A0A0",
     main = "  Normalized generation time 
     in basal and under acetic acid stress")
abline(v= c(log2(0.9), log2(1.1)), col ="black", lty=2)
hist(Analysis_Final_3$AA_GT_RND1_2_MEAN,
     breaks = 100,
     xlim = c(-0.5, 2),
     ylim = c(0, 3000),
     xlab = "Normalized generation time LSC GT",
     ylab= "Strains per bin at 150mM acetic acid",
     main = "",
     col = "#FF33FF")
abline(v= c(log2(0.9), log2(1.1)), col ="black", lty=2)

for(i in 1:nrow(whole_Gene_list_Final)){
  x <- as.character(unique(Smith_Yepg_data$ORF_Category[(Smith_Yepg_data$gene_name %in% whole_Gene_list_Final$LIB_ID[i])]))
  if(length(x)==0){
    whole_Gene_list_Final[i, 8] <- NA
  } else{
    whole_Gene_list_Final[i, 8] <- x
  }
}
colnames(whole_Gene_list_Final)[8] <- "ORF_Category"
whole_Gene_list_Final$ORF_Category <- as.factor(whole_Gene_list_Final$ORF_Category)
#Missing values were obtained from SGD
whole_Gene_list_Final$ORF_Category[which(is.na(whole_Gene_list_Final$ORF_Category))] <- c("Respiratory",
                                                                                          "Other", 
                                                                                          "Essential", 
                                                                                          "Essential", 
                                                                                          "Respiratory", 
                                                                                          "Other", 
                                                                                          "Respiratory", 
                                                                                          "Respiratory", 
                                                                                          "Other", 
                                                                                          "Respiratory", 
                                                                                          "Essential", 
                                                                                          "Essential", 
                                                                                          "Respiratory")
str(whole_Gene_list_Final)

Atc_liq_data <- read.csv("COMPILED_DATA/ATc_liq_titer_data.csv", na.strings = "NaN", header = TRUE)
str(Atc_liq_data)

Atc_liq_cc23 <- Atc_liq_data[which(Atc_liq_data$gRNA_name=="Ctrl-CC23"), ]

uniq_gRNA <- unique(Atc_liq_data$gRNA_name)
uniq_conc <- unique(Atc_liq_data$Atc_concentration)

Atc_liq_data[, 10:15] <- log(Atc_liq_data[, 4:9])

for(i in 1:length(uniq_gRNA)){
  for(j in 1:length(uniq_conc)){
    Atc_liq_data[which(Atc_liq_data$gRNA_name==uniq_gRNA[i]&
                         Atc_liq_data$Atc_concentration==uniq_conc[j]), 16:21] <- 
      Atc_liq_data[which(Atc_liq_data$gRNA_name==uniq_gRNA[i]&
                           Atc_liq_data$Atc_concentration==uniq_conc[j]), 10:15] - 
      Atc_liq_data[which(Atc_liq_data$gRNA_name=="Ctrl-CC23"&
                           Atc_liq_data$Atc_concentration==uniq_conc[j]), 10:15]
  }
}

for(i in 1:nrow(Atc_liq_data)){
  Atc_liq_data[i, 22] <- mean(as.numeric(Atc_liq_data[i, c(16, 19)][which(!is.na(Atc_liq_data[i, c(16, 19)]))]))
  Atc_liq_data[i, 23] <- sd(as.numeric(Atc_liq_data[i, c(16, 19)][which(!is.na(Atc_liq_data[i, c(16, 19)]))]))
  Atc_liq_data[i, 24] <- mean(as.numeric(Atc_liq_data[i, c(17, 20)][which(!is.na(Atc_liq_data[i, c(17, 20)]))]))
  Atc_liq_data[i, 25] <- sd(as.numeric(Atc_liq_data[i, c(17, 20)][which(!is.na(Atc_liq_data[i, c(17, 20)]))]))
  Atc_liq_data[i, 26] <- mean(as.numeric(Atc_liq_data[i, c(18, 21)][which(!is.na(Atc_liq_data[i, c(18, 21)]))]))
  Atc_liq_data[i, 27] <- sd(as.numeric(Atc_liq_data[i, c(18, 21)][which(!is.na(Atc_liq_data[i, c(18, 21)]))]))
}

colnames(Atc_liq_data)[10:15] <- paste0("log_", colnames(Atc_liq_data)[4:9])
colnames(Atc_liq_data)[16:21] <- paste0("LSC_", colnames(Atc_liq_data)[4:9])
colnames(Atc_liq_data)[22:23] <- paste0(c("Mean_", "SD_"), "LSC_Lag")
colnames(Atc_liq_data)[24:25] <- paste0(c("Mean_", "SD_"), "LSC_GT")
colnames(Atc_liq_data)[26:27] <- paste0(c("Mean_", "SD_"), "LSC_Yield")

name_gRNA_atc_titer <- c("ACT1-NRg-5", "ACT1-NRg-8", "Ctrl-CC11", "SEC21-NRg-5", "VPS1-TRg-1")
test <- data.frame()
Atc_titer_subset <- data.frame()
for(i in 1:length(name_gRNA_atc_titer)){
  test <- Atc_liq_data[which(Atc_liq_data$gRNA_name==name_gRNA_atc_titer[i]), ]
  Atc_titer_subset <- rbind(Atc_titer_subset, test)
}

library(ggplot2)
plt1 <- ggplot(Atc_titer_subset, aes(x=Atc_concentration, y=Mean_LSC_Lag, group=gRNA_name, color=gRNA_name)) + 
  geom_pointrange(aes(ymin=Mean_LSC_Lag-SD_LSC_Lag, ymax=Mean_LSC_Lag+SD_LSC_Lag)) +
  labs(title="Normalized Lag phase at different ATc concentration", x="ATc (ug/ml)", y = "LSC Lag phase")+
  theme_classic()+
  scale_color_manual(values=c('#999999','#E69F00', "green3", "red", "black"))+
  scale_x_continuous(breaks = c(0, 1, 2, 3, 5, 7, 10, 15, 25),
                     labels = c("0", "1", "2", "3", "5", "7", "10", "15", "25"),
                     limits = c(0, 26))+
  scale_y_continuous(breaks = c(-1, -0.5, 0, 0.5, 1, 1.5, 2),
                     labels = c("-1", "-0.5", "0", "0.5", "1", "1.5", "2"),
                     limits = c(-1, 2))+
  theme(legend.position="none")
suppressWarnings(print(plt1))

plt2 <- ggplot(Atc_titer_subset, aes(x=Atc_concentration, y=Mean_LSC_GT, group=gRNA_name, color=gRNA_name)) + 
  geom_pointrange(aes(ymin=Mean_LSC_GT-SD_LSC_GT, ymax=Mean_LSC_GT+SD_LSC_GT)) +
  labs(title="Normalized Generation time at different ATc concentration", x="ATc (ug/ml)", y = "LSC GT")+
  theme_classic()+
  scale_color_manual(values=c('#999999','#E69F00', "green3", "red", "black"))+
  scale_x_continuous(breaks = c(0, 1, 2, 3, 5, 7, 10, 15, 25),
                     labels = c("0", "1", "2", "3", "5", "7", "10", "15", "25"),
                     limits = c(0, 26))+
  scale_y_continuous(breaks = c(-1, -0.5, 0, 0.5, 1, 1.5, 2),
                     labels = c("-1", "-0.5", "0", "0.5", "1", "1.5", "2"),
                     limits = c(-1, 2))+
  theme(legend.position="none")
suppressWarnings(print(plt2))

plt3 <- ggplot(Atc_titer_subset, aes(x=Atc_concentration, y=Mean_LSC_Yield, group=gRNA_name, color=gRNA_name)) + 
  geom_pointrange(aes(ymin=Mean_LSC_Yield-SD_LSC_Yield, ymax=Mean_LSC_Yield+SD_LSC_Yield)) +
  labs(title="Normalized Yield at different ATc concentration", x="ATc (ug/ml)", y = "LSC Yield")+
  theme_classic()+
  scale_color_manual(values=c('#999999','#E69F00', "green3", "red", "black"))+
  scale_x_continuous(breaks = c(0, 1, 2, 3, 5, 7, 10, 15, 25),
                     labels = c("0", "1", "2", "3", "5", "7", "10", "15", "25"),
                     limits = c(0, 26))+
  scale_y_continuous(breaks = c(-2.5, -2, -1.5, -1, -0.5, 0, 0.5),
                     labels = c("-2.5", "-2", "-1.5", "-1", "-0.5", "0", "0.5"),
                     limits = c(-2.5, 0.5))+
  theme(legend.position="bottom")
suppressWarnings(print(plt3))

select_genes <- read.table("COMPILED_DATA/selected_genes.txt", header = FALSE, sep = "\t", as.is = TRUE)

y <- vector(mode = "numeric", length = 0)
for(i in 1:length(select_genes$V1)){
  x <- which(Analysis_Final_3$GENE==select_genes$V1[i])
  y <- c(y, x)
}
select_strains <- Analysis_Final_3$gRNA_name[y]

bot_top_50 <- read.csv("COMPILED_DATA/bottom_top_50.csv", stringsAsFactors = FALSE, header = TRUE)

bot_top_50 <- bot_top_50[-which(bot_top_50$n_CTRL<4 & bot_top_50$CTRL_GT_Mean_all > 0.1), ]
bot_top_50_strains <- bot_top_50$gRNA_name

validation_strains <- as.character(union(bot_top_50_strains, select_strains))

Val_data <- read.csv("COMPILED_DATA/Validation_Bioscreen_data.csv", na.strings = "NaN", header = TRUE)
str(Val_data)

Val_data[, 41:76] <- log(Val_data[, 5:40])
colnames(Val_data)[41:76] <- paste0("log_", colnames(Val_data)[5:40])

Val_data_cc <- Val_data[which(Val_data$gRNA_name=="CC23"), ]

for(j in 1:2){
  #Ctrl_Lag
  Val_data_cc[j, 77] <- mean(as.numeric(Val_data_cc[j, c(41, 47, 53, 59, 65, 71)][which(!is.na(Val_data_cc[j, c(41, 47, 53, 59, 65, 71)]))]))
  #Ctrl_GT
  Val_data_cc[j, 78] <- mean(as.numeric(Val_data_cc[j, c(42, 48, 54, 60, 66, 72)][which(!is.na(Val_data_cc[j, c(42, 48, 54, 60, 66, 72)]))]))
  #Ctrl_Yield
  Val_data_cc[j, 79] <- mean(as.numeric(Val_data_cc[j, c(43, 49, 55, 61, 67, 73)][which(!is.na(Val_data_cc[j, c(43, 49, 55, 61, 67, 73)]))]))
  #AA150_Lag
  Val_data_cc[j, 80] <- mean(as.numeric(Val_data_cc[j, c(44, 50, 56)][which(!is.na(Val_data_cc[j, c(44, 50, 56)]))]))
  #AA150_GT
  Val_data_cc[j, 81] <- mean(as.numeric(Val_data_cc[j, c(45, 51, 57)][which(!is.na(Val_data_cc[j, c(45, 51, 57)]))]))
  #AA150_Yield
  Val_data_cc[j, 82] <- mean(as.numeric(Val_data_cc[j, c(46, 52, 58)][which(!is.na(Val_data_cc[j, c(46, 52, 58)]))]))
  #AA125_Lag
  Val_data_cc[j, 83] <- mean(as.numeric(Val_data_cc[j, c(62, 68, 74)][which(!is.na(Val_data_cc[j, c(62, 68, 74)]))]))
  #AA125_GT
  Val_data_cc[j, 84] <- mean(as.numeric(Val_data_cc[j, c(63, 69, 75)][which(!is.na(Val_data_cc[j, c(63, 69, 75)]))]))
  #AA125_GT
  Val_data_cc[j, 85] <- mean(as.numeric(Val_data_cc[j, c(64, 70, 76)][which(!is.na(Val_data_cc[j, c(64, 70, 76)]))]))
}
colnames(Val_data_cc)[77:79] <- paste0("Mean_Ctrl_", c("Lag", "GT", "Yield"))
colnames(Val_data_cc)[80:82] <- paste0("Mean_AA150_", c("Lag", "GT", "Yield"))
colnames(Val_data_cc)[83:85] <- paste0("Mean_AA125_", c("Lag", "GT", "Yield"))

#Replicate1
for(i in 1:100){
  Val_data[i, 77:82] <- Val_data[i, 41:46]-Val_data_cc["10", 77:82]
}
#Replicate2
for(i in 1:100){
  Val_data[i, 83:88] <- Val_data[i, 47:52]-Val_data_cc["10", 77:82]
}
#Replicate3
for(i in 1:100){
  Val_data[i, 89:94] <- Val_data[i, 53:58]-Val_data_cc["10", 77:82]
}

#Replicate1
for(i in 101:200){
  Val_data[i, 77:82] <- Val_data[i, 41:46]-Val_data_cc["110", 77:82]
}
#Replicate2
for(i in 101:200){
  Val_data[i, 83:88] <- Val_data[i, 47:52]-Val_data_cc["110", 77:82]
}
#Replicate3
for(i in 101:200){
  Val_data[i, 89:94] <- Val_data[i, 53:58]-Val_data_cc["110", 77:82]
}

#Replicate1
for(i in 1:100){
  Val_data[i, 95:100] <- Val_data[i, 59:64]-Val_data_cc["10", c(77:79, 83:85)]
}
#Replicate2
for(i in 1:100){
  Val_data[i, 101:106] <- Val_data[i, 65:70]-Val_data_cc["10", c(77:79, 83:85)]
}
#Replicate3
for(i in 1:100){
  Val_data[i, 107:112] <- Val_data[i, 71:76]-Val_data_cc["10", c(77:79, 83:85)]
}

#Replicate1
for(i in 101:200){
  Val_data[i, 95:100] <- Val_data[i, 59:64]-Val_data_cc["110", c(77:79, 83:85)]
}
#Replicate2
for(i in 101:200){
  Val_data[i, 101:106] <- Val_data[i, 65:70]-Val_data_cc["110", c(77:79, 83:85)]
}
#Replicate3
for(i in 101:200){
  Val_data[i, 107:112] <- Val_data[i, 71:76]-Val_data_cc["110", c(77:79, 83:85)]
}

colnames(Val_data)[77:112] <- paste0("LSC_", colnames(Val_data)[5:40])

#Replicate1_LPI_AA150
Val_data[, 113:115] <- Val_data[, 80:82]-Val_data[, 77:79]
colnames(Val_data)[113:115] <- paste0("LPI_", colnames(Val_data)[8:10])
#Replicate2_LPI_AA150
Val_data[, 116:118] <- Val_data[, 86:88]-Val_data[, 83:85]
colnames(Val_data)[116:118] <- paste0("LPI_", colnames(Val_data)[14:16])
#Replicate3_LPI_AA150
Val_data[, 119:121] <- Val_data[, 92:94]-Val_data[, 89:91]
colnames(Val_data)[119:121] <- paste0("LPI_", colnames(Val_data)[20:22])
#Replicate1_LPI_AA125
Val_data[, 122:124] <- Val_data[, 98:100]-Val_data[, 95:97]
colnames(Val_data)[122:124] <- paste0("LPI_", colnames(Val_data)[26:28])
#Replicate2_LPI_AA125
Val_data[, 125:127] <- Val_data[, 104:106]-Val_data[, 101:103]
colnames(Val_data)[125:127] <- paste0("LPI_", colnames(Val_data)[32:34])
#Replicate3_LPI_AA125
Val_data[, 128:130] <- Val_data[, 110:112]-Val_data[, 107:109]
colnames(Val_data)[128:130] <- paste0("LPI_", colnames(Val_data)[38:40])

for (i in 1:nrow(Val_data)){
  #LPI_AA150_Lag
  Val_data[i, 131] <- mean(as.numeric(Val_data[i, c(113, 116, 119)][which(!is.na(Val_data[i, c(113, 116, 119)]))]))
  Val_data[i, 132] <- sd(as.numeric(Val_data[i, c(113, 116, 119)][which(!is.na(Val_data[i, c(113, 116, 119)]))]))
  Val_data[i, 133] <- length(as.numeric(Val_data[i, c(113, 116, 119)][which(!is.na(Val_data[i, c(113, 116, 119)]))]))
  #LPI_AA150_GT
  Val_data[i, 134] <- mean(as.numeric(Val_data[i, c(114, 117, 120)][which(!is.na(Val_data[i, c(114, 117, 120)]))]))
  Val_data[i, 135] <- sd(as.numeric(Val_data[i, c(114, 117, 120)][which(!is.na(Val_data[i, c(114, 117, 120)]))]))
  Val_data[i, 136] <- length(as.numeric(Val_data[i, c(114, 117, 120)][which(!is.na(Val_data[i, c(114, 117, 120)]))]))
  #LPI_AA150_Yield
  Val_data[i, 137] <- mean(as.numeric(Val_data[i, c(115, 118, 121)][which(!is.na(Val_data[i, c(115, 118, 121)]))]))
  Val_data[i, 138] <- sd(as.numeric(Val_data[i, c(115, 118, 121)][which(!is.na(Val_data[i, c(115, 118, 121)]))]))
  Val_data[i, 139] <- length(as.numeric(Val_data[i, c(115, 118, 121)][which(!is.na(Val_data[i, c(115, 118, 121)]))]))
  #LPI_AA125_Lag
  Val_data[i, 140] <- mean(as.numeric(Val_data[i, c(122, 125, 128)][which(!is.na(Val_data[i, c(122, 125, 128)]))]))
  Val_data[i, 141] <- sd(as.numeric(Val_data[i, c(122, 125, 128)][which(!is.na(Val_data[i, c(122, 125, 128)]))]))
  Val_data[i, 142] <- length(as.numeric(Val_data[i, c(122, 125, 128)][which(!is.na(Val_data[i, c(122, 125, 128)]))]))
  #LPI_AA125_GT
  Val_data[i, 143] <- mean(as.numeric(Val_data[i, c(123, 126, 129)][which(!is.na(Val_data[i, c(123, 126, 129)]))]))
  Val_data[i, 144] <- sd(as.numeric(Val_data[i, c(123, 126, 129)][which(!is.na(Val_data[i, c(123, 126, 129)]))]))
  Val_data[i, 145] <- length(as.numeric(Val_data[i, c(123, 126, 129)][which(!is.na(Val_data[i, c(123, 126, 129)]))]))
  #LPI_AA125_GT
  Val_data[i, 146] <- mean(as.numeric(Val_data[i, c(124, 127, 130)][which(!is.na(Val_data[i, c(124, 127, 130)]))]))
  Val_data[i, 147] <- sd(as.numeric(Val_data[i, c(124, 127, 130)][which(!is.na(Val_data[i, c(124, 127, 130)]))]))
  Val_data[i, 148] <- length(as.numeric(Val_data[i, c(124, 127, 130)][which(!is.na(Val_data[i, c(124, 127, 130)]))]))
}
#Assigning the column names
colnames(Val_data)[131:133] <- paste0(c("Mean_", "SD_", "N_"), "LPI_AA150_Lag")
colnames(Val_data)[134:136] <- paste0(c("Mean_", "SD_", "N_"), "LPI_AA150_GT")
colnames(Val_data)[137:139] <- paste0(c("Mean_", "SD_", "N_"), "LPI_AA150_Yield")
colnames(Val_data)[140:142] <- paste0(c("Mean_", "SD_", "N_"), "LPI_AA125_Lag")
colnames(Val_data)[143:145] <- paste0(c("Mean_", "SD_", "N_"), "LPI_AA125_GT")
colnames(Val_data)[146:148] <- paste0(c("Mean_", "SD_", "N_"), "LPI_AA125_Yield")

dCtrl_strains <- c("CC23", "CC14", "CC2", "CC28", "CC30", "CC32", "CC34")
dCTRL_rows <- vector(mode = "integer", length = 0)
test <- data.frame()
Val_data_dCTRL <- data.frame()
for(i in 1:length(dCtrl_strains)){
  test <- Val_data[which(Val_data$gRNA_name==dCtrl_strains[i]), ]
  dCTRL_rows <- c(dCTRL_rows, which(Val_data$gRNA_name==dCtrl_strains[i]))
  Val_data_dCTRL <- rbind(Val_data_dCTRL, test)
}

m1 <- vector(mode = "numeric", length = 0)
m2 <- vector(mode = "numeric", length = 0)
test1 <- data.frame()
Val_data_dCTRL_F <- data.frame()
Val_data_dCTRL_F[1:7, 1:3] <- Val_data_dCTRL[c("10", "19", "28", "46", "55", "64", "73"), 2:4]
for (i in 1:length(dCtrl_strains)){
  test1 <- Val_data_dCTRL[which(Val_data_dCTRL$gRNA_name==dCtrl_strains[i]), ]
  #LPI_AA150_Lag
  m1 <- as.numeric(test1[1, c(113, 116, 119)][which(!is.na(test1[1, c(113, 116, 119)]))])
  m2 <- as.numeric(test1[2, c(113, 116, 119)][which(!is.na(test1[2, c(113, 116, 119)]))])
  Val_data_dCTRL_F[i, 4] <- mean(c(m1, m2))
  Val_data_dCTRL_F[i, 5] <- sd(c(m1, m2))
  Val_data_dCTRL_F[i, 6] <- length(c(m1, m2))
  #LPI_AA150_GT
  m1 <- as.numeric(test1[1, c(114, 117, 120)][which(!is.na(test1[1, c(114, 117, 120)]))])
  m2 <- as.numeric(test1[2, c(114, 117, 120)][which(!is.na(test1[2, c(114, 117, 120)]))])
  Val_data_dCTRL_F[i, 7] <- mean(c(m1, m2))
  Val_data_dCTRL_F[i, 8] <- sd(c(m1, m2))
  Val_data_dCTRL_F[i, 9] <- length(c(m1, m2))
  #LPI_AA150_Yield
  m1 <- as.numeric(test1[1, c(115, 118, 121)][which(!is.na(test1[1, c(115, 118, 121)]))])
  m2 <- as.numeric(test1[2, c(115, 118, 121)][which(!is.na(test1[2, c(115, 118, 121)]))])
  Val_data_dCTRL_F[i, 10] <- mean(c(m1, m2))
  Val_data_dCTRL_F[i, 11] <- sd(c(m1, m2))
  Val_data_dCTRL_F[i, 12] <- length(c(m1, m2))
  #LPI_AA125_Lag
  m1 <- as.numeric(test1[1, c(122, 125, 128)][which(!is.na(test1[1, c(122, 125, 128)]))])
  m2 <- as.numeric(test1[2, c(122, 125, 128)][which(!is.na(test1[2, c(122, 125, 128)]))])
  Val_data_dCTRL_F[i, 13] <- mean(c(m1, m2))
  Val_data_dCTRL_F[i, 14] <- sd(c(m1, m2))
  Val_data_dCTRL_F[i, 15] <- length(c(m1, m2))
  #LPI_AA125_GT
  m1 <- as.numeric(test1[1, c(123, 126, 129)][which(!is.na(test1[1, c(123, 126, 129)]))])
  m2 <- as.numeric(test1[2, c(123, 126, 129)][which(!is.na(test1[2, c(123, 126, 129)]))])
  Val_data_dCTRL_F[i, 16] <- mean(c(m1, m2))
  Val_data_dCTRL_F[i, 17] <- sd(c(m1, m2))
  Val_data_dCTRL_F[i, 18] <- length(c(m1, m2))
  #LPI_AA125_GT
  m1 <- as.numeric(test1[1, c(124, 127, 130)][which(!is.na(test1[1, c(124, 127, 130)]))])
  m2 <- as.numeric(test1[2, c(124, 127, 130)][which(!is.na(test1[2, c(124, 127, 130)]))])
  Val_data_dCTRL_F[i, 19] <- mean(c(m1, m2))
  Val_data_dCTRL_F[i, 20] <- sd(c(m1, m2))
  Val_data_dCTRL_F[i, 21] <- length(c(m1, m2))
}
colnames(Val_data_dCTRL_F)[4:6] <- paste0(c("Mean_", "SD_", "N_"), "LPI_AA150_Lag")
colnames(Val_data_dCTRL_F)[7:9] <- paste0(c("Mean_", "SD_", "N_"), "LPI_AA150_GT")
colnames(Val_data_dCTRL_F)[10:12] <- paste0(c("Mean_", "SD_", "N_"), "LPI_AA150_Yield")
colnames(Val_data_dCTRL_F)[13:15] <- paste0(c("Mean_", "SD_", "N_"), "LPI_AA125_Lag")
colnames(Val_data_dCTRL_F)[16:18] <- paste0(c("Mean_", "SD_", "N_"), "LPI_AA125_GT")
colnames(Val_data_dCTRL_F)[19:21] <- paste0(c("Mean_", "SD_", "N_"), "LPI_AA125_Yield")

Val_data_curated <- Val_data[-dCTRL_rows, ]

Val_data_curated <- Val_data_curated[-which(Val_data_curated$gRNA_name=="BLANK"), ]

Val_data_column_trimmed <- Val_data_curated[, c(2:4, 131:148)]

Validation_LPI_all <- rbind(Val_data_column_trimmed, Val_data_dCTRL_F)
rownames(Validation_LPI_all) <- Validation_LPI_all$gRNA_name
str(Validation_LPI_all)

Val_whole_data_dCTRL <- Val_data[dCTRL_rows, ]

Val_dCTRL_lag_125 <- c(as.numeric(Val_whole_data_dCTRL[, 122]), as.numeric(Val_whole_data_dCTRL[, 125]), as.numeric(Val_whole_data_dCTRL[, 128]))
Val_dCTRL_GT_125 <- c(as.numeric(Val_whole_data_dCTRL[, 123]), as.numeric(Val_whole_data_dCTRL[, 126]), as.numeric(Val_whole_data_dCTRL[, 129]))
Val_dCTRL_Yield_125 <- c(as.numeric(Val_whole_data_dCTRL[, 124]), as.numeric(Val_whole_data_dCTRL[, 127]), as.numeric(Val_whole_data_dCTRL[, 130]))

for(i in 1:nrow(Val_data_curated)){
  test_lag <- t(Val_data_curated[i, c(122, 125, 128)])
  test_GT <- t(Val_data_curated[i, c(123, 126, 129)])
  test_Yield <- t(Val_data_curated[i, c(124, 127, 130)])
  x1 <- sum(!is.na(test_lag[, 1]))
  x2 <- sum(!is.na(test_GT[, 1]))
  x3 <- sum(!is.na(test_Yield[, 1]))
  if(x1>1){
    P.value_lag_125<- t.test(Val_dCTRL_lag_125, test_lag[which(!is.na(test_lag[, 1]))])
    Val_data_curated[i, 149] <- P.value_lag_125$p.value
  } else {
    Val_data_curated[i, 149] <- NA
  }
  if(x2>1){
    P.value_GT_125<- t.test(Val_dCTRL_GT_125, test_GT[which(!is.na(test_GT[, 1]))])
    Val_data_curated[i, 150] <- P.value_GT_125$p.value
  } else {
    Val_data_curated[i, 150] <- NA
  }
  if(x3>1){
    P.value_Yield_125<- t.test(Val_dCTRL_Yield_125, test_Yield[which(!is.na(test_Yield[, 1]))])
    Val_data_curated[i, 151] <- P.value_Yield_125$p.value
  } else {
    Val_data_curated[i, 151] <- NA
  }
}
colnames(Val_data_curated)[149:151] <- c("P.value_lag_125", "P.value_GT_125", "P.value_Yield_125")

Val_data_curated[which(!is.na(Val_data_curated$P.value_lag_125)), 152] <- p.adjust(Val_data_curated$P.value_lag_125[which(!is.na(Val_data_curated$P.value_lag_125))], 
                                                                                   method = "BH", 
                                                                                   n = length(Val_data_curated$P.value_lag_125[which(!is.na(Val_data_curated$P.value_lag_125))]))
Val_data_curated[which(!is.na(Val_data_curated$P.value_GT_125)), 153] <- p.adjust(Val_data_curated$P.value_GT_125[which(!is.na(Val_data_curated$P.value_GT_125))], 
                                                                                  method = "BH", 
                                                                                  n = length(Val_data_curated$P.value_GT_125[which(!is.na(Val_data_curated$P.value_GT_125))]))
Val_data_curated[which(!is.na(Val_data_curated$P.value_Yield_125)), 154] <- p.adjust(Val_data_curated$P.value_Yield_125[which(!is.na(Val_data_curated$P.value_Yield_125))], 
                                                                                     method = "BH", 
                                                                                     n = length(Val_data_curated$P.value_Yield_125[which(!is.na(Val_data_curated$P.value_Yield_125))]))
colnames(Val_data_curated)[152:154] <- c("P.adj_lag_125", "P.adj_GT_125", "P.adj_Yield_125")
rownames(Val_data_curated) <- Val_data_curated$gRNA_name
str(Val_data_curated)

validation_strains_data <- Analysis_Final_3[validation_strains, ]

dCTRL_data_scan_o_matic <- Analysis_Final_3[(Analysis_Final_3$gRNA_name %in% dCtrl_strains), ]

validation_strains_data <- rbind(validation_strains_data, dCTRL_data_scan_o_matic)

#setting the row order similar to that of the bioscreen dataset Validation_LPI_all
validation_strains_data <- validation_strains_data[rownames(Validation_LPI_all), ]

Validation_new_df <- validation_strains_data[, c(1:8, 96:97, 87, 35:40, 98:99, 89, 100, 102, 104:107, 58, 60, 74:79, 84, 86)]
Validation_new_df <- cbind(Validation_new_df, Validation_LPI_all[, 4:21])

bot50 <- bot_top_50[1:48, ]
Top50 <- bot_top_50[49:98, ]

#dCTRL strains (1)
Validation_new_df[(Validation_new_df$gRNA_name %in% dCtrl_strains), 55] <- 1
#top 50 acetic acid tolerant strains (2),
Validation_new_df[(Validation_new_df$gRNA_name %in% Top50$gRNA_name), 55] <- 2
#most 50 acetic acid sensitive strains (3),
Validation_new_df[(Validation_new_df$gRNA_name %in% bot50$gRNA_name), 55] <- 3
#Other candidates (4)
Validation_new_df[which(is.na(Validation_new_df$V55)), 55] <- 4
#Changing the column name
colnames(Validation_new_df)[55] <- "Strain_category"
str(Validation_new_df)

plot(Validation_new_df$LPI_GT_Mean_all, Validation_new_df$Mean_LPI_AA125_GT,
     xlim = c(-0.5, 2),
     ylim = c(-0.5, 2),
     pch = 16, 
     cex = 0.7, 
     col = "black",
     xlab = "LPI_GT_Scan-O-Matic at 150mM acetic acid",
     ylab = "LPI_GT_Bioscreen at 125mM acetic acid")
points(Validation_new_df$LPI_GT_Mean_all[which(Validation_new_df$Strain_category==1)], 
       Validation_new_df$Mean_LPI_AA125_GT[which(Validation_new_df$Strain_category==1)], 
       pch = 16, 
       cex = 0.8, 
       col = "green")
points(Validation_new_df$LPI_GT_Mean_all[which(Validation_new_df$Strain_category==2)], 
       Validation_new_df$Mean_LPI_AA125_GT[which(Validation_new_df$Strain_category==2)], 
       pch = 16, 
       cex = 0.8, 
       col = "blue")
points(Validation_new_df$LPI_GT_Mean_all[which(Validation_new_df$Strain_category==3)], 
       Validation_new_df$Mean_LPI_AA125_GT[which(Validation_new_df$Strain_category==3)], 
       pch = 16, 
       cex = 0.8, 
       col = "red")
abline(lm(Validation_new_df$Mean_LPI_AA125_GT ~ Validation_new_df$LPI_GT_Mean_all))
text(Validation_new_df$LPI_GT_Mean_all, 
     Validation_new_df$Mean_LPI_AA125_GT,
     labels=Validation_new_df$GENE, 
     cex= 0.1, 
     pos = 2)
stat_GT_125 <- lm(Validation_new_df$Mean_LPI_AA125_GT ~ Validation_new_df$LPI_GT_Mean_all)
summary(stat_GT_125)

print("At 125mM Acetic Acid")
print("lag vs GT")
cor(Validation_new_df$Mean_LPI_AA125_Lag, 
    Validation_new_df$Mean_LPI_AA125_GT,  
    method = "pearson", 
    use = "complete.obs")
print("lag vs Yield")
cor(Validation_new_df$Mean_LPI_AA125_Lag, 
    Validation_new_df$Mean_LPI_AA125_Yield,  
    method = "pearson", 
    use = "complete.obs")
print("GT vs Yield")
cor(Validation_new_df$Mean_LPI_AA125_GT, 
    Validation_new_df$Mean_LPI_AA125_Yield,  
    method = "pearson", 
    use = "complete.obs")
print("At 150mM Acetic Acid")  
print("lag vs GT")
cor(Validation_new_df$Mean_LPI_AA150_Lag, 
    Validation_new_df$Mean_LPI_AA150_GT,  
    method = "pearson", 
    use = "complete.obs")
print("lag vs Yield")
cor(Validation_new_df$Mean_LPI_AA150_Lag, 
    Validation_new_df$Mean_LPI_AA150_Yield,  
    method = "pearson", 
    use = "complete.obs")
print("GT vs Yield")
cor(Validation_new_df$Mean_LPI_AA150_GT, 
    Validation_new_df$Mean_LPI_AA150_Yield,  
    method = "pearson", 
    use = "complete.obs")

#Make a color palette
colfunc5<-colorRampPalette(c("goldenrod4", "goldenrod", "white", "turquoise", "turquoise4"))
plot(rep(1,100), col=colfunc5(100), pch=19,cex=2)

brk1 <- c(seq(-0.5, -0.05, length.out = 48))

brk2 <- c(seq(-0.04, 0.04, length.out = 5))

brk3 <- c(seq(0.05, 2, length.out = 48))

brk_F <- c(brk1, brk2, brk3)

Validation_new_df <- Validation_new_df[order(Validation_new_df$LPI_GT_Mean_all, decreasing = TRUE), ]

Validation_new_df <- Validation_new_df[c(170:183, 1:169), ]

Validation_new_df[, 56] <- Validation_new_df[, 43]*(-1)
Validation_new_df[, 57] <- Validation_new_df[, 52]*(-1)
colnames(Validation_new_df)[56:57] <- c("Mean_LPI_AA150_Yield(-1)", "Mean_LPI_AA125_Yield(-1)")

library(pheatmap)
pheatmap(as.matrix(Validation_new_df[, c(9, 18, 46, 49, 57, 37, 40, 56)]),
         color = colfunc5(100),
         breaks = brk_F,
         border_color = "white",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         cellwidth = 10,
         cellheight = 10)

Validation_new_df_Grp <- read.csv("COMPILED_DATA/Validation_new_df_with_Groups_by_GO.csv", stringsAsFactors = FALSE, na.strings = )
rownames(Validation_new_df_Grp) <- Validation_new_df_Grp$gRNA_name

library(factoextra)
#Plotting PCA for only Proteasomal genes and control strains
Functional_group2 <- c("GO:0005839", "GO:0008540", "GO:0008541", "dCTRL")
dataset_pca_125mM_Proteasome <- na.omit(Validation_new_df_Grp[which(Validation_new_df_Grp$Group_BY_GO_Terms %in% Functional_group2), c(46, 49, 52, 55, 6, 58)])
res.pca_Proteasome <- prcomp(dataset_pca_125mM_Proteasome[, 1:3], scale = TRUE)

print("GO:0005839 = proteasome core complex")
print("GO:0008540 = proteasome regulatory particle, base subcomplex")
print("GO:0008541 = proteasome regulatory particle, lid subcomplex")
print("dCTRL = Control strains")
fviz_pca_biplot(res.pca_Proteasome,
                col.ind = as.character(dataset_pca_125mM_Proteasome$Group_BY_GO_Terms),
                palette = c("green", "magenta", "blue", "red"),
                repel = T,     # Avoid text overlapping
                addEllipses = T,
                ellipse.type = "confidence",
                #label = "var"
)+
  theme_classic()

Gene_set1_Proteasome <- c("RPN8", "RPN9", "RPN12", "RPT1", "RPT2", "RPT4", "PRE4", "PUP3")

dCTRL_GENES <- c("Ctrl_14", "Ctrl_2",  "Ctrl_23", "Ctrl_28", "Ctrl_30", "Ctrl_32", "Ctrl_34")

barplot_dataset7 <- Validation_new_df_Grp[which(Validation_new_df_Grp$GENE %in% c(Gene_set1_Proteasome, dCTRL_GENES)), c(1, 49, 50)]
colnames(barplot_dataset7)[2:3] <- c("LPI_GT", "SD_GT")
library(reshape)
reshape_barplot_dataset7 <- reshape(data=barplot_dataset7, idvar="gRNA_name",
                                    varying = list(colnames(barplot_dataset7)[2], colnames(barplot_dataset7)[3]),
                                    v.name=c("Mean", "SD"),
                                    times = c("LPI_GT"),
                                    new.row.names = 1:10000,
                                    direction="long")

library(ggplot2)
ggplot(reshape_barplot_dataset7, aes(fill=time, y=Mean, x=gRNA_name)) + 
  geom_bar(position=position_dodge(), stat="identity", color="black", size=0.4, width = 0.6) +
  geom_errorbar( aes(x=gRNA_name, ymin=Mean-SD, ymax=Mean+SD), position=position_dodge(.9), width=0.2, colour="black", alpha=1, size=0.5)+
  scale_fill_manual(values=c("white"))+
  scale_y_continuous(breaks = c(-0.8, -0.5, -0.25, 0, 0.25, 0.5, 0.8),
                     labels = c("-0.8", "-0.5", "-0.25", "0",  "0.25", "0.5", "0.8"),
                     limits = c(-0.35, 0.95))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))

test <- Val_data_curated[as.character(Validation_new_df_Grp[which(Validation_new_df_Grp$GENE %in% Gene_set1_Proteasome), 1]), c(2, 150)]
print("P-value ≤ 0.5")
test[which(test$P.value_GT_125<=0.05), ]

Val_data_lsc_mean <- Val_data_curated[, 2:4]
for(i in 1:nrow(Val_data_curated)){
  Val_data_lsc_mean[i, 4] <- mean(na.omit(as.numeric(Val_data_curated[i, c(78, 84, 90, 96, 102, 108)])))
  Val_data_lsc_mean[i, 5] <- sd(na.omit(as.numeric(Val_data_curated[i, c(78, 84, 90, 96, 102, 108)])))
}
row.names(Val_data_lsc_mean) <- Val_data_lsc_mean$gRNA_name

barplot_dataset8 <- Val_data_lsc_mean[as.character(barplot_dataset7$gRNA_name), c(1, 4:5)]
barplot_dataset8 <- barplot_dataset8[-c(1:7), ]
colnames(barplot_dataset8)[2:3] <- c("LSC_GT", "SD_GT")
library(reshape)
reshape_barplot_dataset8 <- reshape(data=barplot_dataset8, idvar="gRNA_name",
                                    varying = list(colnames(barplot_dataset8)[2], colnames(barplot_dataset8)[3]),
                                    v.name=c("Mean", "SD"),
                                    times = c("LSC_GT"),
                                    new.row.names = 1:10000,
                                    direction="long")

library(ggplot2)
ggplot(reshape_barplot_dataset8, aes(fill=time, y=Mean, x=gRNA_name)) + 
  geom_bar(position=position_dodge(), stat="identity", color="black", size=0.4, width = 0.6) +
  geom_errorbar( aes(x=gRNA_name, ymin=Mean-SD, ymax=Mean+SD), position=position_dodge(.9), width=0.2, colour="black", alpha=1, size=0.5)+
  scale_fill_manual(values=c("white"))+
  scale_y_continuous(breaks = c(-0.8, -0.5, -0.25, 0, 0.25, 0.5, 0.8),
                     labels = c("-0.8", "-0.5", "-0.25", "0",  "0.25", "0.5", "0.8"),
                     limits = c(-0.35, 0.95))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))

LID_strains <- Validation_new_df_Grp$gRNA_name[which(Validation_new_df_Grp$Group_BY_GO_Terms %in% c("GO:0008541"))]
BASE_strains <- Validation_new_df_Grp$gRNA_name[which(Validation_new_df_Grp$Group_BY_GO_Terms %in% c("GO:0008540"))]
CP_strains <- Validation_new_df_Grp$gRNA_name[which(Validation_new_df_Grp$Group_BY_GO_Terms %in% c("GO:0005839"))]

LID_strains_sig <- LID_strains[c(1, 6:9)]
BASE_strains_sig <- c("RPT4-NRg-2")
CP_strains_sig <- CP_strains[c(1:2, 5:7)]

boxplot(as.numeric(as.matrix(Val_whole_data_dCTRL[, c(124, 127, 130)])), 
        as.numeric(as.matrix(Val_data_curated[which(Val_data_curated$gRNA_name %in% LID_strains_sig), c(124, 127, 130)])),
        as.numeric(as.matrix(Val_data_curated[which(Val_data_curated$gRNA_name %in% BASE_strains_sig), c(124, 127, 130)])),
        as.numeric(as.matrix(Val_data_curated[which(Val_data_curated$gRNA_name %in% CP_strains_sig), c(124, 127, 130)])),
        names = c("Control strains", "19S LID", "19S BASE", "20S CP"),
        cex.axis=1.2,
        col = c("green", "red", "blue", "magenta"),
        ylab= "Relative Yield at 125mM Acetic acid")

P_val_dCTRL_LID_Yield <- t.test(as.numeric(as.matrix(Val_whole_data_dCTRL[, c(124, 127, 130)])), 
                                as.numeric(as.matrix(Val_data_curated[which(Val_data_curated$gRNA_name %in% LID_strains_sig), c(124, 127, 130)])))
P_val_dCTRL_LID_Yield$p.value
P_val_dCTRL_BASE_Yield <- t.test(as.numeric(as.matrix(Val_whole_data_dCTRL[, c(124, 127, 130)])), 
                                 as.numeric(as.matrix(Val_data_curated[which(Val_data_curated$gRNA_name %in% BASE_strains_sig), c(124, 127, 130)])))
P_val_dCTRL_BASE_Yield$p.value
P_val_dCTRL_CP_Yield <- t.test(as.numeric(as.matrix(Val_whole_data_dCTRL[, c(124, 127, 130)])), 
                               as.numeric(as.matrix(Val_data_curated[which(Val_data_curated$gRNA_name %in% CP_strains_sig), c(124, 127, 130)])))
P_val_dCTRL_CP_Yield$p.value

boxplot(as.numeric(as.matrix(Val_whole_data_dCTRL[, c(122, 125, 128)])), 
        as.numeric(as.matrix(Val_data_curated[which(Val_data_curated$gRNA_name %in% LID_strains_sig), c(122, 125, 128)])),
        as.numeric(as.matrix(Val_data_curated[which(Val_data_curated$gRNA_name %in% BASE_strains_sig), c(122, 125, 128)])),
        as.numeric(as.matrix(Val_data_curated[which(Val_data_curated$gRNA_name %in% CP_strains_sig), c(122, 125, 128)])),
        names = c("Control strains", "19S LID", "19S BASE", "20S CP"),
        cex.axis=1.2,
        col = c("green", "red", "blue", "magenta"),
        ylab= "Relative Lag phase at 125mM Acetic acid")

P_val_dCTRL_LID_Lag <- t.test(as.matrix(Val_whole_data_dCTRL[, c(122, 125, 128)]), 
                              as.numeric(as.matrix(Val_data_curated[which(Val_data_curated$gRNA_name %in% LID_strains_sig), c(122, 125, 128)])))
P_val_dCTRL_LID_Lag$p.value
P_val_dCTRL_BASE_Lag <- t.test(as.matrix(Val_whole_data_dCTRL[, c(122, 125, 128)]), 
                               as.numeric(as.matrix(Val_data_curated[which(Val_data_curated$gRNA_name %in% BASE_strains_sig), c(122, 125, 128)])))
P_val_dCTRL_BASE_Lag$p.value
P_val_dCTRL_CP_Lag <- t.test(as.matrix(Val_whole_data_dCTRL[, c(122, 125, 128)]), 
                             as.numeric(as.matrix(Val_data_curated[which(Val_data_curated$gRNA_name %in% CP_strains_sig), c(122, 125, 128)])))
P_val_dCTRL_CP_Lag$p.value

