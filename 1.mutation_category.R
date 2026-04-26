#针对57case 和 57control，分析变异位点的种类
#主要区分LGD,D-mis,D-mis-HC，Damaging(LGD+D-mis),Damaging HC(LGD+D-mis+HC)

setwd("/Users/shenjingting/Desktop/face_database/we_WGS/1.analysis")
library(readxl)
library(openxlsx)
library(dplyr)
library(tidyr)

#针对case
#读取变异注释结果
list.files()
case <- read_excel("TOTAL_CASE_58.xlsx", sheet = 1)
names(case)
colnames(case)[1] <- "file"
colnames(case)[2] <- "index"
colnames(case)[10] <- "INFO"
case <- case %>% filter(case$POS != "POS")

## 转换gnomad_genome_eas_af, QUAL和QD为数值型
case$gnomad_genome_eas_af <- as.numeric(case$gnomad_genome_eas_af) #转化为数值，引入NA
summary(case$gnomad_genome_eas_af)
case$QUAL <- as.numeric(case$QUAL) #转化为数值
summary(case$QUAL)
case[is.na(case)] <- 0 #将NA替换为0

#条件1：gnomad_genome_eas_af < 0.0001 过滤稀有变异
#条件2: 要求depth > 15reads，Alt > 8 reads 支持
#条件2：功能分类
## synonymous: 同义突变
## LGD: 无义(终止密码子提前）、移码、剪接位点突变。
## D-mis: CADD_score >=20 AND (SIFT="D" OR "PolyPhen"="D")
## D-mis-HC: MetaSVM22 = "D" or CADD score >=30 
## LGD + D-mis (Damaging)
## LGD + D-mis-HC (Damaging HC)

summary(case$gnomad_genome_eas_af)
summary(case$QUAL)

## QUAL > 30, QD > 2 #高质量变异筛选
filtered_by_qual <- subset(case, case$QUAL >= 30)

#罕见变异(<0.001)
filtered_by_gnomAD <- subset(filtered_by_qual, filtered_by_qual$gnomad_genome_eas_af <= 0.001)
summary(filtered_by_gnomAD$gnomad_genome_eas_af)


#>15 reads 支持
filtered_by_gnomAD <- filtered_by_gnomAD %>%
  separate(INFO, into=c("GT", "AD", "DP", "GQ", "PL"),
           sep = ":", remove=FALSE) %>%
  mutate(DP = as.numeric(DP))

filtered_by_gnomAD <- filtered_by_gnomAD %>%
  separate(AD, into=c("Ref_D", "Alt_D"),
           sep = ",", remove=FALSE) %>%
  mutate(Ref_D = as.numeric(Ref_D),
         Alt_D = as.numeric(Alt_D))

filtered_by_reads <- subset(filtered_by_gnomAD, filtered_by_gnomAD$DP > 15 & filtered_by_gnomAD$Alt_D > 8)
summary(filtered_by_reads$DP)
summary(filtered_by_reads$Alt_D)
write.csv(filtered_by_reads, "1.zoutput/high_confidence_rare_variants.csv")
#共找到10508个variants(高置信度的罕见变异)

## synonymous
category0 <- filtered_by_reads$Consequence == "synonymous_variant"

##LGD
LGD_list <- c("frameshift_variant", "frameshift_variant,splice_region_variant", "frameshift_variant,start_lost,start_retained_variant",
              "frameshift_variant,stop_lost", "start_lost", "stop_lost","stop_retained_variant",
              "stop_gained", "stop_gained,frameshift_variant", "stop_gained,splice_region_variant",
              "splice_acceptor_variant", "splice_donor_variant", "splice_region_variant",
              "splice_acceptor_variant,coding_sequence_variant",
              "splice_donor_5th_base_variant,intron_variant",
              "splice_polypyrimidine_tract_variant,intron_variant"
              )
category1 <- filtered_by_reads$Consequence %in% LGD_list

##D-mis
category2 <-  (filtered_by_reads$CADD_phred >= 20) & 
  (filtered_by_reads$SIFT_pred == "D" | filtered_by_reads$Polyphen2_HVAR_pred == "D")

##D-mis-HC
category3 <- (filtered_by_reads$SIFT_pred == "D" | filtered_by_reads$Polyphen2_HVAR_pred == "D") &
  (filtered_by_reads$MetaSVM_pred == "D" | filtered_by_reads$CADD_phred > 30) 

#LGD + D-mis
category4 <- filtered_by_reads$Consequence %in% LGD_list &  (filtered_by_reads$CADD_phred >= 20) & 
  (filtered_by_reads$SIFT_pred == "D" | filtered_by_reads$Polyphen2_HVAR_pred == "D")

#LGD + D-mis-HC
category5 <- filtered_by_reads$Consequence %in% LGD_list & (filtered_by_reads$SIFT_pred == "D" | filtered_by_reads$Polyphen2_HVAR_pred == "D") &
  (filtered_by_reads$MetaSVM_pred == "D" | filtered_by_reads$CADD_phred > 30) 

##输出结果
syn_case <- filtered_by_reads %>% filter((category0))
LGD_case <- filtered_by_reads %>% filter((category1))
Dmis_case <- filtered_by_reads %>% filter((category2))
DmisHC_case <- filtered_by_reads %>% filter((category3))
Damaging <- filtered_by_reads %>% filter((category4))
Damaging_HC <- filtered_by_reads %>% filter((category5))

#输出表格
write.xlsx(list(Sheet1=syn_case, Sheet2=LGD_case,Sheet3=Dmis_case,Sheet4=DmisHC_case), file="1.zoutput/case_output.xlsx")

##############################华丽丽的分割线####################################################
#针对control
#读取变异注释结果
list.files()
control <- read_excel("TOTAL_CONTROL_57.xlsx", sheet = 1)
names(control)
colnames(control)[1] <- "file"
colnames(control)[2] <- "index"
colnames(control)[10] <- "INFO"
control <- control %>% filter(control$POS != "POS")

## 转换gnomad_genome_eas_af为数值型
control$gnomad_genome_eas_af <- as.numeric(control$gnomad_genome_eas_af) #转化为数值，引入NA
control$QUAL <- as.numeric(control$QUAL)
control[is.na(control)] <- 0 #将NA替换为0

#条件1：gnomad_genome_eas_af < 0.001 过滤稀有变异
#条件2: 要求depth > 15reads，Alt > 8 reads 支持
#条件2：功能分类
## synonymous: 同义突变
## LGD: 无义(终止密码子提前）、移码、剪接位点突变。
## D-mis: CADD_score >=20 AND (SIFT="D" OR "PolyPhen"="D")
## D-mis-HC: MetaSVM22 = "D" or CADD score >=30 
## LGD + D-mis (Damaging)
## LGD + D-mis-HC (Damaging HC)
summary(control$QUAL)
summary(control$gnomad_genome_eas_af)

## QUAL > 30, QD > 2 #高质量变异筛选
filtered_by_qual <- subset(control, control$QUAL >= 30)

#罕见变异
filtered_by_gnomAD <- subset(filtered_by_qual, filtered_by_qual$gnomad_genome_eas_af <= 0.001)
summary(filtered_by_gnomAD$gnomad_genome_eas_af)

#>15 reads 支持
filtered_by_gnomAD <- filtered_by_gnomAD %>%
  separate(INFO, into=c("GT", "AD", "DP", "GQ", "PL"),
           sep = ":", remove=FALSE) %>%
  mutate(DP = as.numeric(DP))

filtered_by_gnomAD <- filtered_by_gnomAD %>%
  separate(AD, into=c("Ref_D", "Alt_D"),
           sep = ",", remove=FALSE) %>%
  mutate(Ref_D = as.numeric(Ref_D),
         Alt_D = as.numeric(Alt_D))

filtered_by_reads <- subset(filtered_by_gnomAD, filtered_by_gnomAD$DP > 15 & filtered_by_gnomAD$Alt_D > 8)
summary(filtered_by_reads$DP)
summary(filtered_by_reads$Alt_D)

#共找到10604个variants(高置信度的罕见变异)

## synonymous
category0 <- filtered_by_reads$Consequence == "synonymous_variant"

##LGD
LGD_list <- c("frameshift_variant", "frameshift_variant,splice_region_variant", "frameshift_variant,start_lost,start_retained_variant",
              "frameshift_variant,stop_lost", "start_lost", "stop_lost","stop_retained_variant",
              "stop_gained", "stop_gained,frameshift_variant", "stop_gained,splice_region_variant",
              "splice_acceptor_variant", "splice_donor_variant", "splice_region_variant",
              "splice_acceptor_variant,coding_sequence_variant",
              "splice_donor_5th_base_variant,intron_variant",
              "splice_polypyrimidine_tract_variant,intron_variant"
)
category1 <- filtered_by_reads$Consequence %in% LGD_list

##D-mis
category2 <-  (filtered_by_reads$CADD_phred >= 20) & 
  (filtered_by_reads$SIFT_pred == "D" | filtered_by_reads$Polyphen2_HVAR_pred == "D")

##D-mis-HC
category3 <- (filtered_by_reads$SIFT_pred == "D" | filtered_by_reads$Polyphen2_HVAR_pred == "D") &
  (filtered_by_reads$MetaSVM_pred == "D" | filtered_by_reads$CADD_phred > 30) 

#LGD + D-mis
category4 <- filtered_by_reads$Consequence %in% LGD_list &  (filtered_by_reads$CADD_phred >= 20) & 
  (filtered_by_reads$SIFT_pred == "D" | filtered_by_reads$Polyphen2_HVAR_pred == "D")

#LGD + D-mis-HC
category5 <- filtered_by_reads$Consequence %in% LGD_list & (filtered_by_reads$SIFT_pred == "D" | filtered_by_reads$Polyphen2_HVAR_pred == "D") &
  (filtered_by_reads$MetaSVM_pred == "D" | filtered_by_reads$CADD_phred > 30) 

##输出结果
syn_control <- filtered_by_reads %>% filter((category0))
LGD_control <- filtered_by_reads %>% filter((category1))
Dmis_control <- filtered_by_reads %>% filter((category2))
DmisHC_control <- filtered_by_reads %>% filter((category3))
Damaging <- filtered_by_reads %>% filter((category4))
Damaging_HC <- filtered_by_reads %>% filter((category5))

#输出表格
write.xlsx(list(Sheet1=syn_control, Sheet2=LGD_control,Sheet3=Dmis_control,Sheet4=DmisHC_control), file="1.zoutput/control_output.xlsx")
