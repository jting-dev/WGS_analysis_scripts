setwd("/Users/shenjingting/Desktop/face_database/we_WGS/1.analysis/")
library(readxl)
library(openxlsx)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

##根据Damaging基因的统计新结果表中的携带基因突变的人数统计，进行fisher精确性检验。
summary <- read_excel("5.zoutput/Damaging基因的统计新结果.xlsx", sheet = 5) #重新定义身份(Dmis-Dmishc)

#############################统计检验函数构建####################################
fisher_analysis <- function(input_df, output_file, method = "BH") {
  # 验证校正方法有效性
  valid_methods <- c("holm", "hochberg", "hommel", "bonferroni", 
                     "BH", "fdr", "BY", "none")
  if (!method %in% valid_methods) {
    stop(paste("无效的校正方法。可选方法：", paste(valid_methods, collapse = ", ")))
  }
  
  # 读取数据并验证格式
  tryCatch({
    df <- input_df
    required_cols <- c("Gene", "Case_Mut", "Control_Mut")
    if (!all(required_cols %in% colnames(df))) {
      stop("Excel文件必须包含以下列：Gene, Case_Mut, Control_Mut")
    }
  }, error = function(e) {
    stop(paste("数据读取错误:", e$message))
  })
  
  #过滤突变人数>=3的基因
  df <- subset(df, df$Case_Mut >= 3)
  print (paste0("保留病例组大于3的变异进行负担分析，", nrow(df)))
  #只保留病例组中大于对照组中的基因
  df_new <- df[df$Case_Mut > df$Control_Mut, ]
  print (paste0("保留病例组大于对照组变异进行负担分析，", nrow(df_new)))
  df = df_new
  # 初始化结果数据框
  results <- data.frame(
    Gene = character(),
    Case_Mut = integer(),
    Control_Mut = integer(),
    Case_LGD = integer(),
    Case_Dmis = integer(),
    Case_Dmis_HC = integer(),
    Control_LGD = integer(),
    Control_Dmis = integer(),
    Control_Dmis_HC = integer(),
    P_value = numeric(),
    Adjusted_P = numeric(),
    Method = character(),
    OR = numeric(),
    Lower_CI = numeric(),
    Upper_CI = numeric(),
    Significant = logical(),
    stringsAsFactors = FALSE
  )
  
  # 执行批量分析
  p_values <- c()
  for (i in 1:nrow(df)) {
    # 提取数据
    gene <- df$Gene[i]
    case_mut <- df$Case_Mut[i]
    control_mut <- df$Control_Mut[i]
    case_lgd <- df$Case_LGD[i]
    case_dmis <- df$Case_Dmis[i]
    case_dmis_hc <- df$Case_Dmis_HC[i]
    control_lgd <- df$Control_LGD[i]
    control_dmis <- df$Control_Dmis[i]
    control_dmis_hc <- df$Control_Dmis_HC[i]
    
    # 数据验证
    if (!is.numeric(case_mut) || !is.numeric(control_mut)) {
      warning(paste("基因", gene, "突变数非数值，跳过处理"))
      next
    }
    if (case_mut < 0 || control_mut < 0) {
      warning(paste("基因", gene, "突变数为负数，跳过处理"))
      next
    }
    
    # 构建2x2列联表
    table <- matrix(
      c(case_mut, 58 - case_mut,
        control_mut, 103 - control_mut),
      nrow = 2, byrow = TRUE
    )
    
    # 执行Fisher精确检验
    tryCatch({
      result <- fisher.test(table, alternative = "greater")
      
      # 保存原始p值
      p_values <- c(p_values, result$p.value)
      
      # 保存结果
      new_row <- data.frame(
        Gene = gene,
        Case_Mut = case_mut,
        Control_Mut = control_mut,
        Case_LGD = case_lgd,
        Case_Dmis = case_dmis,
        Case_Dmis_HC = case_dmis_hc,
        Control_LGD = control_lgd,
        Control_Dmis = control_dmis,
        Control_Dmis_HC = control_dmis_hc,
        P_value = result$p.value,
        Adjusted_P = NA,
        Method = method,
        OR = result$estimate,
        Lower_CI = result$conf.int[1],
        Upper_CI = result$conf.int[2],
        Significant = (result$p.value < 0.05),
        stringsAsFactors = FALSE
      )
      
      results <- bind_rows(results, new_row)
      
    }, error = function(e) {
      warning(paste("基因", gene, "计算错误:", e$message))
    })
  }
  
  # 执行多重检验校正（跳过无p值的情况）
  if(length(p_values) > 1) {
    adj_p <- p.adjust(p_values, method = method)
    results$Adjusted_P <- adj_p[order(match(results$Gene, df$Gene))]
  }
  
  # 更新显著性判断
  results$Significant <- ifelse(results$Adjusted_P < 0.05, TRUE, FALSE)
  
  # 添加显著性标记
  results$Significance <- ifelse(results$Significant, "*", "")
  
  # 保存结果到Excel
  write.xlsx(results, output_file, rowNames = FALSE)
  message("分析完成！结果已保存至 ", output_file)
}

# 使用示例（支持多种校正方法）
# 20260126 修正了5.zoutput/Damaging基因的统计新结果.xlsx CEP290是不是要加1个？


table(summary$Case_Mut > summary$Control_Mut)

fisher_analysis(
  input_df = summary,
  output_file = "6.zoutput/fisher_results_fdr_20260308.xlsx",
  method = "fdr"  # 可选值："holm", "hochberg", "bonferroni", "BH", "BY", "none"
)
##################################################################################
# 使用示例（自动选择检验方法）
#fisher_analysis(
#  input_file = "5.zoutput/Damaging基因的统计新结果.xlsx",
#  output_file = "6.zoutput/fisher_results.xlsx",
#  method = "BH"  # 可选值："holm", "hochberg", "bonferroni", "BH", "BY", "none"
#)

#######################################拓展可视化################################
result_df <- read_excel("6.zoutput/fisher_results_fdr.xlsx", sheet = 2)

ggplot(result_df, aes(x = reorder(Gene, -P_value), y = P_value)) +
    geom_point(size = 1) +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
    scale_color_manual(values = c("black", "red")) +
    labs(
      title = "Fisher exact results",
      x = "Gene",
      y = "P. value",
      color = "significant"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      plot.title = element_text(hjust = 0.5)
    ) 

result_df$Adjusted_P

ggplot(result_df, aes(x = reorder(Gene, -Adjusted_P), y = Adjusted_P)) +
  geom_point( size = 1) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  scale_color_manual(values = c("black", "red")) +
  labs(
    title = "Fisher exact results",
    x = "Gene",
    y = "Adjusted P.value"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  ) 



####################使用精确无条件检验############################################
#精确无条件检验，基于【得分统计量】（score statistic）检验两个比例相等的零假设，比fihser检验更少保守，更容易达到显著。
#当期望频数<5时，卡方检验不适用，Fisher检验虽精确但保守（p值偏大）。
#无条件方法通过考虑所有可能的参数值（而非固定边缘总和），通常比Fisher检验更少保守，更易达到显著。
#更适用于检测小样本队列中，病例组罕见变异的富集性。
#install.packages("exact2x2")
library(exact2x2)
#demo
demo <- uncondExact2x2(10, 58, 4, 100,
                       alternative = "greater",
                       method = "score")
demo



################对不同类型的变异进行精确检验#####################################
fisher_analysis <- function(input_file, output_file, method = "BH") {
  # 验证校正方法有效性
  valid_methods <- c("holm", "hochberg", "hommel", "bonferroni", 
                     "BH", "fdr", "BY", "none")
  if (!method %in% valid_methods) {
    stop(paste("无效的校正方法。可选方法：", paste(valid_methods, collapse = ", ")))
  }
  
  # 读取数据并验证格式
  tryCatch({
    df <- read.xlsx(input_file, sheet = 4)
    required_cols <- c("Gene", "Case_Mut", "Control_Mut")
    if (!all(required_cols %in% colnames(df))) {
      stop("Excel文件必须包含以下列：Gene, Case_Mut, Control_Mut")
    }
  }, error = function(e) {
    stop(paste("数据读取错误:", e$message))
  })
  
  df <- subset(df, df$Case_Mut >= 3)
  
  # 初始化结果数据框
  results <- data.frame(
    Gene = character(),
    Case_Mut = integer(),
    Control_Mut = integer(),
    Case_LGD = integer(),
    Control_LGD = integer(),
    Case_Dmis = integer(),
    Control_Dmis = integer(),
    Case_Dmis_HC = integer(),
    Control_Dmis_HC = integer(),
    P_value = numeric(),
    Adjusted_P = numeric(),
    P_value_LGD = numeric(),
    Adjusted_P_LGD = numeric(),
    P_value_Dmis = numeric(),
    Adjusted_P_Dmis = numeric(),
    P_value_Dmis_HC = numeric(),
    Adjusted_P_Dmis_HC = numeric(),
    Method = character(),
    OR = numeric(),
    Lower_CI = numeric(),
    Upper_CI = numeric(),
    stringsAsFactors = FALSE
  )
  
  # 执行批量分析
  p_values <- c()
  p_values_LGD <- c()
  p_values_Dmis <- c()
  p_values_DmisHC <- c()
  
  for (i in 1:nrow(df)) {
    # 提取数据
    gene <- df$Gene[i]
    case_mut <- df$Case_Mut[i]
    control_mut <- df$Control_Mut[i]
    case_LGD <- df$Case_LGD[i]
    control_LGD <- df$Control_LGD[i]
    case_Dmis <- df$Case_Dmis[i]
    control_Dmis <- df$Control_Dmis[i]
    case_Dmis_HC <- df$Case_Dmis_HC[i]
    control_Dmis_HC <- df$Control_Dmis_HC[i]
    
    # 构建2x2列联表
    table <- matrix(
      c(case_mut, 58 - case_mut,
        control_mut, 103 - control_mut),
      nrow = 2, byrow = TRUE
    )
    
    table_LGD <- matrix(
      c(case_LGD, 58 - case_LGD,
        control_LGD, 103 - control_LGD),
      nrow = 2, byrow = TRUE
    )
    
    table_Dmis <- matrix(
      c(case_Dmis, 58 - case_Dmis,
        control_Dmis, 103 - control_Dmis),
      nrow = 2, byrow = TRUE
    )
    
    table_Dmis_HC <- matrix(
      c(case_Dmis_HC, 58 - case_Dmis_HC,
        control_Dmis_HC, 103 - control_Dmis_HC),
      nrow = 2, byrow = TRUE
    )
    
    # 执行Fisher精确检验
    tryCatch({
      result <- fisher.test(table, alternative = "two.sided")
      result_LGD <- fisher.test(table_LGD, alternative = "two.sided")
      result_Dmis <- fisher.test(table_Dmis, alternative = "two.sided")
      result_DmisHC <- fisher.test(table_Dmis_HC, alternative = "two.sided")
      
      # 保存原始p值
      p_values <- c(p_values, result$p.value)
      p_values_LGD <- c(p_values_LGD, result_LGD$p.value)
      p_values_Dmis <- c(p_values_Dmis, result_Dmis$p.value)
      p_values_DmisHC<- c(p_values_DmisHC, result_DmisHC$p.value)
      
      # 保存结果
      new_row <- data.frame(
        Gene = gene,
        Case_Mut = case_mut,
        Control_Mut = control_mut,
        Case_LGD = case_LGD,
        Control_LGD = control_LGD,
        Case_Dmis = case_Dmis,
        Control_Dmis = control_Dmis,
        Case_Dmis_HC = case_Dmis_HC,
        Control_Dmis_HC = control_Dmis_HC,
        P_value = result$p.value,
        P_value_LGD = result_LGD$p.value,
        P_value_Dmis = result_Dmis$p.value,
        P_value_Dmis_HC = result_DmisHC$p.value,
        Method = method,
        OR = result$estimate,
        Lower_CI = result$conf.int[1],
        Upper_CI = result$conf.int[2],
        stringsAsFactors = FALSE
      )
      
      results <- bind_rows(results, new_row)
      
    }, error = function(e) {
      warning(paste("基因", gene, "计算错误:", e$message))
    })
  }
  
  # 执行多重检验校正（跳过无p值的情况）
  if(length(p_values) > 1) {
    adj_p <- p.adjust(p_values, method = method)
    adj_p_LGD <- p.adjust(p_values_LGD, method = method)
    adj_p_Dmis <- p.adjust(p_values_Dmis, method = method)
    adj_p_DmisHC<- p.adjust(p_values_DmisHC, method = method)
    
    results$Adjusted_P <- adj_p[order(match(results$Gene, df$Gene))]
    results$Adjusted_P_LGD <- adj_p_LGD[order(match(results$Gene, df$Gene))]
    results$Adjusted_P_Dmis <- adj_p_Dmis[order(match(results$Gene, df$Gene))]
    results$Adjusted_P_Dmis_HC <- adj_p_DmisHC[order(match(results$Gene, df$Gene))]
  }
  
  
  # 保存结果到Excel
  write.xlsx(results, output_file, rowNames = FALSE)
  message("分析完成！结果已保存至 ", output_file)
}

# 使用示例（支持多种校正方法）
summary <- read_excel("5.zoutput/Damaging基因的统计新结果.xlsx", sheet = 4)

fisher_analysis(
  input_file = "5.zoutput/Damaging基因的统计新结果.xlsx",
  output_file = "6.zoutput/fisher_results_fdr_category_103.xlsx",
  method = "fdr"  # 可选值："holm", "hochberg", "bonferroni", "BH", "BY", "none"
)

