library(readxl)
library(dplyr)
library(stringr)
library(tidyr)

setwd("/Users/shenjingting/Desktop/face_database/we_WGS/6.analysis_SV/03.1kg_sv/Han_sv_annotation")

infile <- "AnnotSV_Dd1JiHtyhX.xlsx"
raw <- readxl::read_excel(infile, sheet = 1)
names(raw) <- str_trim(names(raw))

need <- c("AnnotSV ID","ACMG class","SV type","Gene name","Left breakpoint annotations",
          "Samples_ID","RANK")
miss <- setdiff(need, names(raw))
if (length(miss) > 0) stop("缺少列: ", paste(miss, collapse = ", "))

dat0 <- raw %>%
  mutate(
    ACMG_Class_num = suppressWarnings(as.integer(str_extract(`ACMG class`, "\\d+"))),
    left_bp = as.character(`Left breakpoint annotations`),
    svtype_text = as.character(`SV type`),
    gene_text = as.character(`Gene name`)
  ) %>%
  filter(ACMG_Class_num %in% c(4L, 5L)) %>%
  filter(!is.na(left_bp) & str_trim(left_bp) != "") %>%
  filter(!str_detect(str_to_lower(left_bp), "high\\s*signal"))

# ---- 1) 从 SV type 字段拆出 SVTYPE / SV_length / Gene_count
dat1 <- dat0 %>%
  mutate(
    SVTYPE    = str_extract(svtype_text, "\\b(DEL|DUP|INS|INV|BND|TRA|CNV)\\b"),
    SV_length = suppressWarnings(as.integer(str_match(svtype_text, "SV_length\\s*:\\s*(\\d+)")[,2])),
    Gene_count= suppressWarnings(as.integer(str_match(svtype_text, "Gene_count\\s*:\\s*(\\d+)")[,2]))
  )

# ---- 2) 抽取 AnnotSV_Rank_Score（在 "ACMG class" 列里那段文本中）
dat2 <- dat1 %>%
  mutate(
    AnnotSV_Rank_Score = suppressWarnings(as.numeric(
      str_match(as.character(`ACMG class`), "AnnotSV_ranking_score\\s*:\\s*([0-9.]+)")[,2]
    ))
  )

# ---- 3) 抽取 Cytoband：优先从任意列文本里抓 "CytoBand : xxx"
# （你的预览里 CytoBand 在同一行块里，但不确定是否在单独列；这里做“整行拼接”提取最稳）
dat3 <- dat2 %>%
  mutate(
    row_all = apply(across(everything(), ~ifelse(is.na(.x), "", as.character(.x))), 1, paste, collapse = " || "),
    Cytoband = str_match(row_all, "CytoBand\\s*:\\s*([0-9A-Za-z\\.\\-]+)")[,2]
  ) %>%
  dplyr::select(-row_all)

# ---- 4) 只从 Gene name 列的 "Gene list :" 长串里取基因，并拆分
gene_long <- dat3 %>%
  mutate(
    gene_list_str = str_match(gene_text, "(?is)Gene\\s*list\\s*:\\s*(.*)")[,2],
    gene_list_str = ifelse(is.na(gene_list_str), "", gene_list_str),
    gene_list_str = str_split(gene_list_str, "\\r?\\n", simplify = TRUE)[,1],
    gene_list_str = str_squish(gene_list_str)
  ) %>%
  filter(gene_list_str != "") %>%
  separate_rows(gene_list_str, sep = "\\s*;\\s*") %>%
  mutate(gene = str_trim(gene_list_str)) %>%
  filter(gene != "" & gene != ".") %>%
  transmute(
    AnnotSV_ID = `AnnotSV ID`,
    Samples_ID,
    ACMG_Class_num,
    AnnotSV_Rank_Score,
    Cytoband,
    SVTYPE, SV_length, Gene_count,
    gene
  ) %>%
  distinct()
gene_long[c("a", "sv_id", "b", "n_samples")] <- str_split_fixed(gene_long$Samples_ID, "_", 4)
colnames(gene_long)
gene_long2 <- dplyr::select(gene_long, -c("a", "b"))
colnames(gene_long2)

print(length(gene_long2$gene))
print(length(unique(gene_long2$gene)))

# wide：每个SV一行，基因合并
gene_wide <- gene_long2 %>%
  group_by(AnnotSV_ID, Samples_ID, ACMG_Class_num, AnnotSV_Rank_Score, Cytoband,
           SVTYPE, SV_length, Gene_count, sv_id, n_samples) %>%
  summarise(genes = paste(sort(unique(gene)), collapse = ";"), .groups = "drop")

#write.csv(gene_long2, "Han_sv_ACMG45_long_gene.csv", row.names = FALSE)
#write.csv(gene_wide, "Han_sv_ACMG45_wide.csv", row.names = FALSE)


# ---- 5) 根据cytoband + svtype统计携带者人数
colnames(gene_long)
cytoband_sv <- paste(gene_long$Cytoband, gene_long$SVTYPE, sep="_")
gene_long["cytoband_sv"] <- cytoband_sv

gene_long2 <- gene_long %>%
  group_by(cytoband_sv) %>%
  mutate(n_samples = n_distinct(Samples_ID)) %>%
  ungroup() %>%
  filter(n_samples < 5)

# wide：每个SV一行，基因合并
gene_wide <- gene_long2 %>%
  group_by(AnnotSV_ID, Samples_ID, ACMG_Class_num, AnnotSV_Rank_Score, Cytoband,
           SVTYPE, SV_length, Gene_count, cytoband_sv, n_samples) %>%
  summarise(genes = paste(sort(unique(gene)), collapse = ";"), .groups = "drop")

#write.csv(gene_long2, "52case_sv_ACMG45_long_gene.csv", row.names = FALSE)
#write.csv(gene_wide, "52case_sv_ACMG45_wide.csv", row.names = FALSE)

#gene_wide <- read.csv("52case_sv_ACMG45_wide.csv")
# ---- 6) 分DUP和DEL, 对携带人数进行统计，并可视化
library(forcats)

colnames(gene_wide)
dat_plot <- gene_wide %>% dplyr::select("ACMG_Class_num", "Cytoband", "SVTYPE", "n_samples")
dat_plot <- unique(dat_plot)
colnames(dat_plot)
dat_plot["SVTYPE_ACMG"] <- paste0(dat_plot$ACMG_Class_num, dat_plot$SVTYPE)
table(dat_plot$SVTYPE_ACMG)
dat_plot <- dat_plot[dat_plot$n_samples > 1, ]
plot_bar <- function(dat, svtype = c("DEL","DUP")){
  svtype <- match.arg(svtype)
  dat %>%
    mutate(
      SVTYPE = toupper(SVTYPE),
      n_samples = as.numeric(n_samples)
    ) %>%
    filter(SVTYPE == svtype) %>%
    mutate(Cytoband = forcats::fct_reorder(Cytoband, n_samples, .desc = TRUE)) %>%
    ggplot(aes(x = Cytoband, y = n_samples, fill=SVTYPE_ACMG)) +
    geom_bar(stat = "identity") + 
    scale_fill_manual(values = c("5DUP"="#d73027", "4DUP"="pink", "5DEL"="#4575b4", "4DEL"="lightblue")) +
    labs(x = "Cytoband", y = "", title = paste0(svtype, " CNVs")) +
    theme_bw(base_size = 16) +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
      panel.grid.major.x = element_blank()
    )
}

p_del <- plot_bar(dat_plot, "DEL")
p_dup <- plot_bar(dat_plot, "DUP")

p_del
p_dup
