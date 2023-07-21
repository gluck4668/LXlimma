
LXlimma <- function(gene_matrix){


all_packages <- data.frame(installed.packages())

pack <- data.frame(c("BiocManager","openxlsx","dplyr"))

bioc_pack <- data.frame(c("limma"))

pack$type <- pack[,1] %in% all_packages$Package

for (i in 1:nrow(pack)){if (!requireNamespace(pack[i,1], quietly=TRUE))
  install.packages(pack[i,1],update = F,ask = F)}
rm(i)

for (i in 1:nrow(bioc_pack)){if (!requireNamespace(bioc_pack[i,1], quietly=TRUE))
  BiocManager::install (bioc_pack[i,1],update = F,ask = F) }

rm(i)

packages <- c(pack[,1],bioc_pack[,1])

for(i in packages){
  library(i, character.only = T)}

rm(i)
#-----------------------------------------------------------------------------

exp_m <- read.xlsx(gene_matrix)

# 去掉NA
exp_m <- na.omit(exp_m)

# 去掉symbol重复的行，并取重复行的平均值
exp_m <- limma::avereps(exp_m[,-1],ID=exp_m[,1]) 

# log2 transformation
qx <- as.numeric(quantile(exp_m, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||  (qx[6]-qx[1] > 50 && qx[2] > 0)

if (LogC) {exp_m[which(exp_m <= 0)] <- NaN
            exp_m <- log2(exp_m)}

# 组别数据
group_list <- colnames(exp_m)
groups <- gsub("\\d+$", "", group_list) # 去掉字符串后面的数字


#构建design 样本分组矩阵
design <- model.matrix(~ 0 + groups)
rownames(design) = group_list
colnames(design) <- levels(factor(groups))

# 明确 model vs normal
#cts <- paste(levels(factor(groups))[1], levels(factor(groups))[2], sep="-")
#contrast.matrix <- makeContrasts(cts, levels = design)
#contrast.matrix

# 进行差异分析：
#fit <- lmFit(exp_m, design)
#fit <- contrasts.fit(fit, contrast.matrix)
#fit <- eBayes(fit)
#allDiff <- topTable(fit, number = Inf)


fit <- lmFit(exp_m, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- paste(levels(factor(groups))[1], levels(factor(groups))[2], sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B",number=Inf)

tT <- mutate(tT,gene=rownames(tT))

exp_df <- data.frame(exp_m)

exp_df <- mutate(exp_df,gene=rownames(exp_df))

df_file <- dplyr::inner_join(exp_df,tT,by="gene")

rownames(df_file) <- df_file$gene

file_save <- df_file[,-grep("gene",colnames(df_file))]

if(!dir.exists("analysis results"))
  dir.create("analysis results")

file_dir <- paste0("analysis results/statistical analysis result (",cts,").xlsx")

write.xlsx(file_save, file_dir, rowNames=T)

print("The statistical analysis result can be found in the folder of <analysis results>")

}





