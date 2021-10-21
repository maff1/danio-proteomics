################################################################################
# preprocess and run LIMMA ## v001 ## MAFF #####################################

rm(list = ls())

library(data.table)
library(limma)
library(Biobase)
library(missRanger)
# ---------------------------------------------------------------------------- #

f_clean_dep <-  function(pathDat, impute = TRUE, method = "addZero") {
  
dfNames <- c("uniprot", "symbol",
  toupper(make.names(
  data.table::fread(file = pathDat, 
                             sep = "\t", 
                             select = c(13:27), 
                             dec = ",", 
                             header = FALSE, 
                             skip = 3,
                             nrows = 1,
                             data.table = FALSE)
    )
  )
)

df <- data.table::fread(file = pathDat, 
                            sep = "\t", 
                            select = c(5:6, 13:27), 
                            dec = ",", 
                            header = FALSE, 
                            skip = 4,
                            data.table = FALSE)

df <- as.data.frame(
  cbind(df[, 1:2],
        apply(df[, 3:ncol(df)], 
              2, 
              as.numeric)
        )
)

df$V5 <- trimws(gsub(" .*|-DECOY", "", df$V5))
df$V6 <- trimws(df$V6)
df$V6 <- ifelse(df$V6 == "", df$V5, df$V6)
colnames(df) <- dfNames
df <- df[-as.numeric(rownames(df[duplicated(df$uniprot),])),]

metaDat <- df[, 1:2]
exprDat <- df[, 3:ncol(df)]

expr <- as.data.frame(t(
  data.frame(
  exprDat,
  row.names = metaDat$uniprot
  )
))

colsNA = data.frame(percNA = colSums(is.na(expr)/nrow(expr)*100))
print(paste0("Excluded ", length(rownames(subset(colsNA, percNA == 100))),
             " proteins with all zero........................................."))
`%nin%` <- Negate(`%in%`)
expr <- expr[, colnames(expr) %nin% rownames(subset(colsNA, percNA == 100))]
metaDat <- metaDat[metaDat$uniprot %nin% rownames(subset(colsNA, percNA == 100)),]
saveRDS(metaDat, "./RESULTS/metaData_tramadol.rds")

# IMPUTATION ----------------------------------------------------------------- #
if(impute == TRUE & method == "addZero") {
  expr[is.na(expr)] <- 0
    } else if(impute == TRUE & method == "randomForest") {
    expr <- data.frame(
      missRanger::missRanger(
      data = expr, formula = .~., seed = 2021, num.trees = 10
    ), row.names = rownames(expr)
  )
  } else if(impute == FALSE & method == "none") {
    expr
  }

saveRDS(expr, paste0("./RESULTS/", make.names(pathDat), ".rds"))

# REFORMAT TO ACCOMMODATE FOR EXPRESSION SET OBJECT -------------------------- #
samples <- data.frame(samples = rownames(expr), 
                      groups = gsub("\\..*", "", 
                                    rownames(expr)),
                      row.names = rownames(expr)
)

toMatch = data.frame(
    uniprot = colnames(expr),
    row.names = colnames(expr)
)
anno = merge(toMatch, metaDat, by = "uniprot", all.x = TRUE, sort = FALSE)
anno = data.frame(
  uniprot = anno$uniprot,
  symbol = anno$symbol,
  row.names = anno$uniprot
)
  
phenoData <- new("AnnotatedDataFrame", data = samples)
featureData <- new("AnnotatedDataFrame", data = anno)
eset <- ExpressionSet(
  assayData = as.matrix(t(expr)), 
  phenoData = phenoData,
  featureData = featureData
)

# quantile normalization ----------------------------------------------------- #
quantile_normalisation <- function(df){
  df_rank <- apply(df,2,rank,ties.method="min")
  df_sorted <- data.frame(apply(df, 2, sort))
  df_mean <- apply(df_sorted, 1, mean)
  
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
  
  df_final <- as.data.frame(apply(df_rank, 2, index_to_mean, my_mean=df_mean))
  rownames(df_final) <- rownames(df)
  return(df_final)
}

neset <- ExpressionSet(
  assayData = as.matrix(
    log2(quantile_normalisation(exprs(eset)) + 1)), 
  phenoData = phenoData,
  featureData = featureData
)

pdf(file = "./RESULTS/quantile.normalisation.pdf")
boxplot(exprs(eset), outline=FALSE, horizontal = TRUE, main = "before quantile-normalisation")
boxplot(exprs(neset), outline=FALSE, horizontal = TRUE, main = "after quantile-normalisation & log2(x+1)")
dev.off()

# sample OTRA100.1 has extreme mean differences within the OTRA100 group
print("removed sample OTRA100.1 due to extreme mean differences within the OTRA100 group")
filterSample <- colnames(neset)[neset@phenoData@data$samples %nin% "OTRA100.1"]
neset <- neset[, filterSample]

# LIMMA DPE ------------------------------------------------------------------ #

design <- model.matrix(~0 + groups, data = neset)
colnames(design) <- gsub("groups", "", colnames(design))
contrast =  makeContrasts(contrasts=c("OTRA0-CTR", "OTRA100-CTR", "TRA0-CTR",
                                      "TRA100-CTR"), levels=design)
fit1 <- lmFit(exprs(neset), design)
fit2 <- contrasts.fit(fit1, contrasts = contrast)
fit3 <- eBayes(fit2, trend = TRUE)

lsContrast = list("OTRA0-CTR", "OTRA100-CTR", "TRA0-CTR","TRA100-CTR")
names(lsContrast) <- c("OTRA0-CTR", "OTRA100-CTR", "TRA0-CTR","TRA100-CTR")
lsRes = lapply(lsContrast, function(cf) topTable(fit3, 
                                                 coef=cf, 
                                                 n=Inf, 
                                                 confint = TRUE, 
                                                 sort.by="none")
)
lsResAnn = lapply(lsRes, function(X)
  cbind.data.frame(
    neset@featureData@data,
    X
  )
)
openxlsx::write.xlsx(lsResAnn, "./RESULTS/limma_DEP_tramadol.xlsx", rowNames = TRUE)
# ---------------------------------------------------------------------------- #

return(neset)
}

res <- f_clean_dep(
  pathDat = "Samples View Report for Dataset 2 OTRA 210813 correct.xls",
  impute = TRUE, method = "addZero"
)
# ---------------------------------------------------------------------------- #