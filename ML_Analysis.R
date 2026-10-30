#Extract the GSE151052 dataset and prepare the training data
rm(list = ls())
load("/home/hzy-sc1/AAAAA_bioinformatics_workflow/microarray-pipeline/GEO_public_data/external_validation_processed_GSE32537_IPF_human_lung_tissue/GSE32537_processed.RData")
ls()
exp[1:4,1:4]
table(Group)
gene_data <- read.csv("hub_genes_before_ML.csv", 
                      header = TRUE,
                      stringsAsFactors = FALSE)  # Ensure that gene names are stored as character strings

#Extract genes by column position (if genes are in the first column)
genes <- gene_data[, 1]  
g <- genes[!is.na(genes) & genes != ""]  # Remove empty values and NA entries
g <- unique(g)  # Remove duplicate genes
# Important: check whether the genes are present in the matrix
g <- intersect(g, rownames(exp))
if(length(g) == 0) stop("Error: none of the genes in the CSV file are present in the expression matrix")

#exp = log2(edgeR::cpm(exp)+1)
exp = exp[g,]
exp = t(exp)

# Add the prefix "GSE151052_" to the row names
rownames(exp) <- paste0("GSE11014732537_", rownames(exp))
# Inspect the modified data
head(rownames(exp))
# Append the Group column from clinical to the row names of exp
rownames(exp) <- paste0(rownames(exp), "_", clinical$Group)
# Inspect the modified data
head(rownames(exp))
# Convert the exp matrix to a data frame
exp <- as.data.frame(exp)
# Create the Group column and populate it based on the row names
exp$Group <- ifelse(grepl("IPF", rownames(exp)), 1, 0)

# Save exp_df as a CSV file
# Export to a CSV file
write.csv(exp, file = "train_data.csv", row.names = TRUE)









#Extract the GSE76925 and GSE47460 datasets and prepare the test set
#GSE76925
#Extract the GSE151052 dataset and prepare the training data
rm(list = ls())
load("/home/hzy-sc1/AAAAA_bioinformatics_workflow/high_throughput_sequencing-pipeline/GEO_public_data/external_validation_processed_RNAseq_GSE150910_IPF_25.7.16/GSE150910_processed.RData")
ls()
exp[1:4,1:4]
table(Group)
gene_data <- read.csv("hub_genes_before_ML.csv", 
                      header = TRUE,
                      stringsAsFactors = FALSE)  # Ensure that gene names are stored as character strings

#Extract genes by column position (if genes are in the first column)
genes <- gene_data[, 1]  
g <- genes[!is.na(genes) & genes != ""]  # Remove empty values and NA entries
g <- unique(g)  # Remove duplicate genes
# Important: check whether the genes are present in the matrix
g <- intersect(g, rownames(exp))
if(length(g) == 0) stop("Error: none of the genes in the CSV file are present in the expression matrix")

#exp = log2(edgeR::cpm(exp)+1)
exp = exp[g,]
exp = t(exp)

# Add the prefix "GSE151052_" to the row names
rownames(exp) <- paste0("GSE150910_", rownames(exp))
# Inspect the modified data
head(rownames(exp))
# Append the Group column from clinical to the row names of exp
rownames(exp) <- paste0(rownames(exp), "_", clinical$Group)
# Inspect the modified data
head(rownames(exp))
# Convert the exp matrix to a data frame
exp <- as.data.frame(exp)
# Create the Group column and populate it based on the row names
exp$Group <- ifelse(grepl("IPF", rownames(exp)), 1, 0)


# 4.3 Finally, replace "_ipf_" or "_control_" with "_GSM" (case-insensitive)
rownames(exp) <- gsub("_(ipf|control)_", "_GSM", rownames(exp), ignore.case = TRUE)
# Save exp_df as a CSV file
# Export to a CSV file
write.csv(exp, file = "test1_data.csv", row.names = TRUE)








#Extract the GSE76925 and GSE47460 datasets and prepare the test set
#GSE47460
#Extract the GSE151052 dataset and prepare the training data
rm(list = ls())
load("/home/hzy-sc1/AAAAA_bioinformatics_workflow/high_throughput_sequencing-pipeline/GEO_public_data/external_validation_processed_RNAseq_GSE213001_IPF_25.7.16/GSE213001_processed.RData")
ls()
exp[1:4,1:4]
library(tinyarray)
# library(edgeR)
# exp = log2(cpm(exp)+1)
table(Group)
gene_data <- read.csv("hub_genes_before_ML.csv", 
                      header = TRUE,
                      stringsAsFactors = FALSE)  # Ensure that gene names are stored as character strings

#Extract genes by column position (if genes are in the first column)
genes <- gene_data[, 1]  
g <- genes[!is.na(genes) & genes != ""]  # Remove empty values and NA entries
g <- unique(g)  # Remove duplicate genes
# Important: check whether the genes are present in the matrix
g <- intersect(g, rownames(exp))
if(length(g) == 0) stop("Error: none of the genes in the CSV file are present in the expression matrix")

exp = exp[g,]
exp = t(exp)

# Add the prefix "GSE151052_" to the row names
rownames(exp) <- paste0("GSE213001_", rownames(exp))
# Inspect the modified data
head(rownames(exp))
# Append the Group column from clinical to the row names of exp
rownames(exp) <- paste0(rownames(exp), "_", clinical$Group)
# Inspect the modified data
head(rownames(exp))
# Convert the exp matrix to a data frame
exp <- as.data.frame(exp)
# Create the Group column and populate it based on the row names
exp$Group <- ifelse(grepl("IPF", rownames(exp)), 1, 0)


# 4.3 Finally, replace "_ipf_" or "_control_" with "_GSM" (case-insensitive)
rownames(exp) <- gsub("_ALF", "_GSM", rownames(exp), ignore.case = TRUE)
# Save exp_df as a CSV file
# Export to a CSV file
write.csv(exp, file = "test2_data.csv", row.names = TRUE)



#Extract the GSE76925 and GSE47460 datasets and prepare the test set
#GSE47460
#Extract the GSE151052 dataset and prepare the training data
rm(list = ls())
load("/home/hzy-sc1/AAAAA_bioinformatics_workflow/high_throughput_sequencing-pipeline/GEO_public_data/external_validation_processed_RNAseq_GSE92592_IPF-25.7.16/GSE92592_processed.RData")
ls()
exp[1:4,1:4]
library(tinyarray)
# library(edgeR)
# exp = log2(cpm(exp)+1)
table(Group)
gene_data <- read.csv("hub_genes_before_ML.csv", 
                      header = TRUE,
                      stringsAsFactors = FALSE)  # Ensure that gene names are stored as character strings

#Extract genes by column position (if genes are in the first column)
genes <- gene_data[, 1]  
g <- genes[!is.na(genes) & genes != ""]  # Remove empty values and NA entries
g <- unique(g)  # Remove duplicate genes
# Important: check whether the genes are present in the matrix
g <- intersect(g, rownames(exp))
if(length(g) == 0) stop("Error: none of the genes in the CSV file are present in the expression matrix")

exp = exp[g,]
exp = t(exp)

# Add the prefix "GSE151052_" to the row names
rownames(exp) <- paste0("GSE92592_", rownames(exp))
# Inspect the modified data
head(rownames(exp))
# Append the Group column from clinical to the row names of exp
rownames(exp) <- paste0(rownames(exp), "_", clinical$Group)
# Inspect the modified data
head(rownames(exp))
# Convert the exp matrix to a data frame
exp <- as.data.frame(exp)
# Create the Group column and populate it based on the row names
exp$Group <- ifelse(grepl("IPF", rownames(exp)), 1, 0)

# Save exp_df as a CSV file
# Export to a CSV file
write.csv(exp, file = "test3_data.csv", row.names = TRUE)





rm(list = ls())
load("/home/hzy-sc1/AAAAA_bioinformatics_workflow/microarray-pipeline/GEO_public_data/external_validation_processed_GSE53845_IPF_human_lung_tissue/GSE53845_processed.RData")
ls()
exp[1:4,1:4]
library(tinyarray)
# library(edgeR)
# exp = log2(cpm(exp)+1)
table(Group)
gene_data <- read.csv("hub_genes_before_ML.csv", 
                      header = TRUE,
                      stringsAsFactors = FALSE)  # Ensure that gene names are stored as character strings

#Extract genes by column position (if genes are in the first column)
genes <- gene_data[, 1]  
g <- genes[!is.na(genes) & genes != ""]  # Remove empty values and NA entries
g <- unique(g)  # Remove duplicate genes
# Important: check whether the genes are present in the matrix
g <- intersect(g, rownames(exp))
if(length(g) == 0) stop("Error: none of the genes in the CSV file are present in the expression matrix")

exp = exp[g,]
exp = t(exp)

# Add the prefix "GSE151052_" to the row names
rownames(exp) <- paste0("GSE53845_", rownames(exp))
# Inspect the modified data
head(rownames(exp))
# Append the Group column from clinical to the row names of exp
rownames(exp) <- paste0(rownames(exp), "_", clinical$Group)
# Inspect the modified data
head(rownames(exp))
# Convert the exp matrix to a data frame
exp <- as.data.frame(exp)
# Create the Group column and populate it based on the row names
exp$Group <- ifelse(grepl("IPF", rownames(exp)), 1, 0)

# Save exp_df as a CSV file
# Export to a CSV file
write.csv(exp, file = "test4_data.csv", row.names = TRUE)






train <- read.table("train_data.csv", header = TRUE, row.names = 1, sep = ",")
test1 <- read.table("test1_data.csv", header = TRUE, row.names = 1, sep = ",")
test2 <- read.table("test2_data.csv", header = TRUE, row.names = 1, sep = ",")
test3 <- read.table("test3_data.csv", header = TRUE, row.names = 1, sep = ",")
test4 <- read.table("test4_data.csv", header = TRUE, row.names = 1, sep = ",")


# Inspect the column names of test1, test2, and train
colnames(test1)
colnames(test2)
colnames(test3)
colnames(test4)
colnames(train)

# Remove extra columns from test2, assuming the column is named "extra_column" (replace as appropriate)
test2 <- test2[, colnames(test1)]

# Remove extra columns from test2, assuming the column is named "extra_column" (replace as appropriate)
test3 <- test3[, colnames(test1)]

# Remove extra columns from test2, assuming the column is named "extra_column" (replace as appropriate)
test4 <- test4[, colnames(test1)]

# Remove extra columns from train, assuming the column is named "extra_column" (replace as appropriate)
train <- train[, colnames(test1)]

# Merge test1 and test2
test <- rbind(test1, test2, test3, test4)

# Inspect the merged data frame
head(test)



write.csv(test, file = "test_merge_data.csv", row.names = TRUE)
write.csv(train, file = "train_merge_data.csv", row.names = TRUE)





###The training and test sets have now been prepared. Proceed to machine learning.
#Load packages
library(openxlsx)
library(seqinr)
library(plyr)
library(randomForestSRC)
library(glmnet)
library(plsRglm)
library(gbm)
library(caret)
library(mboost)
library(e1071)
library(BART)
library(MASS)
library(snowfall)
library(xgboost)
library(ComplexHeatmap)
library(RColorBrewer)
library(pROC)

# remove.packages("randomForestSRC")
# # Use devtools to install a specific version of randomForestSRC from GitHub
# install.packages("devtools")
#devtools::install_version("randomForestSRC", version = "3.3.3", repos = "http://cran.us.r-project.org")

#Set the working directory

source("refer.ML.R")

#Read the training data file
train_data <- read.table("train_merge_data.csv", header = TRUE, row.names = 1, sep = ",")
# Separate features and labels
train_features <- train_data[, -ncol(train_data), drop = FALSE]
train_labels <- train_data[, ncol(train_data),drop = F]


# Read the test data
test_data <- read.table("test_merge_data.csv", header = TRUE, row.names = 1, sep = ",")
test_features <- test_data[, -ncol(test_data), drop = FALSE]
test_labels <- test_data[, ncol(test_data),drop = F]
test_labels$Cohort <- gsub("(.+)\\_(.+)\\_(.+)",'\\1',rownames(test_data))



#Identify genes shared by the training and validation sets
# Ensure that the training and test sets contain the same features
common_features <- intersect(colnames(train_features), colnames(test_features))
train_data <- as.matrix(train_data[, common_features])
test_data <- as.matrix(test_data[, common_features])
train_data = scaleData(train_data, centerFlags=T, scaleFlags=T) 
test_data = scaleData(test_data, cohort = test_labels$Cohort, centerFlags = T, scaleFlags = T)


# ==== Diabetes machine learning analysis workflow ====
#Read the file containing machine learning methods
methods <- read.table("113_ML_methods.txt", header=T, sep="\t", check.names=F)
methods <- methods$x

#Prepare parameters for the machine learning models
classVar = "Group"         #Set the name of the classification variable
min.selected.var = 4      #Threshold for the number of genes
Variable = colnames(train_features)
preTrain.method =  c("Lasso","glmBoost","RF","Stepglm[both]","Stepglm[backward]")
######################Build machine learning models using the training data######################
#Select variables using the first machine learning method in each model combination
preTrain.var <- list()       #Store variables selected by each algorithm
set.seed(seed = 123)         #Set the random seed
# Initialize a list to store timing results
time.list <- list()
for (method in preTrain.method) {
  # Measure the runtime of the current method
  time.taken <- system.time({
    preTrain.var[[method]] <- RunML(
      method = method,
      Train_set = train_data,
      Train_label = train_labels,
      mode = "Variable",
      classVar = classVar
    )
  })
  
  # Store the timing result (elapsed time)
  time.list[[method]] <- time.taken[["elapsed"]]  # Unit: seconds
  
  # Print the runtime of the current method
  cat(sprintf("Method [%s] took %.2f seconds\n", method, time.taken[["elapsed"]]))
}
preTrain.var[["simple"]] <- colnames(train_data)

# Install and load the UpSetR package
# install.packages("UpSetR")
library(UpSetR)

# Extract the gene lists selected by each algorithm
gene_lists <- list(
  Lasso = preTrain.var$Lasso,
  glmBoost = preTrain.var$glmBoost,
  RF = preTrain.var$RF,
  Step_both = preTrain.var$`Stepglm[both]`,
  Step_back = preTrain.var$`Stepglm[backward]`
)

# Create the intersection matrix
upset_data <- fromList(gene_lists)

allcolour=c('#DF0A1F','#1C5BA7','#019E73','#ED621B','#E477C1')


# Plot the UpSet diagram
pdf("gene_intersection_upset.pdf", width=10, height=6)
upset(upset_data, sets = names(gene_lists), 
      order.by = "freq", 
      text.scale = 1.2,
      matrix.color = allcolour,
      mainbar.y.label = "Gene number intersected",
      sets.x.label = "Gene number selected")
dev.off()

# Identify the core genes shared by all algorithms
core_genes <- Reduce(intersect, gene_lists)
write.table(core_genes, "core_genes.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)
## Update the threshold accordingly to prevent unintended removal later
##### 0) Prepare the following two lines before the loop ---------------------------------------------
Variable_core   <- core_genes           # 1) Features used consistently across models
min.selected.var <- length(core_genes)  # Update the threshold (=4)
#--------------------------------------------------------------------------

# Build models using the second machine learning method in each model combination
model <- list()                     # Initialize the results list
set.seed(123)                       # Set a fixed random seed
Train_set_bk <- train_data          # Back up the original training matrix

for (method in methods) {
  cat(match(method, methods), ":", method, "\n")
  method_name <- method
  method      <- strsplit(method, "\\+")[[1]]
  if (length(method) == 1) method <- c("simple", method)
  
  ## ====== 1) Change only this line to use the core genes =========
  Variable   <- Variable_core       # 2) Key replacement
  ## ========================================
  
  train_data <- Train_set_bk[, Variable, drop = FALSE]
  Train_label <- train_labels
  
  model[[method_name]] <- RunML(
    method      = method[2],      # Machine learning method
    Train_set   = train_data,     # Training expression matrix
    Train_label = Train_label,    # Training labels
    mode        = "Model",        # Build the model
    classVar    = classVar)
  
  ## ===== 2) Change the condition to < instead of <= ======
  if (length(ExtractVar(model[[method_name]])) < min.selected.var) {
    model[[method_name]] <- NULL  # 3) Prevent removal of models containing four genes
  }
  ## ========================================
}

# Restore train_data to the original training matrix
train_data <- Train_set_bk; rm(Train_set_bk)

# Save all models
saveRDS(model, "model.MLmodel.rds")


#Calculate classification scores for each sample based on gene expression
model <- readRDS("model.MLmodel.rds")            #Calculate scores using the linear combination function of each machine learning model
methodsValid <- names(model)                     #Identify valid models based on the number of feature genes
#Predict sample risk scores based on gene expression
RS_list <- list()
for (method in methodsValid){
  RS_list[[method]] <- CalPredictScore(fit = model[[method]], new_data = rbind.data.frame(train_data,test_data))
}
riskTab=as.data.frame(t(do.call(rbind, RS_list)))
riskTab=cbind(id=row.names(riskTab), riskTab)
write.table(riskTab, "model.riskMatrix.txt", sep="\t", row.names=F, quote=F)

#Predict sample classes based on gene expression
Class_list <- list()
for (method in methodsValid){
  Class_list[[method]] <- PredictClass(fit = model[[method]], new_data = rbind.data.frame(train_data,test_data))
}
Class_mat <- as.data.frame(t(do.call(rbind, Class_list)))
#Class_mat <- cbind.data.frame(Test_class, Class_mat[rownames(Class_mat),]) # Run this line to merge sample information from the test set
classTab=cbind(id=row.names(Class_mat), Class_mat)
write.table(classTab, "model.classMatrix.txt", sep="\t", row.names=F, quote=F)

#Extract variables selected by each machine learning method
fea_list <- list()
for (method in methodsValid) {
  fea_list[[method]] <- ExtractVar(model[[method]])
}
fea_df <- lapply(model, function(fit){
  data.frame(ExtractVar(fit))
})
fea_df <- do.call(rbind, fea_df)
fea_df$algorithm <- gsub("(.+)\\.(.+$)", "\\1", rownames(fea_df))
colnames(fea_df)[1] <- "features"
write.table(fea_df, file="model.genes.txt", sep = "\t", row.names = F, col.names = T, quote = F)

#Calculate the AUC for each model
AUC_list <- list()
for (method in methodsValid){
  AUC_list[[method]] <- RunEval(fit = model[[method]],      #Machine learning model
                                Test_set = test_data,        #Validation expression data
                                Test_label = test_labels,    #Validation class labels
                                Train_set = train_data,      #Training expression data
                                Train_label = train_labels,  #Training class labels
                                Train_name = "Train",       #Training cohort label
                                cohortVar = "Cohort",       #GEO dataset ID
                                classVar = classVar)        #Classification variable
}
AUC_mat <- do.call(rbind, AUC_list)
aucTab=cbind(Method=row.names(AUC_mat), AUC_mat)
write.table(aucTab, "model.AUCmatrix.txt", sep="\t", row.names=F, quote=F)


##############################Plot the AUC heatmap##############################
#Prepare data for plotting
AUC_mat <- read.table("model.AUCmatrix.txt", header=T, sep="\t", check.names=F, row.names=1, stringsAsFactors=F)

#Rank machine learning models by mean AUC
avg_AUC <- apply(AUC_mat, 1, mean)
avg_AUC <- sort(avg_AUC, decreasing = T)
AUC_mat <- AUC_mat[names(avg_AUC),]
#Identify the best model (highest mean AUC across the training and test sets)
fea_sel <- fea_list[[rownames(AUC_mat)[1]]]
avg_AUC <- as.numeric(format(avg_AUC, digits = 3, nsmall = 3))

#Set colors for heatmap annotations

CohortCol <- brewer.pal(n = ncol(AUC_mat), name = "Set1")
names(CohortCol) <- colnames(AUC_mat)
source("refer.ML.R")

#Create the plot
cellwidth = 1; cellheight = 0.5
hm <- SimpleHeatmap(Cindex_mat = AUC_mat,       #AUC matrix
                    avg_Cindex = avg_AUC,       #Mean AUC
                    CohortCol = CohortCol,      #Colors for the datasets
                    barCol = "#2171b5",       #Color of the bar plot on the right
                    cellwidth = cellwidth, cellheight = cellheight,    #Width and height of each heatmap cell
                    cluster_columns = F, cluster_rows = F)      #Whether to cluster the data

#Export the heatmap
pdf(file="ML_AUC_heatmap.pdf", width=cellwidth * ncol(AUC_mat) + 6, height=cellheight * nrow(AUC_mat) * 0.45)
draw(hm, heatmap_legend_side="right", annotation_legend_side="right")
dev.off()






# Ensure that the required packages have been loaded
library(pROC)
#Read the risk score file
rsFile="model.riskMatrix.txt"      #Risk score matrix file
method="Stepglm[both]+NaiveBayes"               #Select the machine learning method (modify according to the heatmap)
riskRT=read.table(rsFile, header=T, sep="\t", check.names=F, row.names=1)
#Extract dataset IDs
CohortID=gsub("(.*)\\_(.*)\\_(.*)", "\\1", rownames(riskRT))
CohortID=gsub("(.*)\\.(.*)", "\\1", CohortID)
riskRT$Cohort=CohortID

# Ensure that the required packages have been loaded
library(pROC)
# Assume method specifies the column to use; replace it with the actual column name
for (Cohort in unique(riskRT$Cohort)) {
  # Extract sample group information (control and experimental groups)
  rt = riskRT[riskRT$Cohort == Cohort, ]
  # Extract the final component of the row names ("Normal" or another class)
  y = gsub(".*_(.*)$", "\\1", rownames(rt))
  # Convert "Normal" to 0 and all other classes to 1
  y = ifelse(y == "Normal", 0, 1)
  # Check whether y contains two levels (0 and 1)
  if(length(unique(y)) != 2) {
    cat("Error: y must contain exactly two unique values. Check data for Cohort:", Cohort, "\n")
    next
  }
  # Plot the ROC curve
  roc1 = roc(y, as.numeric(rt[, method]))  # Obtain the parameters of the model ROC curve
  ci1 = ci.auc(roc1, method = "bootstrap")  # Obtain the confidence interval for the area under the ROC curve
  ciVec = as.numeric(ci1)
  # Save the results as a PDF
  pdf(file = paste0("ROC.", Cohort, ".pdf"), width = 5, height = 4.75)
  # Plot the ROC curve and display the AUC
  plot(roc1, print.auc = TRUE, col = "#006d2c", legacy.axes = TRUE, main = Cohort)
  # Add the 95% confidence interval to the plot
  text(0.39, 0.43, paste0("95% CI: ", sprintf("%.03f", ciVec[1]), "-", sprintf("%.03f", ciVec[3])), col = "#006d2c")
  # Close the PDF device
  dev.off()
}










#Load packages
library(glmnet)
library(pROC)
expFile="merged_normalized_data.txt"      #Expression data file
geneFile="model.genes.txt"      #Gene list file
#Read the expression data file
rt=read.table("train_merge_data.csv", header = TRUE, row.names = 1, sep = ",")
rt =t(rt)
#Extract sample group information (control and experimental groups)
y = gsub(".*_(.*)$", "\\1", colnames(rt))
# Convert "Normal" to 0 and all other classes to 1
y = ifelse(y == "Normal", 0, 1)
#Read the gene list file
geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)
geneRT <- geneRT$V1[geneRT$V2 == 'Stepglm[both]+NaiveBayes']
#Define plot colors
bioCol=rainbow(length(geneRT), s=0.9, v=0.9)

#Loop through the genes and plot ROC curves
aucText=c()
k=0
for(x in as.vector(geneRT)){
  k=k+1
  #Plot the ROC curve
  roc1=roc(y, as.numeric(rt[x,]))     #Obtain the parameters of the ROC curve
  if(k==1){
    pdf(file="ROC.genes.pdf", width=9, height=9)
    plot(roc1, print.auc=F, col=bioCol[k], legacy.axes=T, main="", lwd=3)
    aucText=c(aucText, paste0(x,", AUC=",sprintf("%.3f",roc1$auc[1])))
  }else{
    plot(roc1, print.auc=F, col=bioCol[k], legacy.axes=T, main="", lwd=3, add=TRUE)
    aucText=c(aucText, paste0(x,", AUC=",sprintf("%.3f",roc1$auc[1])))
  }
}
#Draw the legend and display the area under the ROC curve
legend("bottomright", aucText, lwd=3, bty="n", 
       cex=0.8, col=bioCol[1:(ncol(rt)-1)],
       inset = c(0.05, 0))  # Shift to the right by 5% of the plot width
dev.off()


