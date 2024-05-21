setwd("E:\\药物联用项目\\操作\\相似性计算\\蛋白GO\\go_simility_提取/")
drug = read.table("558drug.txt",sep="\t",header=T,check.names=F)

drug_go_matrix <- matrix(data = 0, nrow = length(unique(drug$id)), ncol = length(unique(drug$id)))
rownames(drug_go_matrix) <- unique(drug$id)
colnames(drug_go_matrix) <- unique(drug$id)

go_simi <- read.table("go_similarity.txt",sep="\t",header=T,check.names=F,row.names = 1)

drug_target = read.table("drug_target.txt",sep="\t",header=T,check.names=F)

drug_go_matrix <- as.data.frame(drug_go_matrix)
for (i in rownames(drug_go_matrix[320:nrow(drug_go_matrix),])) {
  row_extract <- drug_target[which(drug_target$drug_name == i),]
  for (r in colnames(drug_go_matrix)) {
    col_extract <- drug_target[which(drug_target$drug_name == r),]
    llput <- NULL
    for (j in 1:nrow(row_extract)) {
      for (k in 1:nrow(col_extract)) {
        llput <- data_frame(id1 = row_extract[j,"gene_name"],
                            id2 = col_extract[k,"gene_name"],
                            score = go_simi[row_extract[j,"gene_name"],col_extract[k,"gene_name"]])
        llput <- rbind(llput,llput)
      
        
      }
      
    }
    drug_go_matrix[i,r] = max(llput$score)
  }
}

write.table(drug_go_matrix, "drug_go_matrix1.txt", sep="\t", quote=F, row.names=T)
