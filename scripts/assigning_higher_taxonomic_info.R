# Get higher taxonomic information
library(stringr)

tax_export <- read.csv("/Users/kbja10/Documents/Cornell/Research/Oneida/Data_analysis/Datasets/oneida_output_8.17.2021/blstn_nonchim_fmt10_scinames_Oneida_withseq_taxonomic.csv", head=FALSE)
tax_export[,1] <- as.character(tax_export[,1])
names(tax_export) <- c("seqid","pident","evalue","qstart","qend","length","sscinames","scomnames","sseq","staxids","qseqid","commnames",
                       "taxinfo","taxclass")
output <- data.frame(seqid=NA,superkingdom=NA,kingdom=NA,phylum=NA,subphylum=NA,superclass=NA,class=NA,subclass=NA,
                     infraclass=NA,cohort=NA,subcohort=NA,superorder=NA,order=NA,suborder=NA,
                     family=NA,subfamily=NA,genus=NA,species=NA)

for (i in 1:nrow(tax_export)){
temp <- tax_export[i,]
df <- data.frame(taxclass=t(do.call('rbind', strsplit(as.character(temp$taxclass),'|',fixed=TRUE))),
                  taxinfo=t(do.call('rbind', strsplit(as.character(temp$taxinfo),'|',fixed=TRUE))))
colnames(df) <- c("class","value")
output[i,1] <- temp$seqid
if ("superkingdom" %in% df$class) {output[i,2] <- as.character(df[df$class=="superkingdom",2])} else output[i,2] <- NA
if ("kingdom" %in% df$class) {output[i,3] <- as.character(df[df$class=="kingdom",2])} else output[i,3] <- NA
if ("phylum" %in% df$class) {output[i,4] <- as.character(df[df$class=="phylum",2])} else output[i,4] <- NA
if ("subphylum" %in% df$class) {output[i,5] <- as.character(df[df$class=="subphylum",2])} else output[i,5] <- NA
if ("superclass" %in% df$class) {output[i,6] <- as.character(df[df$class=="superclass",2])} else output[i,6] <- NA
if ("class" %in% df$class) {output[i,7] <- as.character(df[df$class=="class",2])} else output[i,7] <- NA
if ("infraclass" %in% df$class) {output[i,8] <- as.character(df[df$class=="infraclass",2])} else output[i,8] <- NA
if ("subclass" %in% df$class) {output[i,9] <- as.character(df[df$class=="subclass",2])} else output[i,9] <- NA
if ("cohort" %in% df$class) {output[i,10] <- as.character(df[df$class=="cohort",2])} else output[i,10] <- NA
if ("subcohort" %in% df$class) {output[i,11] <- as.character(df[df$class=="subcohort",2])} else output[i,11] <- NA
if ("superorder" %in% df$class) {output[i,12] <- as.character(df[df$class=="superorder",2])} else output[i,12] <- NA
if ("order" %in% df$class) {output[i,13] <- as.character(df[df$class=="order",2])} else output[i,13] <- NA
if ("suborder" %in% df$class) {output[i,14] <- as.character(df[df$class=="suborder",2])} else output[i,14] <- NA
if ("family" %in% df$class) {output[i,15] <- as.character(df[df$class=="family",2])} else output[i,15] <- NA
if ("subfamily" %in% df$class) {output[i,16] <- as.character(df[df$class=="subfamily",2])} else output[i,16] <- NA
if ("genus" %in% df$class) {output[i,17] <- as.character(df[df$class=="genus",2])} else output[i,17] <- NA
if ("species" %in% df$class) {output[i,18] <- as.character(df[df$class=="species",2])} else output[i,18] <- NA
}
output <- cbind(tax_export[,1:10], output)
# write.csv(output, file="/Users/kbja10/Documents/Cornell/Research/Oneida/Data_analysis/Datasets/oneida_output_8.17.2021/blstn_nonchim_fmt10_scinames_Oneida_withseq_taxonomic_matched.csv")
