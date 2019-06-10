library(feather)

samp <- scan(file = "data/bam_logFile.txt", "character")

anno <- read_feather("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/human_20190408_Cross_areal/anno.feather")

missing <- setdiff(anno$sample_name.1_label, samp)
extra <- setdiff(samp, anno$sample_name.1_label)

table(anno$roi_label[anno$sample_name.1_label %in% missing])

write.table(missing, file = "output/missing_samp.txt", 
            col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(extra, file = "output/extra_samp.txt", 
            col.names = FALSE, row.names = FALSE, quote = FALSE)
