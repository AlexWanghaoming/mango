# This script is used to parse the results of interproscan

library(GO.db)
library(tidyverse)
# write GO base infomation
godb <- AnnotationDbi::select(GO.db, keys(GO.db), columns(GO.db))
write_tsv(godb, "/Users/alexwang/0data/db/godb.txt")

interpro <- read.table("/Users/alexwang/0data/xulab_current/PH1-fusion.interproscan.tsv", fill = T, sep = "\t", 
                       col.names = paste0("V", 1:15), na.strings = "")

# write gene2GO 
gene2go <- interpro %>% dplyr::select(Gene = V1, GOID = V14) %>%
  mutate(Gene = str_replace(Gene, "\\.\\d+$", "")) %>% na.omit() %>%
  separate(GOID, paste0("X", 1:(max(str_count(.$GOID,"\\|"))+1), seq = ""), sep = "\\|") %>%
  gather(key = "X", value = "GOID", -Gene) %>% dplyr::select(Gene, GOID) %>%
  na.omit() %>% base::unique()

write_tsv(gene2go, "/Users/alexwang/0data/xulab_current/go/gene2go.txt")
go_anno <- gene2go %>% left_join(godb)
write_tsv(go_anno, "/Users/alexwang/0data/xulab_current/go/go_annot.txt")
