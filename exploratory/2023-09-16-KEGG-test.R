if (!requireNamespace("BiocManager", quietly=TRUE))
      + install.packages("BiocManager")
BiocManager::install("KEGGREST", force = T)

library(KEGGREST)
listDatabases()
ko <- keggList("brite")
keggGet("K01097")[[1]]$BRITE[2:5] %>% as.vector()


genome = read_tsv('data/test/299R.tsv', skip=5, show_col_types = FALSE) %>%
      select(locus_tag = `Locus Tag`, gene = Gene, product = Product, ref = DbXrefs)

file_list = list.files('data/test', pattern=".tsv", full.names = TRUE)

g_list = lapply(file_list, read_tsv, skip=5, show_col_types = FALSE)


genome %>% 
      extract(ref, into="type", regex = "(COG:COG[0-9]+|KEGG:K[0-9]+)") %>% 
      na.omit
kegg_genome = genome %>% 
      extract(ref, into="type", regex = "(KEGG:K[0-9]+)") %>% 
      na.omit
cog_genome = genome %>% 
      extract(ref, into="type", regex = "(COG:COG[0-9]+)") %>% 
      na.omit


kegg_list <- kegg_genome %>% 
      select(type) %>% 
      extract(type, 'kegg', regex="(K[0-9]+)") %>% 
      as.data.frame()

# translates KEGG classifier to COG
keggGet(kegg_list$kegg[25])[[1]]$DBLINKS %>% str_extract(pattern = "(COG[0-9]+)") %>% na.omit %>% as.vector

keg_cog <- function(x){
      if (is.null(keggGet(kegg_list$kegg[x])[[1]]$DBLINKS) == FALSE) {
            keggGet(kegg_list$kegg[x])[[1]]$DBLINKS %>% 
                  str_extract(pattern = "(COG[0-9]+)") %>% 
                  na.omit %>% 
                  as.vector
      } else {
            print("NA")
      }
}

cog_list <- vector()

for(i in 1:nrow(kegg_list)){
      cog_list[i] <- keg_cog(i)
}


keg_cog(25)
is.null(keggGet(kegg_list$kegg[74])[[1]]$DBLINKS)

cog_list
