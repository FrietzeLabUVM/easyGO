
sel_go = search_GO_by_term("activation of MAPK")
sel_go.all = get_all_children(sel_go)

gene_name2go = get_gene_name2go(org.db = org.Hs.eg.db::org.Hs.eg.db)

get_GO_gene_names(sel_go.all, gene_name2go = gene_name2go, as.list = TRUE)

library(org.Hs.eg.db)
entrez2go = select(org.Hs.eg.db, keys= keys(org.Hs.eg.db), columns = c("GO"))
entrez2gene_name = select(org.Hs.eg.db, keys= keys(org.Hs.eg.db), columns = c("SYMBOL"))
head(entrez2go)
head(entrez2gene_name)



seqsetvis::ssvFeatureVenn(
  list(
    parents = subset(gene_name2go, GO %in% sel_go)$SYMBOL,
    children = subset(gene_name2go, GO %in% setdiff(sel_go.all, sel_go))$SYMBOL

  )
)



select(GO.db, sel_go, tc)
# select(GO.db, sel_go.children, tc)
select(GO.db, sel_go.all, tc)




library(annotate)

keys(org.Hs.egGENENAME)
mappedkeys(org.Hs.egGENENAME)

keys(org.Hs.egGO)
columns(org.Hs.egGO)
