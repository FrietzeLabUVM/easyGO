#' Title
#'
#' @param go_ids
#' @param include_children
#' @param include_parents
#' @param includes_genes
#' @param org.db
#' @param gene_name2go
#'
#' @return
#' @export
#'
#' @examples
get_GO_info = function(go_ids, include_children = FALSE, include_parents = FALSE, includes_genes = FALSE, org.db = NULL, gene_name2go = NULL){
    go_df = as.data.frame(select(GO.db, go_ids, columns(GO.db)))
    if(include_children){
        all_children = lapply(go_ids, get_all_children.no_parents)
        go_df$CHILDREN = sapply(all_children, paste, collapse = ",")

    }
    if(include_parents){
        all_parents = lapply(go_ids, get_parents)
        go_df$PARENTS = sapply(all_parents, paste, collapse = ",")
    }
    if(includes_genes){
        if(!is.null(org.db)){
            message(paste(
                sep = "\n  ",
                "org.db must be provided. Examples:",
                "org.Hs.eg.db",
                "org.Mm.eg.db",
                "org.Dm.eg.db"
            ))
        }
        if(is.null(gene_name2go)){
            message("making gene_name2go.  Supply result of get_gene_name2go(org.db) to save time.")
            gene_name2go = get_gene_name2go(org.db)
        }
        all_children = lapply(go_ids, get_all_children)
        all_genes = lapply(all_children, get_GO_gene_names, gene_name2go = gene_name2go)
        go_df$GENES = sapply(all_genes, paste, collapse = ",")
    }
    go_df
}

#' Title
#'
#' @param search_term
#'
#' @return
#' @export
#'
#' @examples
search_GO_by_term = function(search_term){
    sel = subset(select(GO.db, keys(GO.db), columns(GO.db)), grepl(search_term, TERM))
    sel_go = sel$GOID
    sel_go
}

#' Title
#'
#' @param children_go
#'
#' @return
#' @export
#'
#' @examples
get_parents = function(children_go){
    parent_result = annotate::getGOParents(children_go)
    parent_go = lapply(parent_result, function(x)x$Parents)
    parent_go = lapply(parent_go, function(x)x[names(x) == "is_a"])
    found_go = unlist(parent_go)
    names(found_go) = NULL
    found_go
}

#' Title
#'
#' @param parent_go
#' @param ko_go
#' @param depth
#' @param as.list
#'
#' @return
#' @export
#'
#' @examples
get_all_children = function(parent_go, ko_go = character(), depth = 1, as.list = FALSE){
    if(as.list){
        names(parent_go) = parent_go
        lapply(parent_go, function(go){
            get_all_children(go, as.list = FALSE)
        })
    }else{
        child_result = annotate::getGOChildren(parent_go)
        children_go = unlist(lapply(child_result, function(x)x$Children))
        found_go = union(parent_go, ko_go)
        if(is.null(children_go)){
            return(found_go)
        }
        children_go.new = setdiff(children_go, found_go)
        get_all_children(children_go.new, found_go, depth = depth + 1)
    }
}

#' Title
#'
#' @param parent_go
#' @param as.list
#'
#' @return
#' @export
#'
#' @examples
get_all_children.no_parents = function(parent_go, as.list = FALSE){
    if(as.list){
        names(parent_go) = parent_go
        lapply(parent_go, function(go){
            get_all_children.no_parents(go, as.list = FALSE)
        })
    }else{
        all_go = get_all_children(parent_go)
        setdiff(all_go, parent_go)
    }
}

#' Title
#'
#' @param org.db
#'
#' @return
#' @export
#'
#' @examples
get_gene_name2go = function(org.db = org.Hs.eg.db::org.Hs.eg.db){
    entrez2go = select(org.db, keys= keys(org.db), columns = c("GO"))
    entrez2gene_name = select(org.db, keys= keys(org.db), columns = c("SYMBOL"))
    gene_name2go = merge(entrez2go, entrez2gene_name, by = 'ENTREZID')
    gene_name2go
}

#' Title
#'
#' @param go_ids
#' @param org.db
#' @param as.list
#' @param gene_name2go
#'
#' @return
#' @export
#'
#' @examples
get_GO_gene_names = function(go_ids, org.db = org.Hs.eg.db::org.Hs.eg.db, as.list = FALSE, gene_name2go = NULL){
    if(is.null(gene_name2go)){
        message("making gene_name2go.  Supply result of get_gene_name2go(org.db) to save time.")
        gene_name2go = get_gene_name2go(org.db)
    }


    go_df = as.data.frame(subset(gene_name2go, GO %in% go_ids))
    if(as.list){
        split(go_df$SYMBOL, go_df$GO)
    }else{
        unique(go_df$SYMBOL)
    }

}


