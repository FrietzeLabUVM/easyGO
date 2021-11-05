#' get_GO_info
#'
#' @param go_ids character vector of GO ids, like GO:1234.
#' @param include_children If TRUE, a column of CHILDREN is added. Values are
#'   comma delimitted.
#' @param include_parents If TRUE, a column of PARENTS is added. Values are
#'   comma delimitted.
#' @param include_genes If TRUE, a column of GENES is added. Values are comma
#'   delimitted. OrgDb must be supplied.
#' @param org.db OrgDb such as org.Hs.eg.db, org.Mm.eg.db, org.Dm.eg.db etc.
#'   Required for include_genes.
#' @param gene_name2go Optional result of get_gene_name2go() for org.db. Much
#'   faster if supplied.
#'
#' @return A data.frame of information about supplied GO ids.
#' @export
#'
#' @examples
#' go_ids = head(search_GO_by_term("MAPK"))
#' get_GO_info(go_ids)
#' # Additional columns can be added though these take a little longer.
#' get_GO_info(go_ids,
#'   include_children = TRUE,
#'   include_parents = TRUE,
#'   include_genes = TRUE,
#'   org.db = org.Hs.eg.db::org.Hs.eg.db)
#' # Speed things up by pre-generating gene_name2go and reusing.
#' gene_name2go = get_gene_name2go(org.Hs.eg.db::org.Hs.eg.db)
#' get_GO_info(go_ids,
#'   include_children = TRUE,
#'   include_parents = TRUE,
#'   include_genes = TRUE,
#'   org.db = org.Hs.eg.db::org.Hs.eg.db,
#'   gene_name2go = gene_name2go)
get_GO_info = function(go_ids, include_children = FALSE, include_parents = FALSE, include_genes = FALSE, org.db = NULL, gene_name2go = NULL){
    go_df = as.data.frame(select(GO.db, go_ids, columns(GO.db)))
    if(include_children){
        all_children = lapply(go_ids, get_all_children.no_parents)
        go_df$CHILDREN = sapply(all_children, paste, collapse = ",")

    }
    if(include_parents){
        all_parents = lapply(go_ids, get_parents)
        go_df$PARENTS = sapply(all_parents, paste, collapse = ",")
    }
    if(include_genes){
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

#' search_GO_by_term
#'
#' @param search_term character term to search for in TERM (or DEFINITION) of GO database.
#' @param search_field character field to search in GO database.
#'
#' @return character vector of GO ids containing search_term.
#' @export
#' @import AnnotationDbi
#'
#' @examples
#' search_GO_by_term("MAPK")
#' search_GO_by_term("MAPK", search_field = "DEFINITION")
search_GO_by_term = function(search_term, search_field = c("TERM", "DEFINITION")[1]){
    stopifnot(search_field %in% c("TERM", "DEFINITION"))
    sel = subset(AnnotationDbi::select(GO.db, keys(GO.db), columns(GO.db)), grepl(search_term, get(search_field)))
    sel_go = sel$GOID
    sel_go
}

#' get_parents
#'
#' @param children_go GO to find parents for.
#'
#' @return GO ids of parent terms.
#' @export
#' @import annotate
#'
#' @examples
#' go_ids = head(search_GO_by_term("MAPK"))
#' get_parents(go_ids)
get_parents = function(children_go){
    parent_result = annotate::getGOParents(children_go)
    parent_go = lapply(parent_result, function(x)x$Parents)
    parent_go = lapply(parent_go, function(x)x[names(x) == "is_a"])
    found_go = unlist(parent_go)
    names(found_go) = NULL
    found_go
}

#' get_all_children
#'
#' @param parent_go GO to find children for.
#' @param ko_go These GO terms can't be returned as children.  Mainly used recursively by function and should not be supplied by user.
#' @param depth Current depth of recursion.  Mainly used recursively by function and should not be supplied by user.
#' @param as.list If TRUE, a list with children per input parent_go is returned instead of unique character vector for all parent_go.
#'
#' @return Either a character vector of children (includes parents). If as.list == TRUE, a list per parent_go.
#' @export
#' @import annotate
#'
#' @examples
#' go_ids = head(search_GO_by_term("MAPK"))
#' parent_go = get_parents(go_ids)
#' get_all_children(parent_go)
#' get_all_children(parent_go, as.list = TRUE)
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

#' get_all_children.no_parents
#'
#' @param parent_go GO to find children for.
#' @param as.list If TRUE, a list with children per input parent_go is returned instead of unique character vector for all parent_go.
#'
#' @return Either a character vector of children (includes parents). If as.list == TRUE, a list per parent_go.
#' @export
#'
#' @examples
#' go_ids = head(search_GO_by_term("MAPK"))
#' parent_go = get_parents(go_ids)
#' get_all_children.no_parents(parent_go)
#' get_all_children.no_parents(parent_go, as.list = TRUE)
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

#' get_gene_name2go
#'
#' This function is used internally by other functions to create a data.frame mapping GO terms to gene_names.  May be called to pre-generate this data.frame and supplied to later functions to speed things up greatly.
#'
#' @param org.db OrgDb such as org.Hs.eg.db, org.Mm.eg.db, org.Dm.eg.db etc.
#'
#' @return data.frame mapping GO to SYMBOL
#' @export
#' @import org.Hs.eg.db
#'
#' @examples
#' gene_name2go = get_gene_name2go(org.Hs.eg.db::org.Hs.eg.db)
get_gene_name2go = function(org.db = org.Hs.eg.db::org.Hs.eg.db){
    entrez2go = select(org.db, keys= keys(org.db), columns = c("GO"))
    entrez2gene_name = select(org.db, keys= keys(org.db), columns = c("SYMBOL"))
    gene_name2go = merge(entrez2go, entrez2gene_name, by = 'ENTREZID')
    gene_name2go
}

#' get_GO_gene_names
#'
#' @param go_ids character vector of GO ids, like GO:1234.
#' @param org.db OrgDb such as org.Hs.eg.db, org.Mm.eg.db, org.Dm.eg.db etc.
#' @param as.list If TRUE, a list with children per input parent_go is returned instead of unique character vector for all parent_go.
#' @param gene_name2go Optional result of get_gene_name2go() for org.db. Much
#'   faster if supplied.
#'
#' @return data.frame with mapping for ENTREZID, GO, EVIDENCE, ONTOLOGY, SYMBOL
#' @export
#' @import org.Hs.eg.db
#'
#' @examples
#' gene_name2go = get_gene_name2go(org.Hs.eg.db::org.Hs.eg.db)
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


