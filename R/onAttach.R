.onAttach <- function(libname, pkgname) {
  packageStartupMessage(paste(
    sep = "\n  ",
    "Attaching easyGO. Use with an AnnotationDbi such as:",
    "org.Hs.eg.db::org.Hs.eg.db",
    "org.Mm.eg.db::org.Mm.eg.db",
    "org.Dm.eg.db::org.Dm.eg.db"
  ))
}
