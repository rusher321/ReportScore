# BiocManager::install("KEGGREST")

#' dawnload_data
#'
#' @param db_dir
#'
#' @return
#' @export
#'
#' @examples
download_data <- function(db_dir = "database/"){
  library(KEGGREST)
  library(tidyverse)

  # database directory
  db_dir = "database/"
  dir.create(db_dir)

  # pathway/module list
  pwList = keggList("pathway")
  pwList = data.frame(Description=pwList)
  pwList = pwList %>%
    rownames_to_column(var = "PathwayID") %>%
    mutate(PathwayID = gsub("path:","",PathwayID))

  mdList = keggList("module")
  mdList = data.frame(Description=mdList)
  mdList = mdList %>%
    rownames_to_column(var = "ModuleID") %>%
    mutate(ModuleID = gsub("md:","",ModuleID))

  # categories
  GetCategories <- function(dat){
    dat$Class <- NA
    for (x in 1:nrow(dat)){
      if(!is.na(dat$Class[x])){next}
      cat(x,"\n")
      ID <- dat[x,1]
      res <- keggGet(ID)[[1]]
      clas <- ifelse("CLASS" %in% names(res), res$CLASS, NA)
      dat$Class[x] <- clas
    }
    dat
  }
  pwList = GetCategories(pwList)
  mdList = GetCategories(mdList)

  # Get KO list
  GetKOList <- function(dat){
    dat$n <- NA
    dat$KO <- NA
    for (x in 1:nrow(dat)){
      if(!is.na(dat$KO[x])){next}
      cat(x,"\n")
      ID <- dat[x,1]
      res <- keggLink("ko",ID)
      res <- gsub("ko:","",res)
      dat$n[x] <- length(res)
      dat$KO[x] <- paste(res,collapse = ",")
    }
    dat
  }
  pwList = GetKOList(pwList)
  mdList = GetKOList(mdList)

  # filter pathway/module without KOs
  pwList = pwList[pwList$n>0,]
  mdList = mdList[mdList$n>0,]

  write.table(pwList,paste0(db_dir,"pathway.annotation.txt"))
  write.table(mdList,paste0(db_dir,"module.annotation.txt"))

}
