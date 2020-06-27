setwd("/var/www/html/server_corins/scripts_corins/")
packrat::init()

# arguments ####
args = commandArgs(trailingOnly=TRUE)
directory <- args[1]
print(directory)
print("oi")

version <- paste0(R.Version()$major, ".", unlist(strsplit(R.Version()$minor, ".", fixed = TRUE))[1])

library(igraph)
library(snow)
library(data.table)
library(parallel)
# functions ####

nodesProperties <- function(dir, nCores = parallel::detectCores()){
  files <- list.files(dir)
  idx <- grepl("edges.txt$", files) # estava dando erro por causa da expressao regular
  files <- files[idx]
  rm(idx)
  cl <- snow::makeSOCKcluster(nCores)
  on.exit(snow::stopCluster(cl))
  propertiesList <- snow::parLapply(cl, files, function(i){
    # undirected edge list
    fileContent <- data.table::fread(paste0(dir, "/", i), stringsAsFactors = FALSE)
    fileContent <- fileContent[, c(1, 3, 4)]
    w <- fileContent$Distance
    if(all(w > 0)){
      g <- igraph::graph.data.frame(d = fileContent[, c(1, 2)], directed = FALSE)
      fileContent <- fileContent[, c(1, 2)]
      colnames(fileContent) <- c("n1", "n2")
      aux <- data.frame(n1 = fileContent$n2, n2 = fileContent$n1)
      # directed edge list
      fileContent <- rbind(fileContent, aux)
      rm(aux)
      df <- table(fileContent$n1)
      df <- as.data.frame(df, stringsAsFactors = FALSE)
      colnames(df) <- c("node", "degree")
      df$aminoAcid <- gsub("^.*:_:(.*)$", "\\1", df$node)
      df$triangles <- vapply(df$node, function(x){
        as.integer(igraph::count_triangles(g, vids = x))}, integer(1))
      df$clusteringCoef <- mapply(function(x, y){
        if (y < 2) {
          return(0)
        } else {
          return(x/(y * (y - 1)/2))
        }
      }, df$triangles, df$degree)
      igraph::E(g)$weight <- w
      betweennessWeighted <- igraph::betweenness(g, v = df$node, directed = FALSE,
                                                 nobigint = TRUE,
                                                 normalized = TRUE)
      df$betweennessWeighted <- sapply(df$node, function(i){
        return (unname(betweennessWeighted[i]))
      })
      df$filename <- i
      return (df)
    }else{
      return (NULL)
    }
  })
  properties <- data.frame()
  invisible(sapply(seq.int(length(files)), function(i){
    if(!is.null(propertiesList[[i]])){
      properties <<- rbind(properties, propertiesList[[i]])
    }
    return (NULL)
  }))
  return (properties)
}

graphsMeans <- function(dir, file, aux){
  fileContent <- data.table::fread(paste0(dir, "/", file),
                                   stringsAsFactors = FALSE)
  fileContent <- fileContent[, c(1, 3)]
  g <- igraph::graph.data.frame(d = fileContent, directed = FALSE)
  df <- data.frame(degree = mean(aux$degree),
                   clusteringCoef = mean(aux$clusteringCoef),
                   betweennessWeighted = mean(aux$betweennessWeighted))
  df$graphAssortativity <- igraph::assortativity_degree(g, directed = FALSE)
  df$filename <- file
  return (df)
}

# output ####

nodesResult <- nodesProperties(directory)
files <- list.files(directory)
idx <- grepl("edges.txt$", files) # estava dando erro por causa da expressao regular
files <- files[idx]
rm(idx)
graphsResult <- data.frame()
for(i in files){
  if(i %in% nodesResult$filename){
    aux <- nodesResult[nodesResult$filename == i,]
    temp <- graphsMeans(directory, i, aux)
    graphsResult <- rbind(graphsResult, temp)
  }
}
rm(temp, aux, files, i)

files <- unique(graphsResult$filename)
lista <- lapply(files, function(i){
  nodesResult[nodesResult$filename == i, -7]
})

setwd(directory) # a pasta PycharmProjects n�o existe
for(i in seq.int(1,length(files))){
  write.table(x = lista[[i]], file = paste0("ResultadoNos_", files[i]), row.names = F, quote = F, sep = "\t")
}

lista <- lapply(files, function(i){
  graphsResult[graphsResult$filename == i, -5]
})
for(i in seq.int(1,length(files))){
  write.table(x = lista[[i]], file = paste0("ResultadoGrafo_", files[i]), row.names = F, quote = F, sep = "\t")
}

#save(nodesResult, file = "nodesResult.RData", compress = "xz")
#save(graphsResult, file = "graphsResult.RData", compress = "xz")
