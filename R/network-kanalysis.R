library(igraph)
library(bipartite)
library(ggplot2)

#' Functions for k-core decompose a bipartite graph and computing some newtork indexes

#' Reads a network interaction matrix from a CSV file
#'
#' Args:
#'   namenetwork: CSV file that contains the interaction matrix
#'   guild_a, guild_b: Identifier of the guild of species of each class. Default "pl" (plant)
#'                     "pol" (pollinator)
#'   directory: directory where newtork CSVs are located
#'
#' Return List:
#'   graph: Newtork as an igraph object
#'   m    : Interaction matrix
#'   num_guild_b : number of species of guild_b
#'   num_guild_a" : number of species of guild_a
#'   names_guild_a : names of nodes of guild_a
#'   names_guild_b : names of species of guild_b

read_network <- function(namenetwork, guild_astr = "pl", guild_bstr = "pol", directory="")
{
  # Reading species names
  namesred <- read.csv(paste0(directory,namenetwork),header=FALSE,stringsAsFactors=FALSE)
  names_guild_a <- namesred[1,2:ncol(namesred)]
  names_guild_b <- namesred[2:nrow(namesred),1]

  #Reading matrix data
  m <- read.csv(paste0(directory,namenetwork),header=TRUE,row.names=1)

  num_guild_a <- ncol(m)
  num_guild_b <- nrow(m)
  g <- graph.empty()
  for (i in 1:num_guild_a){
    g <- g + vertices(paste0(guild_astr,i),color="white",guild_id="a",name_species=names_guild_a[i],id=i)
  }
  for (i in 1:num_guild_b){
    g <- g + vertices(paste0(guild_bstr,i),color="red",guild_id="b",name_species=names_guild_b[i],id=i)
  }


  # for (j in 1:num_guild_a)
  # {
  #   for (i in 1:num_guild_b)
  #   {
  #     if (m[i,j]!=0) {
  #       g <- g + igraph::edges(paste0(guild_bstr,i),paste0(guild_astr,j))
  #     }
  #   }
  # }

  mm <- matrix(unlist(list(m)),nrow=num_guild_b,ncol=num_guild_a)
  listedgesn <- which(mm!=0, arr.ind = T)
  listedgesn <- listedgesn[order(listedgesn[,1],listedgesn[,2]),]
  listedgesn[,1] <- paste0(guild_bstr,listedgesn[,1])
  listedgesn[,2] <- paste0(guild_astr,listedgesn[,2])
  g <- g + graph.edgelist(listedgesn)
  calc_values <- list("graph" = g, "matrix" = m, "num_guild_b" = num_guild_b, "num_guild_a" = num_guild_a,
                      "names_guild_a" = names_guild_a, "names_guild_b"=names_guild_b)
  return(calc_values)

}

#' Network analysis using k core decomposition
#'
#' This function performs the kcore decomposition of a bipartite network. The input is
#' the interaction matrix in a .csv file, with the format of  www.web-of-life.es
#' Species of guild a are distributed by colums, and those of guild b by rows. First
#' colum contains the labels of guild b nodes, and first row, the labels of guild a.
#' If the interaction matrix is binary, the cell of species_a_m,species_b_n will be set to 1.
#' If it is weighted, to a real number different of 0.
#'
#' @param namenetwork is the file name that contains the interaction matrix
#' @param directory where the network file is stored
#' @param guild_a prefix for the guild of nodes stored in rows
#' @param guild_b prefix for the guild of nodes stored in columns
#' @param plot_graphs plot kshell histogram and kamada kawai plots
#' @param only_NODF just compute the NODF measurement of nestedness
#' @return \code{calc_values} a list containing the following objects
#' \itemize{
#'  \item{\code{"graph"}}{ an \code{igraph::graph} object}
#'  \item{\code{"max_core"}}{ maximum k shell index}
#'  \item{\code{"nested_values"}}{ a list containing all the values provided by the \code{bipartite::nested} function, unless \code{only_NODF} set \code{TRUE}}
#'  \item{\code{"num_guild_a"}}{ number of nodes of guild a}
#'  \item{\code{"num_guild_b"}}{ number of nodes of guild b}
#'  \item{\code{"links"}}{ number of network links}
#'  \item{\code{"meandist"}}{ network average kradius}
#'  \item{\code{"meankdegree"}}{ network average kdegree}
#'  \item{\code{"spaths_mat"}}{ matrix with node to node shortest distance paths}
#'  \item{\code{"matrix"}}{ adyacency matrix with nodes of guild a by columns and guild b by rows}
#'  \item{\code{"g_cores"}}{ list with the value of kshell for each node}
#'  \item{\code{"modularity_measure"}}{ value of \code{igraph::modularity} function}
#'  }
#' @export
#' @examples result_analysis <- analyze_network("M_PL_007.csv", directory = "data/", guild_a = "Plant", guild_b = "Pollinator")


analyze_network <- function(namenetwork, directory="", guild_a = "pl", guild_b = "pol", plot_graphs = FALSE, only_NODF = FALSE,
                            weight_direction = "none")
{

  calc_kradius <- function(i)
  {
    kradius <- 0
    kradius_core <- mean(spaths_mat[i,][an$guild_maxcore])
    if (!is.na(kradius_core)){
      kradius <- kradius + kradius_core
    }
    V(an$g)[i]$kradius <- kradius
    V(an$g)[i]$guild <- an$guild
  }

  an <<- new.env()
  # zinit_time <- proc.time()

  nread <- read_network(namenetwork, directory = directory, guild_astr = guild_a, guild_bstr = guild_b)
  an$g <- as.undirected(nread$g)
  m <- nread$matrix
  names_guild_a <- nread$names_guild_a
  names_guild_b <- nread$names_guild_b
  num_guild_b <- nread$num_guild_b
  num_guild_a <- nread$num_guild_a
  edge_matrix <- igraph::get.edges(an$g, E(an$g))
  spaths_mat <- shortest.paths(an$g)
  g_cores <- graph.coreness(an$g)

  wtc <- walktrap.community(an$g)
  #modularity(wtc)
  modularity_measure <- modularity(an$g, membership(wtc))


  if (plot_graphs){
    plot(an$g, vertex.size=8, layout=layout.kamada.kawai)
    hist(g_cores,right=FALSE)
  }
  lcores <- unique(g_cores)
  max_core <- max(lcores)
  min_core <- min(lcores)
  p <- rep(NA, length(lcores))
  for (k in lcores)
  {
    p[k] <- list(names(g_cores[g_cores == k]))
  }
  # Find plants and plos of a core

  plants_k <- list(rep(NA,max_core))
  pols_k <- list(rep(NA,max_core))
  for (i in 1:max_core)
  {
    auxincore <- c(p[[i]][grep(guild_a,as.character(unlist(p[i])))])
    if (length(auxincore)>0){
      plants_k[[i]] <- auxincore
    }
    else{
      plants_k[[i]] <- c(NA)
    }
    auxincore <- c(p[[i]][grep(guild_b,as.character(unlist(p[i])))])
    if (length(auxincore)>0){
      pols_k[[i]] <- auxincore
    }
    else{
      pols_k[[i]] <- c(NA)
    }
  }
  plants_maxcore <- p[[max_core]][grep(guild_a,as.character(unlist(p[max_core])))]
  pols_maxcore <- p[[max_core]][grep(guild_b,as.character(unlist(p[max_core])))]
  # Purge possible nodes of kcore number maximum that are not part of the giant component
  # Only one case detected, when kcoremax == 2, network PL_30

  if (max_core ==2)
  {
    meandiscnodes <- mean(spaths_mat== Inf)
    for (i in plants_maxcore)
      if (mean(spaths_mat[i,] == Inf)>2*meandiscnodes)
        plants_maxcore <- plants_maxcore[plants_maxcore!=i]
      for (i in pols_maxcore)
        if (mean(spaths_mat[i,] == Inf)>2*meandiscnodes)
          pols_maxcore <- pols_maxcore[pols_maxcore!=i]
  }
  V(an$g)$kradius <- NA
  V(an$g)$kcorenum <- NA
  V(an$g)$kdegree <- 0


  V(an$g)$guild <- ""
  E(an$g)$weights <- 1

  for(i in 1:length(E(an$g))){
    E(an$g)$weights[i] <- m[edge_matrix[i,][2]- num_guild_a,edge_matrix[i,][1]]
  }
  for (i in 1:max_core)
  {
    lnod <- p[[i]]
    if (sum(!is.na(lnod))>0){
      for (k in lnod)
        V(an$g)[k]$kcorenum <- i
    }
  }

  V(an$g)$krisk <- 0
  listanodos <- grep(guild_a,V(an$g)$name)
  an$guild <- guild_a
  an$guild_maxcore <- pols_maxcore
  lapply(listanodos, calc_kradius)

  listanodos <- grep(guild_b,V(an$g)$name)
  an$guild <- guild_b
  an$guild_maxcore <- plants_maxcore
  lapply(listanodos, calc_kradius)


  meandist <- mean(V(an$g)$kradius[V(an$g)$kradius != Inf])
  if (only_NODF)
    nested_values<- nested(as.matrix(m), method = "NODF")
  else
    nested_values<- nested(as.matrix(m), "ALL")

  # kdegree computation
  aux_graf <- data.frame(kdegree=V(an$g)$kdegree, kradius=V(an$g)$kradius ,
                         krisk=V(an$g)$krisk, kcorenum=V(an$g)$kcorenum)

  for (l in 1:nrow(edge_matrix))
  {
    polvertex = edge_matrix[l,2]
    plantvertex = edge_matrix[l,1]
    diff_kcorenum = aux_graf$kcorenum[plantvertex] - aux_graf$kcorenum[polvertex]
    aux_graf$kdegree[polvertex] = aux_graf$kdegree[polvertex] + 1/aux_graf$kradius[plantvertex]
    aux_graf$kdegree[plantvertex] = aux_graf$kdegree[plantvertex] + 1/aux_graf$kradius[polvertex]

    if ((aux_graf$kradius[polvertex] != Inf) & (diff_kcorenum < 0))
      if (aux_graf$kcorenum[polvertex] > 1)
        aux_graf$krisk[polvertex] = aux_graf$krisk[polvertex] - diff_kcorenum

    if ((aux_graf$kradius[plantvertex] != Inf) & (diff_kcorenum > 0 ))
      if (aux_graf$kcorenum[plantvertex] > 1)
        aux_graf$krisk[plantvertex] = aux_graf$krisk[plantvertex] + diff_kcorenum
  }

  V(an$g)$kdegree <- aux_graf$kdegree
  V(an$g)$kradius <- aux_graf$kradius
  V(an$g)$krisk <- aux_graf$krisk + 0.01*aux_graf$kcorenum
  V(an$g)$kcorenum <- aux_graf$kcorenum
  meankdegree <- mean(V(an$g)$kdegree)


  calc_values <- list("graph" = an$g, "max_core" = max_core, "nested_values" = nested_values, "num_guild_a" = num_guild_a,
                      "num_guild_b" = num_guild_b, "links" = length(E(an$g)), "meandist" = meandist, "meankdegree" = meankdegree,
                      "spaths_mat" = spaths_mat, "matrix" = as.matrix(m), "g_cores" = g_cores, "modularity_measure" = modularity_measure)
  return(calc_values)
}


#' Get bipartite labels
#'
#' Add guild labels to a bipartite network
#'
#' @param g is the newtork in \code{igraph::graph} format
#' @param strg_guild_a is a the label of the class of guild a nodes
#' @param strg_guild_b is a the label of the class of guild b nodes
#' @param plot_graphs set to FALSE, deprecated, kept for backwards compatibility
#' @export
#' @examples get_bipartite(result_analysis$graph, str_guild_a = "Plant", str_guild_b = "Fungus")

get_bipartite <- function(g, str_guild_a = "Plant", str_guild_b = "Pollinator", plot_graphs = FALSE)
{
  bg <- g
  V(bg)$type <- FALSE
  V(bg)$label <- ""
  for (i in V(bg)$name)
    if (length(grep(str_guild_b,i))>0){
      V(bg)[which(V(bg)$name==i)]$type <- TRUE
      V(bg)[which(V(bg)$name==i)]$label <- strsplit(i,str_guild_b)[[1]][2]
    }
  else
    V(bg)[which(V(bg)$name==i)]$label <- strsplit(i,str_guild_a)[[1]][2]
  V(bg)$name <- V(bg)$label
  return(bg)
}

#result_analysis <- analyze_network("M_SD_008.csv", directory = "data/", guild_a = "Plant", guild_b = "Pollinator", plot_graphs = FALSE)

