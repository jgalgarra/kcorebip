library(igraph)
library(bipartite)
library(ggplot2)

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
#' @examples result_analysis <- analyze_network("M_SD_001.csv", directory = "../data/", plot_graphs = FALSE,sep=",",speciesinheader=TRUE)


analyze_network <- function(namenetwork, directory="", guild_a = "pl", guild_b = "pol", plot_graphs = FALSE, only_NODF = FALSE,
                            weight_direction = "none",sep=",",speciesinheader=TRUE)
{
  if(!exists("an")){
    an <<- new.env()  
  }
  an$sep <- sep
  an$speciesinheader <- speciesinheader
  # K radius is the average distance to nodes of maximum k index of the opposite guild
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
  
  # Read interaction matrix file
  nread <- read_network(namenetwork, directory = directory, guild_astr = guild_a, guild_bstr = guild_b, sep=an$sep,
                        speciesinheader = an$speciesinheader)
  # create empty graph
  #an$g <- as.undirected(nread$g)
  an$g <- as_undirected(nread$g)
  m <- nread$matrix
  names_guild_a <- nread$names_guild_a
  names_guild_b <- nread$names_guild_b
  num_guild_b <- nread$num_guild_b
  num_guild_a <- nread$num_guild_a
  # Get egde matrix
  edge_matrix <- igraph::get.edges(an$g, E(an$g))
  # Compute all shortest paths in network
  #spaths_mat <- shortest.paths(an$g)
  spaths_mat <- distances(an$g)
  # k-core decompose the network
  #g_cores <- graph.coreness(an$g)
  g_cores <- coreness(an$g)
  
  #wtc <- walktrap.community(an$g)
  wtc <- cluster_walktrap(an$g)
  #modularity(wtc)
  modularity_measure <- modularity(an$g, membership(wtc))
  
  # This option to plot graphs is only useful for a preliminary exam
  # becasuse quality is not good enough
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
  # Find plants and pollinatorss of a core
  # In general, you must read 'plants' as 'guildA' and 'pols' as 'guildB' species
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
  # Nodes that are not part of the Giant Component
  V(an$g)$kradius <- NA
  V(an$g)$kcorenum <- NA
  V(an$g)$kdegree <- 0
  
  
  V(an$g)$guild <- ""
  E(an$g)$weights <- 1
  
  # Weight of links
  for(i in 1:length(E(an$g))){
    E(an$g)$weights[i] <- m[edge_matrix[i,][2]- num_guild_a,edge_matrix[i,][1]]
  }
  # Assign k-inedx to each species
  for (i in 1:max_core)
  {
    lnod <- p[[i]]
    if (sum(!is.na(lnod))>0){
      for (k in lnod)
        V(an$g)[k]$kcorenum <- i
    }
  }
  
  # k-risk initialization
  V(an$g)$krisk <- 0
  
  # k-radius computation
  listanodos <- grep(guild_a,V(an$g)$name)
  an$guild <- guild_a
  an$guild_maxcore <- pols_maxcore
  lapply(listanodos, calc_kradius)
  
  listanodos <- grep(guild_b,V(an$g)$name)
  an$guild <- guild_b
  an$guild_maxcore <- plants_maxcore
  lapply(listanodos, calc_kradius)
  
  meandist <- mean(V(an$g)$kradius[V(an$g)$kradius != Inf])
  # Computation of nestedness measures
  if (only_NODF)
    nested_values<- nested(as.matrix(m), method = "NODF")
  else
    nested_values<- nested(as.matrix(m), "ALL")
  
  # kdegree computation
  aux_graf <- data.frame(kdegree=V(an$g)$kdegree, kradius=V(an$g)$kradius ,
                         krisk=V(an$g)$krisk, kcorenum=V(an$g)$kcorenum)
  
  # k-risk computation
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
#' @param g is the network in \code{igraph::graph} format
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

read_network <- function(namenetwork, guild_astr = "pl", guild_bstr = "pol", directory="", sep=",", speciesinheader=TRUE)
{
  # Reading species names
  if (speciesinheader){
    namesred <- read.csv(paste0(directory,namenetwork),header=FALSE,stringsAsFactors=FALSE,sep = sep)
    names_guild_a <- unname(unlist(namesred[1,2:ncol(namesred)])) # JULY 2023 
    names_guild_b <- namesred[2:nrow(namesred),1]
    m <- read.csv(paste0(directory,namenetwork),header=TRUE,row.names=1,sep=an$sep)
  } else {
    m <- read.csv(paste0(directory,namenetwork),header=FALSE,sep=an$sep)
    for (i in 1:ncol(m))
      colnames(m)[i] <- paste0("A",i)
    for (i in 1:nrow(m))
      rownames(m)[i] <- paste0("B",i)
    names_guild_a <- colnames(m)
    names_guild_b <- rownames(m)
  }
  
  # Calc number of species of each guild
  num_guild_a <- ncol(m)
  num_guild_b <- nrow(m)
  # Create an graph object
  g <- make_empty_graph()
  #g <- empty_graph()
  # Add one node for each species and name it
  for (i in 1:num_guild_a){
    g <- g + vertices(paste0(guild_astr,i),color="white",guild_id="a",name_species=names_guild_a[i],id=i)
  }
  for (i in 1:num_guild_b){
    g <- g + vertices(paste0(guild_bstr,i),color="red",guild_id="b",name_species=names_guild_b[i],id=i)
  }
  
  # Adding links to the graph object
  mm <- matrix(unlist(list(m)),nrow=num_guild_b,ncol=num_guild_a)
  listedgesn <- which(mm!=0, arr.ind = T)
  listedgesn <- listedgesn[order(listedgesn[,1],listedgesn[,2]),]
  listedgesn[,1] <- paste0(guild_bstr,listedgesn[,1])
  listedgesn[,2] <- paste0(guild_astr,listedgesn[,2])
  g <- g + graph.edgelist(listedgesn)
  #g <- g + as_edgelist(listedgesn)
  # Return values
  calc_values <- list("graph" = g, "matrix" = m, "num_guild_b" = num_guild_b, "num_guild_a" = num_guild_a,
                      "names_guild_a" = names_guild_a, "names_guild_b"=names_guild_b)
  return(calc_values)
  
}

#' Reads and analyzes a network from a file
#'
#' @param directorystr where the network file is stored
#' @param network_file is the file name that contains the interaction matrix
#' @param label_strguilda name of the guild A of nodes
#' @param label_strguilda name of the guild B of nodes
#' @return \code{calc_vals} a list containing the following objects
#' \itemize{
#'  \item{\code{"result_analysis"}}{ a full network analysis}
#'  \item{\code{"str_guild_a"}}{ list of nodes of Guild A}
#'  \item{\code{"str_guild_b"}}{ list of nodes of Guild B}
#'  \item{\code{"name_guild_a"}}{ name of Guild A}
#'  \item{\code{"network_name"}}{ network name}
#'  \item{\code{"spe"}}{ separator character}
#'  \item{\code{"speciesinheader"}}{ species names in matrix header}
#'  }
#' @export
#' @examples p <- read_and_analyze("../data","M_PL_003.csv","Plant","Pollinator")
read_and_analyze <- function(directorystr,network_file,label_strguilda,label_strguildb,sep=",",speciesinheader=TRUE)
{
  if (!exists("an")){
    an <<- new.env() 
    an$sep <- sep
    an$speciesinheader <- speciesinheader
  }
  str_guild_a <- "pl"
  str_guild_b <- "pol"
  name_guild_a <- "GuildA"
  name_guild_b <- "GuildB"
  network_name <- strsplit(network_file,".csv")[[1]][1]
  slabels <- c("Plant", "Pollinator")
  if (grepl("_SD_",network_name)){
    str_guild_b <- "disp"
    name_guild_b <- "Dispersers"
  }
  
  if (nchar(label_strguilda)>0){
    slabels <- c(label_strguilda, label_strguildb)
    name_guild_a <- label_strguilda
    name_guild_b <- label_strguildb
  }
  
  result_analysis <- analyze_network(network_file, directory = directorystr, guild_a = str_guild_a,
                                     guild_b = str_guild_b, only_NODF = TRUE, sep=an$sep,
                                     speciesinheader = an$speciesinheader)
  
  
  calc_vals <- list("result_analysis" = result_analysis, "str_guild_a" = str_guild_a, "str_guild_b" = str_guild_b,
                    "name_guild_a" = name_guild_a, "name_guild_b" = name_guild_b,
                    "network_name" = network_name)
  return(calc_vals)
}



# EXAMPLE. Copy and paste to test.
#result_analysis <- analyze_network("M_SD_001.csv", directory = "../data/", plot_graphs = FALSE,sep=",",speciesinheader=TRUE)
#result_analysis <- analyze_network("kaka.csv", directory = "../data/", plot_graphs = FALSE,sep=";",speciesinheader=FALSE)