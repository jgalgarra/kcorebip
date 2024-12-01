library(reshape2)

#' Plotting a matrix graph
#'
#' This function plots the matrix graph of a bipartite network. Configuration parameters and
#' results are stored in a global environment called mat. This environment is not destroyed
#' after the function is executed, so the developer can store it in a configuration file, retrieve
#' 
#' @param datadir the name of the file of the interaction matrix
#' @param filename the file with the interaction matrix
#' @param orderkcoremaxby sets order of nodes, by kradius, kdegree or degree
#' @param label_strguilda string labels of guild a
#' @param label_strguildb string labels of guild b
#' @param label_size base label size
#' @param color_guild_a default filling for nodes of guild_a
#' @param color_guild_b default filling for nodes of guild_b
#' @param color_links color of matrix interactions
#' @param flip_matrix to switch row/column guilds, 
#' @param links_weight show links weight as a gradient of color for weighted networks
#' @param show_species_names add species name to number,
#' @param show_title show network name as plot title, 
#' @param show_legend show guilds legend,
#' @param print_to_file print plot to file, 
#' @param plotsdir default print directory
#' @param plot_size in inches for printed plot,
#' @param fscale_height scale plot height,
#' @param fscale_width scale plot width,
#' @param ppi dots per inch
#' @param progress for interactive visualization
#' @export
matrix_graph <-function(datadir,filename,
                         orderby = "kradius",
                         label_strguilda = "Plant",
                         label_strguildb = "Pollinator",
                         label_size = 18,
                         color_guild_a = "#4169E1", 
                         color_guild_b = "#FF1808",
                         color_links = 'grey7',
                         flip_matrix = FALSE, 
                         links_weight = FALSE,
                         show_species_names = TRUE,
                         show_title = TRUE, show_legend = TRUE,
                         print_to_file = FALSE, plotsdir ="plot_results/", 
                         plot_size = 10,
                         fscale_height = 1,
                         fscale_width = 1,
                         ppi = 300,
                         progress = NULL
                         )
{
  plot_m <- function(longData,flip_matrix=FALSE,nname="",ncolor="green",links_weight=FALSE,
                     strA="",strB="",colorA="blue",colorB="red",lsize=8,show_title=TRUE,show_legend=TRUE)
  {
    legends_text <- paste0(
      #"<span style = 'text-align: left; color:black'>Network: ",nname,
      "</span><span style = 'text-align: right; color:white'>....</span><span style = 'text-align: right; color:",colorA,"'>",paste('&#9632;',strA),
      "</span> <span style = 'text-align: right;color:",colorB,"'>",paste('&nbsp; &#9632;',strB),"</span>")
    
    if (!flip_matrix){
      if (!links_weight)
        mplots<- ggplot(longData, aes(x = speciesA, y = speciesB,fill=as.factor(value)))
      else
        mplots<- ggplot(longData, aes(x = speciesA, y = speciesB,fill=(value+1)))
    }
    else{
      if (!links_weight)
        mplots<- ggplot(longData, aes(x = rev(speciesB), y = rev(speciesA),fill=as.factor(value)))
      else
        mplots<- ggplot(longData, aes(x = rev(speciesB), y = rev(speciesA),fill=(value+1)))
    }  
    if (!links_weight)
      mplots <- mplots + geom_tile(color = ncolor,alpha=1) + scale_fill_manual(values=c("white",ncolor))
    else 
      mplots <- mplots+geom_tile(color = ncolor,)+scale_fill_continuous(type = "gradient", trans = "log",  
                                                                        low = "white", high = ncolor, 
                                                                        breaks = round)
    mplots <- mplots + ylab("")+xlab("")+coord_fixed()
    if (!flip_matrix){
      mplots <- mplots +scale_x_discrete(position="top",labels=unique(longData$numA))+
        scale_y_discrete(labels=unique(longData$numB))
    } else {
      mplots <- mplots +scale_x_discrete(position="top",labels=rev(unique(longData$numB)))+
        scale_y_discrete(labels=rev(unique(longData$numA)))
    }
    if (links_weight)
      lposition = "right"
    else
      lposition = "none"
    if (show_legend)
      mplots <- mplots+labs(caption=legends_text)
    if (show_title)
      mplots <- mplots+ggtitle(paste("Network:",mat$network_name))
    mplots <- mplots +  theme_void()+theme(legend.position=lposition,
                                           legend.title=element_blank(),
                                           legend.text = element_text(size = lsize-1),
                                           axis.text.x = element_text(size=lsize,hjust=0,vjust=1,angle=90,color=ifelse(!flip_matrix,colorA,colorB)),
                                           axis.text.y = element_text(size = lsize,hjust=1,vjust=0.5,color=ifelse(!flip_matrix,colorB,colorA)),
                                           plot.title = element_text(size=lsize+3,hjust=0.5),
                                           axis.title.x=element_text(size=18,face="bold",color="grey40",
                                                                     hjust=0.5),
                                           plot.subtitle = ggtext::element_markdown(hjust=0,vjust=0),
                                           plot.margin = unit(c(0,0,0,0), "cm"),
                                           plot.caption = ggtext::element_markdown(size=lsize+2,hjust=1,vjust=0))
    return(mplots)
  }
  
  create_labels <- function(M,nums,i,vnames,show_species=FALSE,flip_matrix=TRUE,guild="A"){
    if (guild=="A")
      mn <- vnames[which(vnames$name==colnames(M)[i]),]$num
    else
      mn <- vnames[which(vnames$name==rownames(M)[i]),]$num
    if (!flip_matrix){
      if (guild=="A")
        label <- (paste0("   ",mn," ",ifelse(show_species,colnames(M)[i],""),"   "))
      else
        label <- (paste0(" ",mn," ",ifelse(show_species,paste0(rownames(M)[i],""),""),"    "))
    } else {
      if (guild=="A")
        label <- (paste0(" ",mn," ",ifelse(show_species,colnames(M)[i],""),"    "))
      else
        label <- (paste0(" ",mn," ",ifelse(show_species,paste0(rownames(M)[i],"  "),"    ")))
    }
    return(label)
  }
  
  

  setkradiusorder <- function(dnodes){
    myord <- rev(order(10000*dnodes$kshell-(100*dnodes$kradius-dnodes$degree)))
    return(dnodes$num[myord])
  }
  
  setkdegreeorder <- function(dnodes){
    myord <- rev(order(100*dnodes$kdegree-dnodes$kradius))
    return(dnodes$num[myord])
  }
  
  setdegreeorder <- function(dnodes){
    myord <- rev(order(100*dnodes$degree-dnodes$kradius))
    return(dnodes$num[myord])
  }
  
  mat_argg <- c(as.list(environment()))
  # Create global environment

  f<-read_and_analyze(datadir,filename,label_strguilda, label_strguildb)
  mat <<- new.env()
  mat$result_analysis <- f["result_analysis"][[1]]
  mat$str_guild_a <- f["str_guild_a"][[1]]
  mat$str_guild_b <- f["str_guild_b"][[1]]
  mat$name_guild_a <- f["name_guild_a"][[1]]
  mat$name_guild_b <- f["name_guild_b"][[1]]
  mat$network_name <- f["network_name"][[1]]
  mat$mat_argg <- mat_argg
  lnodes <- V(mat$result_analysis$graph)
  species_a <- data.frame("num"=seq(1:mat$result_analysis$num_guild_a),"name"=colnames(mat$result_analysis$matrix))
  species_a$kdegree <- lnodes$kdegree[1:mat$result_analysis$num_guild_a]
  species_a$kradius <- lnodes$kradius[1:mat$result_analysis$num_guild_a]
  if (sum(species_a$kradius==Inf)>0)
    species_a[species_a$kradius==Inf,]$kradius=100
  species_a$kshell <- lnodes$kcorenum[1:mat$result_analysis$num_guild_a]
  species_a$num <- seq(1:mat$result_analysis$num_guild_a)
  #species_a$num <- species_a[order(species_a$name),]$num
  species_b <- data.frame("num"=seq(1:mat$result_analysis$num_guild_b),"name"=rownames(mat$result_analysis$matrix))
  species_b$kdegree <- lnodes$kdegree[seq(mat$result_analysis$num_guild_a+1,mat$result_analysis$num_guild_a+mat$result_analysis$num_guild_b)]
  species_b$kradius <- lnodes$kradius[seq(mat$result_analysis$num_guild_a+1,mat$result_analysis$num_guild_a+mat$result_analysis$num_guild_b)]
  if (sum(species_b$kradius==Inf)>0)
    species_b[species_b$kradius==Inf,]$kradius=100
  species_b$kshell <- lnodes$kcorenum[seq(mat$result_analysis$num_guild_a+1,mat$result_analysis$num_guild_a+mat$result_analysis$num_guild_b)]
  species_b$num <- seq(1:mat$result_analysis$num_guild_b)
  M<-mat$result_analysis$matrix
  binary_network = (sum(M>1)==0)
  mat$binary_network <- binary_network
  mat$network_type <- ifelse(mat$binary_network,strings$value("LABEL_ZIGGURAT_INFO_BINARY"),strings$value("LABEL_ZIGGURAT_INFO_WEIGHTED"))
  DegM <- M
  DegM[DegM>1] <- 1
  species_a$degree <- colSums(DegM)
  species_b$degree <- rowSums(DegM)
  # Vector that stores the node display order
  if (orderby == "kradius"){
    ordvectorA <- setkradiusorder(species_a)
    ordvectorB <- rev(setkradiusorder(species_b))
  }
  if (orderby == "kdegree"){
    ordvectorA <- setkdegreeorder(species_a)
    ordvectorB <- rev(setkdegreeorder(species_b))
  }
  if (orderby == "degree"){
    ordvectorA <- setdegreeorder(species_a)
    ordvectorB <- rev(setdegreeorder(species_b))
  }
  species_a <- species_a[ordvectorA,]
  species_b <- species_b[ordvectorB,]
  if (!links_weight)
    M[M>1]=1
  M<-M[ordvectorB,ordvectorA]
  num_b = species_b$num
  num_a = species_a$num
  longData<-melt(M)
  mplots<- ggplot(longData, aes(x = Var2, y = Var1,fill=as.factor(value)))+geom_tile(alpha=1) + scale_fill_manual(values=c("white","red"))
  longData$numA <- 0
  longData$numB <- 0
  names(longData) = c("speciesB","speciesA", "value","numB","numA")
  mplots2<- ggplot(longData, aes(x = speciesA, y = speciesB,fill=as.factor(value)))+geom_tile(alpha=1) + scale_fill_manual(values=c("white","blue"))
  for (i in 1:length(colnames(M)))
    longData[longData$speciesA==colnames(M)[i],]$numA = create_labels(M,num_a,i,species_a,show_species=show_species_names,flip_matrix=flip_matrix,guild="A")
  for (i in 1:length(rownames(M)))
    longData[longData$speciesB==rownames(M)[i],]$numB = create_labels(M,num_b,i,species_b,show_species=show_species_names,flip_matrix=flip_matrix,guild="B")
  lsize <- label_size - round(log10(sqrt(nrow(longData))))
  if ( max(nrow(species_a),nrow(species_b)>25) && (min(nrow(species_a),nrow(species_b))>15)){
    lsize <- lsize-5
  }
  
  p <- plot_m(longData,flip_matrix=flip_matrix,nname=mat$network_name,
              strA=mat$name_guild_a,strB=mat$name_guild_b,links_weight = (links_weight && !binary_network),
              colorA=color_guild_a,colorB=color_guild_b,lsize=lsize,ncolor=color_links,show_title = show_title,
              show_legend=show_legend)
  fscaleheight = fscale_height
  fscalewidth = fscale_width
  plsize = plot_size
  dppi = ppi
  aspect = (nrow(species_a)+6)/(nrow(species_b)+6) # 8 boxes are more or less the species name length
  imw = plsize*fscalewidth*dppi
  imh = plsize*fscaleheight*dppi/aspect
  # User decides to plot the file
  if (print_to_file){
    dir.create(mat$mat_argg$plotsdir, showWarnings = FALSE)
    nfile <- paste0(plotsdir,mat$network_name,"_MATRIX_orderby_",orderby,".png")
    if (!flip_matrix)
      png(nfile,width=imw,height=imh,res=dppi)
    else
      png(nfile,width=imh,height=imw,res=dppi)
    print(p)
    dev.off()    
  }
  if (!is.null(progress)){
    ppi = 300
    dir.create("tmp", showWarnings = FALSE)
    nfile <- paste0("tmp/",mat$network_name,"_MATRIX.png")
    if (!flip_matrix)
      png(nfile,width=imw,height=imh,res=dppi)
    else
      png(nfile,width=imh,height=imw,res=dppi)
    print(p)
    dev.off()
    mat$matrix_file <- nfile
  }
  mat$aspect <- aspect
  mat$plot <- p
  return(mat)
}

debugging = FALSE
if (debugging)
  p <- matrix_graph("../data/","dattilo2014.csv",
                  print_to_file = TRUE, plotsdir ="plot_results/", 
                  orderby = "kradius",ppi=300,
                  flip_matrix = FALSE, links_weight = FALSE,
                  show_species_names = TRUE,
                  show_title = TRUE,
                  show_legend = TRUE)

