library(ggplot2)
library(scales)
library(grid)
library(gridExtra)
library(igraph)

paint_kdegree_kradius <- function(graph, num_guild_a, num_guild_b,
                                  lsize_title , lsize_axis, lsize_legend, lsize_axis_title ,
                                  lsize_legend_title,
                                  showtext = "no",
                                  network_name = "",
                                  NODF = 0, Modularity = 0, MeanKradius = 0, MeanKdegree = 0,
                                  printable_range = printable_labels,
                                  fname_append = "",
                                  progress
                                  )
{
  g <- V(graph)
  nga <- sum(g[1:num_guild_a]$kradius != Inf)
  ngb <- sum(g[as.numeric(num_guild_a+1):length(g)]$kradius != Inf)
  g <- g[g$kradius != Inf]
  dfaux <- data.frame(g$kradius,g$kdegree,g$kcorenum,(g$kdegree/max(g$kdegree))^1.5)
  names(dfaux) <- c("kradius","kdegree","kcorenum","normdegree")
  scale_factor <- 20
  dfaux$kradio <- scale_factor*sqrt(dfaux$kdegree/pi)
  dfaux$posx <- NA
  dfaux$posy <- NA
  dfaux$name <- NA
  dfaux$symbol <- 1
  dfaux$kcol_label <- NA
  dfaux$despl <- 1.2
  dfaux$name <- as.character(g$name)
  signo <- 1
  guarda <- 0.25
  min_radius <- 0
  tot_species <- nrow(dfaux)
  rndmult <- runif(tot_species,guarda,pi-guarda)
  rndmult_a <- sample(seq(guarda,pi-guarda,length.out=nga))
  rndmult_b <- sample(seq(guarda,pi-guarda,length.out=ngb))
  dfaux <- dfaux[dfaux$kradius != Inf,]
  maxcore <- max(dfaux$kcorenum)
  extreme <- ceiling(max(dfaux[dfaux$kradius != Inf,]$kradius))
  num_central <- (nga+ngb)%/%5
  more_central_nodes <- head(dfaux[order(dfaux$kradius),]$name, num_central)
  slice_multiplier <- 4
  rnd_central <- seq(guarda,pi-guarda,length.out = num_central*slice_multiplier)
  pal <-colorRampPalette(c("cadetblue","darkorchid4"))
  #jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  jet.colors <- colorRampPalette(c("blue","red"))
  #vcols <- pal(max(8,maxcore))
  vcols <- jet.colors(max(8,maxcore))
  alpha_level <- 0.3
  k <- 1
  if (printable_range >0){
    tailp <- 3
    sort_radiuss <- dfaux[order(dfaux$kradius),]$name
    printable_points <- c(head(sort_radiuss,printable_range), tail(sort_radiuss,tailp))
  } else
    printable_points = c()
  for (i in 1:tot_species){
    if (!is.null(progress)) progress$inc((3/4)*(1/tot_species), detail=paste0(strings$value("MESSAGE_POLAR_PROGRESS_PROCESING_SPECIE"), " ", i , "..."))
    if (length(which(printable_points == dfaux[i,]$name)) > 0)
      dfaux[i,]$kcol_label <- vcols[dfaux[i,]$kcorenum]
    if (i>nga)
    {
      offset <- pi
      dfaux[i,]$symbol <- 21
      rndmult <- rndmult_b[i-nga]
    }
    else{
      offset <- 0
      dfaux[i,]$symbol <- 20
      rndmult <- rndmult_a[i]
    }
    if (length(which(more_central_nodes == dfaux[i,]$name)) > 0){
      rdnmult <- rnd_central[ceiling(k*runif(1,slice_multiplier,slice_multiplier*1.5)*i%/%length(rnd_central))]
      k <- k + 1
    }
    if (dfaux[i,]$kradius != Inf){
      dfaux[i,]$posy <- dfaux[i,]$kradius
      dfaux[i,]$posx <- rndmult + offset
    }
    else{
      dfaux[i,]$posx <- 0
      dfaux[i,]$posy <- 0
      dfaux[i,]$kdegree <- 0.0001
    }
  }
  dfaux$fillcol <- 1 + maxcore - dfaux$kcorenum
  polar_plot <- ggplot(dfaux, aes(x=posx,y=posy),legendTextFont=c(15, "bold.italic", "red")) +
    scale_size_area(max_size=scale_factor,name="k-degree") +
    scale_colour_manual(values = vcols,name="k-shell") +
    guides(col = guide_legend(override.aes = list(shape = 15, size = 8)),
           shape = guide_legend(override.aes = list(size = 8, colour = "slategray1")),
           kdegree = guide_legend(override.aes = list(shape = 15, size = 8, colour = "slategray1")))
  if (showtext == "yes"){
    polar_plot <- polar_plot+ geom_text(aes(size=0.005+0.1*normdegree,angle=0,colour = factor(kcorenum), label = name), alpha = alpha_level+0.1)+
      scale_shape_identity()
  }
  else{
    polar_plot <- polar_plot + geom_point(aes(size=kdegree, colour = factor(kcorenum), shape = factor(symbol)), alpha = alpha_level) +
      scale_shape_manual(values=c(16,15),name="Guild",labels=slabels ) +
      annotate(geom="text", x=dfaux$posx, y=dfaux$posy, label=dfaux$name, colour = factor(dfaux$kcol_label), size=2*(3+4*dfaux$normdegree), hjust = 1, alpha = 1)
  }
  polar_plot <- polar_plot + coord_polar(start = -pi/2) + labs(x = '', y = '')# + theme(axis.text.x = element_blank())
  polar_plot <- polar_plot + scale_y_continuous(breaks=seq(min_radius,extreme), lim=c(min_radius, extreme),labels=seq(min_radius,extreme) )
  polar_plot <- polar_plot + scale_x_continuous(breaks=seq(0, 2*pi, by=pi/2), lim=c(0,2*pi))
  polar_plot <- polar_plot+ theme_bw() + theme(panel.border = element_blank(),
                                               legend.key = element_blank(),
                                               axis.ticks.y = element_blank(),
                                               axis.ticks.x = element_blank(),
                                               panel.grid.major.x = element_blank(),
                                               panel.grid.minor.x = element_blank(),
                                               axis.text.y = element_blank(),
                                               panel.grid.major.y = element_line(size = 0.5, linetype = 2, alpha("darkolivegreen",0.5)),
                                               panel.grid.minor.y = element_blank(),
                                               panel.border = element_blank(),
                                               axis.text.x = element_blank(),
                                               legend.text = element_text(size=lsize_legend),
                                               legend.title = element_text(size=lsize_legend_title),
                                               plot.title = element_text(size=lsize_title,lineheight=.8, face="bold")
                                        )
  ylab <- seq(0,extreme)
  pylab <- ylab
  pylab[2:length(pylab)] <- pylab[2:length(pylab)]-0.05
  ylab[1] <- "k-radius"
  xlab <- rep(pi,length(pylab))
  dftext <- data.frame(xlab,ylab,pylab)
  dftext$fillcol <- maxcore
  polar_plot <- polar_plot + annotate(geom="text",x=xlab,y=pylab,label=ylab,size=5, color="gray20", lineheight=.8)
  polar_plot <- polar_plot + ggtitle(sprintf("Network %s NODF: %.02f Modularity: %.04f\n\n Avg k-radius: %.02f Avg k-degree: %.02f", network_name, NODF, Modularity, MeanKradius, MeanKdegree)) +
    guides(row = guide_legend(nrow = 1))
  histo_dist <- ggplot(dfaux, aes(kradius)) + geom_histogram(alpha = alpha_level,binwidth=extreme/15,
                                                             color="white",fill = "forestgreen") +
    xlim(0,extreme) +
    theme_bw() +
    theme(panel.border = element_blank(),
          legend.key = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_line(linetype = 2, color="gray90"),
          panel.grid.major.x = element_line(linetype = 2, color="gray90"),
          panel.border = element_blank(),
          axis.line = element_line(colour = "black",
                                   size = 0.5, linetype = "solid"),
          legend.text = element_text(size=lsize_legend),
          plot.title = element_text(size=lsize_title,lineheight=.8, face="bold"),
          axis.text = element_text(face="bold", size=lsize_axis),
          axis.title.y  = element_text(face="bold", size=lsize_axis_title),
          axis.title.x = element_blank()
    )+
    ggtitle("k-radius") + ylab("Species")
#   if (log_histograms)
#     histo_dist <- histo_dist + scale_x_log10()
  histo_core <- ggplot(dfaux, aes(x=kcorenum)) +
                geom_histogram(alpha = alpha_level, binwidth = 1,color="white",fill = "slategray3") + theme(legend.position = "none") +theme_bw() +
    #scale_x_continuous(breaks=seq(0, maxcore, by=1), lim=c(1,maxcore+1)) +
    theme(panel.border = element_blank(),
          legend.key = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_line(linetype = 2, color="gray90"),
          panel.grid.major.x = element_line(linetype = 2, color="gray90"),
          axis.line = element_line(colour = "black",
                                   size = 0.5, linetype = "solid"),
          panel.border = element_blank(),
          legend.text = element_text(size=lsize_legend),
          plot.title = element_text(size=lsize_title,lineheight=.8, face="bold"),
          axis.text.y = element_text(face="bold", size=lsize_axis),
          axis.text.x = element_text(angle = 0, face="bold", size=lsize_axis),
          axis.title.y  = element_text(face="bold", size=lsize_axis_title),
          axis.title.x = element_blank()
    ) +
    ggtitle("k-shell")+ ylab("Species")
  histo_degree <- ggplot(dfaux, aes(kdegree)) + geom_histogram(alpha = alpha_level,binwidth=max(dfaux$kdegree)/15,
                                                               color="white",fill = "grey20") +
    xlim(0,ceiling(max(dfaux$kdegree))) +
    theme_bw() +
    theme(panel.border = element_blank(),
          legend.key = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_line(linetype = 2, color="gray90"),
          panel.grid.major.x = element_line(linetype = 2, color="gray90"),
          axis.line = element_line(colour = "black",
                                   size = 0.5, linetype = "solid"),
          panel.border = element_blank(),
          legend.text = element_text(size=lsize_legend),
          plot.title = element_text(size=lsize_title,lineheight=.8, face="bold"),
          axis.text = element_text(face="bold", size=lsize_axis),
          axis.title.x = element_blank(),
          axis.title.y  = element_text(face="bold", size=lsize_axis_title)
    )+
    ggtitle("k-degree")+ ylab("Species")
  calc_grafs <- list("polar_plot" = polar_plot, "histo_dist" = histo_dist, "histo_core" = histo_core,
                     "histo_degree" = histo_degree)
  return(calc_grafs)
}

#' Plotting a polar graph
#'
#' This function plots the polar graph of a bipartite network and the histograms of kshell
#' kradius and kdegree
#'
#' @param nname the name of the file of the interaction matrix
#' @param directorystr the directory where the \code{nname} file is stored
#' @param plotsdir the directory where the plot is stored
#' @param print_to_file if set to FALSE the plot is displayed in the R session window
#' @param pshowtext auxiliar for interactive apps, do not modify
#' @param show_histograms display the histograms if set TRUE
#' @param glabels guild labels
#' @param gshortened guild shortened labels
#' @param lsize_title title label size
#' @param lsize_axis axis label size
#' @param lsize_legend legend label size
#' @param lsize_axis_title axis title size
#' @param lsize_legend_title legend label size
#' @param file_name_append a label that the user may append to the plot file name for convenience
#' @param auxiliar for interactive apps, do not modify
#' @param printable_labels range of labeled species
#' @export
#' @examples polar_graph("M_PL_007.csv","data/",plotsdir="grafresults/",print_to_file = TRUE)

polar_graph <- function( nname, directorystr = "data/", plotsdir = "plot_results/polar/", print_to_file = FALSE, pshowtext = FALSE,
                         show_histograms = TRUE, glabels = c("Plant", "Pollinator"),
                         gshortened = c("pl","pol"),
                         lsize_title = 22, lsize_axis = 12, lsize_legend = 13,
                         lsize_axis_title = 14, lsize_legend_title = 15,
                         file_name_append = "",
                         progress=NULL, printable_labels = 0)
{
  nname_name <- strsplit(nname,".csv")[[1]][1]
  sguild_a <<- gshortened[1]
  sguild_b <<- gshortened[2]
  slabels <<- glabels
  if (grepl("_SD_",nname)){
    sguild_b = "disp"
    slabels <<- c("Plant", "Disperser")
  }

  if (!is.null(progress)) progress$inc(1/4, detail=strings$value("MESSAGE_POLAR_PROGRESS_ANALYZING_NETWORK"))
  result_analysis <- analyze_network(nname, directory = directorystr, guild_a = sguild_a, guild_b = sguild_b, plot_graphs = FALSE)
  numlinks <- result_analysis$links

  if (print_to_file){
    dir.create(plotsdir, showWarnings = FALSE)
    ppi <- 600
    if (file_name_append != "")
      ftname_append <- paste0("_",file_name_append)
    else
      ftname_append <- file_name_append
    if (show_histograms)
      png(paste0(plotsdir,nname_name,"_polar",ftname_append,".png"), width=12*ppi, height=12*ppi, res=ppi)
    else
      png(paste0(plotsdir,nname_name,"_polar",ftname_append,".png"), width=9*ppi, height=9*ppi, res=ppi)
  }

  r <- paint_kdegree_kradius(result_analysis$graph, result_analysis$num_guild_a,result_analysis$num_guild_b,
                             lsize_title , lsize_axis, lsize_legend, lsize_axis_title , lsize_legend_title,
                             network_name = nname_name, NODF = result_analysis$nested_values["NODF"],
                             Modularity =  result_analysis$modularity_measure,
                             MeanKradius = result_analysis$meandist, MeanKdegree = result_analysis$meankdegree,
                             showtext = pshowtext, fname_append = ftname_append,
                             printable_range = printable_labels, progress
                              )
  # Non interactive mode. Plots stored in file or displayer in R window
  if (is.null(progress))
  {
    if (show_histograms)
      grid.arrange(r["polar_plot"][[1]], nrow=2, heights=c(3/4,1/4),arrangeGrob(r["histo_dist"][[1]], r["histo_degree"][[1]], r["histo_core"][[1]],ncol=3, nrow=1, widths=c(1/3,1/3,1/3)))
    else
      print(r["polar_plot"][[1]])
    if (print_to_file)
      dev.off()
  }
  # Message for interactive apps.
  if (!is.null(progress)) {
    progress$inc(0, detail=strings$value("MESSAGE_POLAR_PROGRESS_DONE"))
    return(r)
  }
}

#polar_graph("pl017-minus6plants.csv","datanetworks2015/",print_to_file=TRUE, lsize_title = 24, lsize_axis = 18, lsize_legend = 18, lsize_axis_title = 18, lsize_legend_title = 20)
#source("jga/network-kanalysis.R"), lsize_title = 24, lsize_axis = 18, lsize_legend = 18, lsize_axis_title = 18, lsize_legend_title = 20)
#polar_graph("M_PL_021.csv","data/",print_to_file=TRUE, lsize_title = 24, lsize_axis = 18, lsize_legend = 18, lsize_axis_title = 18, lsize_legend_title = 20)
#polar_graph("M_SD_007.csv","data/",print_to_file=TRUE, lsize_title = 24, lsize_axis = 18, lsize_legend = 18, lsize_axis_title = 18, lsize_legend_title = 20)
#polar_graph("M_PL_031.csv","data/",print_to_file=TRUE, lsize_title = 24, lsize_axis = 18, lsize_legend = 18, lsize_axis_title = 18, lsize_legend_title = 20)
