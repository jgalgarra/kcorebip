library(ggplot2)
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
                                  ptitle = "",
                                  filln = FALSE,
                                  alphal = 0.5,
                                  nfsal = "",
                                  maxkradius = 0,
                                  progress
                                  )
{
  g <- V(graph)
  nga <- sum(g[1:num_guild_a]$kradius != Inf)
  ngb <- sum(g[as.numeric(num_guild_a+1):length(g)]$kradius != Inf)
  numspecies <- nga + ngb
  g <- g[g$kradius != Inf]
  dfaux <- data.frame(g$kradius,g$kdegree,g$kcorenum,(g$kdegree/max(g$kdegree))^1.5)
  names(dfaux) <- c("kradius","kdegree","kcorenum","normdegree")
  scale_factor <- 20
  dfaux$symbolradius <- scale_factor*sqrt(dfaux$kdegree)
  dfaux$posx <- 0
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
  dfaux <- dfaux[dfaux$symbolradius != Inf,]
  maxcore <- max(dfaux$kcorenum)
  numkcores <- length(unique(dfaux$kcorenum))
  if (maxkradius == 0)
    extreme <- ceiling(max(dfaux[dfaux$symbolradius != Inf,]$kradius))
  else
    extreme <- maxkradius
  num_central <- (nga+ngb)%/%5
  more_central_nodes <- head(dfaux[order(dfaux$kradius),]$name, num_central)
  slice_multiplier <- 4
  rnd_central <- seq(guarda,pi-guarda,length.out = num_central*slice_multiplier)
  pal <-colorRampPalette(c("cadetblue","darkorchid4"))
  jet.colors <- colorRampPalette(c("slategray","gray10"))
  vcols <- jet.colors(maxcore)
  alpha_level <- 1
  k <- 1
  symbol_a <- 22
  symbol_b <- 21
  dfaux$classe <- "A"
  if (printable_range >0){
    sort_radiuss <- dfaux[order(dfaux$kradius),]$name
    printable_points <- head(sort_radiuss,printable_range)
  } else
    printable_points = c()
  dfaux$posx <- 0
  for (i in 1:tot_species){
    if (!is.null(progress)) progress$inc((3/4)*(1/tot_species), detail=paste0(strings$value("MESSAGE_POLAR_PROGRESS_PROCESING_SPECIE"), " ", i , "..."))
    if (length(which(printable_points == dfaux[i,]$name)) > 0)
      dfaux[i,]$kcol_label <- vcols[dfaux[i,]$kcorenum]
    if (i>nga)
    {
      offset <- pi
      dfaux$symbol[i] <- symbol_b
      dfaux$classe[i] <- "B"
    }
    else{
      offset <- 0
      dfaux$symbol[i] <- symbol_a
    }
    if (dfaux$kradius[i] != Inf)
      dfaux$posy[i] <- dfaux[i,]$kradius
    else{
      dfaux$posy[i] <- 0
      dfaux$kdegree[i] <- 0.0001
    }
  }

  # dfaux ordered by class and kradius
  dfaux <- dfaux[rev(order(dfaux$classe,dfaux$kradius)),]
  dfaux$orderincore <- 1
  dfaux$vjust <- 0.5
  dfaux$num <-na.omit(as.numeric(unlist(strsplit(unlist(dfaux$name), "[^0-9]+"))))
  primemove <- 7
  topradius <- max(dfaux$kradius)
  counttop <- sum(dfaux$kradius == topradius)
  indtop <- 0
  indvulg <-c("A"=0,"B"=0)
  dfaux$posx[1] <- primemove%%(pi-guarda)
  saltovert <- 1.5
  shift <- 0
  denom <- pi-guarda
  for (j in 2:nrow(dfaux))
  {
    if (dfaux$kradius[j] == topradius){
      dfaux$posx[j] <- indtop*(pi-guarda)/counttop
      indtop <- indtop + 1
    }
    else{
      shift  <- 0
      dfaux$posx[j] <- ((1+shift)*(indvulg[dfaux$classe[j]]*primemove))%%denom
      if (sum( (dfaux$posx == dfaux$posx[j])& (dfaux$posy == dfaux$posy[j]))>1){
        indvulg[dfaux$classe[j]] <- indvulg[dfaux$classe[j]] + 0.7
        dfaux$posx[j] <- (indvulg[dfaux$classe[j]]*primemove)%%denom
      }
      indvulg[dfaux$classe[j]] <- indvulg[dfaux$classe[j]] + 1
    }
    if (dfaux$kcore[j] == dfaux$kcore[j-1])
      dfaux$orderincore[j] <- dfaux$orderincore[j-1] + 1
    if ((dfaux$kcorenum[j] == 1) | (dfaux$kradius[j] == topradius))
      if ((maxcore>2) & (numkcores > 1))
        dfaux$vjust[j] <- (-saltovert*(1+0.2*(dfaux$kcorenum[j]-1)))
      else
        dfaux$vjust[j] <- 0
  }
  dfaux[dfaux$classe == "B",]$posx <- dfaux[dfaux$classe == "B",]$posx + offset
  dfaux$sizelabel <- 3+min(3,dfaux$kcorenum)

  polar_plot <- ggplot(dfaux, aes(x=posx,y=posy),legendTextFont=c(15, "bold.italic", "red")) +
    scale_size_area(max_size=scale_factor,name="k-degree") +
    scale_colour_manual(values = vcols,name="k-shell") +
    scale_fill_manual(values = vcols)+
    guides(col = guide_legend(override.aes = list(shape = 15, size = 8)),
           shape = guide_legend(override.aes = list(size = 8, colour = "black")))
  if (showtext == "yes"){
    polar_plot <- polar_plot+ geom_text(aes(size=0.005+0.05*normdegree,angle=0,
                                            colour = factor(kcorenum),
                                            label = num), alpha = alphal+0.1)
  }
  else{
    if (filln)
       polar_plot <- polar_plot + geom_point(aes(size=kdegree, colour = factor(kcorenum),
                                            shape = factor(symbol),fill=factor(kcorenum)),
                                            alpha = alphal, stroke = 2)+
                     scale_shape_manual(values=c(symbol_a,symbol_b),name="Guild",labels=rev(slabels))
    else
       polar_plot <- polar_plot + geom_point(aes(size=kdegree, colour = factor(kcorenum),
                                            shape = factor(symbol)),
                                            alpha = alphal, stroke = 2)+
                     scale_shape_manual(values=c(symbol_a,symbol_b),name="Guild",labels=rev(slabels))
  }

    polar_plot <- polar_plot +annotate(geom="text", x=dfaux$posx, y=dfaux$posy, label=dfaux$num,
                                      colour = factor(dfaux$kcol_label),
                                      size = dfaux$sizelabel,vjust = dfaux$vjust,
                                      hjust = 0.5, alpha = alphal, fontface="bold")
  polar_plot <- polar_plot + coord_polar(start = -pi/2) + labs(x = '', y = '')
  polar_plot <- polar_plot + scale_y_continuous(breaks=seq(min_radius,extreme),
                                                lim=c(min_radius, extreme),labels=seq(min_radius,extreme) )
  polar_plot <- polar_plot + guides(size="none", fill="none")
  polar_plot <- polar_plot + scale_x_continuous(breaks=seq(0, 2*pi, by=pi/2), lim=c(0,2*pi))
  polar_plot <- polar_plot+ theme_bw() + theme( axis.ticks.y = element_blank(),
                                               axis.ticks.x = element_blank(),
                                               panel.grid.major.x = element_blank(),
                                               panel.grid.minor.x = element_blank(),
                                               axis.text.y = element_blank(),
                                               panel.grid.major.y = element_line(size = 0.33,
                                                                                 linetype = 3,
                                                                                 alpha("darkolivegreen",0.8)),
                                               panel.grid.minor.y = element_blank(),
                                               panel.border = element_blank(),
                                               axis.text.x = element_blank(),
                                               legend.text = element_text(size=lsize_legend),
                                               legend.title = element_text(size=lsize_legend_title),
                                               plot.title = element_text(size=lsize_title,
                                                                         hjust = 0.5,lineheight=.8, face="bold")

                                        )
  ylab <- seq(0,extreme)
  pylab <- ylab
  pylab[2:length(pylab)] <- pylab[2:length(pylab)]-0.05
  ylab[1] <- "k-radius"
  xlab <- rep(pi-(guarda/2),length(pylab))
  dftext <- data.frame(xlab,ylab,pylab)
  dftext$fillcol <- maxcore
  polar_plot <- polar_plot + annotate(geom="text",x=xlab,y=pylab,label=ylab,size=4, color="gray50",
                                      lineheight=.8, alpha = 0.5)
  if (ptitle)
    polar_plot <- polar_plot + ggtitle(sprintf("Network %s", network_name)) +
    guides(row = guide_legend(nrow = 1))
  histo_kradius <- ggplot(dfaux, aes(kradius)) + geom_histogram(alpha = alpha_level,position='dodge',
                                                             binwidth=(extreme+1)/10,
                                                             color="white",fill = "forestgreen") +
    scale_x_continuous(breaks=seq(0, extreme+1, by=1), labels=seq(0, extreme+1, by=1), lim=c(0,extreme+1)) +
    theme_bw() +
    theme(legend.key = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_line(linetype = 3, color="gray90"),
          panel.grid.major.x = element_blank(),
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

  histo_core <- ggplot(dfaux, aes(x=kcorenum)) +
                geom_histogram(alpha = alpha_level, binwidth = 1,color="white",fill = "slategray3") + theme(legend.position = "none") +theme_bw() +
    scale_x_continuous(breaks=seq(0, maxcore+1, by=1), lim=c(0,maxcore+1)) +
    theme(panel.border = element_blank(),
          legend.key = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_line(linetype = 3, color="gray90"),
          panel.grid.major.x = element_blank(),
          axis.line = element_line(colour = "black",
                                   size = 0.5, linetype = "solid"),
          legend.text = element_text(size=lsize_legend),
          plot.title = element_text(size=lsize_title,lineheight=.8, face="bold"),
          axis.text.y = element_text(face="bold", size=lsize_axis),
          axis.text.x = element_text(angle = 0, face="bold", size=lsize_axis),
          axis.title.y  = element_text(face="bold", size=lsize_axis_title),
          axis.title.x = element_blank()
    ) +
    ggtitle("k-shell")+ ylab("Species")


  salto <- (1+ceiling(max(dfaux$kdegree))) %/% 8
  histo_kdegree <- ggplot(dfaux, aes(kdegree)) +
    geom_histogram(alpha = alpha_level,binwidth=(1+ceiling(max(dfaux$kdegree)))/8,
                   position='dodge',color="white",fill = "grey20") +
    theme_bw() +
    theme(panel.border = element_blank(),
          legend.key = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_line(linetype = 3, color="gray90"),
          panel.grid.major.x = element_blank(),
          axis.line = element_line(colour = "black",
                                   size = 0.5, linetype = "solid"),
          legend.text = element_text(size=lsize_legend),
          plot.title = element_text(size=lsize_title,lineheight=.8, face="bold"),
          axis.text = element_text(face="bold", size=lsize_axis),
          axis.title.x = element_blank(),
          axis.title.y  = element_text(face="bold", size=lsize_axis_title)
    )+
    ggtitle("k-degree")+ ylab("Species")
  calc_grafs <- list("polar_plot" = polar_plot, "histo_kradius" = histo_kradius, "histo_core" = histo_core,
                     "histo_kdegree" = histo_kdegree, "polar_file" = nfsal)
  return(calc_grafs)
}

#' Plotting a polar graph
#'
#' This function plots the polar graph of a bipartite network and the histograms of kshell
#' kradius and kdegree
#'
#' @param datadir the directory where the \code{filename} file is stored
#' @param filename the name of the file of the interaction matrix
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
#' @param print_title show title and network parameters
#' @param printable_labels range of labeled species
#' @param alpha_nodes fill transparency level
#' @param fill_nodes if set to FALSE nodes are transparent
#' @param max_kradius if bigger than 0 sets an upper limit different of the default. Useful for comparisons
#' @param progress auxiliar for interactive apps, do not modify
#' @export
#' @examples polar_graph("data/","M_PL_007.csv",plotsdir="grafresults/",print_to_file = TRUE)

polar_graph <- function( datadir, filename, plotsdir = "plot_results/polar/", print_to_file = FALSE, pshowtext = FALSE,
                         show_histograms = TRUE, glabels = c("Plant", "Pollinator"),
                         gshortened = c("pl","pol"),
                         lsize_title = 22, lsize_axis = 16,
                         lsize_legend = 16,lsize_axis_title = 16, lsize_legend_title = 16,
                         file_name_append = "", print_title = TRUE,
                         printable_labels = 0, fill_nodes = TRUE, alpha_nodes = 0.5,
                         max_kradius = 0, progress=NULL)
{

  # This assignment stores the call parameters in polar_argg as a list. This list is useful
  # to save plotting parameters for a future simulation

  polar_argg <- c(as.list(environment()))
  strip_isolated_nodes <- function()
  {
    lgrados <- igraph::degree(result_analysis$graph)
    if (sum(lgrados == 0) > 0)
      for (k in 1:length(lgrados))
      {
        if (lgrados[k] == 0){
          result_analysis$graph <<- delete_vertices(result_analysis$graph,names(lgrados[k]))
          if ( length(grep(sguild_b,names(lgrados[k]) )) >0 )
            result_analysis$num_guild_b <<- result_analysis$num_guild_b -1
          else
            result_analysis$num_guild_a <<- result_analysis$num_guild_a -1
        }
      }
  }

  red_name <- strsplit(filename,".csv")[[1]][1]
  sguild_a <<- gshortened[1]
  sguild_b <<- gshortened[2]
  slabels <<- glabels
  if (grepl("_SD_",filename) & (gshortened[1]=="pol") &  (gshortened[1]=="pl")){
    sguild_b = "disp"
    slabels <<- c("Plant", "Disperser")
  }

  if (!is.null(progress)) progress$inc(1/4, detail=strings$value("MESSAGE_POLAR_PROGRESS_ANALYZING_NETWORK"))
  result_analysis <- analyze_network(filename, directory = datadir, guild_a = sguild_a, guild_b = sguild_b, only_NODF = TRUE)
  strip_isolated_nodes()
  numlinks <- result_analysis$links
  an$result_analysis <<- result_analysis
  if (an$result_analysis$max_core == 1){
    msg = "Max core is 1. Polar plot only works if max core is bigger than 1"
    if (!is.null(progress))
      progress$inc(1/11, detail=strings$value(msg))
    else
      print(msg)
    return(an)
  }
  fsal=""
  if (print_to_file) {
    dir.create(plotsdir, showWarnings = FALSE)
    ppi <- 300
    if (file_name_append != "")
      ftname_append <- paste0("_",file_name_append)
    else
      ftname_append <- file_name_append
    fsal <- paste0(plotsdir,red_name,"_polar",ftname_append,".png")
    if (show_histograms)
      png(fsal, width=14*ppi, height=14*ppi, res=ppi)
    else
      png(fsal, width=10*ppi, height=10*ppi, res=ppi)
  }
  if (exists("zgg"))
    zgg$polar_file = fsal
  r <- paint_kdegree_kradius(result_analysis$graph, result_analysis$num_guild_a,result_analysis$num_guild_b,
                             lsize_title , lsize_axis, lsize_legend, lsize_axis_title , lsize_legend_title,
                             network_name = red_name, NODF = result_analysis$nested_values["NODF"],
                             Modularity =  result_analysis$modularity_measure,
                             MeanKradius = result_analysis$meandist, MeanKdegree = result_analysis$meankdegree,
                             showtext = pshowtext, fname_append = ftname_append,
                             printable_range = printable_labels, ptitle = print_title,
                             filln = fill_nodes, alphal = alpha_nodes, nfsal = fsal,
                             maxkradius = max_kradius, progress
                              )
    if (show_histograms)
      r$full_plot <- grid.arrange(r["polar_plot"][[1]], nrow=2, heights=c(4/5,1/5),
                   arrangeGrob(r["histo_kradius"][[1]],
                               r["histo_kdegree"][[1]],
                               r["histo_core"][[1]],ncol=3, nrow=1, widths=c(1/3,1/3,1/3)))
    else
      r$full_plot <- r["polar_plot"][[1]]
    if (print_to_file){
      print(r$full_plot)
      dev.off()
      }
  # }
  # Message for interactive apps.
  if (!is.null(progress)) {
    progress$inc(0, detail=strings$value("MESSAGE_POLAR_PROGRESS_DONE"))
  }
  r$polar_argg <- polar_argg
  r$result_analysis <- result_analysis
  r$network_name <- red_name
  return(r)
}

#polar_graph("../data/","M_SD_001.csv",print_to_file=TRUE, lsize_title = 24, lsize_axis = 18, lsize_legend = 18, lsize_axis_title = 18, lsize_legend_title = 20)