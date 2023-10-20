library(grid)
library(gridExtra)
library(bipartite)
library(igraph)
library(ggplot2)
library(rlang)
library(ggtext)

#' Plotting a ziggurat graph
#'
#' This function plots the ziggurat graph of a bipartite network. Configuration parameters and
#' results are stored in a global environment called zgg. This environment is not destroyed
#' after the function is executed, so the developer can store it in a configuration file, retrieve
#' network analysis variables and so on. If a new ziggurat_graph is called, zgg is destroyed and
#' created again for the new plot. Plotting options are explained in the user manual.
#'
#' @param datadir the name of the file of the interaction matrix
#' @param filename the file with the interaction matrix
#' @param print_to_file if set to FALSE the plot is displayed in the R session window
#' @param plotsdir the directory where the plot is stored
#' @param flip_results displays the graph in portrait configuration
#' @param aspect_ratio ziggurat plot default aspect ratio
#' @param alpha_level transparency for ziggurats' filling
#' @param color_guild_a default filling for nodes of guild_a
#' @param color_guild_b default filling for nodes of guild_b
#' @param color_link default links color
#' @param alpha_link link transparency
#' @param size_link width of the links
#' @param displace_y_b relative vertical displacement of guild_b inner ziggurats
#' @param displace_y_a relative vertical displacement of guild_a inner ziggurats
#' @param lsize_kcoremax nodes in kshell max label size
#' @param lsize_zig nodes in inner ziggurats label size
#' @param lsize_kcore1 labels of nodes in kshell 1
#' @param lsize_legend legend label size
#' @param lsize_kcorebox default kshell boxes label size
#' @param labels_color default label colors
#' @param height_box_y_expand expand inner ziggurat rectangles default height by this factor
#' @param kcore2tail_vertical_separation expand vertical of kshell 1 species linked to kshell 2 by this factor
#' @param kcore1tail_disttocore expand vertical separation of kshell 1 species from kshell max (guild_a, guild,b)
#' @param innertail_vertical_separation expand vertical separation of kshell species connected to khsell > 2 & < kshell max
#' @param factor_hop_x expand inner ziggurats horizontal distance
#' @param fattailjumphoriz displace kshell 1 species linked to leftmost kshell max species
#' @param fattailjumpvert idem for vertical position
#' @param coremax_triangle_width_factor expand khsell max rectangles width by this factor
#' @param coremax_triangle_height_factor expand khsell max rectangles height by this factor
#' @param paint_outsiders paint species not connected to giant component
#' @param displace_outside_component displace outsider species (horizontal, vertical)
#' @param outsiders_separation_expand multiply by this factor outsiders' separation
#' @param outsiders_legend_expand displace outsiders legend
#' @param specialistskcore2_horizontal_dist_rootleaf_expand expand horizontal distance of specialist tail root node connected to kshell 2
#' @param specialistskcore2_vertical_dist_rootleaf_expand expand vertical distance of specialist tails connected to kshell 2
#' @param specialists_boxes_separation_count specialist species boxes separation count
#' @param root_specialist_expand expand root specialist distances of tails connected to kshell <> 2
#' @param hide_plot_border hide border around the plot
#' @param rescale_plot_area full plot area rescaling (horizontal, vertical)
#' @param kcore1specialists_leafs_vertical_separation expand vertical separation of specialist tails connected to kshell 1 species
#' @param corebox_border_size width of kshell boxes
#' @param kcore_species_name_display display species names of  shells listed in this vector
#' @param kcore_species_name_break allow new lines in species names of  shells listed in this vector
#' @param shorten_species_name number of characters of species name to display
#' @param exclude_species_number do not include species number in species
#' @param label_strguilda string labels of guild a
#' @param label_strguildb string labels of guild b
#' @param landscape_plot paper landscape configuration
#' @param backg_color plot background color
#' @param show_title show plot title
#' @param use_spline use splines to draw links
#' @param spline_points number of points for each spline
#' @param file_name_append a label that the user may append to the plot file name for convenience
#' @param svg_scale_factor only for interactive apps, do not modify
#' @param weighted_links function to add link weight: 'none', 'log10' , 'ln', 'sqrt'
#' @param square_nodes_size_scale scale nodes area of kcore1 and outsiders
#' @param move_all_SVG_up move up all the SVG plot by this fraction, useful to crop upper white space
#' @param progress only for interactive apps, do not modifiy
#' @export
#' @examples ziggurat_graph("data/","M_PL_001.csv",plotsdir="grafresults/",print_to_file = TRUE)

ziggurat_graph <- function(datadir,filename,
                           paintlinks = TRUE, print_to_file = FALSE, plotsdir ="plot_results/ziggurat/", flip_results = FALSE,
                           aspect_ratio = 1,
                           alpha_level = 0.2, color_guild_a = c("#4169E1","#00008B"), color_guild_b = c("#F08080","#FF0000"),
                           color_link = "slategray3", alpha_link = 0.5, size_link = 0.5,
                           displace_y_b = rep(0,20),
                           displace_y_a = rep(0,20),
                           lsize_kcoremax = 3.5, lsize_zig = 3, lsize_kcore1 = 2.5, lsize_legend = 4, lsize_core_box = 2.5,
                           labels_color = c(),
                           height_box_y_expand = 1, kcore2tail_vertical_separation = 1,  kcore1tail_disttocore = c(1,1),
                           innertail_vertical_separation = 1,
                           factor_hop_x = 1, fattailjumphoriz = c(1,1), fattailjumpvert = c(1,1),
                           coremax_triangle_height_factor = 1, coremax_triangle_width_factor = 1,
                           paint_outsiders = TRUE, displace_outside_component = c(0,0),
                           outsiders_separation_expand = 1, outsiders_legend_expand = 1,
                           specialistskcore2_horizontal_dist_rootleaf_expand = 1,
                           specialistskcore2_vertical_dist_rootleaf_expand = 0, specialists_boxes_separation_count = 1,
                           root_specialist_expand = c(1,1), hide_plot_border = TRUE, rescale_plot_area = c(1,1),
                           kcore1specialists_leafs_vertical_separation = 1, corebox_border_size = 0.2,
                           kcore_species_name_display = c(), kcore_species_name_break = c(),
                           shorten_species_name = 0, exclude_species_number = FALSE, label_strguilda = "",
                           label_strguildb = "", landscape_plot = TRUE,
                           backg_color = "white", show_title = TRUE, use_spline =TRUE, spline_points = 10,
                           file_name_append = "", svg_scale_factor= 10, weighted_links = "none",
                           square_nodes_size_scale = 1, move_all_SVG_up = 0,
                           progress=NULL
)
{
  # This assignment stores the call parameters in ziggurat_argg as a list. This list is useful
  # to save plotting parameters for a future simulation

  ziggurat_argg <- c(as.list(environment()))

  # Create global environment
  zgg <<- new.env()
  if (!is.null(progress)) progress$inc(1/11, detail=strings$value("MESSAGE_ZIGGURAT_PROGRESS_ANALYZING_NETWORK"))
  # Analyze network
  f <- read_and_analyze(datadir,filename,label_strguilda, label_strguildb)
  zgg$result_analysis <- f["result_analysis"][[1]]
  zgg$str_guild_a <- f["str_guild_a"][[1]]
  zgg$str_guild_b <- f["str_guild_b"][[1]]
  zgg$name_guild_a <- f["name_guild_a"][[1]]
  zgg$name_guild_b <- f["name_guild_b"][[1]]
  zgg$network_name <- f["network_name"][[1]]
  # Exit if kcore max == 1
  if (zgg$result_analysis$max_core == 1){
    msg = "Max core is 1. Ziggurat plot only works if max core is bigger than 1"
    if (!is.null(progress))
      progress$inc(1/11, detail=strings$value(msg))
    else
      print(msg)
    return(zgg)
  }
  # Copy input parameters to the zgg environment
  def_configuration(paintlinks, print_to_file, plotsdir, flip_results, aspect_ratio,
                    alpha_level, color_guild_a, color_guild_b,
                    color_link, alpha_link, size_link,
                    displace_y_b, displace_y_a, lsize_kcoremax, lsize_zig, lsize_kcore1,
                    lsize_legend, lsize_core_box, labels_color,
                    height_box_y_expand, kcore2tail_vertical_separation,  kcore1tail_disttocore,
                    innertail_vertical_separation,
                    factor_hop_x, fattailjumphoriz, fattailjumpvert,
                    coremax_triangle_height_factor, coremax_triangle_width_factor,
                    paint_outsiders, displace_outside_component,
                    outsiders_separation_expand, outsiders_legend_expand, specialistskcore2_horizontal_dist_rootleaf_expand,
                    specialistskcore2_vertical_dist_rootleaf_expand , specialists_boxes_separation_count,
                    root_specialist_expand, hide_plot_border, rescale_plot_area,kcore1specialists_leafs_vertical_separation,
                    corebox_border_size, kcore_species_name_display,kcore_species_name_break,shorten_species_name,exclude_species_number,
                    label_strguilda, label_strguildb, landscape_plot, backg_color, show_title,
                    use_spline, spline_points, file_name_append, svg_scale_factor, weighted_links,
                    square_nodes_size_scale, move_all_SVG_up, progress
  )
  # Removes nodes without any tie. This is not usual in input files but happens
  # when performing destruction simulations
  strip_isolated_nodes()
  init_working_values()
  draw_ziggurat_plot(svg_scale_factor, progress)
  # Copy input parameters as a string for reroducibility
  zgg$ziggurat_argg <- ziggurat_argg
  return(zgg)
}


# Labels of square nodes: tails, specialist chains and outsiders
gen_sq_label <- function(nodes, joinchars = "\n", is_guild_a = TRUE)
{
  # If kcore1 nodes name are displayed
  dispname <- is.element(1,zgg$kcore_species_name_display)
  if (dispname)
    if (is_guild_a)
      nspec <- colnames(zgg$result_analysis$matrix)
    else
      nspec <- rownames(zgg$result_analysis$matrix)
  nnodes <- length(nodes)
  lrow <- round(sqrt(nnodes))
  ssal <- ""
  for (i in 1:nnodes)
  {
    if (dispname)
      ssal <- paste(ssal,nspec[as.integer(nodes[i])])
    else
      ssal <- paste(ssal,nodes[i])
    if ((i %% lrow == 0) & (nnodes > 1) & (i<nnodes))
      ssal <- gsub("  "," ",paste(ssal,joinchars))
  }
  return(ssal)
}

# Create the label species. May be complex if user chooses to display
# the binomial name
create_label_species <- function(strent,newline = FALSE){
  strchar <- ifelse(newline,"\n","")
  pieces <- unlist(strsplit(unlist(strent)," "))
  if (is.na(pieces[2]))
    pieces[2] = ""
  if (zgg$shorten_species_name>0){
    pieces[1] = paste0(substr(pieces[1],1,zgg$shorten_species_name),".")
    if (nchar(pieces[2])>2)
      pieces[2] = paste0(substr(pieces[2],1,zgg$shorten_species_name),".")
  }
  if (length(pieces)>2)
    strsal <- paste(pieces[1],strchar,"XX")
  else
    strsal <- paste(pieces[1],strchar,pieces[2])
  return(strsal)
}

# Decides length and rotation of labels
name_species_preprocess <- function (kcore, list_dfs, kcore_species_name_display,
                                     kcore_species_name_break) {
  if (is.element(kcore,kcore_species_name_display)){
    if (!zgg$flip_results)
      kcoremaxlabel_angle <- 90
    else
      kcoremaxlabel_angle <- 0
    labelszig <- rep("",nrow(list_dfs))
    pnewline <- is.element(kcore,kcore_species_name_break)
    for (j in 1:length(list_dfs$name_species)){
      labelszig[j] <- create_label_species(list_dfs$name_species[j],newline=pnewline)
      if (!zgg$exclude_species_number)
        labelszig[j] <- paste(list_dfs$label[j],labelszig[j])
    }
  } else {
    kcoremaxlabel_angle <- 0
    labelszig <- list_dfs$label
  }
  calc_values <- list("kcoremaxlabel_angle" = kcoremaxlabel_angle, "labelszig" = labelszig)
  return(calc_values)
}

# Draws a square, both in ggplot2 and SVG flavours
draw_square<- function(idPrefix, grafo,svg,basex,basey,side,fillcolor,alphasq,labelcolor,
                       langle,hjust,vjust,slabel,lbsize = zgg$labels_size,
                       inverse="no",adjustoxy = "yes", edgescolor="transparent")
{
  x1 <- c(basex)
  x2 <- c(basex+side)
  y1 <- c(basey)
  y2 <- c(basey+side/zgg$aspect_ratio)
  ds <- data.frame(x1, x2, y1, y2, fillcolor)
  signo <- 1
  if (inverse == "yes")
  {
    ds$y1 <- -(ds$y1)
    ds$y2 <- -(ds$y2)
    signo <- -1
  }

  p <- grafo + geom_rect(data=ds, linewidth=0.01,
                         mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2),
                         fill = fillcolor, alpha = alphasq, color="transparent")
  pxx <- x1+0.05*(x2-x1)
  pyy <- signo*(y1+(y2-y1)/2)
  p <- p +annotate(geom="text", x=pxx, y=pyy, label=slabel,
                   colour = labelcolor, size=lbsize, hjust = hjust,
                   vjust = vjust, angle = langle)
  svg$rect(idPrefix=idPrefix, data=ds, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2),
           fill = fillcolor, alpha=alphasq, size=0.5, color="transparent")
  svg$text(idPrefix=idPrefix, data=data.frame(x=x1+(x2-x1)/2,y=pyy),
                                              mapping=aes(x=x, y=y),color=labelcolor, label=slabel, size=lbsize, angle=langle)
  calc_vals <- list("p" = p, "svg" = svg)
  return(calc_vals)
}

# Draw a rectangle. Nodes of inner ziggurats and coremax are rectangles
draw_rectangle<- function(idPrefix,basex,basey,widthx,widthy,grafo,svg,bordercolor,fillcolor,palpha,slabel,inverse="no",sizelabel=3, bordersize =0.5 )
{
  x1 <- c(basex)
  x2 <- c(basex+widthx)
  y1 <- c(basey)
  y2 <- c(basey+widthy)
  ds <- data.frame(x1, x2, y1, y2, fillcolor)
  signo <- 1
  if (inverse == "yes")
  {
    ds$y1 <- -(ds$y1)
    ds$y2 <- -(ds$y2)
    signo <- -1
  }
  if (bordersize > 0)
    p <- grafo + geom_rect(data=ds,
                         mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2),
                         fill = fillcolor, alpha =palpha, color=bordercolor, linewidth = bordersize, linetype = 3)
  else
    p <- grafo
  p <- p +annotate(geom="text", x=x1+(x2-x1)/8, y=signo*(y1+(y2-y1)/2), label=slabel,
                   colour = fillcolor, size=sizelabel, hjust = 0)
  svg$rect(idPrefix=idPrefix, data=ds, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2),
           fill=fillcolor, alpha=palpha, color=bordercolor, size=bordersize, linetype=3)
  svg$text(idPrefix=idPrefix, data=data.frame(x=c(x1+(x2-x1)/8), y=c(signo*(y1+(y2-y1)/2))),
           mapping=aes(x=x, y=y), color=fillcolor, label=slabel, size=sizelabel)

  calc_vals <- list("p" = p, "svg" = svg)
  return(calc_vals)
}

# Draws (optional) enclosing shell boxes
draw_core_box <- function(grafo, svg, kcore)
{
  marginy <- ifelse((zgg$df_cores$num_species_guild_a[2]+zgg$df_cores$num_species_guild_b[2]>4),1.2*zgg$height_y,0.7*zgg$height_y)
  marginx <- 1.5*zgg$lado
  if (kcore<zgg$kcoremax)
  {
    x_inf <- min(zgg$list_dfs_b[[kcore]]$x1,zgg$list_dfs_a[[kcore]]$x1) - marginx
    widthx <- (max(zgg$list_dfs_b[[kcore]]$x2,zgg$list_dfs_a[[kcore]]$x2 ) - x_inf) + marginx
    y_inf <- ifelse(length(zgg$list_dfs_a[[kcore]])>0, min(zgg$list_dfs_a[[kcore]]$y2,zgg$list_dfs_a[[kcore]]$y1) - marginy,
                    min(zgg$list_dfs_b[[kcore]]$y2,zgg$list_dfs_b[[kcore]]$y1) - marginy)
    widthy <- ifelse(length(zgg$list_dfs_b[[kcore]])>0, max(zgg$list_dfs_b[[kcore]]$y2) - y_inf + (1+0.45*zgg$kcoremax)*marginy,
                     widthy <- max(zgg$list_dfs_a[[kcore]]$y1) - y_inf + (1+0.45*zgg$kcoremax)*marginy)
  }
  else{
    x_inf <- min(zgg$list_dfs_b[[kcore]]$x1,zgg$list_dfs_a[[kcore]]$x1) - 3*marginx
    widthx <- (max(zgg$list_dfs_b[[kcore]]$x2,zgg$list_dfs_a[[kcore]]$x2 ) - x_inf) + marginx
    y_inf <- min(zgg$list_dfs_a[[kcore]]$y2,zgg$list_dfs_b[[kcore]]$y2) - marginy
    widthy <- max(zgg$list_dfs_a[[kcore]]$y2) - y_inf + 2*marginy
  }
  divcolor <- zgg$corecols[kcore]
  f <- draw_rectangle(paste0("kcore", kcore),x_inf,y_inf,widthx,widthy,grafo,svg, divcolor,"transparent",0.3,"",inverse="no",sizelabel = zgg$labels_size, bordersize = zgg$corebox_border_size)
  p <- f["p"][[1]]
  svg <- f["svg"][[1]]
  if (kcore == zgg$kcoremax){
    position_x_text <- x_inf+1.5*marginx
    corelabel <- paste0(kcore,"-shell")
  }
  else{
    position_x_text <- x_inf-marginx+widthx/2
    corelabel <- paste0(kcore,"-shell")
  }
  if (!is.null(nrow(zgg$list_dfs_b[[kcore]])))
    position_y_text <- y_inf+widthy - marginy
  else
    position_y_text <- y_inf - marginy
  zgg$max_position_y_text_core <- max(zgg$max_position_y_text_core,position_y_text)
  if (kcore != zgg$kcoremax){
    px <- position_x_text
    py <- position_y_text
    pangle <- 0
    phjust <- 0
  }
  else {
    px <- position_x_text+1.2*marginx/2
    py <- 0
    phjust <- 0.5
    ifelse (zgg$flip_results,pangle<-0,pangle <- 90)
  }
  # Print label only if size is not 0
  if (zgg$lsize_core_box>0)
    p <- p +annotate(geom="text", x=px, y=py, label=corelabel, colour = divcolor,  fontface="italic",
                   size=zgg$lsize_core_box, hjust = phjust, vjust = 0, angle = pangle)
  svg$text(idPrefix=paste0("kcore", kcore), data=data.frame(x=c(px), y=c(py)),
           mapping=aes(x=x, y=y), color=divcolor, label=corelabel, size=zgg$lsize_core_box, angle=pangle)

  calc_vals <- list("p" = p, "svg" = svg, "max_position_y_text_core" = zgg$max_position_y_text_core)
  return(calc_vals)
}

# Computes the weight of links depending on the grouping option chosen by the user
get_link_weights <- function(matrixweight)
{
  if (zgg$weighted_links == "none")
    return(1)
  if (zgg$weighted_links == "log10")
    return(1+log10(matrixweight))
  if (zgg$weighted_links == "ln")
    return(1+log(matrixweight))
  if (zgg$weighted_links == "sqrt")
    return(sqrt(matrixweight))
}

# Draw tail, species or set of species of 1-shell connected to higher k-index species
draw_tail <- function(idPrefix, p,svg,fat_tail,lado,color,sqlabel,basex,basey,gap,
                      lxx2=0,lyy2=0,sqinverse = "no",
                      position = "West", background = "no",
                      first_leaf = "yes", spline = "no",
                      psize = zgg$lsize_kcore1, is_guild_a = TRUE, wlink=1)
{
  adjust <- "yes"
  lvjust <- 0
  lhjust <- 0
  langle <- 0
  ecolor <- "transparent"
  bgcolor <- color
  labelcolor <- ifelse(length(zgg$labels_color)>0,zgg$labels_color[2-as.numeric(is_guild_a)],color)
  palpha <- zgg$alpha_link
  sidex <- lado*(0.5+sqrt(nrow(fat_tail)))
  sidex <- sidex * sqrt(zgg$square_nodes_size_scale)
  paintsidex <- sidex
  signo <- 1
  yy <- abs(basey)
  plxx2 <- lxx2
  plyy2 <- lyy2
  # Lower half plane of the ziggurat
  if (sqinverse=="yes")
    signo <- -1
  # Tails connected eastwards to inner ziggurats
  if (position == "West"){
    xx <- basex-gap
    posxx1 <- xx+sidex
    posyy1 = signo*(yy)+signo*(0.5*sidex/(zgg$aspect_ratio))
  }
  # Fat tails, or group of species linked to the highest kdegree species in max shell
  else if (position == "East"){
    gap <- zgg$hop_x/2
    xx <- basex+gap
    posxx1 <- xx
    posyy1 = signo*(yy+sidex/(2*zgg$aspect_ratio))
  }
  # Tails connected to other nodes of the max shell
  else if ((position == "North") |(position == "South")) {
    xx <- basex-sidex/2
    posxx1 <- xx+sidex/2
    posyy1 = signo*(yy)
  }

  if (background == "no")
  {
    ecolor <- "transparent"
    if (zgg$alpha_level != 1)
      palpha <- max(zgg$alpha_level-0.09,0)

    rot_angle <-30
    if (position == "North")
    {
      langle <- rot_angle
    }
    if (position == "South")
    {
      langle <- -rot_angle
      lvjust <- 0
    }
    if (position == "West"){
      adjust <- "yes"
    }
  }
  if ((zgg$flip_results) & (langle == 0) & (position!="West") )
    langle <- langle + 70
  # Draw square node
  f <- draw_square(idPrefix, p,svg, xx,yy,paintsidex*sqrt(zgg$square_nodes_size_scale),
                   bgcolor,palpha,labelcolor,langle,lhjust,lvjust,
                   slabel=sqlabel,lbsize = psize,inverse = sqinverse,
                   adjustoxy = adjust, edgescolor = ecolor)
  p <- f["p"][[1]]
  svg <- f["svg"][[1]]
  # Add tail link
  if (zgg$paintlinks)
    add_link(xx1=posxx1, xx2 = plxx2,
                   yy1 = posyy1, yy2 = plyy2,
                   slink = zgg$size_link*wlink, clink = c(zgg$color_link),
                   alpha_l = zgg$alpha_link , spline= spline)
  calc_vals <- list("p" = p, "svg" = svg, "sidex" = sidex, "xx" = posxx1, "yy" = posyy1)
  return(calc_vals)
}


draw_edge_tails <- function(p,svg,point_x,point_y,kcoreother,long_tail,list_dfs,color_guild, inverse = "no",
                            vertical = "yes", orientation = "East", revanddrop = "no",
                            pbackground = "yes", joinchars = "\n", tspline = "no", is_guild_a = TRUE, wlink = 1)
{

  if (orientation == "West")
    point_y <- point_y - zgg$gap
  rxx <- point_x
  ryy <- point_y
  zgg$joinstr <- joinchars
  signo <- 1
  if (inverse == "yes")
    signo <- -1
  list_spec <- list_dfs[[kcoreother]]$label
  if (revanddrop == "yes")
    list_spec <- rev(list_spec)[1:length(list_spec)-1]
  if (orientation == "East")
    list_spec <- rev(list_spec)
  llspec <- length(list_spec)
  m <- 0
  separacion <- 0.035*zgg$tot_width
  for (i in list_spec)
  {
    conn_species <- which(long_tail$partner == i)
    if (length(conn_species)>0)
    {
      little_tail <- long_tail[long_tail$partner == i,]
      if ((orientation == "East") | (orientation == "West"))
      {
        data_row <- list_dfs[[kcoreother]][which(list_dfs[[kcoreother]]$label == i),]
        xx2 <- data_row$x2
        yy2 <- (data_row$y1 + data_row$y2)/2
      }
      else# if (orientation = "North")
      {
        rxx <- point_x*zgg$kcore1tail_disttocore[1]
        ryy <- point_y*0.9*zgg$kcore1tail_disttocore[2]
        data_row <- list_dfs[[kcoreother]][which(list_dfs[[kcoreother]]$label == i),]
        xx2 <- (data_row$x2+data_row$x1)/2
        yy2 <- data_row$y2
      }
      if (kcoreother == 2) {
        yy2 <- ifelse (zgg$kcoremax > 2, (data_row$y1 + data_row$y2)/2, data_row$y2)
        point_y <- (zgg$kcore2tail_vertical_separation-1)*sign(point_y)*zgg$height_y + point_y
        ryy <- point_y
      }
      else if (kcoreother<zgg$kcoremax) {
        point_y <- (zgg$innertail_vertical_separation-1)*sign(point_y)*zgg$height_y + point_y
        ryy <- point_y
      }
      tailweight <- 0
      for (h in 1:nrow(little_tail))
        if (is_guild_a)
          tailweight <- tailweight + zgg$result_analysis$matrix[as.numeric(little_tail$partner[h]),
                                                                        as.numeric(little_tail$orph[h])]
        else
          tailweight <- tailweight + zgg$result_analysis$matrix[as.numeric(little_tail$orph[h]),
                                                                as.numeric(little_tail$partner[h])]
      little_tail$weightlink <- get_link_weights(tailweight)
      v<- draw_tail(paste0(ifelse(is_guild_a, "edge-kcore1-a-", "edge-kcore1-b-"), i),
                    p,svg,little_tail,zgg$lado,color_guild[2],
                    gen_sq_label(little_tail$orph,joinchars = " ", is_guild_a = is_guild_a),
                    rxx,ryy,zgg$gap,lxx2 = xx2,
                    lyy2 = yy2, sqinverse = inverse, position = orientation,
                    background = pbackground, spline = tspline, psize = zgg$lsize_kcore1,
                    is_guild_a = is_guild_a, wlink = little_tail$weightlink[1])
      p <- v["p"][[1]]
      svg <- v["svg"][[1]]
      rxx <- v["xx"][[1]]
      ryy <- v["yy"][[1]]
      if (vertical == "yes"){
        salto <- v["sidex"][[1]]/zgg$aspect_ratio
        point_y <- point_y + 1.4*signo*salto
        rxx <- point_x
      }

      # tails connected to kcoremax except first species

      else{
        if (orientation == "West")
          salto <- 0
        else
          salto <- 0.4*v["sidex"][[1]]/zgg$aspect_ratio
        point_x <- point_x - separacion - v["sidex"][[1]]
        point_y <- point_y - 1.4*signo*salto
        ryy <- point_y
        rxx <- point_x
      }
    }
    m <- m +1
  }
  calc_vals <- list("p" = p, "svg" = svg, "lastx" = rxx, "lasty" = ryy)
  return(calc_vals)
}


# Species disconnected of the Giant Component, called 'outsiders'
# Compute coordinates
conf_outsiders <- function(outsiders,basex,basey,sidex,fillcolor,strguild)
{
  x1 <- c()
  x2 <- c()
  y1 <- c()
  y2 <- c()
  r <- c()
  col_row <- c()
  numboxes <- length(outsiders)
  pbasex <- basex
  xstep <- 2*sidex*zgg$outsiders_separation_expand
  xsep <- 2.5*zgg$outsiders_separation_expand
  ysep <- xsep
  for (j in (1:numboxes))
  {
    x1 <- c(x1, pbasex+(j*xsep*xstep))
    x2 <- c(x2, x1[j]+xstep)
    y1 <- c(y1, basey-ysep*xstep/zgg$aspect_ratio)
    y2 <- c(y2, y1[j]-xstep/zgg$aspect_ratio)
    r <- c(r,j)
    col_row <- c(col_row,fillcolor)
  }
  d1 <- data.frame(x1, x2, y1, y2, r, col_row)
  d1$label <- ""
  for (i in 1:length(outsiders))
    d1[i,]$label <- strsplit(outsiders[i],strguild)[[1]][2]
  return(d1)
}

# Draw outsider square nodes
draw_sq_outsiders <- function(idPrefix, p,svg,dfo,paintsidex,alpha_level,lsize,is_guild_a = TRUE)
{
  if (length(zgg$labels_color)>0)
    labelscolor <- rep(zgg$labels_color[2-as.numeric(is_guild_a)],nrow(dfo))
  else
    labelscolor <- dfo$col_row
  p <- p + geom_rect(data=dfo, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2),
                     fill = dfo$col_row, alpha = alpha_level,color="transparent") +
       geom_text(data=dfo, aes(x=(x2+x1)/2, y= (y2+y1)/2), color=labelscolor,
              label = dfo$label, size=lsize, vjust=0.5)
  svg$rect(idPrefix=idPrefix, data=dfo, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill=dfo$col_row, alpha=alpha_level,
           color="transparent", size=0.5)
  svg$text(idPrefix=idPrefix, data=dfo, mapping=aes(x=(x2+x1)/2, y= (y2+y1)/2), color=labelscolor, label=dfo$label, size=lsize)
  calc_vals <- list("p" = p, "svg" = svg)
  return(calc_vals)
}

# Full outsiders management procedure
handle_outsiders <- function(p,svg,outsiders,df_chains) {
  if (length(zgg$outsider)>0){
    paintsidex <- 2*zgg$height_y*zgg$aspect_ratio
    paintsidex <- paintsidex * sqrt(zgg$square_nodes_size_scale)
    zgg$outsiders_a <- zgg$outsider$name[grep(zgg$str_guild_a,zgg$outsider$name)]
    zgg$outsiders_b <- zgg$outsider$name[grep(zgg$str_guild_b,zgg$outsider$name)]
    pox <- -(zgg$hop_x/4)+ zgg$tot_width * (zgg$displace_outside_component[1])
    poy <- min(-zgg$last_ytail_b[!is.na(zgg$last_ytail_b)]-4*zgg$lado,df_chains$y1) * (1+zgg$displace_outside_component[2])
    dfo_a <- conf_outsiders(zgg$outsiders_a,pox,poy,
                            zgg$lado*sqrt(zgg$square_nodes_size_scale),zgg$color_guild_a[2],zgg$str_guild_a)
    guild_sep <- poy-max(1,length(zgg$outsider)/10)*6*zgg$lado*sqrt(zgg$square_nodes_size_scale)*zgg$outsiders_separation_expand/zgg$aspect_ratio
    dfo_b <- conf_outsiders(zgg$outsiders_b,pox,guild_sep,
                            zgg$lado*sqrt(zgg$square_nodes_size_scale),zgg$color_guild_b[2],zgg$str_guild_b)
    f <- draw_sq_outsiders("outsiders-kcore1-a",p,svg,dfo_a,paintsidex,zgg$alpha_level,zgg$lsize_kcore1)
    p <- f["p"][[1]]
    svg <- f["svg"][[1]]
    f <- draw_sq_outsiders("outsiders-kcore1-b",p,svg,dfo_b,paintsidex,zgg$alpha_level,zgg$lsize_kcore1, is_guild_a = FALSE)
    p <- f["p"][[1]]
    svg <- f["svg"][[1]]
    for (j in 1:nrow(dfo_a))
    {
      mtxa <- zgg$mtxlinks[which(zgg$mtxlinks$guild_a == paste0(zgg$str_guild_a,dfo_a[j,]$label)),]
      for (i in 1:nrow(dfo_b))
      {
        if (sum(as.character(mtxa$guild_b) == paste0(zgg$str_guild_b,dfo_b[i,]$label))>0)
        {
          bend_line = "no"
          link <- data.frame(x1=c(dfo_a[j,]$x1 + (dfo_a[j,]$x2-dfo_a[j,]$x1)/2),
                             x2 = c(dfo_b[i,]$x1 +(dfo_b[i,]$x2-dfo_b[i,]$x1)/2),
                             y1 = c(dfo_a[j,]$y2),  y2 = c(dfo_b[i,]$y1) )
          lcolor = "orange"
          tailweight <- get_link_weights(zgg$result_analysis$matrix[as.numeric(dfo_b[i,]$label),
                                                                                       as.numeric(dfo_a[j,]$label)])
          add_link(xx1=link$x1, xx2 = link$x2,
                         yy1 = link$y1, yy2 = link$y2,
                         slink = zgg$size_link*tailweight, clink = c(zgg$color_link),
                         alpha_l = zgg$alpha_link , spline = bend_line)
        }
      }
    }
    margin <- zgg$height_y
    x_inf <- min(dfo_a$x1,dfo_b$x1) - 1.5*margin
    widthx <- max(dfo_a$x2,dfo_b$x2) - x_inf + margin
    y_inf <- min(dfo_a$y2,dfo_b$y2) - 2*margin/zgg$aspect_ratio
    widthy <- max(dfo_a$y2,dfo_b$y2) - y_inf + 2*margin/zgg$aspect_ratio
    divcolor <- "grey50"
    position_x_text <- x_inf+20*zgg$outsiders_legend_expand
    corelabel <- paste("Outside the giant component")
    position_y_text <- y_inf + margin/zgg$aspect_ratio + (0.9+0.2*zgg$outsiders_legend_expand)*widthy
    px <- position_x_text
    py <- position_y_text
    if (zgg$flip_results){
      px <- x_inf + widthx + margin * zgg$outsiders_legend_expand
      py <- y_inf
    }
    p <- p +annotate(geom="text", x=px, y=py, label=corelabel, colour = divcolor,
                     size=zgg$lsize_core_box, hjust = 0, vjust = 0, angle = 0)

    svg$text("corelabel", data=data.frame(x=c(px), y=c(py)), mapping=aes(x=x, y=y), color=divcolor, label=corelabel, size=zgg$lsize_core_box)
  }
  calc_vals <- list("p" = p, "svg" = svg)
  return(calc_vals)
}

# Draw one triangle of the innermost shell
draw_coremax_triangle <- function(basex,topx,basey,topy,numboxes,fillcolor,strlabels,
                                  igraphnet,strguild,orderby = "None")
{
  x1 <- c()
  x2 <- c()
  y1 <- c()
  y2 <- c()
  r <- c()
  kdegree <- c()
  kradius <- c()
  col_row <- c()
  name_species <- c()
  pbasex <- zgg$coremax_triangle_width_factor*( basex - (numboxes %/%8) * abs(topx-basex)/3)
  xstep <- (topx-pbasex)*1/numboxes
  ptopy <- topy * zgg$coremax_triangle_height_factor
  ystep <- (ptopy-basey)*0.7/numboxes
  for (j in (1:numboxes))
  {
    x1 <- c(x1, pbasex+(j-1)*xstep)
    x2 <- c(x2, x1[j]+0.9*xstep)
    y1 <- c(y1, basey)
    y2 <- c(y2, ptopy-(j-1)*ystep)
    r <- c(r,j)
    col_row <- c(col_row,fillcolor[1+j%%2])
    kdegree <- c(kdegree,0)
    kradius <- c(kradius,1)
    name_species <- c(name_species,"")
  }
  d1 <- data.frame(x1, x2, y1, y2, r, col_row, kdegree, kradius, name_species, stringsAsFactors=FALSE)
  d1$label <- strlabels
  for (i in 1:nrow(d1)){
    d1[i,]$kdegree <- igraphnet[paste0(strguild,d1[i,]$label)]$kdegree
    d1[i,]$kradius <- igraphnet[paste0(strguild,d1[i,]$label)]$kradius
    d1[i,]$name_species <- igraphnet[paste0(strguild,d1[i,]$label)]$name_species
  }

  if (orderby == "kradius"){
    ordvector <- order(d1$kradius)
    d1$label <- d1[ordvector,]$label
    d1$kradius <- d1[ordvector,]$kradius
    d1$kdegree <- d1[ordvector,]$kdegree
    d1$name_species <- d1[ordvector,]$name_species
  }
  else if (orderby == "kdegree"){
    ordvector <- rev(order(d1$kdegree))
    d1$label <- d1[ordvector,]$label
    d1$kradius <- d1[ordvector,]$kradius
    d1$kdegree <- d1[ordvector,]$kdegree
    d1$name_species <- d1[ordvector,]$name_species
  }
  return(d1)
}

# This function adds the following information for species of kcore 1 in lists_dfs_x: label, name_species, kdegree, kradius
conf_kcore1_info <- function(strguild)
{
  auxlistdf <- data.frame(x1=NA,x2=NA,y1=NA,y2=NA,r=NA,col_row=NA,kdegree=NA,kradius=NA,name_species=NA,label=NA)
  retlistdf <- data.frame(x1=c(),x2=c(),y1=c(),y2=c(),r=c(),col_row=c(),kdegree=c(),kradius=c(),name_species=c(),label=c())
  num_s <- zgg$df_cores[1,]$num_species_guild_a
  if (strguild == zgg$str_guild_a)
    listspecies = zgg$df_cores[1,]$species_guild_a
  else
    listspecies = zgg$df_cores[1,]$species_guild_b
  for (j in listspecies[[1]])
  {
    ind <- paste0(strguild,j)
    auxlistdf$label <- j
    auxlistdf$kdegree <-  V(zgg$result_analysis$graph)[ind]$kdegree
    auxlistdf$kradius <-  V(zgg$result_analysis$graph)[ind]$kradius
    auxlistdf$name_species <-  V(zgg$result_analysis$graph)[ind]$name_species
    retlistdf <- rbind(retlistdf,auxlistdf)
  }
  return(retlistdf)
}

# Similar to conf_kcore1_info for species outside the giant component
conf_outsiders_info <- function(strguild)
{
  auxlistdf <- data.frame(x1=NA,x2=NA,y1=NA,y2=NA,r=NA,col_row=NA,kdegree=NA,kradius=NA,name_species=NA,label=NA)
  retlistdf <- data.frame(x1=c(),x2=c(),y1=c(),y2=c(),r=c(),col_row=c(),kdegree=c(),kradius=c(),name_species=c(),label=c())
  listspecies <- NULL
  num_s <- zgg$cores[1,]$num_species_guild_a
  if (strguild == zgg$str_guild_a)
      listspecies = zgg$outsiders_a
  else
      listspecies = zgg$outsiders_b
  for (j in listspecies)
  {
    auxlistdf$label <- gsub(strguild,"",j)
    auxlistdf$kdegree <-  V(zgg$result_analysis$graph)[j]$kdegree
    auxlistdf$kradius <-  V(zgg$result_analysis$graph)[j]$kradius
    auxlistdf$name_species <-  V(zgg$result_analysis$graph)[j]$name_species
    retlistdf <- rbind(retlistdf,auxlistdf)
  }
  return(retlistdf)
}

# Configuration of inner ziggurats coordinates
conf_ziggurat <- function(igraphnet, basex,widthx,basey,ystep,numboxes,fillcolor, strlabels, strguild, inverse = "no", edge_tr = "no")
{
  kdeg <- rep(0,length(strlabels))
  kdist <- rep(1,length(strlabels))
  knames <- rep("",length(strlabels))
  d2 <- data.frame(strlabels,kdeg,kdist,knames,stringsAsFactors=FALSE)
  names(d2) <- c("label","kdegree","kradius","name_species")
  for (i in 1:nrow(d2)){
    d2[i,]$kdegree <- igraphnet[paste0(strguild,d2[i,]$label)]$kdegree
    d2[i,]$kradius <- igraphnet[paste0(strguild,d2[i,]$label)]$kradius
    name <- igraphnet[paste0(strguild,d2[i,]$label)]$name_species
    d2[i,]$name_species <- igraphnet[paste0(strguild,d2[i,]$label)]$name_species
  }
  d2 <- d2[order(d2$kradius),]
  yjump <- 0.2*zgg$height_y
  x1 <- c()
  x2 <- c()
  y1 <- c()
  y2 <- c()
  r <- c()
  col_row <- c()
  fmult_hght <- 1
  if (edge_tr == "no"){
    xstep <- widthx/numboxes
    if (numboxes > 3)
    {
      topx <- basex + widthx
    }
    else{
      basex <- basex + 0.2*widthx
      topx <- basex + 0.8*widthx
    }
  }
  else{
    xstep <- 1.2*widthx/numboxes
    jump <- 0.25*min((1+0.05*numboxes),2)
    basex <- basex + jump*widthx
    topx <- basex + 0.9*widthx
  }
  for (j in (1:numboxes))
  {
    x1 <- c(x1, basex+(j-1)* (xstep/2))
    if (edge_tr == "yes")
      x2 <- c(x2, topx-(j-1)*(xstep/8)  )
    else
      x2 <- c(x2, topx-(j-1)*xstep/4)
    y1 <- c(y1, basey-(j-1)*(ystep*fmult_hght+yjump)  )
    y2 <- c(y2, y1[j]+ystep*fmult_hght)
    r <- c(r,j)
    col_row <- c(col_row,fillcolor[1+j%%2])
  }
  label <- d2$label
  kdegree <- d2$kdegree
  kradius <- d2$kradius
  name_species <- as.character(d2$name_species)
  d1 <- data.frame(x1, x2, y1, y2, r, col_row, label, kdegree, kradius, name_species, stringsAsFactors=FALSE)
  if (inverse == "yes")
  {
    d1$y1 <- -(d1$y1)
    d1$y2 <- -(d1$y2)
  }
  return(d1)
}

# Draw one of the inner ziggurats
draw_individual_ziggurat <- function(idPrefix, igraphnet, kc, basex = 0, widthx = 0, basey = 0, ystep = 0, numboxes = 0, strlabels = "", strguild = "",
                          sizelabels = 3, colorb = "", grafo = "", svg, zinverse ="no", edge = "no", angle = 0)
{
  dr <- conf_ziggurat(igraphnet, basex,widthx,basey,ystep,numboxes,colorb, strlabels,
                      strguild, inverse = zinverse, edge_tr = edge)
  nsp <- name_species_preprocess(kc,dr,zgg$kcore_species_name_display,
                                 zgg$kcore_species_name_break)
  auxhjust <- rep(0,nrow(dr))
  for (i in 1:length(auxhjust))
    auxhjust[i] <- ifelse(zinverse=="yes",1-i%%2,i%%2)
  p <- grafo + geom_rect(data=dr, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2),
                         fill = dr$col_row, alpha = zgg$alpha_level,color="transparent")
  svg$rect(idPrefix=idPrefix, data=dr, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill=dr$col_row, alpha=zgg$alpha_level,
           color="transparent", size=0.5)
  labelcolor <- ifelse(length(zgg$labels_color)>0,zgg$labels_color[2-as.numeric(strguild == zgg$str_guild_a)],dr$col_row)
  if (is.element(kc,zgg$kcore_species_name_display)){
    pangle <- ifelse(zgg$flip_results,90,0)
    p <- p + geom_text(data=dr, aes(x=x1, y= (y2+y1)/2), color=labelcolor,
                       label = nsp$labelszig, size=zgg$lsize_zig, vjust=0.3, angle = pangle,
                       hjust = 0)
    svg$text(idPrefix=idPrefix, data=dr, mapping=aes(x=x1, y= (y2+y1)/2), color=labelcolor, label = nsp$labelszig, size=zgg$lsize_zig, angle=pangle)
  } else {
    p <- p + geom_text(data=dr, aes(x=x1+0.3*(x2-x1)/2+0.65*((max(r)-r)%%2)*(x2-x1),
                                y= (y2+y1)/2), color=labelcolor,
                   label = nsp$labelszig, size=zgg$lsize_zig, vjust=0.3)
    svg$text(idPrefix=idPrefix, data=dr, mapping=aes(x=x1+0.3*(x2-x1)/2+0.65*((max(r)-r)%%2)*(x2-x1), y=(y2+y1)/2), color=labelcolor, label=nsp$labelszig, size=zgg$lsize_zig)
  }

  calc_grafs <- list("p" = p, "svg" = svg, "dr" = dr)
  return(calc_grafs)
}

# Find orphans. Orphans are species 1-shell connected to species 1-shell
# They form werid chains AKA chains of speciaists
find_orphans <- function(mtxlinks,orphans,gnet,guild_a="yes")
{
  m <- 0
  orph <- NA
  partner <- NA
  kcore <- NA
  repeated <- NA
  df_orph <- data.frame(orph,partner,kcore,repeated)
  if (!is.null(orphans))
    for (i in orphans)
    {
      if (guild_a == "yes")
      {
        partner <- mtxlinks$guild_b[(mtxlinks$guild_a == paste0(zgg$str_guild_a,i))]
        str_opp <- zgg$str_guild_b
        str_own <- zgg$str_guild_a
      }
      else{
        partner <- mtxlinks$guild_a[(mtxlinks$guild_b == paste0(zgg$str_guild_b,i))]
        str_opp <- zgg$str_guild_a
        str_own <- zgg$str_guild_b
      }
      for (t in 1:length(partner))
      {
        m <- m+1
        df_orph[m,]$orph <- i
        df_orph[m,]$partner <- strsplit(as.character(partner[t]),str_opp)[[1]][2]
        df_orph[m,]$kcore <- zgg$g[as.character(partner[t])]$kcorenum
        if (length(partner)>1)
          df_orph[m,]$repeated = "yes"
        else
          df_orph[m,]$repeated <- "no"
      }
      df_orph <- df_orph[!is.na(df_orph$orph),]
    }
  return(df_orph)
}

# Add one link to the graph
add_link <- function(xx1 = 0,xx2 = 0,yy1 = 0,yy2 = 0,
                      slink = 1,clink = c("gray70"),alpha_l = 0.1, spline = "no")
{
  if (!zgg$use_spline)
    spline  <- "no"
  link <- data.frame(x1=xx1, x2 = xx2, y1 = yy1,  y2 = yy2, weightlink = slink)
  npoints_link <- zgg$spline_points
  col_link <- clink[1]
  if (spline == "no")
    zgg$straight_links <- rbind(zgg$straight_links,link)
  else{
    if (spline == "horizontal"){
      x <- c(link$x1,link$x1+(link$x2-link$x1)*0.05,link$x1+(link$x2-link$x1)*0.75,link$x2)
      y <- c(link$y1,link$y1+(link$y2-link$y1)*0.1,link$y1+(link$y2-link$y1)*0.65,link$y2)
    }
    else if (spline == "diagonal"){
      x <- c(link$x1,link$x1+(link$x2-link$x1)*0.5,link$x2)
      y <- c(link$y1,link$y1+(link$y2-link$y1)*0.55,link$y2)
    }
    else if (spline == "vertical"){
      x <- c(link$x1,link$x1+(link$x2-link$x1)*0.85,link$x2)
      y <- c(link$y1,link$y1+(link$y2-link$y1)*0.9,link$y2)
    }
    else if (spline == "lshaped"){
      x <- c(link$x1,link$x1+(link$x2-link$x1)*0.80,link$x1+(link$x2-link$x1)*0.90,link$x2)
      y <- c(link$y1,link$y1+(link$y2-link$y1)*0.95,link$y1+(link$y2-link$y1)*0.99,link$y2)
    }
    else if (spline == "specialisthorizontal"){
      x <- c(link$x1,link$x1+(link$x2-link$x1)*0.40,link$x1+(link$x2-link$x1)*0.75,link$x2)
      y <- c(link$y1,link$y1+(link$y2-link$y1)*0.6,link$y1+(link$y2-link$y1)*0.70,link$y2)
    }
    else if (spline == "arc"){
      x <- c(link$x1,link$x1+(link$x2-link$x1)*0.9,link$x2)
      y <- c(link$y1,link$y1+(link$y2-link$y1)*0.85,link$y2)
    }
    xout <- seq(min(x),max(x),length.out = npoints_link)
    s1 <- spline(x,y,xout=xout,method='natural')
    ds1 <- as.data.frame(s1)
    zgg$count_bent_links <- zgg$count_bent_links + 1
    ds1$number <- zgg$count_bent_links
    ds1$weightlink <- slink
    zgg$bent_links <- rbind(zgg$bent_links,ds1)

  }
  return(0)
}

# Analysis of the chains of specialists
specialist_analysis <- function(specialists,opposite_specialists,species)
{
  ldf <- specialists[specialists$orph == species,]
  if (max(ldf$kcore)>1)
    return(ldf)
}

# Store speciailsts in initermediate df_store data frame and compute positions
# Very hard, indeed
store_specialist_species <- function (row_orph, df_store, strguild, lado, gap, original_specialists_a, original_specialists_b)
{

  sidex <- lado
  index <- nrow(df_store)+1
  df_store[index,]$kcorepartner <- row_orph$kcore
  separation <- (1+zgg$specialists_boxes_separation_count)*sidex
  tot_specialists <- nrow(original_specialists_a)+nrow(original_specialists_b)
  jumpfactor <- (4-min(3,(tot_specialists%/%10)))
  cgap <- (lado+gap/(5-min(3,(tot_specialists%/%10))))

  if (row_orph$kcore > 1){
    df_store$guild[index] <- as.character(strguild)
    df_store$orph[index] <- row_orph$orph
    df_store$partner[index] <- row_orph$partner
    if (strguild == zgg$str_guild_b)
      data_row <- zgg$list_dfs_a[[row_orph$kcore]][zgg$list_dfs_a[[row_orph$kcore]]$label==row_orph$partner,]
    if (strguild == zgg$str_guild_a)
      data_row <- zgg$list_dfs_b[[row_orph$kcore]][zgg$list_dfs_b[[row_orph$kcore]]$label==row_orph$partner,]
    if (row_orph$kcore == zgg$kcoremax)
    {
      if (strguild == zgg$str_guild_b){
        edge_row <- zgg$list_dfs_a[[zgg$kcoremax]][1,]
        xbase <-  min(zgg$last_xtail_a[[zgg$kcoremax]],edge_row$x1 - 2*gap) - gap
      }
      else{
        edge_row <- zgg$list_dfs_b[[zgg$kcoremax]][1,]
        xbase <-  min(zgg$last_xtail_b[[zgg$kcoremax]],edge_row$x1 - 2*gap)- gap
      }
      df_store$x1[index] <- xbase - 2 * gap

      if (zgg$kcoremax > 2)
        df_store$y1[index] <- max(abs(edge_row$y2),abs(edge_row$y1)) + 6*cgap/(zgg$aspect_ratio)
      else{
        xbase <- 0
        df_store$y1[index] <- max(abs(edge_row$y2),abs(edge_row$y1)) + 6*cgap/(zgg$aspect_ratio)
      }

      if (df_store$guild[index] == zgg$str_guild_a){
        df_store$y1[index] = -abs(df_store$y1[index])
        df_store$y2[index] = -abs(df_store$y2[index])
      }
      repetitions <- sum((df_store$partner == row_orph$partner) & (df_store$guild == strguild))
      if (repetitions > 1){
        df_store$x1[index] <- df_store$x1[index] -  (repetitions-1) * sidex
        df_store$y1[index] <- df_store$y1[index] +  sign(df_store$y1[index])*(repetitions-1) * (3+as.integer(zgg$kcoremax>3)) * sidex/zgg$aspect_ratio
      }
      repetitions_root <- sum((df_store$kcorepartner == zgg$kcoremax) & (df_store$guild == strguild))
      if (repetitions_root > 1){
        df_store$x1[index] <- df_store$x1[index] +  2 * (repetitions_root) * (zgg$kcoremax) * sidex
        df_store$y1[index] <- df_store$y1[index] +  sign(df_store$y1[index])*(1/jumpfactor)*((repetitions_root+ 0.7*index) * 3* sidex/zgg$aspect_ratio)
      }

    } else {
      if (row_orph$kcore == 2){
        xoffset <- 2*zgg$specialistskcore2_horizontal_dist_rootleaf_expand*separation   # Controls the separation of specialists root leaves connected to core 2
        zgg$yoffset <- zgg$specialistskcore2_vertical_dist_rootleaf_expand*separation/zgg$aspect_ratio
      } else{
        xoffset <- 0
        zgg$yoffset <- 0
      }
      if (strguild == zgg$str_guild_b)
      {
        data_row_pic <- zgg$list_dfs_a[[row_orph$kcore]]
        df_store$x1[index]<-  min(data_row_pic$x2) + zgg$factor_hop_x* separation + xoffset
        df_store$y1[index] <- zgg$last_ytail_a[row_orph$kcore] - (1+sqrt(zgg$kcoremax))*sidex/zgg$aspect_ratio - zgg$yoffset
        zgg$last_ytail_a[row_orph$kcore] <- -abs(df_store$y1[index])
      }
      if (strguild == zgg$str_guild_a)
      {
        data_row_pic <- zgg$list_dfs_b[[row_orph$kcore]]
        df_store$x1[index] <- min(data_row_pic$x2) + zgg$factor_hop_x* separation + xoffset
        df_store$y1[index] <- zgg$last_ytail_b[row_orph$kcore] +  (1+sqrt(zgg$kcoremax))*sidex/zgg$aspect_ratio + zgg$yoffset
        zgg$last_ytail_b[row_orph$kcore] <- abs(df_store$y1[index])
      }
    }
    if (row_orph$kcore == zgg$kcoremax)
      df_store$x1[index] <- df_store$x1[index]* zgg$root_specialist_expand[1]
    else if (zgg$root_specialist_expand[1]>1)
      df_store$x1[index] <- df_store$x1[index]* zgg$root_specialist_expand[1]
    else if (row_orph$kcore > 2)
      df_store$x1[index] <- df_store$x1[index]* (9+zgg$root_specialist_expand[1])/10
    df_store$y1[index] <- df_store$y1[index]* zgg$root_specialist_expand[2] * zgg$aspect_ratio
  }
  else {                                          # Branch leaves
    df_store$partner[index] <- row_orph$orph
    df_store$orph[index] <- row_orph$partner
    df_store$guild[index] <- as.character(strguild)
    data_row <- df_store[(df_store$orph == row_orph$orph) & (swap_strguild(strguild) == df_store$guild),]
    repetitions <- sum((df_store$partner == row_orph$orph) & (df_store$guild == strguild))

    if (data_row$kcorepartner == zgg$kcoremax){
      if (zgg$kcoremax > 2){
        df_store$x1[index] <- data_row$x1 - 1.5*separation*(repetitions)
        df_store$y1[index] <- data_row$y1 + sign(data_row$y1)*4*(zgg$specialists_boxes_separation_count*zgg$height_y + (repetitions-1)*sidex/zgg$aspect_ratio)
      }
      else{
        df_store$x1[index] <- data_row$x1 - (4+0.2*sqrt(index))*sidex
        df_store$y1[index] <- data_row$y1+ sign(data_row$y1)*(zgg$specialists_boxes_separation_count*zgg$height_y + as.integer(data_row$partner)*sidex/zgg$aspect_ratio) +
          (2+0.2*sqrt(index))*sign(data_row$y1)*(repetitions-1)*sidex/zgg$aspect_ratio
      }
    }
    else{                                   # specialists connected to root leaf connected to kcoremax
      hjump <- sign(data_row$x1)*zgg$factor_hop_x* separation
      df_store$x1[index] <- data_row$x1 + hjump
      df_store$y1[index] <- data_row$y1 + sign(data_row$y1)*2*(repetitions-1)*sidex/zgg$aspect_ratio
    }
    reps <- sum((df_store$partner == df_store$partner[index]) & (df_store$guild == strguild))

    if (zgg$kcore1specialists_leafs_vertical_separation!=1){
      addjump <- sign(df_store$y1[index])*reps*zgg$kcore1specialists_leafs_vertical_separation*sidex/zgg$aspect_ratio
      df_store$y1[index]<- addjump + df_store$y1[index]
    }

  }
  df_store$x2[index] <- df_store$x1[index] + sidex
  df_store$y2[index] <- df_store$y1[index] + sidex/zgg$aspect_ratio
  if (row_orph$kcore == zgg$kcoremax){
    df_store$xx2[index] <- (data_row$x2+data_row$x1)/2
    df_store$yy2[index] <- data_row$y2
  }
  else
  {
    df_store$yy2[index] <- (data_row$y2+data_row$y1)/2
    if (df_store$kcorepartner[index] > 1){
      df_store$xx2[index] <- data_row$x2
    }
    else {
      df_store$xx2[index] <- data_row$x1 + (data_row$x2>0)*sidex
    }
  }
  return(df_store)
}

# Draw a chain of specialists
draw_specialist_chains <- function(grafo, svg, df_chains, ladosq)
{
  p <- grafo
  df_chains$weightlink <- 0
  for (i in 1:nrow(df_chains))
  {
    sqinverse = "no"
    if (df_chains[i,]$guild == zgg$str_guild_a){
      bgcolor <- zgg$color_guild_a[2]
      is_guild_a <- TRUE
    }
    if (df_chains[i,]$guild == zgg$str_guild_b){
      bgcolor <- zgg$color_guild_b[2]
      is_guild_a <- FALSE
    }
    hjust <- 0
    vjust <- 0
    labelcolor <- ifelse(length(zgg$labels_color)>0,zgg$labels_color[2-as.numeric(is_guild_a)], bgcolor)
    sqlabel = gen_sq_label(df_chains[i,]$orph,is_guild_a = is_guild_a)
    f <- draw_square(paste0(ifelse(is_guild_a, "specialist-chains-kcore1-a-", "specialist-chains-kcore1-b-"), i),p,
                     svg,df_chains[i,]$x1,df_chains[i,]$y1, ladosq,
                     bgcolor,zgg$alpha_level,
                     labelcolor,0,hjust,vjust,
                     sqlabel,
                     lbsize = zgg$lsize_kcore1,inverse = "no")
    p <- f["p"][[1]]
    svg <- f["svg"][[1]]
    if (zgg$paintlinks){

      yy2 = df_chains[i,]$yy2
      xx2 = df_chains[i,]$xx2
      if (df_chains[i,]$x2>0){
        xx1 = df_chains[i,]$x1
      } else {
        xx1 = df_chains[i,]$x2
      }
      if (df_chains[i,]$y2>0){
        yy1 <-  (df_chains[i,]$y1+df_chains[i,]$y2)/2
      } else {
        if (df_chains[i,]$kcorepartner != zgg$kcoremax)
          yy1 <-  df_chains[i,]$y2 - ladosq/(2*zgg$aspect_ratio)
        else
          yy1 <-  df_chains[i,]$y2
      }


      if ( (df_chains[i,]$kcorepartner >1) & (df_chains[i,]$kcorepartner < zgg$kcoremax) )
      {
        splineshape = "lshaped"
      }
      else
        splineshape = "arc"
        #splineshape = "specialisthorizontal"
      #color_link = "green"
      tailweight <- 0
      # for (h in 1:nrow(df_chains))
        if (df_chains$guild[i] == zgg$str_guild_a)
          tailweight <- tailweight + zgg$result_analysis$matrix[as.numeric(df_chains$partner[i]),
                                                                      as.numeric(df_chains$orph[i])]
        else
          tailweight <- tailweight + zgg$result_analysis$matrix[as.numeric(df_chains$orph[i]),
                                                                as.numeric(df_chains$partner[i])]
      df_chains[i,]$weightlink <- get_link_weights(tailweight)
      add_link(xx1=xx1, xx2 = xx2,
                     yy1 = yy1, yy2 = yy2,
                     slink = zgg$size_link*df_chains[i,]$weightlink, clink = c(zgg$color_link),
                     alpha_l = zgg$alpha_link , spline = splineshape)
    }
  }
  calc_vals <- list("p" = p, "svg" = svg)
  return(calc_vals)
}

swap_strguild <- function(strguild)
{
  if (strguild == zgg$str_guild_a)
    strguild = zgg$str_guild_b
  else
    strguild = zgg$str_guild_a
  return(strguild)
}

# Store the root node of a chain of specialists
store_root_leaf <- function(specialists,df_chains,strguild,lado, gap, original_specialists_a, original_specialists_b)
{
  for (i in 1:nrow(specialists))
  {
    if (specialists[i,]$kcore > 1){
      df_chains <- store_specialist_species(specialists[i,], df_chains, strguild, lado , gap, original_specialists_a, original_specialists_b)
      specialists[i,]$drawn = "yes"
    }
  }
  calc_vals <- list("df_chains" = df_chains, "specialists" = specialists)
  return(calc_vals)
}

# Store leafs of a specialist chain
store_branch_leaf <- function(specialists, specialists_opp,df_chains, pstrguild, plado, gap, original_specialists_a, original_specialists_b)
{
  for (i in 1:nrow(specialists))
  {
    if (specialists[i,]$drawn == "no"){
      strguild <- pstrguild
      if (sum( ((df_chains$orph == specialists[i,]$orph) & (df_chains$guild == strguild)) )>0 ){
        strguild <- swap_strguild(pstrguild)
        df_chains <- store_specialist_species(specialists[i,], df_chains, strguild, plado, gap, original_specialists_a, original_specialists_b)
        specialists[i,]$drawn = "yes"
        mirror_specialist <- which((specialists_opp$partner == specialists[i,]$orph) & (specialists_opp$orph == specialists[i,]$partner))
        if (length(mirror_specialist)>0)
          specialists_opp[mirror_specialist , ]$drawn = "yes"
      }
    }
  }
  calc_vals <- list("df_chains" = df_chains, "specialists" = specialists, "specialists_opp" = specialists_opp)
  return(calc_vals)
}

# Draw links among species of same k-index
draw_innercores_tails <- function(p,svg,kc,list_dfs,df_orph,color_guild, inverse="no", is_guild_a = TRUE)
{
  lastx <- 0
  lasty <- 0

  lpoint_x <- 0
  if (length(list_dfs[[kc]])>0)
    if (kc>2)
      lpoint_x <- list_dfs[[kc]][nrow(list_dfs[[kc]]),]$x2
  else
    lpoint_x <- list_dfs[[kc]][nrow(list_dfs[[kc]]),]$x2 + 4*zgg$lado
  if (kc>2)
    lpoint_y <- (list_dfs[[kc]][1,]$y1+list_dfs[[kc]][1,]$y2)/2
  else{
    if (zgg$kcoremax > 2)
      lpoint_y <- (list_dfs[[kc]][1,]$y1+list_dfs[[kc]][1,]$y2)/4
    else{
      zgg$kcore2tail_vertical_separation <- zgg$kcore2tail_vertical_separation + 2
      lpoint_y <- 1.5*max(list_dfs[[2]]$y2)
    }
  }
  long_tail <- df_orph[(df_orph$kcore == kc) & (df_orph$repeated == "no"),]
  if (length(long_tail)>0){
    v<-  draw_edge_tails(p,svg,lpoint_x,lpoint_y,kc,long_tail,list_dfs,color_guild,
                         inverse = inverse, joinchars = zgg$joinstr,pbackground = "no",
                         tspline = "lshaped", is_guild_a = is_guild_a)
    p <- v["p"][[1]]
    svg <- v["svg"][[1]]
    if (length(v["lastx"][[1]])>0)
      lastx <- v["lastx"][[1]]
    if (length(v["lasty"][[1]])>0)
      lasty <-v["lasty"][[1]]
  }
  calc_vals <- list("p" = p, "svg" = svg, "point_x" = lpoint_x, "point_y" = lpoint_y, "lastx" = lastx, "lasty" = lasty)
  return(calc_vals)
}

# Draw lnks amnog inner ziggurats
draw_inner_links <- function(p, svg)
{
    for (kcb in seq(zgg$kcoremax,2))
    {
      for (kc in seq(zgg$kcoremax,2))
      {
        labels_a <- zgg$list_dfs_a[[kc]]$label
        for (j in seq(along=labels_a))
        {
          numberlinksa <- sum(zgg$mtxlinks$guild_a == paste0(zgg$str_guild_a,labels_a[j]) )
          data_a <- zgg$list_dfs_a[[kc]][j,]
          foundlinksa <- 0
          labels_b <- zgg$list_dfs_b[[kcb]]$label
          for (i in seq(along =labels_b))
          {
            if (sum(zgg$mtxlinks$guild_a == paste0(zgg$str_guild_a,labels_a[j]) & zgg$mtxlinks$guild_b == paste0(zgg$str_guild_b,labels_b[i]))>0)
            {
              foundlinksa <- foundlinksa + 1
              data_b <- zgg$list_dfs_b[[kcb]][i,]
              weightlink <- get_link_weights(zgg$result_analysis$matrix[as.numeric(data_b$label),
                                                                        as.numeric(data_a$label)])
              bend_line = "no"
              if (((kc == 2) & (kcb == zgg$kcoremax)) | ((kc == zgg$kcoremax) & (kcb == 2)))
                bend_line = "horizontal"
              if ((kc == zgg$kcoremax) & (kcb == zgg$kcoremax))
              {
                link <- data.frame(x1= data_a$x1 + (data_a$x2-data_a$x1)/2,
                                   x2 = data_b$x1 +(data_b$x2-data_b$x1)/2,
                                   y1 = data_a$y1,  y2 = zgg$list_dfs_b[[kcb]]$y1[i] )
                lcolor = "orange"
                bend_line = "no"
              }
              else if (kc == kcb) {
                link <- data.frame(x1=  data_a$x1,
                                   x2 = data_b$x1,
                                   y1 = data_a$y1,
                                   y2 = data_b$y1 )
                bend_line = "no"
                lcolor = "pink"
              }
              else if (kc > kcb) {
                if (kc == zgg$kcoremax)
                  link <- data.frame(x1= (data_a$x2 + data_a$x1)/2,
                                     x2 = data_b$x1,
                                     y1 = data_a$y2,  y2 = data_b$y1 )
                else{
                  link <- data.frame(x1=  data_a$x2 ,
                                     x2 = data_b$x1,
                                     y1 = data_a$y1,  y2 = data_b$y1 )
                  bend_line = "diagonal"
                }
                lcolor = "green"
              }
              else
              {
                if (kcb == zgg$kcoremax){
                  y_2 <- data_b$y2
                  x_2 <- (data_b$x2 + data_b$x1)/2
                }
                else{
                  y_2 <- data_b$y1
                  x_2 <- data_b$x2

                }
                link <- data.frame(x1= data_a$x1,
                                   x2 = x_2,
                                   y1 = data_a$y1,  y2 = y_2)
                lcolor = "blue"
              }

              add_link(xx1=link$x1, xx2 = link$x2,
                             yy1 = link$y1, yy2 = link$y2,
                             slink = zgg$size_link*weightlink, clink =  c(zgg$color_link),
                             alpha_l = zgg$alpha_link , spline = bend_line)
            }
            if (foundlinksa >= numberlinksa )
              break
          }
        }
      }
    }
    calc_vals <- list("p" = p, "svg" = svg)
    return(calc_vals)
}

# Draw tail connected to highest kdegree node
draw_fat_tail<- function(p,svg,fat_tail,nrows,list_dfs,color_guild,pos_tail_x,pos_tail_y,fattailjhoriz,fattailjvert,fgap,
                         inverse="no", is_guild_a =TRUE)
{
  ppos_tail_x <- pos_tail_x * fattailjhoriz
  pos_tail_y <- (0.25+0.1*sqrt(nrows))*(list_dfs[[zgg$kcoremax]][1,]$y2+list_dfs[[zgg$kcoremax]][1,]$y1)/2
  ppos_tail_y <- pos_tail_y * fattailjvert
  if (nrow(fat_tail)>0)
  {
    plyy2 <- ifelse(inverse == "yes", list_dfs[[zgg$kcoremax]][1,]$y1-3*zgg$lado, list_dfs[[zgg$kcoremax]][1,]$y2-3*zgg$lado)
    v<- draw_tail(ifelse(is_guild_a, "fat-kcore1-a", "fat-kcore1-b"), p,svg,
                  fat_tail,zgg$lado,color_guild,gen_sq_label(fat_tail$orph,is_guild_a = is_guild_a),
                  ppos_tail_x,ppos_tail_y,fgap,
                  lxx2 = list_dfs[[zgg$kcoremax]][1,]$x1,
                  lyy2 = plyy2,
                  sqinverse = inverse, background = "no", psize = zgg$lsize_kcore1,
                  is_guild_a = is_guild_a, wlink = fat_tail$weightlink[1])
    p <- v["p"][[1]]
    svg <- v["svg"][[1]]
  }
  calc_vals <- list("p" = p, "svg" = svg, "pos_tail_x" = zgg$pos_tail_x)
  return(calc_vals)
}

# Management of chains of specialists
handle_specialists <- function(p,svg,specialists_a,specialists_b,lado,gap)
{
  ladosq <- 2 * lado * sqrt(zgg$square_nodes_size_scale)
  specialists_a <- data.frame(c())
  specialists_b <- data.frame(c())
  if (exists("df_orph_a", envir = zgg))
    if (nrow(zgg$df_orph_a)>0)
     {
      specialists_a <-  zgg$df_orph_a[zgg$df_orph_a$repeated== "yes",]
      specialists_a <-  specialists_a[rev(order(specialists_a$orph,specialists_a$kcore)),]
      if (nrow(specialists_a)>0)
        specialists_a$drawn <- "no"
      }
  if (exists("df_orph_b", envir = zgg))
    if (nrow(zgg$df_orph_b)>0)
      {
      specialists_b <-  zgg$df_orph_b[zgg$df_orph_b$repeated== "yes",]
      specialists_b <-  specialists_b[rev(order(specialists_b$orph,specialists_b$kcore)),]
      if (nrow(specialists_b)>0)
        specialists_b$drawn <- "no"
      }

  # Create empty df_chains data frame
  zgg$df_chains <- data.frame(x1 = numeric(0), x2 = numeric(0), y1 = numeric(0), y2 = numeric(0),
                          guild = character(0), orph = integer(0), partner = integer(0),
                          kcorepartner = integer(0), xx2 = numeric(0), yy2 = numeric(0), stringsAsFactors = FALSE )

  if  (( ( nrow(specialists_a)+nrow(specialists_b) )>0)) {
    original_specialists_a <- specialists_a
    original_specialists_b <- specialists_b
    while (((nrow(specialists_a)+nrow(specialists_b))>0))
    {
      if (nrow(specialists_a)>0){
        k <- store_root_leaf(specialists_a, zgg$df_chains, zgg$str_guild_a, ladosq, gap, original_specialists_a, original_specialists_b)
        zgg$df_chains <- k["df_chains"][[1]]
        specialists_a <- k["specialists"][[1]]
      }
      if (nrow(specialists_b)>0){
        k <- store_root_leaf(specialists_b, zgg$df_chains, zgg$str_guild_b, ladosq, gap, original_specialists_a, original_specialists_b)
        zgg$df_chains <- k["df_chains"][[1]]
        specialists_b <- k["specialists"][[1]]
      }
      if (nrow(specialists_a)>0){
        k <- store_branch_leaf(specialists_a, specialists_b, zgg$df_chains, zgg$str_guild_a, ladosq, gap, original_specialists_a, original_specialists_b)
        zgg$df_chains <- k["df_chains"][[1]]
        specialists_a <- k["specialists"][[1]]
        specialists_b <- k["specialists_opp"][[1]]
      }
      if (nrow(specialists_b)>0){
        k <- store_branch_leaf(specialists_b, specialists_a, zgg$df_chains, zgg$str_guild_b, ladosq, gap, original_specialists_a, original_specialists_b)
        zgg$df_chains <- k["df_chains"][[1]]
        specialists_b <- k["specialists"][[1]]
        specialists_a <- k["specialists_opp"][[1]]
      }
      # Now they may be some specialists of core 1 linked to core 1 that were not
      # stored in the previous procedure
      specialists_a <- specialists_a[specialists_a$drawn == "no",]
      specialists_b <- specialists_b[specialists_b$drawn == "no",]
    }
    f <- draw_specialist_chains(p, svg, zgg$df_chains, ladosq)
    p <- f["p"][[1]]
    svg <- f["svg"][[1]]
  }
  calc_vals <- list("p" = p, "svg" = svg, "df_chains" = zgg$df_chains)
  return(calc_vals)
}

# Final annotations
write_annotations <- function(p, svg)
{
  title_text = ""
  if (zgg$show_title)
    title_text <- paste0("Network: ",zgg$network_name)
  # Legend size conversion factor
  lzcf <- 3
  legends_text <- paste0(
    "<span style = 'color:",zgg$color_guild_a[1],"; font-size:",lzcf*zgg$lsize_legend,"pt'>",zgg$name_guild_a,
    "</span> <span style = 'color:",zgg$color_guild_b[1],"; font-size:",lzcf*zgg$lsize_legend,"pt'>",zgg$name_guild_b,"</span>")
  p <- p+ ggtitle(title_text,subtitle=paste("<span align='right' size>",legends_text,"</span>"))
  p <- p + coord_fixed(ratio=zgg$aspect_ratio) +theme_bw() + theme(panel.grid.minor.x = element_blank(),
                                                               panel.grid.minor.y = element_blank(),
                                                               panel.grid.major.x = element_blank(),
                                                               panel.grid.major.y = element_blank(),
                                                               axis.text.x = element_blank(),
                                                               axis.text.y = element_blank(),
                                                               axis.ticks.x=element_blank(),
                                                               axis.ticks.y=element_blank(),
                                                               axis.title.x = element_blank(),
                                                               axis.title.y = element_blank(),                                                               plot.background = element_rect(fill = zgg$backg_color),
                                                               panel.background = element_rect(fill = zgg$backg_color),
                                                               plot.title = element_text(size = lzcf*(zgg$lsize_legend+1),
                                                                                         hjust = 0.5,
                                                                                         face="plain"),
                                                               
                                                               plot.subtitle = ggtext::element_markdown()
                                                                )
  if (zgg$hide_plot_border)
    p <- p + theme(panel.border=element_blank())
  landmark_top <- 1.2*max(zgg$last_ytail_b[!is.na(zgg$last_ytail_b)],1.2*zgg$ymax)*zgg$rescale_plot_area[2]
  mlabel <- "."
  landmark_right <- (zgg$tot_width+2*zgg$hop_x)*zgg$rescale_plot_area[1]
  f <- draw_square("annotation",p,svg,landmark_right,0,1,"transparent",0.5,"transparent",0,0,0,slabel="")
  p <- f["p"][[1]]
  svg <- f["svg"][[1]]
  p <- p +annotate(geom="text", x= landmark_right, y=0, label=mlabel,
                   colour = "red", size=1, hjust = 0, vjust = 0, angle = 0)
  svg$text("annotation", data=data.frame(x=landmark_right, y=0), mapping=aes(x=x, y=y), color="red", label=mlabel, size=1, angle=0)
  landmark_left <- min(zgg$last_xtail_a[[zgg$kcoremax]],zgg$last_xtail_b[[zgg$kcoremax]])-min(zgg$hop_x,0.2*min(zgg$last_xtail_a[[zgg$kcoremax]],zgg$last_xtail_b[[zgg$kcoremax]]))
  landmark_left <- min(landmark_left, zgg$pos_tail_x)*zgg$rescale_plot_area[1]
  p <- p +annotate(geom="text", x=landmark_left, y=0, label=mlabel,
                   colour = "red", size=2, hjust = 0, vjust = 0, angle = 0)
  svg$text("annotation", data=data.frame(x=landmark_left, y=0), mapping=aes(x=x, y=y), color="red", label=mlabel, size=1, angle=0)
  x_span <- landmark_right - landmark_left

  if (!(zgg$flip_results)){
    x_legend <- 0.8*landmark_right#*(1+zgg$displace_legend[1])
    y_legend <- 0.8*landmark_top#*(1+zgg$displace_legend[2])
  } else {
    x_legend <- 0.8*landmark_top#*(1+zgg$displace_legend[2])
    y_legend <- 0.8*landmark_right#*(1+zgg$displace_legend[1])
  }
  # p <- p + annotate(geom="text", x=x_legend,
  #                   y=y_legend,
  #                   label=zgg$name_guild_a,
  #                   colour = zgg$color_guild_a[1], size=zgg$lsize_legend,
  #                   hjust = 1, vjust = 0, angle = 0)
  # p <- p + annotate(geom="text", x=x_legend,
  #                   y=y_legend,
  #                   label= paste(rep(" ",length(zgg$name_guild_a))," ",zgg$name_guild_b),
  #                   colour = zgg$color_guild_b[1], size=zgg$lsize_legend,
  #                   hjust = 0, vjust = 0, angle = 0)
  # p <- p +annotate(geom="text", x=landmark_left,
  #                  y = y_legend,
  #                  label="1-shell",
  #                  colour = zgg$corecols[2], size=zgg$lsize_core_box, hjust = 0, vjust = 0,
  #                  angle = 0,
  #                  fontface="italic")
  # svg$text("core-1", data=data.frame(x=landmark_left, y=y_legend), mapping=aes(x=x, y=y), color=zgg$corecols[2],
  #          label="1-shell",
  #          size=zgg$lsize_core_box, angle=0)
  # svg$text("core-1", data=data.frame(x=(0.8+(0.3*zgg$displace_legend[1]))*x_legend, y=y_legend), mapping=aes(x=x, y=y), size=zgg$lsize_legend,
  #          label=zgg$name_guild_a, color=zgg$color_guild_a[1], angle=0)
  # svg$text("core-1", data=data.frame(x=x_legend, y=y_legend), mapping=aes(x=x, y=y), size=zgg$lsize_legend,
  #          label=zgg$name_guild_b, color=zgg$color_guild_b[1], angle=0)

  calc_vals <- list("p" = p, "svg" = svg)
  return(calc_vals)
}

# Handle specialist chain species
handle_orphans <- function(vg)
{
  zgg$df_orph_a <- data.frame(c())
  zgg$df_orph_b <- data.frame(c())
  zgg$mtxlinks <- data.frame(igraph::get.edgelist(vg))
  names(zgg$mtxlinks) <- c("guild_a","guild_b")
  if (length(grep(zgg$str_guild_b,zgg$mtxlinks[1,1]))>0)
    names(zgg$mtxlinks) <- rev(names(zgg$mtxlinks))
  zgg$orphans_a <- zgg$df_cores$species_guild_a[[1]]
  zgg$orphans_b <- zgg$df_cores$species_guild_b[[1]]
  if (!is.null(zgg$orphans_a))
    if (!is.na(zgg$orphans_a[1]))
      zgg$df_orph_a <- find_orphans(zgg$mtxlinks,zgg$orphans_a,zgg$g,guild_a="yes")
  if (!is.null(zgg$orphans_b))
    if (!is.na(zgg$orphans_b[1]))
      zgg$df_orph_b <- find_orphans(zgg$mtxlinks,zgg$orphans_b,zgg$g,guild_a="no")
  calc_vals <- list("mtxlinks" = zgg$mtxlinks, "orphans_a" = zgg$orphans_a,
                    "orphans_b" = zgg$orphans_b, "df_orph_a" = zgg$df_orph_a, "df_orph_b" = zgg$df_orph_b )
  return(calc_vals)
}

# Draw specialist chains connected to inner ziggurats
draw_inner_orphans <- function(p, svg)
{
  if (zgg$kcoremax >2)
  {
    for (kc in seq(from = zgg$kcoremax-1, to = 2))
    {

      if (exists("df_orph_a", envir = zgg)){
        w <- draw_innercores_tails(p,svg,kc,zgg$list_dfs_b,zgg$df_orph_a,zgg$color_guild_a, inverse="no")
        p <- w["p"][[1]]
        svg <- w["svg"][[1]]
        zgg$last_xtail_b[kc] <- w["lastx"][[1]]
        zgg$last_ytail_b[kc] <-w["lasty"][[1]]
      }

      if (exists("df_orph_b", envir = zgg)){
        w <- draw_innercores_tails(p,svg,kc,zgg$list_dfs_a,zgg$df_orph_b,zgg$color_guild_b, inverse="yes", is_guild_a = FALSE)
        p <- w["p"][[1]]
        svg <- w["svg"][[1]]
        zgg$last_xtail_a[kc] <- w["lastx"][[1]]
        zgg$last_ytail_a[kc] <-w["lasty"][[1]]
      }
    }
  } else {                                  # orphans when kcoremax == 2
    zgg$kcore2tail_vertical_separation <- zgg$kcore2tail_vertical_separation + 2
    if (exists("df_orph_a", envir = zgg)){
      w <- draw_innercores_tails(p,svg,2,zgg$list_dfs_b,zgg$df_orph_a,zgg$color_guild_a, inverse="yes")
      p <- w["p"][[1]]
      svg <- w["svg"][[1]]
      zgg$last_xtail_b[2] <- w["lastx"][[1]]
      zgg$last_ytail_b[2] <-w["lasty"][[1]]
    }

    if (exists("df_orph_b", envir = zgg)){
      w <- draw_innercores_tails(p,svg,2,zgg$list_dfs_a,zgg$df_orph_b,zgg$color_guild_b, inverse="no", is_guild_a = FALSE)
      p <- w["p"][[1]]
      svg <- w["svg"][[1]]
      zgg$last_xtail_a[2] <- w["lastx"][[1]]
      zgg$last_ytail_a[2] <- w["lasty"][[1]]
    }
  }

  calc_vals <- list("p" = p, "svg" = svg)
  return(calc_vals)
}

# Manage fat tails
handle_fat_tails <- function(p, svg)
{
  fat_tail_x <- min(zgg$last_xtail_a[[zgg$kcoremax]],zgg$last_xtail_b[[zgg$kcoremax]],zgg$list_dfs_a[[zgg$kcoremax]][1,]$x1,zgg$list_dfs_b[[zgg$kcoremax]][1,]$y2)
  max_b_kdegree <- zgg$list_dfs_b[[zgg$kcoremax]][which(zgg$list_dfs_b[[zgg$kcoremax]]$kdegree == max(zgg$list_dfs_b[[zgg$kcoremax]]$kdegree)),]$label
  if (exists("df_orph_a", envir = zgg)){
    fat_tail_a <- zgg$df_orph_a[(zgg$df_orph_a$partner == max(max_b_kdegree)) & (zgg$df_orph_a$repeated == "no"),]
    tailweight <- 0
    if (nrow(fat_tail_a)>0) {
    for (h in 1:nrow(fat_tail_a))
      tailweight <- tailweight + zgg$result_analysis$matrix[as.numeric(fat_tail_a$partner[h]),
                                                                  as.numeric(fat_tail_a$orph[h])]
      fat_tail_a$weightlink <- get_link_weights(tailweight)
    }
  }
  if (!exists("fat_tail_a"))
    fat_tail_a <- data.frame(c())
  max_a_kdegree <- zgg$list_dfs_a[[zgg$kcoremax]][which(zgg$list_dfs_a[[zgg$kcoremax]]$kdegree == max(zgg$list_dfs_a[[zgg$kcoremax]]$kdegree)),]$label
  if (exists("df_orph_b", envir = zgg)){
    fat_tail_b <- zgg$df_orph_b[(zgg$df_orph_b$partner == max(max_a_kdegree)) & (zgg$df_orph_b$repeated == "no"),]
    tailweight <- 0
    if (nrow(fat_tail_b)>0) {
      for (h in 1:nrow(fat_tail_b))
        tailweight <- tailweight+zgg$result_analysis$matrix[as.numeric(fat_tail_b$orph[h]),
                                                    as.numeric(fat_tail_b$partner[h])]
      fat_tail_b$weightlink <- get_link_weights(tailweight)
    }
  }
  if (!exists("fat_tail_b"))
    fat_tail_b <- data.frame(c())
  fgap <- 0.7*zgg$hop_x
  zgg$pos_tail_x <- min(zgg$last_xtail_a[[zgg$kcoremax]],zgg$last_xtail_b[[zgg$kcoremax]],zgg$list_dfs_b[[zgg$kcoremax]][1,]$x1-fgap,zgg$list_dfs_a[[zgg$kcoremax]][1,]$x1-fgap)
  nrows_fat <- nrow(fat_tail_b)+nrow(fat_tail_a)
  if ((exists("fat_tail_a")) & (zgg$kcoremax > 2)) {
    f <- draw_fat_tail(p,svg,fat_tail_a,nrows_fat,zgg$list_dfs_b,zgg$color_guild_a[2],
                       zgg$pos_tail_x,pos_tail_y,zgg$fattailjumphoriz[1],
                       zgg$fattailjumpvert[1],fgap,inverse="yes")
    p <- f["p"][[1]]
    svg <- f["svg"][[1]]
  }

  if ((exists("fat_tail_b")) & (zgg$kcoremax > 2)) {
    f <- draw_fat_tail(p,svg,fat_tail_b,nrows_fat,zgg$list_dfs_a,zgg$color_guild_b[2],zgg$pos_tail_x,
                       pos_tail_y,zgg$fattailjumphoriz[2],zgg$fattailjumpvert[2],fgap,inverse="no", is_guild_a = FALSE)
    p <- f["p"][[1]]
    svg <- f["svg"][[1]]
  }
  calc_vals <- list("p" = p, "svg" = svg, "pos_tail_x" = zgg$pos_tail_x)
  return(calc_vals)
}

draw_all_ziggurats <- function(p, svg)
{
  if (zgg$kcoremax > 2)
  {
    for (kc in seq(from = zgg$kcoremax-1, to = 2))
    {
      if (sum(zgg$df_cores$num_species_guild_a[kc],zgg$df_cores$num_species_guild_b[kc])>0)
      {
        zgg$pointer_x <- (zgg$kcoremax-kc)*zgg$hop_x
        zgg$pointer_y <- zgg$pointer_y - sign(zgg$pointer_y)*(1/(0.8*(zgg$kcoremax-2))*(zgg$kcoremax-kc))*zgg$strips_height
        if (zgg$df_cores$num_species_guild_a[kc] > 0){
          despl_pointer_y <- zgg$displace_y_a[kc] * zgg$ymax
          if ((kc == 2) )
          {
            edge_core <- "yes"
            zgg$pointer_y <- max(4,zgg$species_in_core2_a)*zgg$height_y+abs(max(zgg$list_dfs_b[[zgg$kcoremax]]$y2))
          }
          else {
            edge_core <- "no"
            if (zgg$primerkcore)
              zgg$pointer_y <- zgg$ymax
            else if ((zgg$df_cores$num_species_guild_a[kc-1] > 5) & !(zgg$primerkcore)){
              zgg$pointer_y <- zgg$pointer_y - (0.8+0.1*(zgg$kcoremax-kc)/zgg$kcoremax)*zgg$strips_height
            }
          }
          pystep <- zgg$height_y
          if (zgg$kcoremax < 5)
            wzig <- (0.8+((kc)/zgg$kcoremax))* zgg$width_zig
          else
            wzig <- 1.3*zgg$width_zig
          if (kc == 2)
            wzig <- wzig *min(2,(1+0.1*sqrt(max(zgg$df_cores$num_species_guild_a[kc],zgg$df_cores$num_species_guild_b[kc]))))
          zig <-  draw_individual_ziggurat(paste0("kcore", kc, "-a"), zgg$g, kc, basex = zgg$pointer_x, widthx = wzig,
                                basey = zgg$pointer_y + despl_pointer_y, ystep = pystep, strlabels = zgg$df_cores$species_guild_a[[kc]],
                                strguild = zgg$str_guild_a,
                                sizelabels = zgg$lsize_zig, colorb = zgg$color_guild_a, numboxes = zgg$df_cores[kc,]$num_species_guild_a,
                                zinverse = "yes", edge = edge_core, grafo = p, svg = svg)
          p <- zig["p"][[1]]
          svg <- zig["svg"][[1]]
          zgg$list_dfs_a[[kc]] <- zig["dr"][[1]]
          zgg$last_xtail_a[[kc]] <- max(zgg$list_dfs_a[[kc]]$x2)
          zgg$last_ytail_a[[kc]] <- -max(abs(zgg$list_dfs_a[[kc]]$y2))
          zgg$posic_zig_a <- -max(zgg$posic_zig_a, abs(min(zgg$list_dfs_a[[kc]]$y2)))
          zgg$posic_zig <-  max(abs(zgg$posic_zig_a), abs(zgg$posic_zig_b))
          zgg$primerkcore <- FALSE
          x <- c(0)
          y <- c(zgg$ymax)
        }
        if (zgg$df_cores[kc,]$num_species_guild_b>0){
          despl_pointer_y <- zgg$displace_y_b[kc] * zgg$ymax
          if ((kc == 2))
          {
            edge_core <- "yes"
            zgg$pointer_y <- max(4,zgg$species_in_core2_b)*zgg$height_y+abs(max(zgg$list_dfs_a[[zgg$kcoremax]]$y2))
          }
          else {
            edge_core <- "no"
            if (zgg$primerkcore)
              zgg$pointer_y <- zgg$ymax
            else if ((zgg$df_cores$num_species_guild_b[kc-1]>5) & !(zgg$primerkcore)) {
              zgg$pointer_y <- zgg$pointer_y - 0.2*zgg$height_y*zgg$df_cores[kc,]$num_species_guild_b
            }
          }
          pystep <- zgg$height_y
          if (zgg$kcoremax < 5)
            wzig <- (0.8+(kc/zgg$kcoremax))* zgg$width_zig
          else
            wzig <- 1.3*zgg$width_zig
          if (kc == 2)
            wzig <- wzig *min(2,(1+0.1*sqrt(max(zgg$df_cores$num_species_guild_a[kc],zgg$df_cores$num_species_guild_b[kc]))))
          zig <-  draw_individual_ziggurat(paste0("kcore", kc, "-b"), zgg$g, kc, basex = zgg$pointer_x, widthx = wzig,
                                basey = zgg$pointer_y + despl_pointer_y,  ystep = pystep, strlabels = zgg$df_cores$species_guild_b[[kc]],
                                strguild = zgg$str_guild_b, sizelabels = zgg$lsize_zig,
                                colorb = zgg$color_guild_b, numboxes = zgg$df_cores[kc,]$num_species_guild_b,
                                zinverse = "no", edge = edge_core, grafo = p, svg=svg)
          p <- zig["p"][[1]]
          svg <- zig["svg"][[1]]
          zgg$list_dfs_b[[kc]]<- zig["dr"][[1]]
          zgg$last_xtail_b[[kc]] <- max(zgg$list_dfs_b[[kc]]$x2)
          zgg$last_ytail_b[[kc]] <- max(abs(zgg$list_dfs_b[[kc]]$y2))
          zgg$posic_zig_b <- max(zgg$posic_zig_b,max(zgg$list_dfs_b[[kc]]$y2))
          zgg$posic_zig <-  max(abs(zgg$posic_zig_a), abs(zgg$posic_zig_b))
          zgg$primerkcore <- FALSE
        }
      }
    }
  }
  calc_vals <- list("p" = p, "svg" = svg, "posic_zig" = zgg$posic_zig, "list_dfs_a" = zgg$list_dfs_a, "list_dfs_b" = zgg$list_dfs_b,
                    "last_xtail_a" = zgg$last_xtail_a, "last_ytail_a" = zgg$last_ytail_a,
                    "last_xtail_b" = zgg$last_xtail_b, "last_ytail_b" = zgg$last_ytail_b)
  return(calc_vals)
}


draw_maxcore <- function(svg)
{
  kcoremax_label_display <- function (idPrefix, gp,svg,kcoremaxlabel_angle,pdata,plabel,plabelsize,phjust=0, is_guild_a = TRUE) {
    labelcolor <- ifelse(length(zgg$labels_color)>0,zgg$labels_color[2-as.numeric(is_guild_a)], pdata$col_row)
    if (kcoremaxlabel_angle == 0) {
      gp <- gp +  geom_text(data=pdata, aes(x=x1+(x2-x1)/2, y=y1+(y2-y1)/2), label=plabel,
                            color = labelcolor, size=plabelsize, angle = kcoremaxlabel_angle)
      svg$text(idPrefix=idPrefix, data=pdata, mapping=aes(x=x1+(x2-x1)/2, y=y1+(y2-y1)/2), label=plabel, color=labelcolor, size=plabelsize, angle=kcoremaxlabel_angle)
    } else {
      gp <- gp + geom_text(data=pdata, aes(x=x1, y=y1+(y2-y1)/20), label=plabel,
                           color = labelcolor, size=plabelsize, angle = kcoremaxlabel_angle,
                           vjust = 1, hjust = phjust)
      svg$text(idPrefix=idPrefix, data=pdata, mapping=aes(x=x1, y=y1+(y2-y1)/20), label=plabel, color=labelcolor, size=plabelsize, angle=kcoremaxlabel_angle)
    }
    calc_vals <- list("p" = gp, "svg" = svg)
    return(calc_vals)
  }

  zgg$last_ytail_a[zgg$kcoremax]<- zgg$toopy
  zgg$last_xtail_a[zgg$kcoremax]<- zgg$topxa
  zgg$list_dfs_a[[zgg$kcoremax]]<- draw_coremax_triangle(zgg$basex,zgg$topxa,zgg$basey,zgg$toopy,
                                                  zgg$num_a_coremax,zgg$color_guild_a,
                                                  zgg$df_cores$species_guild_a[[zgg$kcoremax]],
                                                  zgg$g, zgg$str_guild_a, orderby = "kdegree")

  nsp <- name_species_preprocess(zgg$kcoremax,zgg$list_dfs_a[[zgg$kcoremax]],zgg$kcore_species_name_display,
                                 zgg$kcore_species_name_break)

  kcoremaxlabel_angle <- nsp$kcoremaxlabel_angle
  labelszig <- nsp$labelszig

  p <- ggplot() +
    scale_x_continuous(name="x") +
    scale_y_continuous(name="y") +
    geom_rect(data=zgg$list_dfs_a[[zgg$kcoremax]], mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill = zgg$list_dfs_a[[zgg$kcoremax]]$col_row,  color="transparent",alpha=zgg$alpha_level)


  svg$rect(paste0("kcore", zgg$kcoremax, "-a"), data=zgg$list_dfs_a[[zgg$kcoremax]],
           mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill=zgg$list_dfs_a[[zgg$kcoremax]]$col_row, alpha=zgg$alpha_level,
           color="transparent", size=0.5)
  f <- kcoremax_label_display(paste0("kcore", zgg$kcoremax, "-a"),p,svg,kcoremaxlabel_angle,zgg$list_dfs_a[[zgg$kcoremax]],labelszig,zgg$lsize_kcoremax)
  p <- f["p"][[1]]
  svg <- f["svg"][[1]]
  num_b_coremax <- zgg$df_cores[zgg$kcoremax,]$num_species_guild_b
  zgg$basey <- - zgg$basey
  zgg$topxb <- zgg$topxa
  zgg$toopy <- - zgg$toopy
  zgg$list_dfs_b[[zgg$kcoremax]] <- draw_coremax_triangle(zgg$basex,zgg$topxb,zgg$basey,zgg$toopy,num_b_coremax,
                                                  zgg$color_guild_b,
                                                  zgg$df_cores$species_guild_b[[zgg$kcoremax]],
                                                  zgg$g, zgg$str_guild_b, orderby = "kdegree")

  zgg$last_ytail_b[zgg$kcoremax]<- zgg$toopy
  zgg$last_xtail_b[zgg$kcoremax]<- zgg$topxb
  nsp <- name_species_preprocess(zgg$kcoremax,zgg$list_dfs_b[[zgg$kcoremax]],zgg$kcore_species_name_display,
                                 zgg$kcore_species_name_break)
  kcoremaxlabel_angle <- nsp$kcoremaxlabel_angle
  labelszig <- nsp$labelszig

  p <- p + geom_rect(data=zgg$list_dfs_b[[zgg$kcoremax]] , mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2),
                     fill = zgg$list_dfs_b[[zgg$kcoremax]]$col_row, color="transparent", alpha=zgg$alpha_level)

  svg$rect(paste0("kcore", zgg$kcoremax, "-b"), zgg$list_dfs_b[[zgg$kcoremax]] , mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2),
           fill=zgg$list_dfs_b[[zgg$kcoremax]]$col_row, alpha=zgg$alpha_level,
           color="transparent", size=0.5)
  f <- kcoremax_label_display(paste0("kcore", zgg$kcoremax, "-b"),p,svg,kcoremaxlabel_angle,zgg$list_dfs_b[[zgg$kcoremax]],labelszig,
                              zgg$lsize_kcoremax, phjust = 1, is_guild_a = FALSE)
  p <- f["p"][[1]]
  svg <- f["svg"][[1]]

  calc_vals <- list("p" = p, "svg" = svg, "basey" = zgg$basey, "topy" = zgg$toopy, "topxa" = zgg$topxa, "topxb" = zgg$topxb,
                    "list_dfs_a" = zgg$list_dfs_a, "list_dfs_b" = zgg$list_dfs_b,
                    "last_xtail_a" = zgg$last_xtail_a, "last_ytail_a" = zgg$last_ytail_a,
                    "last_xtail_b" = zgg$last_xtail_b, "last_ytail_b" = zgg$last_ytail_b)
  return(calc_vals)
}

draw_coremax_tails <- function(p, svg)
{
  long_kcoremax_tail <- FALSE
  leftjump <- (1.2-0.02*nrow(zgg$list_dfs_a[[zgg$kcoremax]])) * zgg$hop_x
  point_x <- zgg$list_dfs_a[[zgg$kcoremax]][nrow(zgg$list_dfs_a[[zgg$kcoremax]]),]$x2 - leftjump
  point_y <- zgg$posic_zig * 0.8
  if (exists("df_orph_a", envir = zgg))
    long_tail_a <- zgg$df_orph_a[(zgg$df_orph_a$kcore == zgg$kcoremax) & (zgg$df_orph_a$repeated == "no"),]
  if ((exists("long_tail_a")) & (zgg$kcoremax > 2))
  {
    if (length(long_tail_a)>5)
      long_kcoremax_tail <- TRUE
    v<-  draw_edge_tails(p,svg,point_x,(point_y+zgg$height_y*zgg$aspect_ratio),zgg$kcoremax,
                         long_tail_a,zgg$list_dfs_b,zgg$color_guild_a, inverse = "yes",
                         vertical = "no", orientation = "South", revanddrop = "yes",
                         pbackground = "no", tspline = "arc", joinchars=zgg$joinstr)
    p <- v["p"][[1]]
    svg <- v["svg"][[1]]
    zgg$last_xtail_b[zgg$kcoremax] <- v["lastx"][[1]]
    if (length(v["lasty"][[1]])>0)
      zgg$last_ytail_b[zgg$kcoremax] <- v["lasty"][[1]]
    else                                                   # only for degenerate networks with all nodes in kcoremax
      zgg$last_ytail_b[zgg$kcoremax] <- zgg$toopy
  }
  leftjump <- (1.2-0.02*nrow(zgg$list_dfs_b[[zgg$kcoremax]]))* zgg$hop_x
  point_x <- zgg$list_dfs_b[[zgg$kcoremax]][nrow(zgg$list_dfs_b[[zgg$kcoremax]]),]$x2 - leftjump
  point_y <- -1* zgg$posic_zig * 0.8
  if (exists("df_orph_b", envir = zgg))
    long_tail_b <- zgg$df_orph_b[(zgg$df_orph_b$kcore == zgg$kcoremax) & (zgg$df_orph_b$repeated == "no"),]
  if ( (exists("long_tail_b")) & (zgg$kcoremax > 2) ){
    if (nrow(long_tail_b)>5)
      long_kcoremax_tail <- TRUE
    tailweight <- 0
    if (nrow(long_tail_b)>0){
      for (h in 1:nrow(long_tail_b))
        tailweight <- tailweight + zgg$result_analysis$matrix[as.numeric(long_tail_b$orph[h]),
                                                                    as.numeric(long_tail_b$partner[h])]
      long_tail_b$weightlink <- get_link_weights(tailweight)
    }

    v <-  draw_edge_tails(p,svg,point_x,point_y,zgg$kcoremax,long_tail_b,zgg$list_dfs_a,zgg$color_guild_b,
                         inverse = "no",
                         vertical = "no", orientation = "North", revanddrop = "yes",
                         pbackground = "no", tspline = "arc", joinchars=zgg$joinstr, is_guild_a = FALSE,
                         wlink = long_tail_b$weightlink[1])
    p <- v["p"][[1]]
    svg <- v["svg"][[1]]
    zgg$last_xtail_a[zgg$kcoremax] <- v["lastx"][[1]]
    if (length(v["lasty"][[1]])>0)
      zgg$last_ytail_b[zgg$kcoremax] <- v["lasty"][[1]]
    else                                                   # only for degenerate networks with all nodes in kcoremax
      zgg$last_ytail_b[zgg$kcoremax] <- zgg$toopy
  }
  calc_vals <- list("p" = p, "svg" = svg,
                    "last_xtail_a" = zgg$last_xtail_a, "last_ytail_a" = zgg$last_ytail_a,
                    "last_xtail_b" = zgg$last_xtail_b, "last_ytail_b" = zgg$last_ytail_b)
  return(calc_vals)
}

display_plot <- function(p, printfile, flip, plwidth=14, plheight=11, ppi = 300, landscape = zgg$label_strguild, fname_append = "")
{
  if (flip)
    p <- p + coord_flip()
  if (printfile){
    if (fname_append != "")
      ftname_append <- paste0("_",fname_append)
    else
      ftname_append <- fname_append
    dir.create(zgg$plotsdir, showWarnings = FALSE)
    if (landscape)
      png(paste0(zgg$plotsdir,"/",zgg$network_name,"_ziggurat",ftname_append,".png"), width=(plwidth*ppi), height=plheight*ppi, res=ppi)
    else
      png(paste0(zgg$plotsdir,"/",zgg$network_name,"_ziggurat",ftname_append,".png"), width=(plheight*ppi), height=plwidth*ppi, res=ppi)
  }
  print(p)
  if (printfile)
    dev.off()
}

strip_isolated_nodes <- function()
{
  lgrados <- igraph::degree(zgg$result_analysis$graph)
  if (sum(lgrados == 0) > 0)
    for (k in 1:length(lgrados))
    {
      if (lgrados[k] == 0){
        zgg$result_analysis$graph <- delete_vertices(zgg$result_analysis$graph,names(lgrados[k]))
        if ( length(grep(zgg$str_guild_b,names(lgrados[k]) )) >0 )
          zgg$result_analysis$num_guild_b <<- zgg$result_analysis$num_guild_b -1
        else
          zgg$result_analysis$num_guild_a <<- zgg$result_analysis$num_guild_a -1
      }
    }
}


read_and_analyze <- function(directorystr,network_file,label_strguilda,label_strguildb)
{

  str_guild_a <- "pl"
  str_guild_b <- "pol"
  name_guild_a <- "Plants"
  name_guild_b <- "Pollinators"
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
                                     guild_b = str_guild_b, only_NODF = TRUE)


  calc_vals <- list("result_analysis" = result_analysis, "str_guild_a" = str_guild_a, "str_guild_b" = str_guild_b,
                    "name_guild_a" = name_guild_a, "name_guild_b" = name_guild_b,
                    "network_name" = network_name)
  return(calc_vals)
}

def_configuration <- function(paintlinks, print_to_file, plotsdir, flip_results, aspect_ratio,
                              alpha_level, color_guild_a, color_guild_b,
                              color_link, alpha_link, size_link,
                              displace_y_b, displace_y_a, lsize_kcoremax, lsize_zig, lsize_kcore1,
                              lsize_legend, lsize_core_box, labels_color,
                              height_box_y_expand, kcore2tail_vertical_separation,  kcore1tail_disttocore,
                              innertail_vertical_separation ,
                              factor_hop_x, fattailjumphoriz, fattailjumpvert,
                              coremax_triangle_height_factor, coremax_triangle_width_factor,
                              paint_outsiders, displace_outside_component,
                              outsiders_separation_expand, outsiders_legend_expand, specialistskcore2_horizontal_dist_rootleaf_expand,
                              specialistskcore2_vertical_dist_rootleaf_expand , specialists_boxes_separation_count,
                              root_specialist_expand,hide_plot_border,rescale_plot_area,kcore1specialists_leafs_vertical_separation,
                              corebox_border_size, kcore_species_name_display,kcore_species_name_break,
                              shorten_species_name,exclude_species_number,
                              label_strguilda, label_strguildb, landscape_plot, backg_color, show_title,
                              use_spline, spline_points, file_name_append, svg_scale_factor, weighted_links, square_nodes_size_scale,
                              move_all_SVG_up, progress
                              )
{
  # ENVIRONMENT CONFIGURATION PARAMETERS
  zgg$paintlinks <- paintlinks
  zgg$print_to_file <- print_to_file
  zgg$plotsdir <- plotsdir
  zgg$flip_results <- flip_results
  zgg$alpha_level <- alpha_level
  zgg$color_guild_a <- color_guild_a
  zgg$color_guild_b <- color_guild_b
  zgg$color_link <- color_link
  zgg$alpha_link <- alpha_link
  zgg$size_link <- size_link
  zgg$displace_y_b <- displace_y_b
  zgg$displace_y_a <- displace_y_a
  zgg$aspect_ratio <- aspect_ratio
  zgg$labels_size <- 3.5
  zgg$lsize_kcoremax <- lsize_kcoremax
  zgg$lsize_zig <- lsize_zig
  zgg$lsize_kcore1 <- lsize_kcore1
  zgg$lsize_legend <- lsize_legend
  zgg$lsize_core_box <- lsize_core_box
  zgg$labels_color <- labels_color
  zgg$height_box_y_expand <- height_box_y_expand
  zgg$kcore2tail_vertical_separation <- kcore2tail_vertical_separation                 # Vertical separation of orphan boxes linked to core 2 in number of heights_y
  zgg$kcore1tail_disttocore <- kcore1tail_disttocore                            # Horizontal & Vertical distances of edge/specialist tails linked to core 1 North & South
  zgg$innertail_vertical_separation <- innertail_vertical_separation                  # Vertical separation of orphan boxes linked to inner cores in number of heights_y
  #zgg$horiz_kcoremax_tails_expand <- horiz_kcoremax_tails_expand                  # horizontal separation of edge tails connected to kcoremax.
  zgg$factor_hop_x <- factor_hop_x
  #zgg$displace_legend <- displace_legend   DEPRECATED
  zgg$fattailjumphoriz <- fattailjumphoriz
  zgg$fattailjumpvert <- fattailjumpvert
  zgg$coremax_triangle_height_factor <- coremax_triangle_height_factor
  zgg$coremax_triangle_width_factor <- coremax_triangle_width_factor
  zgg$paint_outsiders <- paint_outsiders
  zgg$displace_outside_component <- displace_outside_component
  zgg$outsiders_separation_expand <- outsiders_separation_expand
  zgg$outsiders_legend_expand <- outsiders_legend_expand
  zgg$specialistskcore2_horizontal_dist_rootleaf_expand <- specialistskcore2_horizontal_dist_rootleaf_expand        # Controls the distance of specialist root leaves to partner in core 2
  zgg$specialistskcore2_vertical_dist_rootleaf_expand <- specialistskcore2_vertical_dist_rootleaf_expand
  zgg$specialists_boxes_separation_count <- specialists_boxes_separation_count                  # Separation of leaves of a specialist tail
  zgg$root_specialist_expand <- root_specialist_expand
  zgg$hide_plot_border <- hide_plot_border
  zgg$rescale_plot_area <- rescale_plot_area
  zgg$kcore1specialists_leafs_vertical_separation <- kcore1specialists_leafs_vertical_separation
  zgg$corebox_border_size <- corebox_border_size
  zgg$kcore_species_name_display <- kcore_species_name_display
  zgg$kcore_species_name_break <- kcore_species_name_break
  zgg$shorten_species_name <- shorten_species_name
  zgg$exclude_species_number <- exclude_species_number
  zgg$landscape_plot <- landscape_plot
  zgg$backg_color <- backg_color
  zgg$show_title <- show_title
  zgg$use_spline <- use_spline
  zgg$spline_points <- spline_points
  zgg$file_name_append <- file_name_append
  zgg$svg_scale_factor <- svg_scale_factor
  zgg$weighted_links <- weighted_links
  zgg$square_nodes_size_scale <- square_nodes_size_scale
  zgg$move_all_SVG_up <- move_all_SVG_up
  zgg$progress <- progress
}

init_working_values <- function()
{
  zgg$joinstr <- " "
  zgg$max_position_y_text_core <- 0
  zgg$rg <- V(zgg$result_analysis$graph)
  zgg$g <- zgg$rg[zgg$rg$kradius != Inf]
  zgg$outsider <- zgg$rg[zgg$rg$kradius == Inf]
  zgg$outsiders_a <- zgg$outsider$name[grep(zgg$str_guild_a,zgg$outsider$name)]
  zgg$outsiders_b <- zgg$outsider$name[grep(zgg$str_guild_b,zgg$outsider$name)]
  zgg$ind_cores <- rev(sort(unique(zgg$g$kcorenum)))
  zgg$kcoremax <- max(zgg$ind_cores)
  palcores <- colorRampPalette(c("salmon3","salmon4"))
  zgg$corecols <- palcores(zgg$kcoremax)
  zgg$last_xtail_a <- rep(NA,zgg$kcoremax)
  zgg$last_ytail_a <- rep(NA,zgg$kcoremax)
  zgg$last_xtail_b <- rep(NA,zgg$kcoremax)
  zgg$last_ytail_b <- rep(NA,zgg$kcoremax)
  species_guild_a <- rep(NA,zgg$kcoremax)
  species_guild_b <- rep(NA,zgg$kcoremax)
  num_species_guild_a <- rep(NA,zgg$kcoremax)
  num_species_guild_b <- rep(NA,zgg$kcoremax)
  zgg$df_cores <- data.frame(species_guild_a, species_guild_b, num_species_guild_a, num_species_guild_b)
  zgg$list_dfs_a <- list()
  zgg$list_dfs_b <- list()
  zgg$df_cores$num_species_guild_a <- 0
  zgg$df_cores$num_species_guild_b <- 0
  zgg$straight_links <- data.frame(x1=c(), x2 = c(), y1 = c(),  y2 = c(), weightlink = c())
  zgg$bent_links <- data.frame(x=c(), y = c(),  number = c(), weightlink = c() )
  zgg$count_bent_links <- 0
}


draw_ziggurat_plot <- function(svg_scale_factor, progress)
{
  if (!is.null(progress)) progress$inc(1/11, detail=strings$value("MESSAGE_ZIGGURAT_PROGRESS_PROCESSING_NODES"))
  zinit_time <- proc.time()
  for (i in zgg$ind_cores) {
    nodes_in_core_a <- zgg$g[(zgg$g$guild == zgg$str_guild_a)&(zgg$g$kcorenum == i)]$name
    nodes_in_core_b <- zgg$g[(zgg$g$guild == zgg$str_guild_b)&(zgg$g$kcorenum == i)]$name
    zgg$df_cores[i,]$species_guild_a <- list(unlist(lapply(nodes_in_core_a, function(x) strsplit(x,zgg$str_guild_a)[[1]][[2]])))
    zgg$df_cores[i,]$num_species_guild_a <- length(nodes_in_core_a)
    zgg$df_cores[i,]$species_guild_b <- list(unlist(lapply(nodes_in_core_b, function(x) strsplit(x,zgg$str_guild_b)[[1]][[2]])))
    zgg$df_cores[i,]$num_species_guild_b <- length(nodes_in_core_b)
  }
  zgg$num_a_coremax <- zgg$df_cores[zgg$kcoremax,]$num_species_guild_a
  base_width <- 2000
  zgg$ymax <- 2*base_width/zgg$aspect_ratio
  zgg$tot_width <- zgg$ymax * (1+0.25*(zgg$kcoremax - 2))
  zgg$species_in_core2_a <- sum(zgg$df_cores[2,]$num_species_guild_a)
  zgg$species_in_core2_b <- sum(zgg$df_cores[2,]$num_species_guild_b)
  zgg$species_in_almond_a <- sum(zgg$df_cores[2:(zgg$kcoremax-1),]$num_species_guild_a)
  zgg$species_in_almond_b <- sum(zgg$df_cores[2:(zgg$kcoremax-1),]$num_species_guild_b)
  zgg$height_y <- zgg$ymax/max(1.3,(1.3*max(zgg$species_in_almond_a,zgg$species_in_almond_b)))
  maxincore2 <- max(zgg$species_in_core2_a,zgg$species_in_core2_b)
  if (zgg$kcoremax < 4)
    if (zgg$species_in_core2_a+zgg$species_in_core2_b < 6)
      zgg$height_y <- (0.08)*zgg$ymax
  zgg$yoffset <- zgg$height_y*maxincore2*zgg$height_box_y_expand
  fmult <- (zgg$ymax+zgg$yoffset)/zgg$ymax
  zgg$ymax <- zgg$ymax + zgg$yoffset
  zgg$tot_width <- zgg$tot_width*fmult
  zgg$height_y <- zgg$height_y * fmult * zgg$height_box_y_expand
  zgg$yoffset <- zgg$height_y*maxincore2
  zgg$ymax <- zgg$ymax * (1+0.1*zgg$height_box_y_expand)
  for (i in seq(3,zgg$kcoremax-1)){
    zgg$displace_y_a[i] <- zgg$displace_y_a[i] + zgg$coremax_triangle_height_factor*zgg$species_in_core2_a*zgg$height_y/zgg$ymax
    zgg$displace_y_b[i] <- zgg$displace_y_b[i] + zgg$coremax_triangle_height_factor*zgg$species_in_core2_b*zgg$height_y/zgg$ymax
  }

  zgg$hop_x <- zgg$factor_hop_x*(zgg$tot_width)/max(1,(zgg$kcoremax-2))
  zgg$lado <- min(0.05*zgg$tot_width,zgg$height_y * zgg$aspect_ratio)
  zgg$basey <- (0.1+0.1*length(zgg$df_cores[zgg$kcoremax,]$num_species_guild_a))*zgg$ymax
  wcormax <- 1.2*zgg$hop_x*zgg$coremax_triangle_width_factor
  zgg$topxa <- 0.65*zgg$hop_x
  zgg$basex <- zgg$topxa - wcormax
  zgg$posic_zig <- 0
  zgg$posic_zig_a <- 0
  zgg$posic_zig_b <- 0
  zgg$toopy <- 0.3*zgg$ymax+zgg$basey
  zgg$strips_height <- 0.6*(zgg$ymax-zgg$yoffset)/max(1,(zgg$kcoremax-2))
  # Draw max core triangles
  svg <-SVG(svg_scale_factor)
  if (!is.null(progress)) progress$inc(1/11, detail=strings$value("MESSAGE_ZIGGURAT_PROGRESS_DRAWING_MAXCORE"))
  f <- draw_maxcore(svg)
  p <- f["p"][[1]]
  svg <- f["svg"][[1]]

  zgg$posic_zig <- f["posic_zig"][[1]]
  zgg$list_dfs_a <- f["list_dfs_a"][[1]]
  zgg$list_dfs_b <- f["list_dfs_b"][[1]]
  zgg$last_xtail_a <- f["last_xtail_a"][[1]]
  zgg$last_ytail_a <- f["last_ytail_a"][[1]]
  zgg$last_xtail_b <- f["last_xtail_b"][[1]]
  zgg$last_ytail_b <- f["last_ytail_b"][[1]]
  zgg$topy <- f["topy"][[1]]
  zgg$topxa <- f["topxa"][[1]]
  zgg$topxb <- f["topxb"][[1]]
  # Draw inner almond ziggurats
  zgg$pointer_x <- max(zgg$topxa, zgg$topxb)+zgg$hop_x
  zgg$pointer_y <- zgg$ymax+zgg$height_y*max(zgg$df_cores$num_species_guild_a[zgg$kcoremax-1],zgg$df_cores$num_species_guild_b[zgg$kcoremax-1])
  zgg$width_zig <- 0.3*zgg$hop_x
  zgg$primerkcore <- TRUE
  if (!is.null(progress)) progress$inc(1/11, detail=strings$value("MESSAGE_ZIGGURAT_PROGRESS_DRAWING_ZIGGURAT"))
  f <- draw_all_ziggurats(p, svg)
  p <- f["p"][[1]]
  svg <- f["svg"][[1]]

  zgg$posic_zig <- f["posic_zig"][[1]]
  zgg$list_dfs_a <- f["list_dfs_a"][[1]]
  zgg$list_dfs_b <- f["list_dfs_b"][[1]]
  zgg$last_xtail_a <- f["last_xtail_a"][[1]]
  zgg$last_ytail_a <- f["last_ytail_a"][[1]]
  zgg$last_xtail_b <- f["last_xtail_b"][[1]]
  zgg$last_ytail_b <- f["last_ytail_b"][[1]]
  # Add kcore1 information
  if (!is.null(zgg$df_cores[1,])){
    if (zgg$df_cores[1,]$num_species_guild_a > 0)
      zgg$list_dfs_a[[1]] <- conf_kcore1_info(zgg$str_guild_a)
    if (zgg$df_cores[1,]$num_species_guild_b > 0)
      zgg$list_dfs_b[[1]] <- conf_kcore1_info(zgg$str_guild_b)
  }
  # Add outsiders information
  info_out_a <- conf_outsiders_info(zgg$str_guild_a)
  info_out_b <- conf_outsiders_info(zgg$str_guild_b)
  zgg$list_dfs_a[[1]] <- rbind(zgg$list_dfs_a[[1]], info_out_a)
  zgg$list_dfs_b[[1]] <- rbind(zgg$list_dfs_b[[1]], info_out_b)
  # Draw core boxex
  for (i in seq(zgg$kcoremax,2))
    if ((length(zgg$list_dfs_a[[i]])+length(zgg$list_dfs_b[[i]]))>0){
      v <- draw_core_box(p, svg, i)
      p <- v[["p"]]
      svg <- v[["svg"]]
      zgg$max_position_y_text_core <- v[["max_position_y_text_core"]]
    }

  # Hanlde orphans, species outside the ziggurat
  if (!is.null(progress)) progress$inc(1/11, detail=strings$value("MESSAGE_ZIGGURAT_PROGRESS_DRAWING_ORPHANS"))
  f <- handle_orphans(zgg$result_analysis$graph)
  zgg$mtxlinks <- f["mtxlinks"][[1]]
  zgg$orphans_a <- f["orphans_a"][[1]]
  zgg$orphans_b <- f["orphans_b"][[1]]
  zgg$df_orph_a <- f["df_orph_a"][[1]]
  zgg$df_orph_b <- f["df_orph_b"][[1]]
  # Species of core 1 linked to max core (except the most generalist)
  zgg$gap <-  4*zgg$height_y
  if (!is.null(progress)) progress$inc(1/11, detail=strings$value("MESSAGE_ZIGGURAT_PROGRESS_DRAWING_MAXCORE_TAILS"))
  f <- draw_coremax_tails(p, svg)
  p <- f["p"][[1]]
  svg <- f["svg"][[1]]

  zgg$last_xtail_a <- f["last_xtail_a"][[1]]
  zgg$last_ytail_a <- f["last_ytail_a"][[1]]
  zgg$last_xtail_b <- f["last_xtail_b"][[1]]
  zgg$last_ytail_b <- f["last_ytail_b"][[1]]
  # Fat tails - nodes of core 1 linked to most generalist of opposite guild. Left side of panel
  z <- handle_fat_tails(p, svg)
  p <- z["p"][[1]]
  svg <- z["svg"][[1]]
  zgg$pos_tail_x <- z["pos_tail_x"][[1]]
  # Nodes of core 1 linked to species in cores kcoremax-1 to core 2.
  if (!is.null(progress)) progress$inc(1/11, detail=strings$value("MESSAGE_ZIGGURAT_PROGRESS_DRAWING_INNER_ORPHANS"))
  z <- draw_inner_orphans(p, svg)
  p <- z["p"][[1]]
  svg <- z["svg"][[1]]

  # Draw inner links
  if (!is.null(progress)) progress$inc(1/11, detail=strings$value("MESSAGE_ZIGGURAT_PROGRESS_DRAWING_INNER_LINKS"))
  if (zgg$paintlinks) {
    z <- draw_inner_links(p, svg)
    p <- z["p"][[1]]
    svg <- z["svg"][[1]]
  }

  # specialists management
  if (!is.null(progress)) progress$inc(1/11, detail=strings$value("MESSAGE_ZIGGURAT_PROGRESS_DRAWING_specialistS"))
  v <- handle_specialists(p,svg,specialists_a,specialists_b,zgg$lado,zgg$gap)
  p <- v["p"][[1]]
  svg <- v["svg"][[1]]

  zgg$df_chains <- v["df_chains"][[1]]
  # Specied outside the giant component
  if (!is.null(progress)) progress$inc(1/11, detail=strings$value("MESSAGE_ZIGGURAT_PROGRESS_DRAWING_OUTSIDERS"))
  if (zgg$paint_outsiders) {
    v <- handle_outsiders(p,svg,outsiders,zgg$df_chains)
    p <- v["p"][[1]]
    svg <- v["svg"][[1]]
  }

  # Legend, title and final annotations
  v <- write_annotations(p, svg)
  p <- v["p"][[1]]
  svg <- v["svg"][[1]]

  # Plot straight links
  if (!is.null(progress)) progress$inc(1/11, detail=strings$value("MESSAGE_ZIGGURAT_PROGRESS_DRAWING_LINKS"))
  if (zgg$paintlinks){
    if (nrow(zgg$straight_links)>0) {
      p <- p+ geom_segment(data=zgg$straight_links, aes(x=x1, y=y1, xend=x2, yend=y2),
                           linewidth=zgg$straight_links$weightlink, color=zgg$color_link ,alpha=zgg$alpha_link)
      factormult <- 0.1*svg_scale_factor
      svg$segment(idPrefix="link", data=zgg$straight_links, mapping=aes(x=x1, y=y1, xend=x2, yend=y2),
                  alpha=zgg$alpha_link, color=zgg$color_link,
                  #size=factormult*zgg$straight_links$weightlink)
                  size=zgg$straight_links$weightlink)
    }
    if (nrow(zgg$bent_links)>0) {
      p <- p + geom_path(data =zgg$bent_links,aes(x,y,group=number), linewidth=zgg$bent_links$weightlink,
                       color=zgg$color_link ,alpha=zgg$alpha_link)
      svg$path(idPrefix="link", data=zgg$bent_links, mapping=aes(x, y, group=number), alpha=zgg$alpha_link,
                      color=zgg$color_link,
               #size=0.1*svg_scale_factor*zgg$bent_links$weightlink)
               size=zgg$bent_links$weightlink)
    }
  }
  if (is.null(progress))
    display_plot(p,zgg$print_to_file,zgg$flip_results, landscape = zgg$landscape_plot, fname_append = zgg$file_name_append)


  # Stores results
  zgg$plot  <- p
  zgg$svg   <- svg

  if (!is.null(progress)) progress$inc(0, detail=strings$value("MESSAGE_ZIGGURAT_PROGRESS_DONE"))

  return(zgg)
}

# ziggurat_graph("../data/","example-try-this-first.csv")
