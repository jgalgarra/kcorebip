debugging = FALSE
if (debugging){
   source("network-kanalysis.R")
   source("SVG.R")
}

#' Plotting a bipartite graph
#'
#' This function plots the ziggurat graph of a bipartite network. Configuration parameters and
#' results are stored in a global environment called bpp. This environment is not destroyed
#' after the function is executed, so the developer can store it in a configuration file, retrieve
#' network analysis variables and so on. If a new bipartite_graph is called, bpp is destroyed and
#' created again for the new plot. Plotting options are explained in the user manual.
#'
#' @param datadir the name of the file of the interaction matrix
#' @param filename the file with the interaction matrix
#' @param print_to_file if set to FALSE the plot is displayed in the R session window
#' @param plotsdir the directory where the plot is stored
#' @param orderkcoremaxby sets order of kcoremax nodes, by kradius or kdegree
#' @param style bipartite representation style: legacy, kcoreorder, chilopod
#' @param guild_gap_increase controls the disptance between guild rows
#' @param flip_results displays the graph in portrait configuration
#' @param aspect_ratio ziggurat plot default aspect ratio
#' @param alpha_level transparency for ziggurats' filling
#' @param color_guild_a default filling for nodes of guild_a
#' @param color_guild_b default filling for nodes of guild_b
#' @param color_link default links color
#' @param alpha_link link transparency
#' @param size_link width of the links
#' @param lsize_kcoremax nodes in kshell max label size
#' @param lsize_zig nodes in inner ziggurats label size
#' @param lsize_kcore1 labels of nodes in kshell 1
#' @param lsize_legend legend label size
#' @param labels_color default label colors
#' @param factor_hop_x expand inner ziggurats horizontal distance
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
#' @param show_legend show plot legend position
#' @param use_spline use splines to draw links
#' @param spline_points number of points for each spline
#' @param file_name_append a label that the user may append to the plot file name for convenience
#' @param svg_scale_factor only for interactive apps, do not modify
#' @param weighted_links function to add link weight: 'none', 'log10' , 'ln', 'sqrt'
#' @param square_nodes_size_scale scale nodes area of kcore1 and outsiders
#' @param move_all_SVG_up move up all the SVG plot by this fraction, useful to crop upper white space
#' @param move_all_SVG_right move rightwards all the SVG plot by this fraction, useful to crop right white space
#' @param progress only for interactive apps, do not modifiy
#' @export
#' 
#' @examples bipartite_graph("data/","M_PL_001.csv",plotsdir="grafresults/",print_to_file = TRUE)

bipartite_graph <- function(datadir,filename,
                           paintlinks = TRUE, print_to_file = FALSE, plotsdir ="plot_results/ziggurat/", 
                           orderkcoremaxby = "kradius", style="legacy", guild_gap_increase = 1, 
                           flip_results = FALSE, aspect_ratio = 1,
                           alpha_level = 0.2, color_guild_a = c("#4169E1","#00008B"), color_guild_b = c("#F08080","#FF0000"),
                           color_link = "slategray3", alpha_link = 0.5, size_link = 0.5,
                           lsize_kcoremax = 3.5, lsize_zig = 3, lsize_kcore1 = 2.5, lsize_legend = 4, lsize_core_box = 2.5,
                           labels_color = c(),
                           factor_hop_x = 1, 
                           paint_outsiders = TRUE, displace_outside_component = c(0,0),
                           outsiders_separation_expand = 1, outsiders_legend_expand = 1,
                           specialistskcore2_horizontal_dist_rootleaf_expand = 1,
                           specialistskcore2_vertical_dist_rootleaf_expand = 0, specialists_boxes_separation_count = 1,
                           root_specialist_expand = c(1,1), hide_plot_border = TRUE, rescale_plot_area = c(1,1),
                           kcore1specialists_leafs_vertical_separation = 1, corebox_border_size = 0.2,
                           kcore_species_name_display = c(), kcore_species_name_break = c(),
                           shorten_species_name = 0, exclude_species_number = FALSE, label_strguilda = "",
                           label_strguildb = "", landscape_plot = TRUE,
                           backg_color = "white", show_title = TRUE, show_legend = 'TOP', use_spline =TRUE, spline_points = 10,
                           file_name_append = "", svg_scale_factor= 10, weighted_links = "none",
                           square_nodes_size_scale = 1, move_all_SVG_up = 0, move_all_SVG_right = 0,
                           progress=NULL
)
{
  # This assignment stores the call parameters in ziggurat_argg as a list. This list is useful
  # to save plotting parameters for a future simulation


  bipartite_argg <- c(as.list(environment()))
  # Create global environment
  bpp <<- new.env()
  fsvgtext <<- 5
  if (!is.null(progress)) 
    progress$inc(1/11, detail=strings$value("MESSAGE_ZIGGURAT_PROGRESS_ANALYZING_NETWORK"))
  # Analyze network
  f <- kcorebip:::read_and_analyze(datadir,filename,label_strguilda, label_strguildb)
  bpp$result_analysis <- f["result_analysis"][[1]]
  bpp$str_guild_a <- f["str_guild_a"][[1]]
  bpp$str_guild_b <- f["str_guild_b"][[1]]
  bpp$name_guild_a <- f["name_guild_a"][[1]]
  bpp$name_guild_b <- f["name_guild_b"][[1]]
  bpp$network_name <- f["network_name"][[1]]
  bpp$mtxlinks <- data.frame(igraph::as_edgelist(bpp$result_analysis$graph))
  names(bpp$mtxlinks) <- c("guild_a","guild_b")
  
  # Exit if kcore max == 1
  if (bpp$result_analysis$max_core == 1){
    msg = "Max core is 1. Ziggurat plot only works if max core is bigger than 1"
    if (!is.null(progress))
      progress$inc(1/11, detail=strings$value(msg))
    else
      print(msg)
    return(bpp)
  }
  # Copy input parameters to the bpp environment
  def_configuration_bip(paintlinks, print_to_file, plotsdir, orderkcoremaxby, style,
                    guild_gap_increase, flip_results, aspect_ratio,
                    alpha_level, color_guild_a, color_guild_b,
                    color_link, alpha_link, size_link,
                    lsize_kcoremax, lsize_zig, lsize_kcore1,
                    lsize_legend, lsize_core_box, labels_color,
                    factor_hop_x,
                    paint_outsiders, displace_outside_component,
                    outsiders_separation_expand, outsiders_legend_expand, specialistskcore2_horizontal_dist_rootleaf_expand,
                    specialistskcore2_vertical_dist_rootleaf_expand , specialists_boxes_separation_count,
                    root_specialist_expand, hide_plot_border, rescale_plot_area,kcore1specialists_leafs_vertical_separation,
                    corebox_border_size, kcore_species_name_display,kcore_species_name_break,shorten_species_name,exclude_species_number,
                    label_strguilda, label_strguildb, landscape_plot, backg_color, show_title, show_legend,
                    use_spline, spline_points, file_name_append, svg_scale_factor, weighted_links,
                    square_nodes_size_scale, move_all_SVG_up, move_all_SVG_right, progress
                    )
  # Removes nodes without any tie. This is not usual in input files but happens
  # when performing destruction simulations
  strip_isolated_nodes(bpp)
  init_working_values(bpp)
  draw_bipartite_plot(svg_scale_factor, progress)
  bpp$bipartite_argg <- bipartite_argg
  return(bpp)
  
}

gen_vert_label <- function(nodes, joinchars = "\n")
{
  nnodes <- length(nodes)
  if (nnodes == 1)
    return(nodes)
  ssal <- ""
  for (i in 1:(nnodes-1))
      ssal <- paste0(ssal,nodes[i],ifelse(i %% 2 == 0, "\n",joinchars))
  ssal <- paste0(ssal,nodes[nnodes])
  ssal <- gsub("  "," ",ssal)
  return(ssal)
}


# Draw tail, species or set of species of 1-shell connected to higher k-index species
draw_tail_bip <- function(idPrefix, p,svg,fat_tail,lado,color,sqlabel,basex,basey,gap,
                      lxx2=0,lyy2=0,sqinverse = "no",
                      position = "West", background = "no",
                      first_leaf = "yes", spline = "no",
                      psize = bpp$lsize_kcore1, is_guild_a = TRUE, wlink=1, 
                      style = "ziggurat", lvp = 0)
{
 
  adjust <- "yes"
  lvjust <- 0
  lhjust <- 0
  langle <- 0
  ecolor <- "transparent"
  bgcolor <- color
  labelcolor <- ifelse(length(bpp$labels_color)>0,bpp$labels_color[2-as.numeric(is_guild_a)],color)
  palpha <- bpp$alpha_link
  sidex <- lado
  paintsidex <- sidex
  signo <- 1
  yy <- abs(basey)
  plxx2 <- lxx2
  plyy2 <- lyy2
  # Lower half plane of the ziggurat
  if (sqinverse=="yes")
    signo <- -1
  # Tails connected eastwards to inner ziggurats. Fat tails
  if (position == "West"){
    adjust = "yes"
    lhjust <- ifelse(bpp$flip_results, 0, 0.5)
    lvjust = ifelse(bpp$flip_results, 1, 0.5)
    if (style=="chilopod"){
      xx <- basex-gap
      posxx1 <- xx+bpp$xstep
      posyy1 = plyy2
    } else {
      xx <- basex-gap
      posxx1 <- xx+sidex
      posyy1 = signo*(yy)+signo*(0.5*gap/(bpp$aspect_ratio))
    }
  }
  # Group of species linked to the highest kdegree species in max shell
  else if (position == "East"){
      gap <- bpp$hop_x/2
      xx <- basex+gap
      posxx1 <- xx
      posyy1 = signo*(yy+sidex/(2*bpp$aspect_ratio))
  }
  # Tails connected to other nodes of the max shell
  else if ((position == "North") |(position == "South")) {
    xx <- basex
    posxx1 <- xx
    posyy1 = signo*(yy)
    if (style=="chilopod"){
      adjust = "yes"
      lvjust = ifelse(bpp$flip_results, 1, 0)
      lhjust = ifelse(bpp$flip_results, 0.5, 0)
    }
  }

  if (background == "no")
  {
    ecolor <- "transparent"
    if (bpp$alpha_level != 1)
      palpha <- max(bpp$alpha_level-0.09,0)
      else if (position == "North")
      {
        langle <- rot_angle
      }
      else if (position == "South")
      {
        langle <- -rot_angle
      }
      else if (position == "West"){
        adjust <- "yes"
      }
  }

  f <- draw_square(idPrefix, p,svg, xx,yy,
                   ifelse(bpp$style=="chilopod",lado,
                   paintsidex*sqrt(bpp$square_nodes_size_scale)),
                   bgcolor,palpha,labelcolor,langle,lhjust,lvjust,
                   slabel=sqlabel,lbsize = 0.8*bpp$lsize_kcoremax,
                   SVGtextfactor=fsvgtext, inverse = sqinverse,
                   adjustoxy = adjust, edgescolor = ecolor)
  p <- f["p"][[1]]
  svg <- f["svg"][[1]]
  # Add tail link
  if (bpp$paintlinks){
    if ((position == "North") |(position == "South"))
      posxx1 = posxx1+sidex/2
    add_link(xx1=posxx1, xx2 = plxx2,
                   yy1 = posyy1, yy2 = plyy2,
                   slink = bpp$size_link*wlink, clink = c(bpp$color_link),
                   alpha_l = bpp$alpha_link, myenv=bpp)
  }
  calc_vals <- list("p" = p, "svg" = svg, "sidex" = sidex, "xx" = posxx1, "yy" = posyy1)
  return(calc_vals)
}

draw_edge_tails_bip <- function(p,svg,point_x,point_y,kcoreother,long_tail,list_dfs,color_guild, inverse = "no",
                            vertical = "yes", orientation = "South", revanddrop = "no",
                            pbackground = "yes", joinchars = "\n", tspline = "no", 
                            is_guild_a = TRUE, wlink = 1)
{

  rxx <- point_x
  ryy <- point_y
  bpp$joinstr <- joinchars
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
  separacion <- 0.035*bpp$tot_width
  last_vertical_position <- 0
  for (i in list_spec)
  {
    conn_species <- which(long_tail$partner == i)
    if (length(conn_species)>0)
    {
      little_tail <- long_tail[long_tail$partner == i,]
      
      data_row <- list_dfs[[kcoreother]][which(list_dfs[[kcoreother]]$label == i),]
      xx2 <- (data_row$x2+data_row$x1)/2
      rxx <- data_row$x1
      yy2 <- data_row$y2
      yfactor <- min(5,nrow(little_tail))
      if (nrow(little_tail>1) & last_vertical_position == yfactor)
        yfactor <- yfactor+1
      ryy <- data_row$y2 + ifelse(yy2>0, yfactor*bpp$xstep, -yfactor*bpp$xstep )
      tailweight <- 0
      for (h in 1:nrow(little_tail))
        if (is_guild_a)
          tailweight <- tailweight + bpp$result_analysis$matrix[as.numeric(little_tail$partner[h]),
                                                                        as.numeric(little_tail$orph[h])]
        else
          tailweight <- tailweight + bpp$result_analysis$matrix[as.numeric(little_tail$orph[h]),
                                                                as.numeric(little_tail$partner[h])]
      little_tail$weightlink <- get_link_weights(tailweight, myenv=bpp)

      v<- draw_tail_bip(paste0(ifelse(is_guild_a, "edge-kcore1-a-", "edge-kcore1-b-"), i),
                    p,svg,little_tail,0.95*bpp$xstep,color_guild[2],
                    gen_vert_label(little_tail$orph,joinchars = " "),
                    rxx,ryy,bpp$gap,lxx2 = xx2,
                    lyy2 = yy2, sqinverse = inverse, position = orientation,
                    background = pbackground, spline = tspline, psize = bpp$lsize_kcore1,
                    is_guild_a = is_guild_a, wlink = little_tail$weightlink[1],style=bpp$style,
                    lvp = last_vertical_position)
      p <- v["p"][[1]]
      svg <- v["svg"][[1]]
      rxx <- v["xx"][[1]]
      ryy <- v["yy"][[1]]
      bpp$landmark_top <- max(bpp$landmark_top,ryy)
      bpp$landmark_bottom <- min(bpp$landmark_bottom,-ryy)
      last_vertical_position <- ifelse(yfactor==1, 0, yfactor)
      if (vertical == "yes"){
        salto <- v["sidex"][[1]]/bpp$aspect_ratio
        point_y <- point_y + 1.4*signo*salto
        rxx <- point_x
      }

      # tails connected to kcoremax except first species
      else{
        if (orientation == "West")
          salto <- 0
        else
          salto <- 0.4*v["sidex"][[1]]/bpp$aspect_ratio
        point_x <- point_x #- separacion - v["sidex"][[1]]
        point_y <- point_y #- 1.4*signo*salto
        ryy <- point_y
        rxx <- point_x
      }
    }
    else
      last_vertical_position <- 0
    m <- m +1
  }

  calc_vals <- list("p" = p, "svg" = svg, "lastx" = rxx, "lasty" = ryy)
  return(calc_vals)
}

# Draw the main guild line
draw_parallel_guilds <- function(basex,topx,basey,topy,numboxes,nnodes,fillcolor,strlabels,
                                 igraphnet,strguild,orderby = "kradius",style="legacy",guild="A")
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
  pbasex <- bpp$coremax_triangle_width_factor*( basex - (nnodes %/%8) * abs(topx-basex)/3)
  xstep <- bpp$square_nodes_size_scale*(topx-pbasex)/max(12,nnodes)
  bpp$xstep <- xstep
  vertsep <- 4
  if ((xstep>2000) && (guild=="A")) {
    vertsep <- 5
    bpp$lsize_kcoremax <-  bpp$lsize_kcoremax + 1
  }
  if ((xstep<1000) && (guild=="A")) {
    vertsep <- 2
    bpp$lsize_kcoremax <-  max(2,bpp$lsize_kcoremax - 1.5)
  }
  

  ptopy <- vertsep*basey+ifelse(basey>0,1,-1)*xstep
  bpp$landmark_bottom <- min(bpp$landmark_bottom,-ptopy)
  bpp$landmark_top <- max(bpp$landmark_top,ptopy)
  ystep <- 0
  
  for (j in (1:numboxes))
  {
    x1 <- c(x1, pbasex+(j-1)*xstep)
    x2 <- c(x2, x1[j]+0.9*bpp$xstep)
    y1 <- c(y1, vertsep*basey)
    y2 <- c(y2, ptopy)
    r <- c(r,j)
    col_row <- c(col_row,fillcolor[1+j%%2])
    kdegree <- c(kdegree,0)
    kradius <- c(kradius,1)
    name_species <- c(name_species,"")
  }
  d1 <- data.frame(x1, x2, y1, y2, r, col_row, kdegree, kradius, name_species, stringsAsFactors=FALSE)
  d1$label <- strlabels
  d1$kcorelabel = 0
  # Remove empty nodes
  d1 <- d1[d1$label!="EMPTY",]
  
  if (style == "legacy")
    nodelabels <- strlabels
  else if ((style == "kcoreorder")||(style == "chilopod")) {
    for (i in 1:nrow(d1)){
      s <- strsplit(d1$label[i],"shell")
      d1$label[i] <- s[[1]][1]
      d1$kcore[i] <- s[[1]][2]
    }
  }
  d1 <- d1[!(d1$label %in% c("","NA")),]         # Remove empty kshell separator cells 
  d1 <- d1[!is.na(d1$x1),] 
  for (i in 1:nrow(d1)){
    d1[i,]$kdegree <- igraphnet[paste0(strguild,d1[i,]$label)]$kdegree
    d1[i,]$kradius <- igraphnet[paste0(strguild,d1[i,]$label)]$kradius
    d1[i,]$kcorelabel <- igraphnet[paste0(strguild,d1[i,]$label)]$kcorenum
    d1[i,]$name_species <- igraphnet[paste0(strguild,d1[i,]$label)]$name_species
  }

  if (style == "legacy"){
    mlinks <- bpp$result_analysis$matrix > 0
    if (guild=="A")
      degrees <- colSums(mlinks)
    else if (guild=="B")
      degrees <- rowSums(mlinks)
    d1$degree <- 0
    for (i in 1:nrow(d1))
      if(sum(gsub("\\."," ",names(degrees))==d1$name_species[i])>0)
        d1$degree[i]=degrees[which(gsub("\\."," ",names(degrees))==d1$name_species[i])]
    
    ordvector <- rev(order(d1$degree))
    d1$label <- d1[ordvector,]$label
    d1$kradius <- d1[ordvector,]$kradius
    d1$kdegree <- d1[ordvector,]$kdegree
    d1$kcorelabel <- d1[ordvector,]$kcorelabel
    d1$name_species <- d1[ordvector,]$name_species
    d1$degree <- d1[ordvector,]$degree
  } else if ((style=="kcoreorder") || (style=="chilopod")){
    subscol <- which(names(d1)=="kdegree"):ncol(d1)
    shells <- sort(unique(d1$kcore))
    for (k in shells){
      d2 <- d1[d1$kcore == k,]
      if (nrow(d2)>1){
        d2 <-d2[rev(order(1000*d2$kdegree-d2$kradius)),]
        d1[d1$kcore==k,][,subscol] <- d2[,subscol]
      }
      
    }
  }
  # Fixed aspect ratio of bipartite plot
  bpp$landmark_right <- max(bpp$landmark_right,max(d1$x2))
  bpp$tot_width <- max(max(d1$x2)+xstep,bpp$tot_width)
  bpp$tot_height <- (9/16)*bpp$tot_width
  bpp$landmark_top <- max(bpp$landmark_top,max(d1$y2)+xstep)

  return(d1)
}

# Analysis of the chains of specialists

swap_strguild <- function(strguild)
{
  if (strguild == bpp$str_guild_a)
    strguild = bpp$str_guild_b
  else
    strguild = bpp$str_guild_a
  return(strguild)
}

# Draw links among species of same k-index
draw_innercores_tails_bip <- function(p,svg,kc,list_dfs,df_orph,color_guild, inverse="no", is_guild_a = TRUE)
{
  lastx <- 0
  lasty <- 0

  lpoint_x <- 0
  if (length(list_dfs[[kc]])>0)
    if (kc>2)
      lpoint_x <- list_dfs[[kc]][nrow(list_dfs[[kc]]),]$x2
  else
    lpoint_x <- list_dfs[[kc]][nrow(list_dfs[[kc]]),]$x2 + 4*bpp$lado
  if (kc>2)
    lpoint_y <- (list_dfs[[kc]][1,]$y1+list_dfs[[kc]][1,]$y2)/2
  else{
    if (bpp$kcoremax > 2)
      lpoint_y <- (list_dfs[[kc]][1,]$y1+list_dfs[[kc]][1,]$y2)/4
    else{
      lpoint_y <- 1.5*max(list_dfs[[2]]$y2)
    }
  }
  long_tail <- df_orph[(df_orph$kcore == kc) & (df_orph$repeated == "no"),]
  if (length(long_tail)>0){
    v<-  draw_edge_tails_bip(p,svg,lpoint_x,lpoint_y,kc,long_tail,list_dfs,color_guild,
                         inverse = inverse, joinchars = bpp$joinstr,pbackground = "no",
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

# Draw tail connected to highest kdegree node
draw_fat_tail_bip<- function(p,svg,fat_tail,nrows,list_dfs,color_guild,pos_tail_x,pos_tail_y,
                             fgap,
                             inverse="no", is_guild_a =TRUE, bipartite = FALSE, gstyle = "ziggurat")
{

  tailgap <- (1+min(4,nrows/8))*bpp$xstep
  ppos_tail_x <- pos_tail_x-tailgap
  pos_tail_y <-list_dfs[[bpp$kcoremax]][1,]$y1
  ppos_tail_y <- pos_tail_y

  if (nrow(fat_tail)>0)
  {
    nodekcoremax <- list_dfs[[bpp$kcoremax]][1,]
    plyy2 <-  (nodekcoremax$y1+nodekcoremax$y2)/2
    v<- draw_tail_bip(ifelse(is_guild_a, "edge-kcore1-a-fat", "edge-kcore1-b-fat"), p,svg,
                  fat_tail,ifelse(bpp$style=="chilopod",bpp$xstep,bpp$lado),
                  color_guild,gen_sq_label(fat_tail$orph,is_guild_a = is_guild_a, myenv=bpp),
                  ppos_tail_x,ppos_tail_y,fgap,
                  lxx2 = list_dfs[[bpp$kcoremax]][1,]$x1,
                  lyy2 = plyy2,
                  sqinverse = inverse, background = "no", psize = bpp$lsize_kcore1,
                  is_guild_a = is_guild_a, wlink = fat_tail$weightlink[1],
                  style=bpp$style)
    p <- v["p"][[1]]
    svg <- v["svg"][[1]]
  }
  calc_vals <- list("p" = p, "svg" = svg, "pos_tail_x" = bpp$pos_tail_x-tailgap)
  return(calc_vals)
}

find_specialists_bip <- function(specialists_a,specialists_b)
{
  specialists_a <- data.frame(c())
  specialists_b <- data.frame(c())
  if (exists("df_orph_a", envir = bpp))
    if (nrow(bpp$df_orph_a)>0)
    {
      specialists_a <-  bpp$df_orph_a[bpp$df_orph_a$repeated== "yes",]
      specialists_a <-  specialists_a[rev(order(specialists_a$orph,specialists_a$kcore)),]
      if (nrow(specialists_a)>0)
        specialists_a$drawn <- "no"
    }
  if (exists("df_orph_b", envir = bpp))
    if (nrow(bpp$df_orph_b)>0)
    {
      specialists_b <-  bpp$df_orph_b[bpp$df_orph_b$repeated== "yes",]
      specialists_b <-  specialists_b[rev(order(specialists_b$orph,specialists_b$kcore)),]
      if (nrow(specialists_b)>0)
        specialists_b$drawn <- "no"
    }
  
  calc_vals <- list("specialists_a" = specialists_a, "specialists_b" = specialists_b)
  return(calc_vals)
}

# Handle specialist chain species
handle_orphans_bip <- function(vg)
{
  bpp$df_orph_a <- data.frame(c())
  bpp$df_orph_b <- data.frame(c())
  if (length(grep(bpp$str_guild_b,bpp$mtxlinks[1,1]))>0)
    names(bpp$mtxlinks) <- rev(names(bpp$mtxlinks))
  bpp$orphans_a <- bpp$df_cores$species_guild_a[[1]]
  bpp$orphans_b <- bpp$df_cores$species_guild_b[[1]]
  if (!is.null(bpp$orphans_a))
    if (!is.na(bpp$orphans_a[1]))
      bpp$df_orph_a <- find_orphans(bpp$mtxlinks,bpp$orphans_a,bpp$g,guild_a="yes",myenv=bpp)
  if (!is.null(bpp$orphans_b))
    if (!is.na(bpp$orphans_b[1]))
      bpp$df_orph_b <- find_orphans(bpp$mtxlinks,bpp$orphans_b,bpp$g,guild_a="no",myenv=bpp)
  calc_vals <- list("mtxlinks" = bpp$mtxlinks, "orphans_a" = bpp$orphans_a,
                    "orphans_b" = bpp$orphans_b, "df_orph_a" = bpp$df_orph_a, "df_orph_b" = bpp$df_orph_b )
  return(calc_vals)
}

# Draw specialist chains connected to inner ziggurats
draw_bipartite_leafs <- function(p, svg)
{
    if (exists("df_orph_a", envir = bpp)){
      w <- draw_innercores_tails_bip(p,svg,2,bpp$list_dfs_b,bpp$df_orph_a,bpp$color_guild_a, inverse="yes")
      p <- w["p"][[1]]
      svg <- w["svg"][[1]]
      bpp$last_xtail_b[2] <- w["lastx"][[1]]
      bpp$last_ytail_b[2] <-w["lasty"][[1]]
    }

    if (exists("df_orph_b", envir = bpp)){
      w <- draw_innercores_tails_bip(p,svg,2,bpp$list_dfs_a,bpp$df_orph_b,bpp$color_guild_b, inverse="no", is_guild_a = FALSE)
      p <- w["p"][[1]]
      svg <- w["svg"][[1]]
      bpp$last_xtail_a[2] <- w["lastx"][[1]]
      bpp$last_ytail_a[2] <- w["lasty"][[1]]
    }
  #}

  calc_vals <- list("p" = p, "svg" = svg)
  return(calc_vals)
}

draw_bipartite_leafs <- function(p, svg)
{
    if (exists("df_orph_a", envir = bpp)){
      w <- draw_innercores_tails_bip(p,svg,2,bpp$list_dfs_b,bpp$df_orph_a,bpp$color_guild_a, inverse="yes")
      p <- w["p"][[1]]
      svg <- w["svg"][[1]]
      bpp$last_xtail_b[2] <- w["lastx"][[1]]
      bpp$last_ytail_b[2] <-w["lasty"][[1]]
    }
    
    if (exists("df_orph_b", envir = bpp)){
      w <- draw_innercores_tails_bip(p,svg,2,bpp$list_dfs_a,bpp$df_orph_b,bpp$color_guild_b, inverse="no", is_guild_a = FALSE)
      p <- w["p"][[1]]
      svg <- w["svg"][[1]]
      bpp$last_xtail_a[2] <- w["lastx"][[1]]
      bpp$last_ytail_a[2] <- w["lasty"][[1]]
    }
  
  calc_vals <- list("p" = p, "svg" = svg)
  return(calc_vals)
}

# Manage fat tails
handle_fat_tails_bip <- function(p, svg, style = "legacy")
{
  fat_tail_x <- min(bpp$last_xtail_a[[bpp$kcoremax]],
                    bpp$last_xtail_b[[bpp$kcoremax]],
                    bpp$list_dfs_a[[bpp$kcoremax]][1,]$x1,
                    bpp$list_dfs_b[[bpp$kcoremax]][1,]$y2)
  if (bpp$orderkcoremaxby == "kdegree")
    max_b_k <- bpp$list_dfs_b[[bpp$kcoremax]][which(bpp$list_dfs_b[[bpp$kcoremax]]$kdegree == max(bpp$list_dfs_b[[bpp$kcoremax]]$kdegree)),]$label
  if (bpp$orderkcoremaxby == "kradius")
    max_b_k <- bpp$list_dfs_b[[bpp$kcoremax]][which(bpp$list_dfs_b[[bpp$kcoremax]]$kradius == min(bpp$list_dfs_b[[bpp$kcoremax]]$kradius)),]$label

  if (exists("df_orph_a", envir = bpp)){
    fat_tail_a <- bpp$df_orph_a[(bpp$df_orph_a$partner == max(max_b_k[1])) 
                                & (bpp$df_orph_a$repeated == "no"),]
    if (nrow(fat_tail_a)>1)
      bpp$df_orph_a <- bpp$df_orph_a[!(bpp$df_orph_a$orph %in% fat_tail_a$orph),]
    tailweight <- 0
    if (nrow(fat_tail_a)>0) {
    for (h in 1:nrow(fat_tail_a))
      tailweight <- tailweight + bpp$result_analysis$matrix[as.numeric(fat_tail_a$partner[h]),
                                                                  as.numeric(fat_tail_a$orph[h])]
      fat_tail_a$weightlink <- get_link_weights(tailweight, myenv=bpp)
    }
  }
  if (!exists("fat_tail_a"))
    fat_tail_a <- data.frame(c())
  if (bpp$orderkcoremaxby == "kdegree")
    max_a_k <- bpp$list_dfs_a[[bpp$kcoremax]][which(bpp$list_dfs_a[[bpp$kcoremax]]$kdegree == max(bpp$list_dfs_a[[bpp$kcoremax]]$kdegree)),]$label
  if (bpp$orderkcoremaxby == "kradius")
    max_a_k <- bpp$list_dfs_a[[bpp$kcoremax]][which(bpp$list_dfs_a[[bpp$kcoremax]]$kradius == min(bpp$list_dfs_a[[bpp$kcoremax]]$kradius)),]$label
  
  if (exists("df_orph_b", envir = bpp)){
    fat_tail_b <- bpp$df_orph_b[(bpp$df_orph_b$partner == max(max_a_k)) & (bpp$df_orph_b$repeated == "no"),]
    if (nrow(fat_tail_b)>1)
      bpp$df_orph_b <- bpp$df_orph_b[!(bpp$df_orph_b$orph %in% fat_tail_b$orph),]
    
    tailweight <- 0
    if (nrow(fat_tail_b)>0) {
      for (h in 1:nrow(fat_tail_b))
        tailweight <- tailweight+bpp$result_analysis$matrix[as.numeric(fat_tail_b$orph[h]),
                                                    as.numeric(fat_tail_b$partner[h])]
      fat_tail_b$weightlink <- get_link_weights(tailweight, myenv=bpp)
    }
  }
  if (!exists("fat_tail_b"))
    fat_tail_b <- data.frame(c())
  nrows_fat <- nrow(fat_tail_b)+nrow(fat_tail_a)
  if (style=="chilopod")
    fgap <- 0.7*bpp$hop_x
  else
    fgap <- 0.7*bpp$hop_x + (1+sum(nrows_fat>40))*bpp$lado
  bpp$pos_tail_x <- min(bpp$list_dfs_b[[bpp$kcoremax]][1,]$x1-bpp$xstep,
                        bpp$list_dfs_a[[bpp$kcoremax]][1,]$x1-bpp$xstep)
  if (exists("fat_tail_a")) {
    f <- draw_fat_tail_bip(p,svg,fat_tail_a,nrows_fat,bpp$list_dfs_b,bpp$color_guild_a[2],
                       bpp$pos_tail_x,pos_tail_y,fgap,inverse="yes")
    p <- f["p"][[1]]
    svg <- f["svg"][[1]]
  }

  if (exists("fat_tail_b")) {
    f <- draw_fat_tail_bip(p,svg,fat_tail_b,nrows_fat,bpp$list_dfs_a,bpp$color_guild_b[2],bpp$pos_tail_x,
                       pos_tail_y,fgap,
                       inverse="no", is_guild_a = FALSE)
    
    p <- f["p"][[1]]
    svg <- f["svg"][[1]]
  }
  exists_fat_tail <- (nrows_fat > 0)
  bpp$landmark_left <- min(bpp$landmark_left,bpp$pos_tail_x-2*bpp$step)
  calc_vals <- list("p" = p, "svg" = svg, "pos_tail_x" = bpp$pos_tail_x, "exists_fat_tail" = exists_fat_tail)
  return(calc_vals)
}

draw_maxcore_bip <- function(svg)
{
  kcoremax_label_display <- function (idPrefix, gp,svg,kcoremaxlabel_angle,pdata,plabel,plabelsize,phjust=0, is_guild_a = TRUE) 
  {

    labelcolor <- ifelse(length(bpp$labels_color)>0,bpp$labels_color[2-as.numeric(is_guild_a)], pdata$col_row)
    if (kcoremaxlabel_angle == 0) 
    {
        gp <- gp +  geom_text(data=pdata, aes(x=x1+(x2-x1)/2, y=y1+(y2-y1)/2), label=plabel,
                              color = labelcolor, size=plabelsize, angle = kcoremaxlabel_angle)
        svg$text(idPrefix=idPrefix, data=pdata, mapping=aes(x=x1+(x2-x1)/2, y=y1+(y2-y1)/2), 
                 label=plabel, color=labelcolor, size=fsvgtext*plabelsize, angle=kcoremaxlabel_angle)
    } 
    else {
        gp <- gp + geom_text(data=pdata, aes(x=x1, y=y1+(y2-y1)/20), label=plabel,
                             color = labelcolor, size=plabelsize, angle = kcoremaxlabel_angle,
                             vjust = 1, hjust = phjust)
        svg$text(idPrefix=idPrefix, data=pdata, mapping=aes(x=x1, y=y1+(y2-y1)/20), label=plabel, color=labelcolor, 
                 size=fsvgtext*plabelsize, angle=kcoremaxlabel_angle)
            }
    calc_vals <- list("p" = gp, "svg" = svg)
    return(calc_vals)
  }
  
  adjust_lengths <- function(spA,spB,fill="right"){
    diff_lengths <- length(spA) - length(spB)
    if (fill=="right"){
      specA <- c(rep("EMPTY",abs(diff_lengths)*(diff_lengths<0)),spA)
      specB <- c(rep("EMPTY",diff_lengths*(diff_lengths>0)),spB)
    }
    else{
      specA <- c(spA,rep("EMPTY",abs(diff_lengths)*(diff_lengths<0)))
      specB <- c(spB,rep("EMPTY",diff_lengths*(diff_lengths>0)))
    }
    calc_vals <- list("specA" = specA, "specB" = specB)
    return(calc_vals)
  }
  
  paint_rect_core <- function(list_dfs,alpha_level){
    q <- geom_rect(data=list_dfs, 
                   mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), 
                   fill = list_dfs$col_row,  
                   color="transparent",alpha=alpha_level)
    return(q)
  }

  paint_rect_svg <- function(guildstr,list_dfs){
    svg$rect(paste0("kcore", bpp$kcoremax, guildstr), data=list_dfs,
             mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), 
             fill=list_dfs$col_row, alpha=bpp$alpha_level,
             color="transparent", size=0.5)
  }
  
  paint_labels <- function(p,svg,guildstr,list_dfs){

    nsp <- name_species_preprocess(bpp$kcoremax,list_dfs,bpp$kcore_species_name_display,
                                   bpp$kcore_species_name_break)
    labelszig <- nsp$labelszig
    kcoremaxlabel_angle <- nsp$kcoremaxlabel_angle
    p <- p + paint_rect_core(list_dfs,alpha=bpp$alpha_level)
    paint_rect_svg(guildstr,list_dfs)
    f <- kcoremax_label_display(paste0("kcore", bpp$kcoremax, guildstr),p,svg,kcoremaxlabel_angle,
                                list_dfs,labelszig,
                                bpp$lsize_kcoremax, phjust = 1, is_guild_a = guildstr=="-a")
    return(f)
  }
  
  #Species outside the giant component
  outsiders_A <- gsub(bpp$str_guild_a,"",bpp$outsiders_a)
  outsiders_B <- gsub(bpp$str_guild_b,"",bpp$outsiders_b)
  if (bpp$style == "legacy"){
    species_A <- unlist(bpp$df_cores$species_guild_a)
    species_B <- unlist(bpp$df_cores$species_guild_b)
    species_A <- c(species_A,outsiders_A,"EMPTY")
    species_B <- c(species_B,outsiders_B,"EMPTY")
  }
  else if ((bpp$style == "kcoreorder") ||(bpp$style == "chilopod")) {
    species_A <- c()
    species_B <- c()
    if (bpp$style == "kcoreorder"){
      # Kcore-1 on the left side of the plot
      # shell number is passed together with node label
      spA1 <- paste0(bpp$df_cores$species_guild_a[[1]],"shell1")
      spB1 <- paste0(bpp$df_cores$species_guild_b[[1]],"shell1")
      sadj <- adjust_lengths(spA1,spB1)
      species_A <- sadj$specA
      species_B <- sadj$specB
    }
    if (bpp$style == "chilopod"){
      spA1_spec <- c()
      spB1_spec <- c()
      # specialist chains
      sbip <- find_specialists_bip(specialists_a,specialists_b)
      if (nrow(sbip$specialists_a)>0)
        spA1_spec <- paste0(unique(sbip$specialists_a$orph),"shell1")
      if (nrow(sbip$specialists_b)>0)
        spB1_spec <- paste0(unique(sbip$specialists_b$orph),"shell1")
      # shell number is passed together with node label
      sadj <- adjust_lengths(spA1_spec,spB1_spec,fill="right")
      spA1_spec <- sadj$specA
      spB1_spec <- sadj$specB
    }
    for (k in bpp$kcoremax:2) {
      sadj <- adjust_lengths(paste0(rev(unlist(bpp$df_cores[k,]$species_guild_a)),"shell",k),
                             paste0(rev(unlist(bpp$df_cores[k,]$species_guild_b)),"shell",k),
                             fill = "left")  
      species_A <- c(species_A, "EMPTY", sadj$specA) 
      species_B <- c(species_B, "EMPTY", sadj$specB) 
    }
    # Specialist chains 
    if (bpp$style == "chilopod"){
      species_A <- c(species_A, "EMPTY", spA1_spec) 
      species_B <- c(species_B, "EMPTY", spB1_spec)
    }  
    if (length(outsiders_A>0)){
      sadj <- adjust_lengths(outsiders_A,outsiders_B)
      species_A <- c(species_A, "EMPTY", paste0(sadj$specA,"shell",0))
      species_B <- c(species_B, "EMPTY", paste0(sadj$specB,"shell",0))
      species_A[grepl("EMPTY",species_A)] <- "EMPTY"
      species_B[grepl("EMPTY",species_B)] <- "EMPTY"
    }
  }
  nnodes <- max(length(species_A), length(species_B))
  bpp$list_dfs_a[[bpp$kcoremax]]<- draw_parallel_guilds(bpp$basex,bpp$topxa,bpp$basey*bpp$guild_gap_increase,bpp$toopy,
                                                         length(species_A),nnodes,bpp$color_guild_a,
                                                         species_A,
                                                         bpp$rg, bpp$str_guild_a, 
                                                         orderby = "kdegree",
                                                         style=bpp$style,guild="A")
  p <- ggplot() + scale_x_continuous(name="x") + scale_y_continuous(name="y")
  f <- paint_labels(p,svg,"-a",bpp$list_dfs_a[[bpp$kcoremax]])
  p <- f["p"][[1]]
  svg <- f["svg"][[1]]
  num_b_coremax <- bpp$df_cores[bpp$kcoremax,]$num_species_guild_b
  bpp$basey <- - bpp$basey
  bpp$topxb <- bpp$topxa
  bpp$toopy <- - bpp$toopy
  bpp$list_dfs_b[[bpp$kcoremax]] <- draw_parallel_guilds(bpp$basex,bpp$topxb,bpp$basey*bpp$guild_gap_increase,bpp$toopy,
                                                          length(species_B),nnodes,
                                                          bpp$color_guild_b,
                                                          species_B,bpp$rg,
                                                          bpp$str_guild_b,  orderby = "kdegree",
                                                          style=bpp$style,guild="B")
    bpp$landmark_right <- max(bpp$list_dfs_b[[bpp$kcoremax]]$x2,
                            bpp$list_dfs_a[[bpp$kcoremax]]$x2)+bpp$xstep
  f <- paint_labels(p,svg,"-b",bpp$list_dfs_b[[bpp$kcoremax]])
  calc_vals <- list("p" = f["p"][[1]], "svg" = f["svg"][[1]], "basey" = bpp$basey, 
                    "topy" = bpp$toopy, "topxa" = bpp$topxa, "topxb" = bpp$topxb,
                    "list_dfs_a" = bpp$list_dfs_a, "list_dfs_b" = bpp$list_dfs_b,
                    "last_xtail_a" = bpp$last_xtail_a, "last_ytail_a" = bpp$last_ytail_a,
                    "last_xtail_b" = bpp$last_xtail_b, "last_ytail_b" = bpp$last_ytail_b)
  return(calc_vals)
}

draw_maxcore_tails_bip <- function(p, svg)
{
  long_kcoremax_tail <- FALSE
  #leftjump <- (1.2-0.02*nrow(bpp$list_dfs_a[[bpp$kcoremax]])) * bpp$hop_x
  point_x <- bpp$list_dfs_a[[bpp$kcoremax]][nrow(bpp$list_dfs_a[[bpp$kcoremax]]),]$x2 #- leftjump
  point_y <- bpp$xstep
  if (exists("df_orph_a", envir = bpp))
    long_tail_a <- bpp$df_orph_a[(bpp$df_orph_a$repeated == "no"),]
  if ((exists("long_tail_a")) & (bpp$kcoremax > 2))
  {
    if (length(long_tail_a)>5)
      long_kcoremax_tail <- TRUE
    v<-  draw_edge_tails_bip(p,svg,point_x,point_y*bpp$aspect_ratio,bpp$kcoremax,
                         long_tail_a,bpp$list_dfs_b,bpp$color_guild_a, inverse = "yes",
                         vertical = "no", orientation = "South", revanddrop = "yes",
                         pbackground = "no", tspline = "arc", joinchars=bpp$joinstr)
    p <- v["p"][[1]]
    svg <- v["svg"][[1]]
    bpp$last_xtail_b[bpp$kcoremax] <- v["lastx"][[1]]
    if (length(v["lasty"][[1]])>0)
      bpp$last_ytail_b[bpp$kcoremax] <- v["lasty"][[1]]
    else                                                   # only for degenerate networks with all nodes in kcoremax
      bpp$last_ytail_b[bpp$kcoremax] <- bpp$toopy
  }
  leftjump <- (1.2-0.02*nrow(bpp$list_dfs_b[[bpp$kcoremax]]))* bpp$hop_x
  point_x <- bpp$list_dfs_b[[bpp$kcoremax]][nrow(bpp$list_dfs_b[[bpp$kcoremax]]),]$x2
  point_y <- -bpp$xstep
  if (exists("df_orph_b", envir = bpp))
    long_tail_b <- bpp$df_orph_b[(bpp$df_orph_b$repeated == "no"),]
  if ( (exists("long_tail_b")) & (bpp$kcoremax > 2) ){
    if (nrow(long_tail_b)>5)
      long_kcoremax_tail <- TRUE
    tailweight <- 0
    if (nrow(long_tail_b)>0){
      for (h in 1:nrow(long_tail_b))
        tailweight <- tailweight + bpp$result_analysis$matrix[as.numeric(long_tail_b$orph[h]),
                                                                    as.numeric(long_tail_b$partner[h])]
      long_tail_b$weightlink <- get_link_weights(tailweight, myenv=bpp)
    }

    v <-  draw_edge_tails_bip(p,svg,point_x,point_y*bpp$aspect_ratio,
                              bpp$kcoremax,long_tail_b,bpp$list_dfs_a,bpp$color_guild_b,
                         inverse = "no",
                         vertical = "no", orientation = "North", revanddrop = "yes",
                         pbackground = "no", tspline = "arc", joinchars=bpp$joinstr, is_guild_a = FALSE,
                         wlink = long_tail_b$weightlink[1])
    p <- v["p"][[1]]
    svg <- v["svg"][[1]]
    bpp$last_xtail_a[bpp$kcoremax] <- v["lastx"][[1]]
    if (length(v["lasty"][[1]])>0)
      bpp$last_ytail_b[bpp$kcoremax] <- v["lasty"][[1]]
    else                                                   # only for degenerate networks with all nodes in kcoremax
      bpp$last_ytail_b[bpp$kcoremax] <- bpp$toopy
  }
  calc_vals <- list("p" = p, "svg" = svg,
                    "last_xtail_a" = bpp$last_xtail_a, "last_ytail_a" = bpp$last_ytail_a,
                    "last_xtail_b" = bpp$last_xtail_b, "last_ytail_b" = bpp$last_ytail_b)
  return(calc_vals)
}

display_plot_bip <- function(p, printfile,  plwidth=14, ppi = 300, landscape = bpp$label_strguild, fname_append = "")
{
  if (printfile){
    if (fname_append != "")
      ftname_append <- paste0("_",fname_append)
    else
      ftname_append <- fname_append
    dir.create(bpp$plotsdir, showWarnings = FALSE)
    if (landscape)
      png(paste0(bpp$plotsdir,"/",bpp$network_name,"_ziggurat",ftname_append,".png"), width=(plwidth*ppi), height=(9/16)*plwidth*ppi, res=ppi)
    else
      png(paste0(bpp$plotsdir,"/",bpp$network_name,"_ziggurat",ftname_append,".png"), width=(plheight*ppi), height=(9/16)*plwidth*ppi, res=ppi)
  }
  print(p)
  if (printfile)
    dev.off()
}

def_configuration_bip <- function(paintlinks, print_to_file, plotsdir, orderkcoremaxby, style, 
                              guild_gap_increase, flip_results, aspect_ratio,
                              alpha_level, color_guild_a, color_guild_b,
                              color_link, alpha_link, size_link,
                              lsize_kcoremax, lsize_zig, lsize_kcore1,
                              lsize_legend, lsize_core_box, labels_color,
                              factor_hop_x,
                              paint_outsiders, displace_outside_component,
                              outsiders_separation_expand, outsiders_legend_expand, specialistskcore2_horizontal_dist_rootleaf_expand,
                              specialistskcore2_vertical_dist_rootleaf_expand , specialists_boxes_separation_count,
                              root_specialist_expand,hide_plot_border,rescale_plot_area,kcore1specialists_leafs_vertical_separation,
                              corebox_border_size, kcore_species_name_display,kcore_species_name_break,
                              shorten_species_name,exclude_species_number,
                              label_strguilda, label_strguildb, landscape_plot, backg_color, show_title, show_legend,
                              use_spline, spline_points, file_name_append, svg_scale_factor, weighted_links, square_nodes_size_scale,
                              move_all_SVG_up, move_all_SVG_right, progress
                              )
{
  # ENVIRONMENT CONFIGURATION PARAMETERS
  bpp$paintlinks <- paintlinks
  bpp$print_to_file <- print_to_file
  bpp$plotsdir <- plotsdir
  bpp$orderkcoremaxby <- orderkcoremaxby
  bpp$style <- style
  bpp$guild_gap_increase <- guild_gap_increase
  bpp$flip_results <- flip_results
  bpp$alpha_level <- alpha_level
  bpp$color_guild_a <- color_guild_a
  bpp$color_guild_b <- color_guild_b
  bpp$color_link <- color_link
  bpp$alpha_link <- alpha_link
  bpp$size_link <- size_link
  bpp$aspect_ratio <- aspect_ratio
  bpp$labels_size <- 4
  bpp$lsize_kcoremax <- lsize_kcoremax
  bpp$lsize_zig <- lsize_zig
  bpp$lsize_kcore1 <- lsize_kcore1
  bpp$lsize_legend <- lsize_legend
  bpp$lsize_core_box <- lsize_core_box
  bpp$labels_color <- labels_color                           # Horizontal & Vertical distances of edge/specialist tails linked to core 1 North & South
  bpp$factor_hop_x <- factor_hop_x
  bpp$coremax_triangle_height_factor <- 3
  bpp$coremax_triangle_width_factor <- 3
  bpp$paint_outsiders <- paint_outsiders
  bpp$displace_outside_component <- displace_outside_component
  bpp$outsiders_separation_expand <- outsiders_separation_expand
  bpp$outsiders_legend_expand <- outsiders_legend_expand
  bpp$specialistskcore2_horizontal_dist_rootleaf_expand <- specialistskcore2_horizontal_dist_rootleaf_expand        # Controls the distance of specialist root leaves to partner in core 2
  bpp$specialistskcore2_vertical_dist_rootleaf_expand <- specialistskcore2_vertical_dist_rootleaf_expand
  bpp$specialists_boxes_separation_count <- specialists_boxes_separation_count                  # Separation of leaves of a specialist tail
  bpp$root_specialist_expand <- root_specialist_expand
  bpp$hide_plot_border <- hide_plot_border
  bpp$rescale_plot_area <- rescale_plot_area
  bpp$kcore1specialists_leafs_vertical_separation <- kcore1specialists_leafs_vertical_separation
  bpp$corebox_border_size <- corebox_border_size
  bpp$kcore_species_name_display <- kcore_species_name_display
  bpp$kcore_species_name_break <- kcore_species_name_break
  bpp$shorten_species_name <- shorten_species_name
  bpp$exclude_species_number <- exclude_species_number
  bpp$landscape_plot <- landscape_plot
  bpp$backg_color <- backg_color
  bpp$show_title <- show_title
  bpp$show_legend <- show_legend
  bpp$use_spline <- use_spline
  bpp$spline_points <- spline_points
  bpp$file_name_append <- file_name_append
  bpp$svg_scale_factor <- svg_scale_factor
  bpp$weighted_links <- weighted_links
  bpp$square_nodes_size_scale <- square_nodes_size_scale
  bpp$move_all_SVG_up <- move_all_SVG_up
  bpp$move_all_SVG_right <- move_all_SVG_right
  bpp$progress <- progress
}

draw_inner_links_bip <- function(p, svg)
{
  for (kcb in seq(bpp$kcoremax,2))
  {
    for (kc in seq(bpp$kcoremax,2))
    {
      labels_a <- bpp$list_dfs_a[[kc]]$label
      for (j in seq(along=labels_a))
      {
        numberlinksa <- sum(bpp$mtxlinks$guild_a == paste0(bpp$str_guild_a,labels_a[j]) )
        data_a <- bpp$list_dfs_a[[kc]][j,]
        foundlinksa <- 0
        labels_b <- bpp$list_dfs_b[[kcb]]$label
        for (i in seq(along =labels_b))
        {
          if (sum(bpp$mtxlinks$guild_a == paste0(bpp$str_guild_a,labels_a[j]) & bpp$mtxlinks$guild_b == paste0(bpp$str_guild_b,labels_b[i]))>0)
          {
            foundlinksa <- foundlinksa + 1
            data_b <- bpp$list_dfs_b[[kcb]][i,]
            weightlink <- get_link_weights(bpp$result_analysis$matrix[as.numeric(data_b$label),
                                                                      as.numeric(data_a$label)],myenv=bpp)
            bend_line = "no"
            if (((kc == 2) & (kcb == bpp$kcoremax)) | ((kc == bpp$kcoremax) & (kcb == 2)))
              bend_line = "horizontal"
            if ((kc == bpp$kcoremax) & (kcb == bpp$kcoremax))
            {
              link <- data.frame(x1= data_a$x1 + (data_a$x2-data_a$x1)/2,
                                 x2 = data_b$x1 +(data_b$x2-data_b$x1)/2,
                                 y1 = data_a$y1,  y2 = bpp$list_dfs_b[[kcb]]$y1[i] )
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
              if (kc == bpp$kcoremax)
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
              if (kcb == bpp$kcoremax){
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
                     slink = bpp$size_link*weightlink, clink =  c(bpp$color_link),
                     alpha_l = bpp$alpha_link , spline = bend_line, myenv=bpp)
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


draw_bipartite_plot <- function(svg_scale_factor, progress)
{
  if (!is.null(progress)) 
    progress$inc(1/11, detail=strings$value("MESSAGE_ZIGGURAT_PROGRESS_PROCESSING_NODES"))
  zinit_time <- proc.time()
  for (i in bpp$ind_cores) {
    nodes_in_core_a <- bpp$g[(bpp$g$guild == bpp$str_guild_a)&(bpp$g$kcorenum == i)]$name
    nodes_in_core_b <- bpp$g[(bpp$g$guild == bpp$str_guild_b)&(bpp$g$kcorenum == i)]$name
    bpp$df_cores[i,]$species_guild_a <- list(unlist(lapply(nodes_in_core_a, function(x) strsplit(x,bpp$str_guild_a)[[1]][[2]])))
    bpp$df_cores[i,]$num_species_guild_a <- length(nodes_in_core_a)
    bpp$df_cores[i,]$species_guild_b <- list(unlist(lapply(nodes_in_core_b, function(x) strsplit(x,bpp$str_guild_b)[[1]][[2]])))
    bpp$df_cores[i,]$num_species_guild_b <- length(nodes_in_core_b)
  }
  bpp$num_a_coremax <- bpp$df_cores[bpp$kcoremax,]$num_species_guild_a
  base_width <- 2000
  bpp$ymax <- 1.75*base_width/bpp$aspect_ratio
  bpp$tot_width <- bpp$ymax

  bpp$species_in_core2_a <- sum(bpp$df_cores[2,]$num_species_guild_a)
  bpp$species_in_core2_b <- sum(bpp$df_cores[2,]$num_species_guild_b)
  maxincore2 <- max(bpp$species_in_core2_a,bpp$species_in_core2_b)
  bpp$height_y <- bpp$height_y
  bpp$ymax <- bpp$ymax * 1.1
  bpp$hop_x <- bpp$factor_hop_x*(bpp$tot_width)/max(1,(bpp$kcoremax-2))
  bpp$lado <- min(0.05*bpp$tot_width,bpp$height_y * bpp$aspect_ratio)
  bpp$basey <- 2.5*bpp$lado
  wcormax <- 1.2*bpp$hop_x*bpp$coremax_triangle_width_factor
  bpp$topxa <- 0.65*bpp$hop_x
  bpp$basex <- bpp$topxa - wcormax
  bpp$posic_zig <- 0
  bpp$posic_zig_a <- 0
  bpp$posic_zig_b <- 0
  bpp$toopy <- 0.3*bpp$ymax+bpp$basey
  bpp$exists_fat_tail <- FALSE
  bpp$landmark_left <- 0
  bpp$landmark_right <- 0
  
  # Draw max core 
  svg <-SVG(svg_scale_factor, style = bpp$style, 
            nnodes=bpp$result_analysis$num_guild_a+
            bpp$result_analysis$num_guild_b,
            flip_coordinates=bpp$flip_results)
  
  f <- handle_orphans_bip(bpp$result_analysis$graph)
  bpp$mtxlinks <- f["mtxlinks"][[1]]
  bpp$orphans_a <- f["orphans_a"][[1]]
  bpp$orphans_b <- f["orphans_b"][[1]]
  bpp$df_orph_a <- f["df_orph_a"][[1]]
  bpp$df_orph_b <- f["df_orph_b"][[1]]
  
  sbip <- find_specialists_bip(specialists_a,specialists_b)
  
  if (!is.null(progress))
    progress$inc(1/11, detail=strings$value("MESSAGE_ZIGGURAT_PROGRESS_DRAWING_MAXCORE"))
  f <- draw_maxcore_bip(svg)
  p <- f["p"][[1]]
  svg <- f["svg"][[1]]
  if ((bpp$style=="legacy") || (bpp$style=="kcoreorder"))
    bpp$pos_tail_x <- min(bpp$list_dfs_a[[bpp$kcoremax]]$x1)  # Leftmost side of the leftmost node
  
  if (bpp$paintlinks) {
    z <- draw_inner_links_bip(p, svg)
    p <- z["p"][[1]]
    svg <- z["svg"][[1]]
  }

  if (bpp$style=="chilopod"){
    bpp$posic_zig <- f["posic_zig"][[1]] 
    bpp$list_dfs_a <- f["list_dfs_a"][[1]]
    bpp$list_dfs_b <- f["list_dfs_b"][[1]] 
    bpp$last_xtail_a <-  f["last_xtail_a"][[1]] 
    bpp$last_ytail_a <- f["last_ytail_a"][[1]]
    bpp$last_xtail_b <- f["last_xtail_b"][[1]] 
    bpp$last_ytail_b <-  f["last_ytail_b"][[1]]
    bpp$topy <- f["topy"][[1]] 
    bpp$topxa <-f["topxa"][[1]] 
    bpp$topxb <- f["topxb"][[1]] 
    # bpp$posic_zig <- f["posic_zig"][[1]] 
    # bpp$list_dfs_a <- f["list_dfs_a"][[1]]
    # bpp$list_dfs_b <- f["list_dfs_b"][[1]] 
    # bpp$last_xtail_a <-f["last_xtail_a"][[1]] 
    # bpp$last_ytail_a <- f["last_ytail_a"][[1]]
    # bpp$last_xtail_b <- f["last_xtail_b"][[1]] 
    # bpp$last_ytail_b <- f["last_ytail_b"][[1]] 
    
    if (!is.null(bpp$df_cores[1,])){
      if (bpp$df_cores[1,]$num_species_guild_a > 0)
        bpp$list_dfs_a[[1]] <- conf_kcore1_info(bpp$str_guild_a,myenv=bpp)
      if (bpp$df_cores[1,]$num_species_guild_b > 0)
        bpp$list_dfs_b[[1]] <- conf_kcore1_info(bpp$str_guild_b,myenv=bpp)
    }
    
    # Fat tails - nodes of core 1  linked to most generalist of opposite guild. Left side of panel 
    z <-  handle_fat_tails_bip(p, svg,style = bpp$style ) 
    p <- z["p"][[1]] 
    svg <- z["svg"][[1]]
    bpp$pos_tail_x <- z["pos_tail_x"][[1]] 
    bpp$exists_fat_tail <- ifelse(!is.null(z["exists_fat_tail"][[1]]),z["exists_fat_tail"][[1]],FALSE)
    # Hanlde orphans, species outside the ziggurat
    if (!is.null(progress))
      progress$inc(1/11,  detail=strings$value("MESSAGE_ZIGGURAT_PROGRESS_DRAWING_ORPHANS"))
    bpp$gap <- 4*bpp$height_y 
    if (!is.null(progress)) 
      progress$inc(1/11,detail=strings$value("MESSAGE_ZIGGURAT_PROGRESS_DRAWING_MAXCORE_TAILS")) 
    f <- draw_maxcore_tails_bip(p, svg) 
    p <- f["p"][[1]] 
    svg <- f["svg"][[1]]
    
    bpp$last_xtail_a <- f["last_xtail_a"][[1]] 
    bpp$last_ytail_a <-f["last_ytail_a"][[1]] 
    bpp$last_xtail_b <- f["last_xtail_b"][[1]]
    bpp$last_ytail_b <- f["last_ytail_b"][[1]] 
    # Nodes of core 1 linked to species in cores kcoremax-1 to core 2. 
    if (!is.null(progress)) 
      progress$inc(1/11,detail=strings$value("MESSAGE_ZIGGURAT_PROGRESS_DRAWING_INNER_ORPHANS")) 
    z <- draw_bipartite_leafs(p, svg) 
    p <- z["p"][[1]] 
    svg <- z["svg"][[1]]
  }
  #v <- write_annotations_bip(p, svg)
  v <- write_final_annotations(p, svg, bpp$style, myenv=bpp)
  p <- v["p"][[1]]
  svg <- v["svg"][[1]]
  bpp$landmark_left <<- v["landmark_left"][[1]]
  bpp$landmark_right <<- v["landmark_right"][[1]]
  bpp$landmark_top <<- v["landmark_top"][[1]]
  bpp$landmark_bottom <- v["landmark_bottom"][[1]]
  # Plot straight links
  if (!is.null(progress))
    progress$inc(1/11, detail=strings$value("MESSAGE_ZIGGURAT_PROGRESS_DRAWING_LINKS"))
  if (bpp$paintlinks){
    if (nrow(bpp$straight_links)>0) {
      p <- p+ geom_segment(data=bpp$straight_links, aes(x=x1, y=y1, xend=x2, yend=y2),
                           linewidth=bpp$straight_links$weightlink, color=bpp$color_link ,alpha=bpp$alpha_link)
      factormult <- 0.1*svg_scale_factor
      svg$segment(idPrefix="link", data=bpp$straight_links, mapping=aes(x=x1, y=y1, xend=x2, yend=y2),
                  alpha=bpp$alpha_link, color=bpp$color_link, size=bpp$straight_links$weightlink)
    }
    if (nrow(bpp$bent_links)>0) {
      p <- p + geom_path(data =bpp$bent_links,aes(x,y,group=number), linewidth=bpp$bent_links$weightlink,
                         color=bpp$color_link ,alpha=bpp$alpha_link)
      svg$path(idPrefix="link", data=bpp$bent_links, mapping=aes(x, y, group=number), alpha=bpp$alpha_link,
               color=bpp$color_link,size=bpp$bent_links$weightlink)
    }
  }
  if (bpp$flip_results)
    p <- p+coord_flip()+scale_x_reverse()
  
  if (is.null(progress))
     display_plot_bip(p,bpp$print_to_file, landscape = bpp$landscape_plot, fname_append = bpp$file_name_append)
  
  # Stores results
  bpp$plot  <- p
  bpp$svg   <- svg
  # 
  # html<-svg$html()
  # cat(html, file = "tmp.svg")
    return(bpp)

}
if (debugging)
  bipartite_graph("../data/","RA_HP_042.csv",square_nodes_size_scale = 2,show_title = TRUE,
                       style="chilopod",orderkcoremaxby = "kdegree",
                       guild_gap_increase = 1,weighted_links = "none",
                       svg_scale_factor = 1,color_link = "#6d6d6e",
                       hide_plot_border = TRUE)