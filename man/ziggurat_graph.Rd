% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ziggurat_graph.R
\name{ziggurat_graph}
\alias{ziggurat_graph}
\title{Plotting a ziggurat graph}
\usage{
ziggurat_graph(
  datadir,
  filename,
  style = "ziggurat",
  paintlinks = TRUE,
  print_to_file = FALSE,
  plotsdir = "plotresults/",
  orderkcoremaxby = "kradius",
  flip_results = FALSE,
  aspect_ratio = 1,
  alpha_level = 0.2,
  color_guild_a = c("#4169E1", "#00008B"),
  color_guild_b = c("#F08080", "#FF0000"),
  color_link = "slategray3",
  alpha_link = 0.5,
  size_link = 0.5,
  displace_y_b = rep(0, 20),
  displace_y_a = rep(0, 20),
  lsize_kcoremax = 3.5,
  lsize_zig = 3,
  lsize_kcore1 = 2.5,
  lsize_legend = 4,
  lsize_core_box = 2.5,
  labels_color = c(),
  height_box_y_expand = 1,
  kcore2tail_vertical_separation = 1,
  kcore1tail_disttocore = c(1, 1),
  innertail_vertical_separation = 1,
  factor_hop_x = 1,
  fattailjumphoriz = c(1, 1),
  fattailjumpvert = c(1, 1),
  coremax_triangle_height_factor = 1,
  coremax_triangle_width_factor = 1,
  paint_outsiders = TRUE,
  displace_outside_component = c(0, 0),
  outsiders_separation_expand = 1,
  outsiders_legend_expand = 1,
  specialistskcore2_horizontal_dist_rootleaf_expand = 1,
  specialistskcore2_vertical_dist_rootleaf_expand = 0,
  specialists_boxes_separation_count = 1,
  root_specialist_expand = c(1, 1),
  hide_plot_border = TRUE,
  rescale_plot_area = c(1, 1),
  kcore1specialists_leafs_vertical_separation = 1,
  corebox_border_size = 0.2,
  kcore_species_name_display = c(),
  kcore_species_name_break = c(),
  shorten_species_name = 0,
  exclude_species_number = FALSE,
  label_strguilda = "",
  label_strguildb = "",
  landscape_plot = TRUE,
  backg_color = "white",
  show_title = TRUE,
  show_legend = "BOTTOM",
  use_spline = TRUE,
  spline_points = 10,
  file_name_append = "",
  svg_scale_factor = 10,
  weighted_links = "none",
  square_nodes_size_scale = 1,
  move_all_SVG_up = 0,
  move_all_SVG_right = 0,
  progress = NULL
)
}
\arguments{
\item{datadir}{the name of the file of the interaction matrix}

\item{filename}{the file with the interaction matrix}

\item{style}{equals 'ziggurat', do not modify}

\item{print_to_file}{if set to FALSE the plot is displayed in the R session window}

\item{plotsdir}{the directory where the plot is stored}

\item{orderkcoremaxby}{sets order of kcoremax nodes, by kradius or kdegree}

\item{flip_results}{displays the graph in portrait configuration}

\item{aspect_ratio}{ziggurat plot default aspect ratio}

\item{alpha_level}{transparency for ziggurats' filling}

\item{color_guild_a}{default filling for nodes of guild_a}

\item{color_guild_b}{default filling for nodes of guild_b}

\item{color_link}{default links color}

\item{alpha_link}{link transparency}

\item{size_link}{width of the links}

\item{displace_y_b}{relative vertical displacement of guild_b inner ziggurats}

\item{displace_y_a}{relative vertical displacement of guild_a inner ziggurats}

\item{lsize_kcoremax}{nodes in kshell max label size}

\item{lsize_zig}{nodes in inner ziggurats label size}

\item{lsize_kcore1}{labels of nodes in kshell 1}

\item{lsize_legend}{legend label size}

\item{labels_color}{default label colors}

\item{height_box_y_expand}{expand inner ziggurat rectangles default height by this factor}

\item{kcore2tail_vertical_separation}{expand vertical of kshell 1 species linked to kshell 2 by this factor}

\item{kcore1tail_disttocore}{expand vertical separation of kshell 1 species from kshell max (guild_a, guild,b)}

\item{innertail_vertical_separation}{expand vertical separation of kshell species connected to khsell > 2 & < kshell max}

\item{factor_hop_x}{expand inner ziggurats horizontal distance}

\item{fattailjumphoriz}{displace kshell 1 species linked to leftmost kshell max species}

\item{fattailjumpvert}{idem for vertical position}

\item{coremax_triangle_height_factor}{expand khsell max rectangles height by this factor}

\item{coremax_triangle_width_factor}{expand khsell max rectangles width by this factor}

\item{paint_outsiders}{paint species not connected to giant component}

\item{displace_outside_component}{displace outsider species (horizontal, vertical)}

\item{outsiders_separation_expand}{multiply by this factor outsiders' separation}

\item{outsiders_legend_expand}{displace outsiders legend}

\item{specialistskcore2_horizontal_dist_rootleaf_expand}{expand horizontal distance of specialist tail root node connected to kshell 2}

\item{specialistskcore2_vertical_dist_rootleaf_expand}{expand vertical distance of specialist tails connected to kshell 2}

\item{specialists_boxes_separation_count}{specialist species boxes separation count}

\item{root_specialist_expand}{expand root specialist distances of tails connected to kshell <> 2}

\item{hide_plot_border}{hide border around the plot}

\item{rescale_plot_area}{full plot area rescaling (horizontal, vertical)}

\item{kcore1specialists_leafs_vertical_separation}{expand vertical separation of specialist tails connected to kshell 1 species}

\item{corebox_border_size}{width of kshell boxes}

\item{kcore_species_name_display}{display species names of  shells listed in this vector}

\item{kcore_species_name_break}{allow new lines in species names of  shells listed in this vector}

\item{shorten_species_name}{number of characters of species name to display}

\item{exclude_species_number}{do not include species number in species}

\item{label_strguilda}{string labels of guild a}

\item{label_strguildb}{string labels of guild b}

\item{landscape_plot}{paper landscape configuration}

\item{backg_color}{plot background color}

\item{show_title}{show plot title}

\item{show_legend}{show plot legend position}

\item{use_spline}{use splines to draw links}

\item{spline_points}{number of points for each spline}

\item{file_name_append}{a label that the user may append to the plot file name for convenience}

\item{svg_scale_factor}{only for interactive apps, do not modify}

\item{weighted_links}{function to add link weight: 'none', 'log10' , 'ln', 'sqrt'}

\item{square_nodes_size_scale}{scale nodes area of kcore1 and outsiders}

\item{move_all_SVG_up}{move up all the SVG plot by this fraction, useful to crop upper white space}

\item{move_all_SVG_right}{move rightwards all the SVG plot by this fraction, useful to crop right white space}

\item{progress}{only for interactive apps, do not modifiy}

\item{lsize_kcorebox}{default kshell boxes label size}
}
\description{
This function plots the ziggurat graph of a bipartite network. Configuration parameters and
results are stored in a global environment called zgg. This environment is not destroyed
after the function is executed, so the developer can store it in a configuration file, retrieve
network analysis variables and so on. If a new ziggurat_graph is called, zgg is destroyed and
created again for the new plot. Plotting options are explained in the user manual.
}
\examples{
ziggurat_graph("data/","M_PL_001.csv",plotsdir="grafresults/",print_to_file = TRUE)
}
