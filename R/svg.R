###############################################################################
#   PFC Universidad Politécnica de Madrid - EUITT
#   Representación gráfica de redes bipartitas basadas en descomposición k-core
#
# Autor         : Juan Manuel García Santi / Javier Garcia Algarra
# Módulo        : svg.R
# Descricpción  : Funciones básicas para la generación de un gráfico en formato
#                 SVG (Scalable Vectors Graphics). Contiene las funciones
#                 necesarias para generar un SVG con rectángulos, rutas y
#                 segmentos, y proporcionar o almacenar el XML correspondiente
#                 al SVG generado
###############################################################################
library(ggplot2)
library(rlang)
#' SVG aux function
#'
#'
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords cats
#' @export
#' @examples
#' SVG()

SVG<-function(scale_factor,style="ziggurat",nnodes=50,flip_coordinates=FALSE) {
  plottype = style
  fontscale = ifelse (style=="ziggurat", (3+(nnodes<50)), 30)
  flinkscale = ifelse (style=="ziggurat", 1, 100)
  # crea el objeto SVG
  this<-list(content=c(""), minx=0, miny=0, maxx=0, maxy=0, scale_factor=scale_factor, font_scale_factor=fontscale)
  
  # guarda el contenido del svg en un fichero
  this$save <- function(fileName, svg) {
    fileConn<-file(fileName)
    header<-"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
    writeLines(paste0(header, this$html()), fileConn)
    close(fileConn)
  }
  
  # devuelve el HTML correspondiente al objeto
  zr <- 1
  zl <- 1
  zu <- 1
  this$html<-function() {
    if (style=='ziggurat'){
      zr <- zgg$move_all_SVG_right
      zl <- zgg$landmark_right
      zu <- zgg$move_all_SVG_up
    }
    else{
      zr <- 0
      zl <- bpp$landmark_right
      zu <- 0
    }
    #redondea el viewBox a la decena mas cercana
    minx  <- (1+zr)*floor(this$minx/10)*10
    
    if (style=='ziggurat'){
      maxx  <- min(ceiling(this$maxx/10)*10, ceiling(1.1*zl/100)*10)
      miny  <- (1-zu)*floor(this$miny/10)*10
      maxy  <- ceiling(this$maxy/10)*10
    } else {
      maxx  <- ceiling(4*zl/100)*10
      miny  <- (bpp$landmark_bottom/10)*10
      if (flip_coordinates)
        miny  <- ((bpp$landmark_bottom-4*bpp$xstep)/10)*10
      maxy  <-(bpp$landmark_top/10)*10
      if (flip_coordinates)
        maxy  <- ceiling(4*this$maxy/10)*10
    }
    fzcale <- 1.05
    tleftx <- minx
    tlefty <- miny
    swidth <- ceiling((maxx-minx)/10)*10
    sheight <- ceiling((maxy-miny)/10)*10
    if (flip_coordinates){
      if (style=='ziggurat'){
        viewBox<-paste0(tleftx, " ", ceiling(adjustleft+1.2*(tleftx-tlefty)/10)*10, " ", 1.2*swidth, " ", 1.2*swidth)
      } else {
        if (style=="chilopod")
          viewBox<-paste0(tleftx, " ", (tleftx-tlefty)+bpp$xstep*1.3, " ", 1.2*swidth, " ", 1.4*swidth)
        else{
          fmulth = 1.2
          fmultw = 1.2
          viewBox<-paste0(tleftx, " ", (tleftx-tlefty)+bpp$xstep*1.3, " ", fmultw*swidth, " ", fmulth*swidth)
        }
      }
      svg0<-paste0("<svg transform='rotate(90)' xmlns=\"http://www.w3.org/2000/svg\" preserveAspectRatio=\"xMidYMid meet\" viewBox=\"", viewBox, "\">\n")
    }
    else{
      if (style!='ziggurat'){
        wv <- (1.1+0.1*(nnodes%/%100))*swidth
        if (style=='chilopod'){
          tleftx <- tleftx-bpp$xstep
          h <-  wv*0.5
          vadjust <- -2*bpp$xstep
        }
        else{
          fdivh = ifelse(min(bpp$num_a_coremax,bpp$num_b_coremax)<=9,3,5)
          h <-  wv/fdivh
          vadjust <- 0
        }
        viewBox<-paste0(1.02*tleftx, " ", (tlefty+vadjust), " ", wv, " ", h)
      }
      else{
        viewBox<-paste0(tleftx, " ", tlefty, " ", swidth, " ", swidth)
      }
      svg0<-paste0("<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"", viewBox, "\">\n")
      #svg0<-paste0("<svg xmlns=\"http://www.w3.org/2000/svg\" transform='scale(",zoom_svgplot,")' viewBox=\"", viewBox, "\">\n")
    }
    
    svg1<-paste0("</svg>")
    return(paste0(svg0, paste0(this$content, collapse=""), svg1, sep=""))
  }
  
  # crea un rectangulo a partir de un conjunto de datos, con parametros similares
  # a ggplot2::geom_rect
  this$rect <- function(idPrefix, data, mapping, fill, alpha, color, size=0, linetype=1) {
    result <- ""
    
    # si solo se ha pasado un color lo utiliza para todos los datos
    if (length(color)==1) {
      color<-rep(color, nrow(data))
    }
    
    # evalua las posiciones
    # cambia de signo las coordenadas y, ya que en SVG el eje y es al contrario de lo que trata R con ggplot
    
    xmin  <- this$round_coords(eval_tidy(mapping$xmin, data)/this$scale_factor)
    xmax  <- this$round_coords(eval_tidy(mapping$xmax, data)/this$scale_factor)
    ymin  <- -this$round_coords(eval_tidy(mapping$ymin, data)/this$scale_factor)
    ymax  <- -this$round_coords(eval_tidy(mapping$ymax, data)/this$scale_factor)
    
    # itera para cada rectangulo
    for (i in 1:nrow(data)) {
      rect2<-this$rect2(id=paste0(plottype,idPrefix, "-", i, "-rect"), xmin=xmin[i], xmax=xmax[i], ymin=ymin[i], ymax=ymax[i], fill=fill[i], alpha=alpha, color=color[i], size=size, linetype=linetype)
      result<-paste0(result, rect2)
    }
    
    # incorpora el resultado al contenido del SVG
    this$content<<-cbind(this$content,c(result))
  }
  
  # funcion axiliar para crear un rectangulo
  this$rect2 <- function(id, xmin, xmax, ymin, ymax, fill, alpha, color, size, linetype) {
    result <- ""
    
    # traduce el "no-relleno"
    if (fill=="transparent") {
      fill<-"none"
    }
    
    # actualiza las coordenadas del viewBox si es necesario
    
    if (xmin<this$minx) this$minx<<-xmin
    if (xmin>this$maxx) this$maxx<<-xmin
    if (xmax<this$minx) this$minx<<-xmax
    if (xmax>this$maxx) this$maxx<<-xmax
    if (ymin<this$miny) this$miny<<-ymin
    if (ymin>this$maxy) this$maxy<<-ymin
    if (ymax<this$miny) this$miny<<-ymax
    if (ymax>this$maxy) this$maxy<<-ymax
    
    
    # dibuja el rectangulo
    result <- paste0(result, "<rect id=\"", id, "\" ")
    result <- paste0(result, "style=\"", "fill:", fill, ";fill-opacity:", alpha, "\" ")
    result <- paste0(result, "stroke=\"", color, "\" ")
    if (linetype>0 && linetype<7) {
      result <- paste0(result, "stroke-dasharray=\"", this$stroke_dasharray(linetype), "\" ")
    }
    result <- paste0(result, "stroke-width=\"", size, "\" ")
    
    result <- paste0(result, "x=\"", min(xmin, xmax), "\" ")
    result <- paste0(result, "y=\"", min(ymin, ymax), "\" ")
    result <- paste0(result, "width=\"", abs(xmax-xmin), "\" ")
    result <- paste0(result, "height=\"", abs(ymax-ymin), "\" ")
    result <- paste0(result, "/>\n")
    
    return(result)
  }
  
  # crea un texto a partir de un conjunto de datos, con parametros similares
  # a ggplot2::geom_text
  this$text <- function(idPrefix, data, mapping, label, color, size, angle=0) {
    if (flip_coordinates)
      angle = angle + 90
    result <- ""
    
    # si solo se ha pasado un color lo utiliza para todos los datos
    if (length(color)==1) {
      color<-rep(color, nrow(data))
    }
    
    # evalua las posiciones
    # cambia de signo las coordenadas y, ya que en SVG el eje y es al contrario de lo que trata R con ggplot
    
    
    x <- this$round_coords(eval_tidy(mapping$x, data)/this$scale_factor)
    y <- -this$round_coords(eval_tidy(mapping$y, data)/this$scale_factor)
    # itera para cada texto
    for (i in 1:nrow(data)) {
      text2<-this$text2(id=paste0(plottype,idPrefix, "-", i, "-text"), 
                        x=x[i], y=y[i], label=label[i], color[i], size, angle)
      result<-paste0(result, text2)
    }
    
    # incorpora el resultado al contenido del SVG
    this$content<<-cbind(this$content,c(result))
  }
  
  # funcion auxiliar para la creacion de texto
  # divide el texto generado en tantos tspan como saltos de linea tenga la etiqueta recibida
  this$text2 <- function(id, x, y, label, color, size, angle) {
    result<-""
    
    # actualiza las coordenadas del viewBox si es necesario
    # anyade el tamanyo correspondiente a toda la longitud del texto
    len  <- nchar(label)
    minx <- x-len*size
    maxx <- x+len*size
    miny <- y-len*size
    maxy <- y+len*size
    
    if (minx<this$minx) this$minx<<-minx
    if (maxx>this$maxx) this$maxx<<-maxx
    if (miny<this$miny) this$miny<<-miny
    if (maxy>this$maxy) this$maxy<<-maxy
    
    # agrupacion de texto
    result <- paste0(result, "<text id=\"", id, "\"", " ")
    result <- paste0(result, "y=\"", y, "\"", " ")
    if (angle!=0) {
      # cambia de signo el angulo, ya que se interpreta distinto que en ggplot
      result <- paste0(result, "transform=\"rotate(", -angle , " ", x, " ", y , ")\" ")
    }
    if (len < 4){
      halignstr = "middle"
      valignstr = "middle"
    }
    else {
      halignstr = "start"
      valignstr = "middle"
    }
    result <- paste0(result, "style=\"text-anchor:",halignstr,";dominant-baseline:",valignstr,
                     ";font-family:Arial;font-size:", size*this$font_scale_factor, "px;fill:", color, "\"")
    result <- paste0(result, ">\n")
    # tspan
    first   <- TRUE
    dy      <- 1.4*size*this$font_scale_factor
    labels  <- strsplit(label, "\n")[[1]]
    if (length(labels)>0) {
      for (i in 1:length(labels)) {
        if (nchar(labels[i])>0) {
          result  <- paste0(result, "<tspan ")
          result  <- paste0(result, "x=\"", x, "\"", " ")
          result  <- paste0(result, "dy=\"", ifelse(first, 0, dy), "\">")
          result  <- paste0(result, labels[i])
          result  <- paste0(result, "</tspan>\n")
          first   <- FALSE
        }
      }
    }
    
    # fin de la agrupacion de texto
    result <- paste0(result, "</text>\n")
    
    return(result)
  }
  
  # crea un segmento a partir de un conjunto de datos, con parametros similares
  # a ggplot2::geom_segment
  this$segment <- function(idPrefix, data, mapping, alpha, color, size=0, linetype=1) {
    result <- ""
    
    # si solo se ha pasado un color lo utiliza para todos los datos
    if (length(color)==1) {
      color<-rep(color, nrow(data))
    }
    
    # evalua las posiciones
    # cambia de signo las coordenadas y, ya que en SVG el eje y es al contrario de lo que trata R con ggplot

    x     <- this$round_coords(eval_tidy(mapping$x, data)/this$scale_factor)
    xend  <- this$round_coords(eval_tidy(mapping$xend, data)/this$scale_factor)
    y     <- -this$round_coords(eval_tidy(mapping$y, data)/this$scale_factor)
    yend  <- -this$round_coords(eval_tidy(mapping$yend, data)/this$scale_factor)
    # itera para cada segmento
    for (i in 1:nrow(data)) {
      segment2<-this$segment2(id=paste0(plottype,idPrefix, "-", i, "-segment"),
                              x=x[i], xend=xend[i], y=y[i], yend=yend[i],
                              alpha=alpha, color=color[i], 
                              size=flinkscale*size[i], 
                              linetype=linetype)
      result<-paste0(result, segment2)
    }
    
    # incorpora el resultado al contenido del SVG
    this$content<<-cbind(this$content,c(result))
  }
  
  # funcion axiliar para la creacion de un segmento
  this$segment2 <- function(id, x, xend, y, yend, alpha, color, size, linetype) {
    result <- ""
    
    # actualiza las coordenadas del viewBox si es necesario
    if (x<this$minx)    this$minx<<-x
    if (x>this$maxx)    this$maxx<<-x
    if (xend<this$minx) this$minx<<-xend
    if (xend>this$maxx) this$maxx<<-xend
    if (y<this$miny)    this$miny<<-y
    if (y>this$maxy)    this$maxy<<-y
    if (yend<this$miny) this$miny<<-yend
    if (yend>this$maxy) this$maxy<<-yend
    
    result <- paste0(result, "<g id=\"", id , "\" ")
    result <- paste0(result, "fill=\"none\" ")
    result <- paste0(result, "stroke=\"", color , "\" ")
    if (linetype>0 && linetype<7) {
      result <- paste0(result, "stroke-dasharray=\"", this$stroke_dasharray(linetype), "\" ")
    }
    result <- paste0(result, "stroke-width=\"", size , "\" ")
    
    
    result <- paste0(result, "stroke-opacity=\"", alpha , "\"")
    result <- paste0(result, ">\n")
    
    result <- paste0(result, "<path d=\"")
    result <- paste0(result, "M", x, " ", y, " ")
    result <- paste0(result, "L", xend, " ", yend, "\"")
    result <- paste0(result, "/>\n")
    result <- paste0(result, "</g>\n")
    return(result)
  }
  
  # crea una ruta a partir de un conjunto de datos, con parametros similares
  # a ggplot2::geom_path
  this$path <- function(idPrefix, data, mapping, alpha, color, size=0, linetype=1) {
    result <- ""
    # si solo se ha pasado un color lo utiliza para todos los datos
    if (length(color)==1) {
      color<-rep(color, nrow(data))
    }
    
    # evalua los grupos
    group <- eval_tidy(mapping$group, data)
    # itera para cada grupo de rutas
    # cambia de signo las coordenadas y, ya que en SVG el eje y es al contrario de lo que trata R con ggplot
    
    for (i in unique(group)) {
      g <- data[data[quo_text(mapping$group)]==i,]
      x <- this$round_coords(g[,c(quo_text(mapping$x))]/this$scale_factor)
      y <- -this$round_coords(g[,c(quo_text(mapping$y))]/this$scale_factor)
      path2<-this$path2(id=paste0(plottype,idPrefix, "-", i, "-path"), x=x, y=y, alpha=alpha, color=color[i], size=size, linetype=linetype)
      result<-paste0(result, path2)
    }
    
    # incorpora el resultado al contenido del SVG
    this$content<<-cbind(this$content,c(result))
    
  }
  
  # funcion auxiliar para la creacion de una ruta
  this$path2 <- function(id, x, y, alpha, color, size, linetype) {
    result <- ""
    
    # actualiza las coordenadas del viewBox si es necesario
    minx<-min(x)
    maxx<-max(x)
    miny<-min(y)
    maxy<-max(y)
    if (minx<this$minx) this$minx<<-minx
    if (minx>this$maxx) this$maxx<<-minx
    if (maxx<this$minx) this$minx<<-maxx
    if (maxx>this$maxx) this$maxx<<-maxx
    if (miny<this$miny) this$miny<<-miny
    if (miny>this$maxy) this$maxy<<-miny
    if (maxy<this$miny) this$miny<<-maxy
    if (maxy>this$maxy) this$maxy<<-maxy
    
    # crea la ruta
    result <- paste0(result, "<g id=\"", id , "\" ")
    result <- paste0(result, "fill=\"none\" ")
    result <- paste0(result, "stroke=\"", color , "\" ")
    if (linetype>0 && linetype<7) {
      result <- paste0(result, "stroke-dasharray=\"", this$stroke_dasharray(linetype), "\" ")
    }
    result <- paste0(result, "stroke-width=\"", size , "\" ")
    result <- paste0(result, "stroke-opacity=\"", alpha , "\"")
    result <- paste0(result, ">\n")
    
    result <- paste0(result, "<path d=\"")
    result <- paste0(result, "M", x[1], " ", y[1], " ")
    for (i in 2:length(x)) {
      result <- paste0(result, " L", x[i], " ", y[i])
    }
    result <- paste0(result, "\"/>\n")
    result <- paste0(result, "</g>\n")
    
    return(result)
  }
  
  # funcion auxiliar para especificar el tipo de linea a utilizar, a
  # partir de los mismos valores de linetype usados en ggplot2
  #   0:blank, 1:solid, 2:dashed, 3:dotted, 4:dotdash, 5:longdash, 6:twodash
  this$stroke_dasharray <- function(linetype) {
    linetypes<-list("1"=c(0), "2"=c(4,4), "3"=c(1,1), "4"=c(1,1,4,1), "5"=c(6,1), "6"=c(1,1,4,1))
    result<-paste(linetypes[[as.character(linetype)]], collapse=" ")
    return(result)
  }
  
  # funcion para redondeo de coordenadas y no arrastrar
  # todos los decimales al svg
  this$round_coords<- function (values) {
    return(round(values, digits=2))
  }
  
  class(this)<-c("svg")
  return(this)
}

# prueba de SVG
test <- function(items=2) {
  d<- data.frame(
    x1 = sample(600, items, replace = TRUE),
    x2 = sample(600, items, replace = TRUE),
    y1 = sample(600, items, replace = TRUE),
    y2 = sample(600, items, replace = TRUE),
    fill = rep(c("red", "green", "blue"), items)[1:items]
  )
  
  s1<-SVG(1)
  s1$rect(data=d, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill=d$fill, alpha=0.1, color=d$fill, size=1, linetype=2)
  s1$text(data=d, mapping=aes(x=x1, y=y1), label=paste0("(",d$x1,",",d$y1,")"), color="blue", size=9, angle=0)
  s1$text(data=d, mapping=aes(x=x1, y=y2), label=paste0("(",d$x1,",",d$y2,")"), color="blue", size=9, angle=0)
  s1$text(data=d, mapping=aes(x=x2, y=y2), label=paste0("(",d$x2,",",d$y2,")"), color="blue", size=9, angle=0)
  s1$text(data=d, mapping=aes(x=x2, y=y1), label=paste0("(",d$x2,",",d$y1,")"), color="blue", size=9, angle=0)
  s1$rect(data=d, mapping=aes(xmin=y1, xmax=y2, ymin=x1, ymax=x2), fill=d$fill, alpha=0.1, color=d$fill, size=1, linetype=3)
  s1$text(data=d, mapping=aes(x=y1, y=x1), label=paste0("(",d$y1,",",d$x1,")"), color="blue", size=9, angle=45)
  s1$text(data=d, mapping=aes(x=y1, y=x2), label=paste0("(",d$y1,",",d$x2,")"), color="blue", size=9, angle=45)
  s1$text(data=d, mapping=aes(x=y2, y=x2), label=paste0("(",d$y2,",",d$x2,")"), color="blue", size=9, angle=45)
  s1$text(data=d, mapping=aes(x=y2, y=x1), label=paste0("(",d$y2,",",d$x1,")"), color="blue", size=9, angle=45)
  s1$save("C:\\Temp\\kk.svg")
}

#test()