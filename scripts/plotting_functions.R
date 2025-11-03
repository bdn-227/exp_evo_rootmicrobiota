library(ape)
library(ggplot2)
library(reshape)
library(scales)
library(grid)
library(gridExtra)
library(stringr)
library(tidyverse)


colAddAlpha <- function(col, alpha=1) { 
    
    x <- as.vector(col2rgb(col)) / 255
    
    return(rgb(x[1], x[2], x[3], alpha))

}

plot_btree <- function(tree, output.file="Rplots.pdf", cex.labels=1, label.offset=5.3,
			           type="fan", use.edge.length=T, edge.width=1, edge.color="black",
                       align.labels=T, color.labels="black", color.alg.lines="dark grey",
                       tips.bg="black", tips.col="black", tips.shape=16, tips.cex=1,
                       barplots.height=1, barplots.width, barplots.col, barplots=F, barplots.offset=0.025,
                       margins=c(.1, .1, .1, .1), height=8.27, width=11.7) {
	if (grepl(".tiff", output.file)) {

        tiff(file=output.file, height=height, width=width)

    } else {

        pdf(file=output.file, height=height, width=width, useDingbats=F)

    }

    if (align.labels) {
        tip.color <- "transparent"
    } else {
        tip.color <- color.labels
        color.labels <- "transparent"
    }

    par(mai=margins)

    plot.phylo(tree, type=type, tip.color=tip.color, use.edge.length=use.edge.length,
               label.offset=label.offset, show.node.label=F, edge.width=edge.width,
               cex=cex.labels, edge.color=edge.color)

	lastPP <- get("last_plot.phylo", envir=.PlotPhyloEnv)
	XX <- lastPP$xx
	YY <- lastPP$yy
    Ntip <- lastPP$Ntip
    xx.tips <- XX[1:Ntip]
    yy.tips <- YY[1:Ntip]

    # align tip labels
    
	if(type=="fan") {

    	angle <- atan2(yy.tips, xx.tips)
    	# xx.labels <- label.offset * cos(angle)
    	# yy.labels <- label.offset * sin(angle)

        xx.labels <- xx.tips + label.offset * cos(angle)
    	yy.labels <- yy.tips + label.offset * sin(angle)

		par(xpd=NA)
		sapply(1:length(xx.labels), function(i) lines(c(xx.tips[i], xx.labels[i]),
                                                      c(yy.tips[i], yy.labels[i]),
                                                      lty=3, col=color.alg.lines,
                                                      cex=3))
		
		s <- xx.labels < 0
		angle <- angle * 180 / pi
		angle[s] <- angle[s] + 180
		adj <- rep(0, length(tree$tip.label))
		adj[which(xx.labels < 0)] <- 1
	
	    # sapply(1:length(tree$tip.label), function(i) text(x=xx.labels[i], y=yy.labels[i],
        #                                                   labels=tree$tip.label[i], 
        #                                                   srt=angle[i], offset=1, cex=cex.labels,
        #                                                   adj=adj[i], col=color.labels))

        # TODO: plot N rings 
        # plot ring
        if (barplots) {
            for (j in 1:dim(barplots.height)[2]) {
            for (i in 1:dim(barplots.height)[1]) {
      
                w <- barplots.width
                dx <- xx.labels[i]
                dy <- yy.labels[i]
                theta <- atan(dy/dx)
                x1<- dx + (w/2) * cos(pi/2 - theta)
                y1<- dy - (w/2) * sin(pi/2 - theta)
                x2<- dx - (w/2) * cos(pi/2 - theta)
                y2<- dy + (w/2) * sin(pi/2 - theta)
                s <- if(dx > 0) .1 else -.1
                x3<- s * barplots.height[i, j] * cos(theta) + x2
                y3<- s * barplots.height[i, j] * sin(theta) + y2
                x4<- s * barplots.height[i, j] * cos(theta) + x1
                y4<- s * barplots.height[i, j] * sin(theta) + y1
            
                polygon(c(x1, x2, x3, x4), c(y1, y2, y3, y4),
                        col=barplots.col[i, j], border="transparent", lwd=.1)

                offset <- 0.1
                xx.labels[i] <- s * offset * cos(theta) + xx.labels[i]
                yy.labels[i] <- s * offset * sin(theta) + yy.labels[i]
                   
            }
            }
        }

	} else if(type=="phylogram") {

        xx.labels=xx.tips + label.offset
        yy.labels=yy.tips + label.offset
        
        sapply(1:length(xx.labels), function(i) lines(c(xx.tips[i], max(xx.tips) + label.offset),
                                                      c(yy.tips[i], yy.tips[i]),
                                                      lty=1, col=color.alg.lines, lwd=edge.width-0.5))


        sapply(1:length(tree$tip.label), function(i) text(x=max(xx.labels)+label.offset, y=yy.labels[i],
                                                          labels=tree$tip.label[i], 
                                                          offset=1, cex=cex.labels,
                                                          col=color.labels, adj = c(0,0)))

        if (barplots) {
    
            for (i in 1:dim(barplots.height)[2]) {
            for (j in 1:dim(barplots.height)[1]) {
    
                w <- barplots.width
                x1 <- x2 <- max(xx.labels) + barplots.offset * i + label.offset
                x3 <- x4 <- max(xx.labels) + barplots.offset * i + label.offset + as.numeric(barplots.height[j, i])
                y1 <- y4 <- yy.tips[j] - w
                y2 <- y3 <- yy.tips[j] + w
                
                if (barplots.height[j, i] != 0) {
    
    	        par(xpd=NA)
                polygon(c(x1, x2, x3, x4), c(y1, y2, y3, y4),
                        col=barplots.col[j, i], border="black", lwd=.5)
                
                }
    
            }
            }

        }

    } else if(type=="unrooted"){

        sapply(1:length(tree$tip.label), function(i) text(x=xx.tips[i] + label.offset, y=yy.tips[i],
                                                          labels=tree$tip.label[i], 
                                                          offset=1, cex=cex.labels,
                                                          col=color.labels))


    }

    # plot leaf shapes

    for (i in 1:length(xx.tips)) {

        points(x=xx.tips[i], y=yy.tips[i],
               col=tips.col[i], pch=tips.shape[i],
               cex=tips.cex[i], bg=tips.bg[i], lwd=0.25)

    }  

	dev.off()

}

# colors

alpha <- .7
c_yellow <-          rgb(255 / 255, 255 / 255,   0 / 255, alpha)
c_blue <-            rgb(  0 / 255, 000 / 255, 255 / 255, alpha)
c_orange <-          rgb(255 / 255,  69 / 255,   0 / 255, alpha)
c_green <-           rgb(  50/ 255, 220 / 255,  50 / 255, alpha)
c_dark_green <-      rgb( 50 / 255, 200 / 255, 100 / 255, alpha)
c_very_dark_green <- rgb( 50 / 255, 150 / 255, 100 / 255, alpha)
c_sea_green <-       rgb( 46 / 255, 129 / 255,  90 / 255, alpha)
c_black <-           rgb(  0 / 255,   0 / 255,   0 / 255, alpha)
c_grey <-            rgb(180 / 255, 180 / 255,  180 / 255, alpha)
c_dark_brown <-      rgb(101 / 255,  67 / 255,  33 / 255, alpha)
c_red <-             rgb(200 / 255,   0 / 255,   0 / 255, alpha)
c_dark_red <-        rgb(255 / 255, 130 / 255,   0 / 255, alpha)

# ggplot2 theme

main_theme <- theme(panel.background=element_blank(),
              panel.grid=element_blank(),
              axis.line.x=element_line(color="black"),
              axis.line.y=element_line(color="black"),
              axis.ticks=element_line(color="black"),
              axis.text=element_text(colour="black", size=14),
              axis.title = element_text(size = 15),
              legend.text = element_text(size = 14),
              legend.title=element_text(size=15),
              legend.background=element_blank(),
              legend.key=element_blank(),
              text=element_text(family="sans"),
              strip.text.x = element_text(size = 16), 
    	      strip.background = element_rect(fill="white"))


facette_theme <- main_theme + 
                 theme(plot.title = element_text(hjust = 0.5), 
                        panel.border = element_rect(color = "black", fill = NA, size = 1),
                        strip.background = element_rect(fill="black"),
                        strip.text = element_text(colour = 'white'),
                        strip.text.y = element_text(size = 14))

facette_theme_white <- theme(panel.background=element_blank(),
                                panel.grid=element_blank(),
                                axis.line.x=element_line(color="black"),
                                axis.line.y=element_line(color="black"),
                                axis.ticks=element_line(color="black"),
                                axis.text=element_text(colour="black", size=14),
                                axis.title = element_text(size = 15),
                                legend.text = element_text(size = 14),
                                legend.title=element_text(size=15),
                                legend.background=element_blank(),
                                legend.key=element_blank(),
                                text=element_text(family="sans"),
                                strip.text.x = element_text(size = 16), 
                                strip.background = element_rect(fill="white")) + 
                        theme(plot.title = element_text(hjust = 0.5), 
                                panel.border = element_rect(color = "white", fill = NA, size = 1),
                                strip.background = element_rect(fill="white"),
                                strip.text = element_text(colour = 'black'),
                                strip.text.y = element_text(size = 14))

tree_theme <- theme(
                    panel.background=element_blank(),
                    panel.grid=element_blank(),
                    axis.line.x=element_line(color="black"),
                    axis.line.y=element_blank(),
                    axis.ticks=element_line(color="black"),
                    axis.text=element_text(colour="black", size=14),
                    axis.title = element_text(size = 15),
                    legend.text = element_text(size = 14),
                    legend.title=element_text(size=15),
                    legend.background=element_blank(),
                    legend.key=element_blank(),
                    text=element_text(family="sans"),
                    strip.text.x = element_text(size = 16), 
                    strip.background = element_rect(fill="white")
                    )


colors_host <- data.frame(group=c("col0",    "gifu",  "soil", "collection"),
                     color=c("#f8756b", "#00b8e3", "#654321", "black"))
                     
colors_syncom <- data.frame(group=c("at",    "lj",  "soil", "axenic"),
                     color=c("#f8756b", "#00b8e3", "#654321", "black"))
                     


# read the tables
color_family = read.delim("mapping/colors_family.tsv", header = TRUE, sep = "\t")
taxo_table = readr::read_tsv("mapping/taxonomy.tsv", show_col_types = FALSE)
tax_levels = c("family", "order", "class", "phylum", "kingdom")
cols = c("strain", tax_levels)
taxo_table = taxo_table[cols]
taxo_table$strain_family = glue("{taxo_table$strain} - {taxo_table$family}")
colnames(taxo_table) = c("strain", tax_levels, "strain_family")



# get color_df
color_strain_df = taxo_table[, c("strain", "family")] %>% 
                left_join(color_family, by = c("family" = "Family"))
color_strain = color_strain_df$Color
names(color_strain) = color_strain_df$strain
color_strain = color_strain[!is.na(color_strain)]


# adjust the names
color_strain_family = color_strain
names(color_strain_family) = taxo_table$strain_family[match(names(color_strain_family), taxo_table$strain)]


# color family
color_fam = color_family %>% 
                    filter(Family != "unkown") %>% 
                    select(Family, Color) %>% 
                    deframe()



c25 = c(
      "dodgerblue2", 
      "#E31A1C", # red
      "green4",
      "#6A3D9A", # purple
      "#FF7F00", # orange
      "black", 
      "gold1",
      "skyblue2", 
      "#FB9A99", # lt pink
      "palegreen2",
      "#CAB2D6", # lt purple
      "#FDBF6F", # lt orange
      "gray70", 
      "khaki2",
      "maroon", 
      "orchid1", 
      "deeppink1", 
      "blue1", 
      "steelblue4",
      "darkturquoise", 
      "green1", 
      "yellow4", 
      "yellow3",
      "darkorange4", 
      "brown"
)




save_fig = function(p, f_name, width=7, height=7){

  #save pdf
  f = paste0(f_name, ".pdf")
  ggsave(f, p, width = width, height = height, limitsize = F)

  f = paste0(f_name, ".svg")
  ggsave(f, p, width = width, height = height, limitsize = F)

  #save png
  f = paste0(f_name, ".png")
  ggsave(f, p, width = width, height = height, limitsize = F, dpi=1200)
}






