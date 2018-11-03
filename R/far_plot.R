require(broom)
require(glue)
require(rasterVis)
require(scales) # for scales::squish()
require(grid) # for grid::unit()

# Theme
## Setting basic theme options for plot with ggplot2
theme_base <- theme(axis.line=element_blank(),
										axis.text.x=element_blank(),
										axis.text.y=element_blank(),
										axis.ticks=element_blank(),
										axis.title.x=element_blank(),
										axis.title.y=element_blank(),
										legend.position="none",
										plot.margin=grid::unit(c(0,0,0,0),"null"),
										panel.spacing=grid::unit(c(0,0,0,0),"null"),
										plot.background=element_rect(fill="transparent"),
										panel.background=element_rect(fill="transparent"),
										panel.grid.major=element_blank(),
										panel.grid.minor=element_blank(),
										panel.border=element_blank())

# rho_plot()
rho_plot <- function(input_raster, input_vector, output_file, quantiles_legend=c(0.025,0.975), ...) {
	# Rho limits for legend
	rho_quantiles <- quantile(values(input_raster),quantiles_legend,na.rm=TRUE) 
	rho_bound <- max(sqrt(rho_quantiles^2))
	rho_limits <- c(-rho_bound,rho_bound)
	# Call to ggplot
	p <- rasterVis::gplot(input_raster) +
		geom_raster(aes(fill=value)) +
		scale_fill_gradientn(colours=c("forestgreen","yellow","red"),na.value="transparent",
												 limits=rho_limits, oob=scales::squish) +
		geom_polygon(data=broom::tidy(input_vector), aes(x=long, y=lat, group=id), colour="black", fill="transparent", size=0.3) +
		theme_bw() + theme_base + coord_fixed()
	# Save plot
	ggsave(output_file, ...)
	# Return message an plot
	cat(glue("Results plotted to \"{output_file}\"\n"))
	return(p)
}

# diff_plot()
diff_plot <- function(input_raster, input_vector, output_file, maxpixels=50000, ext=NULL, rect=NULL, ...) {
	# Crop raster
	if (!is.null(ext)) {
		input_raster <- crop(input_raster,ext)
		e <- extent(ext)
		xlim = c(xmin(e),xmax(e))
		ylim = c(ymin(e),ymax(e))
	} else {
		xlim = NULL
		ylim = NULL
	}
	# Call to ggplot
	p <- rasterVis::levelplot(input_raster, maxpixels=maxpixels) +
		geom_raster(aes(fill=factor(value))) +
		scale_fill_manual(values=c("red","forestgreen","darkblue","lightblue"), na.value="transparent") +
		{if (!is.null(rect))
			geom_rect(data=rect, inherit.aes=FALSE, aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,group=id),
								fill="transparent", colour="black", size=0.3)
		} +
		geom_polygon(data=broom::tidy(input_vector), inherit.aes=FALSE, aes(x=long, y=lat, group=id),
								 colour="black", fill="transparent", size=0.3) +
		theme_bw() + theme_base + coord_fixed(ratio=1, xlim, ylim, expand=FALSE)
	# Save plot
	ggsave(output_file, p, ...)
	# Return message and plot
	cat(glue("Results plotted to \"{output_file}\""))
	return(p)
}

# End