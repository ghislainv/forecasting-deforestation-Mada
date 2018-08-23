require(broom)
require(glue)
require(rasterVis)
require(scales)
require(grid)

# Theme
## Setting basic theme options for plot with ggplot2
theme_base <- theme(axis.line=element_blank(),
										axis.text.x=element_blank(),
										axis.text.y=element_blank(),
										axis.ticks=element_blank(),
										axis.title.x=element_blank(),
										axis.title.y=element_blank(),
										legend.position="none",
										plot.margin=unit(c(0,0,0,0),"null"),
										panel.spacing=unit(c(0,0,0,0),"null"),
										plot.background=element_rect(fill="transparent"),
										panel.background=element_rect(fill="transparent"),
										panel.grid.major=element_blank(),
										panel.grid.minor=element_blank(),
										panel.border=element_blank())

# Zooms
zoom1 <- list(xmin=346000,xmax=439000,ymin=7387000,ymax=7480000)
zoom2 <- list(xmin=793000,xmax=886000,ymin=7815000,ymax=7908000)
Extent.zoom1 <- paste(zoom1$xmin,zoom1$ymin,zoom1$xmax,zoom1$ymax)
Extent.zoom2 <- paste(zoom2$xmin,zoom2$ymin,zoom2$xmax,zoom2$ymax)

# rho_plot()
rho_plot <- function(input_raster, input_vector, output_file, quantiles_legend=c(0.025,0.975), ...) {
	# Rho limits for legend
	rho_quantiles <- quantile(values(input_raster),quantiles_legend,na.rm=TRUE) 
	rho_bound <- max(sqrt(rho_quantiles^2))
	rho_limits <- c(-rho_bound,rho_bound)
	# Call to ggplot
	rasterVis::gplot(input_raster) +
		geom_raster(aes(fill=value)) +
		scale_fill_gradientn(colours=c("forestgreen","yellow","red"),na.value="transparent",
												 limits=rho_limits, oob=scales::squish) +
		geom_polygon(data=broom::tidy(input_vector), aes(x=long, y=lat, group=id), colour="black", fill="transparent", size=0.3) +
		theme_bw() + theme_base + coord_equal()
		#coord_equal(xlim=c(300000,1100000),ylim=c(7165000,8685000))
	# Save plot
	ggsave(output_file, ...)
	# Return message
	return(glue("Results plotted to \"{output_file}\""))
}
