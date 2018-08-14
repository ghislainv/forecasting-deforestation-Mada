require(rgdal)

# Madagascar borders
mada <- readOGR(dsn="data/mada",layer="mada38s")
mada.df <- tidy(mada)

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
