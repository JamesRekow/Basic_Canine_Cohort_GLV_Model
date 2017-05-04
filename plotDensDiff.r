#  James Rekow

plotDensDiff = function(doc1, doc2, titleText = NULL){
  
  #  ARGS: doc1 - first input doc. Regions where doc2 has a higher density than doc1 will be
  #               colored blue
  #        doc2 - second input doc. Regions where doc2 has a higher density than doc1 will be
  #               colored red.
  #        title - optional input character string to use as the title of the output plot
  #
  #  RETURNS: a plot of the difference between the 2D densities of doc1 and doc2. The plotted
  #           density is that of density(doc2) - density(doc1)
  
  #  import libraries
  library(MASS)
  library(reshape2)
  library(ggplot2)
  
  #  store default title text if none is supplied
  if(is.null(titleText)){
    titleText = "Difference in Density of OD Points"
  } #  end if
  
  #  calculate the common overlap and dissimilarity range
  overRng = range(c(doc1[[1]], doc2[[1]]))
  dissRng = range(c(doc1[[2]], doc2[[2]]))
  
  #  calculate the 2D density estimate over the common range
  d1 = kde2d(doc1[[1]], doc1[[2]], lims = c(overRng, dissRng), n = 200)
  d2 = kde2d(doc2[[1]], doc2[[2]], lims = c(overRng, dissRng), n = 200)
  
  #  calculate the difference between the 2D density estimates
  diff12 = d1
  diff12$z = d2$z - d1$z
  
  #  melt data into long format
  #  first, add row and column names (x and y grid values) to the z-value matrix
  rownames(diff12$z) = diff12$x
  colnames(diff12$z) = diff12$y
  
  #  now melt it to long format
  diff12m = melt(diff12$z, id.var = rownames(diff12))
  names(diff12m) = c("overlap", "dissimilarity", "densDiff")
  lims = c(-1, 1) * 10
  
  diff12m$densDiff[diff12m$densDiff < lims[1]] = lims[1]
  diff12m$densDiff[diff12m$densDiff > lims[2]] = lims[2]

  #  plot difference between the two densities
  p = ggplot(diff12m, aes(overlap, dissimilarity, fill = densDiff)) +
    geom_raster() +
    #stat_contour(aes(col = ..level..), binwidth = 0.1) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = lims) +
    #scale_colour_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = lims) +
    coord_cartesian(xlim = overRng, ylim = dissRng, expand = c(0, 0)) #+
    #guides(colour = FALSE) +
    # ggtitle(titleText)
    #labs(title = titleText) +
    #theme(plot.title = element_text(hjust = 0.5))
  
  #  NOTE: the commented out lines above can be added in to modify the display, but are
  #  currently left out because they clutter the plot or add unnecessary visual detail that
  #  bloats the size of the resulting image file
  
  return(p)
  
} #  end plotDensDiff function
