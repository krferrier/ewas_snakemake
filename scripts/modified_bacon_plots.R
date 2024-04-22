setGeneric("ggtraces", function(object, burnin=TRUE, index=1){ standardGeneric("ggtraces") })
setMethod( "ggtraces", "Bacon", function(object, burnin=TRUE, index=1){
  if(burnin)
    gstraces <- object@traces[,,index]
  else
    gstraces <- object@traces[-c(1:object@nburnin),,index]
  
  thetahat_df <- as.data.frame(estimates(object)[index,]) %>% 
    rownames_to_column("variable") %>%
    rename(value = "estimates(object)[index, ]")

  gstraces_df <- as.data.frame(gstraces)
  gstraces_df$iteration <- 1:nrow(gstraces_df)
  gstraces_melted <- melt(gstraces_df, id.vars = "iteration") %>%
    mutate(variable = as.character(variable))

  g <- ggplot(gstraces_melted, aes(x = iteration, y = value)) +
    geom_line() +
    geom_hline(data = thetahat_df, aes(yintercept=value), color = "red") +
    facet_wrap(~variable, scales = "free_y", strip.position = "left") +
    scale_x_continuous(labels = c("",1000, "", 3000, "", 5000)) +
    theme_cowplot(font_size = 12) +
    theme(strip.background = element_blank(),
          strip.placement = "outside",
          plot.margin = margin(0.15, 0.25, 0.15, 0.15, "in")) +
    xlab("Iteration") +
    ylab("Trace")
  g
  if (!burnin) {
    g <- g + xlim(c((object@nburnin+1), object@niter))
  }
  
  return(g)
})

setGeneric("ggposteriors", function(object,
                                  thetas = c("sigma.0", "p.0"), index = 1,
                                  alphas=c(0.95, 0.9, 0.75), xlab="", ylab="", ...){
  standardGeneric("ggposteriors")
})
setMethod("ggposteriors", "Bacon", function(object, thetas, index, alphas, xlab, ylab, ...){
  
  if(any(!(thetas %in% colnames(object@traces[,, index]))))
    stop("'thetas' should be two of: ", paste(colnames(object@traces[,, index]), collapse=", "), "!")
  
  gstraces <- object@traces[-c(1:object@nburnin), thetas, index]
  
  if(xlab=="") xlab <- thetas[1]
  if(ylab=="") ylab <- thetas[2]
  
  df <- data.frame(x = gstraces[,1], y = gstraces[,2])
  est_df <- data.frame(x = estimates(object)[index, thetas[1]], y = estimates(object)[index, thetas[2]])
  
  # Plot using ggplot
  p <- ggplot(df, aes(x = x, y = y)) +
    geom_point(shape = 20) +
    stat_ellipse(level = 0.95, col = "blue") +
    stat_ellipse(level = 0.9, col = "blue") +
    stat_ellipse(level = 0.75, col = "blue") +
    labs(x = xlab,
         y = ylab) + 
    ggtitle(paste("median at:", round(estimates(object)[thetas], 3))) +
    geom_point(data = est_df, aes(x = x, y = y), color = "red", shape = 17, size = 4)
  
  return(p)
})


setGeneric("ggfit", function(object, index=1, ...){ standardGeneric("ggfit") })
setMethod("ggfit", "Bacon", function(object, index, col="grey75", border="grey75", ...){
  gg_plotnormmix(tstat(object, corrected=FALSE)[, index], estimates(object)[index, ], ...)
})

gg_plotnormmix <- function(x, theta, ...){
  x <- data.frame(x = x)
  theta <- data.frame(y = theta)
  fit <- ggplot(x, aes(x=x , y = after_stat(density))) +
    geom_histogram(fill = "grey", color="black", binwidth = 1.5) +
    geom_line(aes(x=x, y =dnorm(x, mean(x), sd(x))), lwd=1) +
    geom_line(aes(x=x, y=theta["p.0",]*dnorm(x, theta["mu.0",], theta["sigma.0",])), 
              color="red",
              lwd = 1.5) +
    geom_line(aes(x=x, y=theta["p.1",]*dnorm(x, theta["mu.1",], theta["sigma.1",])), 
              color="green",
              lwd=1.5) +
    geom_line(aes(x=x, y=theta["p.2",]*dnorm(x, theta["mu.2",], theta["sigma.2",])), 
              color="blue",
              lwd=1.5) +
    theme_cowplot(font_size = 12) +
    xlab("test statistics") +
    ylab("Density")
  
  return(fit)
} 
