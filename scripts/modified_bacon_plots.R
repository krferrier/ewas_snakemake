library(tibble)
library(cowplot)
traces_ggplot <- function(object, burnin=TRUE, index=1){
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
          strip.placement = "outside") +
    xlab("Iteration") +
    ylab("Trace")
  g
  if (!burnin) {
    g <- g + xlim(c((object@nburnin+1), object@niter))
  }
  
  print(g)
}

traces_ggplot(bc) + labs(title="title")
