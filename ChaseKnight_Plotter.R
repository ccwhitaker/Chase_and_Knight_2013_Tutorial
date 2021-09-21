# Plotting function - We will simplify the process of making our charts to a
# single function call. This will save us from the normally repetitive process
# of making figures.
plotter <- function(contr                                # control
                    , Tx                                 # treatment
                    , x
                    , metric = c("rich", "ENS_PIE")      
                    , meth = c("accum", "effect")
                    , title){

# We want to differentiate between our plots. Here we will be plotting our
# accumulation plots, for both richness and ENS_PIE
  if(meth == "accum"){
    plot(x = log(x), y = if(metric == "rich"){contr$sac.richness} else{contr$ENS_PIE}
         , ylim = c(0, 50)
         , col = "blue"
         , pch = 17
         , xlab = "Log(Number of Quadrats)"
         , ylab = ifelse(metric == "rich"
                         , "Species Richness"
                         , "ENS_PIE")
         , main = title
         , las = 1)
    points(log(x)
           , y = if(metric == "rich"){Tx$sac.richness}else{Tx$ENS_PIE}
           , col = "red", pch = 15)

    legend("bottomright"
       , legend=c("Control", "Treatment")
       , col=c("blue", "red")
       , pch = c(17, 15)
       , cex = 0.6
       )

  }
# Now we will plot the log effect size of our treatments on the simulated
# communities.
  if(meth == "effect"){
    effect_size = if(metric == "rich"){log(contr$sac.richness) - log(Tx$sac.richness)}
    else{log(contr$ENS_PIE) - log(Tx$ENS_PIE)}
    plot(log(x)
         , effect_size
         , ylim = c(-0.5, 3)
         , col = "black"
         , pch = 17
         , xlab = "Log(Number of Quadrats)"
         , ylab = "Log Ratio Effect Size"
         , main = title
         , las = 1)
  }
}
