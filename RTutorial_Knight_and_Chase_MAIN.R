# In this tutorial, we will use the package "mobsim" to generate biodiversity at
# scale and replicate a similar methodology to that used in the Knight and Chase
# (2013) paper!

# We will start things off with a toy model that gives us a sense
# of what we will be working through, then we will write two functions that
# will give us the tools to repeatedly simulate communities under different
# treatments. A couple of bulky functions that support the simulation and
# visualization but are not the focus of the tutorial have been placed in a
# supplemental script that we can pull from.
# Now with all that out of the way, let's get started!

# Install packages: "mobsim", "mobr", "vegan"
# install.package("mobsim")
# install.package("mobr")
# install.package("vegan")

# Load the required packages:
library(mobsim)
library(mobr)
library(vegan)

# There are some bulky functions that we can ignore for now, but will use later.
source("./ChaseKnight_Cumulative_ENS_PIE.R")
source("./ChaseKnight_Plotter.R")


# In the next few lines of code, we will run our toy simulation to get a feel
# for how the functions work. This is just a test run to get us used to
# simulating community data, so the numbers we are using are small (a simulated
# species pool of 10 species, drawn from a community of 320 individuals). We
# will repeatedly simulate a larger community later on.

# First we will simulate a SAD (Species Abundance Distribution) for our
# simulated community. The SAD follows a log-normal distribution, a continuous
# curve modelling species' relative abundances. To fit a community dataset to
# this curve, we convert a continuous distribution to a discrete dataset
# (integer number of species with integer abundances; you can't have fractions
# of an individual!)


# The coefficient of variation for species abundance (cv_abund) determines the
# size of the "spread" of species from high abundance to low abundance. A low
# number corresponds to a more even community, while a higher number corresponds to a more uneven community.

sad_lnorm <- sim_sad(s_pool = 10 # number of species
                     , n_sim = 320 # number of individuals
                     , sad_type = "lnorm" # theoretical distribution of relative abundance
                     , sad_coef = list(cv_abund = .5)) # parameters for theoretical distribution
# Our 'cv_abund' here of '.5' represents a moderately even community.


# 'sim_sad()' simply returns species abundances. To add biological realism, we
# will also generate locations for each individual


# Let's add some coordinates to our SAD to simulate a spatial distribution
# pattern. We will demonstrate both a random and an aggregated spatial
# distribution. This will allow us to sample our community similar to how a
# "real" community would have individuals spread out over a physical landscape,
# our simulated community will have individuals spread out over a coordinate
# plane, in either a randomized or an aggregated distribution.

# Here we will add spatially random coordinates to our SAD, in a grid space of 
# 2 x 2
practice_sad_coords <- sim_poisson_coords(sad_lnorm
                                          , xrange = c(0, 2)
                                          , yrange = c(0, 2))

# Alternately, here we will introduce intraspecific aggregation for a different
# resultant community landscape. Often when we study communities in the field,
# we find species aggregate with other like individuals, which has profound
# implications for measuring biodiversity at scale as the paper explored.
practice_sad_coords_agg <- sim_thomas_coords(sad_lnorm
                                             , xrange = c(0, 2)
                                             , yrange = c(0, 2))


# Let's take a look at what our simulated communities looks like!
par(mfrow = c(2, 2))
plot(practice_sad_coords, axes = FALSE) # The non-aggregated community
axis(1, at = 0:2)
axis(2, at = 0:2, las = 1)
plot(practice_sad_coords_agg, axes = FALSE) # The aggregated community
axis(1, at = 0:2)
axis(2, at = 0:2, las = 1)



# Now we want to divide our community into a grid and sample quadrats

# mobsim_sample_quadrats returns a data.frame. We will sample the data in 0.25
# x 0.25 quadrats, then display a plot of the sampling, and sample in a grid
# format. NOTE: This sampling method does not perform the nested sampling that
# the paper mentions.

practice_sample_data <- sample_quadrats(practice_sad_coords
                                        , n_quadrats = 64
                                        , quadrat_area = .0625
                                        , plot = TRUE
                                        , method = "grid"
                                        , x0 = 0
                                        , y0 = 0
                                        , delta_x = .25
                                        , delta_y = .25)

# Let's do the same again for our aggregated community.
practice_sample_data_agg <- sample_quadrats(practice_sad_coords_agg
                                            , n_quadrats = 64
                                            , quadrat_area = .0625
                                            , plot = TRUE
                                            , method = "grid"
                                            , x0 = 0
                                            , y0 = 0
                                            , delta_x = .25
                                            , delta_y = .25)



# Our resulting data.frames, practice_sample_data and practice_sample_data_agg,
# have two variables, $spec_dat and $xy_dat. $spec_dat is a matrix that has the
# abundance of individuals of each species in each sampled quadrat. $xy_dat is
# our coordinate plane, connecting coordinates to site (quadrat) number.
# We can see this by looking at the head of the matrices.
head(practice_sample_data$spec_dat)
head(practice_sample_data_agg$spec_dat) # Notice the numerical "clumping" of species in the aggregated community
head(practice_sample_data$xy_dat) # This should be the same between both communities, as we divided them up the same way.




# Now let's make a Species Accumulation Curve from the sampled quadrats we just simulated.
practice_sac <- specaccum(practice_sample_data$spec_dat, method = "collector")
practice_sac_agg <- specaccum(practice_sample_data_agg$spec_dat, method = "collector")
par(mfrow = c(1, 1))
# And here we will again plot the results.
x <- 1:64
plot(x
     , practice_sac$richness
     , col = "blue"
     , pch = 17
     , xlab = "Quadrats Sampled"
     , ylab = "Species Richness"
     , main = "Toy Simulation\nSpecies Accumulation Curves"
     , las = 1
     , ylim = c(0, 10))
points(x, practice_sac_agg$richness
       , col = "red"
       , pch = 15)
legend("bottomright"
       , inset = .01
       , legend=c("Control", "Treatment")
       , col=c("blue", "red")
       , pch = c(17, 15)
       , box.lty  = 0)

# Why are are these "curves" so "blocky"?
# Notice the shift in x-intercept of the treatment's SAC, as well as the
# shallower curve. Do these changes match what we would expect from the paper's
# predictions?


###############################################
## Checkpoint 1 ###############################
###############################################


# Now let's turn what we just did into something a little more representative of
# the paper's methods:

# We will write a function that takes parameters to easily repeat the simulation
# MANY times under different initial conditions (treatments). This will allow us
# to remove variation that we are not interested in.

treatment <- function(s_pool_input = 45
                      , n_sim_input = 32000
                      , cv_abund_input = 0.5
                      , agg = FALSE){
  

  # Here we are inputting our parameters into our simulated SAD
  sad_lnorm <- sim_sad(s_pool = s_pool_input
                       , n_sim = n_sim_input
                       , sad_type = "lnorm"
                       , sad_coef = list(cv_abund = cv_abund_input))

  if(agg == TRUE){
    sad_coords <- sim_thomas_coords(abund_vec = sad_lnorm
                                    , xrange = c(0, 20)
                                    , yrange = c(0, 20))


  }else{

    sad_coords <- sim_poisson_coords(abund_vec = sad_lnorm
                                     , xrange = c(0, 20)
                                     , yrange = c(0, 20))


  }

  sample_data <- sample_quadrats(sad_coords,
                                 n_quadrats = 400,
                                 quadrat_area = 1,
                                 plot = FALSE,
                                 method = "grid",
                                 x0 = 0,
                                 y0 = 0,
                                 delta_x = 1,
                                 delta_y = 1)

  return(sample_data)

}


###############################################
## Checkpoint 2 ###############################
###############################################


# We've now written a function that makes it easy to quickly simulate
# communities with different initial conditions, Our next step will be to write
# a function that allows us to easily set our initial conditions and run that
# simulation many times.

run_loop <- function(s_pool_input_treatment = 45
                     , n_sim_input_treatment = 32000
                     , cv_abund_input_treatment = 0.5
                     , agg_treatment = FALSE
                     , loop = 1){

  sample_sim_plots <- treatment(s_pool_input = s_pool_input_treatment
                                , n_sim_input = n_sim_input_treatment
                                , cv_abund_input = cv_abund_input_treatment
                                , agg = agg_treatment)

  sac <- specaccum(sample_sim_plots$spec_dat, method = "collector")
  ENS_PIE <- ens_PIE(sample_sim_plots$spec_dat)

  for(i in 2:loop){

    sample_sim_plots <- treatment(s_pool_input = s_pool_input_treatment
                                  , n_sim_input = n_sim_input_treatment
                                  , cv_abund_input = cv_abund_input_treatment
                                  , agg = agg_treatment)

    sac$richness <- sac$richness + specaccum(sample_sim_plots$spec_dat, method = "collector")$richness
    ENS_PIE <- ENS_PIE + ens_PIE(sample_sim_plots$spec_dat)


  }

  sac$richness <- sac$richness / (loop) # take average, richness already summed across iterations
  ENS_PIE <- ENS_PIE / (loop) # take average

  results <- data.frame(sac$richness, ENS_PIE)

  return(results)
}

# Let's now run our simulation! We are also going to change how many times we
# run our control and treatment. Notice the difference in time this takes to
# complete.


# Baseline - Pre-treatment
control_45 <- run_loop(loop = 10)

# Treatment - Introduce interspecies aggregation into the community
agg_45 <- run_loop(agg_treatment = TRUE, loop = 10)


# You can also try running the simulation at a much higher number of replicates
# this will take a long time though, so be warned. On my computer I estimate it
# takes about 20 minutes.
# control_45 <- run_loop(loop = 100)
# agg_45 <- run_loop(agg_treatment = TRUE, loop = 100)



###############################################
## Checkpoint 3 ###############################
###############################################
# Now to top it all off, let's turn our data into some nice figures! We will
# use the second imported file, to turn our data into nicely laid out plots.
# After importing that second function we will create a group of plots that
# investigates the effect of intraspecies aggregation. We will look at the SAC,
# ENS_PIE, and Log-Ratio Effect Size (log(control)-log(treatment)) of the
# of the treatment on both richness and the ENS_PIE.
source("./ChaseKnight_Plotter.R")

x <- 1:400


# Plot group
par(mfrow=c(2, 2))
plotter(control_45
        , agg_45
        , x = x
        , meth = "accum"
        , metric = "rich"
        , title = "Intraspecific Aggregation SAC")

plotter(control_45
        , agg_45
        , x = x
        , meth = "effect"
        , metric = "rich"
        , title = "Intraspecific Aggregation\nEffect Size on Richness")

plotter(control_45
        , agg_45
        , x = x
        , meth = "accum"
        , metric = "ENS_PIE"
        , title = "Intraspecific Aggregation ENS_PIE")

plotter(control_45
        , agg_45
        , x = x
        , meth = "effect"
        , metric = "ENS_PIE"
        , title ="Intraspecific Aggregation\nEffect Size on ENS_PIE")

###############################################
## Checkpoint 4 ###############################
###############################################
# I have commented out the following code, but it is to see the effects of
# treatments on raw richness as compared to ENS_PIE. If you are curious to 
# examine this at greater depth, uncomment the following code.

# # Treatment - Reduce Species Abundance Distribution and add rare species
# # Again, I set it to 10 replicates, but you can run it as many times as
# # you would like if you are willing to wait. This will be the case for all
# # other treatments going forward.
# medSAD_45 <- run_loop(cv_abund_input_treatment = 5, loop = 10)
# par(mfrow=c(2, 2))
# plotter(control_45
#                    , medSAD_45
#                    , x = x
#                    , meth = "accum"
#                    , metric = "rich"
#                    , title = "Reduced SAD SAC")
# 
# plotter(control_45
#                      , medSAD_45
#                      , x = x
#                      , meth = "effect"
#                      , metric = "rich"
#                      , title = "Reduced SAD\nEffect Size on Richness")
# 
# plotter(control_45
#                      , medSAD_45
#                      , x = x
#                      , meth = "accum"
#                      , metric = "ENS_PIE"
#                      , title = "Reduced SAD ENS_PIE")
# 
# plotter(control_45
#                      , medSAD_45
#                      , x = x
#                      , meth = "effect"
#                      , metric = "ENS_PIE"
#                      , title ="Reduced SAD\nEffect Size on ENS_PIE")


# # Treatment - Dramatically reduce SAD, make few common species and many rare
# lowSAD_45 <- run_loop(cv_abund_input_treatment = 10, loop = 10)
# par(mfrow=c(2, 2))
# plotter(control_45
#                    , lowSAD_45
#                    , x = x
#                    , meth = "accum"
#                    , metric = "rich"
#                    , title = "Highly Reduced SAD SAC")
# 
# plotter(control_45
#                      , lowSAD_45
#                      , x = x
#                      , meth = "effect"
#                      , metric = "rich"
#                      , title = "Highy Reduced SAD\nEffect Size on Richness")
# 
# plotter(control_45
#                      , lowSAD_45
#                      , x = x
#                      , meth = "accum"
#                      , metric = "ENS_PIE"
#                      , title = "Highly Reduced SAD ENS_PIE")
# 
# plotter(control_45
#                      , lowSAD_45
#                      , x = x
#                      , meth = "effect"
#                      , metric = "ENS_PIE"
#                      , title ="Highly Reduced SAD\nEffect Size on ENS_PIE")
# 
# # Treatment - Reduce density by halving number of individuals
# density_45 <- run_loop(n_sim_input_treatment = 16000, loop = 10)
# par(mfrow=c(2, 2))
# plotter(control_45
#         , density_45
#         , x = x
#         , meth = "accum"
#         , metric = "rich"
#         , title = "Reduced Density SAC")
# 
# plotter(control_45
#         , density_45
#         , x = x
#         , meth = "effect"
#         , metric = "rich"
#         , title = "Reduced Density\nEffect Size on Richness")
# 
# plotter(control_45
#         , density_45
#         , x = x
#         , meth = "accum"
#         , metric = "ENS_PIE"
#         , title = "Reduced Density ENS_PIE")
# 
# plotter(control_45
#         , density_45
#         , x = x
#         , meth = "effect"
#         , metric = "ENS_PIE"
#         , title ="Reduced Density\nEffect Size on ENS_PIE")


# As we approach the end of this tutorial, try creating your own treatment using
# the "run_loop()" function we wrote. Before running the "plotter()" function,
# make some predictions about the changes in the Richness and ENS_PIE curves
# and effect size on richness and ENS_PIE that will result your treatment.


# Use the next line to write your treatment.
treatment_x <- run_loop(
                     , loop = 10)

# Now write what your predictions are for the shape of the SAC, ENS_PIE and the 
# corresponding log-ratio effect size curves in the comment space below:
#
#
#
#

# Now run the following lines of code to test your predictions.
par(mfrow=c(2, 2))
plotter(control_45
        , treatment_x
        , x = x
        , meth = "accum"
        , metric = "rich"
        , title = "Treatment X SAC")

plotter(control_45
        , treatment_x
        , x = x
        , meth = "effect"
        , metric = "rich"
        , title = "Treatment X\nEffect Size on Richness")

plotter(control_45
        , treatment_x
        , x = x
        , meth = "accum"
        , metric = "ENS_PIE"
        , title = "Treatment X ENS_PIE")

plotter(control_45
        , treatment_x
        , x = x
        , meth = "effect"
        , metric = "ENS_PIE"
        , title ="Treatment X\nEffect Size on ENS_PIE")


# You have now written functions to model communities under a variety of
# treatments, and then run those simulations numerous times.
# I hope you found this tutorial both engaging and valuable. Happy coding!