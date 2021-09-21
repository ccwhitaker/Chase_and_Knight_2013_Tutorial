# Chase_and_Knight_2013_Tutorial
A tutorial I wrote for my "Frontiers in Biodiversity Measurement" course, based on the 2013 paper “Scale-dependent effect sizes of ecological drivers on biodiversity: why standardised sampling is not enough” by Jonathan Chase and Tiffany Knight.

“RTutorial_Knight_and_Chase_Main.R” is the main file for the tutorial.
 
The tutorial as it is presented will first guide the user through a “toy” simulation so that you get a basic feel for how the functions work and to ease in to exploring the topics presented in the paper. This involves:

    • Generating two simulated community data sets, and assigning them spatial coordinates on a 2-dimensional plane
    • Dividing the communities into quadrats and obtaining a census
    • We conclude the toy simulation with an examination of the two communities’ respective SAC’s

Then we build on this with a couple of functions that make running replicates of the simulation, and adjustments to the treatments, easy. I found that I needed some backend code to make the tutorial more closely match my vision, but also felt that that code would detract from the final presentation. I placed these accessory functions in separate files ("./ChaseKnight_Cumulative_ENS_PIE.R", and "./ChaseKnight_Plotter.R"). They have appropriate documentation if you are curious to investigate them, but they are not meant to be consumed as part of the tutorial.

Similarly, I had originally intended to replicate all four treatments that Chase and Knight had included in their paper. However, that proved quite ambitious with the time-limited nature of the tutorial and I decided it would best illustrate the core concept of the paper to focus on aggregated vs non-aggregated treatments (ie. a community treatment with highly scale dependent effects). The code to replicate these treatments was excluded from the tutorial in class due to time constraints, but I have re-added it and at the end. Following this, I invite the user to test different treatments and predict a result.

My goal was to present how the ENSPIE metric reflects differences in biodiversity between communities with respect to scale. You can see this well enough with the comparison between the aggregated community and the control. ENS_PIE is clearly responding to a treatment that is fundamentally scale dependent. To really get a comprehensive picture though, you have to compare different treatments, which as mentioned before was cut from the class tutorial for time constraints. I have gone ahead and included the code to run through the other treatments which will allow you to gain the insight from this comparison.
