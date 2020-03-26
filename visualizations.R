# Meta-Analitic Visualizations 

apatheme=theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line(),
        text=element_text(size=20),
        legend.position='bottom')

# Forest Plot
plot_Forest <- function(REModel){
  # Decrease margins so the full space is used
  par(mar=c(2.5,4,1,2.5), cex = .9, font = 1)
  
  # Create the forest plot
  # We need to costumize our table depending on the effect sizes and the thigs we want to investigate
  forest(REModel, 
         xlim = c(-2.5, 1.8),                 # adjust horizontal plot region limits
         order ="obs",                        # order by effect size
         addfit = T,                          # add standard summary polygon
         annotate = T,                        # remove annotations
         width = 0,                           # width of annotations
         efac = .55,                          # height of CI bars
         pch = 19,                            # changing point symbol to filled circle
         col = "gray40",                      # change color of points/CIs
         clim = c(-1 ,1),                     # set absolute limits for CIs
         cex.lab = 1,                         # increase size of x-axis title
         cex.axis = 1,                        # increase size of x-axis labels
         cex = .85,                           # set font expansion factor
         lty = c("solid",                     # CI line type
                 "solid",                     # credibility interval line type
                 "solid"),                    # horizontal line type
         xlab = "",                           # label X axis
         mlab = "RE Model :ST,  p=.077, I^2=20.99",    # label summary estimate
         showweights=F,                       # include weights given to effects in model fitting
         steps = 5)                             # number of tick marks for the x-axis
  
  # Switch to bold font
  par(cex = .9, font = 2)
}