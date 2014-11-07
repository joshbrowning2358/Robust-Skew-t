library(shiny)

# Define UI for application that plots random distributions 
shinyUI(pageWithSidebar(

  # Application title
  headerPanel("Skew-t Distribution"),

  # Sidebar with a slider input for number of observations
  sidebarPanel(
     sliderInput("xi_1", "xi_1", value=0, min=-20, max=20)
    ,sliderInput("xi_2", "xi_2", value=0, min=-20, max=20)
    ,sliderInput("omega_11", "omega_11", value=1, min=0, max=20)
    ,sliderInput("omega_12", "omega_12", value=0, min=-20, max=20)
    ,sliderInput("omega_22", "omega_22", value=1, min=0, max=20)
    ,sliderInput("alpha_1", "alpha_1", value=0, min=-20, max=20)
    ,sliderInput("alpha_2", "alpha_2", value=0, min=-20, max=20)
    ,sliderInput("logNu", "log(nu)", value=1, min=0, max=10)
    ,sliderInput("k", "k", value=5, min=1, max=10)
    ,sliderInput("gridPts", "Number of grid points:", value=50, min=5, max=500)
    ,sliderInput("xRng", "Range of observation values", value=c(-10,10), min=-50, max=50)
  ),

  # Show a plot of the generated distribution
  mainPanel(
    tabsetPanel(
      tabPanel("Influence Functions", plotOutput("influence", height="600px"))
#     ,tabPanel("Confidence Regions", plotOutput("confRegion", height="600px"))
#     ,tabPanel("Change of Variance Functions", plotOutput("cvf", height="600px"))
#     ,tabPanel("Skewed Portion of Likelihood", plotOutput("skewed_ll", height="600px"))
    )
  )
))