library(shiny)
library(ggplot2)
library(sn)
library(numDeriv)
library(gridExtra)

shinyServer(function(input, output) {

  output$influence <- renderPlot({
    xVals = seq(input$xRng[1], input$xRng[2], length.out=input$gridPts)
    xi_1 = input$xi_1
    xi_2 = input$xi_2
    omega_11 = input$omega_11
    omega_12 = input$omega_12
    omega_22 = input$omega_22
    alpha_1 = input$alpha_1
    alpha_2 = input$alpha_2
    nu = exp(input$logNu)

    xi_1IF = sapply( xVals, function(x){try(grad(func=function(xi_1){-log(dmst(x, xi=c(xi_1,xi_2)
           ,Omega=matrix(c(omega_11, omega_12, omega_12, omega_22), ncol=2)
           ,alpha=c(alpha_1, alpha_2)
           ,nu=nu ) )}, x=xi_1 ))})
    xi_1IF = as.numeric(xi_1IF)
    xi_2IF = sapply( xVals, function(x){try(grad(func=function(xi_2){-log(dmst(x, xi=c(xi_1,xi_2)
           ,Omega=matrix(c(omega_11, omega_12, omega_12, omega_22), ncol=2)
           ,alpha=c(alpha_1, alpha_2)
           ,nu=nu ) )}, x=xi_2 ) )})
    xi_2IF = as.numeric(xi_2IF)
    omega_11IF = sapply( xVals, function(x){try(grad(func=function(omega_11){-log(dmst(x, xi=c(xi_1,xi_2)
           ,Omega=matrix(c(omega_11, omega_12, omega_12, omega_22), ncol=2)
           ,alpha=c(alpha_1, alpha_2)
           ,nu=nu ) )}, x=omega_11 ) )})
    omega_11IF = as.numeric(omega_11IF)
    omega_12IF = sapply( xVals, function(x){try(grad(func=function(omega_12){-log(dmst(x, xi=c(xi_1,xi_2)
           ,Omega=matrix(c(omega_11, omega_12, omega_12, omega_22), ncol=2)
           ,alpha=c(alpha_1, alpha_2)
           ,nu=nu ) )}, x=omega_12 ) )})
    omega_12IF = as.numeric(omega_12IF)
    omega_22IF = sapply( xVals, function(x){try(grad(func=function(omega_22){-log(dmst(x, xi=c(xi_1,xi_2)
           ,Omega=matrix(c(omega_11, omega_12, omega_12, omega_22), ncol=2)
           ,alpha=c(alpha_1, alpha_2)
           ,nu=nu ) )}, x=omega_22 ) )})
    omega_22IF = as.numeric(omega_22IF)
    alpha_1IF = sapply( xVals, function(x){try(grad(func=function(alpha_1){-log(dmst(x, xi=c(xi_1,xi_2)
           ,Omega=matrix(c(omega_11, omega_12, omega_12, omega_22), ncol=2)
           ,alpha=c(alpha_1, alpha_2)
           ,nu=nu ) )}, x=alpha_1 ) )})
    alpha_1IF = as.numeric(alpha_1IF)
    alpha_2IF = sapply( xVals, function(x){try(grad(func=function(alpha_2){-log(dmst(x, xi=c(xi_1,xi_2)
           ,Omega=matrix(c(omega_11, omega_12, omega_12, omega_22), ncol=2)
           ,alpha=c(alpha_1, alpha_2)
           ,nu=nu ) )}, x=alpha_2 ) )})
    alpha_2IF = as.numeric(alpha_2IF)
    nuIF = sapply( xVals, function(x){try(grad(func=function(nu){-log(dmst(x, xi=c(xi_1,xi_2)
           ,Omega=matrix(c(omega_11, omega_12, omega_12, omega_22), ncol=2)
           ,alpha=c(alpha_1, alpha_2)
           ,nu=nu ) )}, x=nu ))})
    nuIF = as.numeric(nuIF)
    grid.arrange( qplot(xVals, -xi_1IF, geom="line") + labs(title="xi_1", x="x", y="Influence Function")
                 ,qplot(xVals, -xi_2IF, geom="line") + labs(title="xi_2", x="x", y="Influence Function")
                 ,qplot(xVals, -omega_11IF, geom="line") + labs(title="omega_11", x="x", y="Influence Function")
                 ,qplot(xVals, -omega_12IF, geom="line") + labs(title="omega_12", x="x", y="Influence Function")
                 ,qplot(xVals, -omega_22IF, geom="line") + labs(title="omega_22", x="x", y="Influence Function")
                 ,qplot(xVals, -alpha_1IF, geom="line") + labs(title="alpha_1", x="x", y="Influence Function")
                 ,qplot(xVals, -alpha_2IF, geom="line") + labs(title="alpha_2", x="x", y="Influence Function")
                 ,qplot(xVals, -nuIF, geom="line") + labs(title="nu", x="x", y="Influence Function")
                 ,ncol=2, nrow=4
    )
  })

  output$confRegion <- renderPlot({
    xi_1 = input$xi_1
    xi_2 = input$xi_2
    omega_11 = input$omega_11
    omega_12 = input$omega_12
    omega_22 = input$omega_22
    alpha_1 = input$alpha_1
    alpha_2 = input$alpha_2
    nu = exp(input$logNu)
    dist <- makeSECdistr(dp=list(xi=c(xi_1,xi_2)
                    ,Omega=matrix(c(omega_11, omega_12, omega_12, omega_22), nrow=2)
                    ,alpha=c(alpha_1,alpha_2)
                    ,df=nu)
            ,family="ST")
    plot(dist)
  })

  output$cvf <- renderPlot({
    
  })

  output$skewed_ll <- renderPlot({
    
  })


})