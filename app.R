#' Shiny App to do PERT 3-point estimates
#' 
#' 
library("shiny")
library("shinyMatrix")
if(!exists("PoolOpinions", mode="function")) 
    source("threepoint_core.R")

#' Constant to define into how many pieces the distribution will be split. 
#' The more pieces, the more precision.
Precision <- 100000

# Define UI for app that draws a histogram ----
ui <- fluidPage(
    # App title ----
    titlePanel("Three-Point Estimation"),
    
    tabsetPanel(
        id = "tabs",
        # Single Panel  ----
        tabPanel("Single", "Three-point estimation for a single person",
                 sidebarLayout(
                     sidebarPanel(
                         
                         sliderInput(inputId = "Optimistic0",
                                     label = "Optimistic:",
                                     min = 1,
                                     max = 200,
                                     value = 15),
                         sliderInput(inputId = "Typical0",
                                     label = "Typical:",
                                     min = 1,
                                     max = 200,
                                     value = 25,
                                     step = 1),
                         sliderInput(inputId = "Pessimistic0",
                                     label = "Pessimistic:",
                                     min = 1,
                                     max = 200,
                                     value = 50)
                     ), #sidebarPanel
                     
                     # Main panel for displaying outputs ----
                     mainPanel(
                         plotOutput(outputId = "SinglePlot"),
                         textOutput('safeError')
                     )
                 ),
        ),
        # Team Panel  ----
        tabPanel("Team", "Combined team estimates",
                 sidebarLayout(
                     sidebarPanel(
                         matrixInput("TeamMatrix",
                                     rows = list(names = TRUE, editableNames = TRUE, extend = TRUE),
                                     cols = list(names = TRUE, editableNames = FALSE),
                                     value = matrix(c(11,12,8,9,13,13,10,10,17,16,13,15), 4, 3,
                                                    dimnames = list(c("Peter", "Jane", "Paul", "Susan"),
                                                                    c("Optimal", "Typical", "Pessimistic"))),
                                     #value = matrix(1:6, 2, 3,
                                     #                dimnames = list(c("MA1", "MA2"),
                                     #                                c("Opt", "Typ", "Pess"))),
                                     class = "numeric"
                         )
                     ),
                     # Main panel for displaying outputs ----
                     mainPanel(
                         plotOutput(outputId = "MultiPlot"),
                         #verbatimTextOutput("CombiStats"),
                         #textOutput(outputId = "CombiAverage"),
                         #textOutput(outputId = "CombiStdDev"),
                         #textOutput(outputId = "CombiVar"),
                         #textOutput(outputId = "CombiSkewness"),
                         #textOutput(outputId = "CombiKurtosis"),
                         plotOutput(outputId = "CombiPlot")
                     )
                 )
        )
    )#tabsetPanel
)



# Server ----
server <- function(input, output) {
    
    # Single PERT distribution ----
    # with optimistic, typical and pessimistic value
    output$SinglePlot <- renderPlot({
        shiny::validate(
            need(input$Typical0 > input$Optimistic0,     
                 "'Typical' must be larger than 'Optimistic'"),
            need(input$Pessimistic0 > input$Typical0, 
                 "'Pesimistic' must be larger than 'Typical'")
        )
        
        X <- seq(input$Optimistic0, input$Pessimistic0, length.out=Precision)
        # plot PERT distribution
        plot(X, dpert(X, min=input$Optimistic0, mode=input$Typical0, max=input$Pessimistic0),
             type = "l", 
             main = "Estimated Distribution", xlab = "Estimate", ylab = "Subjective Probability")
        
        # Add the three point estimates (opt., typical, pess.)
        points(c(input$Optimistic0, input$Typical0, input$Pessimistic0), c(0,0,0), 
               pch = 17, col = "red", cex = 3)
        
        UnitX <- range(X) / 100
        #unitY #TODO
        
        # Explicit 3-point estimate using the integral
        Average <- qpert(0.5, min=input$Optimistic0, mode=input$Typical0, max=input$Pessimistic0)
        abline(v = Average, col = "blue")
        text(Average-0.75, y = 0, labels = round(Average, 2), cex = 1, col = "blue")
        
        # Weighted 3-point estimate
        Average <- (input$Optimistic0 + 4L*input$Typical0 + input$Pessimistic0) / 6
        abline(v = Average)
        text(Average+0.5, y = 0.01, labels = round(Average, 2), cex = 1, col = "black")
    })
    
    
    # Team Data Set ----
    # Generate the data set required for the Team tab 
    MultiPertData <- reactive({
        Values <- na.omit( input$TeamMatrix )
        Estimators <- nrow(Values)
        shiny::validate(
            need(Estimators > 0, "I need at least one triple of estimates to show something")
        )
        # Values on x axis
        X <- seq(min(Values), max(Values), length.out = Precision)
        
        
        # Compute all PERT functions
        for(Line in 1:nrow(Values)) {
            dfunc <- function(x) dpert(x, min  = Values[Line, "Optimal"], 
                                       mode = Values[Line, "Typical"], 
                                       max  = Values[Line, "Pessimistic"])
            D[Line] <- AbscontDistribution(d = dfunc, withgaps = FALSE)#, withStand = TRUE)
        }
        rownames(Data) <- rownames(Values)
        
        # 
        # PooledData <- PoolOpinions(Data, Weights = Weight)
        MixedDistr <- UnivarMixingDistribution(Norm(3, .25), Norm(1, .5))
        #plot(MixedDistr)
        
        # Add 50% threshold of pooled distribution
        Threshold50 <- Area50P(X, PooledData)
        
        
        list(
            Estimates = Values,
            PertData = Data,
            PooledPertData = PooledData,
            Threshold50 = Threshold50,
            XValues = X,         # x axis
            Weights = Weight
        )
    })
    
    
    # Multiple Pert distributions ----
    # Plot a chart with several Pert distributions in it 
    output$MultiPlot <- renderPlot({
        PertData  <- MultiPertData()$PertData
        Estimates <- MultiPertData()$Estimates
        X <- MultiPertData()$Xvalues
        Weights <- MultiPertData()$Weights
        
        # Determine appropriate size of y-axis
        MaxY <- max(PertData)
        # Determine values on x-axis
        X <- seq( min(Estimates), max(Estimates), length.out = length(PertData[1,]) )
        
        
        # plot first PERT distribution
        plot(X, PertData[1,],
             type = "l", ylim = c(0, MaxY),
             main = "Individual Distributions", xlab = "Estimate", ylab = "Subjective Probability")
        
        # plot additional distributions
        if(nrow(PertData) > 1) {
            for(Line in 2:nrow(PertData))
                lines(X, PertData[Line,], col="black", lty=Line)
        }
        
        legend("topright", legend=rownames( Estimates ), lty=1:nrow(Estimates), cex=0.8)
    })
    
    
    
    # Pooled distribution ----
    # Plot a chart after distributions have been pooled together
    output$CombiPlot <- renderPlot({
        Pooled <- MultiPertData()$PooledPertData
        Estimates <- MultiPertData()$Estimates
        
        # Get values for x axis
        X <- seq( min(Estimates), max(Estimates), length.out = length(Pooled) )
        
        plot(X, Pooled,
             type = "l", 
             main = "Combined Distribution", xlab = "Estimate", ylab = "Subjective Probability")
        
        abline(v = MultiPertData()$Threshold50, col = "blue")
        text(MultiPertData()$Threshold50, y = 0, labels = round(MultiPertData()$Threshold50, 2), cex = 1, col = "blue")
    })
    
}


#
# Run the App ----
shinyApp(ui, server)