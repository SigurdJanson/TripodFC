#' Shiny App to do PERT 3-point estimates
#' 
#' 
library(shiny)
library(shinyMatrix)
library(ggplot2)
if(!exists("PoolOpinions", mode="function")) 
    source("threepoint_core.R")

#' Constant to define into how many pieces the distribution will be split. 
#' The more pieces, the more precision.
Precision <- 50000

# Define UI for app that draws a histogram ----
ui <- fluidPage(
    titlePanel("Tripod FC"),
    
    tabsetPanel(
        id = "tabs",
        # Single Panel  ----
        tabPanel("Single", h3("Three-Point Estimate of a Single Expert"),
                 sidebarLayout(
                     sidebarPanel(
                         sliderInput(inputId = "Optimistic0",
                                     label = "Optimistic",
                                     min = 1, max = 200,
                                     value = 15),
                         sliderInput(inputId = "Typical0",
                                     label = "Typical",
                                     min = 1, max = 200,
                                     value = 25,
                                     step = 1),
                         sliderInput(inputId = "Pessimistic0",
                                     label = "Pessimistic",
                                     min = 1, max = 200,
                                     value = 50),
                         wellPanel(
                             h4("Distribution Stats"),
                             tableOutput("SinglePertStats")
                         )
                     ), #sidebarPanel
                     
                     # Main panel for displaying outputs ----
                     mainPanel(
                         plotOutput(outputId = "SinglePlot")
                     )
                 ),
        ),
        # Team Panel  ----
        tabPanel(
            "Team", h3("Combined Team Estimates"),
             sidebarLayout(
                 sidebarPanel(
                     matrixInput("TeamMatrix",
                                 rows = list(names = TRUE, editableNames = TRUE, extend = TRUE),
                                 cols = list(names = TRUE, editableNames = FALSE),
                                 value = matrix(c(11,12,8,9,13,13,10,10,17,16,13,15), 4, 3,
                                                dimnames = list(c("Peter", "Jane", "Paul", "Susan"),
                                                                c("Optimal", "Typical", "Pessimistic"))),
                                 class = "numeric"
                     ),
                     wellPanel(
                         h4("Stats for Mixed Distribution"),
                         tableOutput("CombiPertStats")
                     )
                 ),
                 # Main panel for displaying outputs ----
                 mainPanel(
                     plotOutput(outputId = "MultiPlot"),
                     plotOutput(outputId = "CombiPlot")
                 )
             )
        )
    )#tabsetPanel
)



# Server ----
server <- function(input, output) {
    
    # Single PERT distribution ----
    
    # STATS
    output$SinglePertStats <- renderTable({
        req(input$Typical0 > input$Optimistic0, 
            input$Pessimistic0 > input$Typical0)
        
        Mean <- (input$Optimistic0 + 4*input$Typical0 + input$Pessimistic0) / 6
        # Median
        Median <- qBetaPert(0.5, min=input$Optimistic0, mode=input$Typical0, max=input$Pessimistic0)
        Var    <- (Mean - input$Optimistic0) * (input$Pessimistic0 - Mean) / 7
        StdDev <- sqrt(Var)
        SixthSigma <- (input$Pessimistic0 - input$Optimistic0) / 6
        
        data.frame(
            Parameter = c("Mean", "Median", "Variance", "StdDev", "1/6th sigma"),
            Value     = c(Mean, Median, Var, StdDev, SixthSigma)
        )
    })
    
    
    # PLOT
    output$SinglePlot <- renderPlot({
        # PRECONDITIONS
        shiny::validate(
            need(input$Typical0 > input$Optimistic0,     
                 "'Typical' must be larger than 'Optimistic'"),
            need(input$Pessimistic0 > input$Typical0, 
                 "'Pesimistic' must be larger than 'Typical'")
        )
        
        X <- seq(input$Optimistic0, input$Pessimistic0, length.out=Precision)
        Data <- data.frame(
            X = X, 
            Y = dBetaPert(X, min=input$Optimistic0, mode=input$Typical0, max=input$Pessimistic0)
        )
        
        # plot PERT distribution
        p <- ggplot(Data, aes(x = X, y = Y)) +
                geom_line() +
                annotate("point", x = input$Optimistic0, y = 0, size = 5, colour = "blue", alpha = 1/5) +
                annotate("point", x = input$Typical0, y = max(Data$Y), size = 5, colour = "blue", alpha = 1/5) +
                annotate("point", x = input$Pessimistic0, y = 0, size = 5, colour = "blue", alpha = 1/5) +
                ggtitle("Estimated Distribution") +
                scale_x_continuous(name="Estimate") + scale_y_continuous(name="Subjective Probability")
        
        #-print(paste("max(Data$Y)", max(Data$Y), " - X: ", Data$X[which.max(Data$Y)]))

        # Explicit 3-point estimate using the integral
        Avg50p <- qBetaPert(0.5, min=input$Optimistic0, mode=input$Typical0, max=input$Pessimistic0)
        AvgMean <- (input$Optimistic0 + 4L*input$Typical0 + input$Pessimistic0) / 6

        # Plot median, i.e. Avg50p
        Alignment <- ifelse(Avg50p >= AvgMean, "left", "right")
        p <- p + annotate("point", x = Avg50p, y = 0, size = 5, colour = "red", alpha = 1/2) +
                 geom_vline(aes(xintercept = Avg50p), color = "red", alpha = 1/5) +
                 annotate("label", hjust = Alignment, colour = "red", alpha = 5/10,
                          x = Avg50p, y = max(Data$Y)/20,
                          label = format(Avg50p, digits = 2L, nsmall = 2L))
        
        # Plot weighted 3-point estimate, i.e. mean
        Alignment <- ifelse(Avg50p < AvgMean, "left", "right")
        p <- p + annotate("point", x = AvgMean, y = 0, size = 5, colour = "green", alpha = 1/2) +
                 geom_vline(aes(xintercept = AvgMean), color = "green", alpha = 1/5) +
                 annotate("label", hjust = Alignment, colour = "green", alpha = 5/10,
                          x = AvgMean, y = max(Data$Y)/20, 
                          label = format(AvgMean, digits = 2L, nsmall = 2L))
        
        p
    })
    
    
    # Team Distributions ----
    MultiPertValues <- reactive({
        Values <- na.omit( input$TeamMatrix )
        Estimators <- nrow(Values)
        if (Estimators == 0)
            return(NULL)
        else
            return(Values)
    })
    
    
    # Generate the data set required for the Team tab 
    MultiPertData <- reactive({
        Values <- req(MultiPertValues())
        
        X = seq(min(Values), max(Values), length.out=Precision)
        Data <- data.frame()
        
        for (Line in 1:nrow(Values)) {
            Y <- dBetaPert(X, min = Values[Line, 1], mode=Values[Line, 2], max=Values[Line, 3])
            GroupName <- ifelse(is.null(rownames(Values)[Line]), Line, rownames(Values)[Line])
            Group <- rep(rownames(Values)[Line], length(X))
            Data  <- rbind(Data, data.frame(X = X, Y = Y, Group = Group))
        }
        Data$Group <- factor(Data$Group)
        return(Data)
    })
    
    # Plot a chart with several Pert distributions in it 
    output$MultiPlot <- renderPlot({
        Data   <- MultiPertData()
        Values <- MultiPertValues()
        
        shiny::validate(
          need(Values > 0, "I need at least one triple of estimates to show something"),
          need(Data, "No data found")
        )
        
        # plot combined PERT distribution
        p <- ggplot(Data, aes(x = X, y = Y, colour = Group)) +
                    geom_line() +
                    ggtitle("Individual Distributions") +
                    scale_x_continuous(name="Estimate") + 
                    scale_y_continuous(name="Subjective Probability") +
                    theme(legend.position="top") + 
                    scale_colour_discrete(name="Expert")

        for (Line in 1:nrow(Values)) {
            Name <- rownames(Values)[Line]
            PosX <- Values[Line, 2] + ((Values[Line, 3] - Values[Line, 1]) * 0.1 )
            PosY <- dBetaPert(Values[Line, 2], min = Values[Line, 1], mode=Values[Line, 2], max=Values[Line, 3])
            p <- p + annotate("label", x = PosX, y = PosY, label = Name, alpha = 3/10)
        }
        
        p
    })
    
    
    
    # Combined = Pooled distribution ----
    CombinedPertData <- reactive({
        Values <- MultiPertValues()
        Data   <- MultiPertData()
        X <- seq(min(Data$X), max(Data$X), length.out = nrow(Data) / nlevels(Data$Group))
        Y <- rowsum(Data$Y, Data$X) / nlevels(Data$Group)
        Data <- data.frame(X = X, Y = Y) 
        
        return(Data)
    })
    
    
    output$CombiPertStats <- renderTable({
        Data <- CombinedPertData()
        
        Mean   <- AvgOfDensities(Data$X, Data$Y)
        Median <- Area50P(Data$X, Data$Y)
        Var    <- VarOfDensities(Data$X, Data$Y)
        StdDev <- sqrt(Var)
        SixthSigma <- (max(Data$X) - min(Data$X)) / 6
        
        data.frame(
            Parameter = c("Mean", "Median", "Variance", "StdDev", "1/6th sigma"),
            Value     = c(Mean, Median, Var, StdDev, SixthSigma)
        )
    })
    
    
    
    # Plot a chart after distributions have been pooled together
    output$CombiPlot <- renderPlot({
        Values <- req(MultiPertValues())
        Data   <- req(CombinedPertData())
        
        # Compute 50% "integral", i.e. quantile
        Avg50p  <- Area50P(Data$X, Data$Y)
        # Mean
        AvgMean <- AvgOfDensities(Data$X, Data$Y)
        # Var
        Variance <- VarOfDensities(Data$X, Data$Y)
        
        # plot PERT distribution
        p <- ggplot(Data, aes(x = X, y = Y)) +
            geom_line() +
            annotate("point", x = min(Values), y = 0, size = 5, colour = "blue", alpha = 1/5) +
            annotate("point", x = max(Values), y = 0, size = 5, colour = "blue", alpha = 1/5) +
            ggtitle("Mixed Distribution") +
            scale_x_continuous(name="Estimate") + scale_y_continuous(name="Subjective Probability")
        
        # Add median marker, i.e. Area50p
        Alignment <- ifelse(Avg50p > AvgMean, "left", "right")
        p <- p +
            annotate("point", x = Avg50p, y = 0, size = 5, colour = "red", alpha = 1/5) +
            annotate("label", x = Avg50p, y = max(Data$Y)/20, hjust = Alignment, 
                     size = 5, colour = "red", alpha = 1/5, 
                     label = format(Avg50p, digits = 2L, nsmall = 2L))
            
        Alignment <- ifelse(Avg50p <= AvgMean, "left", "right")
        p <- p +
            annotate("point", x = AvgMean, y = 0, size = 5, colour = "green", alpha = 1/5) +
            annotate("label", x = AvgMean, y = max(Data$Y)/20, hjust = Alignment,
                     size = 5, colour = "green", alpha = 1/5, 
                     label = format(AvgMean, digits = 2L, nsmall = 2L))
            
        return(p)
    })
}


#
# Run the App ----
shinyApp(ui, server)