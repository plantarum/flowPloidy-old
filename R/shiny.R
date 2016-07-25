#' @importFrom shiny fluidPage nearPoints reactive radioButtons actionButton plotOutput reactive eventReactive shinyApp

flowShiny <- function(flowList){
  .fhI <- 1
  .fhList <- flowList
  runApp(shinyApp(ui = ui, server = server))
  return(.fhList)
}

ui <- fluidPage(
  titlePanel("flowPloidy"),
  sidebarLayout(
    sidebarPanel(
      htmlOutput("flowNumber", align = "center"),
      fluidRow(
        column(6, actionButton("nxt", label = "Next")),
        column(6, actionButton("prev", label = "Previous"))),
      tags$hr(),
      radioButtons(inputId = "peakPicker",
                   label = "Move peak:", 
                   choices = list("A" = "A", "B" = "B"),
                   selected = "A", inline = TRUE),
      tags$hr(),
      actionButton("exit", label = "Return to R")),
    mainPanel(
      plotOutput("init", click = "pointPicker")))
  ## plotOutput("nls")
)

server <- function(input, output){
  nxtVal <- 0
  prevVal <- 0
  updateVal <- 0
  
  fhInitPlot <- reactive({
    xPt <- nearPoints(.fhList[[fhCurrent()]]$data, input$pointPicker, "x",
                      "intensity", threshold = 25, maxpoints = 1)
    if(nrow(xPt) > 0){
      if(input$peakPicker == "A"){
        .fhList[[fhCurrent()]] <<- sliderPeaks(.fhList[[fhCurrent()]],
                                        xPt[1,1], 
                                        .fhList[[fhCurrent()]]$init$Mb) 
        .fhList[[fhCurrent()]] <<- fhAnalyze(.fhList[[fhCurrent()]])
      } else {
        .fhList[[fhCurrent()]] <<- sliderPeaks(.fhList[[fhCurrent()]],
                                        .fhList[[fhCurrent()]]$init$Ma,
                                        xPt[1,1])
        .fhList[[fhCurrent()]] <<- fhAnalyze(.fhList[[fhCurrent()]])
      }
    }
    .fhList[[fhCurrent()]]
  })

  fhCurrent <- reactive({
    if(input$nxt > nxtVal){
      nxtVal <<- input$nxt
      if(.fhI < length(.fhList))
        .fhI <<- .fhI + 1
    }

    if(input$prev > prevVal){
      prevVal <<- input$prev      
      if(.fhI > 1)
        .fhI <<- .fhI - 1
    }

    ## if(input$update > updateVal){
    ##   updateVal <<- input$update      
    ##   .fhList[[fhCurrent()]] <<- fhAnalyze(.fhList[[fhCurrent()]])
    ## }

    .fhI
  })
  
  observe({
    if(input$exit > 0)
      stopApp()
  })

  fhNLS <- eventReactive(input$update, {
    .fhList[[fhCurrent()]] <<- fhAnalyze(.fhList[[fhCurrent()]])
  })

  output$plot_clickedpoints <- renderTable({
    res <- nearPoints(.fhList[[fhCurrent()]]$data, input$pointPicker, "x",
                      "intensity", threshold = 25, maxpoints = 1)
    if (nrow(res) == 0)
      return()
    res
  })
  
  output$click_info <- renderPrint({
    cat("input$pointPicker:\n")
    str(input$pointPicker)
  })
  
  output$init <- renderPlot({
    plot(fhInitPlot(), init = TRUE, nls = TRUE, comps = TRUE)
  })

  output$flowNumber <- renderText({
    paste("File", tags$b(fhCurrent()), "of", length(.fhList))
  })
  ## output$nls <- renderPlot({
  ##   plot(fhNLS())
  ## })

}
    

## shinyApp(ui = ui, server = server)     
    
sliderPeaks <- function(fh, peakA, peakB){
  pA <- fh$data[round(peakA, 0), ]
  pB <- fh$data[round(peakB, 0), ]

  fh$init$Ma <- pA[, "x"]
  fh$init$Sa <- fh$init$Ma / 20
  fh$init$a1 <- pA[, "intensity"] * fh$init$Sa / 0.45

  if((pA[, "x"] * 2) <= max(fh$data[ ,"x"])){
    fh$init$a2 <- fh$data[fh$init$Ma * 2, "intensity"] *
      fh$init$Sa * 2 / 0.4
    fh$comps$fA2 <- fA2
  } else {
    fh$init$a2 <- NULL
    fh$comps$fA2 <- NULL
  }
  
  fh$init$Mb <- pB[, "x"]
  fh$init$Sb <- fh$init$Mb / 20
  fh$init$b1 <- pB[, "intensity"] * fh$init$Sb / 0.45

  if((pB[, "x"] * 2) <= max(fh$data[ ,"x"])){
    fh$init$b2 <- fh$data[fh$init$Mb * 2, "intensity"] *
      fh$init$Sb * 2 / 0.4
    fh$comps$fB2 <- fB2
  } else {
    fh$init$b2 <- NULL
    fh$comps$fB2 <- NULL
  }

  fh <- makeModel(fh)

  fh$nls <- NULL
  return(fh)
}
