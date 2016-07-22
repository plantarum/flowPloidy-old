ui <- fluidPage(
  radioButtons(inputId = "peakPicker", label = "Move peak:",
               choices = list("A" = "A", "B" = "B"), selected = "A"), 
  actionButton(inputId = "update", label = "Update"),
  plotOutput("init", click = "pointPicker"),
  plotOutput("nls")
)

server <- function(input, output){
  fhReact <- reactive({
    xPt <- nearPoints(.fh$data, input$pointPicker, "x", "intensity",
                      threshold = 25, maxpoints = 1)
    if(nrow(xPt) > 0){
      if(input$peakPicker == "A"){
        .fh <<- sliderPeaks(.fh, xPt[1,1], .fh$init$Mb)
      } else {
        .fh <<- sliderPeaks(.fh, .fh$init$Ma, xPt[1,1])
      }
    }
    .fh
  })

  fhNLS <- eventReactive(input$update, {
    .fh <<- fhAnalyze(.fh)
  })

  output$plot_clickedpoints <- renderTable({
    res <- nearPoints(.fh$data, input$pointPicker, "x", "intensity",
                      threshold = 25, maxpoints = 1)
    if (nrow(res) == 0)
      return()
    res
  })
  
  output$click_info <- renderPrint({
    cat("input$pointPicker:\n")
    str(input$pointPicker)
  })
  
  output$init <- renderPlot({
    plot(fhReact(), init = TRUE)
  })

  output$nls <- renderPlot({
    plot(fhNLS())
  })

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
