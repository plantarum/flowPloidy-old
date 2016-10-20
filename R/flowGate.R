gatePlot <- function(fh, x, y, ...){
  chans <- viewFlowChannels(fh)
  if(! x %in% chans)
    stop("Selected flow channel: ", x, " does not exist!")
  if(! y %in% chans)
      stop("Selected flow channel: ", y, " does not exist!")

  plot(exprs(fhRaw(fh))[, c(x, y)], pch = 16, col = "#05050510")
}

setGate <- function(fh, gate){
  fhGate(fh) <- gate
  fh <- setBins(fh, fhBins(fh1))
  fh
}

gateFlowHist <- function(fh){
  x <- "SS.INT.LIN"
  y <- "FS.TOF.LIN"
  
  ui <- fluidPage(
    fluidRow(
      column(width = 2,
             selectInput('xcol', 'X Variable', viewFlowChannels(fh)),
             selectInput('ycol', 'Y Variable', viewFlowChannels(fh))),
      column(width = 4,
             plotOutput("plot1", height = 300,
                                        # Equivalent to: click = clickOpts(id = "plot_click")
                        click = "plot1_click",
                        brush = brushOpts(
                          id = "plot1_brush"
                        )
                        )
             ),
      column(width = 4,
             plotOutput("plot2", height = 300))
    ),
    fluidRow(
      column(width = 6,
             h4("Points near click"),
             verbatimTextOutput("click_info")
             )
    )
  )

  server <- function(input, output) {
    ##browser()                           
    output$plot1 <- renderPlot({
      message(length(fhGate(fh)))
      message(sum(fhGate(fh)))
      gatePlot(fh, input$xcol, input$ycol)
    })

    output$plot2 <- renderPlot({
      plot(fhHistPlot(), nls = FALSE, init = FALSE, comps = FALSE)
    })
    
    output$click_info <- renderPrint({
      nearPoints(data.frame(exprs(fhRaw(fh))), xvar = input$xcol, yvar =
      input$ycol, input$plot1_click, addDist = TRUE) 
    })

    fhHistPlot <- reactive({
      bp <- brushedPoints(data.frame(exprs(fhRaw(fh))),
                                   xvar = input$xcol, yvar = input$ycol,
                                   input$plot1_brush,
                                   allRows = TRUE)$selected_
      fh <<- setGate(fh, bp)
      fh
    })
    
  }
  runApp(shinyApp(ui = ui, server = server))
}
