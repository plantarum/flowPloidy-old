#' @importFrom shiny selectInput

plotGate <- function(fh, x, y, ...){
  chans <- viewFlowChannels(fh)
  if(! x %in% chans)
    stop("Selected flow channel: ", x, " does not exist!")
  if(! y %in% chans)
      stop("Selected flow channel: ", y, " does not exist!")

  plot(exprs(fhRaw(fh))[, c(x, y)], pch = 16, col = "#05050510")
}

plotResid <- function(fh, ...){
  plot(fhHistData(fh)$gateResid, type = 'n', main = fhFile(fh),
       sub = "Gate Resdiduals", xlab = fhChannel(fh), ...)
  polygon(x = c(fhHistData(fh)$xx, max(fhHistData(fh)$xx) + 1),
          y = c(fhHistData(fh)$gateResid, 0),
          col = "lightgray", border = NA)
}

setGate <- function(fh, gate){
  fhGate(fh) <- gate
  fh <- setBins(fh, fhBins(fh))
  fh
}

gateFlowHist <- function(fh){
  x <- "SS.INT.LIN"
  y <- "FS.TOF.LIN"
  
  ui <- fluidPage(
    fluidRow(
      column(width = 2,
             selectInput('xcol', 'X Variable', viewFlowChannels(fh),
                         selected = viewFlowChannels(fh)[1]),
             selectInput('ycol', 'Y Variable', viewFlowChannels(fh),
                         selected = viewFlowChannels(fh)[2])),
      column(width = 4,
             plotOutput("gatePlot", height = 300,
                                        # Equivalent to: click = clickOpts(id = "plot_click")
                        click = "gatePlot_click",
                        brush = brushOpts(
                          id = "gatePlot_brush"
                        )
                        )
             ),
      column(width = 4,
             plotOutput("gatedData", height = 300))
    ),
    fluidRow(
      column(width = 6,
             h4("Points near click"),
             verbatimTextOutput("click_info")
             ),
      column(width = 6,
             plotOutput("gateResiduals"))
    )
  )

  server <- function(input, output) {
    ##browser()                           
    output$gatePlot <- renderPlot({
      message(length(fhGate(fh)))
      message(sum(fhGate(fh)))
      plotGate(fh, input$xcol, input$ycol)
    })

    output$gateResiduals <- renderPlot({
      message("Residual plot: ", length(fhGate(fh)))
      message(sum(fhGate(fh)))
      tmp <- input$xcol                 # dummy to trigger redraw
      tmp <- input$gatePlot_brush
      plotResid(fh)
    })

    output$gatedData <- renderPlot({
      plot(fhHistPlot(), nls = FALSE, init = FALSE, comps = FALSE)
    })
    
    output$click_info <- renderPrint({
      nearPoints(data.frame(exprs(fhRaw(fh))), xvar = input$xcol, yvar =
      input$ycol, input$gatePlot_click, addDist = TRUE) 
    })

    fhHistPlot <- reactive({
      bp <- brushedPoints(data.frame(exprs(fhRaw(fh))),
                                   xvar = input$xcol, yvar = input$ycol,
                                   input$gatePlot_brush,
                                   allRows = TRUE)$selected_
      fh <<- setGate(fh, bp)
      fh
    })
    
  }
  runApp(shinyApp(ui = ui, server = server))
}
