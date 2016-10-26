#' @importFrom shiny selectInput

plotGate <- function(fh, x, y, ...){
  chans <- viewFlowChannels(fh)
  if(! x %in% chans)
    stop("Selected flow channel: ", x, " does not exist!")
  if(! y %in% chans)
      stop("Selected flow channel: ", y, " does not exist!")

  dat <- exprs(fhRaw(fh))[, c(x, y)]
  dat[,y] <- dat[,y] / dat[,x]
  plot(dat, pch = 16, col = "#05050510")
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

isGated <- function(fh){
  ## returns TRUE if the FlowHist data is gated
  sum(fhGate(fh)) != 0
}

gateFlowHist <- function(fh){
  raw <- exprs(fhRaw(fh))
  chan1 <- viewFlowChannels(fh)[1]
  chan2 <- viewFlowChannels(fh)[2]
  
  initGateData <- data.frame(x = raw[, chan1],
                             y = raw[, chan2] / raw[, chan1]) 

  ui <- fluidPage(
    fluidRow(
      column(width = 2,
             selectInput('xcol', 'X Variable', viewFlowChannels(fh),
                         selected = chan1),
             selectInput('ycol', 'Y Variable', viewFlowChannels(fh),
                         selected = chan2),
             sliderInput("yrange", "yrange", min = 0,
                         max = max(initGateData$y),
                         value = c(0, max(initGateData$y)),
                         step = NULL, round = FALSE,
                         ticks = TRUE, animate = FALSE,
                         width = NULL, sep = ",", pre = NULL, post = NULL,
                         dragRange = FALSE)
             ),
      
      column(width = 4,
             plotOutput("gatePlot", height = 300,
                        clickOpts(id = "plot_click"),
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
    gateData <- reactive({
      chan1 <- viewFlowChannels(fh)[input$xcol]
      chan2 <- viewFlowChannels(fh)[input$ycol]
      
      data.frame(x = raw[, chan1],
                 y = raw[, chan2] / raw[, chan1])
    })

    output$gatePlot <- renderPlot({
      plotGate(gateData(), ylim = input$yrange)
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
      nearPoints(data.frame(exprs(fhRaw(fh))), xvar = input$xcol,
                 yvar = input$ycol, input$gatePlot_click, addDist = TRUE)
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
