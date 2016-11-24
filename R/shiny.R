#' @importFrom shiny fluidPage nearPoints reactive radioButtons
#'   actionButton plotOutput reactive eventReactive shinyApp titlePanel
#'   sidebarLayout sidebarPanel htmlOutput fluidRow tags mainPanel
#'   renderPrint renderTable renderPlot renderText column observe runApp
#'   stopApp wellPanel updateRadioButtons HTML numericInput sliderInput
#'   brushOpts 
NULL

#' @importFrom utils str
NULL

#' Visually assess histogram fits
#'
#' Visually assess histogram fits
#' @title browseFlowHist
#' @param flowList either a \code{\link{FlowHist}} object, or a list of
#'   \code{\link{FlowHist}} objects 
#' @param debug boolean, turns on debugging messages
#' @return Returns the list of \code{\link{FlowHist}} objects, updated by any
#'   changes made in the GUI.
#' @author Tyler Smith
#' @examples
#' library(flowPloidyData)
#' batch1 <- batchFlowHist(flowPloidyFiles, channel = "FL3.INT.LIN")
#' \dontrun{
#' batch1 <- browseFlowHist(batch1)
#' }
#' @export
browseFlowHist <- function(flowList, debug = FALSE){
  if(class(flowList) == "FlowHist"){
    flowList <- list(flowList)
    names(flowList) <- fhFile(flowList[[1]])
  }
  
  .fhI <- 1
  .fhList <- flowList

  initialLinearity <- fhLinearity(.fhList[[1]])
  if(initialLinearity == "fixed")
    initialLinearity <- "ON"
  else
    initialLinearity <- "OFF"

  if(debug) message("init Linearity: ", initialLinearity)
  
  initialDebris <- fhDebris(.fhList[[1]])

  initialSamples <- fhSamples(.fhList[[1]])
  
  if(debug) message("init Debris: ", initialDebris)

  raw <- exprs(fhRaw(.fhList[[.fhI]]))
  chan1 <- fhChannel(.fhList[[.fhI]])
  chan2 <- viewFlowChannels(.fhList[[.fhI]])[2]      # why 2? sometimes, by
                                        # chance, the 
                                        # second column is the SS value.
  
  initGateData <- data.frame(x = raw[, chan1],
                             y = raw[, chan2] / raw[, chan1]) 
  
  ui <- fluidPage(
    tags$head(
           tags$style(HTML("
      .col-sm-3 {
          max-width: 300px;
      }
    "))),
    fluidRow(
      column(width = 3,
             wellPanel(
               fluidRow(
                 column(5, htmlOutput("flowNumber", align = "center")),
                 column(3, actionButton("prev", label = "Prev")),
                 column(4, actionButton("nxt", label = "Next"))),
               tags$hr(),
               fluidRow(
                 column(4, 
                        numericInput('sampSelect', 'Samples',
                                     initialSamples, min = 2, max = 3)),
                 column(8,
                        radioButtons(inputId = "peakPicker",
                                     label = "Move peak:", 
                                     choices = list("A" = "A",
                                                    "B" = "B",
                                                    "C" = "C"), 
                                     selected = "A", inline = TRUE))),
               tags$hr(),
               radioButtons(inputId = "linearity",
                            label = "Linearity",
                            choices = list("Fixed" = "ON",
                                           "Variable" = "OFF"), 
                            inline = TRUE, selected = initialLinearity),
               tags$hr(),
               radioButtons(inputId = "debris",
                            label = "Debris Model",
                            choices = list("MC" = "MC", "SC" = "SC"),  
                            inline = TRUE, selected = initialDebris),
               tags$br(),
               actionButton("exit", label = "Return to R")
             )),
      column(width = 9,
             plotOutput("init", click = "pointPicker"))
    ),
    fluidRow(
      column(width = 3,
             wellPanel(
               selectInput('xcol', 'X Variable',
                           viewFlowChannels(.fhList[[.fhI]]), 
                         selected = chan1),
               selectInput('ycol', 'Y Variable',
                           viewFlowChannels(.fhList[[.fhI]]), 
                         selected = chan2),
               sliderInput("yrange", "Zoom", min = 0, ticks = FALSE,
                           step =
                             max(4, ceiling(log(max(initGateData$y))))/20, 
                           max = max(4, ceiling(log(max(initGateData$y)))),
                           value = 0, dragRange = FALSE))),
      column(3,
             plotOutput("gatePlot",
                        click = "gatePlot_click",
                        brush = brushOpts(
                          id = "gatePlot_brush",
                          resetOnNew = TRUE
                        )
                        )
             ),
      column(3,
             plotOutput("gatedData")
             ),
      column(3,
             plotOutput("gateResiduals")))

  )
    
  server <- function(input, output, session){
    nxtVal <- 0
    prevVal <- 0
    prefix <- ""
    
    fhPlot <- reactive({
      if(debug){
        message(prefix, "fhPlot ",
                        environmentName(environment())) 
        prefix <<- paste(prefix, " ", sep = "")}

      tmp <- .fhList[[fhCurrent()]]
      if(debug) message(prefix, "fh@linearity = ", fhLinearity(tmp))
      if(debug) message(prefix, "button value = ", input$linearity)
      xPt <- nearPoints(fhHistData(.fhList[[fhCurrent()]]),
                        input$pointPicker, "xx", "intensity",
                        threshold = 25, maxpoints = 1)
      if(debug) message(prefix, "xPt set")
      if(nrow(xPt) > 0){
        if(debug) message(prefix, "peak Picker")
        if(input$peakPicker == "A"){
          .fhList[[fhCurrent()]] <<-
            selectPeaks(.fhList[[fhCurrent()]],
                        xPt[1,1],
                        fhInit(.fhList[[fhCurrent()]])$Mb,
                        fhInit(.fhList[[fhCurrent()]])$Mc) 
          .fhList[[fhCurrent()]] <<- fhAnalyze(.fhList[[fhCurrent()]])
        } else if(input$peakPicker == "B"){
          .fhList[[fhCurrent()]] <<-
            selectPeaks(.fhList[[fhCurrent()]],
                        fhInit(.fhList[[fhCurrent()]])$Ma,
                        xPt[1,1],
                        fhInit(.fhList[[fhCurrent()]])$Mc)
          .fhList[[fhCurrent()]] <<- fhAnalyze(.fhList[[fhCurrent()]])
        } else if(input$peakPicker == "C"){
          .fhList[[fhCurrent()]] <<-
            selectPeaks(.fhList[[fhCurrent()]],
                        fhInit(.fhList[[fhCurrent()]])$Ma,
                        fhInit(.fhList[[fhCurrent()]])$Mb,
                        xPt[1,1])
          .fhList[[fhCurrent()]] <<- fhAnalyze(.fhList[[fhCurrent()]])
        }
      }

      if(debug) message(prefix, "checking linearity for element ", .fhI)
      if(input$linearity == "ON" &&
         fhLinearity(.fhList[[fhCurrent()]]) != "fixed")
      {
        if(debug) message(prefix, "fixing linearity")
        .fhList[[fhCurrent()]] <<-
          updateFlowHist(.fhList[[fhCurrent()]],
                         linearity = "fixed", analyze = TRUE)
      }
      else 
        if(input$linearity == "OFF" &&
           fhLinearity(.fhList[[fhCurrent()]]) != "variable") { 
          if(debug) message(prefix, "modeling linearity")
          .fhList[[fhCurrent()]] <<-
            updateFlowHist(.fhList[[fhCurrent()]], 
                           linearity = "variable", analyze = TRUE)
        }

      if(debug) message(prefix, "input$debris: ", input$debris)
      if(input$debris == "SC" &&
         fhDebris(.fhList[[fhCurrent()]]) != "SC")
      {
        if(debug) message(prefix, "switching to SC")
        .fhList[[fhCurrent()]] <<-
          updateFlowHist(.fhList[[fhCurrent()]],
                         debris = "SC", analyze = TRUE)
      }
      else 
        if(input$debris == "MC" &&
           fhDebris(.fhList[[fhCurrent()]]) != "MC") { 
          if(debug) message(prefix, "switching to MC")
          .fhList[[fhCurrent()]] <<-
            updateFlowHist(.fhList[[fhCurrent()]], 
                           debris = "MC", analyze = TRUE)
        }

      if(debug){
        prefix <<- substring(prefix, 2)
        message(prefix, "returning from fhPlot")
      }

      if(input$sampSelect == 2 &&
         fhSamples(.fhList[[fhCurrent()]]) != 2)
      {
        if(debug) message(prefix, "switching to 2 samples")
        .fhList[[fhCurrent()]] <<-
          updateFlowHist(.fhList[[fhCurrent()]],
                         samples = 2, analyze = TRUE)
      }
      else 
        if(input$sampSelect == 3 &&
           fhSamples(.fhList[[fhCurrent()]]) != 3) { 
          if(debug) message(prefix, "switching to 3 samples")
          .fhList[[fhCurrent()]] <<-
            updateFlowHist(.fhList[[fhCurrent()]], 
                           samples = 3, analyze = TRUE)
        }

      dat <- gateData()
      bp <- brushedPoints(dat, xvar = names(dat)[1], yvar = names(dat)[2],
                          input$gatePlot_brush, allRows = TRUE)$selected_
      if(sum(bp) > 0){
        .fhList[[fhCurrent()]] <<- setGate(.fhList[[fhCurrent()]], bp)
      }
      
      .fhList[[fhCurrent()]]
    })

    fhCurrent <- reactive({
      if(debug) {
        message(prefix, "fhCurrent, starting at ", .fhI)
        prefix <<- paste(prefix, " ", sep = "")}
      if(input$nxt > nxtVal){
        if(debug) message(prefix, "moving forwards to ", .fhI + 1)
        nxtVal <<- input$nxt
        if(.fhI < length(.fhList))
          .fhI <<- .fhI + 1
      }

      if(input$prev > prevVal){
        if(debug){
          message(prefix, "moving backwards to ", .fhI - 1)
        }
        prevVal <<- input$prev
        if(.fhI > 1)
          .fhI <<- .fhI - 1
      }
      
      if(fhLinearity(.fhList[[.fhI]]) == "fixed"){
        if(debug) message(prefix, "turning button ON")
        linVal <- "ON"
      } else {
        if(debug) message(prefix, "turning button OFF")
        linVal <- "OFF"
      }

      updateRadioButtons(session, "linearity",
                         selected = linVal)

      updateRadioButtons(session, "debris",
                         selected = fhDebris(.fhList[[.fhI]]))

      updateNumericInput(session, "sampSelect",
                         value = fhSamples(.fhList[[.fhI]])) 

      if(debug){
              prefix <<- substring(prefix, 2)
              message(prefix, "returning from fhCurrent")
      }
      .fhI
    })
    
    observe({
      if(input$exit > 0){
        stopApp()
      }
      
    })

    output$init <- renderPlot({
      if(debug){
        message(prefix, "renderPlot")
        prefix <<- paste(prefix, " ", sep = "")
      }
      plot(fhPlot(), init = TRUE, nls = TRUE, comps = TRUE)
      if(debug){
        prefix <<- substring(prefix, 2)
        message(prefix, "returning from renderPlot")
      }
    })

    output$flowNumber <- renderText({
      if(debug) message(prefix, "flowNumber")
      paste("File", tags$b(fhCurrent()), "of", length(.fhList))
    })

    observe({
      dat <- gateData()
      updateSliderInput(session, "yrange", 
                        step = max(4, ceiling(log(max(dat[,2]))))/20,
                        max = max(4, ceiling(log(max(dat[,2])))),
                        value = 0, min = 0)
    })

    ##browser()                           
    gateData <- reactive({
      chan1 <- input$xcol
      chan2 <- input$ycol

      raw <<- exprs(fhRaw(.fhList[[fhCurrent()]]))
      
      df <- data.frame(x = raw[, chan1],
                       y = raw[, chan2] / raw[, chan1])
      names(df) <- c(chan1, paste(chan1, chan2, sep = "/"))
      df
    })

    output$gatePlot <- renderPlot({
      dat <- gateData()
      plot(dat, ylim = c(0, exp(log(max(dat[, 2])) - input$yrange)),
           pch = 16, col = "#05050510") 
    })

    output$gateResiduals <- renderPlot({
      plotResid(fhHistPlot(), main = "Gate Residuals")
    })

    output$gatedData <- renderPlot({
      plot(fhHistPlot(), nls = FALSE, init = FALSE, comps = FALSE)
    })
    
    fhHistPlot <- reactive({
      dat <- gateData()
      bp <- brushedPoints(dat, xvar = names(dat)[1], yvar = names(dat)[2],
                          input$gatePlot_brush, allRows = TRUE)$selected_
      .fhList[[fhCurrent()]] <<- setGate(.fhList[[fhCurrent()]], bp)
      .fhList[[fhCurrent()]]
    })

  }
  runApp(shinyApp(ui = ui, server = server))
  return(.fhList)
}

selectPeaks <- function(fh, peakA, peakB, peakC){
  pA <- fhHistData(fh)[round(peakA, 0), c("xx", "intensity")]
  if(is.numeric(peakB))                 
    pB <- fhHistData(fh)[round(peakB, 0), c("xx", "intensity")]
  if(is.numeric(peakC))                 
    pC <- fhHistData(fh)[round(peakC, 0), c("xx", "intensity")]
  
  fh <- resetFlowHist(fh)

  if(is.numeric(peakC))
    newPeaks <- as.matrix(rbind(pA, pB, pC))
  else if(is.numeric(peakB))
    newPeaks <- as.matrix(rbind(pA, pB))
  else
    newPeaks <- as.matrix(rbind(pA))
  
  colnames(newPeaks) <- c("mean", "height")

  fhPeaks(fh) <- newPeaks
  
  fh <- addComponents(fh)
  fh <- setLimits(fh)
  fh <- makeModel(fh)
  fh <- getInit(fh)

  return(fh)
}
