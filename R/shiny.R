#' @importFrom shiny fluidPage nearPoints reactive radioButtons
#'   actionButton plotOutput reactive eventReactive shinyApp titlePanel
#'   sidebarLayout sidebarPanel htmlOutput fluidRow tags mainPanel
#'   renderPrint renderTable renderPlot renderText column observe runApp
#'   stopApp wellPanel
NULL

#' @importFrom utils str
NULL

#' Visually assess histogram fits
#'
#' Visually assess histogram fits
#' @title browseFlowHist
#' @param flowList a list of \code{FlowHist} objects
#' @return Returns the list of \code{FlowHist} objects, updated by any
#'   changes made in the GUI.
#' @author Tyler Smith
#' @examples
#' library(flowPloidyData)
#' batch1 <- batchFlowHist(flowPloidyFiles, channel = "FL3.INT.LIN")
#' \dontrun{
#' batch1 <- browseFlowHist(batch1)
#' }
#' @export
browseFlowHist <- function(flowList){
  .fhI <- 1
  .fhList <- flowList

  initialLinearity <- .fhList[[1]]@linearity
  if(initialLinearity == "fixed")
    initialLinearity <- "ON"
  else
    initialLinearity <- "OFF"

  initialDebris <- .fhList[[1]]@opts
  if("SC" %in% initialDebris)
    initialDebris <- "SC"
  else
    initialDebris <- "MC"

  ui <- fluidPage(
    ## sidebarLayout(
      fluidRow(
        column(3,
               wellPanel(
               htmlOutput("flowNumber", align = "center"),
               fluidRow(
                 column(5, actionButton("prev", label = "Previous")),
                 column(5, actionButton("nxt", label = "Next"))))),
        ## tags$hr(),
        column(5,
               fluidRow(
                 column(4, 
                        wellPanel(radioButtons(inputId = "peakPicker",
                                               label = "Move peak:", 
                                               choices = list("A" = "A",
                                                              "B" = "B"), 
                                               selected = "A", inline =
                                                                 TRUE))), 
                 column(4,
                        wellPanel(radioButtons(inputId = "linearity",
                                               label = "Fixed linearity?",
                                               choices =
                                                 list("ON" = "ON",
                                                      "OFF" = "OFF"), 
                                               inline = TRUE,
                                               selected =
                                                 initialLinearity))),
                 column(4,
                        wellPanel(radioButtons(inputId = "debris",
                                               label = "Debris Model",
                                               choices =
                                                 list("MC" = "MC", "SC" =
                                                                     "SC"),  
                                               inline = TRUE,
                                               selected = initialDebris))))),
        ## tags$hr(),
        column(2, actionButton("exit", label = "Return to R"))),
    fluidRow(
      plotOutput("init", click = "pointPicker"))
  )

  server <- function(input, output, session){
    nxtVal <- 0
    prevVal <- 0
    prefix <- ""
    
    fhInitPlot <- reactive({
      message(prefix, "fhInitPlot ", environmentName(environment()))
      prefix <<- paste(prefix, " ", sep = "")
      tmp <- .fhList[[fhCurrent()]]
      message(prefix, "fh@linearity = ", tmp@linearity)
      message(prefix, "button value = ", input$linearity)
      xPt <- nearPoints(.fhList[[fhCurrent()]]@histData,
                        input$pointPicker, "xx", "intensity",
                        threshold = 25, maxpoints = 1)
      message(prefix, "xPt set")
      if(nrow(xPt) > 0){
        message(prefix, "peak Picker")
        if(input$peakPicker == "A"){
          .fhList[[fhCurrent()]] <<-
            selectPeaks(.fhList[[fhCurrent()]], xPt[1,1],
                        .fhList[[fhCurrent()]]@init$Mb) 
          .fhList[[fhCurrent()]] <<- fhAnalyze(.fhList[[fhCurrent()]])
        } else {
          .fhList[[fhCurrent()]] <<-
            selectPeaks(.fhList[[fhCurrent()]],
                        .fhList[[fhCurrent()]]@init$Ma, xPt[1,1])
          .fhList[[fhCurrent()]] <<- fhAnalyze(.fhList[[fhCurrent()]])
        }
      }

      message(prefix, "checking linearity for element ", .fhI)
      if(input$linearity == "ON" &&
         .fhList[[fhCurrent()]]@linearity != "fixed")
      {
        message(prefix, "fixing linearity")
        .fhList[[fhCurrent()]] <<-
          updateFlowHist(.fhList[[fhCurrent()]],
                         linearity = "fixed", analyze = TRUE)
      }
      else 
        if(input$linearity == "OFF" &&
           .fhList[[fhCurrent()]]@linearity != "variable") { 
          message(prefix, "modeling linearity")
          .fhList[[fhCurrent()]] <<-
            updateFlowHist(.fhList[[fhCurrent()]], 
                           linearity = "variable", analyze = TRUE)
        }
      prefix <<- substring(prefix, 2)
      message(prefix, "returning from fhInitPlot")

      if(input$debris == "SC" &&
         ! "SC" %in% .fhList[[fhCurrent()]]@opts)
      {
        message(prefix, "switchint to SC")
        .fhList[[fhCurrent()]] <<-
          updateFlowHist(.fhList[[fhCurrent()]],
                         opts = list("SC"), analyze = TRUE)
      }
      else 
        if(input$debris == "MC" &&
           "SC" %in% .fhList[[fhCurrent()]]@opts) { 
          message(prefix, "switching to MC")
          .fhList[[fhCurrent()]] <<-
            updateFlowHist(.fhList[[fhCurrent()]], 
                           opts = list(), analyze = TRUE)
        }
      prefix <<- substring(prefix, 2)
      message(prefix, "returning from fhInitPlot")
      
      .fhList[[fhCurrent()]]
    })

    fhCurrent <- reactive({
      message(prefix, "fhCurrent, starting at ", .fhI)
      prefix <<- paste(prefix, " ", sep = "")
      if(input$nxt > nxtVal){
        message(prefix, "moving forwards to ", .fhI + 1)
        nxtVal <<- input$nxt
        if(.fhI < length(.fhList))
          .fhI <<- .fhI + 1

        message(prefix, "updating linval forward")
        if(.fhList[[.fhI]]@linearity == "fixed"){
          message(prefix, "turning button ON")
          linVal <- "ON"
        } else {
          message(prefix, "turning button OFF")
          linVal <- "OFF"
        }
        
        updateRadioButtons(session, "linearity",
                           selected = linVal)

      }

      if(input$prev > prevVal){
        message(prefix, "moving backwards to ", .fhI - 1)
        prevVal <<- input$prev      
        if(.fhI > 1)
          .fhI <<- .fhI - 1

        message(prefix, "updating linval backwards")
        if(.fhList[[.fhI]]@linearity == "fixed"){
          message(prefix, "turning button ON")
          linVal <- "ON"
        } else {
          message(prefix, "turning button OFF")
          linVal <- "OFF"
        }
        updateRadioButtons(session, "linearity",
                           selected = linVal)
      }

      if("SC" %in% .fhList[[.fhI]]@opts){
        message(prefix, "Switching to SC")
        debVal <- "SC"
      } else {
        message(prefix, "Switching to MC")
        debVal <- "MC"
      }
      
      updateRadioButtons(session, "debris",
                         selected = debVal)
      
      prefix <<- substring(prefix, 2)
      message(prefix, "returning from fhCurrent")
      .fhI
    })
    
    observe({
      if(input$exit > 0){
        stopApp()
      }
      
    })

    output$init <- renderPlot({
      message(prefix, "renderPlot")
      prefix <<- paste(prefix, " ", sep = "")
      plot(fhInitPlot(), init = TRUE, nls = TRUE, comps = TRUE)
      prefix <<- substring(prefix, 2)
      message(prefix, "returning from renderPlot")
    })

    output$flowNumber <- renderText({
      message(prefix, "flowNumber")
      paste("File", tags$b(fhCurrent()), "of", length(.fhList))
    })

  }
  runApp(shinyApp(ui = ui, server = server))
  return(.fhList)
}

selectPeaks <- function(fh, peakA, peakB){
  pA <- fh@histData[round(peakA, 0), c("xx", "intensity")]
  if(is.numeric(peakB))                 
    pB <- fh@histData[round(peakB, 0), c("xx", "intensity")]
  
  fh <- resetFlowHist(fh)

  if(is.numeric(peakB))
    fh@peaks <- as.matrix(rbind(pA, pB))
  else
    fh@peaks <- as.matrix(rbind(pA))
  
  colnames(fh@peaks) <- c("mean", "height")

  fh <- addComponents(fh)
  fh <- makeModel(fh)
  fh <- getInit(fh)

  return(fh)
}
