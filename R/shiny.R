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
  
  ui <- fluidPage(
    ## sidebarLayout(
      fluidRow(
        column(4,
               wellPanel(
               htmlOutput("flowNumber", align = "center"),
               fluidRow(
                 column(5, actionButton("prev", label = "Previous")),
                 column(5, actionButton("nxt", label = "Next"))))),
        ## tags$hr(),
        column(3,
               wellPanel(radioButtons(inputId = "peakPicker",
                                      label = "Move peak:", 
                                      choices = list("A" = "A", "B" = "B"),
                                      selected = "A", inline = TRUE))),
        ## tags$hr(),
        column(3, actionButton("exit", label = "Return to R"))),
    fluidRow(
      plotOutput("init", click = "pointPicker"))
  )

  server <- function(input, output){
    nxtVal <- 0
    prevVal <- 0
    updateVal <- 0
    
    fhInitPlot <- reactive({
      xPt <- nearPoints(.fhList[[fhCurrent()]]@histData,
                        input$pointPicker, "x", "intensity",
                        threshold = 25, maxpoints = 1)
      if(nrow(xPt) > 0){
        if(input$peakPicker == "A"){
          .fhList[[fhCurrent()]] <<-
            sliderPeaks(.fhList[[fhCurrent()]], xPt[1,1],
                        .fhList[[fhCurrent()]]@init$Mb) 
          .fhList[[fhCurrent()]] <<- fhAnalyze(.fhList[[fhCurrent()]])
        } else {
          .fhList[[fhCurrent()]] <<-
            sliderPeaks(.fhList[[fhCurrent()]],
                        .fhList[[fhCurrent()]]@init$Ma, xPt[1,1])
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

      .fhI
    })
    
    observe({
      if(input$exit > 0){
        stopApp()
      }
      
    })

    fhNLS <- eventReactive(input$update, {
      .fhList[[fhCurrent()]] <<- fhAnalyze(.fhList[[fhCurrent()]])
    })

    output$plot_clickedpoints <- renderTable({
      res <- nearPoints(.fhList[[fhCurrent()]]@histData,
                        input$pointPicker, "x", "intensity",
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
      plot(fhInitPlot(), init = TRUE, nls = TRUE, comps = TRUE)
    })

    output$flowNumber <- renderText({
      paste("File", tags$b(fhCurrent()), "of", length(.fhList))
    })

  }
  runApp(shinyApp(ui = ui, server = server))
  return(.fhList)
}

## shinyApp(ui = ui, server = server)     
   
sliderPeaks <- function(fh, peakA, peakB){
  pA <- fh@histData[round(peakA, 0), c("x", "intensity")]
  pB <- fh@histData[round(peakB, 0), c("x", "intensity")]

  fh <- resetFlowHist(fh)
  
  fh@peaks <- as.matrix(rbind(pA, pB))
  colnames(fh@peaks) <- c("mean", "height")

  fh <- addComponents(fh)
  fh <- makeModel(fh)
  fh <- getInit(fh)

  return(fh)
}
