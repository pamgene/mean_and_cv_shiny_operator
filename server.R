library(shiny)
library(shinyjs)
library(tercen)
library(dplyr)
library(tidyr)
library(ggplot2)

############################################
#### This part should not be modified
getCtx <- function(session) {
  # retreive url query parameters provided by tercen
  query <- parseQueryString(session$clientData$url_search)
  token <- query[["token"]]
  taskId <- query[["taskId"]]
  
  # create a Tercen context object using the token
  ctx <- tercenCtx(taskId = taskId, authToken = token)
  return(ctx)
}
####
############################################

shinyServer(function(input, output, session) {
  
  output$body <- renderUI({
    fluidPage(
      shinyjs::useShinyjs(),
      tags$script(HTML('setInterval(function(){ $("#hiddenButton").click(); }, 1000*30);')),
      tags$footer(shinyjs::hidden(actionButton(inputId = "hiddenButton", label = "hidden"))),
      titlePanel("CV plots"),
      sidebarLayout(
        sidebarPanel(
          wellPanel(
            checkboxInput("dofit", "Show error model fit", value = TRUE),
            sliderInput("phigh", label = h5("Quantile for high signal spots"),min = 0.8, max = 1, value = 0.98),
            sliderInput("plow" , label = h5("Quantile for low signal spots"), min = 0, max = 0.2, value = 0.05)
          ),
          wellPanel(
            selectInput("plottype", label = h4("Select plot type"),
                        choices = list("CV plot" = 1, "SNR plot" = 2, "SD plot" = 3), selected = 1),
            checkboxInput("collapse", "Collapse panes", value = FALSE),
            checkboxInput("logx" , "Logarithmic x-axis", value = FALSE),
            textInput("xmin",    label = "x-axis lower limit", value = "0"),
            textInput("xmax",    label=  "x-axis upper limit", value = "auto"),
            textInput('ymin',    label = "y-axis lower limit", value = "0"),
            textInput('ymax',    label = "y-axis upper limit", value = "auto")
          )
        ),
        mainPanel(plotOutput("cvplot"),
                  tableOutput("fitresult")
        )
      )
    )
  })
  
  getCVData <- reactive({
    computeCVData(getData(session))
  })
  
  doFit <- reactive({
    result <- getCVData()
    if(!input$collapse){
      result = result %>% group_by(pane) %>% do(cvmodel(., pLow = input$plow, pHigh = input$phigh))
    } else {
      result = cvmodel(result, pLow = input$plow, pHigh = input$phigh)
    }
    result
  })
  
  cfgPlot <-reactive({
    
    suppressWarnings({
      #if (input$dofit){
      aResult <- doFit()
      #}
      xLim = as.numeric(c(input$xmin, input$xmax))
      yLim = as.numeric(c(input$ymin, input$ymax))
      
      
      if(input$plottype == 1){
        aPlot = cvPlot(aResult, showFit = input$dofit,
                       xLog = input$logx,
                       xLim = xLim,
                       yLim = yLim,
                       collapse = input$collapse)
        
      } else if (input$plottype == 2){
        aPlot = snrPlot(aResult, showFit = input$dofit,
                        xLog = input$logx,
                        xLim = xLim,
                        yLim = yLim,
                        collapse = input$collapse)
      } else if (input$plottype == 3){
        aPlot = sdPlot(aResult, showFit = input$dofit,
                       xLog = input$logx,
                       xLim = xLim,
                       yLim = yLim,
                       collapse = input$collapse)
      }
      aPlot = aPlot + theme_bw()
    })
  })
  
  output$cvplot = renderPlot({
    suppressWarnings({
      aPlot = cfgPlot()
    })
    return(aPlot)
  })
  
  output$fitresult = renderTable({
    if(input$dofit){
      aResult = doFit()
      if(!input$collapse){
        aTable = aResult %>% group_by(pane) %>%
          summarise(.,std.low = sqrt(identity1(ssq0)) , CV.high = sqrt(identity1(ssq1)))
      } else {
        aTable  = data.frame(pane = "All panes combined", std.low = aResult$ssq0[1], CV.high = sqrt(aResult$ssq1[1]))
      }
    } else {
      aTable = data.frame()
    }
    if (dim(aTable)[2] == 3)
    {
      aTable = data.frame(aTable, snr.high = -10*log10(aTable$CV.high))
      colnames(aTable) = c("Pane", "Std low signals", "CV high signals", "SNR high signals (dB)")
      return(aTable)
    }
  })
  
})

getData <- function(session) {
  ctx <- getCtx(session)
  ctx %>% 
    select(.ci, .ri, .y) %>%
    mutate(.color = ctx$select(ctx$colors[[1]]) %>% pull()) 
}

#' Computes stats per cell.
#'@return a data frame with various per cell stats as columns.
computeCVData <- function(df, value='.y', color='.color', groupBy=c('.ri','.ci')) {
  
  data <- data.frame(value=df[[value]], Color=df[[color]], df[groupBy])
  suppressWarnings({
    result = data %>% group_by_(.dots=groupBy) %>%
      summarise(m = mean(value),
                stdev = sd(value),
                cv = stdev/m,
                snr_db = 10 * log10(abs(m)/stdev),
                nreps = as.double(length(value)),
                lvar = var(log(value)),
                pane = identity1(Color))
  })
  result
}

pooledVarEst<-function(s2, n){
  est = sum(s2 * (n-1), na.rm = TRUE ) / (sum(n) - length(n))
  return(est)
}

identity1 = function(.).[1]

cvmodel <- function(aResult, pLow = 0.05, pHigh = 0.95, maxIter = 25) {
  
  bLow 	= aResult$m <= quantile(aResult$m, pLow);
  bHigh 	= aResult$m >= quantile(aResult$m, pHigh);
  var = aResult$stdev^2
  
  bModel = bLow|bHigh
  
  ssq0 = NaN
  ssq1 = NaN
  pres = NaN
  iter = 0
  tryCatch({
    while(any(bModel)){
      ssq0 = pooledVarEst( var[bLow], aResult$nreps[bLow]);
      lssq1 = median(aResult$lvar[bHigh], na.rm = TRUE)
      ssq1 = exp(lssq1) - 1
      #calculate presence values for all spots
      m0 = aResult$m;
      m0[m0<0] = 0
      pres = ( sqrt(ssq1) * m0) /(sqrt(ssq1)*m0 + sqrt(ssq0))
      
      bLow = pres < pLow
      bHigh = pres > pHigh
      
      if (all(!bLow) | all(!bHigh)){
        ssq0 = NaN
        ssq1 = NaN
        break
      }
      if ( all( (bLow|bHigh) == bModel) ){
        break
      } else {
        bModel = bLow|bHigh
        iter = iter + 1
      }
      if (iter > maxIter){
        break
      }
    }
  }
  , error = function(e)e
  )
  if (!is.nan(ssq0) & !is.nan(ssq1)){
    cvFit = sqrt( ( (ssq1)*(aResult$m^2) + (ssq0))/(aResult$m^2) )
    snrFit = -10 *log10(cvFit)
    sdFit = sqrt(ssq1 * aResult$m^2 + ssq0)
  } else {
    cvFit = NaN
    snrFit = NaN
    sdFit = NaN
  }
  
  data.frame(aResult, 
             ssq0   = ssq0, 
             ssq1   = ssq1, 
             cvFit  = cvFit, 
             snrFit = snrFit, 
             sdFit  = sdFit, 
             bHigh  = bHigh,  
             presence = pres)
}

cvPlot <- function(aFrame, xLim = c(NA,NA),xLog = FALSE, yLim = c(0, NA), showFit = TRUE, collapse = FALSE){
  suppressWarnings({
    
    if (is.na(yLim[2])) yLim[2] = 0.5
    
    if (showFit){
      aPlot = ggplot(aFrame, aes(x = m, y = cv, colour = pane, shape = bHigh) ) + geom_point()  + geom_line(aes(x = m, y = cvFit) , colour = "red")
    } else {
      aPlot = ggplot(aFrame, aes(x = m, y = cv)) + geom_point(colour ="blue")
    }
    if (!collapse){
      aPlot = aPlot + facet_wrap(~pane)
    }
    aPlot = aPlot + ylim(yLim)
    aPlot = aPlot + xlim(xLim)
    aPlot = aPlot + scale_colour_brewer(type = "qual")
    if(xLog){
      if (!any(is.na(xLim)) & all(xLim > 0)){
        aPlot = aPlot + scale_x_log10(limits = xLim)
      } else {
        aPlot = aPlot + scale_x_log10()
      }
    }
  })
  return(aPlot)
}

snrPlot = function(aFrame, xLim = c(NA,NA),xLog = FALSE, yLim = c(0, NA), showFit = TRUE, collapse = FALSE){
  suppressWarnings({
    
    if (is.na(yLim[2])) yLim[2] = 20
    
    if (showFit){
      aPlot = ggplot(aFrame, aes(x = m, y = snr_db, colour = pane, shape = bHigh) ) + geom_point()  + geom_line(aes(x = m, y = snrFit) , colour = "red")
    } else {
      aPlot = ggplot(aFrame, aes(x = m, y = snr_db ,colour = pane)) + geom_point()
    }
    if (!collapse){
      aPlot = aPlot + facet_wrap(~pane)
    }
    aPlot = aPlot + ylim(yLim)
    aPlot = aPlot + xlim(xLim)
    aPlot = aPlot + scale_colour_brewer(type = "qual") + ylab("Signal to Noise ratio (dB)")
    if(xLog){
      if (!any(is.na(xLim)) & all(xLim > 0)){
        aPlot = aPlot + scale_x_log10(limits = xLim)
      } else {
        aPlot = aPlot + scale_x_log10()
      }
    }
  })
  return(aPlot)
}

sdPlot = function(aFrame, xLim = c(NA,NA),xLog = FALSE, yLim = c(0, NA), showFit = TRUE, collapse = FALSE){
  suppressWarnings({
    
    if (is.na(yLim[2])) yLim[2] = 0.25* max(aFrame$m)
    
    if (showFit){
      aPlot = ggplot(aFrame, aes(x = m, y = stdev, colour = pane, shape = bHigh) ) + geom_point()  + geom_line(aes(x = m, y = sdFit) , colour = "red")
    } else {
      aPlot = ggplot(aFrame, aes(x = m, y = stdev ,colour = pane)) + geom_point()
    }
    if (!collapse){
      aPlot = aPlot + facet_wrap(~pane)
    }
    aPlot = aPlot + ylim(yLim)
    aPlot = aPlot + xlim(xLim)
    aPlot = aPlot + scale_colour_brewer(type = "qual") + ylab("Standard Deviation")
    if(xLog){
      if (!any(is.na(xLim)) & all(xLim > 0)){
        aPlot = aPlot + scale_x_log10(limits = xLim)
      } else {
        aPlot = aPlot + scale_x_log10()
      }
    }
  })
  return(aPlot)
}
