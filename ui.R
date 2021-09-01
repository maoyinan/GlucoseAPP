library(shiny)
library(rsconnect)

source('functionsShiny.R')
load('variables')

ui <- fluidPage(

  titlePanel(title=div(img(src="bonbon.png", height = 70),
                       "CGM data visualization tool")),

  sidebarLayout(
    sidebarPanel(
      # textOutput("temp"),

      radioButtons(inputId= 'select_blockLab',
                   label = 'Select scheme to divide CGM into intervals:',
                   choices = c('3 months'='3m',
                               '2 weeks'='2wk',
                               '24 hours'='24h')),
      selectInput(inputId = "select_id",
                  label = "Select subject ID to visualise:",
                  choices = unique(datID$ID)),

      selectInput(inputId = "select_block",
                   label = 'Select interval ID:',
                  choices = unique(datID$block))
    ),

    # beginning of right side


    mainPanel(

      tabsetPanel(
        tabPanel("Notes",
                 includeHTML("notes.html")),

        tabPanel("Demographics",
                 fluidRow(
                   plotOutput("demsumm", height = 250*nrows.plot(length(varsDem)))
                 )),
        tabPanel("Measurements",
                 fluidRow(
                   plotOutput("meassumm", height = 200*nrows.plot(length(varsMeas)))
                 )),
        tabPanel("CGM",
                 fluidRow(
                   h4('CGM over time with highlighted peaks'),
                   plotOutput("timeseries", height = 300),
                   h4('Magnified to selected range'),
                   sliderInput(inputId= "select_range",
                               label = "Range of days:",
                               min = 0, max = 100, value = c(0, 100), step=1, ticks=F),
                   plotOutput("timeseries1", height = 300),
                   h4('Daily fluctuation impression'),
                   plotOutput('tsDay',height = 300),
                   h4('Wavelet power spectrum'),
                   imageOutput('wvpower', height=400, width=800)
                 )),
        tabPanel("Features",
                 fluidRow(
                   plotOutput("featsumm", height = 200*nrows.plot(length(varsFeat)))
                 )),
        # tabPanel("Data of CGM features",
        #          h4("Summary glucose data for all patients"),
        #          downloadButton('downloadData', 'Download Data', class="btn-xs btn-info"),
        #          div(dataTableOutput("tb_feat"), style = "font-size:80%")),
        tabPanel("Correlation",
                 fluidRow(
                   selectInput(inputId= "select_measure",
                               label = "Select body measurement to compare with CGM features",
                               choices = titlesMeas, selected= titlesMeas[1]),
                   plotOutput("measfeat",
                              width = 220*ncols.plot(length(varsFeat)),
                              height = 250*nrows.plot(length(varsFeat)))
                 ))
      )
    )
  )
)
