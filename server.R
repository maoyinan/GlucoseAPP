library(shiny)
library(rsconnect)

source('functionsShiny.R')
load('variables')

server <- function(input, output, session) {

  #------- Initialize the Memory ----------
  selected_vals = reactiveValues(
    block_lab= '3m',
    id = 3, bl = '00000390000004000001')

  #------ Whenever any of the inputs are changed, it only modifies the memory----
  observe({
    req(input$select_blockLab,input$select_id,input$select_block)

    selected_vals$block_lab <- input$select_blockLab
    selected_vals$id <- input$select_id
    selected_vals$bl <- input$select_block
  })

  #------ Update all UI elements using the values stored in memory ------

  observe({
    updateSelectInput(session, "select_id",
                      choices = datID %>%
                        filter(label==selected_vals$block_lab) %>%
                        select(ID) %>%
                        unique %>% unname %>% unlist %>% as.character)
  })

  observe({
    dfID <- datID %>%
      filter(label==selected_vals$block_lab)
    allBlocks <- dfID %>% filter(ID==selected_vals$id) %>% select(block) %>%
      unname %>% unlist %>% as.character

    updateSelectInput(session, "select_block",
                      label= sprintf("Select interval out of %d:", length(allBlocks)),
                      choices= allBlocks)
  })

  observe({
    path <- switch(selected_vals$block_lab,
                   '2wk' = 'https://www.dropbox.com/s/myii27r52ukyiqm/processed2wk?dl=1',
                   '3m' = 'https://www.dropbox.com/s/rpuyntip2a4tzum/processed3m?dl=1',
                   '24h' = 'https://www.dropbox.com/s/ifkubqjgkqzpibb/processed24h?dl=1')
    load(url(path))
    id <- selected_vals$id
    bl <- selected_vals$bl
    datiT <-
      datB %>%
      filter(ID==id& block==bl)
    idxID <- which(datMeas$ID==id&datMeas$block==bl)
    if(length(idxID))
      tb_peaks <- lsPeaks[[idxID]]

    # input update
    idDays <- datiT$day
    updateSliderInput(session, "select_range", min=min(idDays), max=max(idDays)+1, value= c(0,1))

    # define outputs
    output$demsumm <- renderPlot({
      plot.summary(datDem, demOut,vars=varsDem, titles=titlesDem, id=id)
    })
    output$meassumm <- renderPlot({
      plot.summary(datMeas, measOut,vars=varsMeas, titles=titlesMeas, id=id, bl=bl)
    })
    output$timeseries <- renderPlot({
      plot.peaks(datiT,tb_peaks)
    })
    output$timeseries1 <- renderPlot({
      plot.peaks(datiT,tb_peaks,input$select_range)
    })
    output$tsDay <- renderPlot({
      plot.day(datiT,d=10)
    })
    output$wvpower <- renderImage({
      list(src = file.path('www',sprintf('power_ID%s.png',id)),
           contentType = 'image/png',
           width = 1600,
           height = 400,
           alt = "No image rendered")
    },deleteFile = F)

    # output$tb_feat  <-
    #   renderDataTable(datFeat)

    output$featsumm <- renderPlot({
      plot.summary.density(datFeat, vars=varsFeat, titles=titlesFeat, id=id, bl=bl)
    })
    output$measfeat <- renderPlot({
      plot.meas.feat(datMeasFeat, input$select_measure, id=id, bl=bl)
    })

  })
}
