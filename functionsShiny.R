library(gam)
library(grid)
library(ggpmisc)
library(tidyverse)
filter <- dplyr::filter
library(lubridate)
library(WaveletComp)

# vague matching tables by ID and between time window
group.find <- function(ID,vals, tb){
  x <- mapply(function(a,b) (vals>=a)&(ID==b),
              tb$date,tb$ID)
  if(sum(x)==0) return(which(ID==tb$ID & vals<=tb$date)[1])
  else return(last(which(x)))
}

mage.locate <- function (x, n = 1)
{# mean amplitude of glycemic excursion
  # iglu::mage # (Service, F. J. & Nelson, R. L. (1980))
  x %>%
    mutate(abs_diff_mean = abs(record - mean(record, na.rm = TRUE)),
           mage_loc= abs_diff_mean >
             (n * sd(record, na.rm = TRUE)))
}

# separate time series into smaller blocks
block.separate <- function(block_lab,block_len, datT){
  dat1 <- datT %>%
    select(block=!!sprintf('block%s',block_lab), everything())
  dat1 %>%
    group_by(block) %>%
    summarise(t1=max(time),t2=min(time),
              spread=as.numeric(difftime(t1,t2, units='days'))) %>%
    filter(spread>=block_len*.95) %>%
    select(block) %>%
    unlist -> blocks

  dat1 %>%
    group_by(block) %>%
    group_modify(~ {
      mage.locate(.x)
    }) %>%
    mutate(hr= as.numeric(difftime(time,min(time),units = 'hours')),
           day= floor(as.numeric(difftime(time,min(time),units = 'days')))) %>%
    ungroup %>%
    filter(block%in%blocks) %>%
    select(ID,block,time,hr,day,tod,record,mage_loc, abs_diff_mean)
}

# aggregate into peaks and valleysd ----------------------------------------------------------

peakdet <- function(v, delta, x = NULL)
{
  maxtab <- NULL
  mintab <- NULL

  if (is.null(x))
  {
    x <- seq_along(v)
  }

  if (length(v) != length(x))
  {
    stop("Input vectors v and x must have the same length")
  }

  if (!is.numeric(delta))
  {
    stop("Input argument delta must be numeric")
  }

  if (delta <= 0)
  {
    stop("Input argument delta must be positive")
  }

  mn <- Inf
  mx <- -Inf

  mnpos <- NA
  mxpos <- NA

  lookformax <- TRUE

  for(i in seq_along(v))
  {
    this <- v[i]

    if (this > mx)
    {
      mx <- this
      mxpos <- x[i]
    }

    if (this < mn)
    {
      mn <- this
      mnpos <- x[i]
    }

    if (lookformax)
    {
      if (this < mx - delta)
      {
        maxtab <- rbind(maxtab, data.frame(
          pos = mxpos,
          val = mx,
          index = i))

        mn <- this
        mnpos <- x[i]

        lookformax <- FALSE
      }
    }
    else
    {
      if (this > mn + delta)
      {
        mintab <- rbind(mintab, data.frame(
          pos = mnpos,
          val = mn,
          index = i))

        mx <- this
        mxpos <- x[i]

        lookformax <- TRUE
      }
    }
  }

  list(maxtab = maxtab, mintab = mintab)
}

summary_peak <- function(dati, delta){
  dat_peaks <- peakdet(dati$record,delta=delta, x=dati$hr)
  if(is.null(dat_peaks$maxtab)) tb_peaks <- NULL
  else if (nrow(dat_peaks$maxtab)<3) tb_peaks <- NULL
  else tb_peaks <- bind_cols(dat_peaks$maxtab,flag='peak') %>%
      bind_rows(bind_cols(dat_peaks$mintab,flag='valley')) %>%
      arrange(pos) %>%
      mutate(index0 = lag(index,1,default=1)) %>%
      rowwise() %>%
      mutate(gap= max(diff(dati$hr[index0:index]))) %>%
      ungroup() %>%
      mutate(
        pos_bef= lag(pos),
        pos_aft= lead(pos),
        val_bef= lag(val),
        val_aft= lead(val),
        diff_hour_bef= c(NA,diff(pos)),
        diff_hour_aft= c(diff(pos),NA),
        change_bef= c(NA,diff(val)),
        change_aft= c(diff(val),NA),
        rate_bef= change_bef/diff_hour_bef,
        rate_aft=change_aft/diff_hour_aft) %>%
      filter(gap<=1) %>%
      na.omit

  tb_peaks
}

aggregate_features <- function(dati,varsFeat,tb_peaks, tb_fit){
  feat <- rep(NA,length(varsFeat))
  ## build feature matrix
  feat[1] <- max(dati$record, na.rm=T)
  feat[2] <- min(dati$record, na.rm=T)
  feat[3] <- mean(dati$record, na.rm=T)
  feat[4] <- sd(dati$record, na.rm=T)
  feat[5] <- dati %>%
    filter(mage_loc) %>%
    summarise(mean(abs_diff_mean)) %>%
    unlist
  feat[6] <- mean(dati$record>13.9, na.rm = T)
  feat[7] <- mean(dati$record>10, na.rm = T)
  feat[8] <- mean(between(dati$record,3.9,10), na.rm = T)
  feat[9] <- mean(dati$record<3.9, na.rm = T)
  feat[10] <- mean(dati$record<3, na.rm = T)


  ## calculate peak related features
  if(!is.null(tb_peaks)) feat[11:14] <- tb_peaks %>%
    filter(flag=='peak') %>%
    select(change_bef:rate_aft) %>%
    summarise_all(mean, na.rm=T) %>%
    unlist

  ## time of day spline fit
  feat[15:(15+23)] <- tb_fit[[2]] %>%
    filter(x%in%0:23) %>%
    select(y) %>%
    unlist - mean(dati$record, na.rm=T)
  feat
}

# summarise patient level data in histogram (continuous) or bar plot (categorical) -------
# highlight bar where selected patient id belongs
dem.summary <- function(dat, vars){
  out <- list()

  for(i in seq_along(vars)){
    var <- vars[i]
    obs <- dat %>%
      select(.dots=var) %>%
      unlist()

    if(is.character(obs)) {
      x_labels <- names(table(obs))
      x_at <- seq_along(x_labels)
      x_rg <- c(.5, max(x_at)+.5)
      bar_x <- x_at
      bar_y <- table(obs)
    }
    else {
      temp <- hist(obs,plot=F)
      x_at = x_labels <- temp$breaks
      x_rg <- range(x_at)
      bar_x <- temp$mids
      bar_y <- temp$counts
    }
    y_rg <- c(0,max(bar_y))
    bar_wid <- diff(bar_x[1:2])
    out[[i]] <- list(
      x_at=x_at,
      x_labels=x_labels,
      x_rg=x_rg,
      bar_x=bar_x,
      y_rg=y_rg,
      bar_y=bar_y,
      bar_wid=bar_wid
    )
  }
  names(out) <- vars
  out
}
# dem_out <- dem.summary(dat=datDem,vars=colnames(datDem)[-1])
ncols.plot <- function(nVars) min(ceiling(sqrt(nVars)),5)
nrows.plot <- function(nVars) ceiling(nVars/ncols.plot(nVars))
plot.summary <- function(dat, dem_out, vars, titles, id=NULL, bl=NULL){
  # plot the histogram of selected variable.
  # highlight region of selected patient ID
  color_bar <- c(fill='lightskyblue3',col='white')
  color_hl <- c(fill= 'tomato',col='white')

  ncolPlot <- ncols.plot(length(vars))
  pushViewport(plotViewport(c(0,1,1,0)))
  pushViewport(viewport(
    layout=grid.layout(nrow=nrows.plot(length(vars)), ncol=ncolPlot)))
  for(k in seq_along(vars)){
    var <- vars[k]

    pushViewport(viewport(layout.pos.row=(k-1)%/%ncolPlot+1,layout.pos.col = ifelse(k%%ncolPlot==0, ncolPlot, k%%ncolPlot)))
    x_at=dem_out[[var]]$x_at
    x_labels=dem_out[[var]]$x_labels
    x_rg=dem_out[[var]]$x_rg
    bar_x=dem_out[[var]]$bar_x
    y_rg=dem_out[[var]]$y_rg
    bar_y=dem_out[[var]]$bar_y
    bar_wid=dem_out[[var]]$bar_wid

    # highlight bar
    if(!is.null(id)){
      id <- as.numeric(id)
      dat %>%
        filter(ID==id) -> dat1
      if(!is.null(bl)){
        dat1 %>%
          filter(block==bl) -> dat1
      }
      dat1 %>%
        select(.dots=var) %>%
        unlist() -> y0
      if(is.character(y0)){
        hl_at <- which(y0==x_labels)
      }
      else{
        hl_at <- which(y0<=x_at[-1])[1]
      }
    }


    # plot
    pushViewport(plotViewport(c(4,2,1,1),
                              xscale=x_rg, yscale=y_rg))
    grid_bar(bar_y, bar_x, bar_wid, color_bar)
    if(!is.null(id)){
      grid_bar(bar_y[hl_at],bar_x[hl_at], bar_wid, color_hl)
      # grid.text(label = sprintf('ID = %d', id),  x = bar_x[hl_at], y = bar_y[hl_at], default.units = 'native', just = 'top', gp=gpar(col='black'))
    }
    grid.lines(x=c(0,1),y=0,default.units = 'npc')
    grid.xaxis(at = x_at, label = x_labels, gp=gpar(cex=.8))
    grid.text(label = titles[k],  x = .5, y = unit(-2,'lines'))
    grid.yaxis(at = pretty(y_rg), label = pretty(y_rg), gp=gpar(cex=.8))

    popViewport()
    popViewport()

  }

  popViewport()
  popViewport()

}
# density plot only for numberic variable summary
plot.summary.density <- function(dat, vars, titles, id=NULL, bl=NULL){

  color_hl <- c(fill= 'tomato',col='white')

  ncolPlot <- ncols.plot(length(vars))
  pushViewport(plotViewport(c(0,1,1,0)))
  pushViewport(viewport(
    layout=grid.layout(nrow=nrows.plot(length(vars)), ncol=ncolPlot)))
  for(k in seq_along(vars)){
    var <- vars[k]
    xi <- dat %>% select(.dots=var) %>% unlist %>% na.omit
    dens <- density(xi)
    x <- dens$x
    y <- dens$y

    pushViewport(viewport(layout.pos.row=(k-1)%/%ncolPlot+1,layout.pos.col = ifelse(k%%ncolPlot==0, ncolPlot, k%%ncolPlot)))

    # highlight selected
    if(!is.null(id)){
      id <- as.numeric(id)

      dat %>%
        filter(ID==id) -> dat1
      if(!is.null(bl)){
        dat1 %>%
          filter(block==bl) -> dat1
      }
      dat1 %>%
        select(.dots=var) %>%
        unlist() -> x0
      y0 <- y[which(x0<=x)[1]]
    }


    # plot
    pushViewport(plotViewport(c(3,2,1,1),
                              xscale=range(x), yscale=range(y)*c(1,1.05)))
    grid.rect()
    grid.lines(x,y, default.units = 'native')
    if(!is.null(id)){
      grid.points(x0,y0,pch=21, gp=gpar(fill=color_hl[1], col=color_hl[2]))
      # grid.text(label = sprintf('ID = %d', id),  x =x0, y = y0, default.units = 'native', just = 'bottom', gp=gpar(col='black'))
    }
    grid.xaxis(at = axis.trunc(x), label = axis.trunc(x), gp=gpar(cex=.8))
    grid.text(label = titles[k],  x = .5, y = unit(-2,'lines'))
    grid.yaxis(at =axis.trunc(y), label = axis.trunc(y), gp=gpar(cex=.8))

    popViewport()
    popViewport()

  }

  popViewport()
  popViewport()

}

grid_bar <- function(bar_y, bar_x, bar_wid, color_bar){
  # plot bars with heights given
  nbar <- length(bar_x)
  grid.rect(x=bar_x, y= rep(0, nbar),just='bottom',
            width=rep(bar_wid,nbar), height=bar_y,
            default.units='native',
            gp=gpar(fill=color_bar['fill'], col=color_bar['col']))
}


# plot time series of CGM --------------------------------------------------------------------

day.fit <- function(dati,d){
  xi1 <-
    dati %>%
    select(tod) %>%
    unlist
  yi <-
    dati %>%
    select(record) %>%
    unlist
  d1 = data.frame(x=xi1,y=yi)
  d2 = data.frame(x=seq(0,24,1/240))
  fit = gam(as.formula(sprintf('y~s(x,df=%d)',d)),data=d1)
  pre = predict.Gam(fit,newdata = d2)
  d2$y <- pre
  list(d1,d2)
}
plot.day <- function(dati,d){
  temp <- day.fit(dati,d)
  d1 <- temp[[1]]
  d2 <- temp[[2]]

  # plot daily fluctuations with time of day variable tod
  pushViewport(plotViewport(c(4,4,1,1),
                            xscale=c(0,24),yscale=range(d1$y)))
  grid.rect()

  for(k in 1:23) grid.lines(unit(c(k,k),'native'),c(0,1),gp=gpar(col=gray(0.8)))
  # for(k in c(0,3,6,9,12,15,18,21,24))grid.lines(unit(c(k,k),'native'),c(0,1),gp=gpar(lwd=1))

  grid.points(d1$x,d1$y,gp=gpar(cex=0.33,col=rgb(0,0,0,0.2)),pch=16)
  grid.lines(d2$x,d2$y,default.units = 'native')
  grid.lines(c(0,1,1,0,0),c(0,0,1,1,0))
  grid.xaxis(at=c(0,3,6,9,12,15,18,21,24))
  grid.yaxis()
  grid.text('Time of day (hour)',y=unit(-2.5,'lines'))
  grid.text('Glucose (mmol/L)',x=unit(-3,'lines'),rot=90)
  popViewport()
}

# grid plot  --------------------------------------------------

clock.label <- function(t, tod0){
  ticks <- unique(floor(t))
  idx <- which((ticks+tod0)%%6==0)
  list(at=ticks[idx],
       label=((ticks+tod0)%%24)[idx])
}
plot.peaks <- function(dati,tb_peaks, rg=NULL) {
  if(!is.null(rg)){
    dati <- dati %>%
      filter(between(hr, 24*rg[1], 24*rg[2]))
    tb_peaks <- tb_peaks %>%
      filter(between(pos, 24*rg[1], 24*rg[2]))
  }
  xi <- dati$hr
  yi <- dati$record

  peak_colors <- c('tomato', 'lightskyblue3')

  pushViewport(plotViewport(c(5,4,3,1),xscale=range(xi),yscale=c(0,max(yi))))
  tbp <- tb_peaks %>% filter(flag=='peak')
  for(k in seq(nrow(tbp))){
    grid.rect(
      x=max(min(xi),tbp$pos_bef[k]),
      width = tbp$pos[k]-max(min(xi),tbp$pos_bef[k]),
      just = 'left',
      default.units = 'native',
      gp = gpar(col=NA,fill='tomato', alpha=.2))
  }
  grid.lines(xi,yi,default.units = 'native')
  grid.points(tbp %>% select(pos) %>% unlist,
              tbp %>% select(val) %>% unlist,
              gp=gpar(cex=1,fill=peak_colors[1]),pch=21) # peak
  temp <- tbp %>% select(pos_bef) %>% unlist
  idx <- temp>=min(xi)
  temp1 <- tbp %>% select(val_bef) %>% unlist
  grid.points(temp[idx],
              temp1[idx],
              gp=gpar(cex=1,fill=peak_colors[2]),pch=21) # valley

  grid.lines(c(0,1,1,0,0),c(0,0,1,1,0))
  grid.text('Time of day (hour)',y=unit(1,'npc')+unit(2.5,'lines'))
  tod0 <- floor(dati$tod[1])
  lab_tod <- clock.label(xi,tod0)
  grid.xaxis(at= lab_tod$at, label = lab_tod$label, main=F)

  xa <- pretty(xi/24)
  xa <- xa[xa<=max(xi/24)&xa>=min(xi/24)]
  grid.xaxis(at= xa*24, label=xa)
  grid.yaxis()
  grid.text('Time (day)',y=unit(-2.5,'lines'))
  grid.text('Glucose (mmol/L)',x=unit(-3,'lines'),rot=90)
  popViewport()
  pushViewport(plotViewport(c(5,4,3,1)))
  grid.legend(labels=c('Peak', 'Valley'),
              gp=gpar(fill= peak_colors),
              pch= 21, nrow=1, ncol=length(peak_colors),
              vp= viewport(y=unit(-4,'lines'), height=unit(2,'lines')))
  popViewport()
}


axis.trunc <- function(x){
  x_axis <- pretty(x)
  x_axis <- x_axis[between(x_axis,min(x,na.rm=T),max(x,na.rm=T))]
  x_axis
}
range.bigger <- function(x){
  a <- min(x,na.rm = T)
  b <- max(x,na.rm = T)
  dt <- (b-a)*.05
  c(a-dt,b+dt)
}
plot.meas.feat <- function(datMeasFeat, varMeasCho='BMI', id=NULL, bl=NULL){
  varsFeat <- datMeasFeat %>% select(Feature) %>% unique %>% na.omit %>% unlist

  color_hl <- c('tomato')

  ncolPlot <- ncols.plot(length(varsFeat))

  pushViewport(plotViewport(c(0,0,1,0)))
  pushViewport(viewport(
    layout=grid.layout(nrow=nrows.plot(length(varsFeat)), ncol=ncolPlot)))
  for(k in seq_along(varsFeat)){
    var <- varsFeat[k]
    dat <- datMeasFeat %>%
      filter(Measurement==varMeasCho&Feature==var)

    x <- dat$xvalue
    y <- dat$yvalue
    x_axis <- axis.trunc(x)
    y_axis <- axis.trunc(y)

    pushViewport(viewport(layout.pos.row=(k-1)%/%ncolPlot+1,layout.pos.col = ifelse(k%%ncolPlot==0, ncolPlot, k%%ncolPlot)))

    # highlight bar
    if(!is.null(id)&!is.null(bl)){
      dat %>%
        filter(ID==id&block==bl) %>%
        select(xvalue,yvalue) %>%
        unlist() -> xy0
    }


    pushViewport(plotViewport(c(4,2,1,1),
                              xscale=range.bigger(x), yscale=range.bigger(y)))
    grid.rect()
    grid.points(x=x,y=y)
    if(!is.null(id)){
      grid.points(x = xy0[1], y = xy0[2], default.units = 'native', pch=21, gp=gpar(fill=color_hl))
    }
    grid.xaxis(at = x_axis, label = x_axis, gp=gpar(cex=.8))
    grid.text(label = var,  x = .5, y = unit(-2.5,'lines'))
    grid.yaxis(at = y_axis, label = y_axis, gp=gpar(cex=.8))

    popViewport()
    popViewport()

  }

  popViewport()
  popViewport()

}
