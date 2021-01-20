#' plot peaks that are DNA adducts.
#'
#' @param pk_group numeric index
#' @param mz mz null
#' @param bw bandwidth
#'
#' @export
#'
plot_peak_group <- function(pk_group = NULL,mz = NULL, bw = 0.35){
   
   library('plotly')
   con <- RSQLite::dbConnect(RSQLite::SQLite(),db_name)
   e = e2 = ht = NULL
   if(!is.null(pk_group)){
      e <- 	RSQLite::dbGetQuery(con, paste("SELECT * FROM assigned_peak_groups WHERE pk_group = ",pk_group))
      e2 <- 	RSQLite::dbGetQuery(con, paste("SELECT * FROM peak_group_data WHERE pk_group = ",pk_group))
      ht = RSQLite::dbGetQuery(con, paste("SELECT * FROM hit_table WHERE pk_group = ",pk_group))
      #title = RSQLite::dbGetQuery(con, paste("SELECT * FROM peak_group_data WHERE pk_group = ",pk_group))
      
      e = e[order(e$scan),]
   }
   
   
   
   
   RSQLite::dbDisconnect(con)
   
   mat = cbind(e$scan,e$rt_min,e$p_int, e$nl_int)
   mat = mat[order(mat[,1]),]
   
   mat = as.matrix(mat)
   
   rtdelt = mean(diff(mat[,2]))
   scandelt = min(diff(e$scan))
   
   
   mat = rbind(cbind(e$scan[1]-scandelt,e$rt_min[1]-rtdelt,0,0),
               mat,
               cbind(e$scan[nrow(e)]+scandelt,e$rt_min[nrow(e)]+rtdelt,0,0))
   
   mat = rbind(cbind(mat[1,1]-scandelt,mat[1,2]-rtdelt,0,0),
               mat,
               cbind(mat[nrow(mat),1]+scandelt,mat[nrow(mat),2]+rtdelt,0,0))
   mat = rbind(cbind(mat[1,1]-scandelt,mat[1,2]-rtdelt,0,0),
               mat,
               cbind(mat[nrow(mat),1]+scandelt,mat[nrow(mat),2]+rtdelt,0,0))
   
   
   smooth_data.p <- smooth.spline(mat[,2:3],spar = bw)
   smooth_data.b <- smooth.spline(mat[,c(2,4)],spar = bw)
   
   smat.p = cbind(smooth_data.p$x, smooth_data.p$y)
   smat.b = cbind(smooth_data.b$x, smooth_data.b$y)
   
   smat.p[which(smat.p[,2] < 0),2] = 0
   smat.b[which(smat.b[,2] < 0),2] = 0
   
   
   r = round(cor(smat.p[,2], smat.b[,2]),3)
   
   dat = data.frame('x' = round(mat[,2],2),
                    'scan' = mat[,1],
                    'yp' = round(smat.p[,2],0),
                    'yb' = round(smat.b[,2],0))
   

   
   
   if(!is.na(ht$ID[1])){
      if(length(ht$ID) > 1){
         ID = paste(ht$ID, collapse = ',')
      }else{
         ID = ht$ID
      }
      tt = paste(ID,'\n',e2$mz_ms1,'@',e2$rt_ms1)
   }else{
      tt = paste(e2$mz_ms1,'@',e2$rt_ms1)
   }
   
   m <- list(
      l = 100,
      r = 100,
      b = 50,
      t = 80,
      pad =0
   )
   
   mzr = wSIMCity::getMassTolRange(e2$mz_ms1,15)
   #mzr = wSIMCity::getMassTolRange(e2$mz_ms1+1.003355,15)
#e$rt_min[ which(e$mz > mzr[1] & e$mz < mzr[2])]# & e$rt_min > (e2$rt_ms1 - 0.5) & e$rt_min < (e2$rt_ms1+0.5))]
   #M0 peaks
   fig = plotly::plot_ly(dat, 
                         x = ~x, 
                         y = ~yp, 
                         type = 'scatter', 
                         mode = 'lines', 
                         name = 'precursor',
                         line = list( shape = 'spline', smoothing = bw,color = 'rgb(22, 96, 167)', width = 2)) 
   fig = plotly::add_trace(fig, y = ~yb, type = 'scatter',
                           name = 'ms2',
                           mode = 'lines',
                           line = list( shape = 'spline', smoothing = bw,color = 'rgb(205, 12, 24)', width = 2))
   

   
   
   
   fig = fig %>% plotly::layout(
      title = list(text = tt,font = list(family = 'Times New Roman', size = 12, color = "#7f7f7f")),
      xaxis = list(title = "rt (min)",font = list(family ='Times New Roman', size = 7, color = "#7f7f7f")),
      yaxis = list(title = "Intensity",font = list(family ='Times New Roman', size = 7, color = "#7f7f7f"),exponentformat = "e"),
      showlegend = FALSE, margin = m)
   #fig
   
   
   
   #############################################
   ##### Fig 2 raw sticks - peaks
   minx = min(dat$x)
   maxx = max(dat$x)
   
   fig2 = plotly::plotly_empty(dat, 
                         x = ~x, 
                         y = ~yp, 
                         #xaxis = list(range = c(minx,maxx)),
                         type = 'scatter', 
                         mode = 'none')#, 
                         #name = 'precursor')#,
                         #line = list( shape = 'spline', smoothing = bw,color = 'rgb(22, 96, 167)', width = 2)) 
   #fig2
   
   fig2 <-   plotly::add_segments(
       #Line Vertical
       fig2,
       x = ~x,
       y = ~rep(0,times = length(x)),
       
       xend = ~x,
       yend = ~yp,
       
       data = dat,
       
       line = list(color = 'rgb(22, 96, 167)'),

       inherit=T
      )
   
   fig2 <-   plotly::add_segments(
      #Line Vertical
      fig2,
      x = ~x,
      y = ~rep(0,times = length(x)),
      xend = ~x,
      yend = ~-1*yb,
      data = dat,
      line = list(color = 'rgb(205, 12, 24)'),
      inherit=T
   )
   
   #fig2
   
   fig2 = fig2 %>% plotly::layout(
      
      title = list(text = 'MS scans detected',
                   font = list(family = 'Times New Roman', size = 18, color = "#7f7f7f")),
     
      xaxis = list(title = "rt (min)",
                   font = list(family ='Times New Roman', size = 7, color = "#7f7f7f"), 
                   showgrid = T,
                   zeroline = F,
                   showline = F
                   ),
      
      yaxis = list(zeroline = T,showgrid = T,
                   title = "Intensity",
                   font = list(family ='Times New Roman', size = 7, color = "#7f7f7f"),
                   
                   exponentformat = "e"),
      
      showlegend = FALSE,
      
      margin = m)
   

   #############################################
   ### Plot the ms1 peak
   mz = e2$mz_ms1
   wmx = which.max((  lapply(spectra[e$scan], function(x) getPeak(x,mz,7)$y)))
   
   spectrum = spectra[[e$scan[wmx]]]
   
   fig3 = plotPeak.plotly(spectrum, mz,7,col='rgb(10, 10, 256)', title = 'MS1 isotopes')
   
   ax = list(title = 'm/z', showgrid= T)
   ay = list(zeroline = T, showgrid = T)
   
   fig3 <- fig3 %>% layout( xaxis = ax, yaxis = ay, margin = m)

   
   #############################################
   ### Plot the ms2 peak
   mz = e2$mz_ms2
   wmx = which.max((  lapply(spectra[(e$scan+1)], function(x) getPeak(x,mz,7)$y)))
   
   spectrum = spectra[[(e$scan+1)[wmx]]]
   
   fig4 = plotPeak.plotly(spectrum, mz,7,col='rgb(256, 10, 10)', title = 'MS2 isotopes')
   
   ax = list(title = 'm/z', showgrid= T)
   ay = list(zeroline = T, showgrid = T)
   
   fig4 <- fig4 %>% layout( xaxis = ax, yaxis = ay ,margin = m )
   
   
   ##############################################
   
   sfig = plotly::subplot(fig, fig2,fig3,fig4,nrows = 2, margin = 0.075)
   sfig = sfig %>% plotly::layout(title = tt )
   sfig
   
   wd = getwd()
   fileHandle <- paste(wd,'/',plot_dir,'/',sample_name,'_',pk_group,".html",sep="")
   htmlwidgets::saveWidget(sfig, fileHandle, selfcontained = T)
   
   
}


#' Create tables for hit lists
#'
#' @export
#'
createTables <- function(){
   
   library(DT)
   con <- RSQLite::dbConnect(RSQLite::SQLite(),db_name)
   
   ht = RSQLite::dbGetQuery(con, paste("SELECT 
   ID,
   pk_group,
   score_total as 'score',
   f1,
   f2,
   mz_ms1 as mz,
   mz_ms2 as mz2,
   ppm as 'delta_ppm',
   round(int_ms1,0) as int,
   round(int_ms2,0) as int,
   ratio,
   rho as 'peak correlation',
   fwhm1 as 'fwhm: ms1',
   fwhm2 as 'fwhm: ms2',
   
   rt_ms1 as rt,
   n_pk as '# scans'
                                       
   FROM hit_table WHERE n_pk > 4 AND rho > 0.5 and abs(ratio) < 1"))
   RSQLite::dbDisconnect(con)
   
   for(g in ht$pk_group){
      cat(paste(g,'\r'))
      plot_peak_group(g)
   }
   
   
   fileHandle <- paste('../',sample_name,'/plots/',sample_name,'_',ht$pk_group,".html",sep="")
   ht_href = paste('<a href="',
                   fileHandle,
                   '" target=\"popup\" onclick=\"window.open(\'',
                   fileHandle,
                   ',\'name\',\'width=200,height=200\')\">',                   
                   ht$pk_group,'</a>', sep = '')
   
   #fileHandle
   
   ht$pk_group = ht_href
   
   
   
   
   
   dt = DT::datatable(ht,escape = F,# extensions = 'Responsive',
                  
                  caption = paste('Table of hits', sample_name),
                  filter = 'top', 
                  options = list(
                     autoWidth = F,
                     pageLength = 20, 
                     scrollX = T, 
                     scrollY = T,
                     width =1200,
                     height = 800,
                     dom='t'
                     # fixedColumns = list(leftColumns = 3
                  ))#)
   
   wd = getwd()
   fileHandle <- paste(wd,'/',sample_dir,'/',sample_name,".html",sep="")
   htmlwidgets::saveWidget(dt, fileHandle, selfcontained = T)
   
   
}

getSmoothPeak = function(x,y, bw){
   if(length(x) > 3){
      ss = smooth.spline(x,y, spar = bw)
      ss = spline( ss$x, ss$y,n=length(ss$x)*50)
      ss$y[which(ss$y < 0)] = 0
      
      hm = max(ss$y) / 2
      wt = which(ss$y > hm)
      
      return(list(peak = ss, x1 = min(wt), x2 = max(wt), fwhm = ss$x[max(wt)] - ss$x[min(wt)]))
      
   }else{
      
      return(list(peak = NA, x1 = NA, x2 = NA, fwhm = NA))
   }
}

#plot_peak_group(52, bw = 0.35)

