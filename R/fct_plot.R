#' Get the peak
#'
#' @param spectrum matrix of mz and intensities
#' @param mz target mz value
#' @param ppm tolerance window in ppm
#'
#' @return list of mz and intensities
#' @export
getPeak = function(spectrum, mz, ppm){
   
   mr = getMassTolRange(mz,ppm)
   
   w = which(spectrum[,1] > mr[1] & spectrum[,1] < mr[2])
   
   if(length(w)>0){ #TODO 1 or 0 was 1?
      
      dM = spectrum[w,1]-mz
      
      w = w[which.min(dM)]
      
      return(list(x = spectrum[w,1], y= spectrum[w,2]))
      
   }else{
      
      return(list(x = NA, y= NA))
   }
   
}


#' Plots the spectrum of precursor and aglycones
#'
#' @param spectrum matrix of mz and intensities
#' @param mz target mz value
#' @param ppm tolerance window in ppm
#' @param main character vector for the plot title
#' @param col integer color value
#' @export
#'
plotPeak = function(spectrum, mz, ppm, main = NULL,col = 1){
   
   pk0 = getPeak(spectrum,mz,ppm)
   
   pk1 = getPeak(spectrum, mz + 1.003355,ppm)
   
   pk2 = getPeak(spectrum, mz + 2*1.003355,ppm)
   
   plot(NA, xlim = c(mz-0.5, mz+2.5), ylim = c(0,pk0$y*1.1), type = 'n', bty = 'n', main = main,
        xlab = 'm/z', ylab = 'intensity', cex.axis = 0.8,las= 2)
   
   abline(h=0)
   
   lines(spectrum, type = 'h',col = 'grey50')
   
   lines(pk0$x, pk0$y, type = 'h', col=col, lwd = 2)
   
   if(!is.na(pk1$y)){
      
      if(pk1$y < pk0$y*0.5){
         
         lines(pk1$x, pk1$y, type = 'h', col=col, lwd = 2)
         
         if(!is.na(pk2$y)){
            
            if(pk2$y < pk1$y*0.5){
               
               lines(pk2$x, pk2$y, type = 'h', col=col, lwd = 2)
            }
         }
      }
   }
}

#' Plots the spectrum of precursor and aglycones using plotly
#'
#' @param spectrum matrix of mz and intensities
#' @param mz target mz value
#' @param ppm tolerance window in ppm
#' @param main character vector for the plot title
#' @param col integer color value
#' @param titl character the title of the plot
#' @export
#'
plotPeak.plotly = function(spectrum, mz, ppm, main = NULL,col = 1, title = NULL){
   
   
   pk0 = getPeak(spectrum,mz,ppm)
   pk0 = data.frame(x = pk0$x, y = pk0$y)
   
   pk1 = getPeak(spectrum, mz + 1.003355,ppm)
   pk1 = data.frame(x = pk1$x, y = pk1$y)
   
   pk2 = getPeak(spectrum, mz + 2*1.003355,ppm)
   pk2 = data.frame(x = pk1$x, y = pk1$y)
   
   
   spectrum.df = data.frame(x=spectrum[,1], y = spectrum[,2])
   ##########################
   fig = plotly::plotly_empty(spectrum.df, 
                               x = ~x, 
                               y = ~y, 
                               type = 'scatter', 
                               mode = 'none')#, 
   #name = 'precursor')#,
   #line = list( shape = 'spline', smoothing = bw,color = 'rgb(22, 96, 167)', width = 2)) 
   fig <-   plotly::add_segments(
      #Line Vertical
      fig,
      x = ~x,
      y = ~rep(0,times = length(x)),
      xend = ~x,
      yend = ~y,
      data = spectrum.df,
      inherit=T,
      line = list( shape = 'spline', smoothing = bw,color = 'rgb(10, 10, 10)', width = 1, opacity = 0.95)
   )
   
   
   fig <-   plotly::add_segments(
      #Line Vertical
      fig,
      x = ~x,
      y = ~rep(0,times = length(x)),
      xend = ~x,
      yend = ~y,
      data = pk0,
      inherit=T,
      line = list( shape = 'spline', smoothing = bw,color = 'rgb(256, 10, 10)',  opacity = 0.75, width = 2)
   )

   
 
   
   if(!is.na(pk1$y)){
      
      if(pk1$y < pk0$y*0.5){
         
         fig <-   plotly::add_segments(
            #Line Vertical
            fig,
            x = ~x,
            y = ~rep(0,times = length(x)),
            xend = ~x,
            yend = ~y,
            data = pk1,
            inherit=T,
            line = list( shape = 'spline', smoothing = bw,color = 'rgb(256, 10, 10)',  opacity = 0.5, width = 3)
         )
         
         if(!is.na(pk2$y)){
            
            if(pk2$y < pk1$y*0.5){
               
               fig <-   plotly::add_segments(
                  #Line Vertical
                  fig,
                  x = ~x,
                  y = ~rep(0,times = length(x)),
                  xend = ~x,
                  yend = ~y,
                  data = pk2,
                  inherit=T,
                  line = list( shape = 'spline', smoothing = bw,color = 'rgb(256, 10, 10)',  opacity = 0.5, width = 3)
               )
            }
         }
      }
   }
   
   #title = title
   ##########################
   fig = fig %>% plotly::layout(
      title = list(text = title,font = list(family = 'Times New Roman', size = 12, color = "#7f7f7f")),
      xaxis = list(range = list(mz- 0.25, mz + 2.5),title = "m/z",font = list(family ='Times New Roman', size = 7, color = "#7f7f7f")),
      yaxis = list(range = list(0, pk0$y*1.2), title = "Intensity",font = list(family ='Times New Roman', size = 7, color = "#7f7f7f"),exponentformat = "e"),
      showlegend = FALSE, margin = m)
   
   
   return(fig)  
 
  
}



#' Find peaks from a vector of intensities
#'
#' @param x vector of intensities
#' @param m minimum distance between peaks
#'
#' @return vector of peak indices
#' @export
find_peaks <- function (x, m = 3){
   
   shape <- diff(sign(diff(x, na.pad = FALSE)))
   
   pks <- sapply(which(shape < 0), FUN = function(i){
      z <- i - m + 1
      z <- ifelse(z > 0, z, 1)
      w <- i + m + 1
      w <- ifelse(w < length(x), w, length(x))
      if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
   })
   
   pks <- unlist(pks)
   
   pks
}
