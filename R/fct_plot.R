#' Get an extracted ion chromatogram (XIC)
#' 
#' Helper function to extract an extracted ion chromatogram using mzR.
#'
#' @param spectra list of mzR spectra objects
#' @param header header extracted from mzML using mzR
#' @param mz numeric containing the mz value
#' @param ppm numeric the ppm tolerance window to extract the mz value
#' @param rt numeric containing the retention time in seconds
#' @param rt_tol numeric containing the retention time range to extract the XIC

#'
#' @return matrix of smoothed retention times (min) and intensities
#' @export
#'
#'
getXIC <- function(spectra,header,mz,ppm,rt,rt_tol, smooth = TRUE,sp = NULL){
  
  ## Filter scans by RT store for extraction of mz (likely more efficient)
  
  rtRange <- getRTRangeByRT(rt*60,rt_tol*60)   #+/- 1.5 min, mzR reads in seconds
  
  w <- which(header$retentionTime > rtRange[1] & header$retentionTime < rtRange[2])
  
  subsetOfSpectra <- spectra[w]
  
  rtVals <- header$retentionTime[w]
  ## Now extract the range of needed mz and intensities
  
  mzRange  <- getMassTolRange(mz,ppm)  # I like 6ppm in hi res orbitrap data that I'm working with
  

  
  i<-1## index for retention time
  out <- matrix(nrow = length(w), ncol = 2)
  out[,1] = rtVals/60
  out[,2] = 0
  
  for(spectrum in subsetOfSpectra){  
    
    w <- which(spectrum[,1] > mzRange[1] & spectrum[,1] < mzRange[2])
    
    if(length(w) > 0){
      
      spectrum <- spectrum[w,]
      
      if(length(w) > 1){
        
        # Handles multiple mz within tol, in thiscase extract the closest by mass.   
        # Hopefully intensity matches expected modality
        w <- which(abs(spectrum[,1]-mz)==min(abs(spectrum[,1]-mz)))
        
        spectrum <- spectrum[which(abs(spectrum[,1]-mz)==min(abs(spectrum[,1]-mz))),]  
        out[i,2] <- spectrum[2]
        
        
      }
      
    
    }
    
    i<-i+1 
  }
  
  mat <- out
  
  if(smooth){
    
    smooth_data <- smooth.spline(out[,1],out[,2],spar = sp)
    mat <- as.matrix(cbind(smooth_data$x, smooth_data$y))
    
    mat[which(mat[,2]<0),2] <- 0
    
    #plot(mat,type="l")
  }
  #dev(new)
  #plot(mat,type="l")
  
  mat
  
}





#' Plot an extracted ion chromatogram
#' 
#' Extracts a smoothed chormatogram for a mz value and retention time.
#'
#' @param mz numeric containing the mz value
#' @param ppm numeric the ppm tolerance window to extract the mz value
#' @param rt numeric containing the retention time in seconds
#' @param rt_tol numeric containing the retention time range to extract the XIC
#' @param ms1_mzml_file character vector listing the ms1 mzML file
#' @param ms2_mzml_file character vector listing the ms1 mzML file
#' @param plot_dir character vector listing the directory to plot in
#'
#' @export
#'

plotXIC <- function(ms1_mz,ms2_mz,ppm,rt,rt_tol,ms1_mzml_file,ms2_mzml_file,smooth = TRUE, sp = NULL, plot_dir=NULL){
  require(mzR)
  ms1_raw_data <- mzR::openMSfile(ms1_mzml_file)
  ms2_raw_data <- mzR::openMSfile(ms2_mzml_file)
  
  ms1_spectra <- mzR::spectra(ms1_raw_data)
  ms2_spectra <- mzR::spectra(ms2_raw_data)
  
  ms1_header <- mzR::header(ms1_raw_data)
  ms2_header <- mzR::header(ms2_raw_data) 
  
  ms1_XIC <- getXIC(ms1_spectra,ms1_header,mz = ms1_mz,ppm=ppm,rt=rt,rt_tol=rt_tol,smooth = smooth, sp = sp)
  ms2_XIC <- getXIC(ms2_spectra,ms2_header,mz = ms2_mz,ppm=ppm,rt=rt,rt_tol=rt_tol, smooth = smooth, sp = sp)  
  
  max_intensity <- max(
    max_ms1_intensity <- ms1_XIC[which( ms1_XIC[,2] == max(ms1_XIC[,2])),2],
    max_ms2_intensity <- ms2_XIC[which( ms2_XIC[,2] == max(ms2_XIC[,2])),2] )
  

  fileHandle <- paste(round(ms1_mz,4),"@",round(rt,2),".png",sep="")
  par(mar = c(4,4,2,2))
  png(width = 300,height = 300,filename = fileHandle)
  
  dev.neew()
  plot(ms1_XIC, type="l", col= rgb(0.1,0.1,.9,0.8),lwd=2, bty="n",xlab = "", ylab = "", ylim = c(0,max_intensity*1.3),main = sub(".png","",fileHandle))
  
  lines(ms2_XIC,col=rgb(0.9,0.1,.1,0.8), lwd= 2)
  
  dev.off()
  
  
  fileHandle <- paste(round(ms1_mz,4),"@",round(rt,2),".pdf",sep="")
  
  pdf(width = 2.5,height = 2.5,file = fileHandle, pointsize = 7)
  
  plot(ms1_XIC, type="l", col= rgb(0.1,0.1,.9,0.8),lwd=1, bty="n",xlab = "", ylab = "", ylim = c(0,max_intensity*1.3),main = sub(".pdf","",fileHandle))
  
  lines(ms2_XIC,col=rgb(0.9,0.1,.1,0.8), lwd= 1)
  
  dev.off()
  dev.off()
  
  
  
  
  
}



#' Get the list of sim windows
#' 
#' Retrieves the scan windows for plotting functions.
#'
#' @return a character vector of scan windows
#' @export 
#'
getSIMWindows = function(){
  
  list(simStart = unique(scandef[,4]),
       simEnd = unique(scandef[,5]),
       
       srch = c(paste("WSIM",unique(scandef[,4]),unique(scandef[,5]), sep="_"),
                paste("NL",unique(scandef[,4]),unique(scandef[,5]), sep="_"))
  )
  
}

#' getRTRangebyRT
#'
#' @param RT 
#' @param tol 
#'
#' @return
#' @export
#'
#' 
getRTRangeByRT <-function(RT,tol){    
  c(RT-tol, RT+tol)
}



