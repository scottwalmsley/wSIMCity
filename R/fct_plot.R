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
    n = spline(smooth)
    
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

plotXIC <- function(ms1_mz,ms2_mz,ppm,rt,rt_tol,ms1_mzml_file,ms2_mzml_file,smooth = TRUE, sp = NULL, sampleDir=NULL, title = NULL){
  require(mzR)
  ms1_raw_data <- mzR::openMSfile(ms1_mzml_file)
  ms2_raw_data <- mzR::openMSfile(ms2_mzml_file)
  
  ms1_spectra <- mzR::spectra(ms1_raw_data)
  ms2_spectra <- mzR::spectra(ms2_raw_data)
  
  ms1_header <- mzR::header(ms1_raw_data)
  ms2_header <- mzR::header(ms2_raw_data) 
  
  ms1_XIC <- getXIC(ms1_spectra,ms1_header,mz = ms1_mz,ppm=ppm,rt=rt,rt_tol=rt_tol, smooth = smooth, sp = sp)
  ms2_XIC <- getXIC(ms2_spectra,ms2_header,mz = ms2_mz,ppm=ppm,rt=rt,rt_tol=rt_tol, smooth = smooth, sp = sp)  
  
  max_intensity <- max(
    max_ms1_intensity <- ms1_XIC[which( ms1_XIC[,2] == max(ms1_XIC[,2])),2],
    max_ms2_intensity <- ms2_XIC[which( ms2_XIC[,2] == max(ms2_XIC[,2])),2] )
  
  
  if(!is.null(title) & !is.null(sampleDir)){
    fileHandle <- paste(sampleDir,"/plots/",title,"_",round(ms1_mz,4),"@",round(rt,2),".png",sep="")
  }else{
    fileHandle <- paste(sampleDir,"/plots/",round(ms1_mz,4),"@",round(rt,2),".png",sep="")
  }
  
  par(mar = c(4,4,2,2))
  png(width = 300,height = 300,filename = fileHandle)
  
  
  plot(ms1_XIC, type="l", col= rgb(0.1,0.1,.9,0.8),lwd=2, bty="n",xlab = "", ylab = "", ylim = c(0,max_intensity*1.3),main = sub(".png","",fileHandle))
  
  lines(ms2_XIC,col=rgb(0.9,0.1,.1,0.8), lwd= 2)
  
  dev.off()
  
  
  if(!is.null(title) & !is.null(sampleDir)){
    fileHandle <- paste(sampleDir,"/plots/",title,"_",round(ms1_mz,4),"@",round(rt,2),".pdf",sep="")
  }else{
    fileHandle <- paste(sampleDir,"/plots/",round(ms1_mz,4),"@",round(rt,2),".pdf",sep="")
  }
  
  pdf(width = 2.5,height = 2.5,file = fileHandle, pointsize = 7)
  
  plot(ms1_XIC, type="l", col= rgb(0.1,0.1,.9,0.8),lwd=1, bty="n",xlab = "", ylab = "", ylim = c(0,max_intensity*1.3),main = sub(".pdf","",fileHandle))
  
  lines(ms2_XIC,col=rgb(0.9,0.1,.1,0.8), lwd= 1)
  
  dev.off()
  dev.off()
  
  
}



#' Plot XIC of known ions
#' 
#' Generates an EIC of a known compound.
#'
#' @param scandef_file 
#' @param sampleDir 
#' @param ppm 
#' @param rt_tol 
#' @export
#'
plotKnowns <- function(scandef_file, sampleDir,ppm = 6, rt_tol = 1.5){
  
  data("knowns.db")
  
  w <- which(knowns.db$Type =="M")
  
  ms1_mz <- knowns.db$MZ[w]
  ms2_mz <- knowns.db$MZ[-w]
  
  ms1_nm <- knowns.db$Adduct[w]  
  ms2_nm <- knowns.db$Adduct[-w]
  
  rt <- knowns.db$RT[w]
  
  scandef <- read.delim(scandef_file)
  windows <- getSIMWindows(scandef)
  
  
  lf <- list.files(path = sampleDir, pattern = "mzML", full.names = TRUE, recursive = TRUE)
  
  for(i in 1:length(ms1_mz)){
    
    w = which(windows$simEnd > ms1_mz[i] & windows$simStart < ms1_mz[i])
    
    if(length(w)>1){
      
      w <- w[1] 
      
    }
    
    if(length(w)>0){  
      
      plotXIC(ms1_mz[i],ms2_mz[i],ppm=ppm,rt=rt[i],rt_tol=rt_tol,ms1_mzml_file = lf[w+10],ms2_mzml_file = lf[w],smooth = FALSE, title = ms1_nm[i], sampleDir = sampleDir)
      
    }
    
    
  }
  
}


#' Get the list of sim windows
#' 
#' Retrieves the scan windows for plotting functions.
#' @param scandef
#' @return a character vector of scan windows
#' @export 
#'
getSIMWindows = function(scandef){
  
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



#' getPlots
#'
#' @param searchResultList 
#' @param scandef_file 
#' @param sampleDir 
#' @param min_score 
#' @param min_intensity 
#' @param ppm 
#' @param rt_tol 
#'
#' @export
#'
#' @examples
getPlots <- function(searchResultList,scandef_file,sampleDir,min_score = 0.95, min_intensity = 1000, ppm = 6,rt_tol = 1.5 , min_ratio = 10){
  
  require(mzR)
  
  score_feature <- unlist(lapply(searchResultList, function(x)x$best_candidate$score_feature ))
  mz <- unlist(lapply(searchResultList, function(x) x$search$mz))
  intensity_1 <- unlist(lapply(searchResultList, function(x) x$search$intensity))
  intensity_2 <- unlist(lapply(searchResultList, function(x) x$best_candidate$intensity))
  ratio <- unlist(lapply(searchResultList,function(x) x$best_candidate$ratio_intensity))
  
  w <-  which(score_feature > min_score & mz > 330 & intensity_1 > min_intensity & intensity_2 > min_intensity & ratio > (1/min_ratio) & ratio < min_ratio)
  
  sub_list <- searchResultList[w]
  
  ms1_mz <- unlist(lapply(sub_list, function(x) x$search$mz))
  
  rt <- unlist(lapply(sub_list, function(x) x$search$rt))
  
  ms2_mz <- unlist(lapply(sub_list, function(x) x$best_candidate$mz))
  
  scandef <- read.delim(scandef_file)
  
  windows <- getSIMWindows(scandef)
  
  
  lf <- list.files(path = sampleDir, pattern = "mzML", full.names = TRUE, recursive = TRUE)
  
  plot_dir <- paste(sampleDir,"/plots",sep="")
  
  i<-1
  
  for(mz in ms1_mz){
    
    w = which(windows$simEnd > mz & windows$simStart < mz)
    
    if(length(w)>1){
      
      w <- w[1] # get the window with the highest sim range
      
    }
    
    if(length(w)>0){  
      
      plotXIC(mz,ms2_mz[i],ppm = ppm, rt = rt[i], rt_tol = rt_tol, ms1_mzml_file = lf[w+10],ms2_mzml_file = lf[w],smooth = F, sampleDir = sampleDir)
      
    }
    
    i <- i+1
  }
  
}



#' Detect a peak in an XIC
#'
#' @param XIC matrix of the XIC
#' @param rt numeric rt of the peak in minutes
#'
#' @return matrix of the detected peak
#' @export
#'
detectPeak <- function(XIC,rt){
  
  n = spline(XIC[,1],XIC[,2])
  n$y[which(n$y<0)] = 0  
  
  #plot(n, type= "l", ylim = c(-30000,60000))
  
  der = diff(n$y)
  der = c(0,der)
  rder = round(der)
  der2 = diff(rder)
  der2 = c(0,der2)
  der3 = diff(der2)
  der3 = c(0,der3)
  
  
  mat  = matrix(nrow = length(n$x),ncol=3)
  mat[,1] = sign(rder)
  mat[,2] = sign(der2)
  mat[,3] = sign(der3)
  mat[which(mat==0)] = -1
  mat = cbind(mat,n$x,n$y)
  
  
  
  j=k=1
  last = -1
  peaks = list()
  newPk = NULL
  
  #which.min(abs(newPk[,1] - rt))
  
  
  for(i in 1:nrow(mat)){
    
    vec = mat[i,]
    
    if(last == -1 & sum(vec[1:2])==2){
      if(!is.null(newPk)){
        w.max = which.max(newPk[,2])
        peaks[[k]] = list("XIC" =n , "area" =sum(newPk[,2],na.rm = T), max.intensity = newPk[w.max,2], "peak" = newPk,"mode_x" = newPk[w.max,1] )
        newPk = NULL
        k=k+1
      }
      if(is.null(newPk)){
        newPk = vec[4:5]
      }
     
    }
    if(!is.null(newPk)){
     
      newPk = rbind(newPk,vec[4:5])
    }
    
    j=j+1
    
    last = vec[1]	
    
    
  }
  rt.pk = sapply(peaks, function(x) x$mode_x)
  if(length(peaks) ==0)
    return(NULL)
  
  w = which.min(abs(rt.pk-rt))
  peaks[[w]]
}



#' plots and XIC and detects the peak in it
#'
#' @param ms1_mz numeric ms1 mz
#' @param ms2_mz numeric ms2 mz
#' @param ppm numeric ppm tolerance for extraciton
#' @param rt numeric retention time in minutes
#' @param rt_tol numeric retention time window in minutes
#' @param smooth boolean T/F to apply a smoothing function
#' @param sp numeric the spar for the smoother, typically 0.12
#' @param sampleDir character vector of the sample directory containing the mzML files
#' @param title character vector of the title to use in the plot
#'
#' @export
#'
plotPeak <- function(ms1_mz,ms2_mz,ppm,rt,rt_tol,smooth = FALSE, sp = NULL, sampleDir=NULL, file = NULL){
  require(mzR)
  
  scandef <- read.delim(scandef_file)
  
  windows <- getSIMWindows(scandef)
  w = which(windows$simEnd > ms1_mz & windows$simStart < ms1_mz)
  if(length(w)>1){
    
    w <- w[1] # get the window with the highest sim range
    
  }
  
  
  lf <- list.files(path = sampleDir, pattern = "mzML", full.names = TRUE, recursive = TRUE)
  
  g <- grep(paste(windows$simStart[w],windows$simEnd[w],sep="_") ,lf)
  lf <- lf[g]
  
  g <- grep(paste("WSIM",windows$simStart[w],windows$simEnd[w],sep="_"),lf)
  ms1_mzml_file <- lf[g] 
  ms2_mzml_file <-lf[-g]
  
  ms1_raw_data <- mzR::openMSfile(ms1_mzml_file)
  ms2_raw_data <- mzR::openMSfile(ms2_mzml_file)
  
  ms1_spectra <- mzR::spectra(ms1_raw_data)
  ms2_spectra <- mzR::spectra(ms2_raw_data)
  
  ms1_header <- mzR::header(ms1_raw_data)
  ms2_header <- mzR::header(ms2_raw_data) 
  
  ms1_XIC <- getXIC(spectra = ms1_spectra,
                    header = ms1_header,
                    mz = ms1_mz,ppm=ppm,rt=rt,rt_tol=rt_tol, smooth = smooth, sp = sp)
  ms2_XIC <- getXIC(spectra= ms2_spectra,
                    header = ms2_header,
                    mz = ms2_mz,ppm=ppm,rt=rt,rt_tol=rt_tol, smooth = smooth, sp = sp)  
  
  
  pk_ms1 = detectPeak(ms1_XIC,rt)
  pk_ms2 = detectPeak(ms2_XIC,rt)
  
  if(is.null(pk_ms1) | is.null(pk_ms2)){
    return(NULL)
    #return(0)
  }
  
  max.int = c(pk_ms1$max.intensity,pk_ms2$max.intensity)
  
  max.int = max.int[which.max(max.int)]
  main.txt = paste(round(ms1_mz,4),"@",round(rt,2),sep="")
  
  
  if(!is.null(file)){
    pdf(width = 4,height = 4,file = file, pointsize = 8)
  }
  
  par(mar = c(4,4,1,1))
  plot(pk_ms1$XIC, type="l", col=4, lwd=2, ylim = c(0,max.int*1.1), bty="n",
       main = main.txt,
       xlab = expression(paste(italic("t")["R"]," (min)")), ylab = "", las= 2)
  lines(pk_ms2$XIC, type="l", col=2, lwd=2)
  
  lines(pk_ms2$peak,type = "h", lwd = 2,col=rgb(0.9,0.0,0.0,0.23))
  lines(pk_ms1$peak,type = "h", lwd = 2,col=rgb(0.1,0.1,0.9,0.23))
  
  sp2 = spline(pk_ms2$peak, n = 100)
  sp1 = spline(pk_ms1$peak, n = 100)
  
  y1 = sp1$y
  y2 = sp2$y
  
  cr = cor(y2,y1, method = "spearman")
  
  text(min(pk_ms1$XIC$x )+rt_tol*0.1,max.int*0.8,
       paste("Area MS1: ",round(pk_ms1$area),"\n",
             "Area MS2: ", round(pk_ms2$area),"\n",
             "RA: ", round(pk_ms2$area/pk_ms1$area,2),"\n",
             "RI: ", round(pk_ms2$max.intensity/pk_ms1$max.intensity,2),"\n",
             "cor: ", round(cr,2),sep="")
       ,cex = 1,pos= 4)
  
  if(!is.null(file)){
    dev.off()
  }
  
}
