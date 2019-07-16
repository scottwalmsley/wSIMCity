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
getXIC <- function(spectra,header,mz,ppm=NULL,mzdec=NULL, rt,rt_tol,rt_offset = NULL, smooth = TRUE,sp = NULL){
  set.seed(123456)
  ## Filter scans by RT store for extraction of mz (likely more efficient)
  
  rtRange <- getRTRangeByRT(rt*60,rt_tol*60)   #+/- 1.5 min, mzR reads in seconds
  
  w <- which(header$retentionTime > rtRange[1] & header$retentionTime < rtRange[2])
  
  subsetOfSpectra <- spectra[w]
  
  rtVals <- header$retentionTime[w]
  ## Now extract the range of needed mz and intensities
  
  if(!is.null(ppm)){
    mzRange  <- getMassTolRange(mz,ppm)  # I like 6ppm in hi res orbitrap data that I'm working with
  }
  if(!is.null(mzdec)){
    mzRange <- c(mz-mzdec, mz+mzdec)
  }
  
  
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




#' Extract the XIC from raw data
#'
#' @param spectra object mzR spectra
#' @param header object msR header infor
#' @param mz numeric desired mz to extract
#' @param ppm numeric tolerance window in ppm (7 is standard in our system)
#' @param mzdec 
#' @param rtmin numeric minimum rt
#' @param rtmax numeric maximum rt
#' @param smooth boolean 
#' @param sp numeric span value
#'
#' @return matrix
#' @export
getRawXIC <- function(spectra,header,mz,ppm=NULL,mzdec=NULL, rtmin,rtmax,smooth = FALSE,sp = NULL){
  
  # Filter scans by RT store for extraction of mz (likely more efficient)
  w <- which(header$retentionTime > rtmin*60 & header$retentionTime < rtmax*60)
  
  subsetOfSpectra <- spectra[w]
  
  rtVals <- header$retentionTime[w]
  
  ## Now extract the range of needed mz and intensities
  if(!is.null(ppm)){
    mzRange  <- getMassTolRange(mz,ppm)  # I like 6-7ppm in hi res orbitrap data that I'm working with
  }
  if(!is.null(mzdec)){
    mzRange <- c(mz-mzdec, mz+mzdec)
  }
  
  
  i<-1## index for retention time
  out <- matrix(nrow = length(w), ncol = 3)
  out[,2] = rtVals/60
  out[,3] = 0
  out[,1] = NA
  
  for(spectrum in subsetOfSpectra){  
    
    w <- which(spectrum[,1] > mzRange[1] & spectrum[,1] < mzRange[2])
    
    if(length(w) > 0){
      
      spectrum <- spectrum[w,]
      
      if(length(w) > 1){
        
        # Handles multiple mz within tol, in this case extract the closest by mass.   
        # Hopefully intensity matches expected modality
        w.min <- which.min(abs(spectrum[,1]-mz))
        w.max <- which.max(abs(spectrum[,1]-mz))
        mz.min <- spectrum[,]
        spectrum <- spectrum[w.min,]  
        
        
      }
      out[i,1] <- spectrum[1]
      out[i,3] <- spectrum[2]
      
    }
    
    i<-i+1 
  }
  
  mat <- out
  
  if(smooth){
    
    smooth_data <- smooth.spline(out[,2],out[,3],spar = sp)
    #n = spline(smooth)
    
    mat <- as.matrix(cbind(mat[,1],smooth_data$x, smooth_data$y))
    
    mat[which(mat[,3]<0),3] <- 0
    
  }
  
  mat[which(is.na(mat[,1])),1] = mean(mat[,1],na.rm=T)
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
#' @param scandef data.frame containing the SIM windows
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
#' @param RT numeric retention time
#' @param tol numeric tolerance
#'
#' @return vector containing limits of RT range
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
##@examples
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
  
  n = spline(XIC[,2],XIC[,3])
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





#' Plot peaks from raw data
#'
#' @param ms1_mz 
#' @param ms2_mz 
#' @param ppm 
#' @param rtmin 
#' @param rtmax 
#' @param smooth 
#' @param sp 
#' @param sampleDir 
#' @param file 
#' 
#' @export
#'
# @examples
plotRAWPeak <- function(ms1_mz,ms2_mz,ppm,rtmin,rtmax,smooth = FALSE, sp = NULL, sampleDir=NULL, file = NULL){
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
  
  ms1_XIC <- getRAWXIC(spectra = ms1_spectra,
                       header = ms1_header,
                       mz = ms1_mz,ppm=ppm,rtmin=rtmin,rtmax=rtmax, smooth = smooth, sp = sp)
  ms2_XIC <- getRAWXIC(spectra= ms2_spectra,
                    header = ms2_header,
                    mz = ms2_mz,ppm=ppm,rtmin=rtmin,rtmax=rtmax, smooth = smooth, sp = sp)  
  
  
  #pk_ms1 = detectPeak(ms1_XIC,rt)
  #pk_ms2 = detectPeak(ms2_XIC,rt)
  
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
  plot(ms1_XIC, type="l", col=4, lwd=2, ylim = c(0,max.int*1.1), bty="n",
       main = main.txt,
       xlab = expression(paste(italic("t")["R"]," (min)")), ylab = "", las= 2)
  lines(ms2_XIC, type="l", col=2, lwd=2)
  
  
  
  if(!is.null(file)){
    dev.off()
  }
  
}



#' Get the peaks from an XIC pair
#'
#' @param ms1_XIC XIC matrix from the ms1 data
#' @param ms2_XIC XIC matrix from the ms2 data
#' @param adduct_mass the adduct mz value.  Dont forget the "-" if a loss.
#' @param rt numeric retention time value
#' @param sp numeric span for the smoother
#' @param plot boolean plot or no plot
#'
#' @return dataframe
#' @export
#'
getPeaks = function(ms1_XIC,ms2_XIC,adduct_mass,rt,sp = 0.4, plot = TRUE, filename = NULL){
  set.seed(123456)
  if(colSums(ms1_XIC)[3]==0 | colSums(ms2_XIC)[3] == 0 ){
    return(NULL)
  }
  
  # Detects the peaks at the given retention time
  pk_ms1 <- detectPeakAtRT(ms1_XIC,rt)
  pk_ms2 <- detectPeakAtRT(ms2_XIC,rt)
  
  # 
  pkl_1 <- dim(pk_ms1$peak)[1]
  pkl_2 <- dim(pk_ms2$peak)[1]
  if(is.null(pkl_1 ) | is.null(pkl_2)){
    return(NULL)
  }
  
  
  
  
  # not used yet
  mx1 <- which(pk_ms1$peak[,3] == pk_ms1$max.intensity)
  mx2 <- which(pk_ms2$peak[,3] == pk_ms2$max.intensity)
  
  if(colSums(pk_ms1$peak)[2]==0 | colSums(pk_ms2$peak)[2] == 0 ){
    return(NULL)
  }
  
  # for plotting
  max.int <- max(pk_ms1$max.intensity,pk_ms2$max.intensity)
  par(mar = c(4,4,1,1))
  
  # Compute the correlations using trc:: package!
  options(warn = -1)
  trc_pk <- tryCatch(trc::trc_cor_test(pk_ms1$peak[,3],pk_ms2$peak[,3],nperm = 100)$measure, error =  function(err) NA)
  options(warn = 0)
  
  
  # compute deltas for scoring
  dRT <- diff(c(pk_ms1$mode_x,pk_ms2$mode_x))
  score_rt = mean(similarityScore_gauss(dRT,0.2)) # hard coded tolerance ~ 0.2 to ensure data are within 2-3 scans of each other, subject to change.
  
  
  dM <-  mean(pk_ms2$peak[,1]) - (mean(pk_ms1$peak[,1]) + adduct_mass)
  
  dM_ppm = dM/(mean(pk_ms1$peak[,1]) + adduct_mass) * 1e6
  score_mz = mean(similarityScore_gauss(dM,0.007))  # daltons, tolerance is 0.007, since in our system peak centroiding flickers the data within ~0.003-0.005 Da.
  
  # feature score
  score = (0.5 * score_mz) + (0.5 * score_rt)
  
  
  # recompute the ratio
  new_ratio = pk_ms2$max.intensity/pk_ms1$max.intensity
  
  if(!is.na(trc_pk[1])){
    if(plot==TRUE & score > 0.8 & new_ratio > 0.1 & new_ratio < 10 & trc_pk[1] > 0.8){
      
      svg(filename = filename)
      
      plot(pk_ms1$XIC, type="l", col=4, lwd=2, ylim = c(0,max.int*1.21), bty="n",
           xlab = expression(paste(italic("t")["R"]," (min)")), ylab = "", las= 1)
      lines(pk_ms2$XIC, type="l", col=2, lwd=2)
      
      lines(pk_ms1$peak[,2:3],type="h", col=rgb(0.1,0.1,0.8,0.2),lwd=3)
      lines(pk_ms2$peak[,2:3],type="h", col=rgb(0.8,0.1,0.1,0.2),lwd=3)
      
      
      text(min(pk_ms1$XIC$x ),max.int*0.9,
           paste("mz (wsim,nl): ",round(mean(pk_ms1$peak[,1]),4),", ",round(mean(pk_ms2$peak[,1]),4),"\n",
                 "deltas (ppm,rt): ",round(dM_ppm,2),", ",round(dRT,2),"\n", 
                 "Areas: ",round(pk_ms1$area),", ",round(pk_ms2$area),"\n",
                 "Ratios (area,int): ", round(pk_ms2$area/pk_ms1$area,2),", ",round(pk_ms2$max.intensity/pk_ms1$max.intensity,2),"\n",
                 "Cor (rho): ",round(trc_pk[1],2),"\n",
                 "Scores (feat,final): ",round(score,3),", ",round(score * trc_pk[1],3),
                 sep="")
           ,cex = 0.75,pos= 4)
      
      dev.off()
      
      png(filename = sub(".svg",".png",filename))
      
      plot(pk_ms1$XIC, type="l", col=4, lwd=2, ylim = c(0,max.int*1.21), bty="n",
           xlab = expression(paste(italic("t")["R"]," (min)")), ylab = "", las= 1)
      lines(pk_ms2$XIC, type="l", col=2, lwd=2)
      
      lines(pk_ms1$peak[,2:3],type="h", col=rgb(0.1,0.1,0.8,0.2),lwd=3)
      lines(pk_ms2$peak[,2:3],type="h", col=rgb(0.8,0.1,0.1,0.2),lwd=3)
      text(min(pk_ms1$XIC$x ),max.int*0.9,
           paste("mz (wsim,nl): ",round(mean(pk_ms1$peak[,1]),4),", ",round(mean(pk_ms2$peak[,1]),4),"\n",
                 "deltas (ppm,rt): ",round(dM_ppm,2),", ",round(dRT,2),"\n", 
                 "Areas: ",round(pk_ms1$area),", ",round(pk_ms2$area),"\n",
                 "Ratios (area,int): ", round(pk_ms2$area/pk_ms1$area,2),", ",round(pk_ms2$max.intensity/pk_ms1$max.intensity,2),"\n",
                 "Cor (rho): ",round(trc_pk[1],2),"\n",
                 "Scores (feat,final): ",round(score,3),", ",round(score * trc_pk[1],3),
                 sep="")
           ,cex = 0.75,pos= 4)
      dev.off()
      
    }
  
  }
  data.frame("ms1_intensity" = pk_ms1$max.intensity, 
             "ms2_intensity" = pk_ms2$max.intensity, 
             "ms1.mz" = mean(pk_ms1$peak[,1]),
             "ms2.mz" = mean(pk_ms2$peak[,1]),
             "ms1.RT" =  mean(pk_ms1$mode_x),
             "ms2.RT" =  mean(pk_ms2$mode_x),
             "recomputed_ratio" = new_ratio,
             "corr" = trc_pk[1],
             "score_mz" = score_mz,
             "score_rt" = score_rt,
             "score" = score,
             "final_score" = score * trc_pk[1])
  
}


#' Extract XIC data and compute correlations on peaks.
#'
#' @param sampleDir character vector pointing to the sample directory containing the peak list.
#' @param scandef_file character vector pointing to the scan definition file.
#'
#' @return list of dataframes containing the correlation and scores for precursor NL pairs.
#' @export
#'
getXIC_cor <- function(sampleDir, scandef_file, plot = TRUE){
  
  options(warn = -1)
  
  cmd <- paste("del ",sampleDir,"/plots/*", sep = "")
  
  system(cmd)
  
  cat(paste("Processing", sampleDir,"\n"))
  
  set.seed(123456)
  
  # read the search results to obtain the list of precursor ions.
  # for your case, you can alternative ly read the csv file  you produce and then obtain the mz and rt values.
  lf <- list.files(path = sampleDir, pattern = "filtered_SearchResults", recursive = TRUE, full.names = TRUE)
  
  
  d = read.csv(file = lf[1], header=T, row.names = 1)
  d <- d[order(d$mz_MS1),]
  
  # I have a minimum filter in place.
  d<- d[which(d$mz_MS1>330 &  d$mz_MS1 < 634),]
  
  
  # The scan file is a necessary evil to get around lack of header info.  
  scandef = read.delim(scandef_file)
  
  
  # Gets the scan window information in order to extract the correct data from the files. See the segmentmzml function.
  win = getSIMWindows(scandef)
  win_min = win$simStart
  win_max = win$simEnd
  win_srch = win$srch
  
  print("Loading mzML raw data for computing correlations....")
  
  # Obtain the correct mzML files
  lf = list.files(path = sampleDir, pattern = "mzML", recursive = TRUE, full.names = TRUE)
  
  # mzML files are segregated into ms1 and ms2 data. In the case of your 0-25 method, this works for those scan types as well.
  lf_1 = grep("/WSIM",lf)
  lf_2 = grep("/NL",lf)
  
  # open the file connections and pre-load them into memory for later access.
  raw_ms1 <- lapply(lf[lf_1],function(x) mzR::openMSfile(x))
  raw_ms2 <- lapply(lf[lf_2],function(x) mzR::openMSfile(x))
  
  # Assign peaks to scan ranges
  file_idx <- NULL#unlist(lapply(d$mz_MS1,function(x) {w = which( win_min < x & win_max > x); if(length(w)>1){w = w[1]};w}))
  
  # Index each precuror tothe correct scan range mzml file set.
  for(mz in d$mz_MS1){
    
    w = which( win_min < mz & win_max > mz)
    file_idx <- c(file_idx,w[1])
    #cat(paste(w,win_min[w[1]],win_max[w[1]]),mz,"\n")
  }
  
  
  # holds the results
  corr = list()
  
  # for each scan range and for each mz, get the raw data XICs (ms1, and ms2)
  for(idx in unique(file_idx)){
    
    # Load spectra
    spectra_ms1 <- mzR::spectra(raw_ms1[[idx]])
    spectra_ms2 <- mzR::spectra(raw_ms2[[idx]])
    
    # get header infor
    header_ms1 <- mzR::header(raw_ms1[[idx]])
    header_ms2 <- mzR::header(raw_ms2[[idx]])
    
    # w_idx: which file to open
    w_idx <- which(file_idx == idx)
    
    
    # Subset the data from the original data file.
    ms1_RT <- d$RT_MS1[w_idx]
    ms1_mz <- d$mz_MS1[w_idx]
    ms2_mz <- d$mz_MS2[w_idx]
    ratio <- d$ratio_int[w_idx]
    adduct <- as.character(d$adduct[w_idx])
    adduct_mass <- d$adduct_mass[w_idx]
    
    # Hard code this at ~0.3, seems to avoid false negatives. 
    sp = 0.4
    
    # extract ms1 data for all the peaks in the list indexed to the mzml files.
    raw_XIC_ms1 = lapply(1:length(w_idx), function(x) wSIMCity::getRawXIC(spectra = spectra_ms1,
                                                                          header = header_ms1,
                                                                          mz = ms1_mz[x], 
                                                                          ppm = 7, 
                                                                          rtmin = ms1_RT[x]-1.5, 
                                                                          rtmax = ms1_RT[x]+1.5, 
                                                                          smooth = TRUE, 
                                                                          sp = sp))
    
    # extract ms2 data....
    raw_XIC_ms2 = lapply(1:length(w_idx), function(x) wSIMCity::getRawXIC(spectra = spectra_ms2,
                                                                          header = header_ms2,
                                                                          mz = ms2_mz[x], 
                                                                          ppm = 7, 
                                                                          rtmin = ms1_RT[x]-1.5, 
                                                                          rtmax = ms1_RT[x]+1.5,  
                                                                          smooth = TRUE, 
                                                                          sp = sp))
    
    
    # detect the peaks in the data and then return the list of dataframe data
    for( k in (1:length(w_idx))){
      
      if(plot){
        filename <- paste(sampleDir,"/plots/",adduct[k],"_", round(ms1_mz[k],4),"@",round(ms1_RT[k],3),".svg",sep="") 
      }
      #cat(paste("idx:",idx,"k:",k,ms1_mz[k],adduct_mass[k]),"\n")
      corr[[k]] =   
        getPeaks(	raw_XIC_ms1[[k]],raw_XIC_ms2[[k]],adduct_mass[k] ,rt = ms1_RT[k], plot = plot,file =  filename)
      
    }
    
    
    
  }
  outfile = NULL
  for(l in corr){
    if(!is.null(l)){
      outfile <- rbind(outfile,l)
    }
  }
  
  
  fh <- paste(sampleDir,"/finalSearchResults.csv",sep="")
  write.csv(file = fh,outfile)
  options(warn = 0)
  corr
}



#' Detect peaks in an XIC
#'
#' @param XIC matrix XIC of 3 columns (zmz, rt, intensity)
#' @param rt numeric target retention time
#'
#' @return list
#' @export
#'
detectPeakAtRT <- function(XIC,rt){
  set.seed(123456)
  # The smoother
  n = spline(XIC[,2],XIC[,3])
  n$y[which(n$y<0)] = 0  
  
  # need imputed mz values
  nmz = spline(XIC[,2],XIC[,1])
  nmz$y[which(n$y<0)] = 0  
  
  
  #Temp function to aid background subtractions, will migrate from pracma package later.
  #pv = quantmod::findValleys(n$y)-1
  
  # Background subtraction
  #bk = spline(n$x[pv],n$y[pv], n = length(n$x),xmin = min(n$x), xmax = max(n$x))
  #n$y = n$y-bk$y
  #n$y[which(n$y<0)] = 0  
  
  
  # Peak selection and boundary recognition
  pk = pracma::findpeaks(n$y,threshold = 500)
  if(is.null(pk))
    return(NULL)
  
  
  # Select the peak at the indicated RT determined using the feature finding software.
  w = which.min(abs(n$x-rt))
  
  # Get peak boundaries
  pk_idx = pk[which.min(abs(pk[,2]-w)),2]
  pk_start = pk[which(pk[,2]==pk_idx),3]
  pk_end = pk[which(pk[,2]==pk_idx),4]
  
  #cat(paste(rt,pk_idx,pk_start,pk_end,"\n"))
  # 3col matrix: mz, rt, intensity  
  
  
  pk <- as.matrix(cbind(nmz$y[pk_start:pk_end],n$x[pk_start:pk_end],n$y[pk_start:pk_end]))
  
  
  
  
  # max intensityu
  max.int <- pk[which.max(pk[,3]),3]
  
  peaks = list("XIC" = n , 
               "area" =sum(n$y[pk_start:pk_end],na.rm = T), 
               max.intensity = max.int, 
               "mode_x" = n$x[which.max(n$y)], # the retention time
               "rt_start"= n$x[pk_start],
               "rt_end" = n$x[pk_end],
               "peak" = pk)
  
  
  if(length(peaks) ==0)
    return(NULL)
  
  peaks
}



