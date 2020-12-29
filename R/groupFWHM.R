#' Omit NA values
#'
#' @param x vector
#'
#' @return vector
#' @export
naOmit <- function(x) {
   return (x[!is.na(x)]);
}




#' Group peaks by retention time
#'
#' @param sigma number of standard deviation arround the mean (6 = 2 x 3 left and right)
#' @param perfwhm percentage of width of peak at fwhm
#'
#' @export
groupFWHM <- function( sigma=6, perfwhm=0.6) {
   # grouping after retentiontime 
   # sigma - number of standard deviation arround the mean (6 = 2 x 3 left and right)
   # perfwhm - 0.3;
   
   con <- RSQLite::dbConnect(RSQLite::SQLite(),db_name)
   peakmat <- 	RSQLite::dbGetQuery(con, "SELECT * FROM peak_group_data")
   RSQLite::dbDisconnect(con)
   
   
   sample    <- 1
   pspectra  <- list();
   psSamples <- NA;
   
   cat("Start grouping after retention time.\n")
   
   
   
   maxo    <- peakmat$int_ms2#[, 'into']; #max intensities of all peaks
   maxo    <- cbind(peakmat$pk_group,maxo);
   i=1
   while(length(maxo)> 0){
      iint   <- which.max(maxo[,2]);
      rtmed  <- peakmat[iint, "rt_ms2"]; #highest peak in whole spectra
      rt.min <- peakmat[iint, "rtmin"];
      rt.max <- peakmat[iint, "rtmax"]; #begin and end of the highest peak
      hwhm   <- ((rt.max - rt.min) / sigma * 2.35 * perfwhm) / 2; #fwhm of the highest peak
      #all other peaks whose retensiontimes are in the fwhm of the highest peak
      irt    <- which(peakmat[, 'rt_ms2'] > (rtmed - hwhm) & peakmat[, 'rt_ms2'] < (rtmed + hwhm)) 
      if(length(irt)>0){
         #if peaks are found
         idx <- maxo[irt,1];
         pspectra[[length(pspectra)+1]] <- idx#maxo[idx,1]; #create groups
         maxo <- maxo[-irt, ,drop=FALSE]; #set itensities of peaks to NA, due to not to be found in the next cycle
         peakmat <- peakmat[-irt, ,drop=FALSE];
         i=i+1
      }else{
         # break
         idx <- maxo[iint,1];
         cat("Warning: Feature ",idx," looks odd for at least one peak. Please check afterwards.\n");
         pspectra[[length(pspectra)+1]] <- idx; #create groups
         maxo       <- maxo[-iint, ,drop=FALSE]; #set itensities of peaks to NA, due to not to be found in the next cycle
         peakmat  <- peakmat[-iint, ,drop=FALSE];
         i=i+1
      }
      
   }
   con <- RSQLite::dbConnect(RSQLite::SQLite(),db_name)
   peakmat <- 	RSQLite::dbGetQuery(con, "SELECT * FROM peak_group_data")
   RSQLite::dbDisconnect(con)
   
   pspectra.idx = array(dim = nrow(peakmat))
   
   for(i  in 1:length(pspectra)){
      pk.idx = pspectra[[i]]
      pspectra.idx[match(pk.idx, peakmat$pk_group)] = i
      
   }
   
   
   con <- RSQLite::dbConnect(RSQLite::SQLite(),db_name)
   
   RSQLite::dbWriteTable(con,"peak_group_data",cbind('pspec' = pspectra.idx,peakmat), overwrite=TRUE)
   
   RSQLite::dbDisconnect(con)
   
   
   
   
}

