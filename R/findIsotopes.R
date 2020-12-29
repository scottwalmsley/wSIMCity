
#' Create a matrix for isotope abundanc computation and id
#'
#' @param maxiso max number of isotopes
#'
#' @return matrix
#' @export
calcIsotopeMatrix <- function(maxiso=4){
   
   if(!is.numeric(maxiso)){
      stop("Parameter maxiso is not numeric!\n")  
   } else if(maxiso < 1 | maxiso > 8){
      stop(paste("Parameter maxiso must between 1 and 8. ",
                 "Otherwise use your own IsotopeMatrix.\n"),sep="")
   }
   
   isotopeMatrix <- matrix(NA, 8, 4);
   colnames(isotopeMatrix) <- c("mzmin", "mzmax", "intmin", "intmax")
   
   isotopeMatrix[1, ] <- c(1.000, 1.0040, 1.0, 150)
   isotopeMatrix[2, ] <- c(0.997, 1.0040, 0.01, 200)
   isotopeMatrix[3, ] <- c(1.000, 1.0040, 0.001, 200)
   isotopeMatrix[4, ] <- c(1.000, 1.0040, 0.0001, 200)
   isotopeMatrix[5, ] <- c(1.000, 1.0040, 0.00001, 200)
   isotopeMatrix[6, ] <- c(1.000, 1.0040, 0.000001, 200)
   isotopeMatrix[7, ] <- c(1.000, 1.0040, 0.0000001, 200)
   isotopeMatrix[8, ] <- c(1.000, 1.0040, 0.00000001, 200)  
   
   return(isotopeMatrix[1:maxiso, , drop=FALSE])
   
}

#' Find isotopic peaks in peak data
#'
#' @param mslevel the ms level to search
#' @param maxiso the maximum number of isotopes
#'
#' @export
findIsotopes <- function(mslevel = 1, maxiso = 3){
   
   
   con <- RSQLite::dbConnect(RSQLite::SQLite(),db_name)
   dat <- 	RSQLite::dbGetQuery(con, "SELECT * FROM peak_group_data ")
   RSQLite::dbDisconnect(con)
   
   npeaks.global <- 0; #Counter for % bar
   npspectra <- length(unique(dat$pspec));
   
   ncl <- nrow(dat);
   
   
   if(mslevel == 1){
      imz  <- as.matrix(dat[, "mn_mz_ms1", drop=FALSE]);
      irt  <- as.matrix(dat[, "rt_ms1", drop=FALSE]);
      mint <- as.matrix(dat[, "int_ms1", drop=FALSE]);      
   }
   if(mslevel == 2){
      imz  <- as.matrix(dat[, "mn_mz_ms2", drop=FALSE]);
      irt  <- as.matrix(dat[, "rt_ms1", drop=FALSE]);
      mint <- as.matrix(dat[, "int_ms2", drop=FALSE]);      
   }
   
   isomatrix <- matrix(ncol=5, nrow=0);
   colnames(isomatrix) <- c("mpeak", "isopeak", "iso", "charge", "intrinsic")
   
   for( i in seq(along = unique(dat$pspec))){
      sub = dat[which(dat$pspec == i),]
      ipeak = sub$pk_group
      if(length(ipeak) > 1){
         #peak mass and intensity for pseudospectrum
         #mz  <- imz[ipeak];
         #int <- mint[ipeak, , drop=FALSE];
         isomatrix <-  findIsotopesPspec(isomatrix, sub, ppmtol = )    
         #print(dim(isomatrix))
      }
      
      
   }
   
   #clean isotopes
   if(is.null(nrow(isomatrix))) {
      isomatrix = matrix(isomatrix, byrow=F, ncol=length(isomatrix)) 
   }
   
   
   #check if every isotope has only one annotation
   if(length(idx.duplicated <- which(duplicated(isomatrix[, 2]))) > 0){
      peak.idx <- unique(isomatrix[idx.duplicated, 2]);
      for( i in 1:length(peak.idx)){
         #peak.idx has two or more annotated charge
         #select the charge with the higher cardinality
         peak <- peak.idx[i];
         peak.mono.idx <- which(isomatrix[,2] == peak)
         if(length(peak.mono.idx) < 2){
            #peak has already been deleted
            next;
         }
         peak.mono <- isomatrix[peak.mono.idx,1]
         #which charges we have
         charges.list   <- isomatrix[peak.mono.idx, 4];
         tmp <- cbind(peak.mono,charges.list);
         charges.length <- apply(tmp,1, function(x,isomatrix) { 
            length(which(isomatrix[, 1] == x[1] & isomatrix[,4] == x[2])) }, 
            isomatrix);
         idx <- which(charges.length == max(charges.length));
         if(length(idx) == 1){
            #max is unique
            isomatrix <- isomatrix[-which(isomatrix[, 1] %in% peak.mono[-idx] & isomatrix[, 4] %in% charges.list[-idx]),, drop=FALSE]
         }else{
            #select this one, which lower charge
            idx <- which.min(charges.list[idx]);
            isomatrix <- isomatrix[-which(isomatrix[, 1] %in% peak.mono[-idx] & isomatrix[, 4] %in% charges.list[-idx]),, drop=FALSE]
         }
      }
   }
   
   #check if every isotope in one isotope grp, have the same charge
   if(length(idx.duplicated <- which(duplicated(paste(isomatrix[, 1], isomatrix[, 3])))) > 0){
      #at least one pair of peakindex and number of isotopic peak is identical
      peak.idx <- unique(isomatrix[idx.duplicated,1]);
      for( i in 1:length(peak.idx)){
         #peak.idx has two or more annotated charge
         #select the charge with the higher cardinality
         peak <- peak.idx[i];
         #which charges we have
         charges.list   <- unique(isomatrix[which(isomatrix[, 1] == peak), 4]);
         #how many isotopes have been found, which this charges
         charges.length <- sapply(charges.list, function(x,isomatrix,peak) { length(which(isomatrix[, 1] == peak & isomatrix[, 4] == x)) },isomatrix,peak);
         #select the charge which the highest cardinality
         idx <- which(charges.length == max(charges.length));
         if(length(idx) == 1){
            #max is unique
            isomatrix <- isomatrix[-which(isomatrix[, 1] == peak & isomatrix[, 4] %in% charges.list[-idx]),, drop=FALSE]
         }else{
            #select this one, which lower charge
            idx <- which.min(charges.list[idx]);
            isomatrix <- isomatrix[-which(isomatrix[, 1] == peak & isomatrix[, 4] %in% charges.list[-idx]),, drop=FALSE]
         }
      }
   }
   
   #Combine isotope cluster, if they overlap
   index2remove <- c();
   
   if(length(idx.duplicated <- which(isomatrix[, 1] %in% isomatrix[, 2]))>0){
      for(i in 1:length(idx.duplicated)){
         index <-  which(isomatrix[, 2] == isomatrix[idx.duplicated[i], 1])
         index2 <- sapply(index, function(x, isomatrix) which(isomatrix[, 1] == isomatrix[x, 1] & isomatrix[,3] == 1),isomatrix)
         if(length(index2) == 0){
            index2remove <- c(index2remove,idx.duplicated[i])
         }
         max.index <- which.max(isomatrix[index,4]);
         isomatrix[idx.duplicated[i], 1] <- isomatrix[index[max.index], 1];
         isomatrix[idx.duplicated[i], 3] <- isomatrix[index[max.index], 3]+1;
      }
   }
   
   if(length(index <- which(isomatrix[,"iso"] > maxiso)) > 0){
      index2remove <- c(index2remove, index)
   }
   
   if(length(index2remove) > 0){
      isomatrix <- isomatrix[-index2remove,, drop=FALSE];
   }
   
   isomatrix <- isomatrix[order(isomatrix[,1]),,drop=FALSE]
   
   
   #isotopes 
   
   ## label peak groups in SQLdb
   
   isotopes <- matrix(ncol = 3,nrow = nrow(dat) )
   
   for( i in  unique(isomatrix[,1])){
      
      isocluster <- isomatrix[which(isomatrix[,1]==i),]
      if(length(which(isomatrix[,1] == i))>1){
         isotopes[which(dat$pk_group == isocluster[1,'mpeak']),1] = 0
         isotopes[which(dat$pk_group == isocluster[1,'mpeak']),2] = isocluster[1,'charge']
         isotopes[which(dat$pk_group == isocluster[1,'mpeak']),3] = i
      }
      if(length(which(isomatrix[,1] == i))==1){
         isotopes[which(dat$pk_group == isocluster['mpeak']),1] = 0
         isotopes[which(dat$pk_group == isocluster['mpeak']),2] = isocluster['charge']
         isotopes[which(dat$pk_group == isocluster['mpeak']),3] = i
      }
      
      
      if(!is.null(dim(isocluster))){
         for( j in 1:nrow(isocluster)){
            iso <- isocluster[j,] 
            w = which(dat$pk_group == iso['isopeak'])
            
            isotopes[w,1] = 1
            isotopes[w,2] = iso['charge']
            isotopes[w,3] = i
            
         }
      }else{
         w = which(dat$pk_group == isocluster['isopeak'])
      
         isotopes[w,1] = 1
         isotopes[w,2] = isocluster['charge']
         isotopes[w,3] = i
      }
      
   }
   
   if(mslevel==1){
      colnames(isotopes) = c('ion_type_ms1','z_ms1','pk_group_parent_ms1')
      
      
   }
   if(mslevel==2){
      colnames(isotopes) = c('ion_type_ms2','z_ms2','pk_group_parent_ms2')
      
   }
   isotopes <- as.data.frame(isotopes)
   
   con <- RSQLite::dbConnect(RSQLite::SQLite(),db_name)
   
   RSQLite::dbWriteTable(con,"peak_group_data",cbind(dat,isotopes), overwrite=TRUE)
   
   RSQLite::dbDisconnect(con)
   
}



