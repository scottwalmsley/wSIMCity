#' Find isotopes in psuedo spectra
#'
#' @param isomatrix the isomatrix
#' @param dat the data frame
#' @param ppmtol the ppm tolerance
#' @param maxiso max # isotopologues
#' @param maxcharge max allowable charge
#' @param mzabs max mz deviation in  the search
#' @param minfrac min sample fraction (not used)
#'
#' @return matrix
#' @export
findIsotopesPspec <- function(isomatrix,dat,  ppmtol = 15, maxiso = 4, maxcharge = 2, mzabs = 0.005 , minfrac = 0.5){
   ## isomatrix - isotope annotations (5 column matrix)
   ## mz - m/z vector, contains all m/z values from specific pseudospectrum
   ## int - int vector, see above
   ## maxiso - how many isotopic peaks are allowed
   ## maxcharge - maximum allowed charge
   ## devppm - scaled ppm error
   ## mzabs - absolut error in m/z
   
   ## matrix with all important informationen
   #isomatrix <- matrix(ncol=5, nrow=0);
   calc.isomatrix = calcIsotopeMatrix(maxiso)
   filter = FALSE
   
   
   #con <- RSQLite::dbConnect(RSQLite::SQLite(),db_name)
   #dat <- 	RSQLite::dbGetQuery(con, "SELECT * FROM peak_group_data where pspec == 5")
   #RSQLite::dbDisconnect(con)

   spectra <- matrix(c(dat$mn_mz_ms1, dat$pk_group), ncol=2)
   int     <- as.matrix(dat$int_ms1)
   int     <- int[order(spectra[, 1]), , drop=FALSE]
   mz <- as.matrix( dat$mn_mz_ms1)
   mz <- mz[order(spectra[, 1]), , drop=FALSE]
   
   spectra <- spectra[order(spectra[, 1]), ];    
   cnt     <- nrow(spectra);
   
   ## calculate error
   devppm = ppmtol / 1000000
   error.ppm <- devppm * mz;
   
   ## for every peak in pseudospectrum
   for ( j in 1:(length(mz) - 1)){
      ## create distance matrix
      MI <- spectra[j:cnt, 1] - spectra[j, 1];
      ## Sum up all possible/allowed isotope distances + error(ppm of peak mz and mzabs)
      max.index <- max(which(MI < (sum(calc.isomatrix[1:maxiso, "mzmax"]) + error.ppm[j] + mzabs )))
      ## check if one peaks falls into isotope window
      if(max.index == 1){
         ## no promising candidate found, move on
         next;
      }
      
      ## IM - isotope matrix (column diffs(min,max) per charge, row num. isotope)
      ## STN: isn't it rather:
      ## IM - isotope matrix (column diffs(min,max) per ISOTOPE, row num. charge state)
      IM <- t(sapply(1:maxcharge,function(x){
         mzmin <- (calc.isomatrix[, "mzmin"]) / x;
         mzmax <- (calc.isomatrix[, "mzmax"]) / x;
         error <- (error.ppm[j]+mzabs) / x
         res   <- c(0,0);
         for(k in 1:length(mzmin)){
            res <- c(res, mzmin[k]+res[2*k-1], mzmax[k]+res[2*k])
         }
         res[seq(1,length(res),by=2)] <- res[seq(1,length(res),by=2)]-error
         res[seq(2,length(res),by=2)] <- res[seq(2,length(res),by=2)]+error
         return (res[-c(1:2)])
      } ))
      
      ## Sort IM to fix bug, with high ppm and mzabs values 
      ## TODO: Find better solution and give feedback to user!
      IM <- t(apply(IM,1,sort))
      
      ## find peaks, which m/z value is in isotope interval
      hits <- t(apply(IM, 1, function(x){ findInterval(MI[1:max.index], x)}))
      rownames(hits) <- c(1:nrow(hits))
      colnames(hits) <- c(1:ncol(hits))
      hits[which(hits==0)] <-NA
      hits <- hits[, -1, drop=FALSE]
      hits.iso <- hits%/%2 + 1;
      
      
      ## check occurence of first isotopic peak
      for(iso in 1:min(maxiso,ncol(hits.iso))){
         hit <- apply(hits.iso,1, function(x) any(naOmit(x)==iso))
         hit[which(is.na(hit))] <- TRUE
         if(all(hit))
            break;
         hits.iso[!hit,] <- t(apply(hits.iso[!hit,,drop=FALSE],1, function(x) {
            if(!all(is.na(x))){
               ini <- which(x > iso)
               
               ## Here the following condition was previously:
               ## if(!is.infinite(ini) && length(ini) > 0){
               ##
               ## The fix for issue #44 assumes the follwoing:
               ## "There is at least one hit" Not sure why
               ## ini as return value of which() would contain inf at all
               
               if(!is.infinite(ini)[1] & length(ini) > 0){
                  x[min(ini):ncol(hits.iso)] <- NA  
               }
            }
            x
         }))
      }
      
      ## set NA to 0
      hits[which(is.na(hits.iso))] <- 0
      ## check if any isotope is found
      hit <- apply(hits, 1, function(x) sum(x)>0)
      ## drop nonhits  
      hits <- hits[hit, , drop=FALSE]
      
      ## if no first isotopic peaks exists, next
      if(nrow(hits) == 0){
         next;
      }
      
      
      ## getting max. isotope cluster length
      ## TODO: unique or not????
      isohits   <- lapply(1:nrow(hits), function(x) which(hits[x, ] %% 2 !=0))
      isolength <- sapply(isohits, length)
      
      ## Check if any result is found
      if(all(isolength==0)){
         next;
      }
      
      ## itensity checks
      ## candidate.matrix
      ## first column - how often succeded the isotope intensity test
      ## second column - how often could a isotope int test be performed
      candidate.matrix <- matrix(0, nrow=length(isohits), ncol=max(isolength)*2);
      
      for(iso in 1:length(isohits)){
         for(candidate in 1:length(isohits[[iso]])){
            for(sample.index in c(1:ncol(int))){
               charge <- as.numeric(row.names(hits)[iso])
               int.c12 <- int[j, sample.index]
               isotopePeak <- hits[iso,isohits[[iso]][candidate]]%/%2 + 1;
               if(isotopePeak == 1){
                  ## first isotopic peak, check C13 rule
                  int.c13 <- int[isohits[[iso]][candidate]+j, sample.index];
                  int.available <- all(!is.na(c(int.c12, int.c13)))
                  if (int.available){
                     theo.mass <- spectra[j, 1] * charge; #theoretical mass
                     numC      <- abs(round(theo.mass / 12)); #max. number of C in molecule
                     inten.max <- int.c12 * numC * 0.011; #highest possible intensity
                     inten.min <- int.c12 * 1    * 0.011; #lowest possible intensity
                     if((int.c13 < inten.max && int.c13 > inten.min) || !filter){
                        candidate.matrix[iso,candidate * 2 - 1] <- candidate.matrix[iso,candidate * 2 - 1] + 1
                        candidate.matrix[iso,candidate * 2 ] <- candidate.matrix[iso,candidate * 2] + 1
                     }else{
                        candidate.matrix[iso,candidate * 2 ] <- candidate.matrix[iso,candidate * 2] + 1
                     }
                  } else {
                     ## todo
                  } 
               } else {
                  ## x isotopic peak
                  int.cx <- int[isohits[[iso]][candidate]+j, sample.index];
                  int.available <- all(!is.na(c(int.c12, int.cx)))
                  if (int.available) {
                     intrange <- c((int.c12 * calc.isomatrix[isotopePeak,"intmin"]/100),
                                   (int.c12 * calc.isomatrix[isotopePeak,"intmax"]/100))
                     ## filter Cx isotopic peaks muss be smaller than c12
                     if(int.cx < intrange[2] && int.cx > intrange[1]){
                        candidate.matrix[iso,candidate * 2 - 1] <- candidate.matrix[iso,candidate * 2 - 1] + 1
                        candidate.matrix[iso,candidate * 2 ] <- candidate.matrix[iso,candidate * 2] + 1                        
                     }else{
                        candidate.matrix[iso,candidate * 2 ] <- candidate.matrix[iso,candidate * 2] + 1
                     }
                  } else {
                     candidate.matrix[iso,candidate * 2 ] <- candidate.matrix[iso,candidate * 2] + 1
                  }#end int.available
               }#end if first isotopic peak
            }#for loop samples
         }#for loop candidate
      }#for loop isohits
      
      ## calculate ratios
      candidate.ratio <- candidate.matrix[, seq(from=1, to=ncol(candidate.matrix),
                                                by=2)] / candidate.matrix[, seq(from=2, 
                                                                                to=ncol(candidate.matrix), by=2)];
      if(is.null(dim(candidate.ratio))){
         candidate.ratio <- matrix(candidate.ratio, nrow=nrow(candidate.matrix))
      }
      if(any(is.nan(candidate.ratio))){
         candidate.ratio[which(is.nan(candidate.ratio))] <- 0;
      }
      
      ## decision between multiple charges or peaks
      for(charge in 1:nrow(candidate.matrix)){
         if(any(duplicated(hits[charge, isohits[[charge]]]))){
            ## One isotope peaks has more than one candidate
            ## check if problem is still consistent
            for(iso in unique(hits[charge, isohits[[charge]]])){
               if(length(index <- which(hits[charge, isohits[[charge]]]==iso))== 1){
                  ## now duplicates next
                  next;
               }else{
                  ## find best
                  index2 <- which.max(candidate.ratio[charge, index]);
                  save.ratio <- candidate.ratio[charge, index[index2]]
                  candidate.ratio[charge,index] <- 0
                  candidate.ratio[charge,index[index2]] <- save.ratio
                  index <- index[-index2]
                  isohits[[charge]] <- isohits[[charge]][-index]
               }
            }
         }#end if
         
         for(isotope in 1:ncol(candidate.ratio)){
            if(candidate.ratio[charge, isotope] >= minfrac){
               isomatrix <- rbind(isomatrix, 
                                  c(spectra[j, 2],
                                    spectra[isohits[[charge]][isotope]+j, 2], 
                                    isotope, as.numeric(row.names(hits)[charge]), 0))
            } else{
               break;
            }
         }
      }# for(charge in 1:nrow(candidate.matrix)){
   }# end for ( j in 1:(length(mz) - 1)){
   
   return(isomatrix)
   
   
   
}
