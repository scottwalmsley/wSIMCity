#' search_filter
#'
#' @param filter character vector of the Orbitrap filter string
#' @param type 'SIM' for true wide SIM DIA data or any other value for pseudo MS1 0 volt data \(true MS2\)
#'
#' @return data.frame contraining segment search results
#' @export
#'
search_filter <- function(filter,type = 'SIM'){

   peak.idx <- dat <- dRT <- RT <- dM <- dMppm <- dint <- cMZ <- pMZ <- tint <- pint <- scan <- filters <-  NULL


   if(length(filterString[-g]) != 0){
      cat(paste("Processing acquisiton window #",filteridx,"/",length(filterString[-g]),", for filter ID: ",filter, "\r", sep =""))
   }else{
      cat(paste("Processing acquisiton window #",filteridx,"/",length(filterString),", for filter ID: ",filter, "\r", sep =""))
   }
   # Get scans from selected filter
   m <- match(header$filterString,filter)

   idx <- which(!is.na(m))

   rt <- header$retentionTime[which(!is.na(m))]

   # Iterate over the scans
   rtidx <- 1

   for(scanindex in idx){


      ms1_spectrum = spectra[[scanindex]]
      ms1_spectrum = 	ms1_spectrum[which(ms1_spectrum[,1] > scandef[scandef_idx,4] & ms1_spectrum[,1] < scandef[scandef_idx,5]),]

      if(type == 'SIM'){
         ms2_spectrum = spectra[[scanindex+1]]
      }else{
         ms2_spectrum = spectra[[scanindex]]
      }
      ms2_spectrum = ms2_spectrum[which(ms2_spectrum[,1] > mzmin +delta_search_mass & ms2_spectrum[,1] < scandef[scandef_idx,5]+delta_search_mass),] # 229 minimum feasable mz for nucleoside - adduct

      if(class(ms1_spectrum)== "matrix" & class(ms2_spectrum)== "matrix"){

         if( nrow(ms1_spectrum)>0 & nrow(ms2_spectrum)>0){

            # defines a peak index array
            pk.idx = seq(1,length.out = nrow(ms1_spectrum))

            # get the theoretical neutral loss mz for each ms1 mz
            th_nl_mz = round(ms1_spectrum[,1]+delta_search_mass,4)  #l = ms1

            # get each parent mz for the scan
            p_mz = round(ms1_spectrum[,1],4)

            # parent intensity
            p_int = ms1_spectrum[,2]

            #MS2 intensities
            th_nl_int = ms2_spectrum[,2]

            #Vectors to hold the values
            ms1scan = ms2scan = rtmat = pmat = pintmat = array(dim = nrow(ms2_spectrum)) #l - ms1

            # matrices to hold the values
            mat = mzmat = ppmmat = intmat = dintmat = tintmat = matrix(nrow = nrow(ms2_spectrum), ncol = length(th_nl_mz))

            # iterate over each neutral loss mz (theoretical)  dim = nrow ms1_spectrum
            for(i in 1:length(th_nl_mz)){  #MS1 length

               #i=397 #541 = ms2 match

               ms1scan[i] = scanindex
               if(type == 'SIM'){
                  ms2scan[i] = scanindex + 1
               }else{
                  ms2scan[i] = scanindex
               }
               rtmat[i] = 	rt[rtidx]  #

               #dRTmat[i] =

               mzmat[,i] = round(ms2_spectrum[,1],4)

               #dRTmat[i]
               pmat[i] = round(p_mz[i],4)#

               mat[,i] = round(ms2_spectrum[,1],4) - round(th_nl_mz[i],4)
               ppmmat[,i] = round(mat[,i]/th_nl_mz[i] * 1e6,2)

               intmat[,i] = ms2_spectrum[,2]
               dintmat[,i] = ms2_spectrum[,2]/th_nl_int[i]

               #tol = 20

               w = which(abs(ppmmat[,i]) < ppm_tol) #541

               if(length(w)>0){

                  scan = c(scan, rep(ms1scan[i], times = length(w))) #


                  RT = c(RT,rep(rt[rtidx], times = length(w)))   #
                  #dRT = c(dMppm,dRTmat[w,i]) #

                  cMZ = c(cMZ,mzmat[w,i])                          #
                  pMZ = c(pMZ,rep(p_mz[i], times = length(w)))     #

                  dM = c(dM,round(mat[w,i],4)) #
                  dMppm = c(dMppm,ppmmat[w,i]) #


                  dint = c(dint,dintmat[w,i]) #
                  tint = c(tint,intmat[w,i])  #

                  pint = c(pint, rep(p_int[i], times = length(w)))
                  peak.idx = c(peak.idx, rep(pk.idx[i],times = length(w)))
                  filters = c(filters,rep(filter,times = length(w)))


               }
            }
         }
      }
      rtidx=rtidx+1


   }

   assign('filteridx', filteridx+1,.GlobalEnv) #needed?
   assign('scandef_idx', scandef_idx+2,.GlobalEnv)

   dat = data.frame(peak.idx,'filter' = filters, scan, 'rt_min' = round(RT/60 ,4), 'mz' = pMZ, 'nlMZ' = cMZ, dM,'ppm'= dMppm,'p_int' = pint, 'nl_int' =tint)

}
