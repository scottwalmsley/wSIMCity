#' compute group statistics
#'
#' @export
#'
compute_group_data <- function(){
   
   con <- RSQLite::dbConnect(RSQLite::SQLite(),db_name)
   dat <- 	RSQLite::dbGetQuery(con, "SELECT * FROM assigned_peak_groups")
   RSQLite::dbDisconnect(con)
   
   mx.rt.ms1 = mx.rt.ms2 = mx.scan.ms1 = mx.scan.ms2 = mx.int.ms1 = mx.int.ms2 = mn.mz1 = mn.mz2 = mx.mz1 = mx.mz2 = c.g = n.pk = filter =  array(dim = length(unique(dat$pk_group)))
   #mz1.dev = mz2.dev = ppm = dM  = corr = score.rt = score.mz = score.total = ion_type = ion_type.ms2 = isotope_ratio.ms1 = isotope_ratio.ms2 = array(dim = length(unique(dat$pk_group)))
   mz1.dev = mz2.dev = ppm = dM  = corr = score.rt = score.mz = score.total = 
      rtmin = rtmax = fwhm1 = fwhm2 = array(dim = length(unique(dat$pk_group)))
   
   i=1
   for(g in unique(dat$pk_group)){
      #g = 637
      if(g %% 1000 == 0)
         cat(paste(g, "\r"))
      
      sub = dat[which(dat$pk_group == g),]
      
      w = which(sub$p_int == max(sub$p_int) )  # maximum precursor intensity index
      
      
      if(length(w > 0)){
         if(length(w) > 1){
            bw = which.min(abs(sub$ppm[w]))
            w = w[bw]
         }
         
         filter[i] =as.character(sub$filter[w])
         
         mx.scan.ms1[i] = sub$scan[w]
         mx.int.ms1[i] = sub$p_int[w]
         mx.rt.ms1[i] = sub$rt_min[w]
         mn.mz1[i] = round(mean(sub$mz),4)
         mx.mz1[i] = sub$mz[w]
         mz1.dev[i] = mean(abs(round( ((sub$mz - mean(sub$mz))/mean(sub$mz)*1e6),2)   )) 
         
         rtmin[i] = min(sub$rt_min, na.rm=T)
         rtmax[i] = max(sub$rt_min, na.rm=T)
         
         xmax <- sub$rt_min[sub$p_int==max(sub$p_int)]
         hm = sub$p_int[w]/2
         
         hw = which(sub$p_int > hm)
         
         #sub$rt_min[sub$rt_min < xmax & sub$rt_min[sub$rt_min < xmax]]
         
         # x1 <- sub$rt_min[sub$rt_min < xmax][which.min(abs(sub$p_int[w][sub$rt_min < xmax]-sub$p_int[w]/2))]
         # x2 <- sub$rt_min[sub$rt_min > xmax][which.min(abs(sub$p_int[w][sub$rt_min > xmax]-sub$p_int[w]/2))]
         
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
            
         
         sm.pk = getSmoothPeak(x = sub$rt_min,y = sub$p_int,bw = 0)
   
         fwhm1[i] = sm.pk$fwhm
         
         
         w = which(sub$nl_int == max(sub$nl_int))
         
         if(length(w) > 1){
            bw = which.min(abs(sub$ppm[w]))
            w = w[bw]
         }
         
         
         
         mx.scan.ms2[i] = sub$scan[w]
         mx.int.ms2[i] = sub$nl_int[w]
         mx.rt.ms2[i] = sub$rt_min[w]
         mn.mz2[i] = round(mean(sub$nlMZ),4)
         mx.mz2[i] = sub$nlMZ[w]
         mz2.dev[i] = mean(abs(round( ((sub$nlMZ - mean(sub$nlMZ))/mean(sub$nlMZ)*1e6),2)   ))
         
         
         sm.pk = getSmoothPeak(sub$rt_min,sub$nl_int,bw = 0.24)
         
         #lines(sm.pk$peak, col = 2)
         #points(sm.pk$peak$x[sm.pk$x1], sm.pk$peak$y[sm.pk$x1], col = 4, pch = 19, cex = 1.5)
         #points(sm.pk$peak$x[sm.pk$x2], sm.pk$peak$y[sm.pk$x2], col = 4, pch = 19, cex = 1.5)
         
         fwhm2[i] = sm.pk$fwhm
         
         
         
         #ion_type.ms2[i] = sub$ion_type_ms2[w]
         #isotope_ratio.ms2[i] = sub$isotope_ratio_ms2[w]
         
         c.g[i] = g
         n.pk[i] = length(which(dat$pk_group == g))
         
         corr[i] = round(cor(sub$p_int , sub$nl_int),2)
         score.mz[i] = round( similarityScore_gauss(sub$dM[w],0.005), 2)
         score.rt[i] = round( similarityScore_gauss(mx.rt.ms1[i] - mx.rt.ms2[i],0.2), 2)
         
         score.total[i] = score.mz[i] * score.rt[i]
         
         ppm[i] = sub$ppm[w]
         dM[i] = sub$dM[w]  ## TODO add dM to peak data
         #return(0)
      }
      i=i+1
      
   }
   
 
   
   
   meta = data.frame(pk_group = c.g,
                     n_pk = n.pk,
                     mz_ms1 = mx.mz1, mz_ms2 = mx.mz2,
                     dM = dM, ppm = ppm,
                     score_mz = score.mz,
                     score_rt = score.rt,
                     score_total = score.total,
                     #ion_type_ms1 = ion_type,
                     #ion_type_ms2 =ion_type.ms2,
                     #isotope_ratio_ms1 = isotope_ratio.ms1,
                     #isotope_ratio_ms2 = isotope_ratio.ms2,
                     #ratio_diff = isotope_ratio.ms1 - isotope_ratio.ms2,
                     rt_ms1 = mx.rt.ms1, rt_ms2 = mx.rt.ms2,
                     rtmin,rtmax,fwhm1,fwhm2,
                     int_ms1 = mx.int.ms1, int_ms2 = mx.int.ms2,
                     rho = corr,ratio = round(log(mx.int.ms2/mx.int.ms1),5),
                     ppm_dev_ms1 = round(mz1.dev,3), ppm_dev_ms2 = round(mz2.dev,3),
                     scan_ms1 =  mx.scan.ms1, scan_ms2 =  mx.scan.ms2,
                     mn_mz_ms1 = mn.mz1,mn_mz_ms2 = mn.mz2,
                     filter)
   
   meta = meta[order(meta$scan_ms1,meta$mn_mz_ms1),]
   meta = meta[which(meta$n_pk>0),]
   con <- RSQLite::dbConnect(RSQLite::SQLite(),db_name)
   
   RSQLite::dbWriteTable(con,"peak_group_data",meta, overwrite = TRUE)
   RSQLite::dbDisconnect(con)
   
   
}
