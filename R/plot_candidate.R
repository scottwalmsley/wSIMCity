#' Plot candidate DNA adduct
#'
#' @param ppm tolerance window in ppm
#'
#' @export
plot_candidate <- function(ppm = 7){

   con <- RSQLite::dbConnect(RSQLite::SQLite(),db_name)

   dat <- 	RSQLite::dbGetQuery(con, "SELECT * FROM hit_table where n_pk > 3 AND ABS(ratio) < 2.3 AND score_total > 0.8 AND rho > 0.8 ")
   
   id.dat <- 	RSQLite::dbGetQuery(con, "SELECT * FROM identified_known_adducts where n_pk > 3 AND score_total > 0.5")
   
   RSQLite::dbDisconnect(con)

   for(i in 1:nrow(dat)){

      con <- RSQLite::dbConnect(RSQLite::SQLite(),db_name)
      query = paste("SELECT * from assigned_peak_groups WHERE pk_group =",dat$pk_group[i])
      e <- 	RSQLite::dbGetQuery(con, query)
      RSQLite::dbDisconnect(con)

      e = e[order(e$rt_min),]

      mat = cbind(e$rt_min,e$p_int, e$nl_int)
      mat = mat[order(mat[,1]),]


      rtdelt = mean(diff(mat[,1]))

      mat = rbind(c(mat[1,1]-rtdelt,0.5*mat[1,2],0.5*mat[1,3]),mat)
      mat = rbind(c(mat[1,1]-rtdelt,0.5*mat[1,2],0.5*mat[1,3]),mat)
      mat = rbind(mat,c(mat[nrow(mat),1]+rtdelt,0.5*mat[nrow(mat),2],0.5*mat[nrow(mat),3]))
      mat = rbind(mat,c(mat[nrow(mat),1]+rtdelt,0.5*mat[nrow(mat),2],0.5*mat[nrow(mat),3]))

      smooth_data.p <- smooth.spline(mat[,1:2],spar = 0.3)
      smooth_data.b <- smooth.spline(mat[,c(1,3)],spar = 0.3)

      smat.p = cbind(smooth_data.p$x, smooth_data.p$y)
      smat.b = cbind(smooth_data.b$x, smooth_data.b$y)

      smat.p[which(smat.p[,2] < 0),2] = 0
      smat.b[which(smat.b[,2] < 0),2] = 0


      r = round(cor(smat.p[,2], smat.b[,2]),3)


      cat(paste(dat$f1[i],dat$f2[i]),"\n")



      fileHandle <- paste(db_name,round(dat$mz_ms1[i],4),"@",round(dat$rt_ms1[i],2),".png",sep="")
      
      #plot_peak_group(pk_group = dat$pk_group[i])



      png(fileHandle, width = 400, height = 600)

      layout(matrix(c(1,1,2,2,3,4), nrow = 3, ncol = 2, byrow =T))


      main = paste(dat$mz_ms1[i], dat$f1[i], dat$f2[i],'\ncor:',r,'\nscore:',dat$score_total[i])


      mx.y = max(c(max(smat.p[,2]), max(smat.b[,2]),max(mat[,2:3])))

      #######PNG
      plot(mat[,1:2], type = "l", ylim = c(0,mx.y*1.1),lty = 2, col=4, main =main, cex.main = 1.25,
           xlab = "RT", ylab = 'Intensity', bty = 'n',las= 2)
      lines(mat[,c(1,3)], type = "l", lty = 2, col=2)
      lines(smat.p, col=4, lwd =1, lty = 1)
      lines(smat.b, col=2, lwd =1, lty = 1)
      abline(h=0)

      plot(e$scan, e$p_int, lwd = 1, ylim = c(0,mx.y*1.1),type = 'h', xlab = 'Scan', ylab = 'Intensity', bty = 'n', col=4)
      lines(e$scan+5, e$nl_int, lwd = 1, type = 'h', col=2, las= 2)
      abline(h=0)


      mz = dat$mz_ms1[i]


      wmx = which.max((  lapply(spectra[e$scan], function(x) getPeak(x,mz,ppm)$y)))

      spectrum = spectra[[e$scan[wmx]]]
      plotPeak(spectrum, mz,ppm,col=4)


      mz = dat$mz_ms2[i]


      idx = e$scan+1
      wmx = which.max(lapply(spectra[idx], function(x) getPeak(x,mz,ppm)$y))


      spectrum = spectra[[ e$scan[wmx]+1]]
      plotPeak(spectrum, mz,ppm,col = 2)


      dev.off()

      ###########PDF

      pdf(width = 5,height = 8,file = sub('.png','.pdf',fileHandle), pointsize = 11)
      par(mfrow=c(4,1))
      layout(matrix(c(1,1,2,2,3,4), nrow = 3, ncol = 2, byrow =T))
      main = paste(dat$mz_ms1[i], dat$f1[i], dat$f2[i], '\n',r
      )

      plot(mat[,1:2], type = "l", ylim = c(0,mx.y*1.1),lty = 2, col=4, main =main, cex.main = 1.25,
           xlab = "RT", ylab = 'Intensity', bty = 'n', cex.axis = 0.8,cex.main = 0.8, las = 2)
      lines(mat[,c(1,3)], type = "l", lty = 2, col=2)
      lines(smat.p, col=4, lwd =1, lty = 1)
      lines(smat.b, col=2, lwd =1, lty = 1)
      abline(h=0)

      plot(e$scan, e$p_int, lwd = 1, ylim = c(0,mx.y*1.1), type = 'h', xlab = 'Scan', ylab = 'Intensity', bty = 'n',
           cex.axis = 0.8, las= 2, main = 'Raw intensities',cex.main = 0.8, col=4)
      lines(e$scan+5, e$nl_int, lwd = 1, type = 'h', col=2)
      abline(h=0)


      mz = dat$mz_ms1[i]
      s.spectra = spectra[e$scan]
      wmx = which.max(lapply(s.spectra, function(x) getPeak(x,mz,ppm)$y ))
      spectrum = spectra[[e$scan[wmx]]]
      plotPeak(spectrum, mz,ppm,main = 'Precursor', col = 4)


      mz = dat$mz_ms2[i]
      wmx = which.max(lapply(spectra[(e$scan+1)], function(x) getPeak(x,mz,ppm)$y))
      spectrum = spectra[[ e$scan[wmx]+1]]

      plotPeak(spectrum, mz,ppm, main = "Aglycone",col=2)


      dev.off()




      }

 





}




