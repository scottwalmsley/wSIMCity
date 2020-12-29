#' count_peak_groups
#'
#' @param mzmin the minimum mz value
#' @param mzmax the maximum mz value
#' @param mzwid the mz slice \(default = 0.1\)
#' @export
count_peak_groups <- function(mzmin, mzmax, mzwid){

   con <- RSQLite::dbConnect(RSQLite::SQLite(),db_name)

   dat <- 	RSQLite::dbGetQuery(con, "SELECT * FROM raw_peaklist")

   RSQLite::dbDisconnect(con)

   cat(paste(dim(dat)))

   cat(paste("\nAssigning detected peaks with similar mass to peak groups and clustering by min rt tolerance of ",rt_tol," seconds....\n"))

   slices = seq(mzmin,mzmax,by = mzwid)

   hasgroup = array(dim = length(slices))

   gp = keepidx = NULL # keepidx

   gpidx = 1


   # Count data in each mz slice
   for(i in 1:(length(slices)-1)){


      w = which(dat$mz > slices[i] & dat$mz <= slices[i+1] + (mzwid /2))  ## beacause of overlap, the same cluster can enter the db twice.

      if(length(w)>0){

         hasgroup[i] = length(unique(w))


         if(length(dat$rt_min[w])>1){

            dv = (dat$mz[w] - mean(dat$mz[w]))/mean(dat$mz[w]) * 1e6

             d = dist(dv, method = 'manhattan')

             hc = hclust(d,method = "centroid")

             ct = cutree(hc,h=7)#4.5# 7 seconds

            for(ui in unique(ct)){

               id = which(ct==ui)

               keepidx= c(keepidx,w[id])

               gp = c(gp,rep(gpidx,length(w[id]) ))

               gpidx = gpidx+1

            }

         }
      }
      if(i %% 1000 == 0){
         cat(paste(i,"\r"))
      }
   }



   ki = group.index = array(dim = length(unique(keepidx)));

   i=1

   for(uk in unique(keepidx)){

      w = which(keepidx == uk)

      if(length(w)> 1){

            w = w[1]

      }

      ki[i] = keepidx[w]

      group.index[i] = gp[w]

      i=i+1

      if(i %% 1000 ==0)
         cat(paste(i,"\r"))

   }

   dat = dat[ki,]

   dat = cbind(pk_group = group.index,dat)

   cat(paste("Counted ", length(unique(gp)), " peak groups\n"))


   # remove duplicate peak groups
   #cat(paste("\nRemoving duplicate peak groups...\n"))

   out = data.frame()

   i=1



   for(g in unique(dat$pk_group)){

      if(i %% 1000 ==0)
         cat(paste(i,"\r"))

      sub = dat[which(dat$pk_group == g),]

      if(length(sub$scan) > length(unique(sub$scan))){

         lss = 1:length(sub$scan)

         for(j in lss){

            s = 	sub$scan[j]

            w = which(sub$scan == s)

            if(length(w) > 1){

               tmp.sub = sub[w,]

               wp = which(abs(tmp.sub$ppm) == min(abs(tmp.sub$ppm)))

               sub = sub[-w[-wp],]

            }


         }

      }

      out = rbind(out, sub)

      i=i+1
   }

   

   con <- RSQLite::dbConnect(RSQLite::SQLite(),db_name)
   RSQLite::dbWriteTable(con,"assigned_peak_groups",out, overwrite=TRUE)
   RSQLite::dbDisconnect(con)
   cat(paste("\n"))

   # Merge peak groups with similar mass
   dat = out
   rm(out)
   upg = master.upg =  unique(dat$pk_group)
   mean.mz = array(dim = length(upg))

   # get mean mz for each group
   cat(paste('\n'))

   i=1

   for(gp in upg){

      w = which(dat$pk_group == gp)

      mean.mz[i] = mean(dat$mz[w])



      if(i %% 1000 == 0){
         cat(paste(i,'\r'))
      }

      i=i+1

   }
   cat(paste('\n'))

   combined_groups = list()

   i = j =  1

   while(i < length(upg)+1){

      gp = upg[i]

      mz = mean.mz[i]

      mr = getMassTolRange(mz,ppm = 7*2)

      w = which(mean.mz > mr[1] & mean.mz < mr[2])

      if(length(w) > 0 ){


         combined_groups[[j]] = upg[w]

         j=j+1

         upg = upg[-w]

         mean.mz = mean.mz[-w]

      }


   }


   # Now iterate the masster list and assign groups to single groups
   for(i in 1:length(combined_groups)){

      gp = combined_groups[[i]]


      if(length(gp)>1){

         w = which.min(gp)

         idx = unlist(lapply(gp,function(x) which(dat$pk_group == x)))

         dat$pk_group[idx] = gp[w]


      }


   }
   con <- RSQLite::dbConnect(RSQLite::SQLite(),db_name)
   RSQLite::dbWriteTable(con,"assigned_peak_groups",dat, overwrite=TRUE)
   RSQLite::dbDisconnect(con)
   cat(paste("\n"))

   #Now look for breaks in peaks
   upg = master.upg =  unique(dat$pk_group)

   mean.mz = array(dim = length(upg))

   # get mean mz for each group
   cat(paste('\n'))

   i=1

   for(gp in upg){

      w = which(dat$pk_group == gp)

      mean.mz[i] = mean(dat$mz[w])


      if(i %% 1000 == 0){
         cat(paste(i,'\r'))
      }

      i=i+1

   }
   cat(paste('\n'))

########################################################################################################################
   ######################
   out = list()

   dc = nrow(scandef)

   i = j =  1

   while(i < length(upg)+1){

      gp = upg[i]

      #if(gp == 701)
       # break
      e = dat[which(dat$pk_group == gp),]

      if(nrow(e) > 1){

         e = e[order(e$scan),]


         deltDC = c(diff(e$scan),0)

         idx = 1

         for(m in 1:length(deltDC)){

            d = deltDC[m]

            if(idx == nrow(e)){

               out[[j]] = e[idx:m,]
               j=j+1
            }

            if(d > dc * 5){

               out[[j]] = e[idx:m,]

               idx = m+1

               j = j+1
            }


         }

      }


      if(i %% 1000 == 0){
         cat(paste(i,'\r'))
      }

      i = i+1
   }#EOL

   

      
      
      
      
      

   # Now remove single peaks again
   f.out = list()
   j=1
   for(e in out){
      if(nrow(e) > 1){
         f.out[[j]] = e
         j=j+1
      }
   }


   # now reindex the peaks
   i=1
   for(e in f.out){

      f.out[[i]]$pk_group = i
      i=i+1
   }
   f.out = do.call(rbind.data.frame, f.out)
   
   ## clean up multiple peak assignments to same neutral loss peak 
   out = list()
   i=1
   for( g in unique(f.out$pk_group)){
       e = f.out[which(f.out$pk_group == g),]
       # iterate each peak look for fuplicate rt values
       dat = data.frame()
       urt = unique(e$rt_min)
       for(rt in urt){
          w = which(e$rt_min == rt)
          if(length(w) > 1){
             w  = w[which.min(abs(e$ppm[w]))]
          }
          dat = rbind(dat,e[w,])
       }
       out[[i]]  = dat
       i=i+1
   }

   f.out = do.call(rbind.data.frame, out)
   
   con <- RSQLite::dbConnect(RSQLite::SQLite(),db_name)

   RSQLite::dbWriteTable(con,"assigned_peak_groups",f.out, overwrite=TRUE)

   RSQLite::dbDisconnect(con)

   cat(paste("\n"))

}
