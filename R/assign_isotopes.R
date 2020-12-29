#' Assign isotopes to a feature set
#'
#' @export
#'
assign_isotopes = function(){



   con <- RSQLite::dbConnect(RSQLite::SQLite(),db_name)
   dat <- 	RSQLite::dbGetQuery(con, "SELECT * FROM peak_group_data")
   RSQLite::dbDisconnect(con)

   filter_list = unique(dat$filter)


   # Make sure data is ordered

   out = data.frame()

   #
   for(filt in filter_list){

      # Get the subset of data matching the filter
      sub = dat[which(dat$filter == filt),]
      sub = sub[order(sub$scan_ms1,sub$mz_ms1),]
      out = rbind(out,sub)

   }

   dat = out
   rm(out)





   ion_type = has_isotope = is_isotope = isotope_ratio = pk_idx_isotope = array(dim = nrow(dat),NA)
   ion_type.ms2 = has_isotope.ms2 = is_isotope.ms2 = isotope_ratio.ms2 = pk_idx_isotope.ms2 = array(dim = nrow(dat),NA)

   k=1
   stop = FALSE

   for(filt in filter_list){

      #filt = 'FTMS + p NSI SIM ms [298.0000-332.0000]'
      # Get the subset of data matching the filter
      sub = dat[which(dat$filter == filt),]

      #Iterate over each unique scan and their mz peaks
      #for(scan in unique(sub$scan_ms1)){

         #pspec = sub[which(sub$scan_ms1==scan),]
         #iterate over each mz value in the sub
        #if(nrow(pspec)>1)
           #break
         start_idx = k
         #for(i in 1:nrow(pspec)){
         for(i in 1:nrow(sub)){
            # get the mass range for a mass's isotope
            mr = wSIMCity::getMassTolRange( (sub$mz_ms1[i]+1.003355),15)

            # check to see if a mass exists in that range
            w = which(sub$mz_ms1 > mr[1] & sub$mz_ms1 < mr[2])

            # Precursor peak
            if(length(w)>0){

               if(length(w)>1){
                  wm = which.min(abs(sub$ppm[w]))
                  w = w[wm]

               }
               ratio  = sub$int_ms1[w] / sub$int_ms1[i]
               #ratio  =  pspec$nl_int[w] / pspec$nl_int[i]
               if(stop){break}


               if(ratio < 1 ){

                  isotope_ratio[k] = ratio

                  has_isotope[k] = 1

                  pk_idx_isotope[k] = dat$pk_group[start_idx+w-1]

                  is_isotope[start_idx+w-1] = 1
               }

            }
            if(stop){break}

            mr = wSIMCity::getMassTolRange((sub$mz_ms2[i]+1.003355), 15  )

            ########################################
            ## product peak -- MS2 data
            w = which(sub$mz_ms2 > mr[1] & sub$mn_mz_ms2 < mr[2])

            if(length(w)>0){
               if(length(w)>1){
                  wm = which.min(abs(sub$ppm[w]))
                  w = w[wm]

               }
               ratio  =  sub$int_ms2[w] / sub$int_ms2[i]
               if(ratio < 1 ){

                  isotope_ratio.ms2[k] = ratio

                  has_isotope.ms2[k] = 1

                  pk_idx_isotope.ms2[k] = dat$pk_group[start_idx+w-1]

                  is_isotope.ms2[start_idx+w-1] = 1
               }

            }



            if(stop){break}

            k=k+1
         }
         if(stop){break}
      #}

      if(stop){break}
   }

   ion_type = array(dim=nrow(dat),NA)
   ion_type[which(is_isotope == 1)] = "[M1+H]"
   ion_type[which(has_isotope == 1)] = "[M0+H]"

   ion_type.ms2 = array(dim=nrow(dat),NA)
   ion_type.ms2[which(is_isotope.ms2 == 1)] = "[M1+H]"
   ion_type.ms2[which(has_isotope.ms2 == 1)] = "[M0+H]"

   dat = cbind(dat,
               ion_type_ms1 = ion_type,
               ion_type_ms2 = ion_type.ms2,
               isotope_ratio,
               isotope_ratio_ms2 = isotope_ratio.ms2,
               m1_idx = pk_idx_isotope,
               b1_idx = pk_idx_isotope.ms2)
   ########
   ## Refinement
   i=1
   out = NULL
   cat(paste("\n"))
   for(grp in unique(dat$pk_group)){

      w = which(dat$pk_group == grp)
      sub = dat[w,]

      g = grep('M',sub$ion_type_ms1)
      g2 = grep('M',sub$ion_type_ms2)



      ##################################
      if(length(g) > 0 & length(g2)>0){
         if(length(g) > 1 & length(g2)==1){
            sub$ion_type_ms1 = sub$ion_type_ms2 = unique(sub$ion_type_ms2[which(!is.na(sub$ion_type_ms2))])

         }
         if(length(g2) > 1 & length(g)==1){
            sub$ion_type_ms2 = sub$ion_type_ms1 = unique(sub$ion_type_ms1[which(!is.na(sub$ion_type_ms1))])
         }

         if(length(g) == length(g2)){
            tb = table(sub$ion_type_ms1)
            tb2 = table(sub$ion_type_ms2)
            if(dim(tb) == 1 & dim(tb2) ==1){
               sub$ion_type_ms1 = sub$ion_type_ms2 =  unique(sub$ion_type_ms2[which(!is.na(sub$ion_type_ms2))])
            }
            if(dim(tb) == 1 & dim(tb2)>1){
               sub$ion_type_ms1 = sub$ion_type_ms2 =  names(tb2[which.max(tb2)])
            }
            if(dim(tb) > 1 & dim(tb2)==1){
               sub$ion_type_ms1 = sub$ion_type_ms2 =  names(tb2[which.max(tb)])
            }
         }
      }

      if(length(g) > 0 & length(g2) ==0){ # case

         if(length(g) == 1){
            sub$ion_type_ms2 = sub$ion_type_ms1 = sub$ion_type_ms1[which(!is.na(sub$ion_type_ms1))]
         }

         if(length(g) > 1){
            tb = table(sub$ion_type_ms1)
            if(dim(tb) == 1){
               sub$ion_type_ms2 = sub$ion_type_ms1 = unique(sub$ion_type_ms1[which(!is.na(sub$ion_type_ms1))])
            }
            #
         }
      }

      if(i %% 500 == 0){
         cat(paste(i,"\r"))
      }
      i=i+1
      out = rbind(out,sub)
   }
   cat(paste("\n"))


   con <- RSQLite::dbConnect(RSQLite::SQLite(),db_name)
   RSQLite::dbWriteTable(con,"peak_group_data",out, overwrite=TRUE)
   RSQLite::dbDisconnect(con)


}
