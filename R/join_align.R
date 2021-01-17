#' Perform retention time aligment using the joing align algorithm
#'
#' @param ref_db the SQLlite reference database
#' @param db.list the list of SQLlite DBs to align
#' @param ppm mz tolerance window in ppm
#' @param rt_tol rt tolerance window in minutes
#' @param N number of samples to align
#'
#' @export
join_align = function(ref_db, db.list,  ppm = 7, rt_tol = 0.5, N){
   
   
   con <- RSQLite::dbConnect(RSQLite::SQLite(),ref_db )
   
   ref_table = RSQLite::dbGetQuery(con,   'SELECT * from reference_table')
   
   RSQLite::dbDisconnect(con)
   
   
   dat = NULL
   
   ####
   i = 1
   
   for(db in db.list){
      
      con <- RSQLite::dbConnect(RSQLite::SQLite(),db )
      
      dat[[i]] <- RSQLite::dbGetQuery(con,  'SELECT * from hit_table')
      
      RSQLite::dbDisconnect(con)
      
      i=i+1
      
   }
   
   
   requireNamespace('tools')
   
   requireNamespace('RSQLite')
   
   ms1_mat = ms2_mat = rtmat = ms1_mzmat = ms2_mzmat = itmat =  fmat = npk = ppmmat = rat = 
      fwhm1 = fwhm2 =  ID = rho = score_mz = score_rt = score_total = matrix(nrow = nrow(ref_table), ncol = N, NA )
   
   IDX = seq(1, nrow(rtmat))
   # Get reference
   #rt_tol = rt_tol / 60
   
   query_mzr = lapply(ref_table$mz_ms1, function(x) getMassTolRange(x,ppm))
   
   query_rtr = lapply(ref_table$rt_ms1, function(x) c(x-rt_tol,x+rt_tol))
   
   for(i in 1:length(query_mzr)){
      
      e = ref_table[i,]
      
      k = 1  # column counter
      
      for(s in 1:N){
         
         #g = grep(nm,l)
         sub = dat[[s]]
         
         w = NULL
         
         w = which(sub$mz_ms1 > query_mzr[[i]][1] & sub$mz_ms1 < query_mzr[[i]][2]
                   
                   & sub$rt_ms1 > query_rtr[[i]][1] & sub$rt_ms1 < query_rtr[[i]][2] )
         
         if(length(w) > 0){
            
            if(length(w)>1){
               
               df = data.frame(abs(sub$ppm1[w]),abs(sub$ppm2[w]))
               
               df = apply(df,1,sum)
               
               ms1_mat[i,s] = sub$int_ms1[w[which.min(df)]]
               
               ms2_mat[i,s] = sub$int_ms2[w[which.min(df)]]
               
               ms1_mzmat[i,s] = sub$mz_ms1[w[which.min(df)]]
               
               ms2_mzmat[i,s] = sub$mz_ms2[w[which.min(df)]]
               
               rtmat[i,s] = sub$rt_ms1[w[which.min(df)]]
               
               itmat[i,s] = sub$ion_type_ms1[w[which.min(df)]]
               
               fmat[i,s] = sub$f1[w[which.min(df)]]
               
               ppmmat[i,s] = sub$ppm[w[which.min(df)]]
               
               npk[i,s] = sub$n_pk[w[which.min(df)]]
               
               fwhm1[i,s] = sub$fwhm1[w[which.min(df)]]
               
               fwhm2[i,s] = sub$fwhm2[w[which.min(df)]]
               
               ID[i,s] = sub$ID[w[which.min(df)]]
               
               rho[i,s] = sub$rho[w[which.min(df)]]
               
               score_mz[i,s] = sub$score_mz[w[which.min(df)]]
               
               score_rt[i,s] = sub$score_rt[w[which.min(df)]]
               
               score_total[i,s] = sub$score_total[w[which.min(df)]]
               
            }else{
               ms1_mat[i,s] = sub$int_ms1[w]
               
               ms2_mat[i,s] = sub$int_ms2[w]
               
               ms1_mzmat[i,s] = sub$mz_ms1[w]
               
               ms2_mzmat[i,s] = sub$mz_ms2[w]
               
               rtmat[i,s] = sub$rt_ms1[w]
               
               itmat[i,s] = sub$ion_type_ms1[w]
               
               fmat[i,s] = sub$f1[w]
               
               ppmmat[i,s] = sub$ppm[w]
               
               npk[i,s] = sub$n_pk[w]
               
               fwhm1[i,s] = sub$fwhm1[w]
               
               fwhm2[i,s] = sub$fwhm2[w]
               
               ID[i,s] = sub$ID[w]
               
               rho[i,s] = sub$rho[w]
               
               score_mz[i,s] = sub$score_mz[w]
               
               score_rt[i,s] = sub$score_rt[w]
               
               score_total[i,s] = sub$score_total[w]
               
            }
            
         }
         
         k=k+1
         
      }
      
      
   }
   
   con <- RSQLite::dbConnect(RSQLite::SQLite(),ref_db )
   
   RSQLite::dbWriteTable(con,'fwhm1', as.data.frame(cbind(IDX,fwhm1)),overwrite=T)
   
   RSQLite::dbWriteTable(con,'fwhm2', as.data.frame(cbind(IDX,fwhm2)),overwrite=T)
   
   RSQLite::dbWriteTable(con,'id', as.data.frame(cbind(IDX,ID)),overwrite=T)
   
   RSQLite::dbWriteTable(con,'rho', as.data.frame(cbind(IDX,rho)),overwrite=T)
   
   RSQLite::dbWriteTable(con,'score_mz', as.data.frame(cbind(IDX,score_mz)),overwrite=T)
   
   RSQLite::dbWriteTable(con,'score_rt', as.data.frame(cbind(IDX,score_rt)),overwrite=T)
   
   RSQLite::dbWriteTable(con,'score_total', as.data.frame(cbind(IDX,score_total)),overwrite=T)
   
   RSQLite::dbWriteTable(con,'ms1_intensities', as.data.frame(cbind(IDX,ms1_mat)),overwrite=T)
   
   RSQLite::dbWriteTable(con,'ms2_intensities', as.data.frame(cbind(IDX,ms2_mat)),overwrite=T)
   
   RSQLite::dbWriteTable(con,'retention_times', as.data.frame(cbind(IDX,rtmat)),overwrite=T)
   
   RSQLite::dbWriteTable(con,'ms1_mz',          as.data.frame(cbind(IDX,ms1_mzmat)),overwrite=T)
   
   RSQLite::dbWriteTable(con,'ms2_mz',          as.data.frame(cbind(IDX,ms2_mzmat)),overwrite=T)
   
   RSQLite::dbWriteTable(con,'ion_types',       as.data.frame(cbind(IDX,itmat)),overwrite=T)
   
   RSQLite::dbWriteTable(con,'ppm_neutral_loss',       as.data.frame(cbind(IDX,ppmmat)),overwrite=T)
   
   RSQLite::dbWriteTable(con,'num_peaks',       as.data.frame(cbind(IDX,npk)),overwrite=T)
   
   RSQLite::dbWriteTable(con,'ms1_formulae',       as.data.frame(cbind(IDX,fmat)),overwrite=T)
   
   RSQLite::dbDisconnect(con)
   
}
