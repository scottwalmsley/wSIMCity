#' Rine table and remove duplicates
#'
#' @param merged_table the original master table
#' @param ppm the mz tolerance window in ppm
#' @param rt_tol the retention time window in minutes
#'
#' @return data frame of the refined merged table
#' @export
refine_table <- function(merged_table, ppm = 10 , rt_tol = 1.5){

   #merged_table    =  o.merged_table

   #o.merged_table  =  merged_table

   mz = merged_table$mz_ms1

   rt = merged_table$rt_ms1

   mzr = lapply(mz, function(x) getMassTolRange(x, ppm))

   rtr = lapply(rt, function(x) c(x-rt_tol, x+rt_tol))

   refined_table = NULL

   l = length(mz)

   i=1

   while(i < (l+1) ){

      w =  which(      mz >  mzr[[i]][1] &
                       mz <  mzr[[i]][2] &
                       rt >  rtr[[i]][1] &
                       rt <  rtr[[i]][2])

      if(length(w) == 1){

         refined_table = rbind(refined_table,merged_table[w,])

         wm = 1

      }

      if(length(w) > 1){

         wm = which.max(merged_table$score_total[w])

      }

      merged_table = merged_table[-w[wm],]

      mz = mz[-w[wm]]

      rt = rt[-w[wm]]

      mzr = mzr[-w[wm]]

      rtr = rtr[-w[wm]]

      l = length(mz)

   }

   refined_table

}
