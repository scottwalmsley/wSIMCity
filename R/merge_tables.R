
#' Merge two sqlite DB tables together for alignment
#'
#' @param X the first data frame
#' @param Y the second data frame
#' @param ppm the mz tolerance for searching for a mass in the ref DB
#' @param rt_tol the rt tolerance for searching for a mass peak in the ref DB in minutes
#' @return data frame of merged tables
#' @export
merge_tables = function(X,Y,ppm, rt_tol){



   mzX = X$mz_ms1

   rtX = X$rt_ms1

   # get list of real features
   mzY = Y$mz_ms1

   rtY = Y$rt_ms1

   query_mzr = lapply(mzX, function(x) getMassTolRange(x,ppm))

   query_rtr = lapply(rtX, function(x) c(x-rt_tol,x+rt_tol))

   k = 1;

   counter= 0

   idx = merged_table = NULL


   for(i in 1:length(mzX)){

      w = which(mzY >  query_mzr[[i]][1] & mzY < query_mzr[[i]][2] & rtY > query_rtr[[i]][1] & rtY < query_rtr[[i]][2])

      e = X[i,]


      if(length(w) > 0){

         idx = c(idx,i)

         q = Y[w,]
         if(nrow(q)>1){
            w = which.min(q$ppm)
            q = q[w,]
         }

         # TODO pick formula

         merged_table = rbind(merged_table,e)

         k=k+1

         # remove from Y to avoid re-referencing and inserting duplicates
         Y = Y[-w,]

         mzY = mzY[-w]

         rtY = rtY[-w]

         counter = counter + 1

      }else{

         if(is.null(merged_table)){

            merged_table = e

         }else{

            merged_table = rbind(merged_table,e)
         }


         k=k+1

      }

   }

   merged_table = rbind(merged_table,Y)

   merged_table = merged_table[order(merged_table$mz_ms1,merged_table$rt_ms1),]

   merged_table
}

