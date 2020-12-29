#' Compute molecular formulae for a list of masses in the database
#'
#' @param DBE the minimum required double bond equivalence

#' @export
#'
compute_formulae <- function(DBE = 5){

   con <- RSQLite::dbConnect(RSQLite::SQLite(),db_name)

   dat <- 	RSQLite::dbGetQuery(con, "SELECT * FROM peak_group_data")

   RSQLite::dbDisconnect(con)

   elements <- Rdisop::initializeElements( c("C","H","N","O"))

   f1 = f2 = DBE1 = DBE2 =  deltaF = ppm = ppm2 = array(dim = nrow(dat), NA)

   PROTON = 1.00728

   options(warn = 0)

   for(i in 1:nrow(dat)){

      e = dat[i,]

      f = NULL
      f = Rdisop::decomposeMass(e$mz_ms1 - PROTON, z=0, elements= elements, ppm = 10)

      val = NULL

      if(!is.null(f)){
         val = getValid(f, DBE = 5, precursor = T)

      }else{
         val = NULL
      }

      f = NULL
      f = Rdisop::decomposeMass(e$mz_ms2 - PROTON, z=0, elements= elements,  ppm = 10)


      val2 = NULL
      if(!is.null(f)){
         val2 = getValid(f, DBE = 4, precursor = F)
      }else{
         val2 = NULL
      }

      if(!is.null(val) & !is.null(val2)){
         if(  length(val$formula > 0) & length(val2$formula) > 0){

            mat =  matrix(nrow = length(val$formula), ncol = length(val2$formula))

            ppmerr = array(dim = length(val$exactmass))

            ppmerr2 = array(dim = length(val2$exactmass))

            rownames(mat) = val$formula

            colnames(mat) = val2$formula

            for(k in 1:length(val$formula)){
               for(j in 1:length(val2$formula)){
                  mat[k,j] = delta_formula(val$formula[k],val2$formula[j])
                  ppmerr[k] = val$exactmass[k]
                  ppmerr2[j] = val2$exactmass[j]
               }
            }

            ppmerr = (ppmerr - (e$mz_ms1 - PROTON)) / (e$mz_ms1  - PROTON) * 1e6

            ppmerr2 = (ppmerr2 - (e$mz_ms2 - PROTON)) / (e$mz_ms2  - PROTON) * 1e6

            w = NULL

            if(!is.null(mat)){

               w = which(mat == min(mat,na.rm = T), arr.ind = TRUE)

            }

            if(nrow(w)>0){

               f1[i] =  paste(val$formula[w[,1]], collapse = ',')
               DBE1[i] = paste(val$DBE[w[,1]], collapse = ',')


               f2[i] =  paste(val2$formula[w[,2]], collapse = ',')
               DBE2[i] = paste(val2$DBE[w[,2]], collapse = ',')


               ppm[i] = paste(round(ppmerr[w[,1]],2), collapse = ',')
               ppm2[i] = paste(round(ppmerr2[w[,2]],2), collapse = ',')


               deltaF[i] = paste(mat[w], collapse = ',')

            }

         }
      }
      if( i %% 1000 == 0){
         cat(paste(i,"\r"))
      }

   }





   best_f1 = best_f2 = best_DBE1 = best_DBE2 = best_ppm1 = best_ppm2 = best_df = array(NA, dim = length(ppm))

   for(i in seq_along(ppm)){
      if(!is.na(ppm[i])){
         s = abs(as.numeric(unlist(strsplit(ppm[i],split = ',')))) + abs(as.numeric(unlist(strsplit(ppm2[i],split = ','))))
         w  = 	which.min(s)
         best_f1[i] = unlist(strsplit(f1[i], split = ','))[w]
         best_f2[i] =  unlist(strsplit(f2[i], split = ','))[w]
         best_df[i] =  as.numeric(unlist(strsplit(deltaF[i], split = ','))[w])
         best_DBE1[i] = as.numeric(unlist(strsplit(DBE1[i], split = ','))[w])
         best_DBE2[i] = as.numeric(unlist(strsplit(DBE2[i], split = ','))[w])
         best_ppm1[i] = as.numeric(unlist(strsplit(ppm[i], split = ','))[w])
         best_ppm2[i] = as.numeric(unlist(strsplit(ppm2[i], split = ','))[w])

      }

   }


   w = unique(c(which(is.na(f1)),which(is.na(f2))))

   df = data.frame(pk_idx = dat$pk_group[-w], f1 = f1[-w],f2 = f2[-w],df = deltaF[-w],ppm1 = ppm[-w], ppm2 = ppm2[-w] ,D1 = DBE1[-w],D2 = DBE2[-w])

   con <- RSQLite::dbConnect(RSQLite::SQLite(),db_name)
   RSQLite::dbWriteTable(con,"all_peak_formulae",df, overwrite = TRUE)
   RSQLite::dbDisconnect(con)


   df = data.frame(pk_idx = dat$pk_group[-w], f1 = best_f1[-w],f2 = best_f2[-w],df = best_df[-w],ppm1 = best_ppm1[-w], ppm2 = best_ppm2[-w] ,D1 = best_DBE1[-w],D2 = best_DBE2[-w])

   con <- RSQLite::dbConnect(RSQLite::SQLite(),db_name)
   RSQLite::dbWriteTable(con,"peak_formulae",df, overwrite = TRUE)
   RSQLite::dbDisconnect(con)


}


#' get valid atomic formulae
#'
#' @param flist list of molecular formulae
#' @param DBE minimum double bond equivalent
#' @param precursor is a precursor
#'
#' @return list containing: forula, DBE, exact_mass
#' @export
getValid <- function(flist, DBE =4, precursor = T){

   l=NULL
   w = which(flist$valid == 'Valid' & flist$DBE > DBE)

   if(length(w) > 0){


      l = list(
         formula = flist$formula[w],
         DBE = flist$DBE[w],
         exactmass = flist$exactmass[w]

      )
   }else{

      l = list(
         formula = NULL,
         DBE = NULL,
         exactmass = NULL
      )
      return(l)
   }


   g = grep('^C',l$formula)
   if(length(g) > 0){

      l = list(
         formula = l$formula[g],
         DBE = l$DBE[g],
         exactmass = l$exactmass[g]

      )
   }else{

      l = list(
         formula = NULL,
         DBE = NULL,
         exactmass = NULL
      )
      return(l)

   }



   atom_counts = NULL;

   for(i in 1:length(l$formula)){

      atom_counts = rbind(atom_counts,get_atom_counts(l$formula[i]))


   }

   w = NULL
   if(precursor){
      w = which(atom_counts[,2] > 8 & atom_counts[,3] > 12 & atom_counts[,4] > 1 & atom_counts[,5] > 2 & atom_counts[,12] > 4 )


   }else{
      w = which(atom_counts[,2] > 3 & atom_counts[,3] > 4 & atom_counts[,4] > 1  & atom_counts[,12] >  3)
   }

   if(length(w) > 0){


      l$formula = l$formula[w]
      l$DBE = l$DBE[w]
      l$exactmass = l$exactmass[w]

   }else{
      l = list(
         formula = NULL,
         DBE = NULL,
         exactmass = NULL
      )
      return(l)
   }


   return(l)

}


#' Numeric values for nucleotides
#'
#' @return
#'
formulaNucleotide <- function(){


   cytosine = c('C' = 9, 'H' = 13, 'N' = 3, 'O' = 4, 'DBE' = 5 )

   thymine  = c('C' = 10, 'H' = 14, 'N' = 2, 'O' = 5, 'DBE' = 5)

   guanine  = c('C' = 10, 'H' = 13, 'N' = 5, 'O' = 4, 'DBE' = 7)

   adenine  = c('C' = 10, 'H' = 13, 'N' = 5, 'O' = 3, 'DBE' = 7)


   a_cytosine = c('C' = 4, 'H' = 5, 'N' = 3, 'O' = 1, 'DBE' = 4 )

   a_thymine  = c('C' = 5, 'H' = 6, 'N' = 2, 'O' = 2,'DBE' = 4 )

   a_guanine  = c('C' = 5, 'H' = 5, 'N' = 5, 'O' = 1,'DBE' = 6 )

   a_adenine  = c('C' = 5, 'H' = 5, 'N' = 5, 'O' = 0,'DBE' = 6 )

}






