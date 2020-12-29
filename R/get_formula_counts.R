#' Retrieve the atomic counts and mass defect for a given formula / mz
#'
#' @param f the molecular formula
#' @param mz the exact mass
#' @param kendrik_formula character vector string of the molecular formulae
#'
#' @return data.frame of results
#' @export
get_atom_counts <- function(f, mz = NULL, kendrik_formula = NULL){

   arr = unlist(strsplit(f,split = ''))

   atom_count = c(0,0,0,0)


   C = 12.00000
   H = 1.007825
   N = 14.003074
   O = 15.994915


   n_atoms = function(atom,arr){
      cnt = NULL
      pos = grep(atom,arr)
      if(  length(pos) > 0  ){


         g = grep('[0-9]', arr)
         w = which(g == pos + 1)
         if(length(w) == 1){
            cnt = paste(cnt, arr[g[w]], sep = '')
            w2 = which(g == pos + 2)
            if(length(w2) == 1){
               cnt = paste(cnt, arr[g[w2]],sep = '')
               w3 =  which(g == pos + 3)
               if(length(w3) == 1){
                  cnt = paste(cnt, arr[g[w3]],sep='')


               }

            }

         }

      }

      if(is.null(cnt)){
         if(length(pos) > 0){
            cnt = 1
         }else{
            cnt = 0
         }
      }else{
         cnt = as.numeric(cnt)
      }

      cnt
   }


   nC = n_atoms('C',arr)

   nH = n_atoms('H', arr)

   nN = n_atoms('N', arr)

   nO = n_atoms('O', arr)

   nominal_mass = 12*nC + 1*nH + 14*nN + 16*nO

   exact_mass = C*nC + H*nH + N*nN + O*nO

   DBE = nC + 1 - (nH/2) + (nN) / 2

   h_rule = (nominal_mass + nH)%%4

   mass_defect  = exact_mass - nominal_mass


   df = list( 'formula' = f,
              'C' = nC,
              'H' = nH,
              'N' = nN,
              'O' = nO,
              'exact_mass' = exact_mass,
              'nominal_mass' = nominal_mass,
              'mass_defect' = mass_defect,

              'kendrik_mass' = NULL,
              'nominal_kendrik_mass' = NULL,
              'kendrik_mass_defect' = NULL,
              'DBE' = DBE,
              'h_rule' = h_rule
   )

   if(!is.null(mz)){
      df$exact_mass = mz - H
      df$mass_defect = round(df$exact_mass)

   }

   if(!is.null(kendrik_formula)){

      arr = unlist(strsplit(kendrik_formula,split = ''))
      atom_count = c(0,0,0,0)

      nC = n_atoms('C',arr)

      nH = n_atoms('H', arr)

      nN = n_atoms('N', arr)

      nO = n_atoms('O', arr)

      nominal_mass = 12*nC + 1*nH + 14*nN + 16*nO

      exact_mass = C*nC + H*nH + N*nN + O*nO

      df$kendrik_mass = df$exact_mass*(nominal_mass/exact_mass)
      df$nominal_kendrik_mass = round(df$kendrik_mass)
      df$kendrik_mass_defect = df$nominal_kendrik_mass - df$kendrik_mass

   }
   df
}


#' Compute the delta formula for two features, the precursor and the aglycone
#'
#' @param f1 the first formula
#' @param f2 the aglycone formula
#'
#' @return numeric delta value
#' @export
delta_formula <- function(f1,f2){

   if(is.na(f1) | is.na(f2)){
      return(NA)

   }

   f1 =  get_atom_counts(f1)

   f2 = get_atom_counts(f2)

   f1  = c(f1$C,f1$H,f1$N,f1$O)

   f2  = c(f2$C,f2$H,f2$N,f2$O)

   delta = as.numeric(f1) - as.numeric(f2)

   delta = delta - c(5,8,0,3) #dR loss

   sum(abs(delta))

}





