#' Create a hit table of results in the SQLite data base.
#'
#' @export
#'
create_hit_table <- function(){


   QRY = 'SELECT pk_group as pk_idx,
filter,
rt_ms1,
n_pk,
mz_ms1,
mz_ms2,
ppm,
ratio,
rho,
fwhm1,
fwhm2,
score_total,
int_ms1,
int_ms2,
ion_type_ms1,
ion_type_ms2,
peak_formulae.f1,
peak_formulae.f2,
peak_formulae.df,
peak_formulae.ppm1,
peak_formulae.ppm2
FROM peak_group_data
LEFT JOIN peak_formulae ON
peak_group_data.pk_group = peak_formulae.pk_idx
WHERE rho > 0.7
AND score_total > 0.7
AND ABS(ratio) < 2.38
AND df like \'0%\'
ORDER BY mz_ms1'


   con <- RSQLite::dbConnect(RSQLite::SQLite(),db_name)

   dat <- 	RSQLite::dbGetQuery(con, QRY)
   RSQLite::dbWriteTable(con,"hit_table",dat, overwrite = TRUE)
   RSQLite::dbDisconnect(con)

   con <- RSQLite::dbConnect(RSQLite::SQLite(),db_name)
   t = RSQLite::dbListTables(con)
   w = which(t == 'identified_known_adducts')

   if(length(w)>0){
      RSQLite::dbDisconnect(con)
      QRY = 'SELECT pk_idx as pk_group,
hit_table.filter,
hit_table.rt_ms1,
hit_table.n_pk,
hit_table.mz_ms1,
hit_table.mz_ms2,
hit_table.ppm,
hit_table.ratio,
hit_table.rho,
hit_table.fwhm1,
hit_table.fwhm2,
hit_table.score_total,
hit_table.int_ms1,
hit_table.int_ms2,
hit_table.ion_type_ms1,
hit_table.ion_type_ms2,
hit_table.f1,
hit_table.f2,
hit_table.df,
hit_table.ppm1,
hit_table.ppm2,
	identified_known_adducts.ID
FROM hit_table
LEFT JOIN identified_known_adducts ON
hit_table.pk_idx = identified_known_adducts.pk_group'
      #	identified_known_adducts.ID
      con <- RSQLite::dbConnect(RSQLite::SQLite(),db_name)

      dat <- 	RSQLite::dbGetQuery(con, QRY)
      RSQLite::dbWriteTable(con,"hit_table",dat, overwrite = TRUE)
      RSQLite::dbDisconnect(con)


   }

}
