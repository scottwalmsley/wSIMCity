compute_alignment_stats <- function(tb_name){
  
  con <- RSQLite::dbConnect(RSQLite::SQLite(),ref_db )
  ref_table = RSQLite::dbGetQuery(con,   'SELECT * from reference_table')
  RSQLite::dbDisconnect(con)
  
  RSQLite::dbWriteTable(con,'ms2_intensities', as.data.frame(cbind(IDX,ms2_mat)),overwrite=T)
  
}