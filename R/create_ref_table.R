
#' Creates a reference  / master table of compounds for RT alingments
#'
#' @param table_name the name of the SQLlite DB to output.   Do not add .sqlite to the name
#' @param db.list the list of SQLlite DBs to combine
#' @param ppm the ppm tolerance for matching the masses
#' @param rt_tol the rt toleracne in minumtes for searching for a mass peak in the ref DB
#'
#' @export
create_ref_table = function(table_name,db.list, ppm, rt_tol){

   assign("ref_table", paste(tolower(table_name),".sqlite",sep=""), .GlobalEnv)

   dat = NULL

   # get all the tables
   i = 1

   for(db in db.list){

      con <- RSQLite::dbConnect(RSQLite::SQLite(),db )

      dat[[i]] <- RSQLite::dbGetQuery(con,  'SELECT * from certified_hits') #

      RSQLite::dbDisconnect(con)

      i=i+1

   }



   # start with the largest table and work forward.
   print(nrow(dat[[1]]))

   mt = merge_tables(dat[[1]],dat[[2]], ppm = ppm, rt_tol = rt_tol)
   mt = refine_table(mt, ppm = ppm, rt_tol = rt_tol)
   print(nrow(mt))

   for(i in 3:length(dat)){

      mt = merge_tables(mt,dat[[i]], ppm = ppm,rt_tol = rt_tol)
      mt = refine_table(mt, ppm = ppm, rt_tol = rt_tol)
      print(nrow(mt))

   }


   con <- RSQLite::dbConnect(RSQLite::SQLite(),ref_table )

   RSQLite::dbWriteTable(con,'reference_table',mt,overwrite=T)

   RSQLite::dbDisconnect(con)


}
