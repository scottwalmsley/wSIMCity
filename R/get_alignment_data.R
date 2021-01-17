#' Get alignment data from aligned DNA adducts
#'
#' @param db_name name of the SQL table
#'
#' @return list of all tables (data frames) in the SQL database
#' @export
#'
getAlignmentData <- function(db_name){
  
  con <- RSQLite::dbConnect(RSQLite::SQLite(),paste(db_name,'.sqlite',sep = '') )
  
  tbs = RSQLite::dbListTables(con)
  
  db_data <- list()
  
  for(tb in tbs){
    
    db_data[[paste(tb)]] <- RSQLite::dbGetQuery(con,paste('SELECT * FROM',tb)) 
    
  }
  
  RSQLite::dbDisconnect(con)
  
  db_data
  
}
