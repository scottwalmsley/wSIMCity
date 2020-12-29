#' search_delta_mass
#'
#' @param type the search type, either 'SIM' or other for pequdo SIM \(eg. 0-25V consecutive MS2 cycles\).
#' @param rt_tol the retention time tolerance to search for the neutral loss mass aglycone peak in seconds.  Usually set to the time for 2 duty cycles.
#' @param ppm_tol the mass tolerance in ppm.  Use a wide window, typically 20-50 to aid detection of distribution of mass errors.
#'
#' @export
#'
search_delta_mass <- function(type = 'SIM',rt_tol = 9,ppm_tol = 20){

  cat(paste("\nEntering mass search for delta mass: ",delta_search_mass," Da\n", sep=""))

  peak.idx <- dat <- dRT <- RT <- dM <- dMppm <- dint <- cMZ <- pMZ <- tint <- pint <- scan <- filters <-  NULL

  assign('filterString',unique(header$filterString),.GlobalEnv)

  assign('g',grep("Full ms2", filterString),.GlobalEnv)

  if(length(filterString[-g]) != 0){
    cat(paste("There are ",length(filterString[-g]),"SIM or mass aquisition windows (MS2 - 0 meV) to be processed.\n"))
  }else{
    cat(paste("There are ",length(filterString),"mass aquisition windows to be processed.\n"))
  }


  if(length(filterString[-g]) == 0){
    ms1FilterString = filterString
  }else{
    ms1FilterString = filterString[-g]
  }

  bk = F
  filteridx <<-1
  scandef_idx <<- 1
  t = Sys.time()
  
  dat = lapply(ms1FilterString, function(x) search_filter(x,type = type))
  print(Sys.time() - t)

  dat <- do.call("rbind", dat)

  cat(paste("Finished searching mass list for delta mass.\n"))


  assign("db_name", paste(sample_dir,"/", tolower(sample_name),"-",format(Sys.time(),"%Y%m%d_%H%M"),".sqlite",sep=""), .GlobalEnv)

  con <- RSQLite::dbConnect(RSQLite::SQLite(),db_name)
  RSQLite::dbWriteTable(con,"raw_peaklist",dat)
  RSQLite::dbDisconnect(con)

  cat(paste("\nIntermediate search results located in SQLite DB: ",db_name),"\n")
  cat(paste("This DB will be used in subsequent steps if you keep this workspace loaded.\n"))

 # cat(paste(dim(dat)))



}
