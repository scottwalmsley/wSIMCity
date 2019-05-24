#' Read MSDIAL result file
#'
#' Reads a single .msdial data file and returns a data.frame
#'
#' @param fh string
#'
#' @return data.frame
#' @export
#'
#' @examples
#'
#'
#'
#' #data <- msdial_data(fh = "results.msdial")
#'
msdial_data <- function(fh = NULL){
  
  require(data.table)
  
  cat(paste("   ...processing subfile:",fh,"\n"))
 
  d <- read.delim(fh, header= T)

  mz <- d$Precursor.m.z
  rt <- d$RT.min.

  adduct <- d$AdductIon
  isotope <- d$Isotope

  height <- d$Height
  area <- d$Area

  rowID <- d$PeakID

  name <- d$MetaboliteName
  comment <- d$Comment

  dt <- data.table(peak_index = rowID,
               
               rt = rt,
               mz = mz,
                   
               intensity = height,
               area = area,
             
               adduct_type = adduct,
               isotope = isotope,
               
               comment = comment
             
            )

  dt
}


#' Combine MSDIAL results
#'
#' Combines msdial data from the various SIM windows.
#'
#' @param sample_result_path
#'
#' @return list containing two dataframes, one each for the MS1 and MS2 data.
#' @export
#'
#' @examples
#'
#' #msdial_combineDataFiles("results/sample1/")
#'
msdial_combineDataFiles = function(sample_directory){

  cat(paste("processing msdial results for file:",sample_directory,"\n"))
  lf = list.files(pattern = "msdial$", path = sample_directory, full.names = T, recursive = F)

  nl_files = lf[grep("NL_", lf)]
  wsim_files = lf[grep("WSIM_",lf)]

  wsim_data = lapply(wsim_files,function(x) msdial_data(fh = x) )
  nl_data = lapply(nl_files,function(x) msdial_data(fh = x) )

  wsim_list = wsim_data[[1]]
  nl_list = nl_data[[1]]

 
  for(i in 2:length(wsim_data)){
    
    wsim_list = rbind(wsim_list,wsim_data[[i]])
    nl_list = rbind(nl_list,nl_data[[i]])
  
  }
  index <- seq_len(nrow(wsim_list))
  wsim_list <- cbind("index" = index, wsim_list)

  index <- seq_len(nrow(nl_list))
  nl_list <- cbind("index" = index, nl_list)
  
  
  list("wsim" = wsim_list,"nl" = nl_list)

}



#' getMSDIAL_results
#'
#' getMSDIAL_results is a wrapper function for running msdial_combineDataFiles.
#'
#'
#' @param sample_results_path character vector containing the parent directory of results
#'
#' @return list of data.frames containing wsim and nl feature data produced using MSDIAL
#' @export
#'
#' @examples
#'
#' multipleResults <- getMSDIAL_results("results")
#'
#'
getMSDIAL_results <- function(sample_directories){
  
  lapply(sample_directories, function(x) msdial_combineDataFiles(x))
  
}






