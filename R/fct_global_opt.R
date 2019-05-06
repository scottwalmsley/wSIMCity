

#' setGlobal.path
#'
#' Sets the global environment variables called through readGlobalOpts.
#'
#' @param option character vector
#' @param val character vector read from the global paramters file.
#'
#'
#'
#'
#' @export
#'
#'
setGlobal.path <- function(option,val){

       switch(option,
            mzml = assign("mzml_path",val, envir = .GlobalEnv),
            msdial = assign("msdial_path",val, envir = .GlobalEnv),
            msdial_param_path = assign("msdial_param_path",val, envir = .GlobalEnv),
            msconvert = assign("msconvert_path", val,envir = .GlobalEnv),
            knowns = assign("knowns", val,envir = .GlobalEnv),
            knowns_db = assign("knowns_db", val,envir = .GlobalEnv),
            scandef = assign("scandef_file", val,envir = .GlobalEnv),
         
            results = assign("results.path", val,envir = .GlobalEnv)

        )
}


#' Read the global options file and set global values
#'
#' Reads the global options parameter file and loads the variables into the global environment.
#'
#' @param fh character vector indicating the path to the global options file provided by the user.
#'
#'
#' @export
#'
#' @examples
#'
#' readGlobalOpts("global_params.txt")
readGlobalOpts <- function(fh){

   checkDependencies()

   library(mzR, quietly = T)
   library(foreach, quietly = T)
   library(doParallel,quietly = T)

  assign("nCore",detectCores() ,envir = .GlobalEnv)
   if(nCore > 4){
     nCore = nCore - 2
   }





   opts <- read.delim(file = fh, sep = "\t", header = T,stringsAsFactors = F)
   param <- as.character(opts$Parameter)


   setGlobal.path("mzml", opts[which(param == "analysis_folder"),2])
   setGlobal.path("msdial", opts[which(param == "msdial"),2])
   setGlobal.path("msdial_param_path", opts[which(param == "msdial_params"),2])
   setGlobal.path("msconvert", opts[which(param == "msconvert"),2])
   setGlobal.path("knowns", opts[which(param == "knowns"),2])
   setGlobal.path("scandef", opts[which(param == "scandef"),2])
   setGlobal.path("analysis", opts[which(param == "analysis_folder"),2])
   setGlobal.path("results", opts[which(param == "results_path"),2])


   if(!dir.exists(results_path)){
     dir.create(results_path)
   }


   assign("raw_file_list", list.files(mzml_path, pattern = "*.raw",full.names = T),envir = .GlobalEnv)

   assign("mzml_original",sub(".raw",".mzML",raw_file_list),envir = .GlobalEnv)

   lf = list.files(pattern = ".raw$", path = mzml_path, full.names = F)
   lf = sub(".raw","",lf)
   assign("sample_names", lf,envir = .GlobalEnv)
   rm(lf)

   assign("scandef", read.delim(scandef_file,header=T),envir = .GlobalEnv)
   assign("knowns_db", read.delim(knowns, header=T),envir = .GlobalEnv)
   assign("nl_list", read.delim("nl_list.txt", header=T,comment="#"),envir = .GlobalEnv)

}




#' Check for installed dependencies
#'
#' Checks for installed dependencies and then prompts the user to install if needed.
#' Called by readGlobalOptions().
#'
#'
#' @export
#'
#'
checkDependencies <- function(){
  
  if (!requireNamespace("mzR", quietly = TRUE)) {
    stop("Package 'mzR' needed for this function to work. Please install it.", call. = FALSE)
  }
  if (!requireNamespace("foreach", quietly = TRUE)) {
    stop("Package 'foreach' needed for this function to work. Please install it.", call. = FALSE)
  }
  if (!requireNamespace("doParallel", quietly = TRUE)) {
    stop("Package 'doParallel' needed for this function to work. Please install it.", call. = FALSE)
  }
  if(!(R.version$major >= 3 & as.numeric(R.version$minor)>=5.3)){
    warning("wSIMCity was developed with R version 3.5.3")
  }



}


global_options <- setRefClass(Class = "global_options",
                              fields = list(
                                raw_file_path = "character",
                                msconvert_path = "character",
                                scandef_file = "character",
                                msdial_path = "character",
                                msdial_param_path = "character",
                                raw_file_dir = "character"
                              
                                
                                
                                
                                
                              ))


