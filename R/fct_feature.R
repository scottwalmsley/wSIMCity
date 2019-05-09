
#' Find features using MSDIAL
#'
#' A function which calls the MSDIAL Console application.  The location of the console appplication must be provided in the
#' global parameters file.   There must also be a msdial parameters file included.
#'
#' @param sample.result.path the path to the directory where the results will be stored, read using the global parameters file.
#'
#'
#'
#' @export
#'
#' @examples
#'
#'  #findPeaksMSDIAL("results/sample1/")
#'
findPeaksMSDIAL = function(sample_directory,scandef_file, msdial_path, msdial_param_path){
  require(doParallel)
  
  scandef <- read.delim(scandef_file)
  
  ld = list.dirs(sample_directory, full.names = T,recursive = T)
  
  g = grep("NL_",ld)
  
  g = c(grep("WSIM_",ld),g)
  
  ld = ld[g]
  
  sub_files = paste(scandef[,1],scandef[,4], scandef[,5],sep = "_")
  
  print("Step away, have a coffee..... this is a resource hog!")
  
  #cluster = makeCluster(nCore / 2)
  if(detectCores() > 3){
    cluster = parallel::makeCluster(parallel::detectCores() - 2)
    doParallel::registerDoParallel(cluster)
  }
  else{
    foreach::registerDoSEQ()
  }
  
 
  foreach(sub_analysis_path  = ld,.export = c("msdial_path","msdial_param_path")) %dopar% {
    
    
    print("Step away, have a coffee..... this is a resource hog!")
    
    
    system2(command = msdial_path, args = c("lcmsdda", paste("-i ",sub_analysis_path,"/ ",sep=""),paste("-o",sub_analysis_path),paste("-m",msdial_param_path)))
    Sys.sleep(15)
    
    msdial_file <- list.files(pattern = ".msdial",path = sub_analysis_path)
    
   
    system2(command = "mv",args = c(paste(sub_analysis_path,"/",msdial_file,sep=""), paste(sample_directory,"/",msdial_file,sep=""))  )
    Sys.sleep(15)
  }
  
  stopCluster(cluster)

}





#' Find peaks using MSDIAL, multiple samples
#'
#' Wrapper function to run MSDIAL on multiple samples.
#'
#'
#' @export
#'
#' @examples
#'
#' findPeaksMSDIAL_sample()
#'
findPeaksMSDIALSample <- function(sample_directories,scandef_file, msdial_path, msdial_param_path){
  
  sapply(sample_directories,function(x) findPeaksMSDIAL(x, scandef_file, msdial_path, msdial_param_path))
  
}











