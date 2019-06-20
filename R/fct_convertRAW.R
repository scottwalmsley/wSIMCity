#' Get the list of raw files
#' 
#' Returns a list of .raw files.
#'
#' @param raw_file_dir 
#'
#' @return
#' @export
#'
getRawFileList <- function(raw_file_dir){
  
  list.files(path = raw_file_dir,pattern = ".raw",full.names = T)
}



#' Get the sample names
#' 
#' Gets the sample names indicated using the raw file names.
#'
#' @param raw_file_dir character vector indicating the raw file directory to search.
#'
#' @return vector of sample names
#' @export
#'
getSampleNames <- function(raw_file_dir){
  
  
  sample_names <- sapply(getRawFileList(raw_file_dir),function(x) gsub(pattern = ".raw$","",x))
  sample_names <- sapply(sample_names, function(x) sub(pattern = raw_file_dir, "",x))
  sample_names <- sapply(sample_names, function(x) sub(pattern = "/","",x))
  names(sample_names) = NULL
  return(sample_names)

  
}



#' Make directories for each sample
#'
#' Creates directories for sample of interest indicateg in the global options parameter file.
#'
#'
#' @export
#'
#' @examples
#'
#' #makeSampleDir()
makeSampleDir <-function(result_dir,sample_names, scandef_file){
  
  
  
  l <- lapply(sample_names,function(x)
    dir.create(path = paste(result_dir,"/",x,sep=""))
  )
  
  
  l <- lapply(sample_names, function(x)
    makeSIMDir(scandef_file,paste(result_dir,"/",x,sep=""))
  )
  
  
  l <- lapply(sample_names,function(x)
    dir.create(path = paste(result_dir,"/",x,"/plots",sep=""))
  )
  
  l <- lapply(sample_names,function(x)
    dir.create(path = paste(result_dir,"/",x,"/results",sep=""))
  )
  
}



#' Create SIM directories
#'
#' makeSIMDir creates sim directories for performing independent  MSDIAL results files.
#'
#' @param scandef_file a non empty character vector giving the wide SIM scan definition file.
#' @param out_dir a non empty character vector giving the location used to write results.
#'
#'
#' @export
#'
#' @examples
#' #makeSIMDir(out_path = "results")
makeSIMDir <- function(scandef_file, out_path){

  scandef <- read.delim(scandef_file)
  
  for(i in 1:dim(scandef)[1]){

    scandef_line = paste(scandef[i,1],"_",scandef[i,4],"_",scandef[i,5],sep="")
    
    dir_path = paste(out_path,"/",scandef_line,sep="")

    if(!dir.exists(dir_path)){
   
         dir.create(dir_path)
    
    }
  }
}


#' Get sample directories
#'
#' @param results_path character vector indicating the directory in which to search.
#'
#' @return
#' @export
#'
getSampleDirectories <- function(results_path){
  list.dirs(path = results_path, full.names = TRUE, recursive = FALSE)
}




#' Convert raw data
#'
#' convertRaw converts vendor raw data using msconvert into an mzML file. This mzML file is then used by segmentmzML data to segment the file into SIM m/z ranges.
#'
#' @export
#'
#' @examples
#' #convertRaw()
#'
#'
convertRaw <- function(raw_file_dir, msconvert_path ){
  
  
  raw_file_list <- list.files(pattern = ".raw", path = raw_file_dir, full.names = TRUE)
 
  raw_file_dir
  sapply(raw_file_list, function(x) {

    cat(paste("Converting file", x ),"\n")
    
      system2(msconvert_path, args = c('--simAsSpectra', '--mzML', paste('-o', shQuote(raw_file_dir)), as.character(x)),invisible = FALSE)
  })

}


#' Segment mzML data into SIM windows
#'
#' Segregates mzml file into separate mzml files defined by mslevel and sim window.
#'
#'
#' @export
#'
#' @examples
#' #segmentMzMLDataSample
#'
segmentMzMLDataSample <- function(raw_file_directory, sample_directories, scandef_file = scandef_file){
  
  mzml_files <- list.files(pattern = "mzML",path = raw_file_dir, full.names = TRUE)
  i=1
  for(x in mzml_files){
    segmentMzMLData(single_mzml_file =x, sample_directory_path = sample_directories[i], scandef_file = scandef_file)
    i = i+1

  }
}




#' segmentMzMlData
#'
#' @param mzFile The input mzml file
#' 
#' @param sample_directory_path The path to your outputdata.  The folder must be preexisting.
#'
#'
#' @export
#'
#' @examples
#' #segmentMZMLData("test.mzML", out_path = "results/")
segmentMzMLData <- function(single_mzml_file,sample_directory_path, scandef_file = scandef_file){

  cat(paste("Opening file:",single_mzml_file,"\n"))

  mzml_data <- mzR::openMSfile(single_mzml_file)

  spectra <- mzR::spectra(mzml_data)

  header <- mzR::header(mzml_data)
  
  scandef = read.delim(scandef_file)

  start_scan <- seq(1, dim(scandef)[1], by=1)

  if(parallel::detectCores() > 4){
    cl <- parallel::makeCluster(2)
    doParallel::registerDoParallel(cl)
  }
  else{
    doParallel::registerDoParallel(cl <- makeCluster(1))
  }

  
  require(doParallel, quietly = TRUE)
  
  require(foreach,quietly = TRUE)
  
  foreach(start = start_scan,.packages= "mzR", .export = c("scandef")) %dopar% {

    scanSeq <- seq(start,length(spectra), by=20)
    
    sub_spectra <- spectra[scanSeq]
    
    sub_header <- header[scanSeq,]

    scan_no <- seq(1,length(sub_spectra), by=1)
    
    sub_header$seqNum  <- sub_header$acquisitionNum <- scan_no
    
    if(sub_header$msLevel[1] == 2){
      
      sub_header$msLevel = 1
    
    }
    
    if(sub_header$msLevel[1] == 3){
    
      print("This is not WSIM data, leaving process.")
      #shell(paste("rmdir",subDir))
      return();
    
    }

    rownames(sub_header) = scan_no

    if(scandef[start,1] == "WSIM"){
      mzlow  <- scandef[start,4]
      mzhigh <- scandef[start,5]
    }
    if(scandef[start,1] == "NL"){
      mzlow  <- scandef[start,2]
      mzhigh <- scandef[start,3]
    }

    out_file <- sub(".mzML$","", single_mzml_file)
    
    out_file <- strsplit(split = "/", out_file)[[1]]
    out_file <- out_file[length(out_file)]

    scandef_line <- paste(scandef[start,1],"_",scandef[start,4],"_",scandef[start,5],sep="")

    out_file_path <- paste(sample_directory_path,"/",scandef_line,"/", sep = "")


    fileName <- out_file

    fileName <- paste(out_file_path,scandef_line,"_",
                     fileName,".mzML",sep="")
    print(fileName)
    print("Writing mzML data...")

    mzR::writeMSData(object = sub_spectra,file = fileName, header= sub_header)
    Sys.sleep(5)
  }
  parallel::stopCluster(cl)
  #close(mzml_data)
}



