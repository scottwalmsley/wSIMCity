#' set_scan_def
#'
#' @param header a mzR header object
#' @param type either 'SIM' to denote MS1 as SIM spectra or any other as MS1 data collected in the C Trap.
#' @param sim_width the emz width of the SIM acquisition window.
#' @param overlapMZ optional numeric value describing how much overlap exists between SIM windows.
#'
#' @return dataframe containing scan mz window ranges
#' @export
#'
#'
setScanDef <- function(header, type = 'SIM',sim_width = NULL , overlapMZ = NULL){

  filterString = unique(header$filterString)

  if(type=='SIM'){

    g = grep('SIM', filterString)

    ms1 = filterString[g]

    ms2 = filterString[-g]

  }else{

    ms1 = ms2 = filterString

  }


  ScanType = WindowStart = WindowEnd = AquisitionStart = AquisitionEnd = array(dim = length(filterString))

  if(type=='SIM'){

    ScanType[g] = 'WSIM'

    ScanType[-g] = 'NL'

    }else{

      ScanType = 'stepped'

    }

  srch = lapply( filterString, function(x) strsplit(x,' ')[[1]])
  srch = lapply( srch , function(x) x[length(x)])
  srch = lapply(srch, function(x) sub('\\[','', x))
  srch = lapply(srch, function(x) sub('\\]','', x))
  srch = lapply( srch, function(x) strsplit(x,'-')[[1]])


  WindowStart = unlist(lapply(srch, function(x) as.numeric(x[1]) ))
  WindowEnd = unlist(lapply(srch, function(x) as.numeric(x[2]) ))

  if(type=='SIM'){
    WindowEnd[-g] = WindowEnd[g]
  }


  AcquisitionStart = WindowStart
  if(type=='SIM'){
    AcquisitionStart[-g] = WindowStart[g]
  }

  AcquisitionEnd = WindowEnd
  if(type=='SIM'){
    AcquisitionEnd[-g] = WindowEnd[g]
  }



  scandef = data.frame (
    ScanType = ScanType,
    WindowStart  = WindowStart,
    WindowEnd = WindowEnd,
    AcquisitionStart = AcquisitionStart,
    AcquisitionEnd = AcquisitionEnd

  )

  if(!is.null(sim_width)){
    if(type == 'SIM'){
      scandef$AcquisitionStart[g] = scandef$AcquisitionEnd[g] - sim_width - overlapMZ
      scandef$AcquisitionStart[-g] = scandef$WindowStart[-g]
    }
  }

  nscandef = NULL
  if(type != 'SIM'){
    nl = scandef
    nl$ScanType = 'pNL'

    for( i in 1:nrow(nl)){
      nscandef = rbind(nscandef, rbind(scandef[i,],nl[i,]))
    }
    scandef = nscandef
  }

  rownames(scandef) = NULL

  scandef

}
