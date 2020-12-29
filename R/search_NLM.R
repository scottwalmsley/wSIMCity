#' search_NLM Wrappper function to search for the neutral loss in wide SIM DIA data.
#'
#' @param file the mzml file to search. Must be wide SIM MS1 / MS2 DIA data with SIM segments the same mz range
#' @param type 'SIM' for true wide SIM DIA data or any other value for pseudo MS1 0 volt data \(true MS2\)
#' @param raw the Thermo .raw file
#' @param delta_search_mass decimal value indicating the Da loss or gain.  Must be signed if negative \(-\)
#' @param search_name character value indicating the value used
#' @param scandef_file the name of the scan definition file. Currently not used.
#' @param sim_width the mz width of the wSIM windows
#' @param overlapMZ the mz overlap (if any) between wSIM windows
#' @param ppm_tol the ppm mass tolerance in parts per million
#' @param rt_tol the retention time tolerance for the search in minutes
#' @param mzmin the minimum mz, default = 135
#' @param mzmax the maximum mz, default = 634
#' @param mzwid the mz slice size, default = 0.01 for Orbitrap
#'
#' @export
#'
search_NLM <-  function(file = NULL,type = 'SIM', raw = FALSE, delta_search_mass = NULL, search_name = "",scandef_file = NULL, sim_width = NULL,overlapMZ = NULL,ppm_tol = 20, rt_tol = 9,mzmin = 135,mzmax = 634,mzwid = 0.01){

  requireNamespace('mzR', quietly = TRUE)

  if(is.null(file) || is.null(delta_search_mass)){

    return(cat("\nBoth an input file and a mass difference to search for is required.\n"))

  }

  if(grepl("\\.mzml$", tolower(file)) && isTRUE(raw) ){

    str <- "You have indicated a raw file (\'raw = TRUE\') and a mzml file as input."
    str <- paste(str, "If using a raw file correct the file input name to a .raw.\nIf using a mzml file, set \'raw = FALSE\'.\n",sep = "\n")
    return(cat(str))

  }

  if(is.null(delta_search_mass)){

        str <- paste("\nA delta_search_mass is required for the search. Define delta_search_mass and include a negative sign for a mass loss (eg: \'-116.0474\') \n")

        return(cat(str))

  }

  data("mdnadb")
  
  assign("sample_name", sub(".mzml","", tolower(basename(file))),.GlobalEnv)

  assign('sample_dir', sample_name,envir = .GlobalEnv)

  assign('plot_dir', paste(sample_name,'/plots',sep = ''),envir = .GlobalEnv)
  
  cat(paste("\nOpening and processing mzML file:",file,"\n" ,sep = ""))

  assign("delta_search_mass",delta_search_mass ,envir = .GlobalEnv)

  assign("search_name",search_name ,envir = .GlobalEnv)

  assign("rt_tol", rt_tol,envir = .GlobalEnv)

  assign("ppm_tol", ppm_tol, envir = .GlobalEnv)

  assign("mzml_data",mzR::openMSfile(file),envir = .GlobalEnv)

  assign("spectra",mzR::spectra(mzml_data),envir = .GlobalEnv)

  assign("header",mzR::header(mzml_data),envir = .GlobalEnv)

  mzR::close(mzml_data)

  if(!dir.exists(sample_dir)){

    dir.create(sample_dir)
    dir.create(plot_dir)
    #dir.create()

  }else{

    dir_files = list.files(path = sample_dir)

    if(length(dir_files) > 0 ){

      do.call(what = file.remove, list(list.files(path = sample_dir,full.names = T)))

    }

  }

  assign("scandef", setScanDef(header, type = type),.GlobalEnv)
  #assign("scandef", setScanDef(header, type = type, sim_width = sim_width, overlapMZ = overlapMZ),.GlobalEnv)

  cat(paste("Read mzMLfile. Now starting delta mass search for loss(-) / gain of: ", delta_search_mass," Da +/- ", ppm_tol, "ppm\n",sep = ""))

  search_delta_mass(rt_tol=rt_tol, ppm_tol = ppm_tol)

  cat(paste("\nAsssigning peaks to peak groups in retention time...\n"))

  count_peak_groups(mzmin,mzmax, mzwid)

  compute_group_data()
  
  
  wSIMCity::groupFWHM(sigma = 6, perfwhm = 0.5)
  wSIMCity::findIsotopes(mslevel = 1)
  wSIMCity::findIsotopes(mslevel = 2)
  
  compute_formulae(DBE = 4)

  knwn  = find_known_adducts(search_tol = 10,minscore = 0.7)

  create_hit_table()

  plot_dna_adduct_map(knowns = F)
  
  plot_dna_adduct_map(knowns = T)

  plot_candidate(ppm = 7)
  
  createTables()
  

}





