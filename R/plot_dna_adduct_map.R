#' Plot a DNA adduct heat map
#'
#' @param knowns boolean T / F plot only identified known DNA adducts?
#'
#' @return ggplot plot
#' @export
plot_dna_adduct_map <- function(knowns = F){
   
   #requireNamespace('RColorBrewer')
   requireNamespace('ggplot2')
   requireNamespace('plotly')
   requireNamespace('htmlwidgets')
   requireNamespace('withr')
   
   
   con <- RSQLite::dbConnect(RSQLite::SQLite(),db_name)
   
   if(knowns){
      
      d <- 	RSQLite::dbGetQuery(con, "SELECT * FROM identified_known_adducts where n_pk > 2 AND ABS(ratio) < 2.3 AND score_total > 0.5")#AND ion_type_ms1 like '%M%' ")
      
   }else{
      
      d <- 	RSQLite::dbGetQuery(con, "SELECT * FROM hit_table where (n_pk > 3 AND ABS(ratio) < 2.3 AND score_total > 0.8 AND rho > 0.8) OR ID like '%d%'" )# AND ion_type_ms1 like '%0%'")
   }
   
   RSQLite::dbDisconnect(con)
   
   if(nrow(d) == 0){
      
      return(NULL)
   }
   
   d$ID[which(is.na(d$ID))] = ''
   if(knowns){
      dat = data.frame('idx' = d$pk_group,
                       'ppm' = d$ppm,
                       'rt'=d$rt_ms1,
                       'mz'=d$mz_ms1,
                       'int'=d$int_ms2,
                       'id'=d$ID,
                       'score' = d$score_total,
                       'corr' = d$rho,
                       'info' = paste(d$ID,
                                      '\nms1 intensity: ',round(d$int_ms1,0),
                                      '\nms2 intensity: ',round(d$int_ms2,0),
                                      '\nlog ratio: ',round(d$ratio,2),'\nScore:',d$score_total,' ',' R:',d$rho,sep = ''),
                       stringsAsFactors = F
      )
   }else{
      dat = data.frame('idx' = d$pk_group,
                       'ppm' = d$ppm,
                       'rt'=d$rt_ms1,
                       'mz'=d$mz_ms1,
                       'int'=d$int_ms2,
                       'id'=d$ID,
                       'score' = d$score_total,
                       'corr' = d$rho,
                       'info' = paste(d$ID,
                                      '\nms1 intensity: ',round(d$int_ms1,0),
                                      '\nms2 intensity: ',round(d$int_ms2,0),
                                      '\nlog ratio: ',round(d$ratio,2),'\nScore:',d$score_total,' ',' R:',d$rho,sep = ''),
                       stringsAsFactors = F
      )
   }
   
   
   fig <- plotly::plot_ly(
      
      dat, x = ~rt, y = ~mz, 
      type = 'scatter', mode = 'markers',#size = ~log(int), 
      colors =  colorRampPalette(c('blue','red'))(100),
      color = ~score,
      sizes = c(1,200),
      
      text = ~info,
      hovertemplate = paste(
         '<i>m/z</i>: %{y:.4f}',
         '  rt(min): %{x:.3f}<br>',
         '%{text}',
         '<extra></extra> '),
      marker = list(opacity = ~0.2, sizemode = 'diameter',  size = ~log(int))#,sizeref = sizeref)
   )
   print(fig)
   
   
   
   if(knowns){
      htmlwidgets::saveWidget(fig, paste(sample_dir,'_knowns_map.html',sep = ''))
      
   }else{
      htmlwidgets::saveWidget(fig, paste(sample_dir,'_map.html',sep = ''))
   }
   
   
   
}





