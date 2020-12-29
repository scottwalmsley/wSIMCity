



plot_density <- function(){
   requireNamespace('plotly')
   
   
   con <- RSQLite::dbConnect(RSQLite::SQLite(),db_name)
   dat <- 	RSQLite::dbGetQuery(con, paste("SELECT * FROM assigned_peak_groups"))
   RSQLite::dbDisconnect(con)
   
   
   
   dens_ppm <- density(dat$ppm)
   dens <- density(dat$dM/dat$nlMZ * 1e6)
   
   df <- data.frame(
      x = unlist(lapply(dens, "[[", "x")),
      y = unlist(lapply(dens, "[[", "y")),
      cut = rep(names(dens), each = length(dens[[1]]$x))
   )
   
   
   df = data.frame(x = dens$x, y = dens$y)
   
   fig <- plot_ly(df, x = ~x, y = ~y) 
   fig <- fig %>% add_lines()
   
   fig
   
   
}


plot_mz_histogram <- function(){
   requireNamespace('plotly')
   
   
   con <- RSQLite::dbConnect(RSQLite::SQLite(),db_name)
   dat <- 	RSQLite::dbGetQuery(con, paste("SELECT mz,nlMZ FROM assigned_peak_groups"))
   RSQLite::dbDisconnect(con)
   
   fig <- plot_ly(alpha = 0.6)
   fig <- fig %>% add_histogram(x = ~dat$mz)
   fig <- fig %>% add_histogram(x = ~dat$nlMZ)
   fig <- fig %>% layout(barmode = "overlay")
   
   fig
   
}


#' Get plotly figure for table
#'
#' @return plotly fig
#' @export
#'
get_hit_table <- function(){
   
   con <- RSQLite::dbConnect(RSQLite::SQLite(),db_name)
   dat <- 	RSQLite::dbGetQuery(con, paste("SELECT * FROM identified_known_adducts"))
   RSQLite::dbDisconnect(con)
   
   
   fig <- plotly::plot_ly(
      type = 'table',
      header = list(
         values = c("<b>Adducts</b>", names(dat)),
         align = c('left', rep('center', ncol(dat))),
         line = list(width = 1, color = 'black'),
         fill = list(color = 'rgb(235, 100, 230)'),
         font = list(family = "Arial", size = 14, color = "white")
      ),
      cells = list(
         values = rbind(
            rownames(dat), 
            t(as.matrix(unname(dat))),
            height = 20,width = 40
         ),
         align = c('left', rep('center', ncol(dat))),
         line = list(color = "black", width = 1),
         fill = list(color = c('rgb(235, 193, 238)', 'rgba(228, 222, 249, 0.65)')),
         font = list(family = "Arial", size = 12, color = c("black"))
      ))

   fig
   
   
}