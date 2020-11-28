
#' Make partial dependence plot including confidence interval
#'
#' \code{plot_partial_dependence} makes a PDP plot including confidence interval from fixed-effect covariance
#'
#' @author James Thorson
#' @export
plot_partial_dependence <-
function( x,
      train,
      filename = NULL,
      n_samples = 0,
      width = 4,
      height = 4,
      probs = c(0.025,0.975),
      ... )
{
  # load packages
  library(pdp)

  # Get MLE prediction and append intervals
  Return = NULL
  Return$Partial = partial(fit, prediction_type=1, train=train, origdata=train, ...)

  if( n_samples == 0 ){
    # Make plot
    library(ggplot2)
    if(!is.null(filename)) ThorsonUtilities::save_fig( file=filename, width=width, height=height )
      autoplot(Return$Partial)
    if(!is.null(filename)) dev.off()
  }else{
    yhat_zs = NULL
    for( sI in 1:n_samples ){
      if( sI%%max(1,floor(n_samples/10)) == 0 ){
        message( "  Finished sample ", sI, " of ",n_samples )
      }
      # predict( fit, prediction_type=2, newdata=Data_zk, origdata=Data_zk, seed=sI )
      #Partial = partial(fit, prediction_type=2, train=Data_zk, origdata=Data_zk, pred.var=c("BT","season"), type="regression", seed=sI )
      Partial = partial(fit, prediction_type=2, train=train, origdata=train, seed=sI, ...)
      yhat_zs = cbind( yhat_zs, Partial$yhat )
    }

    # Get MLE prediction and append intervals
    yhat_zz = t(apply( yhat_zs, MARGIN=1, FUN=quantile, probs=probs ))
    Return$Partial = cbind( Return$Partial, yhat_zz )
    Return$yhat_zs = yhat_zs

    # Make plot
    library(FishStatsUtils)
    if(!is.null(filename)) ThorsonUtilities::save_fig( file=filename, width=width, height=height )
      plot_timeseries( x = Return$Partial[,1],
        y = Return$Partial[,'yhat'],
        ybounds = Return$Partial[,ncol(Return$Partial)-1:0],
        bounds_type = "shading",
        bounds_args = list(col = rgb(0,0,0,0.2)),
        fn  =  plot,
        type  =  "l",
        xlab = colnames(Return$Partial)[1],
        ylab = "Partial dependence")
    if(!is.null(filename)) dev.off()
  }

  # Invisible return
  return(invisible(Return))
}
