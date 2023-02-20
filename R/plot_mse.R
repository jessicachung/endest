#' Plot MSE curve
#'
#' Plot the MSE curve of a sample using base graphics.
#'
#' @param x Cycle time estimation list generated from the \code{estimate_cycle_time} function.
#' @param sample_name Sample name of sample to plot corresponding to one of the column names in the original expression matrix.
#' @param ... Other optional arguments to be passed to \code{plot}.
#'
#' @return None
#' @export
#' @importFrom graphics plot abline mtext
#'
#' @examples
#' # Simulate gene expression data for 100 genes and 3 samples
#' entrez_ids <- 9801:9900
#' y <- matrix(rnorm(100*3, 5, 2), nrow=100)
#' colnames(y) <- paste0("sample_", 1:3)
#' results <- estimate_cycle_time(exprs=y, entrez_ids=entrez_ids)
#' plot_mse(results, sample_name="sample_1")
plot_mse <- function(x, sample_name, ...) {
  plot(y=x$mse[,sample_name], x=0:99, type="l",
       xlab="Cycle Time", ylab="Mean Squared Error", ...)
  abline(v=x$estimated_time[sample_name], col="#999999", lty=2)
  mtext(sample_name)
}
