#' Find the elbow point in the curve of variance explained by each successive PC. This can be used to determine the number of PCs to retain.
#'
#' @param variance Numeric vector containing the variance explained by each PC.
#'   Should be monotonic decreasing.
#'
#' @details Find the elbow point in the curve of variance explained by each successive PC. This can be used to determine the number of PCs to retain.
#'
#' @return An integer scalar specifying the number of PCs at the elbow point.
#'
#' @author Aaron Lun
#'
#' @examples
#'   col <- 20
#'   row <- 1000
#'   mat <- matrix(rexp(col*row, rate = 1), ncol = col)
#'
#'   # Adding some structure to make it more interesting.
#'   mat[1:100,1:3] <- mat[1:100,1:3] + 5
#'   mat[1:100+100,3:6] <- mat[1:100+100,3:6] + 5
#'   mat[1:100+200,7:10] <- mat[1:100+200,7:10] + 5
#'   mat[1:100+300,11:15] <- mat[1:100+300,11:15] + 5
#'
#'   p <- pca(mat)
#'   chosen <- findElbowPoint(p$variance)
#'
#'   plot(p$variance)
#'   abline(v=chosen, col="red")
#' 
#' @export
findElbowPoint <- function(variance) {
  if (is.unsorted(-variance)) {
    stop("'variance' should be sorted in decreasing order")
  }

  # Finding distance from each point on the curve to the diagonal.
  dy <- -diff(range(variance))
  dx <- length(variance) - 1
  l2 <- sqrt(dx^2 + dy^2)
  dx <- dx/l2
  dy <- dy/l2

  dy0 <- variance - variance[1]
  dx0 <- seq_along(variance) - 1

  parallel.l2 <- sqrt((dx0 * dx)^2 + (dy0 * dy)^2)
  normal.x <- dx0 - dx * parallel.l2
  normal.y <- dy0 - dy * parallel.l2
  normal.l2 <- sqrt(normal.x^2 + normal.y^2)

  # Picking the maximum normal that lies below the line.
  # If the entire curve is above the line, we just pick the last point.
  below.line <- normal.x < 0 & normal.y < 0
  if (!any(below.line)) {
      length(variance)
  } else {
      which(below.line)[which.max(normal.l2[below.line])]
  }
}
