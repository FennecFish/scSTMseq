#' Plot effects of covariates on topic proportions
#'
#' Visualizes how a covariate influences topic proportions from an
#' \code{estimateEffect}-style object produced for scSTMseq models, with options
#' for point estimates, contrasts between levels, or continuous trends. The
#' method wraps internal helpers (\code{plotPointEstimate()}, \code{plotDifference()},
#' \code{plotContinuous()}) to produce base R plots and returns (invisibly) the
#' data used to draw them.
#'
#' @param x An \code{estimateEffect}-like object containing the fitted regression
#'   of topic prevalence on covariates (typically returned by an
#'   \code{estimateEffect()} function in your package). Must include at least
#'   \code{$topics}, \code{$parameters}, and (optionally) \code{$ref.vec}.
#' @param covariate Character scalar; the name of the covariate whose effect you
#'   want to plot.
#' @param model The fitted scSTM/STM model used to derive metadata when needed
#'   (e.g., to infer levels for \code{method = "difference"} when \code{ref}/\code{alt}
#'   are omitted). Can be \code{NULL} if not required by the chosen method.
#' @param topics Integer or character vector of topic indices/names to plot.
#'   Defaults to \code{x$topics}.
#' @param method One of \code{"pointestimate"}, \code{"difference"}, or
#'   \code{"continuous"}. See Details.
#' @param ref,alt Reference and alternative level for \code{method = "difference"}.
#'   If omitted, and if \code{covariate} appears in \code{x$ref.vec}, the function
#'   will use \code{x$ref.vec[[covariate]]} as \code{ref} and infer a single
#'   \code{alt} level from the model metadata. Otherwise both must be supplied.
#' @param moderator Optional character scalar giving the name of a moderator
#'   variable for conditional effects (e.g., an interaction).
#' @param moderator.value Value of \code{moderator} at which to evaluate the
#'   effect for plotting.
#' @param npoints Number of grid points for \code{method = "continuous"}.
#'   Default \code{100}.
#' @param nsims Number of posterior simulations of coefficients used to form
#'   uncertainty intervals. Default \code{500}.
#' @param ci.level Confidence level for intervals, as a probability in \eqn{(0,1]}.
#'   Default \code{0.95}.
#' @param xlim,ylim Numeric length-2 vectors giving x/y-axis limits, or
#'   \code{NULL} for automatic.
#' @param xlab,ylab Axis labels. If \code{ylab = NULL}, a sensible default is used.
#' @param main Plot title.
#' @param printlegend Logical; print a legend for topics when appropriate.
#'   Default \code{TRUE}.
#' @param labeltype How to label topics in legends/titles; typically
#'   \code{"numbers"} (default).
#' @param n Integer tuning parameter passed to labeling/selection helper(s)
#'   (e.g., number of top words per topic when constructing labels). Default \code{7}.
#' @param frexw Numeric in \eqn{[0,1]}; FREX weighting for topic labels when
#'   applicable. Default \code{0.5}.
#' @param add Logical; add to an existing plot instead of creating a new one.
#'   Default \code{FALSE}.
#' @param linecol Line color. Default \code{"black"}.
#' @param width Numeric; graphical width parameter used by some plot types
#'   (e.g., bar widths). Default \code{25}.
#' @param verbose.labels Logical; if \code{TRUE}, use more verbose topic labels
#'   where available. Default \code{TRUE}.
#' @param family Base graphics font family.
#' @param custom.labels Optional character vector to override automatic topic
#'   labels (must align with \code{topics}).
#' @param omit.plot Logical; if \code{TRUE}, compute and return plot data but do
#'   not draw the plot. Default \code{FALSE}.
#' @param ... Additional arguments passed to the method-specific plotting helper.
#'
#' @details
#' \describe{
#'   \item{\code{method = "pointestimate"}}{Plots point estimates (and confidence
#'   intervals) of topic prevalence as a function of the specified covariate.}
#'   \item{\code{method = "difference"}}{Plots the contrast in topic prevalence
#'   between two levels (\code{ref} vs \code{alt}) of a categorical covariate.
#'   If \code{ref}/\code{alt} are omitted and the covariate was registered in
#'   \code{x$ref.vec}, the function attempts to infer them; otherwise both must be supplied.}
#'   \item{\code{method = "continuous"}}{Plots a smooth relationship for a numeric
#'   covariate over a grid of \code{npoints}, with uncertainty from posterior
#'   simulations of coefficients.}
#' }
#'
#' The function relies on \code{simBetas()} to draw posterior coefficient samples,
#' and on \code{produce_cmatrix()} to build covariate design matrices for the
#' requested contrast or grid. It is compatible with objects whose underlying
#' model has class \code{"scSTM"} and/or \code{"STM"}.
#'
#' @return
#' Invisibly returns the object produced by the method-specific helper
#' (\code{plotPointEstimate()}, \code{plotDifference()}, or \code{plotContinuous()}),
#' typically a list containing the data used to draw the plot (e.g., grid,
#' estimates, and confidence intervals). A plot is drawn as a side effect unless
#' \code{omit.plot = TRUE} or \code{add = TRUE} with no new device opened.
#' @export
#' @method plot estimateEffect
plot.estimateEffect <- function(x, covariate, model=NULL,
                                topics=x$topics,
                                method=c("pointestimate", "difference","continuous"),
                                ref=NULL, alt=NULL,
                                moderator=NULL, moderator.value=NULL,
                                npoints=100, nsims=500, ci.level=.95,
                                xlim=NULL, ylim=NULL, xlab="",ylab=NULL,
                                main="", printlegend=T,
                                labeltype="numbers", n=7, frexw=.5,
                                add=F, linecol="black", width=25,
                                verbose.labels=T, family=NULL,
                                custom.labels=NULL, omit.plot=FALSE,...){

  method <- match.arg(method)
  if(method=="difference" && (is.null(ref) | is.null(alt))) {
    # stop("For method='difference' both ref and alt must be specified.")
    if(! covariate %in% names(x$ref.vec)) stop("For method='difference', either the covariate needs to be an covariate specified in the `estimateEffect`, or both ref and alt must be specified.")
    ref.vec <- x$ref.vec
    ref <- ref.vec[match(covariate, names(ref.vec))]
    metadata <- colData(model$settings$sce)
    alt <- setdiff(metadata[[covariate]], ref)
    if(length(alt) != 1) stop("Both ref and alt must be specified, since more than 2 levels are found in metadata")
  }

  #Produce cdata (data in original form) and
  #cmatrix (data in design matrix form)
  cthis <- produce_cmatrix(prep=x, covariate=covariate, method=method,
                           ref=ref,
                           alt=alt, npoints=npoints,
                           moderator=moderator, moderator.value=moderator.value)
  cdata <- cthis$cdata
  cmat <- cthis$cmatrix
  #Simulate betas
  simbetas <- simBetas(x$parameters, nsims=nsims)

  #Find offset for confidence level
  offset <- (1-ci.level)/2


  #Plot for each method
  if(method=="continuous"){
    toreturn <- plotContinuous(prep=x,covariate=covariate,topics=topics, cdata=cdata, cmat=cmat, simbetas=simbetas,
                   offset=offset,xlab=xlab, ylab=ylab, main=main,
                   xlim=xlim, ylim=ylim, linecol=linecol, add=add,
                   labeltype=labeltype,n=n,custom.labels=custom.labels,model=model,frexw=frexw,printlegend=printlegend,omit.plot=omit.plot,...)
    return(invisible(toreturn))
  }
  if(method=="pointestimate"){
    toreturn <- plotPointEstimate(prep=x,covariate=covariate,topics=topics, cdata=cdata, cmat=cmat, simbetas=simbetas,
                      offset=offset,xlab=xlab, ylab=ylab, main=main,
                      xlim=xlim, ylim=ylim, linecol=linecol, add=add,
                      labeltype=labeltype,n=n,
                                  custom.labels=custom.labels,model=model,frexw=frexw,width=width,
                                  verbose.labels=verbose.labels,omit.plot=omit.plot,...)
    return(invisible(toreturn))
  }
  if(method=="difference"){
    if(missing(ref)) stop("Missing a value for ref. See documentation.")
    if(missing(alt)) stop("Missing a value for alt. See documentation.")
    toreturn <- plotDifference(prep=x,covariate=covariate,topics=topics, cdata=cdata, cmat=cmat, simbetas=simbetas,
                   offset=offset,xlab=xlab, ylab=ylab, main=main,
                   xlim=xlim, ylim=ylim, linecol=linecol, add=add,
                   labeltype=labeltype,n=n,
                   custom.labels=custom.labels, printlegend=printlegend,
                   model=model,frexw=frexw,width=width,
                   ref=ref,
                               alt=alt,verbose.labels=verbose.labels,omit.plot=omit.plot,...)
    return(invisible(toreturn))
  }
}
