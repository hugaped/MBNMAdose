# Functions for plotting in MBNMAdose
# Author: Hugo Pedder
# Date created: 2019-04-18

## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c(".", "studyID", "agent", "dose", "Var1", "value",
                                                        "Parameter", "fupdose", "groupvar", "y",
                                                        "network", "a", "param", "med", "l95", "u95", "value",
                                                        "Estimate", "2.5%", "50%", "97.5%", "treatment"))




#' Calculate position of label with respect to vertex location within a circle
#'
#' Useful for graphs drawn using `igraph` to reposition labels relative to vertices when vertices
#' are laid out in a circle (as is common in network plots). `igraph` interprets position within
#' `vertex.label.degree` as radians, so it is necessary to convert locations into radian values. This
#' is the main role of this function.
#'
#' @param x A numeric vector of positions around a circle, typically sequentially numbered.
#' @param start A number giving the offset from 12 o'clock in radians for the label locations.
#' @param direction Either `1` for clockwise numbering (based on the order of `x`) or `-1` for
#' anti-clockwise.
#'
#' @examples
#' radian.rescale(c(1:10), start=0, direction=1)
#'
#' @noRd
#'
#' @references
#' https://gist.github.com/kjhealy/834774/a4e677401fd6e4c319135dabeaf9894393f9392c
radian.rescale <- function(x, start=0, direction=1) {
  c.rotate <- function(x) (x + start) %% (2 * pi) * direction
  c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
}




#' Get large vector of distinct colours using Rcolorbrewer
#' @noRd
genmaxcols <- function() {

  cols1 <- RColorBrewer::brewer.pal(9, "Set1")
  cols2 <- RColorBrewer::brewer.pal(9, "Pastel1")
  cols <- c(rbind(cols1, cols2))

  cols1 <- RColorBrewer::brewer.pal(8, "Set2")
  cols2 <- RColorBrewer::brewer.pal(8, "Pastel2")
  cols <- c(cols, c(rbind(cols1, cols2)))

  cols <- c(cols, RColorBrewer::brewer.pal(12, "Set3"))

  return(cols)
}





#' Overlays observations as shaded regions on a time-course
#'
#' @inheritParams mbnma.run
#' @inheritParams plot.mbnma.predict
#' @param g An object of `class("ggplot")`
#' @param max.col.scale The maximum rgb numeric value to use for the colour scale
#' @noRd
disp.obs <- function(g, network, predict, col="red", max.col.scale=NULL) {
  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(g, c("gg", "ggplot"), add=argcheck)
  checkmate::assertClass(network, "mbnma.network", add=argcheck)
  checkmate::assertClass(predict, "mbnma.predict", add=argcheck)
  checkmate::assertInt(max.col.scale, lower=1, null.ok = TRUE, add=argcheck)
  checkmate::reportAssertions(argcheck)


  raw.data <- network[["data.ab"]]
  raw.data$agent <- factor(raw.data$agent, labels=network[["agents"]])
  predict.data <- summary(predict)
  predict.data[["count"]] <- NA
  predict.data[["cum.count"]] <- NA

  # Identify counts of doses in raw data
  for (i in 1:nrow(predict.data)) {
    if (i==1) {
      dose.low <- -0.01
    } else if (predict.data$agent[i-1]!=predict.data$agent[i]) {
      dose.low <- -0.01
    } else {
      dose.low <- predict.data$dose[i-1]
    }

    dose.high <- predict.data$dose[i]
    p.agent <- predict.data$agent[i]

    count <- length(raw.data$studyID[as.character(raw.data$agent)==p.agent &
                                     raw.data$dose>dose.low & raw.data$dose<=dose.high])
    cum.count <- length(raw.data$studyID[as.character(raw.data$agent)==p.agent &
                                       raw.data$dose<=dose.high])

    predict.data$count[i] <- count
    predict.data$cum.count[i] <- cum.count
  }

  # Identify if placebo is in data
  if (any(unique(predict.data$agent) %in% c("1", "Placebo"))) {
    plac.count <- predict.data$count[predict.data$agent %in% c("1", "Placebo")]
    message(paste0(plac.count, " placebo arms in the dataset are not shown within the plots"))
    agents <- unique(predict.data$agent)[!(unique(predict.data$agent) %in% c("1", "Placebo"))]
  } else {
    agents <- unique(predict.data$agent)
  }


  # Check max.col.scale
  n.cut <- max(predict.data$count[predict.data$dose!=0])
  if (!is.null(max.col.scale)) {
    if (!is.numeric(max.col.scale)) {
      stop("max.col.scale must be a number greater than the maximum number of observed counts in the plotted treatments")
    }
    if (n.cut > max.col.scale) {
      stop("max.col.scale must be a number greater than the maximum number of observed counts in the plotted treatments")
    }
    n.cut <- max.col.scale
  }

  # Generate colours
  cols <- alpha.scale(n.cut=n.cut, col=col)

  for (agent in seq_along(agents)) {
    subset <- predict.data[predict.data$agent==agents[agent],]
    subset$agent <- factor(subset$agent, labels=levels(g$data$agent)[agents[agent]])

    # Start assuming lowest time = 0
    for (m in 2:nrow(subset)) {
      g <- g + ggplot2::geom_ribbon(data=subset[subset$dose<=subset$dose[m] &
                                                  subset$dose>=subset$dose[m-1]
                                                ,],
                                    ggplot2::aes(x=dose,
                                                 ymin=`2.5%`,
                                                 ymax=`97.5%`),
                                    fill=cols[subset$count[m]+1])
    }

  }

  return(list(g, cols))
}




#' Generates colours with varying degrees of transparency
#'
#' Identical to function in `MBNMAtime` package
#'
#' @param ncut A number indicating the number of different counts in the dataset
#' @param col Colour to use for shading
#' @noRd
alpha.scale <- function(n.cut, col="blue") {
  # Run checks
  checkmate::assertIntegerish(n.cut, lower=0, len=1)

  # Set colour intensities
  if (is.character(col)) {
    if (col=="blue") {
      hue <- c(0,0,200)
    } else if (col=="red") {
      hue <- c(120,0,0)
    } else if (col=="green") {
      hue <- c(0,100,0)
    } else {
      stop("Permitted colours are `blue`, `red`, `green`, or an RGB code")
    }
  } else if (is.numeric(col)) {
    if (length(col)!=3) {
      stop("Specified RGB code must have length 3")
    }
    if (any(col>255) | any(col<0)) {
      stop("RGB values must be between 0 and 255")
    }
    hue <- col
  }

  cut <- (255-0)/(n.cut)
  alpha <- 0
  alpha.vec <- alpha

  for (i in 1:n.cut) {
    alpha <- alpha+cut
    alpha.vec <- append(alpha.vec, alpha)
  }

  hexcol <- vector()
  for (i in 1:(n.cut+1)) {
    hexcol <- append(hexcol, grDevices::rgb(hue[1], hue[2],
                                 hue[3], alpha.vec[i], maxColorValue=255
    ))
  }

  return(hexcol)
}




#' Overlays results of split (treatment-level) NMA on a plot
#'
#' Requires automatically running an NMA
#'
#' @inheritParams plot.mbnma.predict
#' @inheritParams predict.mbnma
#' @inheritParams disp.obs
#' @noRd
overlay.split <- function(g, network, E0=unique(g$data$`50%`[g$data$dose==0]), method="common",
                          likelihood="binomial", link="logit", ...) {

  # Check/assign link and likelihood
  likelink <- check.likelink(network$data.ab, likelihood=likelihood, link=link)
  likelihood <- likelink[["likelihood"]]
  link <- likelink[["link"]]

  splitNMA <- nma.run(network=network, method=method,
                      likelihood=likelihood, link=link, ...)

  splitmat <- splitNMA[["jagsresult"]]$BUGSoutput$sims.list$d

  if (splitNMA[["trt.labs"]][1]!="Placebo_0") {
    warning("Split NMA network does not include placebo:\nresults cannot be displayed over MBNMA predictions")
    return(g)
  }

  trtvec <- splitNMA[["trt.labs"]]
  agentvec <- as.vector(sapply(splitNMA[["trt.labs"]],
                           function(x) strsplit(x, split="_", fixed=TRUE)[[1]][1]))
  dosevec <- as.vector(as.numeric(sapply(splitNMA[["trt.labs"]],
                                     function(x) strsplit(x, split="_", fixed=TRUE)[[1]][2])))
  agentvec <- factor(agentvec, levels=network$agents)

  # Ensure only agents whose predictions are plotted are included
  splitmat <- splitmat[, agentvec %in% g$data$agent]
  dosevec <- dosevec[agentvec %in% g$data$agent]
  trtvec <- trtvec[agentvec %in% g$data$agent]
  agentvec <- agentvec[agentvec %in% g$data$agent]


  # Ensure E0 is correct length to add to splitmat
  if (length(E0)<nrow(splitmat)) {
    E0 <- rep(E0, length.out=nrow(splitmat))
  } else if (length(E0)>nrow(splitmat)) {
    E0 <- replicate(nrow(splitmat), sample(E0))
  }

  splitmat <- splitmat + E0

  quants <- apply(splitmat, MARGIN=2, FUN=function(x) {
    stats::quantile(x, probs=c(0.025,0.5,0.975))
  })

  # Return to natural scale
  quants <- rescale.link(quants, direction="natural", link=link)

  split.df <- data.frame("agent"=agentvec, "dose"=dosevec, "treatment"=trtvec,
                         "2.5%"=quants[1,], "50%"=quants[2,], "97.5%"=quants[3,])

  names(split.df)[4:6] <- c("2.5%", "50%", "97.5%")


  # Add split NMAs to plot
  g <- g + ggplot2::geom_point(data=split.df, ggplot2::aes(x=dose, y=`50%`)) +
    ggplot2::geom_errorbar(data=split.df, ggplot2::aes(x=dose, ymin=`2.5%`, ymax=`97.5%`))


  #### Report split NMA model fit statistics ####
  nma.msg <- vector()
  nma.msg <- c(nma.msg, paste0("Split NMA residual Deviance: ",
                               round(splitNMA$jagsresult$BUGSoutput$median$totresdev, 1)))

  if (method=="random") {
    sd <- splitNMA$jagsresult$BUGSoutput$summary[rownames(splitNMA$jagsresult$BUGSoutput$summary)=="sd",]
    sd <- paste0(signif(sd[5], 3), " (", signif(sd[3], 3), ", ", signif(sd[7], 3), ")")

    nma.msg <- c(nma.msg, paste0("Split NMA between-study SD: ", sd
                                 ))
  }

  message(cat(paste(nma.msg, collapse="\n")))

  return(g)

}







#' Plot deviance contributions from an MBNMA model
#'
#' @param mbnma An S3 object of class `"mbnma"` generated by running
#'   a dose-response MBNMA model
#' @param dev.type *STILL IN DEVELOPMENT FOR MBNMAdose!* Deviances to plot - can be either residual
#' deviances (`"resdev"`, the default) or deviances (`"dev"`)
#' @param plot.type Deviances can be plotted either as scatter points (`"scatter"` - the default)
#' or as boxplots (`"box"`)
#' @param facet A boolean object that indicates whether or not to facet (by agent for `MBNMAdose`
#' and by treatment for `MBNMAtime`)
#' @param ... Arguments to be sent to `ggplot2::geom_point()` or `ggplot2::geom_boxplot`
#' @inheritParams predict.mbnma
#' @inheritParams R2jags::jags
#'
#' @return Generates a plot of deviance contributions and returns a list containing the
#' plot (as an object of `class(c("gg", "ggplot"))`), and a data.frame of posterior mean
#' deviance/residual deviance contributions for each observation.
#'
#' @details
#' Deviances should only be plotted for models that have converged successfully. If deviance
#' contributions have not been monitored in `mbnma$parameters.to.save` then additional
#' iterations will have to be run to get results for these.
#'
#' For `MBNMAtime`, deviance contributions cannot be calculated for models with a multivariate likelihood (i.e.
#' those that account for correlation between observations) because the covariance matrix in these
#' models is treated as unknown (if `rho = "estimate"`) and deviance contributions will be correlated.
#'
#' @examples
#' \donttest{
#' # Using the triptans data
#' network <- mbnma.network(HF2PPITT)
#'
#' # Run an Emax dose-response MBNMA and predict responses
#' emax <- mbnma.emax(network, method="random")
#'
#' # Plot deviances
#' devplot(emax)
#'
#' # Plot deviances using boxplots
#' devplot(emax, plot.type="box")
#'
#' # Plot deviances on a single plot (not facetted by agent)
#' devplot(emax, facet=FALSE, plot.type="box")
#'
#' # A data frame of deviance contributions can be obtained from the object
#' #returned by `devplot`
#' devs <- devplot(emax)
#' head(devs$dev.data)
#'
#' # Other deviance contributions not currently implemented but in future
#' #it will be possible to plot them like so
#' #devplot(emax, dev.type="dev")
#' }
#'
#' @export
devplot <- function(mbnma, plot.type="scatter", facet=TRUE, dev.type="resdev",
                    n.iter=mbnma$BUGSoutput$n.iter, n.thin=mbnma$BUGSoutput$n.thin,
                    ...) {
  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(mbnma, "mbnma", add=argcheck)
  checkmate::assertChoice(dev.type, choices = c("dev", "resdev"), add=argcheck)
  checkmate::assertChoice(plot.type, choices = c("scatter", "box"), add=argcheck)
  checkmate::assertInt(n.iter, lower=1, null.ok = TRUE, add=argcheck)
  checkmate::assertInt(n.thin, lower=1, upper=n.iter, null.ok = TRUE, add=argcheck)
  checkmate::reportAssertions(argcheck)

  if (!is.null(mbnma$model.arg$rho)) {
    stop("Correlation between time points has been modelled using a multivariate likelihood.
         Deviances cannot be calculated for this model.")
  }

  if (!(dev.type %in% mbnma$parameters.to.save)) {
    msg <- paste0("`", dev.type, "` not monitored in mbnma$parameters.to.save.\n",
                  "additional iterations will be run in order to obtain results for `", dev.type, "`")
    message(msg)

    dev.df <- mbnma.update(mbnma, param=dev.type, n.iter=n.iter, n.thin=n.thin)

  } else {

    dev.df <- get.theta.dev(mbnma, param=dev.type)
  }

  # Add studyID in addition to study index from model and sort data frame
  dev.df$studyID <- mbnma$model.arg$jagsdata$studyID
  dev.df <- dplyr::arrange(dev.df, dev.df$study, dev.df$arm)

  # Plots the residual deviances
  if (mbnma$type=="time") {
    xlab <- "Follow-up count"
    facetscale <- "fixed"
  } else if (mbnma$type=="dose") {
    agents <- mbnma$network$agents

    # Remove placebo results if they are present
    if (mbnma$network$agents[1]=="Placebo") {
      dev.df <- dev.df[dev.df$facet!=1,]
      agents <- agents[-1]
    }

    xlab <- "Dose"
    facetscale <- "free_x"
  }

  if ("agents" %in% names(mbnma$network)) {
    dev.df$facet <- factor(dev.df$facet, labels=agents)
  } else if ("treatments" %in% names(mbnma$network)) {
    dev.df$facet <- factor(dev.df$facet, labels=mbnma$treatments)
  }

  if (plot.type=="scatter") {
    g <- ggplot2::ggplot(dev.df, ggplot2::aes(x=fupdose, y=mean), group=groupvar) +
      ggplot2::geom_point(...)
  } else if (plot.type=="box") {
    g <- ggplot2::ggplot(dev.df, ggplot2::aes(x=factor(fupdose), y=mean)) +
      ggplot2::geom_boxplot(...)
  }

  # Add axis labels
  g <- g + ggplot2::xlab(xlab) +
    ggplot2::ylab("Posterior mean") +
    ggplot2::theme_bw()

  if (facet==TRUE) {
    g <- g + ggplot2::facet_wrap(~facet, scales = facetscale) +
      ggplot2::expand_limits(x=0)
  }

  graphics::plot(g)
  return(invisible(list("graph"=g, "dev.data"=dev.df)))
}





#' Extracts fitted values or deviance contributions into a data.frame with indices
#'
#' @inheritParams predict.mbnma
#' @param The parameter for which to extract values - can take either `"theta"`, `"dev"` or `"resdev"`
#' @noRd
get.theta.dev <- function(mbnma, param="theta") {
  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(mbnma, "mbnma", add=argcheck)
  checkmate::assertChoice(param, choices = c("dev", "resdev", "theta"), add=argcheck)
  checkmate::reportAssertions(argcheck)

  dev.df <- as.data.frame(
    mbnma$BUGSoutput$summary[grep(paste0("^", param, "\\["),
                                  rownames(mbnma$BUGSoutput$summary)),]
  )

  if (mbnma$type=="time") {
    # Takes the study, arm and follow-up measurement identifier for each residual deviance point
    id <- matrix(unlist(strsplit(
      gsub(
        "([a-z]+\\[)([0-9]+)(,)([0-9]+)(,)([0-9]+)(\\])",
        "\\2.\\4.\\6", rownames(dev.df)),
      split="\\.")),
      byrow=TRUE,
      ncol=3
    )

    dev.df$fupdose <- as.numeric(id[,3])
    dev.df$groupvar <- paste(as.numeric(id[,1]), as.numeric(id[,2]), sep="_")

    # Treatment as facet
    temp <- replicate(max(dev.df$fupdose), mbnma$model$data()$treatment)
    dev.df$facet <- as.vector(temp)[
      stats::complete.cases(as.vector(temp))
      ]

  } else if (mbnma$type=="dose") {
    # Takes the study, arm and follow-up measurement identifier for each residual deviance point
    id <- matrix(unlist(strsplit(
      gsub(
        "([a-z]+\\[)([0-9]+)(,)([0-9]+)(\\])",
        "\\2.\\4", rownames(dev.df)),
      split="\\.")),
      byrow=TRUE,
      ncol=2
    )

    # Agent as facet
    dev.df$facet <- as.vector(mbnma$model.arg$jagsdata$agent)[
      stats::complete.cases(as.vector(mbnma$model.arg$jagsdata$agent))
      ]

    dev.df$fupdose <- as.vector(mbnma$model.arg$jagsdata$dose)[
      stats::complete.cases(as.vector(mbnma$model.arg$jagsdata$dose))
      ]
    dev.df$groupvar <- as.numeric(id[,1])
  }

  dev.df$study <- as.numeric(id[,1])
  dev.df$arm <- as.numeric(id[,2])


  return(dev.df)
}





#' Plot fitted values from MBNMA model
#'
#' @param mbnma An S3 object of class `"mbnma"` generated by running
#'   a dose-response MBNMA model
#' @param disp.obs A boolean object to indicate whether raw data responses should be
#' plotted as points on the graph
#' @param ... Arguments to be sent to `ggplot2::geom_point()` or `ggplot2::geom_line()`
#' @inheritParams predict.mbnma
#' @inheritParams R2jags::jags
#'
#' @return Generates a plot of fitted values from the MBNMA model and returns a list containing
#' the plot (as an object of `class(c("gg", "ggplot"))`), and a data.frame of posterior mean
#' fitted values for each observation.
#'
#' @details
#' Fitted values should only be plotted for models that have converged successfully.
#' If fitted values (`theta`) have not been monitored in `mbnma$parameters.to.save`
#' then additional iterations will have to be run to get results for these.
#'
#' @examples
#' \donttest{
#' # Using the triptans data
#' network <- mbnma.network(HF2PPITT)
#'
#' # Run an Emax dose-response MBNMA and predict responses
#' emax <- mbnma.emax(network, method="random")
#'
#' # Plot fitted values and observed values
#' fitplot(emax)
#'
#' # Plot fitted values only
#' fitplot(emax, disp.obs=FALSE)
#'
#' # A data frame of fitted values can be obtained from the object
#' #returned by `fitplot`
#' fits <- fitplot(emax)
#' head(fits$fv)
#' }
#'
#' @export
fitplot <- function(mbnma, disp.obs=TRUE,
                    n.iter=mbnma$BUGSoutput$n.iter, n.thin=mbnma$BUGSoutput$n.thin, ...) {
  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(mbnma, "mbnma", add=argcheck)
  checkmate::assertLogical(disp.obs, add=argcheck)
  checkmate::assertInt(n.iter, lower=1, null.ok = TRUE, add=argcheck)
  checkmate::assertInt(n.thin, lower=1, upper=n.iter, null.ok = TRUE, add=argcheck)
  checkmate::reportAssertions(argcheck)

  if (!("theta" %in% mbnma$parameters.to.save)) {
    msg <- paste0("`theta` not monitored in mbnma$parameters.to.save.\n",
                  "additional iterations will be run in order to obtain results")
    message(msg)

    theta.df <- mbnma.update(mbnma, param="theta", n.iter=n.iter, n.thin=n.thin)
  } else {
    theta.df <- get.theta.dev(mbnma, param="theta")
  }

  # Add studyID in addition to study index from model
  theta.df$studyID <- mbnma$model.arg$jagsdata$studyID


  # Obtain raw responses to plot over fitted
  if (mbnma$model.arg$likelihood=="binomial") {
    theta <- mbnma$model$data()$r / mbnma$model$data()$N
  } else if (mbnma$model.arg$likelihood=="poisson") {
    theta <- mbnma$model$data()$r / mbnma$model$data()$E
  } else if (mbnma$model.arg$likelihood=="normal") {
    theta <- mbnma$model$data()$y
  }

  theta <- rescale.link(theta, direction="link", link=mbnma$model.arg$link)
  theta.df$y <- as.vector(theta)[stats::complete.cases(as.vector(theta))]

  # Add E0 for each agent and drop dose==0
  E0 <- mean(theta.df$mean[theta.df$fupdose==0])
  theta.df <- theta.df[theta.df$fupdose!=0,]

  for (i in seq_along(unique(theta.df$study))) {
    subset <- theta.df[theta.df$study==unique(theta.df$study)[i],]
    for (k in seq_along(unique(subset$facet))) {
      row <- subset[subset$facet==unique(subset$facet)[k],][1,]
      row$mean <- E0
      row$y <- NA
      row$fupdose <- 0
      theta.df <- rbind(theta.df, row)
    }
  }

  theta.df <- dplyr::arrange(theta.df, theta.df$study, theta.df$facet, theta.df$fupdose)


  # Axis labels
  if (mbnma$type=="time") {
    xlab <- "Follow-up count"
    facetscale <- "fixed"
  } else if (mbnma$type=="dose") {
    xlab <- "Dose"
    facetscale <- "free_x"
  }
  ylab <- "Response on link scale"

  # Add facet labels
  if ("agents" %in% names(mbnma$network)) {
    if (mbnma$network$agents[1]=="Placebo" & mbnma$network$treatments[1]=="Placebo_0") {
      labs <- mbnma$network$agents[-1]
    } else {labs <- mbnma$network$agents}
    theta.df$facet <- factor(theta.df$facet, labels=labs)
  } else if ("treatments" %in% names(mbnma$network)) {
    labs <- mbnma$network$treatments[-1]
    theta.df$facet <- factor(theta.df$facet, labels=labs)
  }

  # Generate plot
  g <- ggplot2::ggplot(theta.df,
                       ggplot2::aes(x=fupdose, y=mean, group=groupvar)) +
    ggplot2::geom_line(...)

  # Overlay observed responses
  if (disp.obs==TRUE) {
    g <- g + ggplot2::geom_point(ggplot2::aes(y=y), size=1, ...)
  }

  # Add facets
  g <- g + ggplot2::facet_wrap(~facet, scales = facetscale)

  # Add axis labels
  g <- g + ggplot2::xlab(xlab) +
    ggplot2::ylab(ylab) +
    ggplot2::theme_bw()

  suppressWarnings(graphics::plot(g))

  return(invisible(list("graph"=g, "fv"=theta.df)))
}







plot.invisible <- function(...){
  ff <- tempfile()
  grDevices::png(filename=ff)
  res <- graphics::plot(...)
  grDevices::dev.off()
  unlink(ff)
  return(res)
}






#' Plot cumulative ranking curves from MBNMA models
#'
#' @param x An object of class `"mbnma.rank"` generated by `rank.mbnma()`
#' @param sucra A logical object to indicate whether Surface Under Cumulative Ranking Curve (SUCRA)
#' values should be calculated and returned as a data frame. Areas calculated
#' using \code{\link[rgeos]{readWKT}}.
#' @param ... Arguments to be sent to `ggplot::geom_line()`
#' @inheritParams rank.mbnma
#'
#' @return Line plots showing the cumulative ranking probabilities for each agent/class and
#' dose-response parameter in `x`. The object returned is a list which contains the plot
#' (an object of `class(c("gg", "ggplot")`) and a data frame of SUCRA values
#' if `sucra = TRUE`.
#'
#' @examples
#' \donttest{
#' # Using the triptans data
#' network <- mbnma.network(HF2PPITT)
#'
#' # Estimate rankings  from an Emax dose-response MBNMA
#' emax <- mbnma.emax(network, emax="rel", ed50="rel", method="random")
#' ranks <- rank(emax)
#'
#' # Plot cumulative rankings for both dose-response parameters simultaneously
#' # Note that SUCRA values are also returned
#' cumrank(ranks)
#' }
#' @export
cumrank <- function(x, params=NULL, sucra=TRUE, ...) {

  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(x, "mbnma.rank", add=argcheck)
  checkmate::assertCharacter(params, null.ok=TRUE, add=argcheck)
  checkmate::assertLogical(sucra, null.ok=FALSE, add=argcheck)
  checkmate::reportAssertions(argcheck)

  output <- list()

  if (is.null(params)) {
    params <- names(x)
  }

  df <- data.frame()
  for (param in seq_along(params)) {
    if (!params[param] %in% names(x)) {
      stop(paste0(params[param], " is not a ranked parameter in x"))
    }

    cum.mat <- x[[params[param]]]$cum.matrix
    treats <- colnames(cum.mat)

    melt <- reshape2::melt(cum.mat)
    melt$param <- params[param]

    df <- rbind(df, melt)

  }

  df$Parameter <- factor(df$param)

  g <- ggplot2::ggplot(df, ggplot2::aes(x=Var1, y=value, linetype=Parameter, colour=Parameter), ...) +
    ggplot2::geom_line(size=1)

  g <- g + ggplot2::facet_wrap(~factor(Var2)) +
    ggplot2::xlab("Rank (1 = best)") +
    ggplot2::ylab("Cumulative probability") +
    ggplot2::labs(linetype="Parameter", colour="Parameter") +
    ggplot2::theme_bw()

  graphics::plot(g)

  output <- list("cumplot"=g)

  # Calculate AUC
  if (sucra==TRUE) {
    df.auc <- df %>%
      dplyr::group_by(df$Var2, df$param) %>%
      dplyr::do(data.frame(sucra=calcauc(.)))

    df.auc <- dplyr::ungroup(df.auc)

    names(df.auc) <- c("agent", "parameter", "sucra")

    output[["sucra"]] <- df.auc

    print(df.auc)
  }

  return(invisible(output))
}




#' Returns a forest plot of nodesplit results
#'
#' @param ... Optional arguments to be passed to `forestplot::forestplot()`
#' @inheritParams plot.nodesplit
#' @noRd
forest.splits <- function(x, ...) {

  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(x, "nodesplit", add=argcheck)
  checkmate::reportAssertions(argcheck)


  df <- summary.nodesplit(x)

  df$Evidence[df$Evidence %in% c("MBNMA", "NMA")] <- "Overall"

  comps <- lapply(x, FUN=function(y) {paste(y$comparison, collapse=" vs\n")})

  df$Comparison <- unlist(comps)[match(df$Comparison, names(comps))]

  df <- dplyr::arrange(df, Comparison, Evidence)

  names(df)[4:5] <- c("l95", "u95")


  df$p.value <- format(df$p.value)
  df$p.value[df$p.value=="0.000"] <- "<0.001"

  # Remove surplus labels
  df$Comparison[df$Evidence!="Direct"] <- NA
  df$p.value[df$Evidence!="Indirect"] <- NA


  # Add blank between groups
  temp <- df[0,]
  for (i in 1:(nrow(df)/3)) {
    seg <- rbind(df[((i-1)*3+1):(i*3),], rep(NA,ncol(df)), rep(NA,ncol(df)))
    temp <- rbind(temp, seg)
  }
  df <- temp

  # Plot forest plot
  forestplot::forestplot(labeltext=cbind(df$Comparison, df$Evidence, df$p.value),
                         mean=df$Median, lower=df$l95, upper=df$u95,
                         boxsize=0.2, graph.pos=3,
                         xlab="Effect size (95% CrI)", hrzl_lines = TRUE,
                         ...)

}
