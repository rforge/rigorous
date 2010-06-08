##
##
## Copyright (c) 2009, Brandon Whitcher and Volker Schmid
## All rights reserved.
## 
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are
## met:
## 
##     * Redistributions of source code must retain the above copyright
##       notice, this list of conditions and the following disclaimer. 
##     * Redistributions in binary form must reproduce the above
##       copyright notice, this list of conditions and the following
##       disclaimer in the documentation and/or other materials provided
##       with the distribution.
##     * The names of the authors may not be used to endorse or promote
##       products derived from this software without specific prior
##       written permission.
## 
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
## "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
## LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
## A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
## HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
## SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
## LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
## DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
## THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
## (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
## OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
## 
## $Id: plot.R 332 2010-01-29 16:54:07Z bjw34032 $
##

#############################################################################
## image() for class="nifti"
## promptMethods(image, "image-methods.Rd")
#############################################################################

image.nifti <- function(x, z=1, w=1, col=gray(0:64/64),
                        plot.type=c("multiple","single"), zlim=NULL,
                        xlab="", ylab="", axes=FALSE, oma=rep(0,4),
                        mar=rep(0,4), bg="black", ...) {
  ## set dimensions
  X <- nrow(x)
  Y <- ncol(x)
  Z <- nsli(x)
  W <- ntim(x)
  ## check dimensions
  if (X == 0 || Y == 0 || Z == 0)
    stop("size of NIfTI volume is zero, nothing to plot")
  if (z < 1 || z > Z) {
    stop("slice \"z\" out of range")
  }
  ## check for z-limits; use internal by default
  if (is.null(zlim)) {
    zlim <- c(x@"cal_min", x@"cal_max")
    if (max(zlim) == 0) {
      zlim <- c(x@"glmin", x@"glmax")
    }
  }
  breaks <- c(min(x,zlim),
              seq(min(zlim), max(zlim), length=length(col)-1),
              max(x,zlim))
  ## single or multiple images?
  if (plot.type[1] == "multiple") {
    index <- 1:Z
  } else {
    index <- z
  }
  lz <- length(index)
  ## plotting
  oldpar <- par(no.readonly=TRUE)
  par(mfrow=ceiling(rep(sqrt(lz),2)), oma=oma, mar=mar, bg=bg)
  if (is.na(W)) { # three-dimensional array
    for (z in index) {
      graphics::image(1:X, 1:Y, x[,,z], col=col, breaks=breaks,
                      axes=axes, xlab=xlab, ylab=ylab, ...)
    }
  } else { # four-dimensional array
    if (w < 1 || w > W)
      stop("volume \"w\" out of range")
    for (z in index) {
      graphics::image(1:X, 1:Y, x[,,z,w], col=col, breaks=breaks,
                      axes=axes, xlab=xlab, ylab=ylab, ...)
    }
  }
  par(oldpar)
  invisible()
}

setMethod("image", signature(x="nifti"), image.nifti)
setMethod("image", signature(x="anlz"), image.nifti)

#############################################################################
## overlay() for class="nifti"
#############################################################################

overlay.nifti <- function(x, y, z=1, w=1, col.x=gray(0:64/64),
                          col.y=hotmetal(), zlim.x=NULL, zlim.y=NULL,
                          plot.type=c("multiple","single"),
                          xlab="", ylab="", axes=FALSE, oma=rep(0,4),
                          mar=rep(0,4), bg="black", ...) {
  ## both volumes must have the same dimension
  if (! all(dim(x)[1:3] == dim(y)[1:3]))
    stop("dimensions of \"x\" and \"y\" must be equal")
  ## set dimensions
  X <- nrow(x)
  Y <- ncol(x)
  Z <- nsli(x)
  W <- ntim(x)
  ## check dimensions
  if (X == 0 || Y == 0 || Z == 0)
    stop("size of NIfTI volume is zero, nothing to plot")
  ## check for z-limits in x; use internal by default
  if (is.null(zlim.x)) {
    zlim.x <- c(x@"cal_min", x@"cal_max")
    if (max(zlim.x) == 0)
      zlim.x <- c(x@"glmin", x@"glmax")
  }
  breaks.x <- c(min(x,zlim.x,na.rm=TRUE),
                seq(min(zlim.x,na.rm=TRUE), max(zlim.x,na.rm=TRUE),
                    length=length(col.x)-1),
                max(x,zlim.x,na.rm=TRUE))
  ## check for z-limits in y; use internal by default
  if (is.null(zlim.y)) {
    zlim.y <- c(y@"cal_min", y@"cal_max")
    if (max(zlim.y) == 0)
      zlim.y <- c(x@"glmin", x@"glmax")
  }
  if (plot.type[1] == "multiple") {
    index <- 1:Z
  } else {
    index <- z
  }
  lz <- length(index)
  if (z < 1 || z > Z)
    stop("slice \"z\" out of range")
  oldpar <- par(no.readonly=TRUE)
  par(mfrow=ceiling(rep(sqrt(lz),2)), oma=oma, mar=mar, bg=bg)
  if (is.na(W)) { # three-dimensional array
    for (z in index) {
      graphics::image(1:X, 1:Y, x[,,z], col=col.x, breaks=breaks.x,
                      zlim=zlim.x, axes=axes, xlab=xlab, ylab=ylab, ...)
      graphics::image(1:X, 1:Y, y[,,z], col=col.y, zlim=zlim.y, add=TRUE)
    }
  } else { # four-dimensional array
    if (w < 1 || w > W)
      stop("volume \"w\" out of range")
    for (z in index) {
      graphics::image(1:X, 1:Y, x[,,z,w], col=col.x, breaks=breaks.x,
                      zlim=zlim.x, axes=axes, xlab=xlab, ylab=ylab, ...)
      graphics::image(1:X, 1:Y, y[,,z], col=col.y, zlim=zlim.y, add=TRUE)
    }
  }
  par(oldpar)
  invisible()
}

setGeneric("overlay", function(x, y, ...) standardGeneric("overlay"))
setMethod("overlay", signature(x="nifti", y="nifti"), overlay.nifti)
setMethod("overlay", signature(x="anlz", y="anlz"), overlay.nifti)
setMethod("overlay", signature(x="anlz", y="nifti"), overlay.nifti)
setMethod("overlay", signature(x="nifti", y="anlz"), overlay.nifti)
setMethod("overlay", signature(x="array", y="array"), overlay.nifti)
setMethod("overlay", signature(x="array", y="nifti"), overlay.nifti)
setMethod("overlay", signature(x="nifti", y="array"), overlay.nifti)
setMethod("overlay", signature(x="array", y="anlz"), overlay.nifti)
setMethod("overlay", signature(x="anlz", y="array"), overlay.nifti)

#############################################################################
## orthographic() for class="nifti"
#############################################################################

orthographic.nifti <- function(x, xyz=NULL, crosshairs=TRUE,
                               col.crosshairs="red", w=1, zlim=NULL,
                               col=gray(0:64/64), xlab="", ylab="",
                               axes=FALSE, oma=rep(0,4), mar=rep(0,4),
                               bg="black", text=NULL, text.color="white",
                               ...) {
  X <- nrow(x)
  Y <- ncol(x)
  Z <- nsli(x)
  W <- ntim(x)
  ## Center crosshairs if not specified
  if (is.null(xyz))
    xyz <- ceiling(c(X,Y,Z)/2)
  ## check dimensions
  if (X == 0 || Y == 0 || Z == 0)
    stop("size of NIfTI volume is zero, nothing to plot")
  ## check for z-limits in x; use internal by default
  if (is.null(zlim)) {
    zlim <- c(x@"cal_min", x@"cal_max")
    if (diff(zlim) == 0)
      zlim <- c(x@"glmin", x@"glmax")
    if (diff(zlim) == 0)
      zlim <- range(x, na.rm=TRUE)
  }
  breaks <- c(min(x,zlim),
              seq(min(zlim), max(zlim), length=length(col)-1),
              max(x,zlim))
  oldpar <- par(no.readonly=TRUE)
  par(mfrow=c(2,2), oma=oma, mar=mar, bg=bg)
  if (is.na(W)) {
    ## Three-dimensional array
    graphics::image(1:X, 1:Z, x[,xyz[2],], col=col, breaks=breaks,
                    asp=x@pixdim[4]/x@pixdim[2],
                    xlab=ylab, ylab=xlab, axes=axes, ...)
    if (crosshairs)
      abline(h=xyz[3], v=xyz[1], col=col.crosshairs)
    graphics::image(1:Y, 1:Z, x[xyz[1],,], col=col, breaks=breaks,
                    asp=x@pixdim[4]/x@pixdim[3],
                    xlab=xlab, ylab=ylab, axes=axes, ...)
    if (crosshairs)
      abline(h=xyz[3], v=xyz[2], col=col.crosshairs)
    graphics::image(1:X, 1:Y, x[,,xyz[3]], col=col, breaks=breaks,
                    asp=x@pixdim[3]/x@pixdim[2],
                    xlab=xlab, ylab=ylab, axes=axes, ...)
    if (crosshairs)
      abline(h=xyz[2], v=xyz[1], col=col.crosshairs)
  } else {
    ## Four-dimensional array    
    if (w < 1 || w > W)
      stop("volume \"w\" out of range")
    graphics::image(1:X, 1:Z, x[,xyz[2],,w], col=col, breaks=breaks,
                    asp=x@pixdim[4]/x@pixdim[2],
                    xlab=ylab, ylab=xlab, axes=axes, ...)
    if (crosshairs)
      abline(h=xyz[3], v=xyz[1], col=col.crosshairs)
    graphics::image(1:Y, 1:Z, x[xyz[1],,,w], col=col, breaks=breaks,
                    asp=x@pixdim[4]/x@pixdim[3],
                    xlab=xlab, ylab=ylab, axes=axes, ...)
    if (crosshairs)
      abline(h=xyz[3], v=xyz[2], col=col.crosshairs)
    graphics::image(1:X, 1:Y, x[,,xyz[3],w], col=col, breaks=breaks,
                    asp=x@pixdim[3]/x@pixdim[2],
                    xlab=xlab, ylab=ylab, axes=axes, ...)
    if (crosshairs)
      abline(h=xyz[2], v=xyz[1], col=col.crosshairs)
  }
  if (! is.null(text)) {
    ## Add user-supplied text to the "fourth" plot
    graphics::image(1:64, 1:64, matrix(NA, 64, 64))
    text(32, 32, text, col=text.color, cex=2)
  }
  par(oldpar)
  invisible()
}

setGeneric("orthographic", function(x, ...) standardGeneric("orthographic"))
setMethod("orthographic", signature(x="nifti"), orthographic.nifti)
setMethod("orthographic", signature(x="anlz"), orthographic.nifti)
setMethod("orthographic", signature(x="array"),
          function(x) {
            x <- as(x, "nifti")
            orthographic.nifti(x)
          })
