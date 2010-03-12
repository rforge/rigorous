##
## Copyright (c) 2010, Brandon Whitcher
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
##     * Neither the name of Rigorous Analytics Ltd. nor the names of
##       its contributors may be used to endorse or promote products 
##       derived from this software without specific prior written 
##       permission.
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
## $Id: $
##

imageTensor <- function(mat3d, weights=NULL, plot=TRUE, add=FALSE,
                        exclude=1, ...) {
  ##
  ## Adapted from _imagenrgb_ by Dr. Agustin Lobo
  ## Instituto de Ciencias de la Tierra (CSIC)
  ## Lluis Sole Sabaris s/n, 08028 Barcelona SPAIN
  ## tel 34 93409 5410
  ## fax 34 93411 0012
  ## alobo@ija.csic.es 
  ##
  
  ## Displays a (m,n,3) array as an rgb image.
  ## If "plot = FALSE", saves the pseudocolor image as a list.
  
  reclass <- function(matriz, origen, imagen, directo=TRUE) {
    if (!directo) {
      aux <- origen
      origen <- imagen
      imagen <- aux
    }
    ## As suggested by P.B. Ripley:
    M <- match(matriz, origen, 0)
    matriz[M > 0] <- imagen[M]
    return(matriz)
  }

  m <- nrow(mat3d)
  n <- ncol(mat3d)

  ## 1. Omitted.

  ## 2. Generates z vectors from a (m,n,3) array.
  triplets <- matrix(c(mat3d), m*n, 3)
  if (!is.null(weights)) {
    triplets <- c(weights)/max(weights, na.rm=TRUE) * triplets
  }
  triplets[is.na(triplets)] <- 0
  triplets[triplets < 0] <- 0
  
  ## 3. Generates RGB colors:
  cols <- rgb(triplets[,1], triplets[,2], triplets[,3])
  ## Generates vector of unique colors:
  cols.unicos <- unique(cols)
  Nc <- length(cols.unicos)
  ## Assigns an integer code to each unique color and transforms the
  ## char color matrix into an integer matrix:
  cols <- as.numeric(reclass(cols, cols.unicos, 1:Nc), drop=FALSE)
  ## Formats vector of color codes as (m,n) matrix:
  cols <- matrix(cols, m, n)
  
  ## 4. Display or save
  if (plot) {
    image(1:m, 1:n, ifelse(cols != exclude, cols, NA), col=cols.unicos,
          add=add, axes=FALSE, ...)
  } else {
    list(ima=cols, cols=cols.unicos)
  }
}

