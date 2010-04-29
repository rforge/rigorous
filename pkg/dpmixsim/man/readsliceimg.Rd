\name{readsliceimg}
\alias{readsliceimg}
\concept{MRI}
\concept{data}

\title{ Read MRI slice data }
\description{
  \code{readsliceimg} reads MRI and mask data.
}
\usage{
readsliceimg(fbase="t1_pn3_rf0", swap=FALSE)
}
\arguments{
  \item{fbase}{ Indicates the dataset prefix of the MRI dataset to use.
  The prefix applies to data files: \file{\{fbase\}_slice_0092.nii.gz}, and 
  \file{\{fbase\}_slice_0092_mask.nii.gz}.  
  These data files were obtained from the BrainWeb repository of the McConnell
  Brain Imaging Center at the Montreal Neurological Institute.
  BrainWeb anatomical models uses MRI slices of dimension 181x217 pixels.
  The datasets included in the package for demonstration correspond to a
  T1 BrainWeb image for slice number 92, with 3\% noise and 0\% intensity
  non-uniformity. }
  \item{swap}{ logical variable (default = \samp{FALSE}) for choosing the
  right/left data display convention consistent with FSLVIEW }  
}

\details{
  The FSL tools may be used to prepare the MRI data and the mask required
  as data input. The package \pkg{oro.nifti} is used for reading
  gzipped NIFTI files.
}


\value{
  a list containing
  \item{fbase}{ dataset prefix of the dataset used in the analysis }
  \item{niislice}{ slice data at all timepoints }
  \item{mask}{ slice mask }
  \item{nrow}{ number of rows }
  \item{ncol}{ number of columns }
  \item{swap}{ relative orientation used in the data setup }
}

\references{
Brandon Whitcher, Volker Schmid and Andrew Thornton, Package \bold{oro.nifti}:
 Rigorous - NIfTI Input / Output, 2010.

FSL/FEAT Analysis tool, FMRIB Software Library (FSL)
 (\url{www.fmrib.ox.ac.uk/fsl})
}

\author{A. Ferreira da Silva, Universidade Nova de Lisboa,
Faculdade de Ciencias e Tecnologia,\cr \email{afs@fct.unl.pt}.
}

\seealso{
  \code{\link{dpmixsim}}
}

\examples{
\dontrun{
slicedata <- readsliceimg(fbase="t1_pn3_rf0", swap=FALSE)
print(str(slicedata))
}}
\keyword{IO}

