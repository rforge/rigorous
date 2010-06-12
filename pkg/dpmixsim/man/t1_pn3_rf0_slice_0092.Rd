\name{t1_pn3_rf0_slice_0092}
\alias{t1_pn3_rf0_slice_0092.Rd}
\docType{data}
\title{Example of a pre-processed MRI slice from the BrainWeb database}
\description{
  The file \file{t1_pn3_rf0_slice_0092.nii.gz} is a pre-processed image of slice \samp{92} with \samp{3\%} noise extracted from the Brainweb database file \file{t1_icbm_normal_1mm_pn3_rf0\[1\].mnc.gz}.
  BrainWeb simulations are based on an anatomical model of normal brain, which can serve as the ground truth for any analysis procedure.
  BrainWeb datasets and are provided by the McConnell Brain Imaging Center at the Montreal Neurological Institute, \url{http://www.bic.mni.mcgill.ca/}, (see Collins et. al. 1998).
}
\usage{
slicedata <- read.sliceimg(fbase="t1_pn3_rf0", swap = TRUE)
}
\format{
  The file \file{t1_pn3_rf0_slice_0092.ni.gz} is in gzipped NIFTI format.
  The R-package \pkg{oro.nifti} is required to read gzipped NIFTI files.
}
\references{
D.L. Collins, et.al., Design and construction of a realistic digital brain
  phantom, \emph{IEEE Trans. on Medical Imaging} 17~(3) (1998) 463-468.      

S.M. Smith, et. al., Advances in Functional and Structural MR Image Analysis and
 Implementation as FSL, \emph{NeuroImage}, 23(S1):208-219, 2004.

Brandon Whitcher, Volker Schmid and Andrew Thornton, Package \bold{oro.nifti}:
 Rigorous - NIfTI Input / Output, 2010.

}
\keyword{datasets}
