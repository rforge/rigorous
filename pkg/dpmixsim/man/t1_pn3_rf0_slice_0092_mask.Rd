\name{t1_pn3_rf0_slice_0092_mask}
\alias{t1_pn3_rf0_slice_0092_mask.Rd}
\docType{data}
\title{Mask file for MRI slice}
\description{
  The \file{t1_pn3_rf0_slice_0092_mask.nii.gz} defines the mask for \file{t1_pn3_rf0_slice_0092.nii.gz}, as used in the examples.
  The mask used here is an all-brain mask; it just removes non-brain regions, as the result of applying a brain extraction tool to the specified dataset.
  Other masks may be defined to select regions of interest (ROIs).
}
\usage{
x <- premask(slicedata, subsamp=TRUE)
}
\format{
  The file \file{t1_pn3_rf0_slice_0092_mask.nii.gz} is in gzipped NIFTI format.
  The R-package \pkg{oro.nifti} is required to read gzipped NIFTI files.
}
\references{
D.~L. Collins, A.~P. Zijdenbos, V.~Kollokian, J.~G. Sled, N.~J. Kabani, C.~J.
  Holmes, A.~C. Evans, Design and construction of a realistic digital brain
  phantom, \emph{IEEE Trans. on Medical Imaging} 17~(3) (1998) 463--468.      

S.M. Smith, et. al., Advances in Functional and Structural MR Image Analysis and
 Implementation as FSL, \emph{NeuroImage}, 23(S1):208-219, 2004.

Brandon Whitcher, Volker Schmid and Andrew Thornton, Package \bold{oro.nifti}:
 Rigorous - NIfTI Input / Output, 2010.
}
\keyword{datasets}
