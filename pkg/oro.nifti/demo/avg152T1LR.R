mniLR <- readNIfTI(file.path(system.file("nifti", package="oro.nifti"),
                             "avg152T1_LR_nifti"))
image(mniLR)
orthographic(mniLR)
