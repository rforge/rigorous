mniLR <- readNIfTI(file.path(system.file("nifti", package="oro.nifti"),
                             "mniLR_nifti"))
image(mniLR)
orthographic(mniLR)
