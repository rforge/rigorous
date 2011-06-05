mniRL <- readNIfTI(file.path(system.file("nifti", package="oro.nifti"),
                             "mniRL_nifti"))
image(mniRL)
orthographic(mniRL)
