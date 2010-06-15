mniRL <- readNIfTI(file.path(system.file("nifti", package="oro.nifti"),
                             "/avg152T1_RL_nifti"))
image(mniRL)
orthographic(mniRL)
