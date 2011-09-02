qaqc <- function(subject, visit, exam, setup, pixelData=TRUE,
                 domain="/usr/cic/apps/bin/domain", verbose=FALSE, ...) {
  ## Clean out 
  system(paste("rm -rf", file.path(setup$dicom, "*")))
  ## Copy over setup$mailText
  mailText <- setup$mailText

  for (i in 1:length(subject)) {
    if (verbose) {
      cat("## SUBJECT:", subject[i], " VISIT:", visit[i], " EXAM:",
          exam[i], fill=TRUE)
    }
    con <- file(setup$logFile, "a", encoding="UTF-8")
    cat("## SUBJECT", subject[i], file=con, fill=TRUE)
    cat("## VISIT", visit[i], file=con, fill=TRUE)
    cat("## EXAM", exam[i], file=con, fill=TRUE)
    cat("#### DOWNLOAD ####", file=con, fill=TRUE)
    close(con)
    ##
    ## Execute Domain download
    ##
    cmd <- paste(domain, "-download",
                 "-study_id", setup$study,
                 "-subject_id", subject[i],
                 "-visit_id", visit[i],
                 "-exam_index", exam[i],
                 "-folder PI", "-o", setup$dicom, ">>", setup$logFile)
    system(cmd)
    ##
    ## Collect DICOM header information
    ##
    dcm <- dicomSeparate(setup$dicom, verbose=verbose, pixelData=pixelData, ...)
    eval(multipleHeader(dcm, "StudyInstanceUID"))
    eval(multipleHeader(dcm, "ClinicalTrialSiteID"))
    eval(multipleHeader(dcm, "ClinicalTrialProtocolID"))
    eval(multipleHeader(dcm, "ClinicalTrialTimePointID"))
    eval(multipleHeader(dcm, "InstitutionName"))
    eval(multipleHeader(dcm, "PatientID"))
    eval(multipleHeader(dcm, "StudyDate"))
    ##
    ## Initialize HTML output
    ##
    if (! is.null(setup$shared.area)) {
      html.path <- file.path("", "scratch", "oncology", "Clinical_Imaging",
                             "Imaging_Initiative", "QC", "CIC_HTML_Outputs",
                             clinicalTrialProtocolID,
                             clinicalTrialSiteID, patientID,
                             clinicalTrialTimePointID, studyDate)
    } else {
      html.path <- file.path(sub("dicom", "html", setup$dicom),
                             clinicalTrialProtocolID,
                             clinicalTrialSiteID, patientID,
                             clinicalTrialTimePointID, studyDate)
    }
    dir.create(html.path, showWarnings=FALSE, recursive=TRUE)
    ## Clean out HTML
    system(paste("rm -rf", file.path(html.path, "*")))
    if (! is.null(setup$css)) {
      system(paste("cp ", setup$css, " ", html.path, sep=""))
    }
    htmlfile <- openPage(paste(studyInstanceUID, "html", sep="."),
                         html.path, link.css=sub(".*/", "", setup$css),
                         title=paste("Report for", patientID, "on",
                           str2date(studyDate)))
    if (verbose) {
      cat(paste("##", html.path), fill=TRUE)
    }
    cat("##", html.path,
        file=(con <- file(setup$logFile, "a", encoding="UTF-8")), fill=TRUE)
    close(con)
    ## Save all DICOM header information to a CSV file
    csv <- dicomTable(dcm$hdr)
    write.csv(csv, file.path(html.path,
                             paste(studyInstanceUID, "csv", sep=".")))
    rm(csv)
    ##
    if (! is.null(setup$shared.area)) {
      html.path <- sub("/scratch/oncology", setup$shared.area, html.path)
    }
    mailText <-
      paste(mailText,
            file.path(html.path, paste(studyInstanceUID, "html", sep=".")),
            sep=", ")
    hwrite(paste("Report for Subject", patientID, "on", str2date(studyDate)),
           htmlfile, heading=1)
    hwrite(paste("<strong>Study:</strong>", setup$study,
                 "<br>\n",
                 "<strong>Study UID:</strong>", studyInstanceUID,
                 "<br>\n",
                 "<strong>Institution Name:</strong>", institutionName,
                 "<br>\n",
                 "<strong>Visit:</strong>", visit[i],
                 "<br>\n",
                 "<strong>Exam:</strong>", exam[i],
                 "<br>\n"), htmlfile)
    ##
    ## Convert DICOM to NIfTI
    ##
    cat("##### CONVERT #####",
        file=(con <- file(setup$logFile, "a", encoding="UTF-8")), fill=TRUE)
    close(con)
    seriesInstanceUID <- extractHeader(dcm$hdr, "SeriesInstanceUID", FALSE)
    for (uid in unique(seriesInstanceUID)) {
      if (verbose) {
        cat("##  ", uid, fill=TRUE)
      }
      cat("##    ", uid, 
          file=(con <- file(setup$logFile, "a", encoding="UTF-8")), fill=TRUE)
      close(con)
      index <- as.vector(unlist(lapply(dcm$hdr, function(x) uid %in% x$value)))
      Z <- sum(index)
      uid.dcm <- list(hdr=dcm$hdr[index], img=dcm$img[index])
      hwrite(hmakeTag("hr"), htmlfile)
      hwrite(paste("<strong>Series UID:</strong>", uid, "<br>\n"), htmlfile)
      ## EXIT = Multiple SeriesDescription fields
      eval(multipleHeader(uid.dcm, "SeriesDescription", htmlfile))
      ##
      hwrite(paste("<strong>Series Description:</strong>", seriesDescription,
                   "<br>\n"), htmlfile)
      ## EXIT = Multiple Modality fields
      eval(multipleHeader(uid.dcm, "Modality", htmlfile))
      ##
      hwrite(paste("<strong>Modality:</strong>", modality,
                   "<br>\n"), htmlfile)
      ## EXIT = Multiple ImageType fields
      eval(multipleHeader(uid.dcm, "ImageType", htmlfile))
      ## EXIT = LOCALIZER or SCOUT or SECONDARY or DERIVED acquisitions
      eval(compareHeader(uid.dcm, "ImageType",
                         c("loc","localizer","scout","secondary","derived"),
                         htmlfile))
      eval(compareHeader(uid.dcm, "SeriesDescription",
                         c("loc","localizer","scout","secondary","derived"),
                         htmlfile))
      ## EXIT = {CR, DX, PR} modality
      eval(compareHeader(uid.dcm, "Modality", c("cr","dx","pr"), htmlfile))
      ## EXIT = pixel dimensions must match across all slices
      pixelSpacing <- try(header2matrix(extractHeader(uid.dcm$hdr, "PixelSpacing", FALSE), 2))
      error <- geterrmessage()
      if (error == "error") {
        hwrite("Series ignored because of error in PixelSpacing.",
               htmlfile, heading=3)
        next
      }
      if (length(unique(pixelSpacing[,1])) > 1 &&
           length(unique(pixelSpacing[,2])) > 1) {
        hwrite("Pixel dimensions are NOT consistent.", htmlfile, heading=3)
        next
      }
      ##
      ## EXIT = ImagePositionPatient doesn't make sense
      imagePositionPatient <-
        header2matrix(extractHeader(uid.dcm$hdr, "ImagePositionPatient",
                                    FALSE), 3)
      dimnames(imagePositionPatient) <- list(1:Z, c("X","Y","Z"))
      movingDimensions <- apply(imagePositionPatient, 2,
                                function(j) any(diff(j) != 0))
      if (sum(movingDimensions) != 1) {
        hwrite("ImagePositionPatient is moving in multiple dimensions.",
               htmlfile, heading=3)
        ## hwrite(imagePositionPatient, htmlfile, digits=1,
        ##        caption="ImagePositionPatient for all slices")
      }
      rm(imagePositionPatient)
      rescale <- ifelse(modality == "CT", TRUE, FALSE)
      uid.nifti <- dicom2nifti(uid.dcm, datatype=4, rescale=rescale,
                               descrip=c("ProtocolName",
                                 "SeriesDescription", "ImageType"),
                               pixelData=pixelData)
      writeNIfTI(uid.nifti, file.path(html.path, uid))
      ##
      ## Produce diagnostic plots of imaging data
      ##
      if (Z > 1) {
        ## Orthographic figure
        graphOrthographic(uid.nifti, uid, modality, html.path, htmlfile)
        ## Slice "thickness" figure
        graphSliceThickness(uid.nifti, uid, imagePositionPatient,
                            movingDimensions, html.path, htmlfile)
      } else {
        ## 2D slice figure
        graphSingleSlice(uid.nifti, uid, modality, html.path, htmlfile)
      }
      if (! is.null(setup$compare)) {
        ## Verify DICOM header information against study template
        output.table <- verifyHeader(uid.dcm$hdr, setup$compare)
        hwrite(output.table, htmlfile, center=TRUE, br=TRUE,
               row.bgcolor="#ffffaa", col.bgcolor="#ffffaa",
               row.style=list("font-weight:bold"), 
               col.style=list("Field"="font-weight:bold",
                 "Result"=ifelse(output.table[,"Result"] == "Issue", "color:#cc0000", "color:#33cc00")), 
               table.frame="void")
      }
      rm(uid.dcm, uid.nifti)
    }
    closePage(htmlfile)
    rm(dcm)
    system(paste("rm -rf", file.path(setup$dicom, "*")))
  }
  if (! is.null(setup$emailList)) {
    cmd <- paste("echo", mailText,
                 "| tr '/' '\\\\' | tr ',' '\012' | mail -s \"",
                 setup$study, "HTML update\"", setup$emailList)
    system(cmd)
  }
  invisible()
}
