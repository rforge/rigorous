studyParameters <- function(argList,
                            date2process=NULL,
                            emailList=NULL,
                            study=NULL,
                            inPath=file.path("", "cic_data", "Oncology"),
                            outPath=file.path("", "scratch"),
                            shared.area=NULL) {
  ## Create vector of arguments, keeping on those with a "="
  ## print(argList)
  argList <- grep("=", unlist(strsplit(argList, " ")), value=TRUE)
  if (length(argList) < 1) {
    stop("No arguments supplied.\n  Please provide \'--args parameter=value parameter=value ...\'\n  Note, all character strings must be quoted using \"\".\n  Note, all dates must be provided in DD-MM-YYYY.")
  }
  out <- vector("list")
  for(i in 1:length(argList)){
    eval(parse(text=argList[i]))
  }
  ## Dates
  if (! is.null(date2process)) {
    d2p <- as.Date(date2process, format="%d-%m-%Y")
  } else {
    warning("Argument \'date2process\' has not been supplied, using ",
            format(Sys.time(), "%Y-%m-%d"))
    d2p <- Sys.time()
  }
  out$today <- format(d2p, "%Y-%m-%d")
  out$domDate <- format(d2p, "%d.%m.%Y")
  ## Email
  if (! is.null(emailList)) {
    out$emailList <- emailList
  } else {
    out$emailList <- "bjw34032@gsk.com"
  }
  out$mailText <- paste(out$today, ", ", sep="")
  ## Define path names for input
  if (! is.null(study)) {
    out$study <- study
  } else {
    stop("Argument \'study\' has not been supplied.")
  }
  studyPath <- file.path(inPath, toupper(study))
  scripts <- file.path(studyPath, "scripts")
  out$css <- file.path(scripts, paste(tolower(study), "css", sep="."))
  out$compare <- file.path(scripts, paste(tolower(study), "specs.csv", sep="_"))
  ## Define path names for output
  scratch <- file.path(outPath, toupper(study))
  out$dicom <- file.path(scratch, "dicom")
  #out$nifti <- file.path(scratch, "nifti") # Now = html directory
  out$logFile <- file.path(scratch, "logs", paste(out$today, ".log", sep=""))
  out$domFile <- file.path(scratch, "logs", paste(out$today, ".dom", sep=""))
  out$shared.area <- shared.area
  return(out)
}
