library("qaqc")

##First read in the arguments listed at the command line
args=(commandArgs(TRUE))

##args is now a list of character vectors
## First check to see if arguments are passed.
## Then cycle through each element of the list and evaluate the expressions.
if(length(args)==0){
    stop("No arguments supplied.Please provide '--args study=[studyID] date2process=[DD-MM-YYYY] production=[TRUE|FALSE]'")
    ##supply default values
    
}else{
    for(i in 1:length(args)){
         eval(parse(text=args[[i]]))
      

    }
       print(study)
       today = date2process
       #make date2process a date object
       
       print(today)
       ## today = "2010-03-22"
    	 toddot <- sub("-", ".", today) #replace first '-'
    	 toddot <- sub("-", ".", toddot)#replace second '-' 
       ## toddot = "22.03.2010"
       if (production){
       	    print("Running on Production")
       	    server.share <- file.path("", "", "hlwsfshr10", "sa200400214")#is this used? how is it passed to qaqc
	    gsk.com <- "bjw34032@gsk.com,cf984373@gsk.com"
	    mail.text <- paste(today, ", ", sep="")
	    data <- file.path("", "cic_data", "Oncology", study)
	    dicom <- file.path("", "scratch", study, "dicom")
	    nifti <- file.path(data, "nifti")
	    scripts <- file.path(data, "scripts")
	    log <- file.path(data, "logs")
	    logfile <- file.path(log, paste(today, ".log", sep=""))
	    domfile <- file.path(log, paste(today, "_domain.txt", sep=""))
	    source(file.path(data, "scripts", paste(tolower(study), "qaqc.R", sep="_")))
       	}
        else{ print("Running on Development")
        
            server.share <- file.path("", "", "hlwsfshr10", "sa200400214")#is this used? how is it passed to qaqc
            gsk.com <- "bjw34032@gsk.com,cf984373@gsk.com"
            mail.text <- paste(today, ", ", sep="")
            data <- file.path("", "cic_data", "Oncology", study, "test")
            dicom <- file.path("", "scratch", study, "dicom")
            nifti <- file.path(data, "nifti")
            scripts <- file.path(data, "scripts")
            log <- file.path(data, "logs")
            logfile <- file.path(log, paste(today, ".log", sep=""))
            domfile <- file.path(log, paste(today, "_domain.txt", sep="")) 
            source(file.path(data, "scripts", paste(tolower(study), "qaqc.R", sep="_")))
        }
	 
}


known.values <- file.path(scripts, paste(tolower(study), "specs.csv", sep="_"))

cat("#### SEARCH ####", file=(con <- file(logfile, "w", encoding="UTF-8")),
    fill=TRUE)
close(con)
domain <- file.path("", "usr", "cic", "apps", "bin", "domain")
system(paste(domain, "-search -study_id", study, "-folder 'PI' >", domfile))
domfile.names <- c("Date", "Time", "File.Name", "Document.Size",
                   "Document.ID", "Mime.Type")
domList <- strsplit(grep(toddot, readLines(domfile), value=TRUE), "[ ]+")
domStrings <- unlist(lapply(domList, function(z) z[4]))
subjectID <- parseDomain(domStrings, index=3)
if (length(subjectID) > 0) {
  
  visitID <- parseDomain(domStrings, c("/","_"), c(4,1))
  examID <- parseDomain(domStrings, c("/","_"), c(4,3))
  #Call to QAQC may need to be changed to reflect new param list
  qaqc(subjectID, visitID, examID, study, logfile, domfile, dicom,
       nifti, compare=known.values, pixelData=FALSE, email=gsk.com,
       css=file.path(scripts, paste(tolower(study), "css", sep=".")),
       verbose=TRUE, debug=FALSE)
} else {
  cat("##### No data to download from DOMAIN #####",
      file=(con <- file(logfile, "a", encoding="UTF-8")), fill=TRUE)
  close(con)
}

subject=subjectID
visit=visitID
exam=examID
study=study
logfile=logfile
domfile=domfile
domdir=dicom
niftidir=nifti
compare=known.values
pixelData=FALSE
email=gsk.com
css=file.path(scripts, paste(tolower(study), "css", sep="."))
verbose=TRUE
debug=FALSE

rm(subject, visit, domdir, niftidir, compare, email, css, verbose, debug)
