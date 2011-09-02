verifyHeader <- function(hdr, reference, nr=NULL, nc=4,
                         out.names=c("Field", "Value", "Reference Value",
                           "Result")) {
  reftable <- read.csv(reference, header=FALSE, colClasses="character")
  fields <- reftable[,1]
  nr <- ifelse(is.null(nr), nrow(reftable), nr)
  out.table <- matrix("", ifelse(is.null(nr), length(fields), nr), nc)
  dimnames(out.table) <- list(NULL, out.names)
  out.table[,"Field"] <- fields
  for (m in 1:nr) {
    value <- extractHeader(hdr, fields[m], FALSE)[1]
    if (length(grep("Date", fields[m], ignore.case=TRUE)) > 0) {
      out.table[m,"Value"] <- str2date(value)
    } else {
      out.table[m,"Value"] <- value
    }
    refindex <- reftable[,1] %in% fields[m]
    refvalue <- reftable[refindex,2]
    reftype <- reftable[refindex,3]
    out.table[m,"Reference Value"] <- refvalue
    ## Check value against "reference value"
    if (! is.na(refvalue)) {
      if (reftype == "character") {
        refvalue <- unlist(strsplit(refvalue, "|", fixed=TRUE))
        result <- any(value == refvalue)
      } else {
        result <- eval(parse(text=paste(value, refvalue)))
      }
      out.table[m,"Result"] <- ifelse(result, "No Issue", "Issue")
    }
  }
  out.table[is.na(out.table)] <- ""
  as.data.frame(out.table, stringsAsFactors=FALSE)
}

multipleHeader <- function(dcm, dicomHeader, htmlfile=NULL, numeric=FALSE) {
  out <- na.omit(extractHeader(dcm$hdr, dicomHeader, numeric=numeric))
  if (length(unique(out)) > 1) {
    errorText <- paste("Multiple", dicomHeader, "strings found.")
    if (is.null(htmlfile)) {
      stop(errorText)
    } else {
      hwrite(errorText, htmlfile, heading=3)
    }
    expression(next)
  } else {
    eval(substitute(expression(x <- y),
                    list(x=paste(tolower(substring(dicomHeader,1,1)),
                           substring(dicomHeader,2,nchar(dicomHeader)),
                           sep=""),
                         y=unique(out))))
  }
}

compareHeader <- function(dcm, dicomHeader, strings, htmlfile,
                          numeric=FALSE) {
  out <- na.omit(extractHeader(dcm$hdr, dicomHeader, numeric=numeric))
  for (i in 1:length(strings)) {
    if (any(matchHeader(out, strings[i]))) {
      hwrite(paste(toupper(strings[i]), "image(s) found and ignored."),
             htmlfile, heading=3)
      return(expression(next))
    } 
  }
  eval(substitute(expression(x <- y),
                  list(x=paste(tolower(substring(dicomHeader,1,1)),
                         substring(dicomHeader,2,nchar(dicomHeader)),
                         sep=""),
                       y=unique(out))))
}
