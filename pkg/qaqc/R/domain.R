parseDomain <- function(dom.list, character.split="/", index=1) {
  out <- NULL
  if (length(character.split) == 1) {
    out <- unlist(lapply(strsplit(dom.list, character.split),
                         function(x) x[index]))
  }
  if (length(character.split) == 2) {
    temp <- unlist(lapply(strsplit(dom.list, character.split[1]),
                          function(x) x[index[1]]))
    out <- unlist(lapply(strsplit(temp, character.split[2]),
                  function(x) x[index[2]]))
  }
  if (is.null(out)) {
    warning("parseDomain produced a NULL result")
  }
  return(out)
}
