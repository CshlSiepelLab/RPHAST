##' @export
##' @nord
.onAttach <- function(libname, pkgname) {
  if(Sys.info()[[1]] == "Windows") { 
    regDir <- paste(shell("echo %appdata%", intern=TRUE), "\\rtfbs\\", sep="")
  }  else {
    regDir <- paste(system("echo $HOME", intern=TRUE), "/.rphast/", sep="")
  }
  
  Sys.setenv(regDir=regDir)
  
  if(!file.exists(paste(regDir, "registered"))) {
    packageStartupMessage("** If you find RPHAST useful, we would appreciate if you let us know so we can better understand our user base.  Registration is free and can be anonymous.  See ?register.rphast for more details.  Thanks! **")
  }
    
}

##' Register RPHAST
##' @param regName Your Name (Optional)
##' @param regEmail Your Email Address (Optional)
##' @param regInstitution Your Institution (Optinal)
##' @export
register.rphast <- function(regName="", regEmail="", regInstitution="") {

  invisible(read.table(file=paste("http://compgen.bscb.cornell.edu/rphast/register.php?name=", regName, "&email=", regEmail, "&institution=", regInstitution, "&version=", packageDescription("rphast")$Version, "&isValid=TRUE", sep="")))
  regDir <- Sys.getenv("regDir")
  dir.create(regDir, showWarnings=FALSE)
  file.create(paste(regDir, "registered"), showWarnings=FALSE)
}

