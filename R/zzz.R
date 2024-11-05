.onLoad <- function(libname, pkgname) {
  source(system.file("scripts", "KnockoffScreen_AL_vMar16.R", package = pkgname), local = FALSE)
  source(system.file("scripts", "ESGWAS.R", package = pkgname), local = FALSE)
  source(system.file("scripts", "GhostKnockoff.R", package = pkgname), local = FALSE)
}