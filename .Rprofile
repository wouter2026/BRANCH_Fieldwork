# Ensure rmarkdown can find Pandoc when running from VS Code.
if (Sys.getenv("RSTUDIO_PANDOC") == "") {
  pandoc_dir <- file.path(Sys.getenv("LOCALAPPDATA"), "Pandoc")
  pandoc_exe <- file.path(pandoc_dir, "pandoc.exe")
  if (file.exists(pandoc_exe)) {
    Sys.setenv(RSTUDIO_PANDOC = pandoc_dir)
  }
}
