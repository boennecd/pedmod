library(rmarkdown)
library(reshape2)
(function(){
  setwd("vignettes")

  render("pedigree_partitioning.Rmd", output_format = "html_document")
  render("pedmod.Rmd", output_format = "html_document")
  render("pedmod.Rmd", output_format = github_document(
    pandoc_args = "--webtex=https://render.githubusercontent.com/render/math?math="),
    output_file = "README.md")
  file.copy("README.md", file.path("..", "README.md"), overwrite = TRUE)
  unlink("README.md")
  unlink("README.html")

  fs <- list.files(file.path("man", "figures"))
  for(f in fs)
    file.copy(file.path("man", "figures", f),
              file.path("..", "man", "figures", f), overwrite = TRUE)
  unlink("man", recursive = TRUE)
})()
