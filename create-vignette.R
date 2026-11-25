(function(){
  setwd("vignettes")

  rmarkdown::render(input = "pedigree_partitioning.Rmd",
                    output_format = "html_document")
  rmarkdown::render(input = "pedmod.Rmd", output_format = "html_document")
  rmarkdown::render(
    input = "pedmod.Rmd",
    output_format = rmarkdown::github_document(
      math_method = "default", html_preview = FALSE),
    output_file = "README.md")
  file.copy("README.md", file.path("..", "README.md"), overwrite = TRUE)
  unlink("README.md")

  fs <- list.files(file.path("man", "figures"))
  for(f in fs)
    file.copy(file.path("man", "figures", f),
              file.path("..", "man", "figures", f), overwrite = TRUE)
  unlink("man", recursive = TRUE)
})()
