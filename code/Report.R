
#needed as link to read-vis.Rmd

para = list(out = out,mapfile = mapfile,inputGen = inputGen)

if (!require("rmarkdown")) {
  install.packages("rmarkdown")
  library(rmarkdown)
}
if (!require("knitr")) {
  install.packages("knitr")
  library(knitr)
}


render("read-vis.Rmd",params = para)



