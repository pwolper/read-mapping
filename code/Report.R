
print(mapfile)
rmarkdown::render("read-vis.Rmd",params = list(output,inputGenome,mapfile))



