using Weave
#pdf
#cd("/home/moralesmendozar/Dropbox/03_UPenn/classes/2019_fall/00_714/partJesus/PS01_jfv")
#weave(Pkg.dir("714_PS01_01.jmd"), informat="markdown", out_path = :pwd, doctype = "md2pdf")
weave("714_PS01_01.jmd", informat="markdown", out_path = :pwd, doctype = "md2pdf")

#![Policy Function for the Stochastic Grid]("Plots/003_PolicyFunction_d_stochastic_20191201.png"){ width=50% }
#![Policy Function for the Stochastic Grid]("Plots/003_PolicyFunction_d_stochastic_20191201.png"){ width=50% }
#<img src="Plots/003_PolicyFunction_d_stochastic_20191201.png" width="200">


# ```{r, fig.cap="Policy Function for the Stochastic Grid", out.width = '50%'}
# knitr::include_graphics('Plots/003_PolicyFunction_d_stochastic_20191201.png')
# ```
