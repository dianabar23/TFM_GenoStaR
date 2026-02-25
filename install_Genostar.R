#To install devtools, previously install Rtools4.5

#Install GenoStaR from GitHub:
Sys.unsetenv("GITHUB_PAT")
devtools::install_github("GenoStaR-Genomics-Tools/GenoStaR", auth_token = NULL)

library("GenoStaR")
