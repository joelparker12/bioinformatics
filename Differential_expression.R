if(!requireNamespace("BiocManager", quietly = TRUE))
  install.packages('BiocManager')  
BiocManager::install("maSigPro")

library(maSigPro)

load("single.Rdata")

ss.edesign

head(ss.DATA)

ss.example <- maSigPro(ss.DATA, ss.edesign)

head(ss.example$p.vector)

### Potato Expirenment. 
data("data.abiotic")
data("edesign.abiotic")
head(edesign.abiotic)


#### Create regression matrix for the full regression model. 
design <- make.design.matrix(edesign.abiotic, degree = 2) ### check to see why degree=2
##This example has three time points, so we can consider up to a quadratic regression model (degree = 2). Larger number of time points would potentially allow a higher polynomial degree.

## NOw we find singnificant genes. 
fit <- p.vector(data.abiotic, design, Q=.05, MT.adjust = "BH", min.obs = 20)

###COmmands= lecture 20 


#NOw find significant differences
tstep <- T.fit(fit, step.method = "backward")


# Obtain lists of significant genes. 
sigs <- get.siggenes(tstep, rsq = .6, vars = "groups")


suma2Venn(sigs$summary[, c(1:4)])

#### How many significant differentially expressed genes?
sigs$sig.genes$ColdvsControl$g

see.genes(sigs$sig.genes$ColdvsControl,show.fit = T,dis = design$dis)


# Looking at NB data
data(NBdata)
data("NBdesign")

d <- make.design.matrix(NBdesign)

library(MASS)
NBp <- p.vector(NBdata, d, counts = TRUE)
NBt <- T.fit(NBp)

get <- get.siggenes(NBt, vars = "all")
get$summary

see.genes(get$sig.genes, k=4)

