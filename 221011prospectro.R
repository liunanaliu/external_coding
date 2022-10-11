#221004 substitute for peaks
#https://cran.r-project.org/web/packages/prospectr/vignettes/prospectr.html
##remotes::install_github('l-ramirez-lopez/prospectr')
library(prospectr)

##remotes::install_github('philipp-baumann/simplerspec')
##library(simplerspec)
library(tidyverse)


##setwd("B:/ORG/Forschung/AN3/Foo_WanLing_WF/07_Chol-Liposome/FTIR/Lipo_quantification_klea/DPPC_spectra/CSV")
##lf <- dir("B:/ORG/Forschung/AN3/Foo_WanLing_WF/07_Chol-Liposome/FTIR/Lipo_quantification_klea/DPPC_spectra/CSV",
          full.names = T)

###########
#221011 batch processing
#do not run!!
names <- dir()
data <- NULL
data <- list(classnames=NULL, data=NULL, labels=NULL, x=NULL)
for (i in 1:length(names)){
  the_ith <- read.csv(names [i])[,-1]
  data$data <- rbind(data$data, t(the_ith))
  data$labels <- c(data$labels, array(i, NCOL(the_ith)))
  data$classnames <- c(data$classnames, names [i])
  print(i)
}

dim(data$data)


#######modified by na
setwd('/media/l33/2708-A9DC/CSV')
dir <- file.path(getwd())
dir
files <- list.files(pattern = '.csv',dir,recursive = T)
files
names(files)
library(readr)
##sapply(strsplit(files,'_'),function(x)x[c(1:3)])
runs <- sapply(strsplit(files,'-'),function(x)x[1])
names(files) <- runs
files
myfiles <- lapply(files,read.csv)
myfiles
################
#do not run!!
install.packages("chromatographR")
library(chromatographR)
read_chroms(paths = dir,
            format_in = 'chemstation_csv')
#error


#####################
myfiles$DPPC_0.5_1_2021$Intensity
values <- do.call(rbind,lapply(myfiles,function(x) x$Intensity))
values
values <-  `colnames<-`(values, paste0("lambda_", myfiles[[1]]$Wavenumber))
df <- data.frame(file_number=rownames(values),values)
dim(values)
dim(df)

bak <- read.csv('/home/l33/Downloads/All_bak_mean.csv')
library(tidyverse)
##str_split(rownames(values),'_',simplify = T)[,1]
is.data.frame(values)

group <- matrix(unlist(str_split(df$file_number,'_',simplify = T)[,1]))
df$group <- group
df$con <- matrix(unlist(str_split(df$file_number,'_',simplify = T)[,2]))
table(df$group,df$con)
table(bak$group,bak$conc)
###
library(tidyr)
bak1 <- gather(df, 'Wavelength', 'Intensity', 
               -c('file_number','group','con'))
table(bak1$group,bak1$con)
head(bak1)
head(bak)






##spc_list <- read_opus_univ(fnames = lf,
                           extract = c('spc'))

##lf <- read.csv('/home/l33/Downloads/DPPC_0.5_1_2021-10-27T18-48-18.csv')
##lf[,1]
##lf[1:900,1]
##str(lf)


##data("NIRsoil")
##str(NIRsoil$spc)
##colnames(NIRsoil$spc)

lf_subtracted <- as.data.frame(t(lf))
colnames(lf_subtracted) <- lf_subtracted[1,]
##the same frame for df, 
##wavelength as colname, 1 column for intensity 
lf_subtracted <- lf_subtracted[2,]
wav <- as.numeric(colnames(lf_subtracted))
#get the wavelength
lf_rm_bs <- baseline(lf_subtracted,
                          wav)
lf_rm_bs <- as.data.frame(lf_rm_bs)
colnames(lf_rm_bs) <- 'intensity'
lf_rm_bs <- as.data.frame(t(lf_rm_bs))

plot(lf[1:900,1],lf_subtracted,
     t='l')
plot(lf[1:900,1],lf_rm_bs,
     t='l')
######################################
##do not run!!
lf_snv <- standardNormalVariate(lf_subtracted)
plot(lf[1:900,1],lf_snv,
     t='l')
########################################
lf_norm <- blockNorm(X=lf_rm_bs,
                     targetnorm = 1)$Xscaled
plot(wav,lf_norm,
     t='l')
#choose 1 of the two
lf_norm <- blockNorm(X=lf_rm_bs,
                     targetnorm = 1)$Xscaled
plot(wav,lf_norm,
     t='l')


