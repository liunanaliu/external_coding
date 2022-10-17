#20221017na for fcs 
#Load the packages
library(flowCore)
library(flowAI)
library(ggcyto)

#How to get help
??flowCore

#Load a single file
#BiocManager::install("flowWorkspaceData")
dataDir <- system.file("extdata",package = "flowWorkspaceData")

library(flowCore)
file.name <- system.file("extdata","0877408774.B08",
                         package="flowCore")
file.name
X=read.FCS(file.name, transformation=FALSE)
X
Y=read.FCS(paste(dataDir,'/CytoTrol_CytoTrol_1.fcs',sep=''),
           transformation=FALSE)
Y
exprs(X[1:10,])
parameters(X)$name
keyword(X)
names(X)
exprs(X)
each_col(X, median)

#load batch .fcs files
fcs.dir <-system.file("extdata","compdata","data",
                      package="flowCore")

fs <- read.flowSet(path = fcs.dir)
fs
sampleNames(fs)
exprs(fs[[1]][1:12,])
head(exprs(fs[[1]]))

#Compensation
fcsfiles <- list.files(pattern ="CytoTrol",
                       system.file("extdata",
                                   package="flowWorkspaceData"),
                       full = TRUE)
fcsfiles
fs <- read.flowSet(fcsfiles)
x <- fs[[1]]
x
#get spillover matrix
comp_list <- spillover(x)
comp_list

comp <- comp_list[[1]]
x_comp <- compensate(x, comp)##for single fcs
comp <- fsApply(fs, 
                function(x)spillover(x)[[1]], 
                simplify=FALSE) ###for multiple fcs
fs_comp <- compensate(fs, comp)

###do not run!!similar function, different code
spillover(fcsfile)
fcsfile_comp <-compensate(fcsfile, spillover(fcsfile)$SPILL)
fcsfile_comp


###check compensation effects
library(gridExtra)
transList <- estimateLogicle(x,c("V450-A","V545-A"))
p1 <- autoplot(transform(x, transList),"V450-A", "V545-A") +
  ggtitle("Before")
p2 <- autoplot(transform(x_comp,transList), "V450-A", "V545-A")+
  ggtitle("After")
grid.arrange(as.ggplot(p1), as.ggplot(p2),ncol = 2)


#Cleaning
#for single file
x_comp_clean <- flow_auto_qc(x_comp)
x_comp_clean
keyword(x_comp_clean) <- keyword(x)
x_comp_clean[1]

#for a folder with multiplle files
fs_comp_clean <- flow_auto_qc(fs_comp)
fs_comp_clean

#Transformation
x_comp_clean[1]
trans <- estimateLogicle(x_comp_clean, 
                         colnames(x_comp_clean[,5:11]))
x_comp_clean_trans <- transform(x_comp_clean, trans)
x_comp_clean_trans
autoplot(x_comp_clean_trans)

ggcyto(x_comp_clean_trans,aes(x="FSC-A", y="SSC-A"))+
  geom_hex(bins=256)

trans <- estimateLogicle(fs_comp_clean[[1]], 
                         colnames(fs_comp_clean[[1]][,5:11]))
fs_comp_clean_trans <- transform(fs_comp_clean, trans)
gs <- GatingSet(fs_comp_clean_trans)
#work for multiple files

#cell gate based on FSC & CCS to exclude debris
rg1 <- rectangleGate("FSC-A"=c(70000,Inf),
                     filterId = "NoneDebris")
gs_pop_add(gs, rg1, parent="root")
recompute(gs)
gs_get_pop_paths(gs)
ggcyto(fs_comp_clean_trans[[1]],aes(x="FSC-A", y="SSC-A"))+
  geom_hex(bins=256)+
  geom_gate(gs_pop_get_gate(gs, "NoneDebris"))
gs_pop_get_stats(gs)

#Singlet gating
ggcyto(fs_comp_clean_trans[[1]], aes(x = "FSC-H", y = 'FSC-W'))+ 
  geom_hex(bins = 256)
rg2 <- rectangleGate("FSC-H"=c(50000, 100000),"FSC-W"=c(30000, 100000))
gs_pop_add(gs, rg2, parent = "NoneDebris", name = "singlets")
gs_get_pop_paths(gs)
recompute(gs)
ggcyto(fs_comp_clean_trans, aes(x = "FSC-H", y = 'FSC-W'))+ 
  geom_hex(bins = 256)+ geom_gate(gs_pop_get_gate(gs, "singlets"))

#exploring the gatingset
plot(gs)
gs_pop_get_stats(gs)
gs_pop_get_stats(gs, "NoneDebris", "percent")


suppressMessages(gs <-load_gs(list.files(dataDir, 
                                         pattern = "gs_manual",full = TRUE)))

# get stats all nodes
dt <- gs_pop_get_stats(gs) #default is"count"

nodes <- c("CD4","CD8")
gs_pop_get_stats(gs, nodes,"percent")
plot(gs)
autoplot(gs[[1]])

##automatic!!!
autoplot(fs_comp_clean_trans[[1]])
autoplot(fs_comp_clean_trans[[2]])

auto_gs <- GatingSet(fs_comp_clean_trans)
fs_data<- gs_pop_get_data(auto_gs)
noneDebris_gate<- fsApply(fs_data, 
                          function(fr) openCyto:::.flowClust.2d(fr, 
                                                                channels= c("FSC-A","SSC-A")))

gs_pop_add(auto_gs, 
           noneDebris_gate, 
           parent = "root", 
           name="noneDebris_gate")
recompute(auto_gs)
autoplot(auto_gs, x="FSC-A", y="SSC-A", "noneDebris_gate", bins=256)

#Singlet gate
fs_data <- gs_pop_get_data(auto_gs, "noneDebris_gate") #get parent data
singlet_gate <- fsApply(fs_data, 
                        function(fr) openCyto:::.singletGate(fr, channels =c("FSC-A", "FSC-H")))
gs_pop_add(auto_gs, 
           singlet_gate, 
           parent = "noneDebris_gate", 
           name = "singlets")
recompute(auto_gs)
autoplot(auto_gs, x = 'FSC-A', y = 'FSC-H',
         "singlets", bins = 256)

#Quad gate
#error ggcyto(gs,aes(x='FSC-A'))+
  geom_density()+
  geom_gate("singlets")
fs_data <- gs_pop_get_data(auto_gs, "singlets") #get parent data
BGquad_gate <- fsApply(fs_data, 
                       function(fr) openCyto:::.quadGate.tmix(fr, #gFunc="mindensity", min=c(0,0), 
                                                             channels =c('B710-A', 'R780-A'),
                                                             K=3,usePrior="no"))
gs_pop_add(auto_gs, BGquad_gate, 
           parent = "singlets", 
           names = c("1", "2", "3", "4"))
recompute(auto_gs)
gs_get_pop_paths(auto_gs[[1]])
plot(auto_gs)
autoplot(fs_comp_clean_trans, x = 'B710-A', y = 'R780-A')
autoplot(auto_gs, x = 'B710-A', y = 'R780-A', 
         gs_get_pop_paths(auto_gs)[4:7], bins = 256)

#not ideal, so fix plot
colnames(fs_comp_clean[[1]])
fs_comp_clean[[1]]@parameters@data

p<-ggcyto(auto_gs[1:2],
          aes(x = 'B710-A', y = 'R780-A'), 
          subset="singlets", arrange = FALSE)
p<- p + geom_hex(bins=256)
p<- p + geom_gate(gs_get_pop_paths(auto_gs)[4:7]) 
p<- p + geom_stats(gs_get_pop_paths(auto_gs)[4:7])
p<- p + theme(strip.text = element_text(size = 7))
myPars <- ggcyto_par_set(limits = list(y = c(0,5), x = c(0,5)))
p<- p  + myPars
p

#Removing stuff
gs_pop_remove(auto_gs, "singlets")

#statistics
gs_pop_get_stats(auto_gs)
gs_pop_get_stats(auto_gs, "noneDebris_gate", "percent")
gs_pop_get_stats(auto_gs, "noneDebris_gate", type = pop.MFI)

pop.quantiles <- function(fr){
  chnls <- colnames(fr)
  res <- matrixStats::colQuantiles(exprs(fr), probs = 0.75)
  names(res) <- chnls
  res
}
gs_pop_get_stats(auto_gs, gs_get_pop_paths(auto_gs), 
                 type = pop.quantiles)

pop.mean <- function(fr){
  chnls <- colnames(fr)
  res <- colMeans(exprs(fr))
  names(res) <- chnls
  res
}
gs_pop_get_stats(auto_gs, gs_get_pop_paths(auto_gs), 
                 type = pop.mean)




