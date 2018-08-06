library(Matrix)
library(limma)
library(reshape2)
library(argparse)
library(methods)

parser = argparse::ArgumentParser(description="Script to generate specificity scores from sci-ATAC-seq data.")
parser$add_argument('-PV','--pvals', help='List of sites that are differentially accessibile for each cluster and their p-values saved as an RDS file.')
parser$add_argument('-PR','--props', help='Matrix of sites by clusters with proportion of cells in each cluster that are accessible for each site. Also needs to be an RDS file.')
parser$add_argument('-DP','--depth', help='Matrix specifying median read depth for each cluster.')
parser$add_argument('-ST','--txtout', help='Name of output text file for writing specific results.')
parser$add_argument('-SR','--rdsout', help='Name of output RDS file for writing specific results.')
parser$add_argument('-AT','--alltxtout', help='Name of output text file for writing all specificity scores.')
parser$add_argument('-AR','--allrdsout', help='Name of output RDS file for writing all specificity scores.')
args = parser$parse_args()

makeprobsvec<-function(p){
  phat<-p/sum(p)
  phat[is.na(phat)] = 0
  phat
}

shannon.entropy <- function(p) {
  if (min(p) < 0 || sum(p) <=0)
    return(Inf)
  p.norm<-p[p>0]/sum(p)
  -sum( log2(p.norm)*p.norm)
}

JSdistVec <- function(p, q){
  JSdiv <- shannon.entropy((p + q)/2) - (shannon.entropy(p) + 
                                           shannon.entropy(q)) * 0.5
  JSdiv[is.infinite(JSdiv)] <- 1
  JSdiv[JSdiv < 0] <- 0
  JSdist <- sqrt(JSdiv)
  JSdist
}

specificity_scorer = function(normpropmat){
  marker_gene_specificities <- lapply(1:ncol(normpropmat), function(cell_type_i){
    perfect_specificity <- rep(0.0, ncol(normpropmat))
    perfect_specificity[cell_type_i] <- 1.0
    apply(normpropmat, 1, function(x) { 
      if (sum(x) > 0) 1 - JSdistVec(makeprobsvec(x), perfect_specificity)
      else 0
    })
  })
  return(do.call(cbind, marker_gene_specificities))
}

markerlistmaker = function(markering,daps){
  markerlist = list()
  betamarkerlist = list()
  nonmarkerlist = list()
  markerback = list()
  for(i in 1:ncol(markering)){
    commoners = intersect(rownames(daps[[match(colnames(markering)[i],names(daps))]]),rownames(markering))
    markerlist[[i]] = markering[match(commoners,rownames(markering)),i]
    betamarkerlist[[i]] = daps[[match(colnames(markering)[i],names(daps))]][match(commoners,rownames(daps[[match(colnames(markering)[i],names(daps))]])),1]
    nonmarkerlist[[i]] = markering[-match(commoners,rownames(markering)),i]
    names(markerlist)[i] = colnames(markering)[i]
    names(betamarkerlist)[i] = colnames(markering)[i]
    names(nonmarkerlist)[i] = colnames(markering)[i]
  }
  markerback$markerlist = markerlist
  markerback$betamarkerlist = betamarkerlist
  markerback$nonmarkerlist = nonmarkerlist
  return(markerback)
}

spec_floorer = function(trueo,nullo){
  nullpass = c()
  for(i in 1:length(trueo)){
    nullpass[i] = length(which(nullo >= trueo[i]))
  }
  nullfrac = nullpass/(nullpass + 1:length(nullpass))
  return(trueo[max(which(nullfrac <= 0.1))])
}

results_writer = function(markerlist,specfloor,daps,txtout,rdsout){
  wbmat = matrix(,,5)
  colnames(wbmat) = c("specificity_score","locusID","cluster","subset_cluster","cluster_name")
  for(i in 1:length(markerlist))
  {
    #print(i)
    specsites = markerlist[[i]][which(log10(markerlist[[i]]) > specfloor)]
    if(length(specsites) == 0){
      next
    }
    currmatout = matrix(names(specsites))
    currcluster = names(markerlist)[i]
    currmatout = cbind(currmatout,daps[[match(currcluster,names(daps))]][match(names(specsites),rownames(daps[[match(currcluster,names(daps))]])),])
    currmatout = cbind(currmatout,specsites)
    colnames(currmatout)[c(1,5)] = c("locusID","specificity_score")
    currmatout = currmatout[order(currmatout[,5],decreasing=T),]
    currmatout[,2:5] = signif(currmatout[,2:5],4)
    #print(dim(currmatout)[1])
    clusterer = strsplit2(currcluster,"[.]")
    clusternow = gsub("clusters_","",clusterer[1])
    subnow = gsub("cluster_","",clusterer[2])
    assigns = cbind(rep(clusternow,nrow(currmatout)),rep(subnow,nrow(currmatout)))
    assigns = cbind(assigns,rep(currcluster,nrow(currmatout)))
    currmatforwbmat = cbind(currmatout[,c(5,1)],assigns)
    rownames(currmatforwbmat) = NULL
    colnames(currmatforwbmat) = c("specificity_score","locusID","cluster","subset_cluster","cluster_name")
    wbmat = rbind(wbmat,currmatforwbmat)
  }
  wbmat = wbmat[-1,]
  write.table(wbmat,txtout,row.names=F,quote=F,sep="\t")
  saveRDS(wbmat,rdsout)
}

print("Loading files...")
daps = readRDS(args$pvals)
propmat = readRDS(args$props)
depthmat = read.table(args$depth)

print("Normalizing proportions...")
depthnorm = mean(depthmat[,2])/depthmat[,2]
logdepthnorm = mean(log10(depthmat[,2]))/log10(depthmat[,2])
propmatnormbylogdepth = t(t(propmat)*logdepthnorm)

print("Calculating specificity scores...")
marker_specificities_out = specificity_scorer(propmatnormbylogdepth)
markerdup = marker_specificities_out^2
markering = markerdup * propmatnormbylogdepth
rownames(markering) = rownames(propmat)
colnames(markering) = colnames(propmat)
markeringlong = melt(as.matrix(markering))
markeringlong = cbind(strsplit2(markeringlong[,2],"[.]"),markeringlong)
markeringlong = markeringlong[,c(5,3,1,2,4)]
colnames(markeringlong) = c("specificity_score","locusID","cluster","subset_cluster","cluster_name")
markeringlong[,3] = gsub("clusters_","",markeringlong[,3])
markeringlong[,4] = gsub("cluster_","",markeringlong[,4])
print("Writing all scores...")
saveRDS(markeringlong,args$allrdsout)
write.table(markeringlong,args$alltxtout,row.names=F,quote=F,sep="\t")
markerlists = markerlistmaker(markering,daps)

print("Calculating specificity threshold...")
nulldist = log10(as.numeric(unlist(markerlists$nonmarkerlist)))
truedist = log10(as.numeric(unlist(markerlists$markerlist)))
pmods = rep(1,times=length(unlist(markerlists$markerlist)))
pmods[which(as.numeric(unlist(markerlists$betamarkerlist)) < 0)] = -1
baremin = max(log10(as.numeric(unlist(markerlists$markerlist)))[which(pmods < 0)])
trueo = sort(truedist[which(truedist >= baremin & pmods > 0)],decreasing = T)
nullo = sort(nulldist[which(nulldist >= baremin & nulldist >= min(trueo))],decreasing = T)
specfloor = spec_floorer(trueo,nullo)

print("Writing specific results...")
results_writer(markerlists$markerlist,specfloor,daps,args$txtout,args$rdsout)
