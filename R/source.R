#' WilcoxTest
#' wilcox test or pair wilcox test
#' @param data, {data.frame or matrix}, ko profile: row is KO ID, col is sample ID
#' @param grp, {data.frame}, group information, two class
#' @param paired , True or FALSE
#'
#' @return Wilcoxon test result
#'
#' @export
#'
#' @examples
WilcoxTest <- function(data, grp, paired = FALSE){

  grp$lv <- as.factor(grp[,1])
  glv <- as.numeric(grp$lv)-1

  count <- apply(data, 1, function(x){sum(as.numeric(x)>0)})
  data <- data[which(count > 0),]
  name <- rownames(data)

  data1 <- t(scale(t(data)))

  logtest<-function(marker,outcome1){
    model <- glm(outcome1~marker, family=binomial())
    coef(summary(model))["marker", c("Estimate", "Std. Error", "Pr(>|z|)")]
  }
  ####### to compute the OR value
  or1 <- apply(data1, 1, function(x){logtest(as.numeric(x), glv)})
  or1 <- t(or1)
  OR   <- exp(as.numeric(or1[ ,1]))
  lower<- exp(as.numeric(or1[ ,1]) - 1.96*as.numeric(or1[, 2]))
  uper <- exp(as.numeric(or1[, 1]) + 1.96*as.numeric(or1[, 2]))

  or2  <- data.frame(round(OR, 4), round(lower, 4),round(uper, 4),
                     com = paste(round(OR,2),"(",round(lower,2),",",round(uper,2),")",sep=""))
  rownames(or2) <- name

  ######## wilcoxon test
  num <- nrow(data)
  res <- matrix(0, nrow=num, ncol=9)
  rownames(res) <- rownames(data)
  rlv <- rev(levels(grp$lv)) # here to rev the level

  for(i in 1:num){
    a <- as.numeric(data[i,])
    res[i,1] <- wilcox.test(a[glv==1], a[glv==0], paired = paired)$p.value
    r <- rank(a)
    res[i, 2] <- round(mean(r[glv == 1]))
    res[i, 3] <- round(mean(r[glv == 0]))
    res[i, 4] <- ifelse(res[i,3] > res[i,2], rlv[1], rlv[2])
    #### to compute the zscore, if enrich == 0, negtive ; enrich ==1, positive
    if( res[i,4] == rlv[2]) {
      res[i,5] <-  round(qnorm(1-(wilcox.test(a[glv==1], a[glv==0],
                        exact=F, paired=paired)$p.value)/2), 2)
    }else{
      res[i,5] <- round(qnorm((wilcox.test(a[glv==1], a[glv==0],
                        exact=F, paired=paired)$p.value)/2), 2)
    }

    res[i,6] = round(sum(a[glv == 1] > 0)/sum(glv == 1), 4)
    res[i,7] = round(sum(a[glv == 0] > 0)/sum(glv == 0), 4)
    res[i,8] <- round(mean(a[glv == 1]), 10)
    res[i,9] <- round(mean(a[glv == 0]), 10)

  }

  qvalue <- p.adjust(res[,2],method="BH")

  out.dat <- cbind(res, or2, qvalue)

  colnames(out.dat) <- c("p.value", paste0("rankIn", rev(levels(grp$lv))), "enrichedIn",
                         "p.1tail", paste0("occ.In", rev(levels(grp$lv))),
                          paste0("meanIn",rev(levels(grp$lv))),
                         "exp.OR","exp.low","exp.up","com","qvalue(BH)")
  return(out.dat)

}

#' GetCoverage
#' To compute the coverage of module or pathway in the current ko profile
#' @param pro, ko profile
#' @param con, group information
#' @param database , database path
#'
#' @return
#'
#' @export
#'
#' @examples
GetCoverage <- function(pro, con, database){

  ptw <- read.table(paste0(database, "/pathway.annotation.txt"), row.names=1,header=T,check.names=F)
  modul <- read.table(paste0(database, "/module.annotation.txt"), row.names=1,header=T,check.names=F)

  # Group will sort alphabetically
  grp <- as.factor(con[,1])
  gn <- levels(grp)
  glv <- as.numeric(grp) -1

  id0 <- which(con[,1]==gn[1])
  id1 <- which(con[,1]==gn[2])
  p0 <- pro[,rownames(con)[id0]]
  p1 <- pro[,rownames(con[id1,,drop=F])]

  # cutoff for ko coverage in each group : 0.1
  cutoff <- function(x){
    id <- x > 0
    n <- length(x[id])
    r <- n/length(x)
    if(r > 0.1){
      row_id <- 1
    }else{
      row_id <- 0
    }
    return(row_id)
  }

  pp0 <- apply(p0,1,cutoff)
  pp1 <- apply(p1,1,cutoff)

  p0 <- p0[which(pp0==1),]
  p1 <- p1[which(pp1==1),]

  # compute coverage of pathway and module
  cov <- function(x){
    coverage <- matrix(NA, nrow(x),2)
    for(i in 1:nrow(x)) {
      ko <- unlist(strsplit(as.character(x[i,ncol(x)]),split=","))
      n0 <- length(intersect(ko,rownames(p0)))
      n1 <- length(intersect(ko,rownames(p1)))
      c0 <- n0/x[i,ncol(x)-1]
      c1 <- n1/x[i,ncol(x)-1]
      coverage[i,] <- c(c0,c1)
    }
    return (coverage)
  }

  rownames(ptw) <- ptw$PathwayID
  rownames(modul) <- modul$ModuleID
  coverage.path <- cov(ptw)
  coverage.modul <- cov(modul)

  #names
  rownames(coverage.path) <- rownames(ptw)
  rownames(coverage.modul) <- rownames(modul)
  colnames(coverage.path) <- paste0("coverage_",gn)
  colnames(coverage.modul) <- paste0("coverage_",gn)

  coverage.path <- cbind(ptw[rownames(coverage.path),],coverage.path)
  coverage.modul <- cbind(modul[rownames(coverage.modul),],coverage.modul)

  list(pathwayCoverage=coverage.path,moduleCoverage=coverage.modul)
}





#' GetZscore
#' to compute the zscore of wilcoxon test
#' @param dat
#' @param con
#' @param paired
#'
#' @return
#' @export
#'
#' @examples
Getzscore <- function(dat, con, paired = F, hist = F, adjust = T){

  dat <- dat[rowSums(dat)!=0, ]
  glv <- as.factor(con[,1])
  grp <- as.numeric(glv) -1
  grp.f <- grp[!is.na(grp)]
  dat.f <- dat[,!is.na(grp)]
  num <- nrow(dat.f)
  nsample <- c(length(grp.f[grp.f==0]),length(grp.f[grp.f==1]))

  # result...
  res <- matrix(0,nrow=num,ncol=16)
  rownames(res) <- rownames(dat.f)

  # test...
  pairinfo<-as.logical(paired)
  for(i in 1:num){
    a <- as.numeric(dat.f[i,])
    res[i,1] <- wilcox.test(a[grp.f==0], a[grp.f==1], paired=pairinfo, alternative="two.sided")$p.value
    res[i,2] <- wilcox.test(a[grp.f==0], a[grp.f==1], paired=pairinfo, alternative="less")$p.value
    onetail <-  res[i,2]
    res[i,3] <- wilcox.test(a[grp.f==0],a[grp.f==1], paired=pairinfo, alternative="greater")$p.value
    res[i,10] <- mean(a[grp.f==0])
    res[i,11] <- mean(a[grp.f==1])
    res[i,12] <- sd(a[grp.f==0])
    res[i,13] <- sd(a[grp.f==1])
    res[i,14] <- sum(a[grp.f==0]>0)
    res[i,15] <- sum(a[grp.f==1]>0)

    if(res[i,1]<0.05){
      if(onetail < 0.05){
        res[i,16] <- 2  # 0<1
      }else{
        res[i,16] <- 1  # 0>1
      }
    }else{
        res[i,16] <- 0  # non
    }
  }
  #p.adjust-->q.value and get z-score
  dat <- res
  zmax <- 8.20953
  zmin <- -8.20953
  for(i in 1:3){
    x <- dat[,i] # pvalue
    x <- x[!is.na(x)]
    res <- res[!is.na(res[,1]),]
    if(adjust){ # add adjust parameter to decide if to adjust the pvalue
      x.a <- p.adjust(x, method = "BH")
    }else{
      x.a <- x
    }
    res[,3+i] <- x.a
    z.a <- qnorm(1-x) # this is fit for less; twotail +/-qnorm(1-p/2); righttail qnorm(p)
    for(j in 1:length(z.a))
    {
      if((z.a[j]>0) && is.infinite(z.a[j])){
        z.a[j] <- zmax
      }
      if((z.a[j]<0) && is.infinite(z.a[j])){
        z.a[j] <- zmin
      }
    }
    res[,6+i] <- z.a
  }
  for(i in 1:dim(res)[1]){
    a <- res[i,]
    if(a[4] > 0.05){
      a[16] <- 0
    }
    res[i,] <- a
  }

  title <- c("p-value(twosided)","p-value(less)", "p-value(greater)",
             "p.adjust(twosided)","p.adjust(less)","p.adjust(greater)",
             "z-score(twosided)","z-score(less)","z-score(greater)",
             paste0(c("mean@","sd@"),rep(levels(glv),2)),paste0("occ.@",levels(glv)),"enrichment_dir")
  colnames(res) <- title

  #####hist####
  dat <- res
  lab = hist(dat[,1],plot=F,breaks=20)
  Peak <- lab$mids[which.max(lab$density)]
  sd <- sd(dat[,1])
  main = paste("P.value.Peak=",Peak,".Sd=",round(sd,3),sep="")
  outname <- "P.value.hist.pdf"
  #pdf(outname)
  if(hist){
      hist(dat[,1], freq=F,breaks=20,xlab="p.value",main=main)
    }
  #dev.off()
  ####FDR####
  x <- dat[,1]
  x <- x[!is.na(x)]
  x.a <- p.adjust(x,method = "BH")
  y <- cbind(x,x.a)
  r <- order(x)
  y <- y[r,]
  fdr.1 <- y[max(which(y[,1] < 0.01)),]
  fdr.5 <- y[max(which(y[,1] < 0.05)),]
  fdr=data.frame(c(0.01,0.05),c(fdr.1[2],fdr.5[2]))
  colnames(fdr) = c("Significant_level","FDR")
  #outname <- paste(prefix,"FDR.txt",sep=".")
  #write.table(fdr,file=outname,quote=F,sep="\t",row.names=F)

  res
}

#' CoverNumer
#'
#' @param dat, KO profile
#' @param cov, group infromation
#'
#' @return
#' @export
#'
#' @examples
CoverNumber <- function(dat, cov){

  dat <- as.data.frame(dat)
  cov$n_cov <- NA
  cov$r_cov <- NA
  for (x in 1:nrow(cov)){
    ko <- unlist(strsplit(cov$KO[x], ","))
    n <- sum(rownames(dat) %in% ko)
    cov$r_cov[x] <- n/cov$n[x]
    cov$n_cov[x] <- n
  }
  cov[cov$n_cov>0, ]
}

#' PermZscore
#' permutation of szscore value
#' @param dat, zscore result
#' @param cov,
#' @param perm, permuation number
#'
#' @return
#' @export
#'
#' @examples
#'
PermZscore <- function(dat, cov, perm = 1000){

  kolst <- lapply(cov$KO, function(x) unlist(strsplit(x, ",")))
  kolst <- unlist(kolst)
  #kolst <- kolst[duplicated(kolst)]
  kolst <- unique(kolst)
  dat <- dat[rownames(dat) %in% kolst,]

  zcol <- 7:9
  zscore <- dat[, zcol]  #in zscore, the first column will be two.sided, second will be less, third will be greater
  names(zscore) <- rownames(dat)
  ds <- dim(dat)

  zscorenum <- dim(zscore)[2]
  #klist
  klist <- union(cov$n_cov, NULL)
  knum <- length(klist)
  #permute
  res <- matrix(0, nrow=knum, ncol=2*zscorenum+1)

  #start permute
  for(j in 1:knum){
    k <- klist[j]
    b <- rep(0, perm)
    res[j, 1] <- k
    for(t in 1:zscorenum){	#number of selected zscore types
      for(i in 1:perm){
        a <- zscore[sample(ds[1], k), t] ### random select k zscore
        b[i] <- sum(a)/sqrt(k)
      }
      res[j,t*2] <- mean(b)   #mean
      res[j,(t*2+1)] <- sd(b)	#sd
    }
  }
  rownames(res) <- res[,1]
  res[,-1]
}

#' GetReporterScore
#' Computer the report score
#' @param dat
#' @param perm
#' @param cov
#' @param occ, cutoff of coverage on pathway or module
#'
#' @return
#' @export
#'
#' @examples
GetReporterScore <- function(dat, perm, cov, occ = 0.4){

  cov$zscore1 <- NA
  cov$zscore2 <- NA
  cov$zscore3 <- NA
  cov$zscore <- NA

  for(x in 1:nrow(cov)){

    ##### to compute the
    ko <- unlist(strsplit(cov$KO[x],","))
    tmp <- dat[rownames(dat) %in% ko, ,drop=F]
    n_cov <- cov$n_cov[x]
    zscore <- tmp[, grep("z-score", colnames(tmp)), drop=F]
    for(i in 1:3){
      zscore[,i] <- zscore[,i]/sqrt(n_cov)
    }
    zscore <- colSums(zscore)
    occ0 <- tmp[, grep("occ",colnames(tmp)),drop=F] ## here
    occ1 <- sum(occ0[,1] > 0)
    occ2 <- sum(occ0[,2] > 0)
    a <- perm[as.character(n_cov), ] # perm zscore
    for(i in 1:3){
      if(a[i*2] == 0){  # sd if zero
        zscore[i] <- NA
      }else{
        zscore[i] <- (zscore[i]-a[i*2-1])/a[i*2]
      }
    }
    cov[x, (ncol(cov)-3):(ncol(cov)-1)] <- zscore


    if(occ1/cov$n[x] >= occ | occ2/cov$n[x] >= occ){
      cov$zscore[x] <- ifelse(zscore[2] > 0, zscore[2], -zscore[3])
    }
  }
  cov
}

#'  ReporterScore
#'  The command to compute the reportscore
#' @param pr
#' @param grp
#' @param paired
#' @param database
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
ReporterScore <- function(pr, grp, paired = F, database = "./database", occ = 0.1, adjust =T){
  # step1 test
  res.test <- WilcoxTest(pr,grp,paired = paired)

  # step2 coverage
  res.cov <- GetCoverage(pr,grp, database)

  # step3 z-score
  res.z <- Getzscore(pr,grp,paired = paired, adjust)

  # step4 cover number
  cov.pw <- CoverNumber(res.z, res.cov$pathwayCoverage)
  cov.md <- CoverNumber(res.z, res.cov$moduleCoverage)

  # step5 perm
  perm.ptw <- PermZscore(res.z, cov.pw, perm = 1000)
  perm.mod <- PermZscore(res.z, cov.md, perm = 1000)

  # step6 reporter score
  out.ptw <- GetReporterScore(res.z, perm.ptw, cov.pw, occ)
  out.mod <- GetReporterScore(res.z, perm.mod, cov.md, occ)

  list(pathway = out.ptw, module = out.mod)

}
