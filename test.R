#### load data #########
source("R/source.R")
pr <- read.table("testdata/all.KEGG.abun.txt",row.names = 1,header = T,sep = "\t")
grp <- read.table("testdata/sampleinfo.txt",row.names = 1,header = T,sep = "\t")


#### D2
grp_sub <- grp[which(grp$TimePoint == "Day2"), 1, drop=F]
grp_sub$Protein <- ifelse(grp_sub$Protein == "C-Pork",
                          "B-Beef", grp_sub$Protein)
pr_sub <- pr[, rownames(grp_sub)]

res_D2 <- ReporterScore(pr_sub, grp_sub, paired = F, database = "./database", occ = 0.1)

#### D7
grp_sub <- grp[which(grp$TimePoint == "Day7"), 1, drop=F]
grp_sub$Protein <- ifelse(grp_sub$Protein == "C-Pork",
                          "B-Beef", grp_sub$Protein)
pr_sub <- pr[, rownames(grp_sub)]

res_D7 <- ReporterScore(pr_sub, grp_sub, paired = F, database = "./database", occ = 0.1)


#### D14
grp_sub <- grp[which(grp$TimePoint == "Day14"), 1, drop=F]
grp_sub$Protein <- ifelse(grp_sub$Protein == "C-Pork",
                          "B-Beef", grp_sub$Protein)
pr_sub <- pr[, rownames(grp_sub)]

res_D14 <- ReporterScore(pr_sub, grp_sub, paired = F, database = "./database", occ = 0.1)



