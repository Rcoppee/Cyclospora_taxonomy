setwd("/media/adm-loc/CryptoDD/RUNS/bam/vcf/VCF.GZ/")

library(vcfR)
library(poppr)
library(ape)
library(RColorBrewer)
library(tidyr)
library(ggplot2)
library(adegenet)
library(dplyr)
library(igraph)

rubi.VCF <- read.vcfR("final3.vcf") #whole-genome
rubi.VCF <- read.vcfR("final_cyclone.vcf") #cyclone

pop.data <- read.table("Continents.txt", sep = "\t", header = TRUE) #continents
#pop.data <- read.table("lineages.txt", sep = "\t", header = TRUE) #lineages
#pop.data <- read.table("Countries.txt", sep = "\t", header = TRUE) #continents

all(colnames(rubi.VCF@gt)[-1] == pop.data$AccessID)

gl.rubi <- vcfR2genlight(rubi.VCF)

ploidy(gl.rubi) <- 2

pop(gl.rubi) <- pop.data$State

tree <- aboot(gl.rubi, tree = "nj", distance = bitwise.dist, sample = 10, showtree = F, cutoff = 50, quiet = T)

pop_levels <- levels(pop(gl.rubi))
#cols <- c("#f21707", "#006619", "#ede000", "#E9DAED", "#45ad18", "#f99500", "#7be0ef", "#2323ff", "#ff869b", "#875A33", "#29A677", "#D600C3", "#b0296c", "#919100") #countries
cols <- c("#2a8a00","#156297", "#eb8808", "#f21707") #continents
cols <- c("#2a8a00", "#156297", "#eb8808") #continents
#cols <- c("#16B6F5","#F59800", "#1DD122") #lineage
pop_cols   <- cols[seq_along(pop_levels)]

plot.phylo(tree, cex = 0.8, font = 2, adj = 0, tip.color =  cols[pop(gl.rubi)])
nodelabels(tree$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.8,font = 3, xpd = TRUE)
legend(
  "topleft",
  legend = pop_levels,
  col    = pop_cols,
  pch    = 19,
  bty    = "n",
  cex    = 1,
  title  = "Population"
)
axis(side = 1)
title(xlab = "Genetic distance (proportion of loci that are different)")

write.tree(tree, file = "continents_WG_tree.nwk")

rubi.pca <- glPca(gl.rubi, nf = 4)
barplot(100*rubi.pca$eig/sum(rubi.pca$eig), col = heat.colors(50), main="PCA Eigenvalues")
title(ylab="Percent of variance\nexplained", line = 2)
title(xlab="Eigenvalues", line = 1)

rubi.pca.scores <- as.data.frame(rubi.pca$scores)
rubi.pca.scores$pop <- pop(gl.rubi)

#use pc1 AND 2 !!! :)
library(ggplot2)
library(ggprism)
set.seed(9)
p <- ggplot(rubi.pca.scores, aes(x=PC1, y=PC2, colour=pop)) +
  geom_point(size=3, show.legend = F) +
  scale_color_manual(values = cols) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme_bw() +
  theme_prism()

p


pnw.dapc <- dapc(gl.rubi, var.contrib = TRUE, scale = FALSE, n.pca = 5, n.da = 3)

scatter( pnw.dapc, col = cols, cex =2, legend = FALSE, clabel = F, scree.pca = FALSE, scree.da = FALSE )

compoplot(pnw.dapc,col = cols, posi = 'top')

dapc.results <- as.data.frame(pnw.dapc$posterior)
dapc.results$pop <- pop(gl.rubi)
dapc.results$indNames <- rownames(dapc.results)

scatter(pnw.dapc, cell = 0, pch = 18:23, cstar = 0, mstree = TRUE, lwd = 2, lty = 2)

contrib <- loadingplot(pnw.dapc$var.contr, axis = 2, thres = 0.07, lab.jitter = 1)


dapc.results <- pivot_longer(dapc.results, -c(pop, indNames))

colnames(dapc.results) <- c("Original_Pop","Sample","Assigned_Pop","Posterior_membership_probability")


p <- ggplot(dapc.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop))
p <- p + geom_bar(stat='identity') 
p <- p + scale_fill_manual(values = cols) 
p <- p + facet_grid(~Original_Pop, scales = "free")
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
p

p <- ggplot(dapc.results, aes(x = Sample, y = Posterior_membership_probability, fill = Assigned_Pop)) +
  geom_col(width = 1, colour = NA, show.legend = FALSE) +
  scale_fill_manual(values = cols) +
  scale_x_discrete(expand = c(0,0)) +
  facet_grid(~Original_Pop, scales = "free") +
  theme_prism() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8)
  )

p
