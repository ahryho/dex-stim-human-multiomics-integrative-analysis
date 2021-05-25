
library(data.table)
library(ggplot2)
library(ggVennDiagram)

data.dir.pre  <- "~/bio/datasets/kimono/"
prior.cpg.gex <- fread(paste0(data.dir.pre, "input_04/prior_cpg_ensg.csv")) # p =0.05, fc = 0
layer.gex     <- fread(paste0(data.dir.pre, "input_04/gex.csv")) # p = 0.05, fc = 0.1

gex.ensg  <- colnames(layer.gex)[-1]
dnam.ensg <- unique(prior.cpg.gex$Ensemble_ID)

ggVennDiagram(list(DMA = dnam.ensg, DEA = gex.ensg), scaled = T) + theme(legend.position = "none") # 791 overlapped genes
# VennDiagram::venn.diagram(list(DMA = dnam.ensg, DEA = gex.ensg), "gex_dnam_overlap.png")

fit <- euler(list(DMA = dnam.ensg, DEA = gex.ensg))
round(fit$original.value["DMA&DEA"] / length(gex.ensg) * 100, 1)
round(fit$original.value["DMA&DEA"] / length(dnam.ensg) * 100, 1)

plot(fit,
     quantities = list(type = c("counts", "percent")),
     fill = c("lightgrey", "orange"))

length(gex.ensg)     