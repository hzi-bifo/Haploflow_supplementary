#!/usr/bin/env Rscript

# takes a labelled tree as well as the annotation file with reconstructed locations and counts how often each location transmits to a new one

library("ape")
library("plyr")

args = commandArgs(trailingOnly=TRUE)
phylo_file = args[1]
annotation_file = args[2]
#output_file = args[3]

#phylo_file <- "/home/sreimering/repos/VirusTracker/augur_adjustment/cov2020_reconstructed.phy"
#annotation_file <- "/home/sreimering/repos/VirusTracker/augur_adjustment/cov2020_reconstructed.annotation.txt"
#output_file <- "/home/sreimering/repos/VirusTracker/augur_adjustment/cov2020_links.txt"

tree <- read.tree(phylo_file)
annotation <- read.table(annotation_file, sep=" ", header=TRUE, stringsAsFactors=FALSE)

# get all labels in the order of the ids
all_label <- c(tree$tip.label, tree$node.label)

# get all locations
all_locations <- unique(annotation$location)

changes <- list()
num_changes <- 0

# go through all edges in the tree
for (i in 1:nrow(tree$edge)){
  # nodes in the edge matrix are represented by ids; use the id to get the corresponding labels stored in all_label
  origin_label <- all_label[tree$edge[i, 1]]
  destination_label <- all_label[tree$edge[i, 2]]
  
  # get the location assigned to the origin and destination nodes
  origin <- annotation$location[which(annotation$label == origin_label)]
  destination <- annotation$location[which(annotation$label == destination_label)]
 
  local <- c()
  for (j in 1:nchar(origin)) # have the same length by definition
  {
    pre <- substr(origin,j,j)
    post <- substr(destination,j,j)
    if (pre != post)
    {
        local <- c(local, paste(pre, j, post, sep=""))
        num_changes <- num_changes + 1
    }
  }
  if (is.null(local))
  {
    changes[[i]] <- c("")
  }
  else
  {
    changes[[i]] <- local
  }
  names(changes)[i] <- paste(c(origin_label, destination_label), collapse='    ')
}

print(num_changes)
print(changes)
#changes.df = do.call("rbind", lapply(changes, as.data.frame))
#write.table(changes, file = output_file, sep='\t')
