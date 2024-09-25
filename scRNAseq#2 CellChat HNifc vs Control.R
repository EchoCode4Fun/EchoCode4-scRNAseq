#Set up control and HNifc cellchat

library(CellChat)
library(Seurat)
library(patchwork)

# Convert Seurat to CellChat
cellchat_HNifc <- createCellChat(object = HNifc_seurat_object, group.by = "ident")
cellchat_Control <- createCellChat(object = control_seurat_object, group.by = "ident")

# Set the database for communication
CellChatDB <- CellChatDB.mouse # or CellChatDB.mouse if using mouse data
cellchat_HNifc@DB <- CellChatDB
cellchat_Control@DB <- CellChatDB

# Subset the expression data of signaling genes
cellchat_HNifc <- subsetData(cellchat_HNifc)
cellchat_Control <- subsetData(cellchat_Control)

# Identify overexpressed genes and interactions
cellchat_HNifc <- identifyOverExpressedGenes(cellchat_HNifc)
cellchat_Control <- identifyOverExpressedGenes(cellchat_Control)

cellchat_HNifc <- identifyOverExpressedInteractions(cellchat_HNifc)
cellchat_Control <- identifyOverExpressedInteractions(cellchat_Control)

#Compute commun prob
cellchat_HNifc <- computeCommunProb(cellchat_HNifc)
cellchat_Control <- computeCommunProb(cellchat_Control)

# Filter out low probability interactions
cellchat_HNifc <- filterCommunication(cellchat_HNifc, min.cells = 10)
cellchat_Control <- filterCommunication(cellchat_Control, min.cells = 10)

#Infer cellular communication networks 
cellchat_HNifc <- computeCommunProbPathway(cellchat_HNifc)
cellchat_Control <- computeCommunProbPathway(cellchat_Control)

cellchat_HNifc <- aggregateNet(cellchat_HNifc)
cellchat_Control <- aggregateNet(cellchat_Control)


# Compare communication networks between HNifc and Control CellChat

#########The code below is for two datasets with different cell compositions 

library(CellChat)
library(patchwork)

##The goal is to make two datasets (slightly different) contain the same cell composition 
##For vastly different cell groups, some functions do not work # Check and print cell types in HNifc_cellchat
##The idea is to lift up the ones with less compositions to match the one with more 

# Check and print cell types in HNifc_cellchat
HNifc_cell_types <- levels(cellchat_HNifc@idents)
print("Cell types in HNifc_cellchat:")
print(HNifc_cell_types)

# Check and print cell types in Control_cellchat
control_cell_types <- levels(cellchat_Control@idents)
print("Cell types in Control_cellchat:")
print(control_cell_types)

# Define the combined cell types from both objects
group.new <- union(levels(cellchat_HNifc@idents), levels(cellchat_Control@idents))

# Lift up the cellchat objects to match the combined cell types
cellchat_HNifc <- liftCellChat(cellchat_HNifc, group.new)
cellchat_Control <- liftCellChat(cellchat_Control, group.new)

# Merge the CellChat objects
object.list <- list(HNifc = cellchat_HNifc, Control = cellchat_Control)
cellchat <- mergeCellChat(object.list, add.names = c("HNifc", "Control"), cell.prefix = TRUE)

# Print the merged CellChat object
print(cellchat)

# Check and print cell types in HNifc_cellchat
HNifc_cell_types <- levels(cellchat_HNifc@idents)
print("Cell types in HNifc_cellchat:")
print(HNifc_cell_types)

# Check and print cell types in Control_cellchat
control_cell_types <- levels(cellchat_Control@idents)
print("Cell types in Control_cellchat:")
print(control_cell_types)

# Check if the cell types are the same
if (setequal(HNifc_cell_types, control_cell_types)) {
  print("Both CellChat objects have the same cell composition.")
} else {
  print("The CellChat objects have different cell compositions.")
}


#########The code below is for two datasets with the same cell composition

#Compare the total number of interactions and interaction strength
ptm = Sys.time()
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

##Compare the number of interactions and interaction strength among different cell populations
#Form 1: circle plot
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

#Form2: heatmap
gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2

#Form3: circle plot (for the number of interactions and interaction weight)
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

# Get the maximum weight for scaling the interaction weights
weight.max <- getMaxWeight(object.list, attribute = c("idents", "weight"))

# Set up the plotting area
par(mfrow = c(1, 2), xpd = TRUE)

# Loop through each CellChat object to visualize
for (i in 1:length(object.list)) {
  netVisual_circle(
    object.list[[i]]@net$weight,  # Use 'weight' instead of 'count'
    weight.scale = TRUE, 
    label.edge = FALSE, 
    edge.weight.max = weight.max[2], 
    edge.width.max = 12, 
    title.name = paste0("Interaction weights - ", names(object.list)[i])
  )
}


# Form 4: Circle plot (only shown the interactions between certain cell types)

#Form4: circle plot (only shown the interactions between certain cell types)
# Define the cell types of interest
group.cellType <- c(rep("Hematopoietic Stem/Progenitor Cells", 1), rep("Neutrophils", 1), rep("Macrophages", 1))

# Set the factor levels to match the cell types you're interested in
group.cellType <- factor(group.cellType, levels = c("Hematopoietic Stem/Progenitor Cells", "Neutrophils", "Macrophages"))

# Merge interactions for the specified cell types
object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})

# Merge the CellChat objects
cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)

# Get the maximum weight for scaling the interaction weights in the merged data
weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents", "count", "count.merged"))

# Set up the plotting area
par(mfrow = c(1, 2), xpd = TRUE)

# Loop through each CellChat object to visualize
for (i in 1:length(object.list)) {
  netVisual_circle(
    object.list[[i]]@net$count.merged, 
    weight.scale = TRUE, 
    label.edge = TRUE, 
    edge.weight.max = weight.max[3], 
    edge.width.max = 12, 
    title.name = paste0("Number of interactions - ", names(object.list)[i])
  )
}

# Calculate max interaction strength for EACH DATASET
weight.max.strength <- sapply(object.list, function(x) {
  max(x@net$weight.merged)
})

# Loop through each dataset in object.list
for (i in seq_along(object.list)) {
  # Plot interaction strength for the current dataset
  netVisual_circle(
    object.list[[i]]@net$weight.merged, 
    weight.scale = TRUE, 
    label.edge = TRUE, 
    edge.weight.max = weight.max.strength[i], # Use max strength for current dataset
    edge.width.max = 12, 
    title.name = paste0("Interaction strength - ", names(object.list)[i])
  )
}

##Compare the major sources and targets in a 2D space

#1. identify cell populations with significant changes in sending or receiving signals
#scatter plot
# Compute the network centrality scores for each object in the list
object.list <- lapply(object.list, netAnalysis_computeCentrality)

# Proceed with the scatter plot
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count) - diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # Control the dot size in the different datasets
gg <- list()

for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}

# Plot the results
gridExtra::grid.arrange(grobs = gg, ncol = 2)


#2. Identify the signalling changes of specific cell populations
##(the code below does not work...)
# Compute centrality scores for the CellChat object
cellchat <- netAnalysis_computeCentrality(cellchat)

# Check Neutrophils
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Activated Neutrophils", signaling.exclude = "MIF")

# Check Stem cells
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Hematopoietic Stem/Progenitor Cellss", signaling.exclude = "MIF")

# Combine the plots
patchwork::wrap_plots(plots = list(gg1, gg2))


#3. Identify altered signaling with distinct interaction strength
##3.1Compare the overall information flow of each signaling pathway or ligand-receptor pair
gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE)

gg1 + gg2

##3.2 Compare outgoing (or incoming) signaling patterns associated with each cell population

library(ComplexHeatmap)
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

#4.Identify the up-gulated and down-regulated signaling ligand-receptor pairs

##4.1.Identify dysfunctional signaling by comparing the communication probabities
#> Comparing communications on a merged object
#> > print(HNifc_cell_types)
#
#[1] "Mature Neutrophils"                    
#[2] "Activated Neutrophils"                 
#[3] "Hematopoietic Stem/Progenitor Cells"   
#[4] "Immature Neutrophils"                  
#[5] "pre B Cells"                           
#[6] "Megakaryocytes/Platelets"              
#[7] "Dendritic Cells"                       
#[8] "Fibroblasts or Mesenchymal Cells"      
#[9] "Erythroid Cells/Erythrocyte Precursors"
#[10] "Pro-B Cells"                           
#[11] "Macrophages Monocytes"                 
#[12] "Antimicrobial Neutrophils"             
#[13] "Undefined"                             
#[14] "Mast Cells/Progenitors"                
#[15] "Natural Killer Cells"                  
#[16] "Plasma Cells"       

#[1] "Hematopoietic Stem/Progenitor Cells" "Neutrophils"                        
#[3] "Macrophages"                         "Megakaryocytes"                     
#[5] "pro-B cells"                         "Basophils"                          
#[7] "Erythroblasts"                       "Mast Cell Progenitors" 


# activated Neutrophils as senders 
ptm = Sys.time()

netVisual_bubble(cellchat, sources.use = 2, targets.use = c(1:8),  comparison = c(1, 2), angle.x = 45)

#anti-microbial neutrophils as senders
ptm = Sys.time()

netVisual_bubble(cellchat, sources.use = 1, targets.use = c(1:8),  comparison = c(1, 2), angle.x = 45)

#mature neutrophils as senders
ptm = Sys.time()

netVisual_bubble(cellchat, sources.use = 3, targets.use = c(1:8),  comparison = c(1, 2), angle.x = 45)

#Immature neutrophils as senders
ptm = Sys.time()

netVisual_bubble(cellchat, sources.use = 4, targets.use = c(1:8),  comparison = c(1, 2), angle.x = 45)


# Stem cells as senders 
ptm = Sys.time()

netVisual_bubble(cellchat, sources.use = 3, targets.use = c(1:20),  comparison = c(1, 2), angle.x = 45)

# Monocytes as senders 
ptm = Sys.time()

netVisual_bubble(cellchat, sources.use = 11, targets.use = c(1:20),  comparison = c(1, 2), angle.x = 45)



##More specifically for the upregulated and the downregulated 
gg1 <- netVisual_bubble(cellchat, sources.use = 3, targets.use = c(1:16),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in HNifc group", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = 3, targets.use = c(1:16),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in HNifc group", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2


###Visually compare cell-cell communication with three different plots
##Compare the signalling gene expression distribution between different datasets 
#Violin plot
cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("HNifc", "Control")) # set factor level
plotGeneExpression(cellchat, signaling = "CXCL", split.by = "datasets", colors.ggplot = T, type = "violin")

#Circle plots
pathways.show <- c("CXCL") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

#heatmap
pathways.show <- c("CXCL") 
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}

ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))

#chord diagram 
# Chord diagram
pathways.show <- c("CXCL") 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
}


save.image(file = "scRNAseq#2CellChat.RData")


#---------shown below is the original code

#########The code below is for two datasets with the same cell composition

#Compare the total number of interactions and interaction strength
ptm = Sys.time()
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

##Compare the number of interactions and interaction strength among different cell populations
#Form 1: circle plot
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

#Form2: heatmap
gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2

#Form3: circle plot (for the number of interactions and interaction weight)
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

# Get the maximum weight for scaling the interaction weights
weight.max <- getMaxWeight(object.list, attribute = c("idents", "weight"))

# Set up the plotting area
par(mfrow = c(1, 2), xpd = TRUE)

# Loop through each CellChat object to visualize
for (i in 1:length(object.list)) {
  netVisual_circle(
    object.list[[i]]@net$weight,  # Use 'weight' instead of 'count'
    weight.scale = TRUE, 
    label.edge = FALSE, 
    edge.weight.max = weight.max[2], 
    edge.width.max = 12, 
    title.name = paste0("Interaction weights - ", names(object.list)[i])
  )
}


# Form 4: Circle plot (only shown the interactions between certain cell types)

#Form4: circle plot (only shown the interactions between certain cell types)
# Define the cell types of interest
group.cellType <- c(rep("Hematopoietic Stem/Progenitor Cells", 1), rep("Antimicrobial Neutrophils", 1), rep("Macrophages Monocytes", 1))

# Set the factor levels to match the cell types you're interested in
group.cellType <- factor(group.cellType, levels = c("Hematopoietic Stem/Progenitor Cells", "Antimicrobial Neutrophils", "Macrophages Monocytes"))

# Merge interactions for the specified cell types
object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})

# Merge the CellChat objects
cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)

# Get the maximum weight for scaling the interaction weights in the merged data
weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents", "count", "count.merged"))

# Set up the plotting area
par(mfrow = c(1, 2), xpd = TRUE)

# Loop through each CellChat object to visualize
for (i in 1:length(object.list)) {
  netVisual_circle(
    object.list[[i]]@net$count.merged, 
    weight.scale = TRUE, 
    label.edge = TRUE, 
    edge.weight.max = weight.max[3], 
    edge.width.max = 12, 
    title.name = paste0("Number of interactions - ", names(object.list)[i])
  )
}

# Calculate max interaction strength for each dataset
weight.max.strength <- sapply(object.list, function(x) {
  max(x@net$weight.merged)
})

# Plot interaction strength
netVisual_circle(
  object.list[[i]]@net$weight.merged, 
  weight.scale = TRUE, 
  label.edge = TRUE, 
  edge.weight.max = weight.max.strength[3], 
  edge.width.max = 12, 
  title.name = paste0("Interaction strength - ", names(object.list)[i])
)


##Compare the major sources and targets in a 2D space

#1. identify cell populations with significant changes in sending or receiving signals
#scatter plot
# Compute the network centrality scores for each object in the list
object.list <- lapply(object.list, netAnalysis_computeCentrality)

# Proceed with the scatter plot
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count) - diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # Control the dot size in the different datasets
gg <- list()

for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}

# Plot the results
gridExtra::grid.arrange(grobs = gg, ncol = 2)


#2. Identify the signalling changes of specific cell populations
##(the code below does not work...)
# Compute centrality scores for the CellChat object
cellchat <- netAnalysis_computeCentrality(cellchat)

# Check Neutrophils
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Activated Neutrophils", signaling.exclude = "MIF")

# Check Stem cells
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Hematopoietic Stem/Progenitor Cellss", signaling.exclude = "MIF")

# Combine the plots
patchwork::wrap_plots(plots = list(gg1, gg2))


#3. Identify altered signaling with distinct interaction strength
##3.1Compare the overall information flow of each signaling pathway or ligand-receptor pair
gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE)

gg1 + gg2

##3.2 Compare outgoing (or incoming) signaling patterns associated with each cell population

library(ComplexHeatmap)
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

#4.Identify the up-gulated and down-regulated signaling ligand-receptor pairs

##4.1.Identify dysfunctional signaling by comparing the communication probabities
#> Comparing communications on a merged object
#> > print(HNifc_cell_types)
#
#[1] "Mature Neutrophils"                    
#[2] "Activated Neutrophils"                 
#[3] "Hematopoietic Stem/Progenitor Cells"   
#[4] "Immature Neutrophils"                  
#[5] "pre B Cells"                           
#[6] "Megakaryocytes/Platelets"              
#[7] "Dendritic Cells"                       
#[8] "Fibroblasts or Mesenchymal Cells"      
#[9] "Erythroid Cells/Erythrocyte Precursors"
#[10] "Pro-B Cells"                           
#[11] "Macrophages Monocytes"                 
#[12] "Antimicrobial Neutrophils"             
#[13] "Undefined"                             
#[14] "Mast Cells/Progenitors"                
#[15] "Natural Killer Cells"                  
#[16] "Plasma Cells"       

# activated Neutrophils as senders 
ptm = Sys.time()

netVisual_bubble(cellchat, sources.use = 2, targets.use = c(1:20),  comparison = c(1, 2), angle.x = 45)

#anti-microbial neutrophils as senders
ptm = Sys.time()

netVisual_bubble(cellchat, sources.use = 12, targets.use = c(1:20),  comparison = c(1, 2), angle.x = 45)

#mature neutrophils as senders
ptm = Sys.time()

netVisual_bubble(cellchat, sources.use = 1, targets.use = c(1:20),  comparison = c(1, 2), angle.x = 45)

#Immature neutrophils as senders
ptm = Sys.time()

netVisual_bubble(cellchat, sources.use = 4, targets.use = c(1:20),  comparison = c(1, 2), angle.x = 45)


# Stem cells as senders 
ptm = Sys.time()

netVisual_bubble(cellchat, sources.use = 3, targets.use = c(1:20),  comparison = c(1, 2), angle.x = 45)

# Monocytes as senders 
ptm = Sys.time()

netVisual_bubble(cellchat, sources.use = 11, targets.use = c(1:20),  comparison = c(1, 2), angle.x = 45)



##More specifically for the upregulated and the downregulated 
gg1 <- netVisual_bubble(cellchat, sources.use = 3, targets.use = c(1:16),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in HNifc group", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = 3, targets.use = c(1:16),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in HNifc group", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2


###Visually compare cell-cell communication with three different plots
##Compare the signalling gene expression distribution between different datasets 
#Violin plot
cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("HNifc", "Control")) # set factor level
plotGeneExpression(cellchat, signaling = "CXCL", split.by = "datasets", colors.ggplot = T, type = "violin")

#Circle plots
pathways.show <- c("CXCL") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

#heatmap
pathways.show <- c("CXCL") 
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}

ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))

#chord diagram 
# Chord diagram
pathways.show <- c("CXCL") 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
}


save.image(file = "scRNAseq#2CellChat.RData")











#########The code below is for two datasets with the same cell composition

#Compare the total number of interactions and interaction strength
ptm = Sys.time()
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

##Compare the number of interactions and interaction strength among different cell populations
#Form 1: circle plot
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

#Form2: heatmap
gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2

#Form3: circle plot (for the number of interactions and interaction weight)
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

# Get the maximum weight for scaling the interaction weights
weight.max <- getMaxWeight(object.list, attribute = c("idents", "weight"))

# Set up the plotting area
par(mfrow = c(1, 2), xpd = TRUE)

# Loop through each CellChat object to visualize
for (i in 1:length(object.list)) {
  netVisual_circle(
    object.list[[i]]@net$weight,  # Use 'weight' instead of 'count'
    weight.scale = TRUE, 
    label.edge = FALSE, 
    edge.weight.max = weight.max[2], 
    edge.width.max = 12, 
    title.name = paste0("Interaction weights - ", names(object.list)[i])
  )
}


# Form 4: Circle plot (only shown the interactions between certain cell types)

#Form4: circle plot (only shown the interactions between certain cell types)
# Define the cell types of interest
group.cellType <- c(rep("Hematopoietic Stem/Progenitor Cells", 1), rep("Antimicrobial Neutrophils", 1), rep("Macrophages Monocytes", 1))

# Set the factor levels to match the cell types you're interested in
group.cellType <- factor(group.cellType, levels = c("Hematopoietic Stem/Progenitor Cells", "Antimicrobial Neutrophils", "Macrophages Monocytes"))

# Merge interactions for the specified cell types
object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})

# Merge the CellChat objects
cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)

# Get the maximum weight for scaling the interaction weights in the merged data
weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents", "count", "count.merged"))

# Set up the plotting area
par(mfrow = c(1, 2), xpd = TRUE)

# Loop through each CellChat object to visualize
for (i in 1:length(object.list)) {
  netVisual_circle(
    object.list[[i]]@net$count.merged, 
    weight.scale = TRUE, 
    label.edge = TRUE, 
    edge.weight.max = weight.max[3], 
    edge.width.max = 12, 
    title.name = paste0("Number of interactions - ", names(object.list)[i])
  )
}

# Calculate max interaction strength for each dataset
weight.max.strength <- sapply(object.list, function(x) {
  max(x@net$weight.merged)
})

# Plot interaction strength
netVisual_circle(
  object.list[[i]]@net$weight.merged, 
  weight.scale = TRUE, 
  label.edge = TRUE, 
  edge.weight.max = weight.max.strength[3], 
  edge.width.max = 12, 
  title.name = paste0("Interaction strength - ", names(object.list)[i])
)


##Compare the major sources and targets in a 2D space

#1. identify cell populations with significant changes in sending or receiving signals
#scatter plot
# Compute the network centrality scores for each object in the list
object.list <- lapply(object.list, netAnalysis_computeCentrality)

# Proceed with the scatter plot
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count) - diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # Control the dot size in the different datasets
gg <- list()

for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}

# Plot the results
gridExtra::grid.arrange(grobs = gg, ncol = 2)


#2. Identify the signalling changes of specific cell populations
##(the code below does not work...)
# Compute centrality scores for the CellChat object
cellchat <- netAnalysis_computeCentrality(cellchat)

# Check Neutrophils
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Activated Neutrophils", signaling.exclude = "MIF")

# Check Stem cells
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Hematopoietic Stem/Progenitor Cellss", signaling.exclude = "MIF")

# Combine the plots
patchwork::wrap_plots(plots = list(gg1, gg2))


#3. Identify altered signaling with distinct interaction strength
##3.1Compare the overall information flow of each signaling pathway or ligand-receptor pair
gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE)

gg1 + gg2

##3.2 Compare outgoing (or incoming) signaling patterns associated with each cell population

library(ComplexHeatmap)
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

#4.Identify the up-gulated and down-regulated signaling ligand-receptor pairs

##4.1.Identify dysfunctional signaling by comparing the communication probabities
#> Comparing communications on a merged object
#> > print(HNifc_cell_types)
#
#[1] "Mature Neutrophils"                    
#[2] "Activated Neutrophils"                 
#[3] "Hematopoietic Stem/Progenitor Cells"   
#[4] "Immature Neutrophils"                  
#[5] "pre B Cells"                           
#[6] "Megakaryocytes/Platelets"              
#[7] "Dendritic Cells"                       
#[8] "Fibroblasts or Mesenchymal Cells"      
#[9] "Erythroid Cells/Erythrocyte Precursors"
#[10] "Pro-B Cells"                           
#[11] "Macrophages Monocytes"                 
#[12] "Antimicrobial Neutrophils"             
#[13] "Undefined"                             
#[14] "Mast Cells/Progenitors"                
#[15] "Natural Killer Cells"                  
#[16] "Plasma Cells"       

# activated Neutrophils as senders 
ptm = Sys.time()

netVisual_bubble(cellchat, sources.use = 2, targets.use = c(1:20),  comparison = c(1, 2), angle.x = 45)

#anti-microbial neutrophils as senders
ptm = Sys.time()

netVisual_bubble(cellchat, sources.use = 12, targets.use = c(1:20),  comparison = c(1, 2), angle.x = 45)

#mature neutrophils as senders
ptm = Sys.time()

netVisual_bubble(cellchat, sources.use = 1, targets.use = c(1:20),  comparison = c(1, 2), angle.x = 45)

#Immature neutrophils as senders
ptm = Sys.time()

netVisual_bubble(cellchat, sources.use = 4, targets.use = c(1:20),  comparison = c(1, 2), angle.x = 45)


# Stem cells as senders 
ptm = Sys.time()

netVisual_bubble(cellchat, sources.use = 3, targets.use = c(1:20),  comparison = c(1, 2), angle.x = 45)

# Monocytes as senders 
ptm = Sys.time()

netVisual_bubble(cellchat, sources.use = 11, targets.use = c(1:20),  comparison = c(1, 2), angle.x = 45)



##More specifically for the upregulated and the downregulated 
gg1 <- netVisual_bubble(cellchat, sources.use = 3, targets.use = c(1:16),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in HNifc group", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = 3, targets.use = c(1:16),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in HNifc group", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2


###Visually compare cell-cell communication with three different plots
##Compare the signalling gene expression distribution between different datasets 
#Violin plot
cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("HNifc", "Control")) # set factor level
plotGeneExpression(cellchat, signaling = "CXCL", split.by = "datasets", colors.ggplot = T, type = "violin")

#Circle plots
pathways.show <- c("CXCL") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

#heatmap
pathways.show <- c("CXCL") 
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}

ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))

#chord diagram 
# Chord diagram
pathways.show <- c("CXCL") 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
}


save.image(file = "scRNAseq#2CellChat.RData")













# Compare communication networks

#########The code below is for two datasets with different cell compositions 

library(CellChat)
library(patchwork)

##The goal is to make two datasets (slightly different) contain the same cell composition 
##For vastly different cell groups, some functions do not work # Check and print cell types in HNifc_cellchat
##The idea is to lift up the ones with less compositions to match the one with more 

# Check and print cell types in HNifc_cellchat
hnifc_cell_types <- levels(HNifc_cellchat@idents)
print("Cell types in HNifc_cellchat:")
print(hnifc_cell_types)

# Check and print cell types in Control_cellchat
control_cell_types <- levels(control_cellchat@idents)
print("Cell types in Control_cellchat:")
print(control_cell_types)

# Define the combined cell types from both objects
group.new <- union(levels(HNifc_cellchat@idents), levels(control_cellchat@idents))

# Lift up the cellchat objects to match the combined cell types
HNifc_cellchat <- liftCellChat(HNifc_cellchat, group.new)
control_cellchat <- liftCellChat(control_cellchat, group.new)

# Merge the CellChat objects
object.list <- list(HNifc = HNifc_cellchat, Control = control_cellchat)
cellchat <- mergeCellChat(object.list, add.names = c("HNifc", "Control"), cell.prefix = TRUE)

# Print the merged CellChat object
print(cellchat)

# Check and print cell types in HNifc_cellchat
hnifc_cell_types <- levels(HNifc_cellchat@idents)
print("Cell types in HNifc_cellchat:")
print(hnifc_cell_types)

# Check and print cell types in Control_cellchat
control_cell_types <- levels(control_cellchat@idents)
print("Cell types in Control_cellchat:")
print(control_cell_types)

# Check if the cell types are the same
if (setequal(hnifc_cell_types, control_cell_types)) {
  print("Both CellChat objects have the same cell composition.")
} else {
  print("The CellChat objects have different cell compositions.")
}



#Once it's done, move on to visualization 
#Comparative visualization and analysis of cell-cell communication using the lifted objects

# Circle plot
pathways.show <- c("CXCL") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}


#########The code below is for two datasets with the same cell composition

#Compare the total number of interactions and interaction strength
ptm = Sys.time()
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

##Compare the number of interactions and interaction strength among different cell populations
#Form 1: circle plot
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

#Form2: heatmap
gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2

#Form3: circle plot (for the number of interactions and interaction weight)
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

# Get the maximum weight for scaling the interaction weights
weight.max <- getMaxWeight(object.list, attribute = c("idents", "weight"))

# Set up the plotting area
par(mfrow = c(1, 2), xpd = TRUE)

# Loop through each CellChat object to visualize
for (i in 1:length(object.list)) {
  netVisual_circle(
    object.list[[i]]@net$weight,  # Use 'weight' instead of 'count'
    weight.scale = TRUE, 
    label.edge = FALSE, 
    edge.weight.max = weight.max[2], 
    edge.width.max = 12, 
    title.name = paste0("Interaction weights - ", names(object.list)[i])
  )
}


#Form4: circle plot (only shown the interactions between certain cell types)
# Define the cell types of interest
group.cellType <- c(rep("Stem cells", 1), rep("Neutrophils", 1), rep("Monocytes", 1))

# Set the factor levels to match the cell types you're interested in
group.cellType <- factor(group.cellType, levels = c("Stem cells", "Neutrophils", "Monocytes"))

# Merge interactions for the specified cell types
object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})

# Merge the CellChat objects
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

# Get the maximum weight for scaling the interaction weights in the merged data
weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents", "count", "count.merged"))

# Set up the plotting area
par(mfrow = c(1, 2), xpd = TRUE)

# Loop through each CellChat object to visualize
for (i in 1:length(object.list)) {
  netVisual_circle(
    object.list[[i]]@net$count.merged, 
    weight.scale = TRUE, 
    label.edge = TRUE, 
    edge.weight.max = weight.max[3], 
    edge.width.max = 12, 
    title.name = paste0("Number of interactions - ", names(object.list)[i])
  )
}

# Plot interaction strength
netVisual_circle(
  object.list[[i]]@net$weight.merged, 
  weight.scale = TRUE, 
  label.edge = TRUE, 
  edge.weight.max = weight.max.strength[3], 
  edge.width.max = 12, 
  title.name = paste0("Interaction strength - ", names(object.list)[i])
)


##Compare the major sources and targets in a 2D space

#1. identify cell populations with significant changes in sending or receiving signals
#scatter plot
# Compute the network centrality scores for each object in the list
object.list <- lapply(object.list, netAnalysis_computeCentrality)

# Proceed with the scatter plot
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count) - diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # Control the dot size in the different datasets
gg <- list()

for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}

# Plot the results
gridExtra::grid.arrange(grobs = gg, ncol = 2)


#2. Identify the signalling changes of specific cell populations
##(the code below does not work...)
# Compute centrality scores for the CellChat object
cellchat <- netAnalysis_computeCentrality(cellchat)

# Check Neutrophils
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Neutrophils", signaling.exclude = "MIF")

# Check Stem cells
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Stem cells", signaling.exclude = "MIF")

# Combine the plots
patchwork::wrap_plots(plots = list(gg1, gg2))


#3. Identify altered signaling with distinct interaction strength
##3.1Compare the overall information flow of each signaling pathway or ligand-receptor pair
gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE)

gg1 + gg2

##3.2 Compare outgoing (or incoming) signaling patterns associated with each cell population

library(ComplexHeatmap)
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

#4.Identify the up-gulated and down-regulated signaling ligand-receptor pairs

##4.1.Identify dysfunctional signaling by comparing the communication probabities
#> Comparing communications on a merged object
#> > print(hnifc_cell_types)
#[1] "B cells"           "B cells, pro"      "Basophils"        
#[4] "DC"                "Epithelial cells"  "ILC"              
#[7] "Macrophages"       "Monocytes"         "Neutrophils"      
#[10] "NK cells"          "NKT"               "Stem cells"       
#[13] "Stromal cells"     "T cells"           "Endothelial cells"
#[16] "Eosinophils"       "Fibroblasts"       "Microglia"        
#[19] "Tgd"     

# Neutrophils as senders 
ptm = Sys.time()

netVisual_bubble(cellchat, sources.use = 9, targets.use = c(7:16),  comparison = c(1, 2), angle.x = 45)

# Stem cells as senders 
ptm = Sys.time()

netVisual_bubble(cellchat, sources.use = 12, targets.use = c(7:16),  comparison = c(1, 2), angle.x = 45)

# Monocytes as senders 
ptm = Sys.time()

netVisual_bubble(cellchat, sources.use = 8, targets.use = c(7:16),  comparison = c(1, 2), angle.x = 45)



##More specifically for the upregulated and the downregulated 
gg1 <- netVisual_bubble(cellchat, sources.use = 12, targets.use = c(7:16),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in LS", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = 12, targets.use = c(7:16),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in LS", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2


###Visually compare cell-cell communication with three different plots
##Compare the signalling gene expression distribution between different datasets 
#Violin plot
cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("HNifc", "Control")) # set factor level
plotGeneExpression(cellchat, signaling = "CXCL", split.by = "datasets", colors.ggplot = T, type = "violin")

#Circle plots
pathways.show <- c("CXCL") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

#heatmap
pathways.show <- c("CXCL") 
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}

ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))

#chord diagram 
# Chord diagram
pathways.show <- c("CXCL") 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
}


############This code is run after mannual annotation!!

#Set up control and HNifc cellchat

library(CellChat)
library(Seurat)
library(patchwork)

# Convert Seurat to CellChat
cellchat_HNifc <- createCellChat(object = HNifc_seurat_object, group.by = "ident")
cellchat_Control <- createCellChat(object = control_seurat_object, group.by = "ident")

# Set the database for communication
CellChatDB <- CellChatDB.mouse # or CellChatDB.mouse if using mouse data
cellchat_HNifc@DB <- CellChatDB
cellchat_Control@DB <- CellChatDB

# Subset the expression data of signaling genes
cellchat_HNifc <- subsetData(cellchat_HNifc)
cellchat_Control <- subsetData(cellchat_Control)

# Identify overexpressed genes and interactions
cellchat_HNifc <- identifyOverExpressedGenes(cellchat_HNifc)
cellchat_Control <- identifyOverExpressedGenes(cellchat_Control)

cellchat_HNifc <- identifyOverExpressedInteractions(cellchat_HNifc)
cellchat_Control <- identifyOverExpressedInteractions(cellchat_Control)

#Compute commun prob
cellchat_HNifc <- computeCommunProb(cellchat_HNifc)
cellchat_Control <- computeCommunProb(cellchat_Control)

# Filter out low probability interactions
cellchat_HNifc <- filterCommunication(cellchat_HNifc, min.cells = 10)
cellchat_Control <- filterCommunication(cellchat_Control, min.cells = 10)

#Infer cellular communication networks 
cellchat_HNifc <- computeCommunProbPathway(cellchat_HNifc)
cellchat_Control <- computeCommunProbPathway(cellchat_Control)

cellchat_HNifc <- aggregateNet(cellchat_HNifc)
cellchat_Control <- aggregateNet(cellchat_Control)


# Compare communication networks between HNifc and Control CellChat

#########The code below is for two datasets with different cell compositions 

library(CellChat)
library(patchwork)

##The goal is to make two datasets (slightly different) contain the same cell composition 
##For vastly different cell groups, some functions do not work # Check and print cell types in HNifc_cellchat
##The idea is to lift up the ones with less compositions to match the one with more 

# Check and print cell types in HNifc_cellchat
HNifc_cell_types <- levels(cellchat_HNifc@idents)
print("Cell types in HNifc_cellchat:")
print(HNifc_cell_types)

# Check and print cell types in Control_cellchat
control_cell_types <- levels(cellchat_Control@idents)
print("Cell types in Control_cellchat:")
print(control_cell_types)

# Define the combined cell types from both objects
group.new <- union(levels(cellchat_HNifc@idents), levels(cellchat_Control@idents))

# Lift up the cellchat objects to match the combined cell types
cellchat_HNifc <- liftCellChat(cellchat_HNifc, group.new)
cellchat_Control <- liftCellChat(cellchat_Control, group.new)

# Merge the CellChat objects
object.list <- list(HNifc = cellchat_HNifc, Control = cellchat_Control)
cellchat <- mergeCellChat(object.list, add.names = c("HNifc", "Control"), cell.prefix = TRUE)

# Print the merged CellChat object
print(cellchat)

# Check and print cell types in HNifc_cellchat
HNifc_cell_types <- levels(cellchat_HNifc@idents)
print("Cell types in HNifc_cellchat:")
print(HNifc_cell_types)

# Check and print cell types in Control_cellchat
control_cell_types <- levels(cellchat_Control@idents)
print("Cell types in Control_cellchat:")
print(control_cell_types)

# Check if the cell types are the same
if (setequal(HNifc_cell_types, control_cell_types)) {
  print("Both CellChat objects have the same cell composition.")
} else {
  print("The CellChat objects have different cell compositions.")
}


