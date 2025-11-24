
#### Load the program ####
library(ape)

#### Set the working directory to the folder in which the tree is ####
setwd("path/to/you/working_directory")

#### 1. Read the Newick tree file #####
my_tree <- read.tree("your_newick.treefile")

#### 2. Define all your outgroup species ####
# Create a character vector with all the species names that form your outgroup.

outgroup_species <- c("Sp_0") ### If there were more outgroup species, you should add them like this -> c("Sp_0", "Sp_X", ...)

#### 3. Root the tree with the multi-species outgroup ####
# This finds the common ancestor of all species in 'outgroup_species' and roots the tree there.

my_tree <- root(my_tree, outgroup = outgroup_species, resolve.root = TRUE)

# (Optional) You can still define your ingroup this way if it's useful for later steps
ingroup_species <- my_tree$tip.label[!my_tree$tip.label %in% outgroup_species]

#### 4. Define the root node for calibration ####
# In a rooted 'phylo' object, the root node number is always Ntip + 1.

root_node_number <- length(my_tree$tip.label) + 1

#### 5. Define a vector containing the nodes you wish to date (not in a romantic way). A.k.a. the Most Recent Common Ancestors ####
node <- c(
  root_node_number, ### MRCA between ingroup_species and outgroup_species
  getMRCA(my_tree, tip = c("Sp_A","Sp_B")), ### MRCA of Sp_A & Sp_B
  getMRCA(my_tree, tip = c("Sp_C","Sp_D")), ### MRCA of Sp_C & Sp_D
  getMRCA(my_tree, tip = c("Sp_A","Sp_B","Sp_C","Sp_D")) ### MRCA of all ingroup_species
)

#### 6. Define the maximum and minimum age estimates for the nodes established in 5. ####
# In these variables we will define vectors containing the max-min age estimates, following the order established in 5. 

# Define minimum ages in millions of years (Mya)
age.min <- c(
  130,
  6.3,
  4.25,
  23
)


# Define maximum ages in Mya
age.max <- c(
  160,
  13.3,
  8.87,
  36.5
)

# Define typology of age constraints Soft (TRUE) vs Hard (FALSE) calibrations 
soft.bound <- c(
  FALSE,
  FALSE,
  FALSE,
  FALSE
)  

#### 7. We create a table containing all calibration information ####
mycalibration <- data.frame(node,age.min,age.max,soft.bound)

#### 8. We calibrate the tree using `chronos` #####
mytimetree <- chronos(my_tree, lambda = 1, model = "relaxed", calibration = mycalibration, control = chronos.control())

#### 9. Save the calibrated tree in Nexus and Newick formats ####
write.tree(mytimetree, file="your_calibrated_tree.tree")     # Newick format
write.tree(mytimetree, file="your_calibrated_tree.nex")      # Nexus format

#### 10. Save the timetree in pdf format with time axis.

pdf("your_calibrated_tree.pdf", width = 8, height = 6) 
plot(mytimetree)

axisPhylo(side = 1,          # Place axis at the bottom (1 = X-axis)
          backward = TRUE,   # Present (0) on the right, past on the left
          las = 1,           # Horizontal axis labels
          col = "black")     # Axis color

# Add custom ticks (e.g., every 5 Mya)
root_age <- max(node.depth.edgelength(mytimetree))
custom_ticks <- seq(0, root_age, by = 5)  # Change '10' to '5' for ticks every 5 Mya

# Add minor ticks (no labels)
axis(1, at = root_age - custom_ticks, labels = FALSE, tcl = -0.3, col.ticks = "gray50")

# Add labeled major ticks (every 20 Mya, optional)
major_ticks <- seq(0, root_age, by = 10)
axis(1, at = root_age - major_ticks, labels = major_ticks, las = 1, col = "black")


# Add axis label
mtext("Time (Million Years Ago)", 
      side = 1,              # Below the plot (X-axis)
      line = 2.5,            # Adjust position
      cex = 1)               # Label size
dev.off()

