require("rphast")

#' read.newick.tree
cat(c("((hg18:0.142679,(mm9:0.083220,rn4:0.090564):0.269385):0.020666,canFam2:0.193569);",
      "(human, (mouse, rat));",
      sep="\n"), file="test.nh")
read.newick.tree("test.nh")
unlink("test.nh")


#' numnodes.tree
numnodes.tree(c("((hg18:0.142679,(mm9:0.083220,rn4:0.090564):0.269385):
                0.020666,canFam2:0.193569);",
                "(human, (mouse, rat));"))
       
#' prune.tree
trees <- c("((hg18, panTro2), mm9);",
           "((hg18:0.142679,(mm9:0.083220,rn4:0.090564):0.269385):
                0.020666,canFam2:0.193569);")
prune.tree(trees, c("panTro2", "mm9"), all.but=TRUE)
prune.tree(trees, "hg18", all.but=FALSE)


#' name.ancestors
trees <- c("((hg18, panTro2), mm9);",
           "((hg18:0.142679,(mm9:0.083220,rn4:0.090564):0.269385):
                0.020666,canFam2:0.193569);")
name.ancestors(trees)


#' subtree
trees <- c("((hg18, panTro2), mm9);",
           "((hg18:0.142679,(mm9:0.083220,rn4:0.090564):0.269385)
            :0.020666,canFam2:0.193569);")
trees <- name.ancestors(trees)
subtree(trees, c("hg18-panTro2", "mm9-rn4"))

#' rescale.tree
trees <- c("((hg18:1.0, panTro2:2.0):3.0, mm9:4.0);",
           "((hg18:0.142679,(mm9:0.083220,rn4:0.090564):0.269385):
                0.020666,canFam2:0.193569);")
rescale.tree(trees, 0.5)
rescale.tree(trees, c(0.5, 2.0))
trees <- name.ancestors(trees)
rescale.tree(trees, 0.5, c("hg18-panTro2", "hg18-mm9"))
      
#' rename.tree
trees <- c("((hg18:1.0, panTro2:2.0):3.0, mm9:4.0);",
           "((hg18:0.142679,(mm9:0.083220,rn4:0.090564):0.269385):
                0.020666,canFam2:0.193569);")
rename.tree(trees,
            old.names=c("hg18", "panTro2", "mm9", "rn4", "canFam2"),
            new.names=c("human", "chimp", "mouse", "rat", "dog"))

#' fix.semicolon.tree
str <- c("213", "345")
fix.semicolon.tree(str)
str <- c("213;", "345;")
fix.semicolon.tree(str)
str <- c("213", "345;")
fix.semicolon.tree(str)


#' label.branches
trees <- c("((hg18:1.0, panTro2:2.0)hg18-panTro2:3.0, mm9:4.0);",
           "((hg18:0.142679,(mm9:0.083220,rn4:0.090564)mm9-rn4:
             0.269385)hg18-rn4:0.020666,canFam2:0.193569);")
label.branches(trees, c("hg18", "mm9"), "humanAndMouse")


#' label.subtree
trees <- c("((hg18:1.0, panTro2:2.0)hg18-panTro2:3.0, mm9:4.0);",
           "(((hg18:0.01, panTro2:0.01)hg18-panTro2:0.07,
              (mm9:0.083220,rn4:0.090564)mm9-rn4:
             0.269385)hg18-rn4:0.020666,canFam2:0.193569);")
label.subtree(trees, "hg18-panTro2", "human-chimp", include.leading=FALSE)
label.subtree(trees, "hg18-panTro2", "human-chimp", include.leading=TRUE)


#' summary.tree
tree <- "(((hg18:0.01, panTro2:0.01)hg18-panTro2:0.07,
              (mm9:0.083220,rn4:0.090564)mm9-rn4:
             0.269385)hg18-rn4:0.020666,canFam2:0.193569);"
summary.tree(tree)
summary.tree(label.subtree(tree, "mm9-rn4", "rodent", include.leading=TRUE))


rm(list = ls())
gc()
