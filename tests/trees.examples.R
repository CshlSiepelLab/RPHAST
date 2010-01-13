require("rphast")

#' read.tree
cat(c("((hg18:0.142679,(mm9:0.083220,rn4:0.090564):0.269385):0.020666,canFam2:0.193569);",
      "(human, (mouse, rat));",
      sep="\n"), file="test.nh")
read.tree("test.nh")
unlink("test.nh")


#' tree.numnodes
tree.numnodes(c("((hg18:0.142679,(mm9:0.083220,rn4:0.090564):0.269385):
                0.020666,canFam2:0.193569);",
                "(human, (mouse, rat));"))
       
#' tree.prune
trees <- c("((hg18, panTro2), mm9);",
           "((hg18:0.142679,(mm9:0.083220,rn4:0.090564):0.269385):
                0.020666,canFam2:0.193569);")
tree.prune(trees, c("panTro2", "mm9"), all.but=TRUE)
tree.prune(trees, "hg18", all.but=FALSE)


#' tree.name.ancestors
trees <- c("((hg18, panTro2), mm9);",
           "((hg18:0.142679,(mm9:0.083220,rn4:0.090564):0.269385):
                0.020666,canFam2:0.193569);")
tree.name.ancestors(trees)


#' tree.subtree
trees <- c("((hg18, panTro2), mm9);",
           "((hg18:0.142679,(mm9:0.083220,rn4:0.090564):0.269385)
            :0.020666,canFam2:0.193569);")
trees <- tree.name.ancestors(trees)
tree.subtree(trees, c("hg18-panTro2", "mm9-rn4"))

#' tree.scale
q("no")
require("rphast")
trees <- c("((hg18:1.0, panTro2:2.0):3.0, mm9:4.0);",
           "((hg18:0.142679,(mm9:0.083220,rn4:0.090564):0.269385):
                0.020666,canFam2:0.193569);")
tree.scale(trees, 0.5)
tree.scale(trees, c(0.5, 2.0))
trees <- tree.name.ancestors(trees)
tree.scale(trees, 0.5, c("hg18-panTro2", "hg18-mm9"))
      
#' tree.rename
trees <- c("((hg18:1.0, panTro2:2.0):3.0, mm9:4.0);",
           "((hg18:0.142679,(mm9:0.083220,rn4:0.090564):0.269385):
                0.020666,canFam2:0.193569);")
tree.rename(trees,
            old.names=c("hg18", "panTro2", "mm9", "rn4", "canFam2"),
            new.names=c("human", "chimp", "mouse", "rat", "dog"))
