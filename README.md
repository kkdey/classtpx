# classtpx

A Repository for Semi-supervised Topic models, applicable in set ups where for some samples we know which subgroups they come from and we wish to use these subgroups to drive the clusters for the other samples for which no information is available. 

The application of this semi-supervised approach lies mainly in Biological applications. For example, we have RNA-seq data on single cells. Some single cells have been FACS sorted and we know which cell type or cell cycle phase these cells belong to. However, this information is not known for all the cells. We can use the information from the FACS sorted cells to drive the overall clustering and we can ensure that to a great extent, these clusters will actually reflect true cell-type or cell-phase clusters. 

The code in this repository is flexible and allows for fitting the model where no information is available on any of the samples (no cell is FACS sorted) and that would correspond to the unsupervised approach as in `maptpx` due to Matt Taddy (2012).

At the other end, we may have all cells FACS sorted in which case we have a classifier and this we can use for `predict.class_topics` for classification. 
