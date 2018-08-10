library(Seurat)
library(Matrix)


make_seurat = function(atac_matrix, metadata=NULL, pca_coords=NULL, tsne_coords=NULL) {
               if (is.null(pca_coords) && is.null(tsne_coords)) {
			stop('One of pca_coords or tsne_coords (or both) must be set to make the fake seurat object.')
	       }

               if (!is.null(metadata)) {
                   seurat_obj = CreateSeuratObject(atac_matrix, meta.data=metadata, project='DEFAULT', min.cells=0, min.genes=0)
	       } else {
		   seurat_obj = CreateSeuratObject(atac_matrix, project='DEFAULT', min.cells=0, min.genes=0)
	       }

	       if (!is.null(pca_coords)) {
			pca_coords = as.matrix(pca_coords)
               		seurat_obj <- SetDimReduction(object = seurat_obj, reduction.type = "pca", slot = "cell.embeddings", new.data = pca_coords)
               		seurat_obj <- SetDimReduction(object = seurat_obj, reduction.type = "pca", slot = "key", new.data = "pca")
	       }

	       if (!is.null(tsne_coords)) {
			tsne_coords = as.matrix(tsne_coords)
                        seurat_obj = CreateSeuratObject(atac_matrix, meta.data=metadata, project='DEFAULT', min.cells=0, min.genes=0)
                        seurat_obj <- SetDimReduction(object = seurat_obj, reduction.type = "tsne", slot = "cell.embeddings", new.data = tsne_coords)
                        seurat_obj <- SetDimReduction(object = seurat_obj, reduction.type = "tsne", slot = "key", new.data = "tsne")
	       }

               return(seurat_obj)
}


# Then you can cluster using FindClusters in Seurat
