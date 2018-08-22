library(Matrix)
library(Rtsne)
library(irlba)

atac_dim_reduction = function(atac_matrix, site_frequency_threshold=0.03) {
               num_cells_ncounted = Matrix::rowSums(atac_matrix)
               threshold = ncol(atac_matrix) * site_frequency_threshold

               ncounts = atac_matrix[num_cells_ncounted >= threshold,]

               ## Normalize the data with TF-IDF
               nfreqs = t(t(ncounts) / Matrix::colSums(ncounts))
               tf_idf_counts = nfreqs * log(1 + ncol(ncounts) / Matrix::rowSums(ncounts))

               ## Do SVD
               set.seed(0)
               SVDtsne = irlba(tf_idf_counts, 50, 50, maxit=1000)
               d_diagtsne = matrix(0, nrow=length(SVDtsne$d), ncol=length(SVDtsne$d))
               diag(d_diagtsne) = SVDtsne$d
               SVDtsne_vd = t(d_diagtsne %*% t(SVDtsne$v))
               rownames(SVDtsne_vd) = colnames(atac_matrix)
               colnames(SVDtsne_vd) = paste0('pca_', 1:ncol(SVDtsne_vd))

               ## Run TSNE to 2 dimensions
               tsnetfidf = Rtsne(SVDtsne_vd, pca=F, perplexity=30, max_iter=5000)

               tsne_coords = as.data.frame(tsnetfidf$Y)
               colnames(tsne_coords) = c('tsne_1', 'tsne_2')
               rownames(tsne_coords) = colnames(ncounts)

               return(list('tsne_coords'=tsne_coords, 'pca_coords'=as.data.frame(SVDtsne_vd)))
}
