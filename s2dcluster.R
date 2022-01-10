#### clustering functions ####
#### sequencial UMAP dim reduction and density based clustering ####
.umap_dbs_clust<-function(data, maxevts, maxlvl, minpts, clust_options, idxs = 1:(nrow(data))>0, lvl = 1, clust = "cluster", multi_thread = T, fast_sgd = TRUE, ret_model = F, myqueue = NULL, verbose = verbose) {
  knn_needed<-FALSE
  if (!is.matrix(data[idxs,])) {
    if (ret_model){
      return(list(cluster = rep(clust, 1), model = list(cluster = clust, data_columns = length(data), maxevts = maxevts, maxlvl = maxlvl, minpts = minpts, clust_options = clust_options)))
    } else return(rep(clust, 1))
  }
  if (nrow(data[idxs,]) <= minpts) {
    if (ret_model){
      return(list(cluster = rep(clust, nrow(data[idxs,])), model = list(cluster = clust, data_columns = ncol(data), maxevts = maxevts, maxlvl = maxlvl, minpts = minpts, clust_options = clust_options)))
    } else return(rep(clust, nrow(data[idxs,])))
  }
  if (nrow(data[idxs,])>maxevts) {
    datasample<-sample(1:nrow(data[idxs,]), maxevts, replace = FALSE)
    knn_needed<-T
  } else datasample<-1:nrow(data[idxs,])
  method = clust_options$method
  mineps = clust_options$mineps
  mindens = clust_options$mindens
  bw = clust_options$bw
  nbins = clust_options$nbins
  mvpratio = clust_options$mvpratio
  # UMAP modeling
  umap_model<-uwot::umap(data[idxs,][datasample,], scale = FALSE, n_threads = ifelse(multi_thread, RcppParallel::defaultNumThreads(), 1), fast_sgd = fast_sgd, ret_model = T, verbose = verbose)
  embdata<-umap_model$embedding
  # clustering UMAP embedding
  clusters<-switch (method,
                    "dbscan" = {
                      num_el<-nrow(embdata)
                      vdist<-sort(as.vector(FNN::knn.dist(embdata, minpts)))
                      x1<-length(vdist)-num_el
                      eps<-max(mineps, vdist[x1])
                      res<-dbscan::dbscan(embdata, eps = eps, minPts = minpts)
                      if (any(res$cluster == 0)) clusters<-res$cluster + 1L else clusters<-res$cluster
                      clusters
                    },
                    "kdbscan" = {
                      kdbscan(embdata, minpts = max(10,minpts*(nrow(embdata)/nrow(data[idxs,]))), mindens = mindens, ret_model = T)
                    },
                    "hdbscan" = {
                      res<-dbscan::hdbscan(embdata, minPts = minpts)
                      if (any(res$cluster == 0)) clusters<-res$cluster + 1L else clusters<-res$cluster
                      clusters
                    },
                    "sdbscan" = {
                      sdbscan(embdata, minpts = max(10,minpts*(nrow(embdata)/nrow(data[idxs,]))), bw = bw, nbins = nbins, mvpratio = mvpratio, ret_model = T)
                    }
  )
  if (inherits(clusters, "s2dcluster")) {
    cl_model<-clusters
    clusters<-clusters$cluster
  } else cl_model<-NULL
  # cluster not sampled data
  if (length(unique(clusters)) > 1) {
    if (knn_needed) {
      new_embdata<-uwot::umap_transform(data[idxs,][-datasample,], umap_model, n_threads = ifelse(multi_thread, RcppParallel::defaultNumThreads(), 1), verbose = verbose)
      if (is.null(cl_model)){
        knn<-FNN::get.knnx(embdata, new_embdata, k=1)
        newclusters<-clusters[knn$nn.index]
      } else {
        newclusters<-predict.s2dcluster(cl_model, new_embdata)
      }
      all_clusters<-integer(nrow(data[idxs,]))
      all_clusters[datasample]<-clusters
      all_clusters[-datasample]<-newclusters
      clusters<-all_clusters
    }
    new_clusters<-sort(unique(clusters))
    new_clusters<-paste0(clust, "_", new_clusters)
    clusters<-new_clusters[clusters]
  } else {
    clusters<-(rep(clust, nrow(data[idxs,])))
    new_clusters<-clust
  }
  newmodels<-list()
  if ((maxlvl > lvl) && length(new_clusters) > 1) {
    if (!is.null(myqueue)) myqueue$producer$fireAssignReactive("umap_dbs_msg", paste0("Found ", length(new_clusters), " new clusters on ", clust, " (level ", lvl, "). Trying deeper levels."))
    for (c in new_clusters) {
      nidxs<-idxs
      nidxs[idxs]<-nidxs[idxs] & (clusters == c)
      c<-c
      if (multi_thread) {
        futureAssign(paste0("nclust", c), .umap_dbs_clust(data, maxevts, maxlvl, minpts, clust_options, idxs = nidxs, lvl = lvl+1, clust = c, multi_thread = multi_thread, fast_sgd = fast_sgd, ret_model = ret_model, myqueue = myqueue, verbose = verbose), seed = T)
      } else {
        assign(paste0("nclust", c), .umap_dbs_clust(data, maxevts, maxlvl, minpts, clust_options, idxs = nidxs, lvl = lvl+1, clust = c, multi_thread = multi_thread, fast_sgd = fast_sgd, ret_model = ret_model, myqueue = myqueue, verbose = verbose))
      }
    }
    for (c in new_clusters) {
      if (ret_model) {
        tempclust<-get(paste0("nclust", c))
        clusters[clusters == c]<-tempclust$cluster
        if (length(unique(tempclust$cluster)) > 1) newmodels[[c]]<-tempclust$model
      } else clusters[clusters == c]<-get(paste0("nclust", c))
    }
  } else if (!is.null(myqueue)) myqueue$producer$fireAssignReactive("umap_dbs_msg", paste0("Found only 1 cluster on ", clust, " (level ", lvl, ") or no deeper level allowed. Skipping deeper levels on this cluster."))
  if (lvl == 1) if (!is.null(myqueue)) myqueue$producer$fireAssignReactive("umap_dbs_msg", paste0("Found a total of ", length(unique(clusters)), " clusters."))
  if (ret_model) {
    umap_file<-tempfile()
    uwot::save_uwot(umap_model, umap_file)
    thismodel<-list(cluster = clust, data_columns = ncol(data), maxevts = maxevts, maxlvl = maxlvl, minpts = minpts, clust_options = clust_options)
    if (length(unique(clusters)) > 1) thismodel<-c(thismodel, list(umap_file = umap_file, umap_sample = datasample, s2dcluster_model = cl_model))
    if (length(newmodels) > 0) thismodel<-c(thismodel, newmodels)
    return(list(cluster = clusters, model = thismodel))
  } else return(clusters)
}

umap_dbs_clust<-function(data, maxevts, maxlvl, minpts, clust_options, multi_thread = T, fast_sgd = TRUE, ret_model = F, myqueue = NULL, verbose = F) {
  if (ret_model) {
    if (!(clust_options$method %in% c("sdbscan", "kdbscan"))) stop("ret_model needs sdbscan or kdbscan clustering method!", call. = F)
  }
  clusters<-.umap_dbs_clust(data = data, maxevts = maxevts, maxlvl = maxlvl, minpts = minpts, clust_options = clust_options, multi_thread = multi_thread, fast_sgd = fast_sgd, ret_model = ret_model, myqueue = myqueue, verbose = verbose)
  if (ret_model) {
    clusters_model<-clusters$model
    clusters<-clusters$cluster
  }
  ctable<-as.data.frame(table(clusters))
  ctable<-ctable[order(ctable$Freq, decreasing = T),]
  new_clusters<-c(1:nrow(ctable))
  names(new_clusters)<-ctable$clusters
  clusters<-unname(new_clusters[clusters])
  if (ret_model) {
    res<-list(cluster = clusters, clusters=new_clusters, model = clusters_model)
    class(res)<-"UMAP_s2dcluster"
    return(res)
  } else return(clusters)
}

.UMAP_s2dcluster_newdata<-function(x, model, multi_thread = T, verbose = F) {
  umap_model<-uwot::load_uwot(model$umap_file)
  embdata<-uwot::umap_transform(x, umap_model, n_threads = ifelse(multi_thread, RcppParallel::defaultNumThreads(), 1), verbose = verbose)
  s2dmodel<-model$s2dcluster_model
  clust<-predict.s2dcluster(s2dmodel, embdata)
  clust<-paste0(model$cluster, "_", clust)
  new_clusters<-unique(clust)
  for (c in new_clusters) {
    c<-c
    if (is.null(model[[c]])) next
    idxs<-clust == c
    nclust<-.UMAP_s2dcluster_newdata(x[idxs,], model[[c]], multi_thread = multi_thread, verbose = verbose)
    clust[idxs]<-nclust
  }
  return(clust)
}

predict.UMAP_s2dcluster<-function(res, newdata, multi_thread = T, verbose = F) {
  if (!inherits(res, "UMAP_s2dcluster")) stop("res must be a UMAP_s2dcluster result object!", call. = F)
  if (!is.matrix(newdata) || ncol(newdata) != res$model$data_columns || nrow(newdata) < 1) stop(paste0("newdata must be a matrix with ", res$model$data_columns, " columns and >0 rows!"), call. = F)
  cluster<-c()
  if (length(res$clusters) > 1) {
    cluster<-.UMAP_s2dcluster_newdata(newdata, res$model, multi_thread = multi_thread, verbose = verbose)
  } else cluster<-rep("cluster", nrow(newdata))
  clusters<-unname(res$clusters[cluster])
  return(clusters)
}

.save_model<-function(model, modeldir) {
  modelfile<-file.path(modeldir, basename(model$umap_file))
  file.copy(model$umap_file, modelfile, overwrite = T)
  model$umap_file<-modelfile
  for (i in paste0(model$cluster, "_", model$s2dcluster_model$clusters)) {
    if (!is.null(model[[i]])) {
      model[[i]]<-.save_model(model[[i]], modeldir)
    }
  }
  return(model)
}

save_model<-function(res, file) {
  if (!inherits(res, "UMAP_s2dcluster")) stop("res must be a UMAP_s2dcluster result object!", call. = F)
  writeaccess<-tryCatch({
    suppressWarnings(write.table(1, file))
  }, error=function(e) FALSE)
  if (!is.null(writeaccess)) {
    stop("user must have write access rights to file!", call. = F)
  } else {
    file.remove(file)
  }
  dir<-dirname(file)
  filename<-sub("([^.]+)\\.[[:alnum:]]+$", "\\1", basename(file))
  modeldir<-paste0(filename, "_models")
  olddir<-getwd()
  setwd(dir)
  dir.create(modeldir)
  res$model<-.save_model(res$model, modeldir)
  saveRDS(res, file)
  setwd(olddir)
  return(invisible(res))
}

### CUML version
.cuml_umap_dbs_clust<-function(data, maxevts, maxlvl, minpts, clust_options, idxs = 1:(nrow(data))>0, lvl = 1, clust = "cluster", multi_thread = T, ret_model = F, myqueue = NULL, verbose = verbose, cuml = NULL) {
  knn_needed<-FALSE
  if (!is.matrix(data[idxs,])) {
    if (ret_model){
      return(list(cluster = rep(clust, 1), model = list(cluster = clust, data_columns = length(data), maxevts = maxevts, maxlvl = maxlvl, minpts = minpts, clust_options = clust_options)))
    } else return(rep(clust, 1))
  }
  if (nrow(data[idxs,]) <= minpts) {
    if (ret_model){
      return(list(cluster = rep(clust, nrow(data[idxs,])), model = list(cluster = clust, data_columns = ncol(data), maxevts = maxevts, maxlvl = maxlvl, minpts = minpts, clust_options = clust_options)))
    } else return(rep(clust, nrow(data[idxs,])))
  }
  if (nrow(data[idxs,])>maxevts) {
    datasample<-sample(1:nrow(data[idxs,]), maxevts, replace = FALSE)
    knn_needed<-T
  } else datasample<-1:nrow(data[idxs,])
  method = clust_options$method
  mineps = clust_options$mineps
  mindens = clust_options$mindens
  bw = clust_options$bw
  nbins = clust_options$nbins
  mvpratio = clust_options$mvpratio
  # UMAP modeling
  umap_model<-cuml$UMAP(n_epochs = 500L)$fit(data[idxs,][datasample,])
  embdata<-umap_model$embedding_
  # clustering UMAP embedding
  clusters<-switch (method,
                    "dbscan" = {
                      num_el<-nrow(embdata)
                      vdist<-sort(as.vector(FNN::knn.dist(embdata, minpts)))
                      x1<-length(vdist)-num_el
                      eps<-max(mineps, vdist[x1])
                      res<-dbscan::dbscan(embdata, eps = eps, minPts = minpts)
                      if (any(res$cluster == 0)) clusters<-res$cluster + 1L else clusters<-res$cluster
                      clusters
                    },
                    "kdbscan" = {
                      kdbscan(embdata, minpts = max(10,minpts*(nrow(embdata)/nrow(data[idxs,]))), mindens = mindens, ret_model = T)
                    },
                    "hdbscan" = {
                      res<-dbscan::hdbscan(embdata, minPts = minpts)
                      if (any(res$cluster == 0)) clusters<-res$cluster + 1L else clusters<-res$cluster
                      clusters
                    },
                    "sdbscan" = {
                      sdbscan(embdata, minpts = max(10,minpts*(nrow(embdata)/nrow(data[idxs,]))), bw = bw, nbins = nbins, mvpratio = mvpratio, ret_model = T)
                    }
  )
  if (inherits(clusters, "s2dcluster")) {
    cl_model<-clusters
    clusters<-clusters$cluster
  } else cl_model<-NULL
  # cluster not sampled data
  if (length(unique(clusters)) > 1) {
    if (knn_needed) {
      new_embdata<-umap_model$transform(data[idxs,][-datasample,])
      if (is.null(cl_model)){
        knn<-FNN::get.knnx(embdata, new_embdata, k=1)
        newclusters<-clusters[knn$nn.index]
      } else {
        newclusters<-predict.s2dcluster(cl_model, new_embdata)
      }
      all_clusters<-integer(nrow(data[idxs,]))
      all_clusters[datasample]<-clusters
      all_clusters[-datasample]<-newclusters
      clusters<-all_clusters
    }
    new_clusters<-sort(unique(clusters))
    new_clusters<-paste0(clust, "_", new_clusters)
    clusters<-new_clusters[clusters]
  } else {
    clusters<-(rep(clust, nrow(data[idxs,])))
    new_clusters<-clust
  }
  newmodels<-list()
  if ((maxlvl > lvl) && length(new_clusters) > 1) {
    if (!is.null(myqueue)) myqueue$producer$fireAssignReactive("umap_dbs_msg", paste0("Found ", length(new_clusters), " new clusters on ", clust, " (level ", lvl, "). Trying deeper levels."))
    for (c in new_clusters) {
      nidxs<-idxs
      nidxs[idxs]<-nidxs[idxs] & (clusters == c)
      c<-c
      assign(paste0("nclust", c), .cuml_umap_dbs_clust(data, maxevts, maxlvl, minpts, clust_options, idxs = nidxs, lvl = lvl+1, clust = c, ret_model = ret_model, myqueue = myqueue, verbose = verbose, cuml = cuml))
    }
    for (c in new_clusters) {
      if (ret_model) {
        tempclust<-get(paste0("nclust", c))
        clusters[clusters == c]<-tempclust$cluster
        if (length(unique(tempclust$cluster)) > 1) newmodels[[c]]<-tempclust$model
      } else clusters[clusters == c]<-get(paste0("nclust", c))
    }
  } else if (!is.null(myqueue)) myqueue$producer$fireAssignReactive("umap_dbs_msg", paste0("Found only 1 cluster on ", clust, " (level ", lvl, ") or no deeper level allowed. Skipping deeper levels on this cluster."))
  if (lvl == 1) if (!is.null(myqueue)) myqueue$producer$fireAssignReactive("umap_dbs_msg", paste0("Found a total of ", length(unique(clusters)), " clusters."))
  if (ret_model) {
    umap_file<-tempfile()
    uwot::save_uwot(umap_model, umap_file)
    thismodel<-list(cluster = clust, data_columns = ncol(data), maxevts = maxevts, maxlvl = maxlvl, minpts = minpts, clust_options = clust_options)
    if (length(unique(clusters)) > 1) thismodel<-c(thismodel, list(umap_file = umap_file, umap_sample = datasample, s2dcluster_model = cl_model))
    if (length(newmodels) > 0) thismodel<-c(thismodel, newmodels)
    return(list(cluster = clusters, model = thismodel))
  } else return(clusters)
}

cuml_umap_dbs_clust<-function(data, maxevts, maxlvl, minpts, clust_options, ret_model = F, myqueue = NULL, verbose = F) {
  ret_model<-F # hardcode ret_model until implemented!
  if (ret_model) {
    if (!(clust_options$method %in% c("sdbscan", "kdbscan"))) stop("ret_model needs sdbscan or kdbscan clustering method!", call. = F)
  }
  cuml<-reticulate::import("cuml")
  clusters<-.cuml_umap_dbs_clust(data = data, maxevts = maxevts, maxlvl = maxlvl, minpts = minpts, clust_options = clust_options, ret_model = ret_model, myqueue = myqueue, verbose = verbose, cuml = cuml)
  if (ret_model) {
    clusters_model<-clusters$model
    clusters<-clusters$cluster
  }
  ctable<-as.data.frame(table(clusters))
  ctable<-ctable[order(ctable$Freq, decreasing = T),]
  new_clusters<-c(1:nrow(ctable))
  names(new_clusters)<-ctable$clusters
  clusters<-unname(new_clusters[clusters])
  if (ret_model) {
    res<-list(cluster = clusters, clusters=new_clusters, model = clusters_model)
    class(res)<-"UMAP_s2dcluster"
    return(res)
  } else return(clusters)
}

.plot.UMAP_s2dcluster_model<-function(model, data, ...) {
  umap_model<-uwot::load_uwot(model$umap_file)
  if (is.null(data)) {
    embdata<-umap_model$embedding
  } else embdata<-uwot::umap_transform(data, umap_model)
  par(fig=c(0,1,0.1,1), mar=c(2.1,2.1,2.1,1.1))
  plot.s2dcluster(model$s2dcluster_model, embdata, main = model$cluster, xlab="", ylab="", ...)
  clusts<-predict.s2dcluster(model$s2dcluster_model, embdata)
  clusters<-sort(unique(clusts))
  opar <- par(fig=c(0, 1, 0, 0.1), oma=c(0, 0, 0, 0), 
              mar=c(0, 0, 0, 0), new=TRUE)
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend("bottom", legend = paste0(model$cluster, "_", clusters), pt.cex = 1.5,
         col = palette()[clusters], pch = 16, horiz=TRUE, inset = 0.03)
  par(opar)
  for (i in clusters) {
    clust<-paste0(model$cluster, "_", i)
    if (!is.null(model[[clust]])) {
      .plot.UMAP_s2dcluster_model(model[[clust]], data[clusts == i,], ...)
    }
  }
}


plot.UMAP_s2dcluster<-function(res, data = NULL, ...) {
  if (!inherits(res, "UMAP_s2dcluster")) stop("res must be a UMAP_s2dcluster result object!", call. = F)
  if (!is.null(data)) if (!is.matrix(data) || ncol(data) != res$model$data_columns || nrow(data) < 1) stop(paste0("data must be a matrix with ", res$model$data_columns, " columns and >0 rows!"), call. = F)
  .plot.UMAP_s2dcluster_model(res$model, data, ...)
}

#### simple 2d clustering methods for data with sparse dense clusters ####
.kdbscan<-function(x, minpts = 100, mindens = 0.001, level = 1, idxs = 1:nrow(x)>0, cluster = "root", maxlvl = 100, alfa = 0.05, theta = 5, ret_model = F) {
  cut_x<-NULL
  for (.theta in seq(0, 90, theta)) {
    angle<-pi/180*.theta
    if (angle != 0) {
      x2<-x %*% matrix(c(cos(angle), sin(angle), -sin(angle), cos(angle)), ncol=2, byrow = T)
    } else x2<-x
    dens_0<-NULL
    for (axis in 1:2) {
      dens<-density(x2[idxs,axis], bw=.1)
      dens_0<-which(dens$y<mindens)
      xmin_alfa_corrected<-min(x2[idxs,axis]) + diff(range(x2[idxs,axis]))*alfa
      xmax_alfa_corrected<-max(x2[idxs,axis]) - diff(range(x2[idxs,axis]))*alfa
      dens_0<-dens_0[dens$x[dens_0]>xmin_alfa_corrected & dens$x[dens_0]<xmax_alfa_corrected]
      if (length(dens_0) > 1) {
        reps<-diff(c(0L,c(which(dens_0[-1L] != dens_0[-length(dens_0) ] + 1 ),length(dens_0))))
        max_rep<-match(max(reps), reps)
        cut_pos<-dens_0[sum(reps[1:max_rep])-abs(max(reps)/2)]
        if (any(as.integer(table(x2[idxs,axis] > dens$x[cut_pos])) < minpts)) next
        cut_x<-dens$x[cut_pos]
        cut_axis<-axis
        cut_angle<-angle
        break
      } else if (length(dens_0) == 1) {
        if (any(as.integer(table(x2[idxs,axis] > dens$x[dens_0])) < minpts)) next
        cut_x<-dens$x[dens_0]
        cut_axis<-axis
        cut_angle<-angle
        break
      }
    }
    if (!is.null(cut_x)) break
  }
  if (is.null(cut_x)) {
    if (ret_model){
      return(list(cluster = rep(cluster, nrow(x[idxs,])), model = list(cluster="none")))
    } else return(rep(cluster, nrow(x[idxs,])))
  }
  clust<-rep(paste0(cluster, "_1"), nrow(x[idxs,]))
  clust[x2[idxs,cut_axis] > cut_x]<-paste0(cluster, "_2")
  clmodel<-list(cluster=cluster, angle=cut_angle, axis=cut_axis, cut=cut_x)
  if (any(data.frame(table(clust))$Freq < minpts)) {
    if (ret_model){
      return(list(cluster = rep(cluster, nrow(x[idxs,])), model = list(cluster="none")))
    } else return(rep(cluster, nrow(x[idxs,])))
  }
  if (level < maxlvl) {
    new_clusters<-unique(clust)
    nmodel<-list()
    cmodel<-c()
    for (c in new_clusters) {
      c<-c
      if (length(clust[clust == c]) < minpts*2) next
      nidxs<-idxs
      nidxs[idxs]<-nidxs[idxs] & (clust == c)
      nclust<-.kdbscan(x, minpts = minpts, mindens = mindens, level = level+1, idxs = nidxs, cluster = c, maxlvl = maxlvl, theta = theta, ret_model = ret_model)
      if (ret_model) {
        if (nclust$model$cluster != "none") {
          cmodel<-c(cmodel, c)
          nmodel[[c]]<-nclust$model
        }
        nclust<-nclust$cluster
      }
      clust[clust == c]<-nclust
    }
    if (ret_model) {
      clmodel<-c(clmodel, nmodel[cmodel])
    }
  }
  if (ret_model) {
    return(list(cluster = clust, model = clmodel))
  } else return(clust)
  
}

kdbscan<-function(x, minpts = 100, mindens = 0.001, maxlvl = 100, alfa = 0.05, theta = 5, ret_model = F) {
  clusters<-.kdbscan(x, minpts = minpts, mindens = mindens, maxlvl = maxlvl, theta = theta, ret_model = ret_model)
  if (ret_model) {
    clusters_model<-clusters$model
    clusters<-clusters$cluster
  }
  ctable<-as.data.frame(table(clusters))
  ctable<-ctable[order(ctable$Freq, decreasing = T),]
  new_clusters<-c(1:nrow(ctable))
  names(new_clusters)<-ctable$clusters
  clusters<-unname(new_clusters[clusters])
  if (ret_model) {
    res<-list(cluster = clusters, clusters=new_clusters, x = substitute(x), model = clusters_model)
    class(res)<-"s2dcluster"
    return(res)
  } else return(clusters)
}

getdensdata<-function(x2, axis, minpts) {
  dens<-density(x2[,axis])
  x2_sorted<-sort(x2[,axis])
  if (!(x2_sorted[minpts] < x2_sorted[length(x2_sorted)-minpts+1])) return(NULL)
  v = optimize(approxfun(dens$x,dens$y),interval=c(x2_sorted[minpts], x2_sorted[length(x2_sorted)-minpts+1]))
  p1 = optimize(approxfun(dens$x,dens$y),interval=c(x2_sorted[minpts], v$minimum), maximum = T)
  p2 = optimize(approxfun(dens$x,dens$y),interval=c(v$minimum, x2_sorted[length(x2_sorted)-minpts+1]), maximum = T)
  return(c(v$minimum, v$objective, p1$objective, p2$objective))
}

.sdbscan<-function(x, minpts = 100, level = 1, idxs = 1:nrow(x)>0, cluster = "root", maxlvl = 100, alfa = 0.05, bw = 0.25, nbins = 1, theta = 5, mvpratio = 0.01, ret_model = F, plotcuts = F) {
  cut_x<-NULL
  vdata<-c()
  for (.theta in seq(0, 90, theta)) {
    angle<-pi/180*.theta
    if (angle != 0) {
      x2<-x %*% matrix(c(cos(angle), sin(angle), -sin(angle), cos(angle)), ncol=2, byrow = T)
    } else x2<-x
    dens_0<-NULL
    range<-list(x=range(x2[idxs,1]), y=range(x2[idxs,2]))
    bins<-c(as.integer(diff(range$x)/bw), as.integer(diff(range$y)/bw))
    for (axis in 1:2) {
      thisrange<-range[[axis]]
      newx<-as.integer(((x2[idxs,axis]-thisrange[1])/(thisrange[2]-thisrange[1]))*bins[[axis]])
      dens_0<-(min(newx):max(newx))[!((min(newx):max(newx)) %in% newx)]
      dens_0<-dens_0[dens_0>(min(newx)+bins[[axis]]*alfa) & dens_0<(max(newx)-bins[[axis]]*alfa)]
      if (length(dens_0) >= nbins) {
        reps<-diff(c(0L,c(which(dens_0[-1L] != dens_0[-length(dens_0) ] + 1 ),length(dens_0))))
        if (!any(reps >= nbins)) next
        rep_order<-order(reps, decreasing = T)
        cut_pos<-NULL
        for (r in seq_along(rep_order)) {
          max_rep<-rep_order[r]
          if (reps[max_rep] < nbins) break
          cut_pos<-dens_0[sum(reps[1:max_rep])-floor(abs(reps[max_rep]/2))]
          if (!any(as.integer(table(newx > cut_pos)) < minpts)) break
          cut_pos<-NULL
        }
        if (is.null(cut_pos)) next
        cut_x<-(((cut_pos+0.5)/bins[[axis]])*(thisrange[2]-thisrange[1]))+thisrange[1]
        cut_axis<-axis
        cut_angle<-angle
        break
      }
      if (mvpratio > 0) {
        res<-getdensdata(x2[idxs,], axis, minpts)
        if (!is.null(res)) vdata<-c(vdata, angle, axis, res)
      }
    }
    if (!is.null(cut_x)) break
  }
  if (is.null(cut_x)) {
    if (mvpratio > 0 && !is.null(vdata)) {
      vdata<-data.frame(matrix(vdata, byrow = T, ncol=6, dimnames = list(NULL, c("angle", "axis", "newx", "value", "maxval_before", "maxval_after"))))
      v_cuts<-vdata[vdata$value<=vdata$maxval_before*mvpratio & vdata$value<=vdata$maxval_after*mvpratio,]
      v_cuts$pvr1<-v_cuts$maxval_before/v_cuts$value
      v_cuts$pvr2<-v_cuts$maxval_after/v_cuts$value
      v_cuts$pvr_sum<-v_cuts$pvr1+v_cuts$pvr2
      v_cuts<-v_cuts[order(v_cuts$value, decreasing = F),]
      for (v in seq_along(v_cuts[,1])){
        v_cut<-v_cuts[v,]
        x2<-x %*% matrix(c(cos(v_cut$angle), sin(v_cut$angle), -sin(v_cut$angle), cos(v_cut$angle)), ncol=2, byrow = T)
        cut_x<-v_cut$newx
        if (any(table(x2[idxs,v_cut$axis]>cut_x) < minpts)) {
          cut_x<-NULL
          next
        }
        cut_axis<-v_cut$axis
        cut_angle<-v_cut$angle
        break
      } 
      if (is.null(cut_x)) {
        if (ret_model){
          return(list(cluster = rep(cluster, nrow(x[idxs,])), model = list(cluster="none")))
        } else return(rep(cluster, nrow(x[idxs,])))
      }
    } else {
      if (ret_model){
        return(list(cluster = rep(cluster, nrow(x[idxs,])), model = list(cluster="none")))
      } else return(rep(cluster, nrow(x[idxs,])))
    }
  }
  clust<-rep(paste0(cluster, "_1"), nrow(x[idxs,]))
  clust[x2[idxs,cut_axis] > cut_x]<-paste0(cluster, "_2")
  if (plotcuts) sdbscan_hist(x2[idxs,], bw = bw, cut_axis=cut_axis, cut_x=cut_x, cex=2, pch = ".", 
                             main = paste0(cluster, " angle = ", cut_angle*180/pi),
                             cex.main = 0.8,
                             col=factor(clust))
  clmodel<-list(cluster=cluster, angle=cut_angle, axis=cut_axis, cut=cut_x)
  if (any(data.frame(table(clust))$Freq < minpts)) {
    if (ret_model){
      return(list(cluster = rep(cluster, nrow(x[idxs,])), model = list(cluster="none")))
    } else return(rep(cluster, nrow(x[idxs,])))
  }
  if (level < maxlvl) {
    new_clusters<-unique(clust)
    nmodel<-list()
    cmodel<-c()
    for (c in new_clusters) {
      c<-c
      if (length(clust[clust == c]) < minpts*2) next
      nidxs<-idxs
      nidxs[idxs]<-nidxs[idxs] & (clust == c)
      nclust<-.sdbscan(x, minpts = minpts, level = level+1, idxs = nidxs, cluster = c, maxlvl = maxlvl, alfa = alfa, bw = bw, nbins = nbins, theta = theta, mvpratio = mvpratio, ret_model = ret_model, plotcuts = plotcuts)
      if (ret_model) {
        if (nclust$model$cluster != "none") {
          cmodel<-c(cmodel, c)
          nmodel[[c]]<-nclust$model
        }
        nclust<-nclust$cluster
      }
      clust[clust == c]<-nclust
    }
    if (ret_model) {
      clmodel<-c(clmodel, nmodel[cmodel])
    }
  }
  if (ret_model) {
    return(list(cluster = clust, model = clmodel))
  } else return(clust)
}

sdbscan<-function(x, minpts = 100, maxlvl = 100, alfa = 0.05, bw = 0.25, nbins = 1, theta = 5, mvpratio = 0.01, ret_model = F, plotcuts = F) {  ### TODO: check parameters
  clusters<-.sdbscan(x, minpts = minpts, maxlvl = maxlvl, alfa = alfa, bw = bw, nbins = nbins, theta = theta, mvpratio = mvpratio, ret_model = ret_model, plotcuts = plotcuts)
  if (ret_model) {
    clusters_model<-clusters$model
    clusters<-clusters$cluster
  }
  ctable<-as.data.frame(table(clusters))
  ctable<-ctable[order(ctable$Freq, decreasing = T),]
  new_clusters<-c(1:nrow(ctable))
  names(new_clusters)<-ctable$clusters
  clusters<-unname(new_clusters[clusters])
  if (ret_model) {
    res<-list(cluster = clusters, clusters=new_clusters, x = substitute(x), model = clusters_model)
    class(res)<-"s2dcluster"
    return(res)
  } else return(clusters)
}

.plot_cut<-function(x, angle, axis, cut, cuts) {
  xmax = 100
  xmin = -100
  dt<-matrix(c(rep(cut, 2), c(xmin, xmax)), ncol=2)
  if (axis == 2) dt<-dt[,c(2,1)]
  dtr<-dt %*% matrix(c(cos(-angle), sin(-angle), -sin(-angle), cos(-angle)), ncol=2, byrow = T)
  inter<-get_intersections(dtr, cuts)
  inter<-inter[!is.infinite(inter[,1]),]
  poly<-x[chull(x),]
  segs<-poly2seg(poly)
  interx<-get_intersections(dtr, segs, onsegment = TRUE)
  if ((inter[1,1] - inter [2,1] < 1e-5) && (inter[1,1] - inter [2,1] > -1e-5)) {#dtr is vertical!
    inter<-inter[order(inter[,2]),]
    interminus<-matrix(inter[inter[,2]<=min(interx[,2]),], ncol = 2)
    interplus<-matrix(inter[inter[,2]>=max(interx[,2]),], ncol = 2)
    dtr<-rbind(interminus[nrow(interminus),], interplus[1,])
    
  } else {
    inter<-inter[order(inter[,1]),]
    interminus<-matrix(inter[inter[,1]<=min(interx[,1]),], ncol = 2)
    interplus<-matrix(inter[inter[,1]>=max(interx[,1]),], ncol = 2)
    dtr<-rbind(interminus[nrow(interminus),], interplus[1,])
  }
  lines(dtr)
  return(dtr)
}

.plot_model<-function(x, res, model, cuts) {
  x2<-x %*% matrix(c(cos(model$angle), sin(model$angle), -sin(model$angle), cos(model$angle)), ncol=2, byrow = T)
  cut<-.plot_cut(x, model$angle, model$axis, model$cut, cuts)
  new_clusters<-c(paste0(model$cluster, "_1"), paste0(model$cluster, "_2"))
  for (c in 1:2) {
    c<-c
    if (is.null(model[[new_clusters[c]]])) next
    if (c == 1) {
      idxs<-x2[,model$axis]<=model$cut
    } else {
      idxs<-x2[,model$axis]>model$cut
    }
    .plot_model(x[idxs,], res, model[[new_clusters[c]]], rbind(cuts, cut))
  }
}

plot.s2dcluster<-function(res, x = NULL, cex=2, ...) {
  if (!inherits(res, "s2dcluster")) stop("res must be a s2dcluster result object!", call. = F)
  if (is.null(x)) {
    x<-try(eval(res$x), silent = T)
    if (!is.matrix(x) || ncol(x) != 2 || nrow(x) < 1) stop("x must be a matrix with 2 columns and >0 rows!", call. = F)
    col<-res$cluster
  } else {
    if (!is.matrix(x) || ncol(x) != 2 || nrow(x) < 1) stop("x must be a matrix with 2 columns and >0 rows!", call. = F)
    col<-predict.s2dcluster(res, as.matrix(x, ncol=2))
  }
  plot(x, pch=".", cex=cex, col=col, ...)
  if (res$model$cluster != "none") {
    .plot_model(x, res, res$model, make_box(10*(max(abs(c(range(x[,1]), range(x[,2])))))))
  }
}

print.s2dcluster<-function(res) {
  if (!inherits(res, "s2dcluster")) stop("res must be a s2dcluster result object!", call. = F)
  cat("s2dcluster result object.\n")
  cat("Data file:", res$x, "with", length(res$cluster), "points.\n")
  cat("Found", length(res$clusters), "clusters.\n")
  cat("Cluster table:")
  print(table(res$cluster))
  cat("Cluster hierarchy:\n")
  cltable<-data.frame(Cluster = res$clusters[order(names(res$clusters))])
  print(cltable)
}

.s2dcluster_newdata<-function(x, model) {
  if (class(x)[1] == "numeric") x<-matrix(x, ncol = 2)
  x2<-x %*% matrix(c(cos(model$angle), sin(model$angle), -sin(model$angle), cos(model$angle)), ncol=2, byrow = T)
  clust<-rep(paste0(model$cluster, "_1"), nrow(x2))
  clust[x2[,model$axis] > model$cut]<-paste0(model$cluster, "_2")
  new_clusters<-unique(clust)
  for (c in new_clusters) {
    c<-c
    if (is.null(model[[c]])) next
    idxs<-clust == c
    nclust<-.s2dcluster_newdata(x[idxs,], model[[c]])
    clust[idxs]<-nclust
  }
  return(clust)
}

predict.s2dcluster<-function(res, newdata) {
  if (!inherits(res, "s2dcluster")) stop("res must be a s2dcluster result object!", call. = F)
  if (!is.matrix(newdata) || ncol(newdata) != 2 || nrow(newdata) < 1) stop("newdata must be a matrix with 2 columns and >0 rows!", call. = F)
  cluster<-c()
  if (res$model$cluster != "none") {
    cluster<-.s2dcluster_newdata(newdata, res$model)
  } else cluster<-rep("root", nrow(newdata))
  clusters<-unname(res$clusters[cluster])
  return(clusters)
}

substract<-function(x) {
  return(x[1]-x[2])
}

#source: Weisstein, Eric W. "Line-Line Intersection." From MathWorld--A Wolfram Web Resource. https://mathworld.wolfram.com/Line-LineIntersection.html
get_intersection<-function(line, cut, onsegment = FALSE) {
  diffline<-apply(line, 2, substract)
  diffline[abs(diffline)<1e-5]<-0
  diffcut<-apply(cut, 2, substract)
  diffcut[abs(diffcut)<1e-5]<-0
  detall<-det(rbind(diffline, diffcut))
  if (detall == 0) return(c(Inf, Inf))
  
  det1<-det(line)
  det2<-det(cut)
  x<-det(matrix(c(det1, diffline[1], 
                  det2, diffcut[1]), ncol = 2, byrow = T))/detall
  y<-det(matrix(c(det1, diffline[2], 
                  det2, diffcut[2]), ncol = 2, byrow = T))/detall
  if (onsegment) {
    if (!(
      x <= max(line[,1]) + 1e-5 &&
      x <= max(cut[,1]) + 1e-5 &&
      x >= min(line[,1]) - 1e-5 &&
      x >= min(cut[,1]) - 1e-5 &&
      y <= max(line[,2]) + 1e-5 &&
      y <= max(cut[,2]) + 1e-5 &&
      y >= min(line[,2]) - 1e-5 &&
      y >= min(cut[,2]) - 1e-5
    )) return(NULL)
  }
  return(c(x, y))
}

get_intersections<-function(line, cuts, ...) {
  inter<-matrix(NA_real_, ncol = 2, nrow = 0)
  for (i in 1:(nrow(cuts)/2)) {
    inter<-rbind(inter, get_intersection(line, cuts[c(2*i-1, 2*i),], ...))
  }
  return(inter)
}

poly2seg<-function(poly) {
  rows<-nrow(poly)
  segs<-matrix(NA_real_, ncol = 2, nrow = rows*2)
  segs[(1:rows)*2-1,]<-poly
  segs[(1:rows)*2,]<-poly[c(2:rows, 1),]
  return(segs)
}

make_box<-function(x) {
  p1<-c(-x, -x)
  p2<-c(x, -x)
  p3<-c(x, x)
  p4<-c(-x, x)
  return(rbind(p1, p2, p2, p3, p3, p4, p4, p1))
}

sdbscan_hist<-function(x, xlab="", ylab="", bw, cut_axis, cut_x, ...){
  layout(matrix(c(2,0,1,3), ncol=2, byrow=TRUE), widths=c(4/5,1/5), heights=c(1/5,4/5))
  xhist = hist(x[,1], plot=FALSE, breaks = as.integer(diff(range(x[,1]))/bw))
  yhist = hist(x[,2], plot=FALSE, breaks = as.integer(diff(range(x[,2]))/bw))
  oldpar<-par(no.readonly = T)
  par(mar=c(3,3,1,1))
  plot(x, xlab = "", ylab = "", ...)
  if (cut_axis == 1) {
    abline(v = cut_x)
  } else abline(h = cut_x)
  par(mar=c(0,3,1,1))
  barplot(xhist$counts, axes=FALSE, ylim=c(0, max(xhist$counts)), space=0)
  par(mar=c(3,0,1,1))
  barplot(yhist$counts, axes=FALSE, xlim=c(0, max(yhist$counts)), space=0, horiz=TRUE)
  par(oma=c(3,3,0,0))
  mtext(xlab, side=1, line=1, outer=TRUE, adj=0, 
        at=.8 * (mean(x[,1]) - min(x[,1]))/(max(x[,1])-min(x[,1])))
  mtext(ylab, side=2, line=1, outer=TRUE, adj=0, 
        at=(.8 * (mean(x[,2]) - min(x[,2]))/(max(x[,2]) - min(x[,2]))))
  par(oldpar)
}