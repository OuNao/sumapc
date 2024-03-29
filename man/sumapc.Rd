% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sumapc.R
\name{sumapc}
\alias{sumapc}
\title{Sequential UMAP 2d embedding and simple 2d space partitioning clustering.}
\usage{
sumapc(
  data,
  maxevts = 10000L,
  transform_maxevts = 1000000L,
  maxlvl = 3L,
  minpts = 100L,
  clust_options = list(method = "sdbscan", mineps = 1, mindens = 0.1, bw = 0.05, nbins
    = 5, mvpratio = 0.5),
  multi_thread = TRUE,
  fast_sgd = TRUE,
  ret_model = FALSE,
  myqueue = NULL,
  verbose = FALSE,
  use_cuml = FALSE,
  seed = NULL,
  ...
)
}
\arguments{
\item{data}{a data matrix.}

\item{maxevts}{max number of random points to apply UMAP.}

\item{transform_maxevts}{max number of points to apply UMAP transform. The points not used in the initial UMAP model are divided in batches of transform_maxevts points to apply UMAP transform. Lower value = lower memory usage.}

\item{maxlvl}{max sequential levels to process.}

\item{minpts}{min cluster size.}

\item{clust_options}{list; 2d clustering options (see Details).}

\item{multi_thread}{logical. Needs \link[future]{future} package loaded and multisession plan selected to use parallel clustering (not only multithreading in UMAP). Only used if use_cuml = FALSE.}

\item{fast_sgd}{logical; see \link[uwot]{umap}.}

\item{ret_model}{logical; should clustering model be exported.}

\item{myqueue}{shinyQueue object. Multithreaded code use \link[future]{future} package. This object is used to send log messages to the base worker. See \link[ipc]{shinyQueue}.}

\item{verbose}{logical.}

\item{use_cuml}{logical. Enable CUML GPU compute acceleration.}

\item{seed}{Integer or NULL.}

\item{...}{additional parameters passed to UMAP}
}
\value{
A vector of cluster numbers with length = nrow(data)
}
\description{
Sequential UMAP 2d embedding and simple 2d space partitioning clustering.
}
\details{
clust_options must be a named list containing:
 method: one of "dbscan", "hdbscan", "kdbscan", "sdbscan"
 mineps: min eps value passed to dbscan
 mindens: min density passed to kdbscan
 bw: bw passed to sdbscan
 nbins: min number of bins of bw width
 mvpratio: max valley/peak ratio)
 Any other value is passed as is to the clustering method.
}
