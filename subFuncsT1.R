# this is the network contruction using the KEGG signaling pathways
# Usage: 1) loadKeggPathways --> 2) generateKeggNetInfoV1 (for each cell type, with log2(fold change) and gene symbols)
# 3) generate 'all' activated signaling network of the given cell type using "generateKeggActivationNetworkV1"
# 4) generate only the 'downstream signaling' of given receptors of the given cell type using "generateKeggReceptorNetworkV1"

# ------------------
# generateKeggNetworkV1
# ------------------
generateKeggReceptorNetworkV1 <- function(netInfo, Rs, T0){
	# library(igraph)
	# Rs: receptors, like: Rs <- c("EGFR")
	# T0: threshold,e.g., log2(1.25)
	# T0 <- log2(1.25) # threshold with mean-node score >= T0

	if (1 == 1){ # network reconstruction 
		
		xt <- netInfo[[1]]
		v_t_a <- netInfo[[2]]
		v_t_ai <- netInfo[[3]]

		# T0 <- log2(1.25) # threshold with mean-node score >= T0
		T1 <- log2(1.25) # threshold with mean-node score <= -T1
		path_t <- c('Source', 'Target', 'PathwayName')
		path_ti <- c('Source', 'Target', 'PathwayName')
		# ---
		nt <- length(xt)
		vt <- rep(0.0, nt)
		for (j in 1:nt){
			yt <- xt[[j]]
			zt1a <- yt[[1]][1]
			zt1b <- yt[[1]][2]
			zt2 <- yt[[2]]
			p_n <- zt2[1,3]  # pathway name
			#
			lt <- dim(zt2)[1]
			if (lt <=2){
				next
			}
			#
			if (zt1a >= T0 & zt1b > v_t_a[lt]){  # activation
					if (length(which(zt2[,1] %in% Rs | zt2[,1] %in% Rs))>0){ # include the receptors
						path_t <- rbind(path_t, zt2)
					}
			}
			if (zt1a >= T1 & zt1b < v_t_ai[lt]){  # inhibition
					path_ti <- rbind(path_ti, zt2)
			}

		} # end of network construction
	}
	path_t <- unique(path_t[,c(1:2)])
	return(path_t)

}
# ------------------
# generateKeggNetworkV1
# ------------------
generateKeggActivationNetworkV1 <- function(netInfo, T0){
	# library(igraph)
	# Rs: receptors, like: Rs <- c("EGFR")
	# T0: threshold,e.g., log2(1.25)
	# T0 <- log2(1.25) # threshold with mean-node score >= T0

	if (1 == 1){ # network reconstruction 
		
		xt <- netInfo[[1]]
		v_t_a <- netInfo[[2]]
		v_t_ai <- netInfo[[3]]

		# T0 <- log2(1.25) # threshold with mean-node score >= T0
		T1 <- log2(1.25) # threshold with mean-node score <= -T1
		path_t <- c('Source', 'Target', 'PathwayName')
		path_ti <- c('Source', 'Target', 'PathwayName')
		# ---
		nt <- length(xt)
		vt <- rep(0.0, nt)
		for (j in 1:nt){
			yt <- xt[[j]]
			zt1a <- yt[[1]][1]
			zt1b <- yt[[1]][2]
			zt2 <- yt[[2]]
			p_n <- zt2[1,3]  # pathway name
			#
			lt <- dim(zt2)[1]
			if (lt <=2){
				next
			}
			#
			if (zt1a >= T0 & zt1b > v_t_a[lt]){  # activation
					path_t <- rbind(path_t, zt2)
			}
			if (zt1a >= T1 & zt1b < v_t_ai[lt]){  # inhibition
					path_ti <- rbind(path_ti, zt2)
			}
			
		}
		#
		# print('# of Edge in path_t(activation):'); print(dim(unique(path_t)))
		# print('# of Edge in path_ti(inhibition):'); print(dim(unique(path_ti)))
	} # end of network construction
	path_t <- unique(path_t[,c(1:2)])
	return(path_t)

}

#############
### function: loadKeggPathways
#############

loadKeggPathways <- function(){
	# get the KEGG signaling pathways
	if (1 == 1){ # load in KeggPathways
      library(graphite)
      library(igraph)

      kps <- pathways("hsapiens", "kegg")
      names_pathways <- names(kps)
      n_p <- length(kps)
      KeggPathways <- list()
      for (i in 1:n_p){
        pt <- kps[[i]]
        pt <- convertIdentifiers(pt, 'symbol')
        et1 <- pt@mixedEdges
        et2 <- pt@protEdges
        et <- rbind(et1, et2)
        KeggPathways[[i]] <- et
      }
      names(KeggPathways) <- names_pathways
      #
      KeggGenes <- c()
      KeggGraphs <- c()
      k <- 0
      for (i in 1:length(KeggPathways)){
        pt <- KeggPathways[[i]]
        gs <- c()
        xt <- pt[,1]; idxt1 <- regexpr('SYMBOL', xt); idxt1 <- which(idxt1>0); nt1 <- sum(idxt1>0)
        if (nt1>0){gs1 <- pt[,2]; gs1 <- gs1[which(idxt1>0)]; gs <- union(gs, gs1)}
        xt <- pt[,3]; idxt2 <- regexpr('SYMBOL', xt); idxt2 <- which(idxt2>0); nt2 <- sum(idxt2>0)
        if (nt2>0){gs2 <- pt[,4]; gs2 <- gs2[which(idxt2>0)]; gs <- union(gs, gs1)}
        if (max(nt1, nt2)>0){ KeggGenes <- union(KeggGenes, gs)}
      	
      	idxt <- intersect(idxt1, idxt2)
      	nt <- length(idxt)
      	if (nt > 1){
        	e_t <- pt[idxt,c(2,4)]
        	KeggGraphs[[i]] <- graph.edgelist(as.matrix(e_t), directed=T)
        	k <- k + 1
         }
      }
      print(k)

  	}
  	keggInfo <- list()
  	keggInfo[[1]] <- KeggGenes
  	keggInfo[[2]] <- KeggPathways
  	keggInfo[[3]] <- KeggGraphs
  	return(keggInfo)
}
# ------------------
# generateKeggReceptorNetworkV1
# ------------------
generateKeggNetInfoV1 <- function(gSyms, fc_tmp, keggInfo){
	# library(igraph)
	# gSyms: gene symbol
	# fc_tmp: fold change (log2-scale) / or significance score (with mean = 0)
	#
	KeggGenes <- keggInfo[[1]]
	KeggPathways <- keggInfo[[2]]
	KeggGraphs <- keggInfo[[3]]

	Ta <- 3.0
	fc_tmp[fc_tmp >Ta] = Ta  #
	fc_tmp[fc_tmp < -1*Ta] = -1*Ta  #

	# get the genes in KEGG only
	idx_gSym <- which(gSyms %in% KeggGenes)
	gSym1 <- gSyms[idx_gSym]
	fc_tmp1 <- fc_tmp[idx_gSym]
	gSyms <- gSym1
	# take the maximum or minimum value of multiple genes
	st1 <- unique(gSyms)
	nt <- length(st1)
	fc_tmp <- rep(0.0, nt)
	for (i in 1:nt){
		idxt <- which(gSyms %in% st1[i])
		vt <- fc_tmp1[idxt]
		if (length(vt)>1){
			idxt1 <- which(abs(vt) == max(abs(vt)))
			fc_tmp[i] <- vt[idxt1[1]]
		} else {
			fc_tmp[i] <- vt
		}
	}
	gSyms <- st1
	#
	length(gSyms)
	length(fc_tmp)
	max(fc_tmp)
	min(fc_tmp)
	#

	if (1 == 1){ # network reconstruction 
		lPath <- c()
		sigPath <- c()
		sigPath <- getSubPathwayKeggV2a1(fc_tmp, gSyms, KeggPathways, KeggGraphs, lPath)
		pt1 <- sigPath$pathTmp
		lPath <- sigPath$lPath

		# extract the paths
		path_score_1 <- c()
		path_score_2 <- c()
		path_score_v <- list()
		path_score_a <- list()
		for (i in 1:50){
			path_score_v[[i]] <- c(0)
			path_score_a[[i]] <- c(0)
		}
		#
		path_t <- c('Source', 'Target', 'PathwayName')
		path_ti <- c('Source', 'Target', 'PathwayName')
		xt <- pt1
		#
		nt <- length(xt)
		vt <- rep(0.0, nt)
		for (j in 1:nt){
			yt <- xt[[j]]
			zt1a <- yt[[1]][1]
			zt1b <- yt[[1]][2]
			zt2 <- yt[[2]]
			p_n <- zt2[1,3]  # pathway name
			#
			lt <- dim(zt2)[1]
			path_score_v[[lt]] <- c(path_score_v[[lt]], zt1a)
			path_score_a[[lt]] <- c(path_score_a[[lt]], zt1b)
			#
			path_score_1 <- c(path_score_1, zt1a)
			path_score_2 <- c(path_score_2, zt1b)
		}
		#
		v_t_a <- rep(0.25, 50)  # activation threshold for the signaling paths with different lengths
		v_t_ai <- rep(-0.25, 50) # inhibition score threshold 
		#
		s_p_a <- matrix(0, 50,3) # to calculate the Real # of signaling pathways, mean-activation, mean-inhibition
		s_p_v <- matrix(0, 50,2) # vertex scores-mean 
		n_t <- matrix(0, 50, 2)
		# for (i in 1:50){
		C1 <- 1.25
		for (i in 1:50){
			# activation score
			vt1 <- path_score_a[[i]]
			idxt <- which(vt1 >0)
			if (length(idxt) < 2){next}
			vt1a <- vt1[vt1>0]
			# v_t_a[i] <- max(0.1, mean(vt1a)+1.0*sd(vt1a));
			v_t_a[i] <- mean(vt1a) + C1*sd(vt1a)
			n_t[i,1] <- (sum(vt1a >= v_t_a[i]))
			idxt <- which(vt1 <0)
			if (length(idxt) < 2){next}
			vt1 <- path_score_a[[i]]
			vt1a <- vt1[vt1<0]  # inhibition
			# v_t_ai[i] <- min(-0.1, mean(vt1a)-1.0*sd(vt1a)); 
			v_t_ai[i] <- mean(vt1a) - C1*sd(vt1a);
			n_t[i,2] <- (sum(vt1a <= v_t_ai[i]))
			# print(sum(vt1 < -v_t_ai[i]))
			s_p_a[i,] <- c(length(vt1), mean(vt1[vt1 > 0]), mean(vt1[vt1 < 0]))
			# vertex score
			vt1 <- path_score_v[[i]]
			s_p_v[i,] <- c(length(vt1), mean(vt1[vt1 > 0]))
		}
		v_t_a[c(11:50)] <- 5.0  # will not consider the paths with more than 11 nodes
		v_t_ai[c(11:50)] <- -5.0
		#
		netInfo <- list()
		netInfo[[1]] <- pt1;
		netInfo[[2]] <- v_t_a;
		netInfo[[3]] <- v_t_ai;

	} # end of network construction
	return(netInfo)

}
# ------------------
# getSubPathwayKeggV2a1
# ------------------
getSubPathwayKeggV2a1 <- function(fc1, gSyms, KeggPathways, KeggGraphs, lPath){
	#
	# print('test')
	pathway_names <- names(KeggPathways)
	n_pathway <- length(KeggPathways)
	# print('npathway')
	# print(n_pathway)
	pathTmp <- list() # save all the paths individually
	n_pt <- 0  # number of paths
	
	n_pathway <- length(KeggGraphs)
	for (j1 in 1:n_pathway){  # for all signaling pathways
		print(j1)
		p_name <- pathway_names[j1]
		gTmp <- KeggGraphs[[j1]]

		if (is.null(gTmp)){print('no pathway'); next;} # no network/pathway 
		# check if there are TFs and Receptors in this signaling pathway
		nodes_tmp <- V(gTmp)$name
		# z-scores of nodes_tmp
		idxt <- which(gSyms %in% nodes_tmp)
		zScore <- fc1[idxt]
		names(zScore) <- gSyms[idxt]  # for each sample
		# get the source/sink nodes
		# xt1 <- KeggPathways[[j1]]
		# n_From <- xt1[,2]
		# n_To <- xt1[,4]

		xt1 <- get.edgelist(gTmp)
		n_From <- xt1[,1]
		n_To <- xt1[,2]
		n_Source <- setdiff(n_From, n_To)
		n_Sink <- setdiff(n_To, n_From)
		#
		rec_set1 <- n_Source
		tf_set1 <- n_Sink
		n_rec1 <- length(rec_set1)
		n_tf1 <- length(tf_set1)
		#
		if (min(n_rec1, n_tf1) < 1){  # no possible Receptor - TF path
			next
		}
		# print(j1)
	  for (j2 in 1:n_rec1){  # each source point
			# print("j2"); print(j2)
			pf1 <- rec_set1[j2]  # from point
		  	paths <- get.shortest.paths(gTmp, pf1, tf_set1, mode='out')
			# paths <- all_simple_paths(gTmp, pf1, tf_set1, mode='out')
			paths <- paths$vpath
			nPath <- length(paths)
			#
			if (nPath > 0){  # connected
		  	for (k in 1:nPath){  # all paths
					pt <- paths[[k]]
					pt <- pt$name
					nPt <- length(pt)
					if (nPt < 2){  # no path (dis == Inf)
						next
					}
					# for each signaling path
					pathTmp1 <- c('Source', 'Target', 'pathwayName')
					for (l in 1:(nPt-1)){
						vt <- c(pt[l], pt[l+1], p_name)
						pathTmp1 <- rbind(pathTmp1, vt)
					}
					# activation evaluation
					if (nrow(pathTmp1) == 2){
						eTmp1 <- pathTmp1[-1,]
						dim(eTmp1) <- c(1,3)
					} else {
						eTmp1 <- pathTmp1[-1,]
					}
					#
					# st <- getScorePathwayKeggV1a(eTmp1, KeggPathways[[j1]], zScore, 0.0)
					st <- getScorePathwayKeggV1b(eTmp1, KeggPathways[[j1]], zScore, 0.0)
					xt <- list()
					xt[[1]] <- st
					xt[[2]]	<- eTmp1
					n_pt <- (n_pt + 1)
					pathTmp[[n_pt]] <- xt
				}
			}
		}
	}	# end of each signaling pathways?
	# pathways[[j]] <- pathTmp  # pathways of the j-th patient
	res1 <- list(pathTmp=pathTmp, lPath=lPath)
	return(res1)
}

# --------------------
# getScorePathwayKeggV1b
# --------------------
getScorePathwayKeggV1b <- function(edgePath, eType, zScore, alpha1){
	# edgePath: nx2 dimension; each row is an edge
	# eType: edge type, e.g., activation, inhibition,
	# zScore: gene expression of genes, e.g., fold change, z-score, protein expression
	# alpha1: weight of the contribution of the parent node
	vScore <- 0.0  # pathway score initialization
	aScore <- 0.0  # activation score
	# pre-processing
	gSym <- names(zScore)
	nEdge <- nrow(edgePath)
	eType1 <- rep(1.0, nEdge)  # 1.0: activation; -1.0: inhibition
	eScore <- matrix(0.0, nEdge, 2)  # scores of nodes on each edge
	#
	eS <- as.character(eType[,2])  # source nodes
	eT <- as.character(eType[,4])  # target nodes
	for (i in 1:nEdge){
		et <- edgePath[i,]
		# construct the eType1
		idxt <- which((eS %in% et[1]) & (eT %in% et[2]))
		if (length(idxt)>0){
			et1 <- as.character(eType[idxt[1], 6])
			# print(et1)
			if (et1 == 'Process(inhibition)'){
				eType1[i] <- -1.0
			}
		}
		# construct the eScore
		idxt <- which(gSym %in% et[1])
		if (length(idxt)>0){
			eScore[i,1] <- zScore[idxt]
		}
		#
		idxt <- which(gSym %in% et[2])
		if (length(idxt)>0){
			eScore[i,2] <- zScore[idxt]
		}
	}
	# the first node has no parent node
	eType1 <- c(1, eType1)
	# calculate the score of the signaling pathway
	for (i in 1:nEdge){
		vt <- 0.5*abs((eScore[i,1]+eType1[i+1]*eScore[i,2]))  # abs
		vScore <- vScore + vt
		#
		vt <- 0.5*(eScore[i,1]+eScore[i,2])  #
		aScore <- aScore + vt
	}
	vScore <- vScore + 0.5*eType1[i+1]*eScore[i,2] + 0.5*eScore[1,1]  # the last node and first node score
	vScore <- vScore/(nEdge+1)
	#
	aScore <- aScore + 0.5*eScore[i,2] + 0.5*eScore[1,1]  # the last node and first node score
	aScore <- aScore/(nEdge+1)
	# print('testing')
	# if (aScore != vScore) {print(c(aScore, vScore))}
	return(c(vScore,aScore))
}
