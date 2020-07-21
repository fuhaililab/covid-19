## ------------------------------------
#
###
#

#
###
#
## ------------------------------------
# venn_diagram_v1 <- function(){


# }

#
###
#
get_WikiPathways <- function(){
	library(rWikiPathways)
	
	pth_name_all <- listPathwayNames(organism = "Homo sapiens") 
	pth_id_all <- listPathwayIds(organism = "Homo sapiens") 
	pth_info_t1 <- getPathwayInfo(pth_id_all[397])
	pth_t1 <- getPathway(pth_id_all[397])

}

#
###
#
get_sim_GO <- function(GO_up_2){
	#
	library(org.Hs.eg.db)
	library(GO.db)
	library(GOSemSim)
	# GO-GO similarity network
	d <- godata('org.Hs.eg.db', ont="BP", computeIC=FALSE)
	n1 <- dim(GO_up_2)[1]
	e_go_up <- matrix(0.0, n1, n1)
	for (i in 1:n1){
		go_t1 <- GO_up_2[i,1]
		for (j in 1:n1){
			go_t2 <- GO_up_2[j,1]
			e_go_up[i,j] <- goSim(go_t1, go_t2, semData=d)
		}
	}
	colnames(e_go_up) <- as.character(GO_up_2[,1])
	print(dim(e_go_up))
	return(e_go_up)
}


#
###
#
generate_display_gene_go_network <- function(Go_t1, net_t1, symbol_1, fc_1){
	#
	node_t1 <- unique(net_t1[,1])
	node_t2 <- unique(net_t1[,2])
	n_t1 <- length(node_t1)
	n_t2 <- length(node_t2)
	att_t1 <- rep('GO-Term', n_t1)
	att_t2 <- rep('Gene', n_t2)
	#
	w_t1 <- rep(0, n_t1)
	name_t1 <- rep('test', n_t1)
	for (i in 1:n_t1){
		st1 <- node_t1[i]
		idxt <- which(Go_t1[,1] %in% st1)
		w_t1[i] <- Go_t1[idxt,5]
		name_t1[i] <- paste(Go_t1[idxt, 1], Go_t1[idxt, 2], sep='-')
	}
	w_t1 <- -log(as.numeric(w_t1))
	w_t1 <- (w_t1 - min(w_t1))/(max(w_t1) - min(w_t1))
	#
	w_t2 <- rep(1.01, n_t2)
	for (i in 1:n_t2){
		st1 <- node_t2[i]
		idxt <- which(symbol_1 %in% st1)
		if (length(idxt) < 1){ next }
		x_fct <- fc_1[idxt]
		w_t2[i] <- max(x_fct, 1.0/x_fct) 
		if (x_fct > log2(1.0)){
			att_t2[i] <- 'Gene_Up'
		} else if (x_fct < log2(1.0)){
			att_t2[i] <- 'Gene_Down'
		}
	}
	w_t2 <- (w_t2 - min(w_t2))/(max(w_t2) - min(w_t2))
	#
	node_att <- cbind(c(node_t1, node_t2), c(name_t1, node_t2), c(att_t1, att_t2), c(w_t1, w_t2))
	colnames(node_att) <- c('Node', 'Name1', 'Type1', 'Weight')

	return(node_att)
	# f_t1 <- paste(f_t, 'gene_GO_Up_Net.txt', sep='-')
	# colnames(net_t1) <- c('Source', 'Target')
	# write.table(net_t1, f_t1, quote=F, row.names=F, sep='\t')
	
	# f_t1 <- paste(dir_res, 'gene_GO_Up_Node.txt', sep='-')
	# write.table(node_att, f_t1, quote=F, row.names=F, sep='\t')

}
#
###
#
if (1 == 2){ # testing
	# x1 <- rnorm(100)
	# plot(x1)
	x1 <- c(-300:700)/100
	y1 <- exp(-x1^2)
	y2 <- exp(-(x1-3)^2)
	y3 <- 1/(1+exp(-(x1-2)))

	plot(x1,y1, type='l', col="blue", lwd=5)
	lines(x1, y2, type='l', col="red", lwd=5)
	lines(x1, y3, type='l', col="green", lwd=5)

	legend("topleft", c("gene_i", "gene_j", "gene_k"), fill=c("blue", "red", "green"))

}



getFcPvDEseq2 <- function(Exp, idx_control, idx_positive, symbol_1){
	# get fold change/ p-value using DESeq2
	library("DESeq2")

	count_t1 <- as.matrix(Exp[,c(idx_control, idx_positive)])
	# colnames(count_t1) <- c('control1', 'control2', 'control3', 'positive1', 'positive2', 'positive3')
	# col_data <- c('control', 'control', 'control', 'positive', 'positive', 'positive')

	n_st1 <- length(idx_control)
	n_st2 <- length(idx_positive)
	colnames(count_t1) <- c(paste('control', rep(1:n_st1), sep=''), paste('positivve', rep(1:n_st2), sep=''))
	col_data <- c(rep('control', n_st1), rep('positive', n_st2))

	col_data <- data.frame(col_data)
	colnames(col_data) <- 'condition'
	rownames(col_data) <- colnames(count_t1)
	dds <- DESeqDataSetFromMatrix(countData = count_t1, colData = col_data, design = ~condition)
	dds$condition <- factor(dds$condition, levels = c("control","positive"))
	# pre-filtering
	idx_t1 <- rowSums(counts(dds)) >= 30
	dds <- dds[idx_t1,]
	symbol_2 <- symbol_1[idx_t1]
	# differential expression analysis
	dds <- DESeq(dds)
	res <- results(dds) #@listData
	# res <- cbind(symbol_2, res[,c(2,5)])
	# print('new haha')
	res1 <- list()
	res1[[1]] <- res
	res1[[2]] <- symbol_2
	return(res1) # log2 folc change + pvalue

}

#
###
#
getFcPvTtest <- function(Exp1, idx_t1, idx_t2){
	# idx_t1 <- idx_IDH_wt
	# idx_t2 <- idx_IDH_mt
	v1 <- Exp1[, idx_t1]
	v2 <- Exp1[, idx_t2]
	# remove the constants values like 0 
	sd1 <- apply(v1, 1, sd)
	sd2 <- apply(v2, 1, sd)

	idx_1a <- which(sd1 == 0)
	# idx_1b <- which(sd1 >= 0.1*max(sd1))

	idx_2a <- which(sd2 == 0)
	# idx_2b <- which(sd2 >= 0.1*max(sd2))

	idx_3a <- union(idx_1a, idx_2a)
	# idx_3b <- union(idx_1b, idx_2b)

	# idx_4a <- union(idx_3a, idx_3b)
	idx_4a <- idx_3a
	idx_4b <- c(1:length(sd1))[-idx_4a]

	#
	fc1 <- rep(1.0, length(sd1))
	v1 <- Exp1[idx_4b, idx_t1]
	v2 <- Exp1[idx_4b, idx_t2]
	m1 <- rowMeans(v1)
	m2 <- rowMeans(v2)
	fc1[idx_4b] <- m1/m2

	pv1 <- rep(1.0, length(sd1))
	idxt1 <- c(1:dim(v1)[2])
 	idxt2 <- c(1:dim(v2)[2])+length(idxt1)
 	pv_t <- apply(cbind(v1, v2), 1, function(x) t.test(x[idxt1],x[idxt2])$p.value)  # p-val
 	pv1[idx_4b] <- pv_t

 	res <- cbind(fc1, pv1)

	return(res)

}

# sub-functions ---
tfEnrichAnalysisV1a <- function(gs1, all_gene, v_TF, v_Tar){
  #
  gs2 <- setdiff(all_gene, gs1)
  N_gs1 <- length(gs1)
  N_gs2 <- length(gs2)
  #
  gs_TF <- unique(v_TF)
  nt <- length(gs_TF)

  xt <- matrix(1.0, 2, 2)
  pv_t1 <- matrix(1.0, nt, 4) 
  for (i in 1:nt){
    # print(i)
    idx_t1 <- which(v_TF == gs_TF[i])
    gs_t1 <- v_Tar[idx_t1]
    # fisher test to identify the enriched GOs
    xt[1,1] <- length(intersect(gs1, gs_t1))
    xt[1,2] <- length(intersect(gs2, gs_t1))
    xt[2,1] <- N_gs1 - xt[1,1]
    xt[2,2] <- N_gs2 - xt[1,2]
    pv_t1[i,1] <- gs_TF[i]
    pv_t1[i,2] <- length(gs_t1)
    pv_t1[i,3] <- xt[1,1]
    pv_t1[i,4] <- fisher.test(xt, alternative = "greater")$p.value
  }
  return(pv_t1)
}

# sub-functions ---
keggEnrichAnalysisV1a <- function(gs1, all_gene, KeggPathways){
  # convert gene symbol to entrezID 
  gs2 <- setdiff(all_gene, gs1)
  #
  N_gs1 <- length(gs1)
  N_gs2 <- length(gs2)

  xt <- matrix(1.0, 2, 2)
  nt <- length(KeggPathways)
  pv_kegg_1 <- matrix(1.0, nt) 
  for (i in 1:nt){
    # print(i)
    st <- KeggPathways[[i]]
    gs_t1 <- union(st[,2], st[,4])
    # gs_t2 = unlist(mget(gs_t1,org.Hs.egSYMBOL))
    # fisher test to identify the enriched GOs
    xt[1,1] <- length(intersect(gs1, gs_t1))
    xt[1,2] <- length(intersect(gs2, gs_t1))
    xt[2,1] <- N_gs1 - xt[1,1]
    xt[2,2] <- N_gs2 - xt[1,2]
    pv_kegg_1[i] <- fisher.test(xt, alternative = "greater")$p.value
  }
  return(pv_kegg_1)
}
#
###
#
getGeneKeggNetworkV1a <- function(gs0, SPs, KeggPathways){
	#
	nt <- length(SPs)  # signaling-pathways 
	name_pathways <- names(KeggPathways)
	geneKeggNet <- c("Pathway", "Gene") 
	for (i in 1:nt){
		j <- SPs[i]
		st_name <- name_pathways[j]
		st <- KeggPathways[[j]]
    	gs_t1 <- union(st[,2], st[,4])
		#
		gs_t2 <- intersect(gs0, gs_t1)

		nt1 <- length(gs_t2)
		if (nt1 > 0){
			xt1 <- rep(st_name, nt1)
			xt2 <- gs_t2
			geneKeggNet <- rbind(geneKeggNet, cbind(xt1,xt2))
		} 

	}
	return(geneKeggNet[-1,])
	#
} 
#
###
#


# sub-functions ---
goEnrichAnalysisV1a <- function(gs1, all_gene){
  library(org.Hs.eg.db)
  library(GO.db)
  # convert gene symbol to entrezID 
  all_gene <- mapIds(org.Hs.eg.db, all_gene, 'ENTREZID', 'SYMBOL')
  idxt <- which(is.na(all_gene))
  if (length(idxt) > 0){
 	 all_gene <- all_gene[-idxt]
  }
  
  gs1 <- mapIds(org.Hs.eg.db, gs1, 'ENTREZID', 'SYMBOL')
  idxt <- which(is.na(gs1))
  if (length(idxt) > 0){
 	 gs1 <- gs1[-idxt]
  }

  gs2 <- setdiff(all_gene, gs1)
  # check all the GO terms related to these genes
  go_t1 <- mget(gs1, org.Hs.egGO)
  nt <- length(go_t1)
  ids_go_1 <- c()
  for (i in 1:nt){
    xt1 <- go_t1[[i]]
    ids_go_1 <- union(ids_go_1, names(xt1))
  }
  #
  N_gs1 <- length(gs1)
  N_gs2 <- length(gs2)

  xt <- matrix(1.0, 2, 2)
  nt <- length(ids_go_1)
  pv_go_1 <- matrix(1.0, nt, 4) 
  for (i in 1:nt){
    # print(i)
    go_idt <- ids_go_1[i]
    gs_t1 <- get(go_idt, org.Hs.egGO2ALLEGS)
    # gs_t2 = unlist(mget(gs_t1,org.Hs.egSYMBOL))
    # fisher test to identify the enriched GOs
    xt[1,1] <- length(intersect(gs1, gs_t1))
    xt[1,2] <- length(intersect(gs2, gs_t1))
    xt[2,1] <- N_gs1 - xt[1,1]
    xt[2,2] <- N_gs2 - xt[1,2]
    pv_go_1[i,3] <- length(gs_t1)
    pv_go_1[i,1] <- Term(GOTERM[[go_idt]])
    pv_go_1[i,2] <- Ontology(GOTERM[[go_idt]])
    pv_go_1[i,4] <- fisher.test(xt, alternative = "greater")$p.value
  }
  pv_go_1 <- cbind(ids_go_1, pv_go_1)
  #
  colnames(pv_go_1) <- c("go_id", "go_term_name", 'go_ontology', '#ofgenes_go', 'p-value')
  return(pv_go_1)
}
#
###
#
getGeneGoNetworkV1a <- function(gs0, GO1){
	#
	library(org.Hs.eg.db)
  	library(GO.db)
	
	gs1 <- mapIds(org.Hs.eg.db, gs0, 'ENTREZID', 'SYMBOL')
  	idxt <- which(is.na(gs1))
	if (length(idxt) > 0){
		gs1 <- gs1[-idxt]
		gs0 <-gs0[-idxt]
	}

	nt <- length(GO1)
	geneGoNet <- c("GO", "Gene") 
	for (i in 1:nt){
		go_idt <- GO1[i]
		gs_t1 <- get(go_idt, org.Hs.egGO2ALLEGS)
		#
		gs_t2 <- intersect(gs1, gs_t1)
		idxt1 <- which(gs1 %in% gs_t2)
		gs_t3 <- gs0[idxt1]

		nt1 <- length(gs_t2)
		if (nt1 > 0){
			xt1 <- rep(go_idt, nt1)
			xt2 <- gs_t3
			geneGoNet <- rbind(geneGoNet, cbind(xt1,xt2))
		} 

	}
	return(geneGoNet[-1,])
	#
} 
#
###
#
getKneighborNetwork <- function(f_net, k, gs){
	library(igraph)
	g_t <- read.delim(f_net, sep='\t')
	g_t <- cbind(as.character(g_t[,1]), as.character(g_t[,2]))
	g_t <- g_t[-which(g_t[,1] == g_t[,2]),]  # remove the self-connection
	gene_g_t <- union(g_t[,1], g_t[,2])  # 21699 genes;
	gTmp <- graph.edgelist(g_t, directed=F)
	#
	gs <- intersect(gs, gene_g_t)
	nt <- length(gs)
	net_1a <- c("Source", "Target")
	#
	for (i in 1:nt){ # from inhibition nodes to target;  
		paths_0a <- get.shortest.paths(gTmp, gs[i], setdiff(gs, gs[i-1]))
		# paths_0a <- get.shortest.paths(G_t, 'MAP3K14', Vi)
		paths_0b <- paths_0a$vpath
		n_t <- length(paths_0b)
		#
		idxt <- rep(0, n_t)
		for (it in 1:n_t){ # get lengths of each path
		  pt <- names(paths_0b[[it]])
		  idxt[it] <- length(pt)
		}
		idxt1 <- order(idxt, decreasing=F)
		idxt2 <- which(idxt[idxt1] >= 2)
		#
		for (it in 1:k){ # convert path to 2-column table/networ
		  pt <- names(paths_0b[[idxt1[idxt2[it]]]])
		  n_pt <- length(pt)
		  if (n_pt < 2){next}
		  for (l in 1:(n_pt-1)){
		    vt <- c(pt[l], pt[l+1])
		    net_1a <- rbind(net_1a, vt)
		  }
		}
	} # end of i_vi
	net_1a <- unique(net_1a[-1,])
	#
	g_1b <- induced.subgraph(gTmp, union(net_1a[,1], net_1a[,2]))
	net_1b <- get.edgelist(g_1b)
	return(net_1b)
	#
}

#
###
#
generateSigNetwork_gs_v1a <- function(f_bioGrid, gene_symbol_0, fc_t0, gs1a, gs2a){ # 
	#
	# gene_symbol_0 <- symbol_2
	# fc_t0 <- 2.0^fc_2
	# gs1a <- gs_t1c
	# gs2a <- gs_t1d
	#
	g_t <- read.delim(f_bioGrid, sep='\t')
	g_t <- cbind(as.character(g_t[,1]), as.character(g_t[,2]))
	g_t <- g_t[-which(g_t[,1] == g_t[,2]),]  # remove the self-connection
	gene_g_t <- union(g_t[,1], g_t[,2])  # 21699 genes;
	#
	library(igraph)
	ig_t0 <- graph.edgelist(g_t, directed=F)
	ig_t1 <- graph.edgelist(g_t, directed=F)
	#
	# set the edge weights 
	# case - 1: prefer/biase to up-regulated genes 
	idx_t <- which(gene_symbol_0 %in% gene_g_t)
	gene_symbol_2a <- gene_symbol_0[idx_t]
	fc_t2a <- fc_t0[idx_t]

	V0 <- gene_symbol_2a  # common genes 
	n_V0 <- length(V0)

	fc_t0a <- fc_t2a  # fold change values
	idx_t <- order(fc_t0a, decreasing=T)
	fc_t0a <- fc_t0a[idx_t]
	V0 <- V0[idx_t]

	fc_t0a[fc_t0a <= 0.5] <- 0.5  # ignore the down-regulated genes 
	fc_t0a[fc_t0a >= 5.0] <- 5.0  # trucate the extrem values  

	e_all <- E(ig_t1)
	e_head <- names(head_of(ig_t1, e_all))
	e_tail <- names(tail_of(ig_t1, e_all))
	#
	e_w1 <- rep(1.0, length(e_all))  # initialization 
	vht <- rep(1.0, length(e_all))
	vtt <- rep(1.0, length(e_all))
	#

	for (i in 1:n_V0){  # update the edge; assign the fold change to the nodes
		vt <- fc_t0a[i]
		idx_t <- which(e_head %in% V0[i])
		vht[idx_t] <- vt
		idx_t <- which(e_tail %in% V0[i])
		vtt[idx_t] <- vt
	}
	#
	e_w1 <- 2.0/(vht + vtt)
	# set_edge_attr("weight", value = 1:10)
	ig_t1 <- set.edge.attribute(ig_t1, "weight", value=e_w1)
	# gs1 <- gs1[which(gs1 %in% gene_g_t)]; #print(length(gs1))
	
	gs1 <- intersect(gs1a, gene_g_t) # root genes; remove genes are not in the bioGrid databse
	gs2 <- intersect(gs2a, gene_g_t) # remove genes are not in the
	n_stop <- round(0.7*length(gs1)) # only consider the first 75 nodes in each roots

	# loop for the gene set nodes
	g_net2A <- c('Source', 'Target') # initilization
	n_gs1 <- length(gs1)
	V0 <- gs1 # first generate the network among the gs1 genes.

	for (i_gs in 1:n_gs1){ # 
		print(i_gs)
		#
		g_net <- c('Source', 'Target')
		V2 <- gs1[i_gs]
		#
		idx_t <- which(V0 %in% V2)
		if (length(idx_t) >0){
			fc_t1 <- fc_t0[-idx_t]
			V1 <- V0[-idx_t]
		} else {
			V1 <- V0
		}
		# print('Start ---')
		flag1 <- 1
		while (flag1 > 0){
			gs_1 <- V1  
			# gs_1 <- intersect(gs_1, gene_g_t)
			# shortest path is limited to find the common paths 
			# paths_0v <- shortest.paths(ig_t0, V2, gs_1)
			paths_1v <- shortest.paths(ig_t1, V2, gs_1)		
			# v_t1 <- paths_1v/paths_0v
			v_t1 <- paths_1v
			
			idx_t1 <- which(v_t1 == min(v_t1))
			if (length(V2) == 1){			
				paths_0 <- get.shortest.paths(ig_t0, V2, gs_1[idx_t1[1]])
			} else {
				#
				idx_t1a <- ceiling(idx_t1[1]/length(V2))
				idx_t1b <- idx_t1[1] %% length(V2) 
				#
				paths_0 <- get.shortest.paths(ig_t0, V2[idx_t1b], gs_1[idx_t1a])
				# paths_1 <- get.shortest.paths(ig_t1, V2, gs_1)
			}
			paths_0 <- paths_0$vpath
			# paths_1 <- paths_1$vpath
			# choose the best path and add it into the network 
			pt <- names(paths_0[[1]])
			nPt <- length(pt)
			# convert to 2-column table/network 
			g_net_tmp <- c("tmp", "tmp")
			for (l in 1:(nPt-1)){
				vt <- c(pt[l], pt[l+1])
				g_net_tmp <- rbind(g_net_tmp, vt)
			}
			V_tmp <- union(g_net_tmp[,1], g_net_tmp[,2])
			V_tmp <- V_tmp[-which(V_tmp == "tmp")]
			g_net <- rbind(g_net, g_net_tmp[-1,])
			V2 <- union(V2, V_tmp)
			idx_t <- which(V0 %in% V2)
			fc_t1 <- fc_t0[-idx_t]
			V1 <- V0[-idx_t]
			# update network size
			g_net <- unique(g_net)
			# print(dim(g_net))

			flag1 <- (length(V1) >= n_stop)

			fc_V2 <- fc_t2a[which(gene_symbol_2a %in% V2)]
			# mean(fc_gs1)
			v_yt <- median(fc_V2)
			print(v_yt)

		}
		#
		if (length(g_net) > 3){
			g_net1 <- g_net[-1,]
		}
		g_net1 <- unique(g_net1)
		# remove (a,b)-(b,a) edge pairs 
		n_e <- dim(g_net1)[1]
		g_net2 <- g_net1[1,]
		for (i in 2:n_e){
			et <- g_net1[i,]
			net_ta <- rbind(g_net2, et)
			net_tb <- rbind(g_net2, c(et[2], et[1]))
			net_ta <- unique(net_ta)
			net_tb <- unique(net_tb)
			x_ta <- dim(net_ta)[1]
			x_tb <- dim(net_tb)[1]
			if (x_ta <= x_tb){
				g_net2 <- net_ta
			} else {
				g_net2 <- net_tb
			}
		}
		#
		gs_t0a <- union(g_net2[,1], g_net2[,2])
		fc_gs1 <- fc_t2a[which(gene_symbol_2a %in% gs_t0a)]
		# mean(fc_gs1)
		v_yt <- median(fc_gs1)
		print(v_yt)

		print(length(intersect(gs_t0a, gs1a)))

		# merge network starting from different nodes
		g_net2A <- rbind(g_net2A, g_net2)
		# g_net2A <- unique(g_net2A)
		
	} # end of loop gs1 nodes
	g_net2A <- g_net2A[-1,]

	g_net2A_a1 <- g_net2A 

	g_net2A <- g_net2A_a1
	# remove noisy (not much paths with low fold change) nodes 
	v_net2A <- union(g_net2A[,1], g_net2A[,2])
	xt1 <- table(g_net2A)
	node_t1 <- setdiff(v_net2A, gs1a) 
	xt2 <- xt1[names(xt1) %in% node_t1]
	xt2 <- xt2[order(xt2, decreasing=F)]

	n_t1 <- round(length(v_net2A) - length(gs1a)*1.3)

	node_rm <- intersect(names(xt2[xt2 <= round(length(gs1a)*0.5)]), names(xt2[1:n_t1]))
	n_rm <- length(node_rm)
	net_t <- unique(g_net2A)

	node_0a <- v_net2A

	for (i in 1:n_rm){
		node_st <- node_rm[i]
		node_0b <- setdiff(node_0a, node_st)

		e_t0 <- net_t[which(net_t[,1] != node_st & net_t[,2] != node_st),]
		node_0c <- union(e_t0[,1], e_t0[,2])

		if (length(setdiff(node_0b, node_0c)) > 0){next} # cannot remove this node

		node_0a <- node_0b
		v_5a <- union(e_t0[,1], e_t0[,2])
		G_t_5a <- graph.edgelist(e_t0, directed=F)
		dis_5a <- shortest.paths(G_t_5a, v_5a[1], v_5a)
		idx_t <- which(is.infinite(dis_5a))
		if (length(idx_t) == 0){
			print(i)
			net_t <- e_t0
		}
	}
	g_net2B <- net_t
	
	
	g_net2B <- g_net2A
	print('network - 1')
	print(dim(g_net2B))
	v_net2B <- union(g_net2B[,1], g_net2B[,2])
	
	print(length(v_net2B))
	fc_t_net <- fc_t2a[which(gene_symbol_2a %in% v_net2B)]
	mean(fc_t_net)
	median(fc_t_net)

	print(intersect(gs1a, v_net2B))
	# plot(G_t_5a, vertex.size=1, vertex.label.font=2, vertex.label.cex=0.5)

	fc_t_gs1 <- fc_t2a[which(gene_symbol_2a %in% gs1a)]
	mean(fc_t_gs1)
	median(fc_t_gs1)

	# connet the gs1 network to gs2 nodes
	g_net2A <- c('Source', 'Target') # initilization
	V2 <- v_net2B

	gs1 <- setdiff(gs2a, gs1a)
	gs1 <- intersect(gs1, gene_g_t)
	n_gs1 <- length(gs1)
	V0 <- gs1

	if (n_gs1 > 0){
		for (i_gs in 1:1){ # 
			#
			g_net <- c('Source', 'Target')
			#
			idx_t <- which(V0 %in% V2)
			if (length(idx_t) >0){
				fc_t1 <- fc_t0[-idx_t]
				V1 <- V0[-idx_t]
			} else {
				V1 <- V0
			}
			# print('Start ---')
			flag1 <- 1
			while (flag1 > 0){
				gs_1 <- V1  
				# gs_1 <- intersect(gs_1, gene_g_t)
				# shortest path is limited to find the common paths 
				# paths_0v <- shortest.paths(ig_t0, V2, gs_1)
				paths_1v <- shortest.paths(ig_t1, V2, gs_1)		
				# v_t1 <- paths_1v/paths_0v
				v_t1 <- paths_1v
				idx_t1 <- which(v_t1 == min(v_t1))
				if (length(V2) == 1){			
					paths_0 <- get.shortest.paths(ig_t0, V2, gs_1[idx_t1[1]])
				} else {
					#
					idx_t1a <- ceiling(idx_t1[1]/length(V2))
					idx_t1b <- idx_t1[1] %% length(V2) 
					#
					paths_0 <- get.shortest.paths(ig_t0, V2[idx_t1b], gs_1[idx_t1a])
					# paths_1 <- get.shortest.paths(ig_t1, V2, gs_1)
				}
				paths_0 <- paths_0$vpath
				# paths_1 <- paths_1$vpath
				# choose the best path and add it into the network 
				pt <- names(paths_0[[1]])
				nPt <- length(pt)
				# convert to 2-column table/network 
				g_net_tmp <- c("tmp", "tmp")
				for (l in 1:(nPt-1)){
					vt <- c(pt[l], pt[l+1])
					g_net_tmp <- rbind(g_net_tmp, vt)
				}
				V_tmp <- union(g_net_tmp[,1], g_net_tmp[,2])
				V_tmp <- V_tmp[-which(V_tmp == "tmp")]
				g_net <- rbind(g_net, g_net_tmp[-1,])
				V2 <- union(V2, V_tmp)
				idx_t <- which(V0 %in% V2)
				fc_t1 <- fc_t0[-idx_t]
				V1 <- V0[-idx_t]
				# update network size
				g_net <- unique(g_net)
				print(dim(g_net))

				flag1 <- length(V1)

			}
			#
			if (length(g_net) > 3){
				g_net1 <- g_net[-1,]
			}
			g_net1 <- unique(g_net1)
			# remove (a,b)-(b,a) edge pairs 
			n_e <- dim(g_net1)[1]
			g_net2 <- g_net1[1,]
			for (i in 2:n_e){
				et <- g_net1[i,]
				net_ta <- rbind(g_net2, et)
				net_tb <- rbind(g_net2, c(et[2], et[1]))
				net_ta <- unique(net_ta)
				net_tb <- unique(net_tb)
				x_ta <- dim(net_ta)[1]
				x_tb <- dim(net_tb)[1]
				if (x_ta <= x_tb){
					g_net2 <- net_ta
				} else {
					g_net2 <- net_tb
				}
			}
			# merge network starting from different nodes
			g_net2A <- rbind(g_net2A, g_net2)
			# g_net2A <- unique(g_net2A)
			
		} # end of loop gs1 nodes
		g_net2A <- g_net2A[-1,]
		#
	}
	#
	print('network - 2')
	print(dim(g_net2A))
	v_net2A <- union(g_net2A[,1], g_net2A[,2])
	
	print(length(v_net2A))
	fc_t_net <- fc_t2a[which(gene_symbol_2a %in% v_net2A)]
	mean(fc_t_net)
	median(fc_t_net)

	print(length(intersect(gs2a, v_net2A)))

	#
	g_net2C <- rbind(g_net2A, g_net2B)
	g_net1 <- unique(g_net2C)
	# remove (a,b)-(b,a) edge pairs 
	n_e <- dim(g_net1)[1]
	g_net2 <- g_net1[1,]
	for (i in 2:n_e){
		et <- g_net1[i,]
		net_ta <- rbind(g_net2, et)
		net_tb <- rbind(g_net2, c(et[2], et[1]))
		net_ta <- unique(net_ta)
		net_tb <- unique(net_tb)
		x_ta <- dim(net_ta)[1]
		x_tb <- dim(net_tb)[1]
		if (x_ta <= x_tb){
			g_net2 <- net_ta
		} else {
			g_net2 <- net_tb
		}
	}
	# merge network starting from different nodes
	g_net2C <- unique(g_net2)

	print('network - 2')
	print(dim(g_net2C))
	v_net <- union(g_net2C[,1], g_net2C[,2])
	print(length(v_net))
	fc_t_net <- fc_t2a[which(gene_symbol_2a %in% v_net)]
	mean(fc_t_net)
	median(fc_t_net)
	
	print(length(intersect(gs2a, v_net)))

	return(g_net2C)
}
#
###
#
generateSigNetworkV3a <- function(f_bioGrid, gene_symbol_0, fc_t0, gs1, net_score, net_size_max, net_size_min, dir_res1, net_name, wFlag){  # network analysis
	# parameters --- 
	# g_t: network matrix (2-column), like biogrid from gene A to gene B 
	# gene_symbol_0: the gene symbols of gene expression files  
	# fc_1: fold change, or other importance score of genes in the gene_symbol_0
	# net_score: to limit the size of the network by score of the network: sum(v)-sum(e)
	# dir_res1: is the folder to save the results
	# net_name: give the name of the network to save 
	#
	# print('Prepare to Start ---')
	# ---- function starts ------
	library(igraph)
	# f1 <- paste(dir_data0, 'BIOGRID-ALL-3.5.174.mitab.Symbol.txt', sep='/')
	# g_t <- read.delim(f1, sep='\t')
	# g_t <- cbind(as.character(g_t[,1]), as.character(g_t[,2]))
	#
	g_t <- read.delim(f_bioGrid, sep='\t')
	g_t <- cbind(as.character(g_t[,1]), as.character(g_t[,2]))
	g_t <- g_t[-which(g_t[,1] == g_t[,2]),]  # remove the self-connection
	gene_g_t <- union(g_t[,1], g_t[,2])  # 21699 genes;
	#
	library(igraph)
	ig_t0 <- graph.edgelist(g_t, directed=F)
	ig_t1 <- graph.edgelist(g_t, directed=F)
	#
	idx_t <- which(gene_symbol_0 %in% gene_g_t)
	V0 <- gene_symbol_0[idx_t]  # common genes 
	n_V0 <- length(V0)
	fc_t0 <- fc_1[idx_t]  # fold change values
	idx_t <- order(fc_t0, decreasing=T)
	fc_t0 <- fc_t0[idx_t]
	V0 <- V0[idx_t]
	# set the edge weights 
	# case - 1: prefer/biase to up-regulated genes 
	e_all <- E(ig_t1)
	e_head <- names(head_of(ig_t1, e_all))
	e_tail <- names(tail_of(ig_t1, e_all))
	#
	e_w1 <- rep(1.0, length(e_all))  # initialization 
	vht <- rep(1.0, length(e_all))
	vtt <- rep(1.0, length(e_all))
	#
	fc_t0a <- fc_t0 
	fc_t0a[fc_t0a <= 1.0] <- 1.0  # ignore the down-regulated genes 
	for (i in 1:n_V0){  # update the edge; assign the fold change to the nodes
		vt <- fc_t0a[i]
		idx_t <- which(e_head %in% V0[i])
		vht[idx_t] <- vt
		idx_t <- which(e_tail %in% V0[i])
		vtt[idx_t] <- vt
	}
	#
	e_w1 <- 2.0/(vht + vtt)
	# set_edge_attr("weight", value = 1:10)
	ig_t1 <- set.edge.attribute(ig_t1, "weight", value=e_w1)
	#
	gs1 <- gs1[which(gs1 %in% gene_g_t)]; #print(length(gs1))
	#
	###
	# loop for the gene set nodes
	g_net2A <- c('Source', 'Target') # initilization
	n_gs1 <- length(gs1)
	for (i_gs in 1:n_gs1){ # 
		N_e <- 1
		nt0 <- 200  # top 10 nodes to reach based on fold change 
		netScore <- 2.0
		w_node_t <- 1.2
		#
		g_net <- c('Source', 'Target')
		V2 <- gs1[i_gs]
		#
		idx_t <- which(V0 %in% V2)
		if (length(idx_t) >0){
			fc_t1 <- fc_t0[-idx_t]
			V1 <- V0[-idx_t]
		} else {
			V1 <- V0
		}
		# print('Start ---')
		# while (N_e < net_size){  # network size limit; at least > net_size_min
		while ((netScore >= net_score & w_node_t >= 1.10 & N_e < net_size_max) | N_e < net_size_min){
		# for (i in 1:50){
			# print('Started ---')
			# fetch the top nt0 genes from V1 based on fold change 
			gs_1 <- V1[1:nt0]  # gene set 1
			gs_1 <- intersect(gs_1, gene_g_t)
			# shortest path is limited to find the common paths 
			paths_0v <- shortest.paths(ig_t0, V2, gs_1)
			paths_1v <- shortest.paths(ig_t1, V2, gs_1)		
			#
			v_t1 <- paths_1v/paths_0v
			idx_t1 <- which(v_t1 == min(v_t1))
			if (length(V2) == 1){			
				paths_0 <- get.shortest.paths(ig_t0, V2, gs_1[idx_t1[1]])
			} else {
				#
				idx_t1a <- ceiling(idx_t1[1]/length(V2))
				idx_t1b <- idx_t1[1] %% length(V2) 
				#
				paths_0 <- get.shortest.paths(ig_t0, V2[idx_t1b], gs_1[idx_t1a])
				# paths_1 <- get.shortest.paths(ig_t1, V2, gs_1)
			}
			paths_0 <- paths_0$vpath
			# paths_1 <- paths_1$vpath
			# choose the best path and add it into the network 
			pt <- names(paths_0[[1]])
			nPt <- length(pt)
			# convert to 2-column table/network 
			g_net_tmp <- c("tmp", "tmp")
			for (l in 1:(nPt-1)){
				vt <- c(pt[l], pt[l+1])
				g_net_tmp <- rbind(g_net_tmp, vt)
			}
			V_tmp <- union(g_net_tmp[,1], g_net_tmp[,2])
			V_tmp <- V_tmp[-which(V_tmp == "tmp")]
			# V_tmp <- setdiff(V_tmp, V2)
			#
			# idx_t <- which(V0 %in% V_tmp)
			# if (length(idx_t) > 0){
			# 	w_node_t <- fc_t0[idx_t] # multiple probes
			# 	w_node_t <- mean(w_node_t) 
			# 	print('w_node_t')
			# 	print(w_node_t)
			# 	if (w_node_t < 1.01){break} # exit while loop
			# } # check if new nodes has NO information already
			#
			g_net <- rbind(g_net, g_net_tmp[-1,])
			# print(g_net)
			# update V1, V2
			# V2 <- union(g_net[,1], g_net[,2])
			V2 <- union(V2, V_tmp)
			idx_t <- which(V0 %in% V2)
			fc_t1 <- fc_t0[-idx_t]
			V1 <- V0[-idx_t]
			# update network size
			g_net <- unique(g_net)
			N_e <- (dim(g_net)[1] -1)
			# sum fold change 
			idx_t <- which(V0 %in% V2)
			netScore <- fc_t0[idx_t] # multiple probes 
			netScore <- mean(netScore) # average node score
			# print('netScore')
			# print('Network_size:')
			# print(N_e)
		}
		#
		# print(dim(g_net))
		#
		if (length(g_net) > 3){
			g_net1 <- g_net[-1,]
		}
		g_net1 <- unique(g_net1)
		# remove (a,b)-(b,a) edge pairs 
		n_e <- dim(g_net1)[1]
		g_net2 <- g_net1[1,]
		for (i in 2:n_e){
			et <- g_net1[i,]
			net_ta <- rbind(g_net2, et)
			net_tb <- rbind(g_net2, c(et[2], et[1]))
			net_ta <- unique(net_ta)
			net_tb <- unique(net_tb)
			x_ta <- dim(net_ta)[1]
			x_tb <- dim(net_tb)[1]
			if (x_ta <= x_tb){
				g_net2 <- net_ta
			} else {
				g_net2 <- net_tb
			}
		}
		# merge network starting from different nodes
		g_net2A <- rbind(g_net2A, g_net2)
		g_net2A <- unique(g_net2A)
		# print('---')
		# print(dim(g_net2A))
		
	} # end of loop gs1 nodes
	g_net2A <- g_net2A[-1,] 
	#
	# write out the network
	if (wFlag == TRUE){	
	#
		# print('signaling network are being saved!')	
		# colnames(g_net2A) <- c('Source', 'Target')
		f1 <- paste(dir_res1, '/sigNet_', net_name, '_', format(Sys.time(), "%Y-%b-%d-%X"), '.txt', sep='')
		write.table(g_net2A, f1, quote=F, row.names=F, sep='\t')
		# attributes 
		nodes <- union(g_net2A[,1], g_net2A[,2])
		node_w1 <- rep(1.0, length(nodes))
		for (i in 1:length(nodes)){
			idx_t <- which(V0 %in% nodes[i])
			if (length(idx_t) > 0){
				node_w1[i] <- fc_t0[idx_t[1]] # multiple probes 
			}
		}
		node_w2 <- node_w1
		# normalize node weights to [1, 10]
		idx_t <- order(node_w1, decreasing=T)
		v_max <- node_w1[idx_t[5]]
		idx_t <- order(node_w1, decreasing=F)
		v_min <- node_w1[idx_t[5]]
		node_w1[node_w1 > v_max] <- v_max
		node_w1[node_w1 < v_min] <- v_min
		node_w1 <- 1 + 9*(node_w1 - v_min)/(v_max - v_min)
		#
		att_1 <- cbind(nodes, node_w1)
		colnames(att_1) <- c('Nodes', 'Weight')
		# f2 <- paste(dir_res1, 'sigNet_att.txt', sep='/')
		f2 <- gsub('.txt', '_att.txt', f1)
		write.table(att_1, f2, quote=F, row.names=F, sep='\t')
		#
		att_2 <- cbind(nodes, node_w2)
		colnames(att_2) <- c('Nodes', 'Weight')
		# f2 <- paste(dir_res1, 'sigNet_att.txt', sep='/')
		f3 <- gsub('.txt', '_att1.txt', f1)
		write.table(att_2, f3, quote=F, row.names=F, sep='\t')

	}
	# return variables 
	return(g_net2A)
}
#
###
#

getNormalizedLaplacian <- function(A){
    D = colSums(A) # degrees of vertices
    n1 = dim(A)[1]
    D1 = diag(1.0/sqrt(D))
    L = diag(n1) - D1 %*% A %*% D1
    return(L)
 }




