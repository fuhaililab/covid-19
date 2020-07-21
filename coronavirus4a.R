# No 'ACE2' gene in the KEGG pathway
if (1 == 2){
	dir_data0 <- c('~/FHLwustl/Projects/OrangeLi/DataSets/')
	f_bioGrid <- paste(dir_data0, 'BIOGRID-ALL-3.5.174.mitab.Symbol.txt', sep='/')
	f_tarDrugBank <- paste(dir_data0, 'drug_tar_drugBank_all.txt', sep='/')
	dir_project <- c('~/FHLwustl/Projects/covid')
	dir_res <- paste(dir_project, 'res', sep='/')
	# f_Exp <- paste(dir_project, 'GSE147507RawReadCounts.tsv', sep='/')

	f_Exp <- paste(dir_project, 'GSE147507_RawReadCounts_Human.tsv', sep='/')
	f_ppi_1 <- paste(dir_project, 'PPI/41586_2020_2286_MOESM6_ESM.txt', sep='/')


	# sub-funcs
	dir_code <- c('~/FHLwustl/Projects/tcga/code')
	source(paste(dir_code, 'allSubFuncsV0.R', sep='/'))

	if (1 == 2){ # save the workspace data
		#
		f_save_1 <- paste(dir_res, 'covid_v4a.RData', sep='/')
		save.image(f_save_1)
	}

}

if (1 == 2){ # data load and pre-analysis
	Exp <- read.delim(f_Exp, sep='\t')
	exp_symbol_0 <- as.character(Exp[,1])
	Exp <- Exp[,-1]
	sample_name_0 <- names(Exp)

	idx_mock_A549 <- c(7, 8, 9)
	idx_cov2_A549 <- c(10, 11, 12)
	idx_mock_A549_ACE2_s6 <- c(27, 28, 29)
	idx_cov2_A549_ACE2_s6 <- c(30, 31, 32)
	idx_mock_A549_ACE2_s16 <- c(70, 71, 72)
	idx_cov2_A549_ACE2_s16 <- c(73, 74, 75)

	idx_mock_NHBE <- c(1, 2, 3)
	idx_cov2_NHBE <- c(4, 5, 6)

	idx_mock_CALU3_s7 <- c(33, 34, 35)
	idx_cov2_CALU3_s7 <- c(36, 37, 38)

	idx_mock_LungBio_s15 <- c(66, 67)
	idx_cov2_LungBio_s15 <- c(68, 69)


	# fold change/ p-value using DEseq2
	fc_pv_LungBio_s15_cov2_vs_mock <- getFcPvDEseq2(Exp, idx_mock_LungBio_s15, idx_cov2_LungBio_s15, exp_symbol_0)
	fc_pv_A549_cov2_vs_mock <- getFcPvDEseq2(Exp, idx_mock_A549, idx_cov2_A549, exp_symbol_0)
	fc_pv_NHBE_cov2_vs_mock <- getFcPvDEseq2(Exp, idx_mock_NHBE, idx_cov2_NHBE, exp_symbol_0)
	fc_pv_CALU3_s7_cov2_vs_mock <- getFcPvDEseq2(Exp, idx_mock_CALU3_s7, idx_cov2_CALU3_s7, exp_symbol_0)
	fc_pv_A549_ACE2_s6_cov2_vs_mock <- getFcPvDEseq2(Exp, idx_mock_A549_ACE2_s6, idx_cov2_A549_ACE2_s6, exp_symbol_0)
	fc_pv_A549_ACE2_s16_cov2_vs_mock <- getFcPvDEseq2(Exp, idx_mock_A549_ACE2_s16, idx_cov2_A549_ACE2_s16, exp_symbol_0)

	symbol_A549_ACE2_s6 <- fc_pv_A549_ACE2_s6_cov2_vs_mock[[2]]
	res_DEseq2 <- fc_pv_A549_ACE2_s6_cov2_vs_mock[[1]]  
	fc_A549_ACE2_s6_cov2_vs_mock <- res_DEseq2[,2]; 	pv_A549_ACE2_s6_cov2_vs_mock <- res_DEseq2[,5]

	symbol_A549_ACE2_s16 <- fc_pv_A549_ACE2_s16_cov2_vs_mock[[2]]
	res_DEseq2 <- fc_pv_A549_ACE2_s16_cov2_vs_mock[[1]]  
	fc_A549_ACE2_s16_cov2_vs_mock <- res_DEseq2[,2]; 	pv_A549_ACE2_s16_cov2_vs_mock <- res_DEseq2[,5]

	symbol_A549 <- fc_pv_A549_cov2_vs_mock[[2]]
	res_DEseq2 <- fc_pv_A549_cov2_vs_mock[[1]]  
	fc_A549_cov2_vs_mock <- res_DEseq2[,2];	pv_A549_cov2_vs_mock <- res_DEseq2[,5]

	symbol_NHBE <- fc_pv_NHBE_cov2_vs_mock[[2]]
	res_DEseq2 <- fc_pv_NHBE_cov2_vs_mock[[1]]  
	fc_NHBE_cov2_vs_mock <- res_DEseq2[,2]; 	pv_NHBE_cov2_vs_mock <- res_DEseq2[,5]

	symbol_CALU3_s7 <- fc_pv_CALU3_s7_cov2_vs_mock[[2]]
	res_DEseq2 <- fc_pv_CALU3_s7_cov2_vs_mock[[1]]  
	fc_CALU3_s7_cov2_vs_mock <- res_DEseq2[,2]; 	pv_CALU3_s7_cov2_vs_mock <- res_DEseq2[,5]
	
	symbol_LungBio_s15 <- fc_pv_LungBio_s15_cov2_vs_mock[[2]]
	res_DEseq2 <- fc_pv_LungBio_s15_cov2_vs_mock[[1]]  
	fc_LungBio_s15_cov2_vs_mock <- res_DEseq2[,2]; 	pv_LungBio_s15_cov2_vs_mock <- res_DEseq2[,5]

	if (1 == 2){ # DEG comparison
		gs_up_NHBE <- symbol_NHBE[which(fc_NHBE_cov2_vs_mock >= log2(1.25)  & pv_NHBE_cov2_vs_mock <= 0.05)]
		gs_up_CALU3_s7 <- symbol_CALU3_s7[which(fc_CALU3_s7_cov2_vs_mock >= log2(1.5)  & pv_CALU3_s7_cov2_vs_mock <= 0.05)]
		gs_up_LungBio_s15 <- symbol_LungBio_s15[which(fc_LungBio_s15_cov2_vs_mock >= log2(1.05)  & pv_LungBio_s15_cov2_vs_mock <= 0.05)]
		gs_up_A549_ACE2_s6 <- symbol_A549_ACE2_s6[which(fc_A549_ACE2_s6_cov2_vs_mock >= log2(1.5) & pv_A549_ACE2_s6_cov2_vs_mock <= 0.05)] 
		gs_up_A549_ACE2_s16 <- symbol_A549_ACE2_s16[which(fc_A549_ACE2_s16_cov2_vs_mock >= log2(1.5) & pv_A549_ACE2_s16_cov2_vs_mock <= 0.05)] 
		length(gs_up_NHBE)
		length(gs_up_CALU3_s7)
		length(gs_up_LungBio_s15)	
		length(gs_up_A549_ACE2_s6)
		length(gs_up_A549_ACE2_s16)
		length(intersect(gs_up_NHBE, gs_up_CALU3_s7))
		length(intersect(gs_up_NHBE, gs_up_LungBio_s15))
		length(intersect(gs_up_NHBE, gs_up_A549_ACE2_s6))
		length(intersect(gs_up_A549_ACE2_s6, gs_up_A549_ACE2_s16))
		length(intersect(gs_up_A549_ACE2_s6, gs_up_CALU3_s7))

		gs_dn_NHBE <- symbol_NHBE[which(fc_NHBE_cov2_vs_mock <= log2(1.0/1.25)  & pv_NHBE_cov2_vs_mock <= 0.05)]
		gs_dn_A549_ACE2_s6 <- symbol_A549_ACE2_s6[which(fc_A549_ACE2_s6_cov2_vs_mock <= log2(1.0/1.5) & pv_A549_ACE2_s6_cov2_vs_mock <= 0.05)] 
		gs_dn_A549_ACE2_s16 <- symbol_A549_ACE2_s16[which(fc_A549_ACE2_s16_cov2_vs_mock <= log2(1.0/1.5) & pv_A549_ACE2_s16_cov2_vs_mock <= 0.05)] 
		length(gs_dn_NHBE)
		length(gs_dn_A549_ACE2_s6)
		length(gs_dn_A549_ACE2_s16)

	}

}

if (1 == 2){ # KEGG signaling network analysis
	symbol_t1 <- symbol_NHBE
	fc_t1 <- fc_NHBE_cov2_vs_mock
	str_t <- 'NHBE-sd125-1-'
	length(symbol_t1)
	length(fc_t1)

	symbol_t1 <- symbol_A549_ACE2_s6
	fc_t1 <- fc_A549_ACE2_s6_cov2_vs_mock
	str_t <- 'A549_ACE2_s6-sd125-1-'
	length(symbol_t1)
	length(fc_t1)

	symbol_t1 <- symbol_CALU3_s7
	fc_t1 <- fc_CALU3_s7_cov2_vs_mock
	str_t <- 'CALU3_s7-sd125-1-'
	length(symbol_t1)
	length(fc_t1)

	symbol_t1 <- symbol_CALU3_s7
	x <- c(1:length(symbol_t1))
	x1 <- sample(x, length(symbol_t1), replace = FALSE)
	fc_t1 <- fc_CALU3_s7_cov2_vs_mock[x1]
	str_t <- 'Random-1-'
	length(symbol_t1)
	length(fc_t1)

	symbol_t1 <- symbol_CALU3_s7
	x <- c(1:length(symbol_t1))
	x1 <- sample(x, length(symbol_t1), replace = FALSE)
	fc_t1 <- fc_CALU3_s7_cov2_vs_mock[x1]
	str_t <- 'Random-2-'
	length(symbol_t1)
	length(fc_t1)

	# network comparison
	str_t <- 'NHBE-sd125-1-'
	f1 <- paste(dir_res, paste(str_t,'sigNet-Activation-Att.txt',sep=''), sep='/')
	node_info1 <- read.delim(f1, sep='\t')
	f1 <- paste(dir_res, paste(str_t,'sigNet-Activation.txt',sep=''), sep='/')
	net_t1 <- read.delim(f1, sep='\t')

	str_t <- 'A549_ACE2_s6-sd125-1-'
	f1 <- paste(dir_res, paste(str_t,'sigNet-Activation-Att.txt',sep=''), sep='/')
	node_info2 <- read.delim(f1, sep='\t')

	str_t <- 'CALU3_s7-sd125-1-'
	f1 <- paste(dir_res, paste(str_t,'sigNet-Activation-Att.txt',sep=''), sep='/')
	node_info3 <- read.delim(f1, sep='\t')

	node_t1_2 <- intersect(node_info1[,1], node_info2[,1])
	node_t2_3 <- intersect(node_info2[,1], node_info3[,1])
	node_t1_3 <- intersect(node_info1[,1], node_info3[,1])
	node_t1_2_3 <- intersect(node_t1_2, node_t1_3)
	node_t1_23 <- union(node_t1_2, node_t1_3)

	dim(node_info1)
	dim(node_info2)
	dim(node_info3)
	length(node_t1_2)
	length(node_t1_3)
	length(node_t2_3)
	length(node_t1_2_3)

	# sub-network
	str_t1 <- 'NHBE-sd125-1-'
	net_t1a <- net_t1[net_t1[,1] %in% node_t1_2_3 & net_t1[,2] %in% node_t1_2_3, ]
	node_t1a <- union(net_t1a[,1], net_t1a[,2])
	length(node_t1a)
	f1 <- paste(dir_res, paste(str_t1,'sigNet-Activation-sub1.txt',sep=''), sep='/')
	write.table(net_t1a, f1, quote=F, col.names=T, row.names=F, sep='\t')

	net_t1b <- net_t1[net_t1[,1] %in% node_t1_23 & net_t1[,2] %in% node_t1_23, ]
	node_t1b <- union(net_t1a[,1], net_t1a[,2])
	length(node_t1b)
	f1 <- paste(dir_res, paste(str_t1,'sigNet-Activation-sub2.txt',sep=''), sep='/')
	write.table(net_t1b, f1, quote=F, col.names=T, row.names=F, sep='\t')

	xt1 <- fc_t1[which(symbol_t1 %in% node_t1_2_3)]
	mean(xt1)
	median(xt1)
	xt1 <- fc_t1[which(symbol_t1 %in% "SMAD1")]

}

if (1 == 2){ # drugBank drugs inhibiting the signaling pathways
	str_t1 <- 'NHBE-sd125-1-'
	f1 <- paste(dir_res, paste(str_t1,'sigNet-Activation-sub2.txt',sep=''), sep='/')
	net_xt1 <- read.delim(f1, sep='\t') #344 edges
	node_xt1 <-union(net_xt1[,1], net_xt1[,2]) #244 nodes
	#
	if (1 == 2){			# compare KEGG signaling and GO
		gs_go_t1 <- unique(net_super_GO_gene_up_2[,2])
		gs_kegg_t1 <- node_xt1
		gs_go_kegg_t1 <- intersect(gs_go_t1, gs_kegg_t1)
	}

	f_t1 <- paste(dir_data0, 'drug links fda.csv', sep='/')
	d_fda <- read.delim(f_t1, sep=',')
	d_fda <- as.character(d_fda[,2]) # 3556

	f_t1 <- paste(dir_data0, 'DrugBank-Drug-Target.txt', sep='/')
	d_xt1 <- read.delim(f_t1, sep='\t')
	idx_xt1 <- which(d_xt1[,2] %in% node_xt1)
	d_xt2 <- d_xt1[idx_xt1,]
	d_xt3 <- unique(as.character(d_xt2[,1]))
	d_xt4 <- intersect(d_xt3, toupper(d_fda))
	idx_xt2 <- which(d_xt2[,1] %in% d_xt4)
	d_xt2a <- d_xt2[idx_xt2,]
	d_xt2a[,1] <- tolower(d_xt2a[,1])

	table_xt1 <- table(as.character(d_xt2a[,2]))
	if (1 == 2){ #
		t1 <- sort(table_xt1)
		tar_t1 <- names(t1)
		n_t1 <- length(tar_t1)
		x_t1 <- matrix('', n_t1, 100)
		z_t1 <- as.character(d_xt2a[,1])
		z_t2 <- as.character(d_xt2a[,2])
		for (i in 1:n_t1){
			st <- tar_t1[i]
			y_t1 <- which(z_t2 %in% st)
			n_yt1 <- length(y_t1)
			x_t1[i, 1:n_yt1] <- z_t1[y_t1]
		}
		x_t2 <- cbind(tar_t1, x_t1)
		#
		f_t1 <- paste(dir_res, 'target_drug_Tabel_1.txt', sep='/')
		write.table(x_t2, f_t1, quote=F, sep='\t', row.names=F)

		f_t1 <- paste(dir_res, 'target_drug_Tabel_2.txt', sep='/')
		write.table(x_t1, f_t1, quote=F, sep=',', row.names=F)

	}

	# generate drug-target signaling network
	net_xt2 <- rbind(as.matrix(net_xt1[,c(1,2)]), as.matrix(d_xt2a))
	f_t1 <- paste(dir_res, 'Drug_signaling_net-up-Activation.txt', sep='/')
	write.table(net_xt2, f_t1, quote=F, sep='\t', row.names=F)

	# node attributes
	z_t0 <- unique(d_xt2a[,2])
	z_t1 <- tolower(unique(d_xt2a[,1]))
	node_yt1 <- c(node_xt1, z_t1)
	type_yt1 <- c(rep('Gene', length(node_xt1)), rep('Drug', length(z_t1)))
	isTarget <- rep(0, length(type_yt1))
	isTarget[which(node_yt1 %in% z_t0)] <- 1
	node_att_1a <- cbind(node_yt1, type_yt1, isTarget)

	f_t1 <- paste(dir_res, 'Drug_signaling_Node-up-Activation.txt', sep='/')
	write.table(node_att_1a, f_t1, quote=F, sep='\t', row.names=F)


}

if (1 == 2){ # GO analysis
	symbol_t1 <- symbol_NHBE
	fc_t1 <- fc_NHBE_cov2_vs_mock
	pv_t1 <- pv_NHBE_cov2_vs_mock
	str_t <- 'NHBE-GO-'
	length(symbol_t1)
	length(fc_t1)
	length(pv_t1)

	gs_up_1 <- symbol_t1[which(fc_t1 >= log2(1.25)  & pv_t1 <= 0.05)]
	gs_dn_1 <- symbol_t1[which(fc_t1 <= -log2(1.25) & pv_t1 <= 0.05)]
	length(gs_up_1)
	length(gs_dn_1)


	symbol_t1 <- symbol_A549_ACE2_s6
	fc_t1 <- fc_A549_ACE2_s6_cov2_vs_mock
	pv_t1 <- pv_A549_ACE2_s6_cov2_vs_mock
	str_t <- 'A549_ACE2_s6-GO-'
	length(symbol_t1)
	length(fc_t1)
	length(pv_t1)

	symbol_t1 <- symbol_CALU3_s7
	fc_t1 <- fc_CALU3_s7_cov2_vs_mock
	pv_t1 <- pv_CALU3_s7_cov2_vs_mock
	str_t <- 'CALU3_s7-GO-'
	length(symbol_t1)
	length(fc_t1)
	length(pv_t1)

	###
	gs_up_1 <- symbol_t1[which(fc_t1 >= log2(2.0)  & pv_t1 <= 0.05)]
	gs_dn_1 <- symbol_t1[which(fc_t1 <= -log2(2.0) & pv_t1 <= 0.05)]
	length(gs_up_1)
	length(gs_dn_1)

	GO_up_1 <- goEnrichAnalysisV1a(gs_up_1, symbol_t1)
	GO_dn_1 <- goEnrichAnalysisV1a(gs_dn_1, symbol_t1)
	dim(GO_up_1)
	dim(GO_dn_1)

	idx_up1 <- which(as.character(GO_up_1[,3]) == "BP" & as.numeric(GO_up_1[,4]) <= 1000 & as.numeric(GO_up_1[,4]) >= 10 & GO_up_1[,5] <= 0.025)
	idx_dn1 <- which(as.character(GO_dn_1[,3]) == "BP" & as.numeric(GO_dn_1[,4]) <= 1000 & as.numeric(GO_dn_1[,4]) >= 10 & GO_dn_1[,5] <= 0.025)

	GO_up_2 <- GO_up_1[idx_up1,]
	GO_dn_2 <- GO_dn_1[idx_dn1,]

	dim(GO_up_2)
	dim(GO_dn_2)

	flag1 <- paste(str_t, 'fc2', sep='')
	f_t1 <- paste(dir_res, paste(flag1, 'list_GO_up_2.txt', sep='-'), sep='/')
	write.table(GO_up_2, f_t1, quote=F, sep='\t', row.names=F)

	f_t2 <- paste(dir_res, paste(flag1, 'list_GO_dn_2.txt', sep='-'), sep='/')
	write.table(GO_dn_2, f_t2, quote=F, sep='\t', row.names=F)

	if (1 == 2){ # identify the common activated/inhibited GOs
		str_t <- 'NHBE-GO-p25-'
		flag1 <- paste(str_t, 'fc125', sep='')
		f_t1 <- paste(dir_res, paste(flag1, 'list_GO_up_2.txt', sep='-'), sep='/')
		go_t1 <- read.delim(f_t1, sep='\t')

		str_t <- 'A549_ACE2_s6-GO-p25-'
		flag1 <- paste(str_t, 'fc2', sep='')
		f_t1 <- paste(dir_res, paste(flag1, 'list_GO_up_2.txt', sep='-'), sep='/')
		go_t2 <- read.delim(f_t1, sep='\t')

		str_t <- 'CALU3_s7-GO-p25-'
		flag1 <- paste(str_t, 'fc2', sep='')
		f_t1 <- paste(dir_res, paste(flag1, 'list_GO_up_2.txt', sep='-'), sep='/')
		go_t3 <- read.delim(f_t1, sep='\t')

		dim(go_t1)
		dim(go_t2)
		dim(go_t3)

		idx_t1_2 <- which(go_t1[,1] %in% go_t2[,1])
		idx_t1_3 <- which(go_t1[,1] %in% go_t3[,1])
		idx_t2_3 <- which(go_t2[,1] %in% go_t3[,1])
		length(idx_t1_2)
		length(idx_t1_3)
		length(idx_t2_3)

		idx_t1 <- union(idx_t1_2,idx_t1_3)
		go_t1a <- go_t1[idx_t1,]
		flag1 <- paste(str_t, 'fc125', sep='')
		f_t1 <- paste(dir_res, paste(flag1, 'list_GO_up_2_sub2_3.txt', sep='-'), sep='/')
		write.table(go_t1a, f_t1, quote=F, sep='\t', row.names=F)

		# Manual selection of the activated GOs
		f_t1a <- paste(dir_res, paste(flag1, 'list_GO_up_2_sub2_3_manual_selection.txt', sep='-'), sep='/')
		GO_up_2a <- read.delim(f_t1a, sep='\t')
		x_t1 <- GO_up_2a[,3]
		idx_t1 <- which(x_t1 > 0)
		GO_up_2a <- GO_up_2a[-idx_t1,-3]
		dim(GO_up_2a)

	}

	#
	symbol_t1 <- symbol_NHBE
	fc_t1 <- fc_NHBE_cov2_vs_mock
	pv_t1 <- pv_NHBE_cov2_vs_mock
	str_t <- 'NHBE-GO-'
	length(symbol_t1)
	length(fc_t1)
	length(pv_t1)

	gs_up_1 <- symbol_t1[which(fc_t1 >= log2(1.25)  & pv_t1 <= 0.05)]
	gs_dn_1 <- symbol_t1[which(fc_t1 <= -log2(1.25) & pv_t1 <= 0.05)]
	length(gs_up_1)
	length(gs_dn_1)
	# take the maximum or minimum value of multiple probes (the same gene)
	symbol_2 <- unique(symbol_t1)
	nt <- length(symbol_2)
	fc_t2 <- rep(0.0, nt)
	pv_t2 <- rep(1.0, nt)
	for (i in 1:nt){
		idxt <- which(symbol_t1 %in% symbol_2[i])
		vt <- fc_t1[idxt]
		vt1 <- pv_t1[idxt]
		if (length(vt)>1){
			idxt1 <- which(abs(vt) == max(abs(vt)))
			fc_t2[i] <- vt[idxt1[1]]
			pv_t2[i] <- vt1[idxt1[1]]

		} else {
			fc_t2[i] <- vt
			pv_t2[i] <- vt1
		}
	}
	
	# GO-gene network 
	netGeneGO_up_1 <- getGeneGoNetworkV1a(gs_up_1, GO_up_2a[,1])
	# netGeneGO_dn_1 <- getGeneGoNetworkV1a(gs_dn_1, GO_dn_2a[,1])
	
	colnames(netGeneGO_up_1) <- c('Source', 'Target')
	# colnames(netGeneGO_dn_1) <- c('Source', 'Target')

	node_netGeneGO_up_1 <- generate_display_gene_go_network(GO_up_2a, netGeneGO_up_1, symbol_2, fc_t2)
	# node_netGeneGO_dn_1 <- generate_display_gene_go_network(GO_dn_2a, netGeneGO_dn_1, symbol_2, fc_t2)
	# node_netGeneGO_diff_1 <- generate_display_gene_go_network(GO_diff_2, netGeneGO_diff_1, symbol_1, fc_1)

	# GO-clustering (as there are too many GOs for CMAP drugs) to find Super-GOs 
	sim_GO_up_2 <- get_sim_GO(GO_up_2a)
	# sim_GO_dn_2 <- get_sim_GO(GO_dn_2a)
	dim(sim_GO_up_2)
	# dim(sim_GO_dn_2)

	# clusters <- hclust(as.dist(1-sim_GO_up_2))
	# plot(clusters)
	library(apcluster)
	c1 <- apclusterK(sim_GO_up_2, K=5)

	idx_c1 <- c1@idx
	name_t1 <- colnames(sim_GO_up_2)
	v_c1 <- unique(idx_c1)
	n_c1 <- length(v_c1)
	net_super_GO_up_2 <- c("super-GO", "GO")
	net_super_GO_gene_up_2 <- c("super-GO", "GeneSymbol") 
	for (i in 1:n_c1){
		idx_t1 <- which(idx_c1 %in% v_c1[i])
		n_t1 <- length(idx_t1)
		x_t1 <- cbind(rep(paste("super-GO", i, sep='-'), n_t1), name_t1[idx_t1])
		net_super_GO_up_2 <- rbind(net_super_GO_up_2, x_t1)
		#
		gs_t1 <- c()
		for (j in 1:n_t1){
			idx_t2 <- which(netGeneGO_up_1[,1] %in% name_t1[idx_t1[j]])
			gs_t1 <- c(gs_t1, netGeneGO_up_1[idx_t2,2])
		}
		gs_t1 <- unique(gs_t1)
		n_t2 <- length(gs_t1)
		x_t2 <- cbind(rep(paste("super-GO", i, sep='-'), n_t2), gs_t1)
		net_super_GO_gene_up_2 <- rbind(net_super_GO_gene_up_2, x_t2)
	}
	net_sGO_GO_gene_up_2 <- rbind(net_super_GO_up_2[-1,], net_super_GO_gene_up_2[-1,])
	
	colnames(net_sGO_GO_gene_up_2) <- c('Source', 'Target')

	#
	sGO_node <- unique(net_super_GO_up_2[-1,1])
	n_sGO <- length(sGO_node)

	node_t <- cbind(sGO_node, sGO_node, rep('superGO', n_sGO), rep('0.5', n_sGO))
	node_superGO_Gene <- rbind(node_netGeneGO_up_1, node_t)

	if (1 == 2){
		f_t1 <- paste(dir_res, paste(flag1, 'net_superGO_GO_gene_up_2-c5.txt', sep='-'), sep='/')
		write.table(net_sGO_GO_gene_up_2, f_t1, quote=F, sep='\t', row.names=F)

		f_t1 <- paste(dir_res, paste(flag1, 'Table_superGO_geneSymbol_up_2-c5.txt', sep='-'), sep='/')
		write.table(net_super_GO_gene_up_2, f_t1, quote=F, sep='\t', row.names=F)

		f_t1 <- paste(dir_res, paste(flag1, 'node_superGO_GO_gene_up_2-c5.txt', sep='-'), sep='/')
		write.table(node_superGO_Gene, f_t1, quote=F, sep='\t', row.names=F)
	}
}




if (1 == 2){ # bait-prey PPIs
	ppi_1 <- read.delim(f_ppi_1, sep='\t')
	names(ppi_1)
	gs_ppi <- as.character(ppi_1[,3])


	GO_t1 <- goEnrichAnalysisV1a(gs_ppi, union(symbol_A549, gs_ppi))
	dim(GO_t1)
	idx_t1 <- which(as.character(GO_t1[,3]) == "BP" & as.numeric(GO_t1[,4]) <= 1000 & as.numeric(GO_t1[,4]) >= 10 & GO_t1[,5] <= 0.05)
	GO_t2 <- GO_t1[idx_t1,]
	dim(GO_t2)

	if (1 == 2){
		flag1 <- "covid"
		f_t1 <- paste(dir_res, paste(flag1, 'list_GO_ppi_2.txt', sep='-'), sep='/')
		write.table(GO_t2, f_t1, quote=F, sep='\t', row.names=F)

		if (1 == 2){ # re-load the up-regulated GOs 
			GO_t2a <- read.delim(gsub('.txt', '_manual_selection1.txt', f_t1))
			x_t1 <- GO_t2a[,3]
			idx_t1 <- which(x_t1 > 0)
			GO_t2a <- GO_t2a[-idx_t1,-3]
			dim(GO_t2a)	
		}
	}

	netGeneGO_t1 <- getGeneGoNetworkV1a(gs_ppi, GO_t2[,1])
	colnames(netGeneGO_t1) <- c('Source', 'Target')
	length(unique(netGeneGO_t1[,2]))

	symbol_tx <- union(symbol_A549, gs_ppi)
	fc_tx <- rep(1.0, length(symbol_tx))
	node_netGeneGO_t1 <- generate_display_gene_go_network(GO_t2, netGeneGO_t1, symbol_tx, fc_tx)

	# node_netGeneGO_up_1 <- generate_display_gene_go_network(GO_up_2, netGeneGO_up_1, symbol_1, fc_1)

	# super-GOs, clustering analysis
	sim_GO_t1 <- get_sim_GO(GO_t2)
	dis_GO_t1 <- as.dist(1.0 - sim_GO_t1)

	library(cluster)
	library(apcluster)

	silhouette_score <- function(k){
		c1 <- apclusterK(sim_GO_t1, K=k)
		idx_c1 <- c1@idx
		ss <- silhouette(c1@idx, dis_GO_t1)
		mean(ss[, 3])
	}
	k <- 2:50
	avg_sil <- sapply(k, silhouette_score)
	plot(k, type='b', avg_sil, xlab='Number of clusters', ylab='Average Silhouette Scores', frame=FALSE)

	# clusters <- hclust(as.dist(1-sim_GO_up_2))
	# plot(clusters)
	library(apcluster)
	c1 <- apclusterK(sim_GO_t1, K=10)

	idx_c1 <- c1@idx
	name_t1 <- colnames(sim_GO_t1)
	v_c1 <- unique(idx_c1)
	n_c1 <- length(v_c1)
	net_super_GO_t1 <- c("super-GO", "GO")
	net_super_GO_gene_t1 <- c("super-GO", "GeneSymbol") 
	for (i in 1:n_c1){
		idx_t1 <- which(idx_c1 %in% v_c1[i])
		n_t1 <- length(idx_t1)
		x_t1 <- cbind(rep(paste("super-GO", i, sep='-'), n_t1), name_t1[idx_t1])
		net_super_GO_t1 <- rbind(net_super_GO_t1, x_t1)
		#
		gs_t1 <- c()
		for (j in 1:n_t1){
			idx_t2 <- which(netGeneGO_t1[,1] %in% name_t1[idx_t1[j]])
			gs_t1 <- c(gs_t1, netGeneGO_t1[idx_t2,2])
		}
		gs_t1 <- unique(gs_t1)
		n_t2 <- length(gs_t1)
		x_t2 <- cbind(rep(paste("super-GO", i, sep='-'), n_t2), gs_t1)
		net_super_GO_gene_t1 <- rbind(net_super_GO_gene_t1, x_t2)
	}
	net_sGO_GO_gene_t1 <- rbind(net_super_GO_t1[-1,], net_super_GO_gene_t1[-1,])
	colnames(net_sGO_GO_gene_t1) <- c('Source', 'Target')

	sGO_node <- unique(net_super_GO_t1[-1,1])
	n_sGO <- length(sGO_node)
	node_t <- cbind(sGO_node, sGO_node, rep('superGO', n_sGO), rep('0.5', n_sGO))
	node_superGO_Gene_t1 <- rbind(node_netGeneGO_t1, node_t)

	# net_superGO_GO_gene_t1 <- net_sGO_GO_gene_t1 
	# net_superGO_Gene_t1 <- net_super_GO_gene_t1
	# node_superGO_GO_gene_t1 <- node_superGO_Gene_t1
	#
	if (1 == 2){ # write out the results
		flag1 <- 'ppi' 
		f_t1 <- paste(dir_res, paste(flag1, 'net_superGO_GO_gene_dn_2.txt', sep='-'), sep='/')
		write.table(net_sGO_GO_gene_t1, f_t1, quote=F, sep='\t', row.names=F)
		f_t1 <- paste(dir_res, paste(flag1, 'genesets_super_GO_geneSymbol_dn_2.txt', sep='-'), sep='/')
		write.table(net_super_GO_gene_t1, f_t1, quote=F, sep='\t', row.names=F)
		f_t1 <- paste(dir_res, paste(flag1, 'node_super_GO_gene_dn_2.txt', sep='-'), sep='/')
		write.table(node_superGO_Gene_t1, f_t1, quote=F, sep='\t', row.names=F)
	}

}

if (1 == 2){ # kegg signaling 
	if (1 == 2){ # load in KeggPathways
	    library(graphite)
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
	    for (i in 1:length(KeggPathways)){
	      pt <- KeggPathways[[i]]
	      gs <- c()
	      xt <- pt[,1]; idxt <- regexpr('SYMBOL', xt); nt1 <- sum(idxt>0)
	      if (nt1>0){gs1 <- pt[,2]; gs1 <- gs1[which(idxt>0)]; gs <- union(gs, gs1)}
	      xt <- pt[,3]; idxt <- regexpr('SYMBOL', xt); nt2 <- sum(idxt>0)
	      if (nt2>0){gs2 <- pt[,4]; gs2 <- gs2[which(idxt>0)]; gs <- union(gs, gs1)}
	      if (max(nt1, nt2)>0){ KeggGenes <- union(KeggGenes, gs)}
	    }
	}
	#
	gs_kegg_ppi <- intersect(KeggGenes, gs_ppi);
	# infer which signaling pathways are involved?
	k <- 0
	net_gene_pathway <- c("Source", "Target")
	for (i in 1:n_p){
		gs <- c()
		pt <- KeggPathways[[i]]
		xt <- pt[,1]; idxt <- regexpr('SYMBOL', xt); nt1 <- sum(idxt>0)
	    if (nt1>0){gs1 <- pt[,2]; gs1 <- gs1[which(idxt>0)]; gs <- union(gs, gs1)}
	    xt <- pt[,3]; idxt <- regexpr('SYMBOL', xt); nt2 <- sum(idxt>0)
	    if (nt2>0){gs2 <- pt[,4]; gs2 <- gs2[which(idxt>0)]; gs <- union(gs, gs1)}

	    yt <- intersect(gs, gs_ppi)
	    n_yt <- length(yt)
	    if (n_yt >= 4){
	    	print(n_yt); print(names_pathways[i]); k <- k+1
	    	zt <- cbind(rep(names_pathways[i], n_yt), yt)
	    	net_gene_pathway <- rbind(net_gene_pathway, zt)
	    }
	} 

	if (1 == 2){ # write out the results 
		flag1 <- 'ppi' 
		f_t1 <- paste(dir_res, paste(flag1, 'net_gene_pathway1.txt', sep='-'), sep='/')
		write.table(net_gene_pathway, f_t1, quote=F, sep='\t', row.names=F)
	}

}

if (1 == 2){ # GO analysis
	# # A549
	# fc_1 <- fc_A549_cov2_vs_mock
	# pv_1 <- pv_A549_cov2_vs_mock
	# symbol_1 <- symbol_A549
	# flag1 <- 'A549'

	# NHBE
	fc_1 <- fc_NHBE_cov2_vs_mock
	pv_1 <- pv_NHBE_cov2_vs_mock
	symbol_1 <- symbol_NHBE
	flag1 <- 'NHBE'

	# take the maximum or minimum value of multiple probes (the same gene)
	symbol_2 <- unique(symbol_1)
	nt <- length(symbol_2)
	fc_t2 <- rep(0.0, nt)
	pv_t2 <- rep(1.0, nt)
	for (i in 1:nt){
		idxt <- which(symbol_1 %in% symbol_2[i])
		vt <- fc_1[idxt]
		vt1 <- pv_1[idxt]
		if (length(vt)>1){
			idxt1 <- which(abs(vt) == max(abs(vt)))
			fc_t2[i] <- vt[idxt1[1]]
			pv_t2[i] <- vt1[idxt1[1]]

		} else {
			fc_t2[i] <- vt
			pv_t2[i] <- vt1
		}
	}
	#

	gs_up_1 <- symbol_1[which(fc_1 >= log2(1.5)  & pv_1 <= 0.05)]
	gs_dn_1 <- symbol_1[which(fc_1 <= -log2(1.5) & pv_1 <= 0.05)]
	gs_diff_1 <- union(gs_up_1, gs_dn_1)
	length(gs_up_1)
	length(gs_dn_1)
	length(gs_diff_1)

	if (1 == 2){
		gs_t_x1a <- intersect(gs_diff_1, gs_ppi)
		gs_t_x1b <- intersect(symbol_1, gs_ppi) # 310 gene intersect
		fc_t_x1b <- fc_1[which(symbol_1 %in% gs_ppi)]
		mean(fc_t_x1b)
		length(gs_t_x1a)
		length(gs_up_1)
	}


	# idxt1 <- which(exp_symbol_0 %in% gs_up_1[3])
	# idxt2 <- which(symbol_1 %in% gs_up_1[3])
	# Exp[idxt1, c(7:12)]
	# t1 <- mean(as.numeric(Exp[idxt1, c(7:9)]))/mean(as.numeric(Exp[idxt1, c(10:12)]))
	# log2(1/t1)
	# fc_1[idxt2]

	GO_up_1 <- goEnrichAnalysisV1a(gs_up_1, symbol_1)
	GO_dn_1 <- goEnrichAnalysisV1a(gs_dn_1, symbol_1)
	GO_diff_1 <- goEnrichAnalysisV1a(gs_diff_1, symbol_1)

	dim(GO_up_1)
	dim(GO_dn_1)
	dim(GO_diff_1)
	#
	idx_up1 <- which(as.character(GO_up_1[,3]) == "BP" & as.numeric(GO_up_1[,4]) <= 1000 & as.numeric(GO_up_1[,4]) >= 10 & GO_up_1[,5] <= 0.05)
	idx_dn1 <- which(as.character(GO_dn_1[,3]) == "BP" & as.numeric(GO_dn_1[,4]) <= 1000 & as.numeric(GO_dn_1[,4]) >= 10 & GO_dn_1[,5] <= 0.05)
	idx_diff1 <- which(as.character(GO_diff_1[,3]) == "BP" & as.numeric(GO_diff_1[,4]) <= 1000 & as.numeric(GO_diff_1[,4]) >= 10 & GO_diff_1[,5] <= 0.05)

	GO_up_2 <- GO_up_1[idx_up1,]
	GO_dn_2 <- GO_dn_1[idx_dn1,]
	GO_diff_2 <- GO_diff_1[idx_diff1,]

	dim(GO_up_2)
	dim(GO_dn_2)
	dim(GO_diff_2)

	# remove the common GOs in both activated/inhibited GOs
	vt1 <- as.character(GO_up_2[,1])
	vt2 <- as.character(GO_dn_2[,1])
	xt1 <- intersect(vt1, vt2)
	idx_t1a <- which(GO_up_2[,1] %in% xt1)
	idx_t1b <- which(GO_dn_2[,1] %in% xt1)
	GO_up_2b <- GO_up_2[-idx_t1a,]
	GO_dn_2b <- GO_dn_2[-idx_t1b,]

	dim(GO_up_2b)
	dim(GO_dn_2b)

	if (1 == 2){ # manual selection of Up-regulated GOs
		f_t1 <- paste(dir_res, paste(flag1, 'list_GO_up_2.txt', sep='-'), sep='/')
		write.table(GO_up_2b, f_t1, quote=F, sep='\t', row.names=F)

		f_t2 <- paste(dir_res, paste(flag1, 'list_GO_dn_2.txt', sep='-'), sep='/')
		write.table(GO_dn_2b, f_t2, quote=F, sep='\t', row.names=F)

		# re-load the up-regulated GOs 
		GO_up_2a <- read.delim(gsub('.txt', '_manual_selection1.txt', f_t1))
		x_t1 <- GO_up_2a[,3]
		idx_t1 <- which(x_t1 > 0)
		GO_up_2a <- GO_up_2a[-idx_t1,-3]
		dim(GO_up_2a)

		# re-load the dn-regulated GOs 
		GO_dn_2a <- read.delim(gsub('.txt', '_manual_selection1.txt', f_t2))
		x_t1 <- GO_dn_2a[,3]
		idx_t1 <- which(x_t1 > 0)
		GO_dn_2a <- GO_dn_2a[-idx_t1,-3]
		dim(GO_dn_2a)
	}

	# GO-gene network 
	netGeneGO_up_1 <- getGeneGoNetworkV1a(gs_up_1, GO_up_2a[,1])
	netGeneGO_dn_1 <- getGeneGoNetworkV1a(gs_dn_1, GO_dn_2a[,1])
	# netGeneGO_diff_1 <- getGeneGoNetworkV1a(gs_diff_1, GO_diff_2[,1])
	
	colnames(netGeneGO_up_1) <- c('Source', 'Target')
	colnames(netGeneGO_dn_1) <- c('Source', 'Target')
	# colnames(netGeneGO_diff_1) <- c('Source', 'Target')
	#
	node_netGeneGO_up_1 <- generate_display_gene_go_network(GO_up_2a, netGeneGO_up_1, symbol_2, fc_t2)
	node_netGeneGO_dn_1 <- generate_display_gene_go_network(GO_dn_2a, netGeneGO_dn_1, symbol_2, fc_t2)
	# node_netGeneGO_diff_1 <- generate_display_gene_go_network(GO_diff_2, netGeneGO_diff_1, symbol_1, fc_1)

	# GO-clustering (as there are too many GOs for CMAP drugs) to find Super-GOs 
	sim_GO_up_2 <- get_sim_GO(GO_up_2a)
	sim_GO_dn_2 <- get_sim_GO(GO_dn_2a)
	dim(sim_GO_up_2)
	dim(sim_GO_dn_2)

	# clusters <- hclust(as.dist(1-sim_GO_up_2))
	# plot(clusters)
	library(apcluster)
	c1 <- apclusterK(sim_GO_up_2, K=10)

	idx_c1 <- c1@idx
	name_t1 <- colnames(sim_GO_up_2)
	v_c1 <- unique(idx_c1)
	n_c1 <- length(v_c1)
	net_super_GO_up_2 <- c("super-GO", "GO")
	net_super_GO_gene_up_2 <- c("super-GO", "GeneSymbol") 
	for (i in 1:n_c1){
		idx_t1 <- which(idx_c1 %in% v_c1[i])
		n_t1 <- length(idx_t1)
		x_t1 <- cbind(rep(paste("super-GO", i, sep='-'), n_t1), name_t1[idx_t1])
		net_super_GO_up_2 <- rbind(net_super_GO_up_2, x_t1)
		#
		gs_t1 <- c()
		for (j in 1:n_t1){
			idx_t2 <- which(netGeneGO_up_1[,1] %in% name_t1[idx_t1[j]])
			gs_t1 <- c(gs_t1, netGeneGO_up_1[idx_t2,2])
		}
		gs_t1 <- unique(gs_t1)
		n_t2 <- length(gs_t1)
		x_t2 <- cbind(rep(paste("super-GO", i, sep='-'), n_t2), gs_t1)
		net_super_GO_gene_up_2 <- rbind(net_super_GO_gene_up_2, x_t2)
	}
	net_sGO_GO_gene_up_2 <- rbind(net_super_GO_up_2[-1,], net_super_GO_gene_up_2[-1,])
	
	colnames(net_sGO_GO_gene_up_2) <- c('Source', 'Target')


	#
	sGO_node <- unique(net_super_GO_up_2[-1,1])
	n_sGO <- length(sGO_node)

	node_t <- cbind(sGO_node, sGO_node, rep('superGO', n_sGO), rep('0.5', n_sGO))
	node_superGO_Gene <- rbind(node_netGeneGO_up_1, node_t)

	if (1 == 2){



		f_t1 <- paste(dir_res, paste(flag1, 'net_super_GO_GO_gene_up_2.txt', sep='-'), sep='/')
		write.table(net_sGO_GO_gene_up_2, f_t1, quote=F, sep='\t', row.names=F)

		f_t1 <- paste(dir_res, paste(flag1, 'net_super_GO_geneSymbol_up_2.txt', sep='-'), sep='/')
		write.table(net_super_GO_gene_up_2, f_t1, quote=F, sep='\t', row.names=F)

		f_t1 <- paste(dir_res, paste(flag1, 'node_super_GO_gene_up_2.txt', sep='-'), sep='/')
		write.table(node_superGO_Gene, f_t1, quote=F, sep='\t', row.names=F)
	}

	# write out results
	if (1 == 2){
		f_t1 <- paste(dir_res, paste(flag1, 'gene_up_1.txt', sep='-'), sep='/')
		write.table(gs_up_1, f_t1, quote=F, sep='\t', row.names=F)
		
		f_t2 <- paste(dir_res, paste(flag1, 'gene_down_1.txt', sep='-'), sep='/')
		write.table(gs_dn_1, f_t2, quote=F, sep='\t', row.names=F)

		f_t3 <- paste(dir_res, paste(flag1, 'Go_up_1.txt', sep='-'), sep='/')
		write.table(GO_up_2, f_t3, quote=F, sep='\t', row.names=F)
		
		f_t4 <- paste(dir_res, paste(flag1,'Go_down_1.txt', sep='-'), sep='/')
		write.table(GO_dn_2, f_t4, quote=F, sep='\t', row.names=F)

		f_t5 <- paste(dir_res, paste(flag1, 'Go_diff_1.txt', sep='-'), sep='/')
		write.table(GO_diff_2, f_t5, quote=F, sep='\t', row.names=F)

		f_t6 <- paste(dir_res, paste(flag1, 'Go_Gene_Up_Net_1.txt', sep='-'), sep='/')
		write.table(netGeneGO_up_1, f_t6, quote=F, sep='\t', row.names=F)

		f_t7 <- paste(dir_res, paste(flag1, 'Go_Gene_Down_Net_1.txt', sep='-'), sep='/')
		write.table(netGeneGO_dn_1, f_t7, quote=F, sep='\t', row.names=F)

		f_t8 <- paste(dir_res, paste(flag1, 'Go_Gene_Diff_Net_1.txt', sep='-'), sep='/')
		write.table(netGeneGO_diff_1, f_t8, quote=F, sep='\t', row.names=F)

		f_t9 <- paste(dir_res, paste(flag1, 'Go_Gene_Up_Node_1.txt', sep='-'), sep='/')
		write.table(node_netGeneGO_up_1, f_t9, quote=F, sep='\t', row.names=F)

		f_t10 <- paste(dir_res, paste(flag1, 'Go_Gene_Down_Node_1.txt', sep='-'), sep='/')
		write.table(node_netGeneGO_dn_1, f_t10, quote=F, sep='\t', row.names=F)

		f_t11 <- paste(dir_res, paste(flag1, 'Go_Gene_Diff_Node_1.txt', sep='-'), sep='/')
		write.table(node_netGeneGO_diff_1, f_t11, quote=F, sep='\t', row.names=F)	


	}

} # end of GO analysis
#

#
if (1 == 2){ # get drugs from CMAP (2837 drugs/compounds) based on the super-GOs. 
	# FDA drugs
	f_t1 <- paste(dir_data0, 'drug links fda.csv', sep='/')
	d_fda <- read.delim(f_t1, sep=',')
	d_fda <- as.character(d_fda[,2]) # 3556

	f_t1 <- paste(dir_res, 'cmap_res1/new1_s5-1a.gct', sep='/')
	d_t1 <- read.delim(f_t1, sep='\t');
	f_t1 <- paste(dir_res, 'cmap_res1/new1_s5-2a.gct', sep='/')
	d_t2 <- read.delim(f_t1, sep='\t')
	f_t1 <- paste(dir_res, 'cmap_res1/new1_s5-3a.gct', sep='/')
	d_t3 <- read.delim(f_t1, sep='\t')
	f_t1 <- paste(dir_res, 'cmap_res1/new1_s5-4a.gct', sep='/')
	d_t4 <- read.delim(f_t1, sep='\t')
	f_t1 <- paste(dir_res, 'cmap_res1/new1_s5-5a.gct', sep='/')
	d_t5 <- read.delim(f_t1, sep='\t')


	# f_t1 <- paste(dir_res, 'cmap_res1/new1-s5-1.gct', sep='/')
	# d_t1 <- read.delim(f_t1, sep='\t');
	# f_t1 <- paste(dir_res, 'cmap_res1/new1-s5-2.gct', sep='/')
	# d_t2 <- read.delim(f_t1, sep='\t')
	# f_t1 <- paste(dir_res, 'cmap_res1/new1-s5-3.gct', sep='/')
	# d_t3 <- read.delim(f_t1, sep='\t')
	# f_t1 <- paste(dir_res, 'cmap_res1/new1-s5-4.gct', sep='/')
	# d_t4 <- read.delim(f_t1, sep='\t')
	# f_t1 <- paste(dir_res, 'cmap_res1/new1-s5-5.gct', sep='/')
	# d_t5 <- read.delim(f_t1, sep='\t')

	# f_t1 <- paste(dir_res, 'cmap_res1/covid-s6.gct', sep='/')
	# d_t6 <- read.delim(f_t1, sep='\t')
	# f_t1 <- paste(dir_res, 'cmap_res1/covid-s7.gct', sep='/')
	# d_t7 <- read.delim(f_t1, sep='\t')
	# f_t1 <- paste(dir_res, 'cmap_res1/covid-s9.gct', sep='/')
	# d_t9 <- read.delim(f_t1, sep='\t')
	# f_t1 <- paste(dir_res, 'cmap_res1/covid-s10.gct', sep='/')
	# d_t10 <- read.delim(f_t1, sep='\t')

	#
	st1 <- 'metformin'
	# st1 <- "naltrexone"

	d_t0 <- d_t1
	d_t0 <- d_t0[-c(1:5),c(3,10)]
	x0 <- toupper(d_fda)
	x1 <- toupper(d_t0[,1])
	x2 <- intersect(x0, x1)
	idx_t1 <- which(x1 %in% x2)
	d_t0a1 <- d_t0[idx_t1, ]
	which(d_t0a1[,1] %in% st1)

	d_t0 <- d_t2; 	d_t0 <- d_t0[-c(1:5),c(3,10)]
	x1 <- toupper(d_t0[,1]); 	idx_t1 <- which(x1 %in% x2)
	d_t0a2 <- d_t0[idx_t1, ]; 	which(d_t0a2[,1] %in% st1)

	d_t0 <- d_t3; 	d_t0 <- d_t0[-c(1:5),c(3,10)]
	x1 <- toupper(d_t0[,1]); 	idx_t1 <- which(x1 %in% x2)
	d_t0a3 <- d_t0[idx_t1, ]; 	which(d_t0a3[,1] %in% st1)

	d_t0 <- d_t4; 	d_t0 <- d_t0[-c(1:5),c(3,10)]
	x1 <- toupper(d_t0[,1]); 	idx_t1 <- which(x1 %in% x2)
	d_t0a4 <- d_t0[idx_t1, ]; 	which(d_t0a4[,1] %in% st1)

	d_t0 <- d_t5; 	d_t0 <- d_t0[-c(1:5),c(3,10)]
	x1 <- toupper(d_t0[,1]); 	idx_t1 <- which(x1 %in% x2)
	d_t0a5 <- d_t0[idx_t1, ]; 	which(d_t0a5[,1] %in% st1)

	# d_t0 <- d_t6; 	d_t0 <- d_t0[-c(1:5),c(3,10)]
	# x1 <- toupper(d_t0[,1]); 	idx_t1 <- which(x1 %in% x2)
	# d_t0a6 <- d_t0[idx_t1, ]; 	which(d_t0a6[,1] %in% st1)

	# d_t0 <- d_t7; 	d_t0 <- d_t0[-c(1:5),c(3,10)]
	# x1 <- toupper(d_t0[,1]); 	idx_t1 <- which(x1 %in% x2)
	# d_t0a7 <- d_t0[idx_t1, ]; 	which(d_t0a7[,1] %in% st1)

	# d_t0 <- d_t9; 	d_t0 <- d_t0[-c(1:5),c(3,10)]
	# x1 <- toupper(d_t0[,1]); 	idx_t1 <- which(x1 %in% x2)
	# d_t0a9 <- d_t0[idx_t1, ]; 	which(d_t0a9[,1] %in% st1)

	# d_t0 <- d_t10; 	d_t0 <- d_t0[-c(1:5),c(3,10)]
	# x1 <- toupper(d_t0[,1]); 	idx_t1 <- which(x1 %in% x2)
	# d_t0a10 <- d_t0[idx_t1, ]; 	which(d_t0a10[,1] %in% st1)

	#
	# d_m1 <- cbind(as.character(d_t0a1[1:100,1]),as.character(d_t0a2[1:100,1]),as.character(d_t0a3[1:100,1]),as.character(d_t0a4[1:100,1]),as.character(d_t0a5[1:100,1]),
	# 		as.character(d_t0a6[1:100,1]), as.character(d_t0a7[1:100,1]), as.character(d_t0a9[1:100,1]), as.character(d_t0a10[1:100,1]))
	
	N1 <- 100
	d_m1 <- cbind(as.character(d_t0a1[1:N1,1]),as.character(d_t0a2[1:N1,1]),as.character(d_t0a3[1:N1,1]),
		as.character(d_t0a4[1:N1,1]),as.character(d_t0a5[1:N1,1]))

	d_v1 <- as.vector(d_m1)
	d_l1 <- unique(d_v1)
	drug_f1 <- table(d_v1)
	#
	drug_f2 <- drug_f1[order(drug_f1, decreasing=T)]
	drug_go_f1 <- unique(d_v1)

	t1 <- drug_f2[drug_f2>=3] 
	t2 <- names(t1)
	drug_go_f5 <- t2


	#
	if (1 == 2){ # get the other drugs in the top 
		d_m2 <- matrix('NA', dim(d_m1)[1], dim(d_m1)[2])
		i <- 1; xt <- setdiff(d_m1[,i], t2); xt <- unique(xt); nt<-length(xt); d_m2[1:nt,i] <-xt 
		i <- 2; xt <- setdiff(d_m1[,i], t2); xt <- unique(xt); nt<-length(xt); d_m2[1:nt,i] <-xt 
		i <- 3; xt <- setdiff(d_m1[,i], t2); xt <- unique(xt); nt<-length(xt); d_m2[1:nt,i] <-xt 
		i <- 4; xt <- setdiff(d_m1[,i], t2); xt <- unique(xt); nt<-length(xt); d_m2[1:nt,i] <-xt 
		i <- 5; xt <- setdiff(d_m1[,i], t2); xt <- unique(xt); nt<-length(xt); d_m2[1:nt,i] <-xt 
		# i <- 6; xt <- setdiff(d_m1[,i], t2); xt <- unique(xt); nt<-length(xt); d_m2[1:nt,i] <-xt 
		# i <- 7; xt <- setdiff(d_m1[,i], t2); xt <- unique(xt); nt<-length(xt); d_m2[1:nt,i] <-xt 
		# i <- 8; xt <- setdiff(d_m1[,i], t2); xt <- unique(xt); nt<-length(xt); d_m2[1:nt,i] <-xt 
		# i <- 9; xt <- setdiff(d_m1[,i], t2); xt <- unique(xt); nt<-length(xt); d_m2[1:nt,i] <-xt 

		#
		f_t1 <- paste(dir_res, 'Drugs_4_super-GO_1.txt', sep='/')
		write.table(d_m2, f_t1, quote=F, sep='\t', row.names=F)

	}

	#
	if (1 == 2){ # get the drug information
		y_t1 <- as.character(d_t1[,3])
		idx_t1 <- which(y_t1 %in% t2)
		y_t2 <- d_t1[idx_t1,c(3:6)]
		y_t3 <- unique(y_t2)
		# get the unique data
		z_t1 <- unique(y_t3[,1])
		n_t1 <- length(z_t1)
		idx_T1 <- rep(0,n_t1)
		for (i in 1:n_t1){
			idx_t1a <- which(y_t3[,1] %in% z_t1[i])
			idx_T1[i] <- idx_t1a[1]
		}
		y_t3a <- y_t3[idx_T1,]

		#
		f_t1 <- paste(dir_res, 'Drug_high_frequency_info_1.txt', sep='/')
		write.table(y_t3a, f_t1, quote=F, sep='\t', row.names=F)

	}

	if (1 == 2){ # comparision with reported drugs
		f_t1 <- paste(dir_res, 'drug_covid_reported_1.txt', sep='/')
		x_t1 <- read.delim(f_t1, header=F, sep='\t')
		x_t1 <- toupper(as.character(x_t1[[1]]))
		t_2a <- toupper(t2)
		x_t2 <- intersect(x_t1, t_2a)
		#
		f_t1 <- paste(dir_res, 'drug_covirus_reported_1.txt', sep='/')
		x_t1 <- read.delim(f_t1, header=F, sep='\t')
		x_t1 <- toupper(as.character(x_t1[[1]]))
		t_2a <- toupper(t2)
		#
		x_t2 <- intersect(x_t1, t_2a)


	}

	if (1 == 2){ # generate the superGO-GO-Drug Networks
		net_sGO_drug_1 <- c('source', 'target')
		idx_t1 <- which(t2 %in% d_m1[,1]); nt<- length(idx_t1); net_sGO_drug_1 <- rbind(net_sGO_drug_1, cbind(rep('super-GO-1',nt), t2[idx_t1]))
		idx_t2 <- which(t2 %in% d_m1[,2]); nt<- length(idx_t2); net_sGO_drug_1 <- rbind(net_sGO_drug_1, cbind(rep('super-GO-2',nt), t2[idx_t2]))
		idx_t3 <- which(t2 %in% d_m1[,3]); nt<- length(idx_t3); net_sGO_drug_1 <- rbind(net_sGO_drug_1, cbind(rep('super-GO-3',nt), t2[idx_t3]))
		idx_t4 <- which(t2 %in% d_m1[,4]); nt<- length(idx_t4); net_sGO_drug_1 <- rbind(net_sGO_drug_1, cbind(rep('super-GO-4',nt), t2[idx_t4]))
		idx_t5 <- which(t2 %in% d_m1[,5]); nt<- length(idx_t5); net_sGO_drug_1 <- rbind(net_sGO_drug_1, cbind(rep('super-GO-5',nt), t2[idx_t5]))
		# idx_t6 <- which(t2 %in% d_m1[,6]); nt<- length(idx_t6); net_sGO_drug_1 <- rbind(net_sGO_drug_1, cbind(rep('super-GO-6',nt), t2[idx_t6]))
		# idx_t7 <- which(t2 %in% d_m1[,7]); nt<- length(idx_t7); net_sGO_drug_1 <- rbind(net_sGO_drug_1, cbind(rep('super-GO-7',nt), t2[idx_t7]))
		# idx_t8 <- which(t2 %in% d_m1[,8]); nt<- length(idx_t8); net_sGO_drug_1 <- rbind(net_sGO_drug_1, cbind(rep('super-GO-9',nt), t2[idx_t8]))
		# idx_t9 <- which(t2 %in% d_m1[,9]); nt<- length(idx_t9); net_sGO_drug_1 <- rbind(net_sGO_drug_1, cbind(rep('super-GO-10',nt), t2[idx_t9]))

		node_sGO_drug_1 <- c('Name1', 'Type1')
		node1 <- unique(net_sGO_drug_1[-1,1]); nt1 <- length(node1)
		node2 <- unique(net_sGO_drug_1[-1,2]); nt2 <- length(node2)
		node_sGO_drug_1 <- rbind(node_sGO_drug_1, cbind(node1, rep('super-GO', nt1)))
		node_sGO_drug_1 <- rbind(node_sGO_drug_1, cbind(node2, rep('Drug', nt2)))


	}

	if (1 == 2){ # generate the individual superGO-GO-Drug Networks
		net_sGO_drug_2 <- c('source', 'target')
		nt <- 100
		net_sGO_drug_2 <- rbind(net_sGO_drug_2, cbind(rep('super-GO-1',nt), d_m1[1:nt,1]))
		net_sGO_drug_2 <- rbind(net_sGO_drug_2, cbind(rep('super-GO-2',nt), d_m1[1:nt,2]))
		net_sGO_drug_2 <- rbind(net_sGO_drug_2, cbind(rep('super-GO-3',nt), d_m1[1:nt,3]))
		net_sGO_drug_2 <- rbind(net_sGO_drug_2, cbind(rep('super-GO-4',nt), d_m1[1:nt,4]))
		net_sGO_drug_2 <- rbind(net_sGO_drug_2, cbind(rep('super-GO-5',nt), d_m1[1:nt,5]))
		# net_sGO_drug_2 <- rbind(net_sGO_drug_2, cbind(rep('super-GO-6',nt), d_m1[1:nt,6]))
		# net_sGO_drug_2 <- rbind(net_sGO_drug_2, cbind(rep('super-GO-7',nt), d_m1[1:nt,7]))
		# net_sGO_drug_2 <- rbind(net_sGO_drug_2, cbind(rep('super-GO-9',nt), d_m1[1:nt,8]))
		# net_sGO_drug_2 <- rbind(net_sGO_drug_2, cbind(rep('super-GO-10',nt),d_m1[1:nt,9]))

		node_sGO_drug_2 <- c('Name1', 'Type1')
		node1 <- unique(net_sGO_drug_2[-1,1]); nt1 <- length(node1)
		node2 <- unique(net_sGO_drug_2[-1,2]); nt2 <- length(node2)
		node_sGO_drug_2 <- rbind(node_sGO_drug_2, cbind(node1, rep('super-GO', nt1)))
		node_sGO_drug_2 <- rbind(node_sGO_drug_2, cbind(node2, rep('Drug', nt2)))


	}

	# if (1 == 2){ # drugbank target on the kegg signaling network
	# 	# get the kegg signaling network
	# 	f_t1 <- paste(dir_res, 'kegg_signaling_upsigNet-Activation.txt', sep='/')
	# 	net_xt1 <- read.delim(f_t1, sep='\t')
	# 	node_xt1 <-union(net_xt1[,1], net_xt1[,2])

	# 	f_t1 <- paste(dir_data0, 'DrugBank-Drug-Target.txt', sep='/')
	# 	d_xt1 <- read.delim(f_t1, sep='\t')
	# 	idx_xt1 <- which(d_xt1[,2] %in% node_xt1)
	# 	d_xt2 <- d_xt1[idx_xt1,]
	# 	d_xt3 <- unique(as.character(d_xt2[,1]))
	# 	d_xt4 <- intersect(d_xt3, toupper(d_fda))
	# 	idx_xt2 <- which(d_xt2[,1] %in% d_xt4)
	# 	d_xt2a <- d_xt2[idx_xt2,]
	# 	d_xt2a[,1] <- tolower(d_xt2a[,1])

	# 	table_xt1 <- table(as.character(d_xt2a[,2]))
	# 	if (1 == 2){ #

	# 		t1 <- sort(table_xt1)
	# 		tar_t1 <- names(t1)
	# 		n_t1 <- length(tar_t1)
	# 		x_t1 <- matrix('', n_t1, 100)
	# 		z_t1 <- as.character(d_xt2a[,1])
	# 		z_t2 <- as.character(d_xt2a[,2])
	# 		for (i in 1:n_t1){
	# 			st <- tar_t1[i]
	# 			y_t1 <- which(z_t2 %in% st)
	# 			n_yt1 <- length(y_t1)
	# 			x_t1[i, 1:n_yt1] <- z_t1[y_t1]
	# 		}
	# 		x_t2 <- cbind(tar_t1, x_t1)

	# 		#
	# 		f_t1 <- paste(dir_res, 'target_drug_Tabel_1.txt', sep='/')
	# 		write.table(x_t2, f_t1, quote=F, sep='\t', row.names=F)

	# 		f_t1 <- paste(dir_res, 'target_drug_Tabel_2.txt', sep='/')
	# 		write.table(x_t1, f_t1, quote=F, sep=',', row.names=F)

	# 	}

	# 	# generate drug-target signaling network
	# 	net_xt2 <- rbind(as.matrix(net_xt1[,c(1,2)]), as.matrix(d_xt2a))

	# 	f_t1 <- paste(dir_res, 'Drug_signaling_net-up-Activation.txt', sep='/')
	# 	write.table(net_xt2, f_t1, quote=F, sep='\t', row.names=F)

	# 	# node attributes
	# 	z_t0 <- unique(d_xt2a[,2])
	# 	z_t1 <- tolower(unique(d_xt2a[,1]))
	# 	node_yt1 <- c(node_xt1, z_t1)
	# 	type_yt1 <- c(rep('Gene', length(node_xt1)), rep('Drug', length(z_t1)))
	# 	isTarget <- rep(0, length(type_yt1))
	# 	isTarget[which(node_yt1 %in% z_t0)] <- 1
	# 	node_att_1a <- cbind(node_yt1, type_yt1, isTarget)

	# 	f_t1 <- paste(dir_res, 'Drug_signaling_Node-up-Activation.txt', sep='/')
	# 	write.table(node_att_1a, f_t1, quote=F, sep='\t', row.names=F)
	# }

	if (1 == 2){
		f_t1 <- paste(dir_res, paste(flag1, 'net_sGO_drug_2.txt', sep='-'), sep='/')
		write.table(net_sGO_drug_2, f_t1, quote=F, sep='\t', col.names=F, row.names=F)	

		f_t1 <- paste(dir_res, paste(flag1, 'node_sGO_drug_2.txt', sep='-'), sep='/')
		write.table(node_sGO_drug_2, f_t1, quote=F, sep='\t', col.names=F, row.names=F)		


		f_t1 <- paste(dir_res, paste(flag1, 'net_sGO_drug_1.txt', sep='-'), sep='/')
		write.table(net_sGO_drug_1, f_t1, quote=F, sep='\t', col.names=F, row.names=F)	

		f_t1 <- paste(dir_res, paste(flag1, 'node_sGO_drug_1.txt', sep='-'), sep='/')
		write.table(node_sGO_drug_1, f_t1, quote=F, sep='\t', col.names=F, row.names=F)		

		f_t1 <- paste(dir_res, paste(flag1, 'topDrugsFrequency.txt', sep='-'), sep='/')
		write.table(t1, f_t1, quote=F, sep='\t', row.names=F)

		# f_t1 <- paste(dir_res, paste(flag1, 'net_super_GO_geneSymbol_up_2.txt', sep='-'), sep='/')
		# write.table(net_super_GO_gene_up_2, f_t1, quote=F, sep='\t', row.names=F)

		# f_t1 <- paste(dir_res, paste(flag1, 'node_super_GO_gene_up_2.txt', sep='-'), sep='/')
		# write.table(node_superGO_Gene, f_t1, quote=F, sep='\t', row.names=F)
	}

}





if (1 == 2){ # Compare with clinical trials
	f_t1 <- paste(dir_res, 'table_trials - 2020-05-27 05_38_12.csv', sep='/')
	# f_t1 <- paste(dir_res, 'table_trials - 2020-06-29 21_40_47.csv', sep='/')
	drug_ct <- read.delim(f_t1, sep=',') 
	#
	xt1 <- as.character(drug_ct[,1]) #665 unique clinical trials
	xt2 <- as.character(drug_ct[,16]) #276 unique treatments/combination treatments
	xt2u <- unique(xt2)
	xt2u1 <- strsplit(xt2u, ',', fixed=T)
	nt <- length(xt2u1)
	drug_ct1 <- c()
	for (i in 1:nt){
		st <- xt2u1[[i]]
		drug_ct1 <- c(drug_ct1, st)
	}
	drug_ct1 <- unique(drug_ct1)
	drug_ct2 <- gsub(' ', '', drug_ct1, fixed=T)
	drug_ct2 <- unique(drug_ct2) #252 drugs/treatments
	drug_ct2 <- gsub("/-", '/', drug_ct2, fixed=T)

	drug_ct2a <- strsplit(drug_ct2, '/', fixed=T)
	# drug_ct2a <- strsplit(drug_ct2, '+', fixed=T) # for drug combinations
	nt <- length(drug_ct2a)
	drug_ct3 <- c()
	for (i in 1:nt){
		st <- drug_ct2a[[i]]
		drug_ct3 <- c(drug_ct3, st)
	}
	drug_ct3 <- unique(drug_ct3)
	drug_ct4 <- toupper(drug_ct3)
	drug_ct4 <- unique(drug_ct4) #255

	drug_ct5 <- gsub("HCQ", 'Hydroxychloroquine', drug_ct4, fixed=T)
	drug_ct5 <- gsub("AZT", "Zidovudine", drug_ct5, fixed=T)
	drug_ct5 <- gsub("TCZ", "Tocilizumab", drug_ct5, fixed=T)
	drug_ct5 <- gsub("LPV", "Lopinavir", drug_ct5, fixed=T)
	drug_ct5 <- gsub("CHLOROQUINEANALOG", "CHLOROQUINE", drug_ct5, fixed=T)
	drug_ct5 <- c(drug_ct5, "CHLOROQUINE")

	drug_ct5 <- toupper(drug_ct5)
	#
	if (1 == 2){ # directly copy the drug/treatment names from the dashboard
		f_t1 <- paste(dir_res, 'table_trials - 2020-06-29 21_40_47_manual.csv', sep='/')
		drug_ct <- read.delim(f_t1, sep=',') 
		xt2 <- as.character(drug_ct[,1])
		xt2u <- unique(xt2)
		xt2u1 <- strsplit(xt2u, ',', fixed=T)
		nt <- length(xt2u1)
		drug_ct1 <- c()
		for (i in 1:nt){
			st <- xt2u1[[i]]
			drug_ct1 <- c(drug_ct1, st)
		}
		drug_ct1 <- unique(drug_ct1)
		drug_ct2 <- gsub(' ', '', drug_ct1, fixed=T)
		drug_ct2 <- unique(drug_ct2) #252 drugs/treatments
		drug_ct2 <- gsub("/-", '/', drug_ct2, fixed=T)
		drug_ct2 <- gsub("or", '/', drug_ct2, fixed=T)
		drug_ct2 <- gsub("with/without", '/', drug_ct2, fixed=T)
		drug_ct2a <- strsplit(drug_ct2, '/', fixed=T)
		# drug_ct2a <- strsplit(drug_ct2, '+', fixed=T) # for drug combinations
		nt <- length(drug_ct2a)
		drug_ct3 <- c()
		for (i in 1:nt){
			st <- drug_ct2a[[i]]
			drug_ct3 <- c(drug_ct3, st)
		}
		drug_ct3 <- unique(drug_ct3)
		drug_ct4 <- toupper(drug_ct3)
		drug_ct4 <- unique(drug_ct4) #255
		drug_ct5b <- toupper(drug_ct4)
		#
	}
	drug_ct5 <- unique(c(drug_ct5, drug_ct5b)) 

	drug_fda1 <- gsub(' ', '', d_fda, fixed=T)
	drug_fda2 <- toupper(drug_fda1)


	drug_common_t1 <- intersect(drug_ct5, drug_fda2) # 76 fda drugs
	length(drug_common_t1)
	#
	drug_signaling_1 <- toupper(unique(d_xt2a[,1]))
	drug_signaling_2 <- gsub(' ', '', drug_signaling_1, fixed=T)
	drug_common_t2 <- intersect(drug_common_t1,drug_signaling_2)
	length(drug_common_t2) # 30 drugs

	drug_go_f5a <- gsub(' ', '', drug_go_f5, fixed=T)
	drug_go_f5a <- toupper(drug_go_f5a)
	drug_common_t3 <- intersect(drug_common_t1,drug_go_f5a)
	length(drug_common_t3) # 1 drugs

	drug_go_f1a <- gsub(' ', '', drug_go_f1, fixed=T)
	drug_go_f1a <- toupper(drug_go_f1a)
	drug_common_t4 <- intersect(drug_common_t1, drug_go_f1a)
	length(drug_common_t4) # 10 drugs

	drug_common_t5 <- intersect(drug_common_t2, drug_common_t4) # "DEXAMETHASONE" "NAPROXEN

	drug_common_t6 <- intersect(drug_go_f1a, drug_signaling_2)
	print(drug_common_t6)

	length(drug_signaling_2)

	drug_common_all <- intersect(drug_common_t2, drug_common_t4)
	tolower(intersect(drug_common_t2, drug_common_t4))

	#
	if (1 == 2){ # write out the results
		f_t1 <- paste(dir_res, 'fda_clinicaltrials_drugs_1a.txt', sep='/')
		xt1 <- c(tolower(drug_common_t1), rep('NA', 1))
		dim(xt1) <- c(23, 5)
		write.table(xt1, f_t1, quote=F, sep='\t', row.names=F)

		f_t2 <- paste(dir_res, 'fda_clinicaltrials_signaling_drugs_1.txt', sep='/')
		xt1 <- tolower(drug_common_t2); dim(xt1) <- c(1, length(xt1))
		write.table(xt1, f_t2, quote=F, sep=',', row.names=F)

	}


} 
#


