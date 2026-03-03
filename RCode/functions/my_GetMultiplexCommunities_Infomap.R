my_GetMultiplexCommunities_Infomap <-
    function(g.list,
             bin.path = NA,
             isDirected,
             seed = 12345,
             includeSelfLinks = F,
             numTrials = 100,
             twoLevel = T,
             preclusterMultiplex = F,
             addMissingPhysicalNodes = T,
             hardPartitions = F,
             verbose = T,
             addAggregateAnalysis = T,
             multilayerRelaxRate = NA,
             multilayerJSRelaxRate = NA,
             outputPrefix = "multimap") {
        if (is.na(bin.path) || !file.exists(bin.path)) {
            stop(
                "Error! You must provide a valid path to the INFOMAP bin. Likely you will find it in the bin/ folder of muxviz, or you must compile it from source in src/ folder. If this is the case, just unzip the infomap archive and run  make  that will generate executable Infomap. Feel free to move the file where you prefer and provide the full path as an argument to this function."
            )
        }
        
        Layers <- length(g.list)
        
        tmpname <- outputPrefix
        inputFile <- paste0(tmpname, "_infomap.edges")
        if (file.exists(inputFile))
            file.remove(inputFile)
        fileConn <- file(inputFile, open = "at")
        
        cat('1/2 Setting up the algorithms...\n')
        
        mergedEdgelist <- data.frame()
        
        #write in the Infomap multilayer format for edge-colored networks
        #this part is different from the one in GetMultilayerCommunities_Infomap
        writeLines(c("*Intra", "#level node node weight"), fileConn)
        
        for (l in 1:Layers) {
            edges <- igraph::get.edgelist(g.list[[l]])
            weights <- igraph::E(g.list[[l]])$weight
            if (is.null(weights))
                weights <- rep(1, igraph::gsize(g.list[[l]]))
            mergedEdgelist <- rbind(mergedEdgelist,
                                    data.frame(
                                        layer = l,
                                        from = edges[, 1],
                                        to = edges[, 2],
                                        weight = weights
                                    ))
            
            if (!isDirected) {
                #this is because multimap requires both directions specified, even for undirected networks
                mergedEdgelist <- rbind(mergedEdgelist,
                                        data.frame(
                                            layer = l,
                                            from = edges[, 2],
                                            to = edges[, 1],
                                            weight = weights
                                        ))
            }
        }
        utils::write.table(
            mergedEdgelist,
            file = fileConn,
            row.names = F,
            col.names = F,
            quote = F
        )
        close(fileConn)
        
        cat('2/2 Finding communities...\n')
        cat(' + Multiplex network...\n')
        
        exePath <- bin.path
        
        outname <- tmpname
        outdir <- getwd()
        
        #default flags
        exeFlags <- paste(inputFile, outdir)
  #      exeFlags <- paste(exeFlags, "--input-format multilayer")
        exeFlags <- paste(exeFlags, "--clu")
        
        exeFlags <- paste(exeFlags, "--seed", seed)
        exeFlags <- paste(exeFlags, "--num-trials", numTrials)
        
        if (isDirected) {
            exeFlags <- paste(exeFlags, "-f directed")
        } else {
            exeFlags <- paste(exeFlags, "-f undirected")
        }
        
        if (is.na(multilayerRelaxRate) && is.na(multilayerJSRelaxRate)) {
            stop("ERROR! You must specify a non-negative value for the relax rate or the JS relax
           rate.")
        } else if (is.na(multilayerRelaxRate)) {
            # then multilayerJSRelaxRate must be non-negative
            if (multilayerJSRelaxRate >= 0) {
                exeFlags <- paste(
                    exeFlags,
                    "--multilayer-js-relax-rate",
                    multilayerJSRelaxRate
                )
            }
        } else if (is.na(multilayerJSRelaxRate)) {
            # then multilayerRelaxRate must be non-negative
            if (multilayerRelaxRate >= 0) {
                exeFlags <- paste(
                    exeFlags,
                    "--multilayer-relax-rate",
                    multilayerRelaxRate
                )
            }
        } else if (!is.na(multilayerRelaxRate) && !is.na(multilayerJSRelaxRate)) {
            if (multilayerRelaxRate >= 0 && multilayerJSRelaxRate >= 0) {
                stop("ERROR! You must specify either the relax rate or the JS relax rate, not both.")
            }
        }
        
        if (includeSelfLinks) {
            exeFlags <- paste(exeFlags, "--include-self-links")
        }
        
        if (twoLevel) {
            exeFlags <- paste(exeFlags, "--two-level")
        }
        
        if (preclusterMultiplex) {
            exeFlags <- paste(exeFlags, "--pre-cluster-multiplex")
        }
        
        if (addMissingPhysicalNodes) {
            exeFlags <- paste(exeFlags, "--multilayer-add-missing-nodes")
        }
        
        if (hardPartitions) {
            exeFlags <- paste(exeFlags, "--hard-partitions")
        }
        
        if (verbose) {
            exeFlags <- paste(exeFlags, "-vvv")
        }
        
        exeFlags <- paste(exeFlags, "--out-name", outname)
        
        # call infomap
        system(paste(exePath, exeFlags), intern = T)
        
        
        #read output. Here I could redirect the output inside the R environment.. but
        #for compatibility with the rest of the code I prefer to read a file
        communityList <- list()
        
        #import the results (clu and modularity value)
        resultFile <- paste0(outputPrefix, "_expanded.clu")
        wmemb_membership <- utils::read.table(resultFile, header = F, sep = " ")
        
        communityList$membership.multi <- wmemb_membership
        
        #if(!hardPartitions){
        #same columns regardless of this flag
        colnames(communityList$membership.multi) <-
            c("layer", "node", "module", "flow")
        #}
        #reorder, for easier inspection
        communityList$membership.multi <-
            communityList$membership.multi[order(communityList$membership.multi$layer,
                                                 communityList$membership.multi$node), ]
        
        
        resultFile <- paste0(outputPrefix, "_expanded.map")
        wtcod <-
            as.numeric(strsplit(readLines(resultFile, n = 5), " ")[[5]][3])
        
        communityList$codelength.multi <- wtcod
        
        cat(paste("    Code length Multiplex: ", wtcod, "\n"))
        numComms <- max(wmemb_membership$V3)
        cat(paste("    Communities Multiplex: ", numComms, "\n"))
        
        communityList$modules.multi <- numComms
        
        communityList$msize.multi <-
            table(communityList$membership.multi$module)
        
        #depending on flags, Infomap can transform into layer IDs the id of isolated nodes.
        #let's remove those ones
        communityList$membership.multi <-
            communityList$membership.multi[which(communityList$membership.multi$layer <=
                                                     Layers), ]
        
        
        #TODO for the future: calculate modularity of the partition. No direct multiplex way from igraph
        #one possibility is to pass the expanded representation of the network
        #but in case of edgecolored the supradjacency matrix would empty off-diagonal
        #resulting in huge modularity due to layers, not modules..
        #igraph::modularity(x, membership, weights = NULL, ...)
        
        if (addAggregateAnalysis) {
            cat(' + Aggregate network...\n')
            
            #calculate same things for the aggregate using R-igraph infomap
            g.agg <- GetAggregateNetworkFromNetworkList(g.list)
            infocom <-
                igraph::cluster_infomap(g.agg, modularity = TRUE)
            wmemb_membership_aggregate <-
                as.numeric(igraph::membership(infocom))
            wtcod_aggregate <- igraph::code_len(infocom)
            
            communityList$membership.aggr <-
                data.frame(node = 1:length(wmemb_membership_aggregate),
                           module = wmemb_membership_aggregate)
            communityList$codelength.aggr <- wtcod_aggregate
            
            cat(paste("    Code length Aggregate: ", wtcod_aggregate, "\n"))
            numCommsAggr <- max(wmemb_membership_aggregate)
            cat(paste("    Communities Aggregate: ", numCommsAggr, "\n"))
            
            communityList$modules.aggr <- numCommsAggr
            communityList$msize.aggr <-
                table(communityList$membership.aggr$module)
        }
        
        cat('Calculation Completed!\n')
        
        return(communityList)
    }
