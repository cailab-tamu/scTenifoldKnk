#' Plot a KO-centered subnetwork and (optionally) annotate with enrichment
#'
#' @title Plot KO network
#' @description Generate and plot a KO-centered subnetwork from the output of
#' `scTenifoldKnk`. The function selects genes with significant differential
#' regulation, extracts their interactions from the reconstructed WT network,
#' filters edges by weight quantile, and displays the network. When
#' `annotate = TRUE` the function queries enrichment databases and overlays
#' category pies on nodes and a legend of significant terms.
#'
#' @param X List. Output from `scTenifoldKnk` function.
#' @param gKO Character. Gene symbol of simulated knockout gene.
#' @param q Numeric. Edge-weight quantile used to threshold weak edges (default
#'  0.99).
#' @param annotate Logical. If TRUE, perform enrichment annotation and display
#'  pies on nodes (default TRUE).
#' @param nCategories Integer. Maximum number of enrichment categories to show
#'  in the legend when annotation is requested (default 20).
#' @param fdrThreshold Numeric. Adjusted p-value cutoff (FDR) for reporting
#'  enriched terms (default 0.05).
#'
#' @return Invisibly returns NULL. The primary purpose is plotting the network
#'  as a side effect.
#'
#' @examples
#' \dontrun{
#' res <- scTenifoldKnk(countMatrix, gKO = "G100")
#' plotKO(res, gKO = "G100")
#' }
#'
#' @export
#' @import enrichR
#' @importFrom grDevices hcl.colors rgb
#' @importFrom graphics legend par
#' @importFrom stats quantile
#' @importFrom reshape2 melt
plotKO <- function(X, gKO, q = 0.99, annotate = TRUE, nCategories = 20, fdrThreshold = 0.05){
  # gKO <- 'Trem2'
  # q = 0.99
  # annotate = TRUE
  # nCategories = 20
  # fdrThreshold = 0.05
  gList <- unique(c(gKO, X$diffRegulation$gene[X$diffRegulation$distance > 1e-10 & X$diffRegulation$p.adj < 0.05]))
  if(length(gList) > 0){
    sCluster <- as.matrix(X$tensorNetworks$WT[gList,gList])
    koInfo <- sCluster[gKO,]
    gList <- gList[!grepl('^mt-|^Rpl|^Rps',gList, ignore.case = TRUE)]
    sCluster[abs(sCluster) <= quantile(abs(sCluster), q)] <- 0
    sCluster[gKO,] <- koInfo
    diag(sCluster) <- 0
    sCluster <-  reshape2::melt(as.matrix(sCluster))
    colnames(sCluster) <- c('from', 'to', 'W')
    sCluster <- sCluster[sCluster$W != 0,]
    netPlot <- igraph::graph_from_data_frame(sCluster, directed = TRUE)
    dPlot <- igraph::centr_degree(netPlot)$res
    W <- rep(1,nrow(sCluster))
    sG   <- (names(igraph::V(netPlot))[dPlot > 1])[-1]
    W[sCluster$from %in% sG] <- 0.2
    W[sCluster$to %in% sG] <- 0.2
    W[sCluster$from %in% gKO] <- 1
    W[sCluster$from %in% gKO & sCluster$to %in% sG] <- 0.8
    set.seed(1)
    layPlot <- igraph::layout_with_fr(netPlot, weights = W)
    dPlot <- (dPlot/max(dPlot))*20
    if(isTRUE(annotate)){
      enrichFunction <- function(X, fdrThreshold = fdrThreshold){
        E <- enrichr(X, c('KEGG_2019_Human', 'GO_Biological_Process_2018','GO_Cellular_Component_2018', 'GO_Molecular_Function_2018','BioPlanet_2019', 'WikiPathways_2019_Human', 'Reactome_2016'))
        E <- do.call(rbind.data.frame, E)
        E <- E[E$Adjusted.P.value < fdrThreshold,]
        E <- E[order(E$Adjusted.P.value),]
        E$Term <- unlist(lapply(strsplit(E$Term,''), function(X){
          X[1] <- toupper(X[1])
          X <- paste0(X,collapse = '')
          X <- gsub('\\([[:print:]]+\\)|Homo[[:print:]]+|WP[[:digit:]]+','',X)
          X <- gsub("'s",'',X)
          X <- unlist(strsplit(X,','))[1]
          X <- gsub('[[:blank:]]$','',X)
          return(X)
        }))
        selectedSet <- rep(FALSE, nrow(E))
        for(i in seq_len(nrow(E))){
          if(i == 1){
            selectedSet[i] <- TRUE
          } else {
            A <- unique(unlist(strsplit(E[which(selectedSet[seq_len(i)]),'Genes'], ';')))
            B <- unlist(strsplit(E[i,'Genes'], ';'))
            selectedSet[i] <- !all(B %in% A)
          }
        }
        gSets <- table(toupper(E$Term))
        gSets <- names(gSets[gSets > 1])
        for(i in gSets){
          selectedSet[which(toupper(E$Term) %in% i)[-1]] <- FALSE
        }
        E <- E[selectedSet,]
        if(nrow(E) > nCategories){
          E <- E[seq_len(nCategories),]
        }
        return(E)
      }
      E <- enrichFunction(gList, fdrThreshold)
      if(isTRUE(nrow(E) > 0)){
        tPlot <- strsplit(E$Genes, ';')
        pPlot <- matrix(0,nrow = length(igraph::V(netPlot)), ncol = nrow(E))
        rownames(pPlot) <- toupper(names(igraph::V(netPlot)))
        for(i in seq_along(tPlot)){
          pPlot[unlist(tPlot[i]),i] <- 1
        }
        pPlot <- lapply(seq_len(nrow(pPlot)), function(X){as.vector(pPlot[X,])})
        names(pPlot) <- names(igraph::V(netPlot))
        tPlot <- unique(unlist(tPlot))
        eGenes <- toupper(names(igraph::V(netPlot))) %in% tPlot
        vColor <- rgb(195/255, 199/255, 198/255 ,0.3)
        if(nrow(E) == 1){
          pieColors <- list(hcl.colors(5, palette = 'Zissou 1', alpha = 0.7)[5])
        } else {
          pieColors <- list(hcl.colors(nrow(E), palette = 'Zissou 1', alpha = 0.7))
        }
        par(mar=c(4,0,0,0), xpd = TRUE)
        suppressWarnings(plot(netPlot,
                              layout = layPlot,
                              edge.arrow.size=.2,
                              vertex.label.color="black",
                              vertex.shape = ifelse(eGenes,'pie','circle'),
                              vertex.pie = pPlot,
                              vertex.size = 10+dPlot,
                              vertex.pie.color=pieColors,
                              vertex.label.family="Times",
                              vertex.label.font=ifelse(eGenes,2,1),
                              edge.color = ifelse(E(netPlot)$W > 0, 'red', 'blue'),
                              edge.curved = ifelse(W == 0.2, 0, 0.1),
                              vertex.color = vColor,
                              vertex.frame.color = NA))
        sigLevel <- formatC(E$Adjusted.P.value, digits = 2, format = 'g', width = 0, drop0trailing = TRUE)
        gSetNames <- lengths(strsplit(E$Genes, ';'))
        gSetNames <- paste0('(', gSetNames,') ', E$Term, ' FDR = ', sigLevel)
        legend(x = -1.05, y = -1.05, legend = gSetNames, bty = 'n', ncol = 2, cex = 1, col = unlist(pieColors), pch = 16)
      } else {
        par(mar=c(0,0,0,0))
        plot(netPlot,
             layout = layPlot,
             edge.arrow.size=.2,
             vertex.label.color="black",
             vertex.size = 10+dPlot,
             vertex.label.family="Times",
             edge.color = ifelse(E(netPlot)$W > 0, 'red', 'blue'),
             edge.curved = ifelse(W == 0.2, 0, 0.1),
             vertex.color = rgb(0,188/255,1,0.3),
             vertex.frame.color = NA)
      }
    } else {
      par(mar=c(0,0,0,0))
      plot(netPlot,
           layout = layPlot,
           edge.arrow.size=.2,
           vertex.label.color="black",
           vertex.size = 10+dPlot,
           vertex.label.family="Times",
           edge.color = ifelse(E(netPlot)$W > 0, 'red', 'blue'),
           edge.curved = ifelse(W == 0.2, 0, 0.1),
           vertex.color = rgb(0,188/255,1,0.3),
           vertex.frame.color = NA)
    }
  } else {
    netPlot <- matrix(0,1,1)
    rownames(netPlot) <- colnames(netPlot) <- gKO
    netPlot <- igraph::graph_from_adjacency_matrix(netPlot)
    par(mar=c(0,0,0,0))
    plot(netPlot,
         edge.arrow.size=.2,
         vertex.size = 50,
         vertex.label.color="black",
         vertex.label.family="Times",
         vertex.color = rgb(0,188/255,1,0.3),
         vertex.frame.color = NA)
  }
}
