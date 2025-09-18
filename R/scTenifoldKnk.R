#' @export scTenifoldKnk
#' @importFrom Matrix Matrix rowMeans rowSums
#' @importFrom scTenifoldNet makeNetworks tensorDecomposition manifoldAlignment
#' @author Daniel Osorio <dcosorioh@lainz.tecnun.es>
#' @title scTenifoldKNK
#' @description Predict gene perturbations
#' @param countMatrix countMatrix
#' @param gKO gKO
#' @param qc A boolean value (TRUE/FALSE), if TRUE, a quality control is applied over the data.
#' @param qc_minLSize An integer value. Defines the minimum library size required for a cell to be included in the analysis.
#' @param qc_mtThreshold A decimal value between 0 and 1. Defines the maximum ratio of mitochondrial reads (mithocondrial reads / library size) present in a cell to be included in the analysis. It's computed using the symbol genes starting with 'MT-' non-case sensitive.
#' @param nc_nNet An integer value. The number of networks based on principal components regression to generate.
#' @param nc_nCells An integer value. The number of cells to subsample each time to generate a network.
#' @param nc_nComp An integer value. The number of principal components in PCA to generate the networks. Should be greater than 2 and lower than the total number of genes.
#' @param nc_symmetric A boolean value (TRUE/FALSE), if TRUE, the weights matrix returned will be symmetric.
#' @param nc_scaleScores A boolean value (TRUE/FALSE), if TRUE, the weights will be normalized such that the maximum absolute value is 1.
#' @param nc_lambda A continuous value between 0 and 1. Defines the multiplicative value (1-lambda) to be applied over the weaker edge connecting two genes to maximize the adjacency matrix directionality.
#' @param nc_q A decimal value between 0 and 1. Defines the cut-off threshold of top q% relationships to be returned.
#' @param td_K An integer value. Defines the number of rank-one tensors used to approximate the data using CANDECOMP/PARAFAC (CP) Tensor Decomposition.
#' @param td_maxIter An integer value. Defines the maximum number of iterations if error stay above \code{td_maxError}.
#' @param td_maxError A decimal value between 0 and 1. Defines the relative Frobenius norm error tolerance.
#' @param td_nDecimal An integer value indicating the number of decimal places to be used.
#' @param ma_nDim An integer value. Defines the number of dimensions of the low-dimensional feature space to be returned from the non-linear manifold alignment.
#' @param nCores An integer value. Defines the number of cores to be used.
#' @examples
#' # Loading single-cell data
#' scRNAseq <- system.file("single-cell/example.csv",package="scTenifoldKnk")
#' scRNAseq <- read.csv(scRNAseq, row.names = 1)
#'
#' # Running scTenifoldKnk
#' scTenifoldKnk(countMatrix = scRNAseq, gKO = 'G100', qc_minLSize = 0)

# scTenifoldKnk主函数：执行单细胞基因调控网络的虚拟敲除实验
# 该函数用于模拟基因敲除对单细胞基因调控网络的影响
scTenifoldKnk <- function(countMatrix, qc = TRUE, gKO = NULL, qc_mtThreshold = 0.1, qc_minLSize = 1000, nc_lambda = 0, nc_nNet = 10, nc_nCells = 500, nc_nComp = 3,
                          nc_scaleScores = TRUE, nc_symmetric = FALSE, nc_q = 0.9, td_K = 3, td_maxIter = 1000,
                          td_maxError = 1e-05, td_nDecimal = 3, ma_nDim = 2, nCores = parallel::detectCores()){
  
  # 第1步：数据质量控制
  # 如果qc参数为TRUE，则对输入的表达矩阵进行质量控制
  # 包括过滤低质量细胞和线粒体基因比例过高的细胞
  if(isTRUE(qc)){
    # 调用scQC函数进行质量控制
    # mtThreshold: 线粒体基因比例阈值
    # minLSize: 最小文库大小阈值
    countMatrix <- scQC(countMatrix, mtThreshold = qc_mtThreshold, minLSize = qc_minLSize)
  }
  
  # 第2步：基因过滤
  # 根据细胞数量决定基因过滤策略
  if(ncol(countMatrix) > 500){
    # 如果细胞数大于500，保留在至少5%细胞中表达的基因
    countMatrix <- countMatrix[rowMeans(countMatrix != 0) >= 0.05,]
  } else {
    # 如果细胞数小于等于500，保留在至少25个细胞中表达的基因
    countMatrix[rowSums(countMatrix != 0) >= 25,]
  }
  
  # 第3步：构建野生型(WT)基因调控网络
  # 使用scTenifoldNet包的makeNetworks函数构建多个网络
  # 这一步会生成基于主成分回归的多个网络
  #set.seed(1)  # 注释掉的随机种子设置
  WT <- scTenifoldNet::makeNetworks(X = countMatrix, q = nc_q, nNet = nc_nNet, nCells = nc_nCells, 
                                   scaleScores = nc_scaleScores, symmetric = nc_symmetric, 
                                   nComp = nc_nComp, nCores = nCores)
  
  # 第4步：张量分解
  # 使用CANDECOMP/PARAFAC (CP)张量分解方法降噪和整合多个网络
  # 这一步将多个网络整合成一个去噪后的权重平均网络
  #set.seed(1)  # 注释掉的随机种子设置
  WT <- scTenifoldNet::tensorDecomposition(xList = WT, K = td_K, maxError = td_maxError, 
                                          maxIter = td_maxIter, nDecimal = td_nDecimal)
  
  # 第5步：执行虚拟基因敲除
  # 从张量分解结果中提取网络矩阵
  WT <- WT$X
  
  # 应用严格方向性约束，增强网络的方向性
  WT <- strictDirection(WT, lambda = nc_lambda)
  
  # 转换为标准矩阵格式
  WT <- as.matrix(WT)
  
  # 将对角线元素设为0（基因与自身的连接）
  diag(WT) <- 0
  
  # 转置矩阵以符合后续分析需求
  WT <- t(WT)
  
  # 创建敲除(KO)网络：复制野生型网络
  KO <- WT
  
  # 执行虚拟敲除：将目标基因的所有出度边设为0
  # 这模拟了该基因被敲除后无法调控其他基因的情况
  KO[gKO,] <- 0
  
  # 第6步：流形对齐
  # 使用非线性流形对齐方法比较野生型和敲除网络
  # 将高维网络数据映射到低维特征空间进行比较
  MA <- manifoldAlignment(WT, KO, d = ma_nDim, nCores = nCores)
  
  # 第7步：差异调控分析
  # 基于流形对齐结果识别受基因敲除影响的差异调控基因
  # 计算每个基因在两种条件下的距离和统计显著性
  DR <- dRegulation(MA, gKO)
  
  # 第8步：准备输出结果
  # 创建包含所有分析结果的列表
  outputList <- list()
  
  # 保存张量网络结果（野生型和敲除型网络）
  outputList$tensorNetworks$WT <- Matrix(WT)  # 野生型网络
  outputList$tensorNetworks$KO <- Matrix(KO)  # 敲除型网络
  
  # 保存流形对齐结果
  outputList$manifoldAlignment <- MA
  
  # 保存差异调控分析结果
  outputList$diffRegulation <- DR
  
  # 返回完整的分析结果
  return(outputList)
}
