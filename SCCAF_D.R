SCCAF_D = function(param, batch_key='sampleID',span=0.3,python_home = Sys.which("python")) {
    reticulate::use_python(python_home)
    source("/home/feng_shuo/Rscript/benchmark.R")
    source("/home/feng_shuo/deconvolution/SCCAF-D-test/deconvolution1.R")
    source("/home/feng_shuo/deconvolution/SCCAF-D-test/Frame.R")
    source("/home/feng_shuo/deconvolution/Real_bulk/DWLS.R")
    # source('/home/feng_shuo/Rscript/framework1.R')
    #####
    if (!reticulate::py_module_available("SCCAF")) {
        stop("python module SCCAF does not seem to be installed; - try running 'pip install SCCAF'")
    }
    reticulate::source_python("/home/feng_shuo/Rscript/sccaf.py")
    reticulate::source_python("/home/feng_shuo/Rscript//scanpy_workflow.py")
    #####
    if (length(param) != 11) {
        print("Please check that all required parameters are indicated or are correct")
        print("Example usage for bulk deconvolution methods: 'Rscript Master_deconvolution.R baron none bulk TMM all nnls 100 none 1'")
        print("Example usage for single-cell deconvolution methods: 'Rscript Master_deconvolution.R baron none sc TMM TMM MuSiC 100 none 1'")
        stop()
    }
    flag = FALSE
    bulk = param[1]
    dataset = param[2]
    transformation = param[3]
    deconv_type = param[4]
    if (deconv_type == "bulk") {
        normalization = param[5]
        marker_strategy = param[6]
    }
    else if (deconv_type == "sc") {
        normalization_scC = param[5]
        normalization_scT = param[6]
    }
    else {
        print("Please enter a valid deconvolution framework")
        stop()
    }
    method = param[7]
    number_cells = round(as.numeric(param[8]), digits = -2)
    to_remove = param[9]
    num_cores = min(as.numeric(param[10]), parallel::detectCores() - 
        1)
    if (param[11] == "T") {
        NormTrans = TRUE
    }
    else {
        NormTrans = FALSE
    }
    ###
    pathway=getwd()
    ####
    dataset1=strsplit(dataset,'\\.')[[1]][1]
    ####
    X2_1 = readRDS(dataset)
    sceasy::convertFormat(X2_1, from = "seurat", to = "anndata",main_layer = "counts",outFile = paste(dataset1,'.h5',sep=''))
    ###
    ad=paste(pathway,'/',dataset1,'.h5',sep='')
    #####selection cells--scanpy
    # batch_key=batch_key
    # span=span
    scanpy_workflow(ad,pathway,batch_key=batch_key,span=span)
    ###selection top100 cells for each cell type
    reference = selection_cells(pathway,dataset1,dataset)
    ###
    print(table(reference$cellType))
    ####
    ####read data
    X1 = read_bulk(bulk)
    X2 = read_data(paste(pathway,'/',dataset1,'_',"sccaf-reference",".rds",sep=""))
    ####
    #####
    # path = paste(getwd(), "/sccaf-d-results", sep = "")
    #####
    if (FALSE) {
        X2 <- QC(X2)
    }
    to_keep = intersect(rownames(X1$data), rownames(X2$data))
    print(paste0("number of intersect features: ", length(to_keep)))
    ####
    X1$data = X1$data[to_keep, ]
    X2$data = X2$data[to_keep, ]
    Xtrain = prepare_train(X2$data, X2$original_cell_names)
    pDataC = X2$pData
    P <- as.data.frame(matrix(0.1, nrow = length(unique(colnames(X2$data))), 
        ncol = dim(X1$data)[2]))
    rownames(P) <- unique(colnames(X2$data))
    colnames(P) <- colnames(X1$data)
    Xtest <- list(T = X1$data, P = P)
    #####
    return(Framework(deconv_type, NormTrans, Xtest, Xtrain, normalization, 
        normalization_scT, normalization_scC, transformation, 
        marker_strategy, to_remove, Xtest$P, method, pDataC))
}
