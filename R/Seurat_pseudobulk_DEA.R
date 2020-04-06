# Function to separate single-cell data by a factor of interest and create 
# pseudo-bulk datasets by combining matrices for each factor
# followed by DEA with DESeq2   

# currently supports only single-factor designs
  
bulkDEA = function(object, #seurat object
                   metaData, #metadata mapping samples to factors. Rownames are sample names.
                   design, #design string for DESeq
                   splitCells = T, #whether or not to split by cell type (or whatever)
                   cellFactor = NULL, #name of the cell-type column
                   organism, #in ensembl format, e.g. mmusculus
                   geneVar, #type of gene ID used, in biomart format e.g. ensembl_gene_id
                   outDir, #directory to output results
                   outPrefix, #file prefix for results files
                   cores = 1) #cores to run in parallel
{
   require(Seurat)
   require(DESeq2)
   require(biomaRt)
   require(future)
   library(BiocParallel)
   register(BPPARAM = MulticoreParam(cores))
   dateString = format(Sys.Date(), "%y%m%d")
   if(!dir.exists(outDir))
      dir.create(outDir)
   setwd(outDir)
   
   compFactor = colnames(metaData)[1] # the factor we will be comparing by
   
   #make mart for conversion
   mart = useMart("ensembl", paste0(organism, "_gene_ensembl"))
   conv = getBM(attributes = c("ensembl_gene_id", "entrezgene_id", "external_gene_name"), mart = mart)
   
   # keep track of number of samples per condition
   completeMetadata = metaData
   allSamples = sort(unique(substr(rownames(object@meta.data), 1, (regexpr("\\_[ATCG]", rownames(object@meta.data)) - 1))))
   sampleData = matrix(nrow = 0, ncol = (length(allSamples) + 1)) # keeps track of number of cells per condition
   sampleDataCols = c("cellType", paste0("Ncells_", allSamples))
   colnames(sampleData) = sampleDataCols
   sampleData = as.data.frame(sampleData)
   
   bulkSampleData = matrix(nrow = 0, ncol = (length(unique(metaData[ , 1])) + 1)) #keep track of number of biological replicates per condition
   bulkSampleDataCols = c("cellType", paste0("Nsamples_", sort(unique(metaData[ , 1]))))
   colnames(bulkSampleData) = bulkSampleDataCols
   bulkSampleData = as.data.frame(bulkSampleData)
   
   
   cellTypes = c("AllCells", levels(factor(object@meta.data[ , cellFactor])))
   for(cell in cellTypes)
   {
      message(cell)
      metaData = completeMetadata #reset for each round
      
      cellSub = object
      if(cell != "AllCells" && splitCells) #for all cells, keep whole dataset
         cellSub = cellSub[ , cellSub@meta.data[ , cellFactor] == cell]
      
      #now get joined counts matrices
      counts = as.matrix(cellSub@assays$RNA@counts) # raw counts
      sd = cell
      for(sample in allSamples) #get number of cells for each (zero therefore means not present)
      {
            sd = c(sd, sum(grepl(sample, colnames(cellSub))))
      }
      samples = allSamples[as.numeric(sd[2:length(sd)]) > 50] # need at least 50 cells to get a real picture
      sampleData$cellType = as.character(sampleData$cellType)
      sampleData = rbind(sampleData, sd)
      colnames(sampleData) = sampleDataCols #necessary on first pass (but harmless afterward)
      
      # make sure comparison is even worth it
      metaData = subset(metaData, rownames(metaData) %in% samples)
      if(length(unique(metaData[, 1])) < 2) # skip if we can't make any comparisons
      {
         warning("Skipped -- not enough samples")
         next
      }      
      #remove factors with fewer than 2 samples
      bulkN = table(metaData[, 1])
      bulkN = data.frame(sample = names(bulkN), N = as.numeric(bulkN))
      bulkN = bulkN[order(bulkN$sample), ]
      keep = bulkN[bulkN$N >= 2, ]$sample
      if(length(keep) > 1)
      {
         cellSub = cellSub[ , cellSub@meta.data[ , compFactor] %in% keep]
      } else
      {
         warning("Skipped -- not enough samples")
         next #no need to continue with analysis in this case
      }
      
      #update sample Ns
      bulkN$N[bulkN$N < 2] = NA
      bulkSampleData$cellType = as.character(bulkSampleData$cellType) #necessary on first pass (but harmless afterward)
      bulkSampleData = rbind(bulkSampleData, c(cell, bulkN$N))
      colnames(bulkSampleData) = bulkSampleDataCols #necessary on first pass (but harmless afterward)
      
      joinedMat = matrix(data = NA, nrow = nrow(counts), ncol = length(samples))
      rownames(joinedMat) = rownames(counts)
      colnames(joinedMat) = samples
      for(sample in samples)
      {
         sampleMat = counts[ , grepl(sample, colnames(counts))]
         sampleCounts = rowSums(sampleMat)
         joinedMat[, sample] = sampleCounts
      }
      
      #now for DEA
      stopifnot(class(design) == "formula") #calling formula in function causes massive problems. Avoid by doing externally!
      joinedMat = joinedMat[, match(rownames(metaData), colnames(joinedMat))]
      des = DESeqDataSetFromMatrix(countData = joinedMat, colData = metaData, design = design)
      saveRDS(des, paste(dateString, outPrefix, cell, "subset", "des.rds", sep = "_"))
      dge = DESeq(des, parallel = T)
      
      allFactors = sort(unique(as.character(metaData[, 1])), decreasing = T) # reverse order so we get older/younger for all comps
      allComps = combn(x = allFactors, m = 2, simplify = F) # every possible combination of ages; list of older, younger for each comp
      allComps = lapply(allComps, as.character) #otherwise we get factor levels!
      for(comp in allComps)
      {
         res = as.data.frame(results(dge, contrast = c(colnames(metaData)[1], comp[1], comp[2]), alpha = 0.05))
         res = rownames_to_column(res, var = geneVar)
         if(geneVar == "ensembl_gene_id")
         {
            res$ensembl_gene_id = ifelse(grepl("\\.", res$ensembl_gene_id),
                                         yes = substr(res$ensembl_gene_id, 1, (regexpr("\\.", res$ensembl_gene_id) - 1)),
                                         no = res$ensembl_gene_id) #remove version numbers
         }
         res = merge(res, conv, all.x = T)
         
         compName = paste(comp[1], "vs", comp[2], sep = "_")
         write.csv(res, paste(dateString, outPrefix, cell, "subset", compName, "dge.csv", sep = "_"))
      }
   }
   
   #output sample statistics
   write.csv(sampleData, paste(dateString, outPrefix, "sampleData.csv", sep = "_"))
   write.csv(bulkSampleData, paste(dateString, outPrefix, "bulkSampleData.csv", sep = "_"))
}
      
         
         