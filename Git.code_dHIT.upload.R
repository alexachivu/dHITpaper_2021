#Code to upload on git: --> dHIT paper

### Packages used: rtracklayer, GenomicRanges, BRGenomics, DEseq2, ggplot.

#1 .Function: Imports BigWigs into R as GRanges object
#Parameters:
#           BW_plus.path  = path to positive/plus strand PROseq
#           BW_minus.path = path to negative/minus strand PROseq
#           coord.list    = GRanges object, list of coordinates to read PROseq within 
#Output: Granges of PROseq reads within the given coordinates list 

Import.BigWigs = function(BW_plus.path, BW_minus.path, coord.list = c() ){
    
    #Import PROseq BigWig files
        PRO_plus  = import.bw(BW_plus.path, which = coord.list);
        PRO_minus = import.bw(BW_minus.path, which = coord.list);

    #Set the strand
        strand(PRO_plus) = "+";
        strand(PRO_minus) = "-";
    
    #Make all PROseq score positive
        PRO_minus$score = (abs(PRO_minus$score));
        PRO_plus$score  = (abs(PRO_plus$score));

    
    #Merge PROseq plus and minus strands & convert to BRG
        PROseq = makeGRangesBRG( append( PRO_plus, PRO_minus ) );
        
    
    return(PROseq)
}

#2 .Function: Counts reads that fall within a given set of coordinates
#Parameters:
#           BW_plus.path  = path to positive/plus strand PROseq
#           BW_minus.path = path to negative/minus strand PROseq
#           coord.list    = GRanges object, list of coordinates to read PROseq within (all entries need to be resized to the same windth)
#           binsize       = integer, bin size to count signal within coord.list
#           field         = field corresponding to signal/score in PROseq GRanges (column in the PROseq GRanges object imported prom paths)
#           FUN           = function to use for binning the reads
#Output: Granges of PROseq reads within the given coordinates list, binned in bins of width equal to the binsize  
Count.reads = function( BW_plus.path, BW_minus.path, coord.list, binsize = 10, field = "score", FUN = "sum"){
    options(warn=-1)
    #Import PROseq BigWig files
        PRO_plus = import.bw(BW_plus.path, which = gene.list);
        PRO_minus = import.bw(BW_minus.path, which = gene.list);

    #Set the strand
        strand(PRO_plus) = "+";
        strand(PRO_minus) = "-";
    
    #Merge PROseq plus and minus strands & convert to BRG
        PROseq = makeGRangesBRG( append( PRO_plus, PRO_minus ) );
        PROseq$Norm.score = abs(PROseq$score);    
    
    #Count signal around TREs 
        PROseq.cts = getCountsByPositions(PROseq, gene.list, binsize = binsize, FUN = FUN, field = field);

    return(PROseq.cts)
    



#3 .Function: Runs DEseq2 on a list of PROseq paths
#Parameters:
#           PROseq.list  = list of PROseq paths (example of PROseq list: PROseq.list <- GRangesList(A_rep1 = Import.BigWigs("PROseq1_plus.bw", "PROseq1_minus.bw"),
#                                                                                                   A_rep2 = Import.BigWigs("PROseq2_plus.bw", "PROseq2_minus.bw"),)
#                                                )
#           ranges       = ranges of data to run DEseq2 on
#           spikeins.list = list of spike-in or scaling factors (if left empty, DEseq will run its intrinsic normalization of the data)
#           ncores        = number of cores to use for computation
#           field         = field to use for getting PROseq signal 
#Output: DEseq output of results
Run.DEseq2 = function( PROseq.list, ranges, spikeins.list = c(), ncores = 1, field = "score" ){       
        
    #Filter data before DEseq analysis
    dds = getDESeqDataSet( PROseq.list, ranges, gene_names = ranges$Gn.name, ncores = ncores, field = field);
    
    #Run DEseq2
    res.dds = getDESeqResults(dds, "B", "A", sizeFactors = spikeins.list, alpha = 0.05);

    #Add gene information to DEseq results
    res.dds_genes = merge(as.data.frame(res.dds), as.data.frame(ranges), by="row.names", all=TRUE);
    res.dds_genes = na.omit(as.data.frame(res.dds_genes)) # --> genes w/o NAs

    #Prepare the function output
    out = list(res.dds, res.dds_genes);

    return(out)
}

#4 .Function: Makes custom MA plot for given DEseq output data
#Parameters:
#           output_Run.DEseq2  = DEseq2 output (can be used with the Run.DEseq2 output)
#           Fig.type    = figure type to write the plot to (can be png, svg, bit, etc)
#           fig.width   = figure width 
#           fig.height  = figure height
#Output: MA plot with density map   
Plot.MAs = function(output_Run.DEseq2, file.type = "image/png", pl.height = 4, pl.width = 4){
    
    options(jupyter.plot_mimetypes = file.type, repr.plot.width = pl.width, repr.plot.height = pl.height);

    #Plot a generic MA graph
    plotMA( na.omit(output_Run.DEseq2[[1]]), alpha = 0.05 )
    
    #Plot a more complex MA, with heatmap
    ggplot( as.data.frame( na.omit(output_Run.DEseq2[[1]]) ), mapping = aes(x=baseMean, y=log2FoldChange) ) +
    geom_point() +
    scale_x_log10() +
    scale_y_continuous() +
    geom_pointdensity(adjust = 0.8, size = 1.5) +
        scale_color_viridis()+
    xlab("baseMean") +
    ylab("log2FC") +
    theme_classic() + 
    geom_hline(yintercept = 0, alpha=0.5, linetype="dashed", col = "red", size = 1)
    
    
}


#5. Heatmap function
LinearHeatmap = function( CountMatrix, nRows, FUN='mean', ... ) {
    CountMatrix = CountMatrix;
    cksz = floor(nrow(CountMatrix)/nRows);
    myit = iter(CountMatrix, by='row', chunksize=cksz);
    linMat = foreach( chunk = myit, .combine='rbind' ) %do% {
        apply(chunk, 2, FUN, ...);
    }
    return( linMat[1:nRows,] );
}

#6 .Function: Takes in a matrix of counts and plots a heatmap
#Parameters:
#           Signal.cts  = matrix of counts (can use output from Count.reads)
#           n.lines     = number of lines to plot in the heatmap (the smaller the number, the lower the resolution)
#           color       = heatmap desireg max color
#           Fig.type    = figure type to write the plot to (can be png, svg, bit, etc)
#           fig.width   = figure width 
#           fig.height  = figure height
#           Fig.title   = figure title
#           data.scale  = adjust outliers in data (all values > data.scale will be reset to data.scale)
#           l.max       = color limits for scale_fill_gradient2
#Output: Granges of PROseq reads within the given coordinates list, binned in bins of width equal to the binsize  
Plot.heatmap = function(Signal.cts, n.lines = 100, color = "#FF0000", Fig.type = "image/png", fig.width = 4, fig.height = 3, Fig.title = "", data.scale = 0.1, l.max = 0.1){
    
    #Figure characteristics
    options(jupyter.plot_mimetypes = Fig.type, repr.plot.width = fig.width, repr.plot.height = fig.height );
    theme_set(
            theme_classic() +

            theme(
                    plot.title = element_text(color="black", size=12, face="bold"),
                    axis.title.x = element_text(color="black", size=10, face="bold"),
                    axis.title.y = element_text(color="black", size=10, face="bold")
                    ) +
            theme(axis.title.y = element_text(margin = margin(t = 0, r = 30, b = 0, l = 0))) +
            theme(axis.title.x = element_text(margin = margin(t = 30, r = 0, b = 0, l = 0))) +
            theme(axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            axis.text.x = element_text(face="bold", color="black", size=14),
            axis.line = element_line(colour = "white")) +
            theme(plot.title = element_text(hjust = 0.5))
    );


    #############################################################################################################
    
    #Working with the data
    #Remove rows with no counts
    hmap = Signal.cts[rowMeans(Signal.cts)>0,];
    #Row-normalize signals
    hmap = hmap / apply(hmap, 1, max);

    linmap = LinearHeatmap(hmap, n.lines);
    
    #Set limits for the data
    linmap[linmap == "NaN"] = 0;
    linmap[linmap > data.scale] = data.scale;


    rownames(linmap) = NULL;
    df.0.BW <- reshape2::melt(linmap, varnames = c("TSS", "position"))

    #gg_plot = ggplot(df.0.BW, aes(x = position, y = TSS, fill = log10(1+10*value))) + 
    gg_plot = ggplot(df.0.BW, aes(x = position, y = TSS, fill = log10(1+value))) + 

              stat_density2d() + 
              geom_tile() +
              ggtitle(Fig.title) +
                        xlab("Distance from TSS (kb)") +
                        ylab(" ") +
              scale_fill_gradient2(
                  low = muted("red"),
                  mid = "white",
                  high = muted("blue"),
                  #midpoint = 0.01, 
                  #low="white", high = color,
                                    na.value = "grey50", guide = "colourbar",
                                   aesthetics = "fill", limits = c(0, l.max) ) 
    
    
    return(gg_plot)
    
}




#7.  Compute linear regression
#Define regions that background regions for the ChIP of interest (example for H3K27ac below).
# The rationale for this analysis follows the mathematical models presented in: Bonhoure et al., 2014 
#("Quantifying ChIP-seq data: a spiking method providing an internal reference for sample-to-sample normalization")



## I. `Remove regions that may contain H3K27ac signal` :

        To achieve this, I masked ed all coordinates corresponding to H3K27ac ENCODE broad peaks from hg19 (+ the adjacent regions)
               
               bedtools complement -i wgEncodeBroadHistoneK562H3k27acStdPk.broadPeak.sorted  -g ../hg19/hg19.sort.chromInfo > hg19.noENCODEpeaks.K27ac.bed



## II. `Split the untranscribed chromatin (hg19) into 5kb non-overlapping windows`:

                
                *bedtools makewindows -b hg19.noENCODEpeaks.K27ac.bed -w 5000 > hg19.noENCODE.5kb.nonOverlspWin.bed
                
            
       The output "hg19.noENCODE.5kb.nonOverlspWin.bed" can be found in:
       
       ../ChIPseq_Trp/NormUsing.ReadsIn.Over.ReadsOutsidePeaks_April2020/K27ac
      
<br>      
      
      - Output looks like:
                    chr1    0       5000
                    chr1    5000    10000
                    chr1    10000   11869    #--> Certain windows are < 5kb. 
                    chr1    31109   34554    # --> To address this, I will import this file into R and remove all winsows < 5kb
                    chr1    36081   41081    # --> A total of 284216 windows were removed
                    chr1    41081   46081
                    chr1    46081   51081
                    chr1    51081   52473
                    chr1    54936   59936
                    chr1    59936   62948
                    
<br> 

        The name of the file containing only untranscribed chromatin (w/o H3K27ac) binned in 5kb non-overlapping windows is: hg19.noENCODE.5kb.nonOverlspWin.fixed.bed

<br>
  
## III. `Count number of ChIP-seq reads falling in each 5kb window`:
       I took H3K27ac ChIP-seq as example and overlapped each ChIPseq BAM with hg19.5kb.nonOverlspWin.bed, counting signal in each bin.
      
                  ex: bedtools coverage -a hg19.noENCODE.5kb.nonOverlspWin.fixed.bed -b K27ac.bam* > K27ac.coverage5kb.nonOverlapWin.bed
                  
                  
                   bedtools coverage -a hg19.noENCODE.5kb.nonOverlspWin.fixed.bed -b ./mergedMasked_K27ac_0h_br1.bam > K27ac_0h_br1.coverage5kb.noENCODE.bed &          

<br>       
<br>

       - Output looks like:
       [agc94@cbsudanko WCE]$ head K27ac_0h_br1.coverage5kb.noENCODE.bed
                chr     start   end    width  strand  depth  nrBases_at_depth    width.of.window   proportion.Of.B.at.depth
                
                chr1    0       5000    5001    *       0            0                 5000                0.0000000
                chr1    5000    10000   5001    *       0            0                 5000                0.0000000
                chr1    36081   41081   5001    *       0            0                 5000                0.0000000
                chr1    41081   46081   5001    *       0            0                 5000                0.0000000
                chr1    46081   51081   5001    *       0            0                 5000                0.0000000
                chr1    54936   59936   5001    *       12           456               5000                0.0912000
                chr1    63887   68887   5001    *       12           407               5000                0.0814000
                chr1    70008   75008   5001    *       8            252               5000                0.0504000,
                
               where:
                     A feture                 = bed file of coordinates 
                     B feature                = BAM files
                     depth                    = nr reads in each window
                     nrBases_at_depth         = nr bases in each window
                     width.of.window          = width of each window
                     proportion.Of.B.at.depth = nrBases_at_depth/5000
       
           
<br>
<br>

# ***lot ChIPseq as a function of WCE signal in these background regions. Get linear regression coeficients to use for ChIPseq normalization.***





#8. Function: Count reads in binned windows (used for CUT&RUN)
#Parameters:
#           mark.BAM = GRanges object, Rdata of BAM file (each entry contains start and end coordinates for each read in the BAM)
#           regions  = regions of interest to count reads within
#           spike    = spike-in/ scaling factor (total number of spike-ins in a given sample) 
#           binsize  = integer, bin size to count signal within regions
#           field    = field corresponding to signal/score in PROseq GRanges (column in the PROseq GRanges object imported prom paths)
#Output: matrix of spike-in normalized reads within given regions 
Count.reads.InBins = function(mark.BAM, regions, spike = 1, binsize = 1, field = "score"){
    
    #Resize each read from BAM to 1bp (center position of each read)
    mark.BAM.res = resize(mark.BAM, width = 1, fix = "center");
    #Compute coverage
    mark.BAM.res = GRanges(coverage(mark.BAM.res, weight="count"));
    mark.BAM.res = makeGRangesBRG(mark.BAM.res[mark.BAM.res$score > 0]);

    #Count DHS signal within each CTCF peak
    mark.BAM.score = getCountsByPositions(mark.BAM.res, regions, binsize = binsize, field="score");
    
    return(mark.BAM.score/spike)
}



#9. Function: Count and normalize reads in binned windows (used for ChIPseq)
#Parameters:
#           ChIP.Rdata.path = path to the ChIPseq Rdata containing GRanges object of BAM file (each entry contains start and end coordinates for each read in the BAM)
#           WCE.Rdata.path  = path to the WCE Rdata containing GRanges object of BAM file (each entry contains start and end coordinates for each read in the BAM)
#           coord.path      = path to regions of interest to count reads within
#           up              = bases upstream from the TSS of coord.path
#           down            = bases downstream from the TSS of coord.path
#           binsize         = integer, bin size to count signal within regions
#           field           = field corresponding to signal/score in PROseq GRanges (column in the PROseq GRanges object imported prom paths)
#           a,b             = linear regression coeficients (obtained as decribed in 7.)
#           spike           = spike-in/ scaling factor (total number of spike-ins in a given sample) 
#Output: matrix of normalized ChIPseq reads within given regions (coord.path)
Norm.ChIP = function(ChIP.Rdata.path, WCE.Rdata.path, coord.path, up = 2000, down = 2000, binsize = 100, field = "score", 
                    a = 0.87326, b = 6.87326, spike = 3002){
    
    #Read in coordinates of interest
        promoters = read.table(coord.path, header = T);
    #Convert them to GRanges
        promoters = makeGRangesFromDataFrame(promoters.maxTSS);
    #Resize them to 4kb: each entries consists of the base-start positions of each maxTSS
        promoters = promoters(promoters, upstream=up, downstream=down);


    #Load ChIPseq BAM file
    load(ChIP.Rdata.path);
    ChIP = reads;
    rm(reads)

    #Load Whole cell extract (WCE) BAM file
    load(WCE.Rdata.path);
    WCE = reads;
    rm(reads)

    #Size select the data for mono-nucleosome fragments
    ChIP = ChIP[width(ChIP) >= 100 & width(ChIP) <= 200];
    WCE  = WCE[width(WCE) >= 100 & width(WCE) <= 200];


    #Count reads within given regions on interest 
    ChIP.perBin = Count.reads.InBins(ChIP, promoters, spike = 1,  binsize = binsize, field = field);
    WCE.perBin  = Count.reads.InBins(WCE,  promoters, spike = 1,  binsize = binsize, field = field);

    #Calculate theotic value fo ChIP signal from the linear regressions computed using (Lnear.reg function)
    ChIP.theoretic.perBin = (a * WCE.perBin) +  b;

    #Calculate residuals between the real number of reads and the theoretical one
    res = ChIP.perBin - ChIP.theoretic.perBin;

    #Replace all residuals < 0 with 0
    res[res < 0] = 0;

    #Subsample residual matrices to compute CI and apply the spike-in normalization 
    #(To normalize, divide by the total number of spike-ins in that given sample)
    res.subs = metaSubsampleMatrix(res/spike);

    return(res.subs)
}




#10. Function: Find ChIP-seq reads that fall within peaks 
#Parameters:
#           peaks.GRanges    = GRanges object containing coordinated of ChIPseq peaks called with macs2
#                              *This function requires the file encoding peaks coordinated (peaks.GRanges) to include a unique identified column with the name of each peak
#           ChIP.GRanges     = GRanges object containing ChIPseq signal (signal needs to have already been adjusted for differences in background)
#           spike            = spike-in/ scaling factor (total number of spike-ins in a given sample) 
#           sumReads.colName = column name to append with the number of reads in each peak
#Output: matrix of normalized reads
GetSignal.inPeaks = function(peaks.GRanges, ChIP.GRanges, spike = 1, sumReads.colName="New.column"){
        
        #Spike-in norm the BAM file
        ChIP.GRanges$score = ChIP.GRanges$score * 1/spike;
    
        hits = findOverlaps(ChIP.GRanges, peaks.GRanges);
        ChIP.GRanges$transcript_id = "nothing";
        ChIP.GRanges$transcript_id[hits@from] = as.character(peaks.GRanges$transcript_id[hits@to]);
        
        #Get number of reads per peak
        aggregated.peaks = aggregate(ChIP.GRanges$score ~ ChIP.GRanges$transcript_id, FUN="sum");
        colnames(aggregated.peaks) = c("transcript_id", "score");

       #Perform left-join to get the number of reads per peak 
        mrg.data2b <- merge(aggregated.peaks, as.data.frame(peaks.GRanges), by="transcript_id", all = TRUE);
   
        #If there were no reads for certain genes/peaks, the merge function will fill score = NA
        #So, I have replaced these NA values with 0
        mrg.data2b[is.na(mrg.data2b$score), ]$score  = 0;
        
        #Change name of transcript_id column to sumReads.colName
        names(mrg.data2b)[1]<-paste(sumReads.colName);
    
    return(mrg.data2b)
    
}
