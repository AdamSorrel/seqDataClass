# seqDataClass
Data class to enable easy loading of sequencing output from DADA2, Deblur, Qiime2 and other platforms into a multiindex pandas data frame connecting taxonomy, metadata and ASV/OTU table. 

This file is functioning, but still in development before it will be hopefully published. It can now automatically load and parse data from DADA2, Qiime2, Deblur and Mothur. Once successfully loaded, it serves as a wrapper for relative abundance stacked barplot convenience function using plotly, as well as several specific multiindex pandas operations that would otherwise be more complicated, but are practical for data analysis. This data class has also entirely open access to the multiindexed data frame for user convenience offering all pandas functionality. 

Comments are adhering to standards that should allow automatic documentation production once this program is finalized. 
