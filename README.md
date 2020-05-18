# MONTAGE
Mosaic CNV Detection Tool  
MONTAGE: MOsaicNumberTAGEasily  
perl Montage.pl list_BafLrr_SignalFiles.txt  
A BafLrr SignalFile looks like this (column order is flexible):  
Name    Chr     Position        sample.B Allele Freq    sample.Log R Ratio  
rs1000000       12      126890980       0.9949  0.03895  
rs1000002       3       183635768       1       -0.3053  
# User Guide  
Input: Single Sample BAF/LRR Signal Files  
To derive from GenomeStudio Table use kcolumn.pl  
Name Chr Position sample1.B Allele Freq sample1.Log R Ratio … sampleN.B Allele Freq sampleN.Log R Ratio  
rs1000000 12 126890980 0.9949 0.03895 … 0.0129 -0.08263  
…  
To derive from GenomeStudio FinalReport use FinalReportToBafLrr.pl  
Name Chr Position B Allele Freq Log R Ratio sampleID  
rs1000000 12 126890980 0.9949 0.03895 sample1  
…  
rs1000000 12 126890980 0.0129 -0.08263 sampleN  
…  
Install: As long as you have standard bash and perl, you should not have to worry about dependencies or admin rights. If on Windows, you can use Cygwin.  
wget https://github.com/CAG-CNV/MONTAGE/archive/master.zip  
Unzip master.zip  
Parallel Computing: Recommended, there are no between sample dependencies in this Mosaic CNV calling.  
Output:  
chr6:57000000-62000000        numsnp=4      length=5,000,001   state2,cn=1 6092196049_R03C02.int.csv startsnp=6_57000000 endsnp=6_62000000  
Use visualize_cnv.pl to plot the BAF/LRR underlying the Mosaic CNV call.  
Use ParseCNV to conduct CNV association studies  
