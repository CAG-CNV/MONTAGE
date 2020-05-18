use warnings;
#use strict;
@ARGV>=1 or die "Usage: $0 <input-file>\n the input file contains one BAF LRR signalfile int.csv path per line
";
#use Getopt::Long;
#GetOptions('a'=>\ my $a);
open (FAM, $ARGV[0]) or die;
open (LOG, ">$ARGV[0]"."_Montage.log");
while (<FAM>) {
	s/[\r\n]+$//;
	#Check Column Order
	open(BafLrrFile, $_);
			$BafLrrHeader=<BafLrrFile>;
                        $BafLrrHeader =~ s/[\r\n]+$//;
                        @BafLrrHeaderLine=split(/\t/,$BafLrrHeader);
                        for $z (0 .. @BafLrrHeaderLine-1)
                        {
                                if ($BafLrrHeaderLine[$z] eq 'Name' or $BafLrrHeaderLine[$z] eq 'SNP' or $BafLrrHeaderLine[$z] eq 'SNP ID' or $BafLrrHeaderLine[$z] eq 'ProbeID')
                                {
                                        $name_index = $z+1;#awk is 1 based not 0 based
                                }
                                elsif ($BafLrrHeaderLine[$z] eq 'Chr' or $BafLrrHeaderLine[$z] eq 'Chromosome')
                                {
                                        $chr_index = $z+1;
                                }
                                elsif ($BafLrrHeaderLine[$z] eq 'Position')
                                {
                                        $pos_index = $z+1;
                                }
                                elsif ($BafLrrHeaderLine[$z] =~ m/(.*)\.?Log R Ratio$/ or $BafLrrHeaderLine[$z] =~ m/(.*)\.?LRR$/)
                                {
                                        $logr_index = $z+1;
                                }
                                elsif ($BafLrrHeaderLine[$z] =~ m/(.*)\.?B Allele Freq(uency)?$/ or $BafLrrHeaderLine[$z] =~ m/(.*)\.?BAF$/)
                                {
                                        $b_index = $z+1;
                                }

                        }
                        print LOG "Indexes: $name_index\t$chr_index\t$pos_index\t$logr_index\t$b_index\n";

        $command = "sed 's/\r//' $_ | grep -v REMOVE | sort -k$chr_index,$chr_index -k$pos_index,$pos_index"."n | awk '{if(\$$pos_index>0)print}' | awk -v OFS=\"\t\" -v ORS=\"\" 'BEGIN{window=1000000;slide=1000000;} {lastChr=chr; chr=\$$chr_index; lastMod=mod; mod=\$$pos_index%window; if(\$$pos_index<=window&&\$$chr_index==lastChr){count++}
else{sum2-=\$$logr_index;if(\$$b_index>=0.4&&\$$b_index<=0.6){sumMid--}else if(\$$b_index>=0.1&&\$$b_index<0.4){sumLow--}else if(\$$b_index>0.6&&\$$b_index<=0.9){sumHig--}else{sumHom--};sumBaf-=\$$b_index;sumsqBaf-=((\$$b_index)^2);if(\$$logr_index<-3){sumHomDel--}} 

if(\$$b_index>=0.4&&\$$b_index<=0.6){sumMid++}else if(\$$b_index>=0.1&&\$$b_index<0.4){sumLow++}else if(\$$b_index>0.6&&\$$b_index<=0.9){sumHig++}else{sumHom++};sum2+=\$$logr_index;if(\$$b_index>=0.1&&\$$b_index<=0.9){sumBaf+=\$$b_index;sumsqBaf+=(\$$b_index)^2};array[mod]=\$$b_index;array2[mod]=\$$logr_index;diff=mod-lastMod;if(\$$logr_index<-3){sumHomDel++}}

diff<0{print \"$_\",\$2\".\";printf \"\%.0f\",(\$$pos_index-slide)/slide;print \"000000\\t\"sumMid+0,sumLow+0,sumHig+0,sumHom+0,sum2/count; if(sumMid+sumLow+sumHig>0&&((sumsqBaf-sumBaf^2/( sumMid+sumLow+sumHig))/( sumMid+sumLow+sumHig))>0){print \"\\t\"sqrt((sumsqBaf-sumBaf^2/( sumMid+sumLow+sumHig))/( sumMid+sumLow+sumHig))} else{print \"\\t0\"}; print sumHomDel+0,\"\\n\";

sum=0;sum2=\$$logr_index;count=1;if(\$$chr_index==lastChr){window+=slide}else{window=1000000};sumMid=0;sumLow=0;sumHig=0;sumHom=0;sumHomDel=0;if(\$$b_index>=0.4&&\$$b_index<=0.6){sumMid++}else if(\$$b_index>=0.1&&\$$b_index<0.4){sumLow++}else if(\$$b_index>0.6&&\$$b_index<=0.9){sumHig++}else{sumHom++};sumBaf=\$$b_index;sumsqBaf=(\$$b_index)^2;if(\$$logr_index<-3){sumHomDel++}}' > Count.5LowHigh_Sliding.txt";
            print LOG "$command\n";
	    `$command`;
	   $command="BafSdThresh=\$(sed 's/\\ /\\t/g' Count.5LowHigh_Sliding.txt | egrep -v 'X|Y' | ./datamash q3 8 iqr 8 | sed 's/\\t/+1.5*/' | bc -l);LrrAvgThresh=\$(sed 's/\\ /\\t/g' Count.5LowHigh_Sliding.txt | egrep -v 'X|Y' | ./datamash q1 7 iqr 7 | sed 's/\\t/-1.5*/' | bc -l); echo \$BafSdThresh \$LrrAvgThresh; awk -v b=\"\$BafSdThresh\" -v l=\"\$LrrAvgThresh\" '{if(\$8>b&&\$7<l&&\$9<3)print}' Count.5LowHigh_Sliding.txt | egrep -v 'X|Y' > Candidates";
	    print LOG "$command\n";
            $o=`$command`;print LOG "$o";
	    $command = "awk '{print \"chr\"\$2,\"numsnp=1 length=1 state2,cn=1\",\$1,\"startsnp=\"\$2,\"endsnp=\"\$2}' Candidates | awk -F'[. ]' '{print \$1\":\"\$2\"-\"\$2+1000000,\$3,\$4,\$5,\$6\".\"\$7\".\"\$8,\$9\"_\"\$2,\$11\"_\"\$2+1000000}' > Candidates.rawcnv";
	     print LOG "$command\n";
            `$command`;
	    $command="BafSdThresh=\$(sed 's/\\ /\\t/g' Count.5LowHigh_Sliding.txt | egrep -v 'X|Y' | ./datamash q3 8 iqr 8 | sed 's/\\t/+1.5*/' | bc -l);LrrAvgThresh=\$(sed 's/\\ /\\t/g' Count.5LowHigh_Sliding.txt | egrep -v 'X|Y' | ./datamash q3 7 iqr 7 | sed 's/\\t/+1.5*/' | bc -l); echo \$BafSdThresh \$LrrAvgThresh; awk -v b=\"\$BafSdThresh\" -v l=\"\$LrrAvgThresh\" '{if(\$8>b&&\$7>l&&\$9<3)print}' Count.5LowHigh_Sliding.txt | egrep -v 'X|Y' > Candidates";
	    print LOG "$command\n";
            $o=`$command`;print LOG "$o";
	    $command = "awk '{print \"chr\"\$2,\"numsnp=1 length=1 state5,cn=3\",\$1,\"startsnp=\"\$2,\"endsnp=\"\$2}' Candidates | awk -F'[. ]' '{print \$1\":\"\$2\"-\"\$2+1000000,\$3,\$4,\$5,\$6\".\"\$7\".\"\$8,\$9\"_\"\$2,\$11\"_\"\$2+1000000}' >> Candidates.rawcnv";
	     print LOG "$command\n";
            `$command`;
		$command = "echo -e 'Name\\tChr\\tPos' > Header; awk '{ORS=\"\";temp=\$2;gsub(/\\./,\"_\",\$2);print \$2\"\\t\";gsub(/_/,\"\\t\",\$2);print \$2\"\\n\";gsub(/.*\\t/,\"\",\$2);gsub(/\\..*/,\"\",temp);print temp\"_\"\$2+1000000\"\\t\"temp\"\\t\"\$2+1000000\"\\n\";}' Count.5LowHigh_Sliding.txt | sort -u > Content; cat Header Content > signalfile";
		print LOG "$command\n";
            `$command`;
		$command = "perl ./clean_cnv.pl combineseg Candidates.rawcnv --signalfile signalfile >> $ARGV[0]"."_Montage.cleancnv 2>clean_cnv_messages";
	    print LOG "$command\n";
            `$command`;
		$command = "rm Header Content signalfile clean_cnv_messages Count.5LowHigh_Sliding.txt Candidates Candidates.rawcnv";
		print LOG "$command\n";
		`$command`;
    	 
	}
