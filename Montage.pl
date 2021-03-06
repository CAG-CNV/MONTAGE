use warnings;
#use strict;
@ARGV>=1 or die "Usage: $0 <input-file>\n the input file contains one BAF LRR signalfile int.csv path per line
optional --signalposfile map.txt needs to be specified if only SNP Name BAF LRR included in signalfiles without Chromosome and Position.
";
%MonthPastStart = ("0","Jan","1","Feb","2","Mar","3","Apr","4","May","5","Jun","6","Jul","7","Aug","8","Sep","9","Oct","10","Nov","11","Dec");
@timeData = localtime(time);
$Year = $timeData[5] + 1900;
if ($timeData[2] > 12){$timeData[2] -= 12; $AmOrPm = "PM";}
elsif ($timeData[2] eq 12){$AmOrPm = "PM";}else{$AmOrPm = "AM";}
if($timeData[2]<10) {$timeData[2]="0".$timeData[2];}
if($timeData[1]<10) {$timeData[1]="0".$timeData[1];}
if($timeData[0]<10) {$timeData[0]="0".$timeData[0];}
print "Run Started ".$timeData[2].":".$timeData[1].":".$timeData[0]." $AmOrPm on ".$MonthPastStart{$timeData[4]}." ".$timeData[3]." ".$Year."\n";
use Getopt::Long;
GetOptions('signalposfile=s'=>\$signalposfile);
if($signalposfile)
{open(MAP, $signalposfile) or die;
$mapHeaderLine = <MAP>;
$mapHeaderLine =~ s/[\r\n]+$//;
@mapHeader=split(/\t/,$mapHeaderLine);
for $z (0 .. @mapHeader-1)
{
	if ($mapHeader[$z] eq 'Name' or $mapHeader[$z] eq 'SNP' or $mapHeader[$z] eq 'SNP ID' or $mapHeader[$z] eq 'ProbeID')
        {
		$map_name_index = $z+1;#awk is 1 based not 0 based
	}
	elsif ($mapHeader[$z] eq 'Chr' or $mapHeader[$z] eq 'Chromosome')
        {
                #$map_chr_index = $z+1;
        }
        elsif ($mapHeader[$z] eq 'Position')
        {
        	#$map_pos_index = $z+1;
	}
}
}
open (FAM, $ARGV[0]) or die;
open (LOG, ">$ARGV[0]"."_Montage.log");
$filelist_or_signalHeaderLine = <FAM>;
$filelist_or_signalHeaderLine =~ s/[\r\n]+$//;
@filelist_or_signalHeader=split(/\t/,$filelist_or_signalHeaderLine);
if($#filelist_or_signalHeader > 1)
{
	close(FAM);
	open (FILELIST, ">$ARGV[0]"."_MontageInputFiles.txt");
	for($countInputs=0;$countInputs<=$#ARGV;$countInputs++)
	{
		print FILELIST "$ARGV[$countInputs]\n";
	}
	close(FILELIST);
	open(FAM, "$ARGV[0]"."_MontageInputFiles.txt");
}
else
{
	close(FAM);
	open (FAM, $ARGV[0]);
}
while ($filename=<FAM>) {
	$filename=~s/[\r\n]+$//;
	#Check Column Order
	open(BafLrrFile, $filename);
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
			if($chr_index && $pos_index)
			{
                        	print LOG "Indexes: $name_index\t$chr_index\t$pos_index\t$logr_index\t$b_index\n";
			}
			elsif($signalposfile)
			{
				print LOG "Indexes: $name_index\t$logr_index\t$b_index\n";
			}
			else
			{
				print "ERROR: No Chromosome and Position fields in input file $_. Please specify --signalposfile map.txt";
			}
if($signalposfile)
{
	$c="sort -k $map_name_index,$map_name_index $signalposfile >a";`$c`;
	$c="sort -k $name_index,$name_index $filename >b";`$c`;
	$command = "join a b -1 $map_name_index -2 $name_index -t \$'\\t' > $ARGV[0]"."tempInput.txt";
	print LOG "$command\n";
        $o=`$command`;print LOG "$o";
	$_=$ARGV[0]."tempInput.txt";

	        	open(BafLrrFile, $_);
                        $BafLrrHeader=<BafLrrFile>;
                        $BafLrrHeader =~ s/[\r\n]+$//;
                        @BafLrrHeaderLine=split(/\t/,$BafLrrHeader);
                        for $z (0 .. @BafLrrHeaderLine-1)
                        {
                                if ($BafLrrHeaderLine[$z] eq 'Name' or $BafLrrHeaderLine[$z] eq 'SNP' or $BafLrrHeaderLine[$z] eq 'SNP ID' or $BafLrrHeaderLine[$z] eq 'ProbeID')
                                {
                                        $name_index = $z+1;
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
                        if($chr_index && $pos_index)
                        {
                                print LOG "Indexes: $name_index\t$chr_index\t$pos_index\t$logr_index\t$b_index\n";
                        }
			else
			{
				print "ERROR: Insufficient fields in signalfile combined with signalposfile\n";
			}

}
else
{
	$_=$filename;
}
        $command = "chmod u+x datamash; median=`sed 's/\\r//' $_ | grep -v REMOVE | grep -v Position | awk '{print \$$logr_index}' | fgrep -v NaN | ./datamash median 1`; echo \$median; grep -v Position $_ | sort -C -k$chr_index,$chr_index -k$pos_index,$pos_index"."n; sorted=\$?; sed 's/\\r//' $_ | grep -v REMOVE | grep -v Position | if [ \$sorted == 0 ] ; then tee; else sort -k$chr_index,$chr_index -k$pos_index,$pos_index"."n; fi | awk '{if(\$$pos_index>0&&NF>=5)print}' | awk -v median=\$median '{if(\$$logr_index!=\"NaN\"){\$$logr_index=\$$logr_index-median}print}' | awk -v OFS=\"\t\" -v ORS=\"\" 'BEGIN{window=1000000;slide=1000000;} {lastChr=chr; chr=\$$chr_index; lastMod=mod; mod=\$$pos_index%window; if(\$$pos_index<=window&&\$$chr_index==lastChr){count++}
else{sum2-=\$$logr_index;if(\$$b_index>=0.4&&\$$b_index<=0.6){sumMid--}else if(\$$b_index>=0.1&&\$$b_index<0.4){sumLow--}else if(\$$b_index>0.6&&\$$b_index<=0.9){sumHig--}else{sumHom--};sumBaf-=\$$b_index;sumsqBaf-=((\$$b_index)^2);if(\$$logr_index<-3){sumHomDel--}} 

if(\$$b_index>=0.4&&\$$b_index<=0.6){sumMid++}else if(\$$b_index>=0.1&&\$$b_index<0.4){sumLow++; if(start == \"\"){start=\$$pos_index}; if(end<\$$pos_index&&\$$pos_index<(start+slide)){end=\$$pos_index}}else if(\$$b_index>0.6&&\$$b_index<=0.9){sumHig++; if(start == \"\"){start=\$$pos_index}; if(end<\$$pos_index&&\$$pos_index<(start+slide)){end=\$$pos_index}}else{sumHom++};sum2+=\$$logr_index;if(\$$b_index>=0.1&&\$$b_index<=0.9){sumBaf+=\$$b_index;sumsqBaf+=(\$$b_index)^2};array[mod]=\$$b_index;array2[mod]=\$$logr_index;diff=mod-lastMod;if(\$$logr_index<-3){sumHomDel++}}

diff<0{print \"$filename\",\$2\"#\";printf \"\%.0f\",(\$$pos_index-slide)/slide;print \"000000\\t\"sumMid+0,sumLow+0,sumHig+0,sumHom+0,sum2/count; if(sumMid+sumLow+sumHig>0&&((sumsqBaf-sumBaf^2/( sumMid+sumLow+sumHig))/( sumMid+sumLow+sumHig))>0){print \"\\t\"sqrt((sumsqBaf-sumBaf^2/( sumMid+sumLow+sumHig))/( sumMid+sumLow+sumHig))\"\\t\"} else{print \"\\t0\\t\"}; print sumHomDel+0,start+0,end+0\"\\n\";

sum=0;sum2=\$$logr_index;count=1;if(\$$chr_index==lastChr){window+=slide}else{window=1000000};sumMid=0;sumLow=0;sumHig=0;sumHom=0;sumHomDel=0;start=\"\";end=\"\";if(\$$b_index>=0.4&&\$$b_index<=0.6){sumMid++}else if(\$$b_index>=0.1&&\$$b_index<0.4){sumLow++; if(start == \"\"){start=\$$pos_index} end=\$$pos_index}else if(\$$b_index>0.6&&\$$b_index<=0.9){sumHig++; if(start == \"\"){start=\$$pos_index} end=\$$pos_index}else{sumHom++};sumBaf=\$$b_index;sumsqBaf=(\$$b_index)^2;if(\$$logr_index<-3){sumHomDel++}}' | awk '{if(\$10>0)print}' > $ARGV[0]"."Count.5LowHigh_Sliding.txt";
            print LOG "$command\n";
	    `$command`;
	   $command="chmod u+x datamash; BafSdThresh=\$(sed 's/\\ /\\t/g' $ARGV[0]"."Count.5LowHigh_Sliding.txt | egrep -v 'X|Y' | ./datamash q3 8 iqr 8 | sed 's/\\t/+1.5*/' | bc -l);LrrAvgThresh=\$(sed 's/\\ /\\t/g' $ARGV[0]"."Count.5LowHigh_Sliding.txt | egrep -v 'X|Y' | ./datamash q1 7 iqr 7 | sed 's/\\t/-1.5*/' | bc -l); echo \$BafSdThresh \$LrrAvgThresh; awk -v b=\"\$BafSdThresh\" -v l=\"\$LrrAvgThresh\" '{if(\$8>b&&\$7<l&&\$9<3)print}' $ARGV[0]"."Count.5LowHigh_Sliding.txt | egrep -v 'X|Y' > $ARGV[0]"."Candidates";
	    print LOG "$command\n";
            $o=`$command`;print LOG "$o";
	    $command = "awk '{ORS=\"\";print \"chr\"\$2,\$10,\$11,\"numsnp=1 length=1 state2,cn=1\";if(0.43922455<\$8&&\$8<=0.55321495){print\".1\"};
if(0.34220825<\$8&&\$8<=0.43922455){print\".2\"};
if(0.26099545<\$8&&\$8<=0.34220825){print\".3\"};
if(0.19441555<\$8&&\$8<=0.26099545){print\".4\"};
if(0.14129795<\$8&&\$8<=0.19441555){print\".5\"};
if(0.10047205<\$8&&\$8<=0.14129795){print\".6\"};
if(0.07076725<\$8&&\$8<=0.10047205){print\".7\"};
if(0.05101295<\$8&&\$8<=0.07076725){print\".8\"};
if(0.04003855<\$8&&\$8<=0.05101295){print\".9\"};
print\" \"\$1,\"startsnp=\"\$2,\"endsnp=\"\$2\"\\n\"}' $ARGV[0]"."Candidates | awk -F'[# ]' '{print \$1\":\"\$3\"-\"\$4,\$5,\$6,\$7,\$8,\$9\"_\"\$3,\$11\"_\"\$4}' > $ARGV[0]"."Candidates.rawcnv";
	     print LOG "$command\n";
            `$command`;
	    $command="BafSdThresh=\$(sed 's/\\ /\\t/g' $ARGV[0]"."Count.5LowHigh_Sliding.txt | egrep -v 'X|Y' | ./datamash q3 8 iqr 8 | sed 's/\\t/+1.5*/' | bc -l);LrrAvgThresh=\$(sed 's/\\ /\\t/g' $ARGV[0]"."Count.5LowHigh_Sliding.txt | egrep -v 'X|Y' | ./datamash q3 7 iqr 7 | sed 's/\\t/+1.5*/' | bc -l); echo \$BafSdThresh \$LrrAvgThresh; awk -v b=\"\$BafSdThresh\" -v l=\"\$LrrAvgThresh\" '{if(\$8>b&&\$7>l&&\$9<3)print}' $ARGV[0]"."Count.5LowHigh_Sliding.txt | egrep -v 'X|Y' > $ARGV[0]"."Candidates";
	    print LOG "$command\n";
            $o=`$command`;print LOG "$o";
	    $command = "awk '{ORS=\"\";print \"chr\"\$2,\$10,\$11,\"numsnp=1 length=1 state5,cn=3\";if(0.03667345<\$8&&\$8<=0.03974705){print\".1\"};
if(0.03974705<\$8&&\$8<=0.04808875){print\".2\"};
if(0.04808875<\$8&&\$8<=0.06052795){print\".3\"};
if(0.06052795<\$8&&\$8<=0.07589405){print\".4\"};
if(0.07589405<\$8&&\$8<=0.09301645){print\".5\"};
if(0.09301645<\$8&&\$8<=0.11072455){print\".6\"};
if(0.11072455<\$8&&\$8<=0.12784775){print\".7\"};
if(0.12784775<\$8&&\$8<=0.14321545){print\".8\"};
if(0.14321545<\$8&&\$8<=0.15565705){print\".9\"};
print\" \"\$1,\"startsnp=\"\$2,\"endsnp=\"\$2\"\\n\"}' $ARGV[0]"."Candidates | awk -F'[# ]' '{print \$1\":\"\$3\"-\"\$4,\$5,\$6,\$7,\$8,\$9\"_\"\$3,\$11\"_\"\$4}' >> $ARGV[0]"."Candidates.rawcnv";
	     print LOG "$command\n";
            `$command`;
	    $command="BafSdThresh=\$(sed 's/\\ /\\t/g' $ARGV[0]"."Count.5LowHigh_Sliding.txt | egrep -v 'X|Y' | ./datamash q3 8 iqr 8 | sed 's/\\t/+1.5*/' | bc -l);LrrAvgThreshLow=\$(sed 's/\\ /\\t/g' $ARGV[0]"."Count.5LowHigh_Sliding.txt | egrep -v 'X|Y' | ./datamash q1 7 iqr 7 | sed 's/\\t/-1.5*/' | bc -l);LrrAvgThreshHigh=\$(sed 's/\\ /\\t/g' $ARGV[0]"."Count.5LowHigh_Sliding.txt | egrep -v 'X|Y' | ./datamash q3 7 iqr 7 | sed 's/\\t/+1.5*/' | bc -l); echo \$BafSdThresh \$LrrAvgThreshLow \$LrrAvgThreshHigh; awk -v b=\"\$BafSdThresh\" -v ll=\"\$LrrAvgThreshLow\" -v lh=\"\$LrrAvgThreshHigh\" '{if(\$8>b&&\$7>=ll&&\$7<=lh&&\$9<3)print}' $ARGV[0]"."Count.5LowHigh_Sliding.txt | egrep -v 'X|Y' > $ARGV[0]"."Candidates";
	    print LOG "$command\n";
            $o=`$command`;print LOG "$o";
	    $command = "awk '{print \"chr\"\$2,\$10,\$11,\"numsnp=1 length=1 state4,cn=2\",\$1,\"startsnp=\"\$2,\"endsnp=\"\$2}' $ARGV[0]"."Candidates | awk -F'[# ]' '{print \$1\":\"\$3\"-\"\$4,\$5,\$6,\$7,\$8,\$9\"_\"\$3,\$11\"_\"\$4}' >> $ARGV[0]"."Candidates.rawcnv";
	     print LOG "$command\n";
            `$command`;
		$command = "echo -e 'Name\\tChr\\tPos' > $ARGV[0]"."Header; awk '{ORS=\"\";temp=\$2;gsub(/\\#/,\"_\",\$2);gsub(/_/,\"\\t\",\$2);gsub(/.*\\t/,\"\",\$2);gsub(/\\#.*/,\"\",temp);print temp\"_\"\$10\"\\t\"temp\"\\t\"\$10\"\\n\";print temp\"_\"\$11\"\\t\"temp\"\\t\"\$11\"\\n\";}' $ARGV[0]"."Count.5LowHigh_Sliding.txt | sort -u > $ARGV[0]"."Content; cat $ARGV[0]"."Header $ARGV[0]"."Content > $ARGV[0]"."signalfile";
		print LOG "$command\n";
            `$command`;
		$command = "perl ./clean_cnv.pl combineseg $ARGV[0]"."Candidates.rawcnv --signalfile $ARGV[0]"."signalfile 2>$ARGV[0]"."clean_cnv_messages | egrep -wv 'numsnp=1|numsnp=2' | sed 's/,cn=3/,cn=2/' >> $ARGV[0]"."_Montage.cleancnv";
	    print LOG "$command\n";
            `$command`;
		$command = "rm $ARGV[0]"."Header $ARGV[0]"."Content $ARGV[0]"."signalfile $ARGV[0]"."clean_cnv_messages $ARGV[0]"."Count.5LowHigh_Sliding.txt $ARGV[0]"."Candidates $ARGV[0]"."Candidates.rawcnv";
		print LOG "$command\n";
		`$command`;
    	 
	}
%MonthPastStart = ("0","Jan","1","Feb","2","Mar","3","Apr","4","May","5","Jun","6","Jul","7","Aug","8","Sep","9","Oct","10","Nov","11","Dec");
@timeData = localtime(time);
$Year = $timeData[5] + 1900;
if ($timeData[2] > 12){$timeData[2] -= 12; $AmOrPm = "PM";}
elsif ($timeData[2] eq 12){$AmOrPm = "PM";}else{$AmOrPm = "AM";}
if($timeData[2]<10) {$timeData[2]="0".$timeData[2];}
if($timeData[1]<10) {$timeData[1]="0".$timeData[1];}
if($timeData[0]<10) {$timeData[0]="0".$timeData[0];}
print "Run Ended ".$timeData[2].":".$timeData[1].":".$timeData[0]." $AmOrPm on ".$MonthPastStart{$timeData[4]}." ".$timeData[3]." ".$Year."\n";
