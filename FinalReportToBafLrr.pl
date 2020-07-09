#!/usr/bin/perl
# Joe Glessner


unless( open( REPORTFILE, $ARGV[0]))
    {
        print "Final Report File not found! $ARGV[0]\n";
        exit;
    }

unless( open( EXCLUDEFILE, $ARGV[1]))
    {
        print "Exclude File not found! $ARGV[1] (fine if you have no samples to exclude)\n";
    }
while ($line = <EXCLUDEFILE>)
{
	chomp($line);
	$line =~ s/[\r\n]+$//;
	@line_array = split(/\t|\,/,$line);
	$Exclude{$line_array[0]}=1;
}
$exclude=0;

while (!($line =~ /SNP Name/))	#Take out header
{	
	$line = <REPORTFILE>;
}
	chomp($line);
	$line =~ s/[\r\n]+$//;
	@header_line = split(/\t|\,/,$line);
	for($i=0;$i<=$#header_line;$i++)
	{
		if($header_line[$i] =~ /SNP Name/)
		{
			$SNP_Name_Index=$i;
		}
		elsif($header_line[$i] =~ /Sample/)
		{
			$Sample_Index=$i;
		}
		elsif($header_line[$i] =~ /B Allele Freq/)
                {
                        $BAF_Index=$i;
                }
		elsif($header_line[$i] =~ /Log R Ratio/)
                {
                        $LRR_Index=$i;
                }
		else
		{
			print "WARNING: $header_line[$i] will not be used\n";
		}
		#print "$i\n";
	}
$i=0;
while (defined ($line = <REPORTFILE>))
{
	chomp($line);
	$line =~ s/[\r\n]+$//;
	@Vals=split(/\t|\,/,$line);
	if($lastID ne $Vals[$Sample_Index])
	{
		close RESULTFILE;
		if(exists($Exclude{$Vals[$Sample_Index]}))
		{
			$exclude++;
		}
		else
		{
		open (RESULTFILE, ">$Vals[$Sample_Index].baflrr");
		$i++;
		print RESULTFILE "Name\t".$Vals[$Sample_Index].".Log R Ratio\t".$Vals[$Sample_Index].".B Allele Freq\n"."$Vals[$SNP_Name_Index]\t$Vals[$LRR_Index]\t$Vals[$BAF_Index]\n";
		}
		
	}
	else
	{
		if(!exists($Exclude{$Vals[$Sample_Index]}))
		{
		print RESULTFILE "$Vals[$SNP_Name_Index]\t$Vals[$LRR_Index]\t$Vals[$BAF_Index]\n";
		}
	}
	$lastID=$Vals[$Sample_Index];
	close DEFFILE;
}
print "$i Samples Processed Final Report To BafLrr and $exclude excluded\n";
close REPORTFILE;
close RESULTFILE;
done

