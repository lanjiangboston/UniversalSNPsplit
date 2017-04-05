#!/usr/bin/perl -w
# Author: Lan Jiang <lan.jiang.boston@gmail.com>;
# Created: 21 Nov 2016;
# Version: 0.01;
# Version: 0.02; try to speed up
#use warnings;
#use strict;
##### Goal: extact allele supporting reads based the SNP information  #####

##### take how much time? It mainly dependent on the row number of pileup file, usually, 212933 lines of pileup files take about 8 hours.
##### how to check the progress? the percentage of done is equal to #line in .runlog file/ #line in .pileup files.
##### what is ..errlog file? 
##### is any difference between SE and PE reads ?



#tophatbam2allelecounttable_v3.sh
#bam2allelecounttable
#pay attention to -Q 0 in 'samtools mpileup' and --readQ 0 in 'python /netscr/lanjiang/bin/pileup2ac.py'

#cd /nfs3/lanjiang/AsymmetricDHS/RNAseq2hybrid_morula_kdm6b/Tophatbam2AlleleCountTable
#sname=hybrid_morula_Kdm6bMut
#bam=/nfs3/lanjiang/AsymmetricDHS/RNAseq2hybrid_morula_kdm6b/tophat_bam/hybrid_morula_Kdm6bMut/accepted_hits.bam    #30222900 reads
#exontagSNPbed=/nfs3/lanjiang/SNP_track/mm9/AKR_PWKPhJ/exon.tagSNP.AKR_PWKPhJ.20161107.simple.bed                   #301855   SNPs
#awk '{print $1"\t"$2"\t"$2+1}' $exontagSNPbed | samtools view -b -q 10 -L - $bam > $sname.SNP.unique.bam           #generate bam file contains reads unique and cover SNP. less than 5min. got 5565261 reads
#samtools index $sname.SNP.unique.bam
#cut -f1,2 $exontagSNPbed | samtools mpileup -O -Q 0 -l - $sname.SNP.unique.bam  > $sname.pileup                    #143646 SNPs
#perl /nfs3/lanjiang/software/Extract_AlleleSupportReads_from_Sam_20161210.pl $sname $exontagSNPbed $sname.SNP.unique.bam  $sname.pileup 1>1 2>2#
#/nfs3/lanjiang/software/subread-1.5.1-source/bin/featureCounts  -O -T 1 -t exon -g gene_id -a /nfs/Shared/genome/mm9/Annotation/Genes/genes.gtf -o $sname.split.genotype1.ReadsCount.txt $sname.split.genotype1.sam
#/nfs3/lanjiang/software/subread-1.5.1-source/bin/featureCounts  -O -T 1 -t exon -g gene_id -a /nfs/Shared/genome/mm9/Annotation/Genes/genes.gtf -o $sname.split.genotype2.ReadsCount.txt $sname.split.genotype2.sam


my $sname       = shift @ARGV;  #$sname
my $SNP_file    = shift @ARGV;  #$exontagSNPbed
my $BAM_file    = shift @ARGV;  #$bam
my $PILEUP_file = shift @ARGV;  #$sname.pileup


open SNP, "<$SNP_file" or die"can't open:$!";    
#format example:
#head /nfs3/lanjiang/SNP_track/mm9/B6_JF1/exon.tagSNP.JF1.vs.B6.20161118.simple.bed 
#chr1    3205922 3205922 G       A
#chr1    3206100 3206100 A       G
my %snp1;
my %snp2;
while(<SNP>){
	chomp;
        my @line = split /\t/, $_;
	    my $loc1 = join("_",$line[0],$line[1],$line[3]); #chr1    3205922 3205922 G       A
	    my $loc2 = join("_",$line[0],$line[1],$line[4]); 
	    $snp1{$loc1} =1;
	    $snp2{$loc2} =1; 

	}
close SNP;


open FILEOUT1,">$sname.split.genotype1.sam" or die"can't open:$!";
open FILEOUT2,">$sname.split.genotype2.sam" or die"can't open:$!";
open FILEOUT3,">$sname.split.stat.txt" or die"can't open:$!";
open FILEOUT4,">$sname.runlog" or die"can't open:$!";
open FILEOUT5,">$sname.errlog" or die"can't open:$!";

my %genotype1_SNP;
my %genotype2_SNP;
my %sam_file;
my %whichSNP_covered_by_reads;
open PILEUP_FILE, "<$PILEUP_file" or die"can't open:$!";  
#format example:
#chr1    4481127 N       3       gga     ^[W     16,6,2
#chr1    4481157 N       5       ggaaa   _TXZ^   46,36,32,19,10
#chr1    4481582 N       7       gaaaaaa VbZbHZ_ 47,31,31,24,23,23,23

while(<PILEUP_FILE>){
    chomp;
    my @pileup_line = split /\t/, $_;                                                #chr1    4481127 N       3       gga     ^[W     16,6,2
    my @Base_postion = split //, $pileup_line[4];
    my @SNP_postion = split /\,/, $pileup_line[6];
    my $SNP_chr = $pileup_line[0];
    my $SNP_pos = $pileup_line[1];
    #my @sam4oneSNP=`samtools view $BAM_file $SNP_chr:$SNP_pos-$SNP_pos`;
    my @sam4oneSNP=`samtools view $BAM_file $SNP_chr:$SNP_pos-$SNP_pos | grep -v "YT:Z:UP"`; #update 20170119, filter reads mapped un concordantly. Becuase when pileup file generated, those reads are filtered. So if when do samtools view, not filter those reads, will cause $#sam2oneSNP+1 != $pileup_line[3] 

    if($#sam4oneSNP+1 == $pileup_line[3])
    {
        print FILEOUT4 $#sam4oneSNP+1, "\n" ;
        for(my $i=0; $i<=$#sam4oneSNP; $i++)
        {
            my $sam4oneread = $sam4oneSNP[$i];                                               
            my @sam_line = split /\t/, $sam4oneread;                                     #1458_1065_1988_F3       0       chr1    34332083        50      50M     *       0       0       AACTTCCACAACTCCCTCCAGGACTTCATCAACTGGCTTACCCAGGCTGA      Aaabccbabba`_acbb`\^_^_[[_VPWWY^[Z\\[^QJY]ZSVNK^bA     XA:i:0  MD:Z:50 NM:i:0  CM:i:0  XS:A:+  NH:i:1
            my $read_name = $sam_line[0];
            my $read_seq = $sam_line[9];                                                 #assume in bam/sam files, the sequence already transform to be in waston strand.
            my $base_support_by_read = substr($read_seq, $SNP_postion[$i]-1, 1); 
            #print $read_name, "\t", $base_support_by_read, "\n";
            if($Base_postion[$i]=~/[ATCGatcg]/)  # to avoid spliced reads covered SNP, 
            {
                $tag=join("_", $SNP_chr, $SNP_pos, $base_support_by_read);                   #
                $genotype1_SNP{$read_name}+=0;
                $genotype2_SNP{$read_name}+=0;
                $genotype1_SNP{$read_name}++ if(exists $snp1{$tag});
                $genotype2_SNP{$read_name}++ if(exists $snp2{$tag});
                $sam_file{$read_name} = $sam4oneread;
            }
            #$whichSNP_covered_by_reads{$read_name}.="$tag\;" if(exists $snp1{$tag});
            #$whichSNP_covered_by_reads{$read_name}.="$tag\;" if(exists $snp2{$tag});
        }
    }
    else
    {
        print FILEOUT5 $_, "\n" ;
    }
}

close PILEUP_FILE;



foreach my $ID(keys %sam_file)
{
    print FILEOUT1 $sam_file{$ID} if( $genotype1_SNP{$ID} > $genotype2_SNP{$ID} );
    print FILEOUT2 $sam_file{$ID} if( $genotype1_SNP{$ID} < $genotype2_SNP{$ID} );
    print FILEOUT3 "$ID\tgenotype1_SNP\t$genotype1_SNP{$ID}\tgenotype2_SNP\t$genotype2_SNP{$ID}\n";
    #print FILEOUT3 "$ID\tgenotype1_SNP\t$genotype1_SNP{$ID}\tgenotype2_SNP\t$genotype2_SNP{$ID}\t$whichSNP_covered_by_reads{$read_name}\n";
}

close FILEOUT1;
close FILEOUT2;
close FILEOUT3;
close FILEOUT4;
close FILEOUT5;



