#!/usr/bin/perl

use strict;
use warnings;
use Tabix;		# requires module load tabix!
use Data::Dumper;
use List::Util qw(sum);

#input files crest .predSV.txt file and collapsed bam file
my $crest = $ARGV[0];
my $bam = $ARGV[1];
#my $normal_bam = $ARGV[2];

my $scriptPath = "/.mounts/labs/PCSI/production/phoenix-pipe";		# path to findSupportingPairs.pl

my $smallDelSize = 1000;			# if type is deletion flag if break point are less than $dist bases apart and fewer than 2 clips on either side
my $smallDelClipMin = 3;

my $itxdist = 300;		# if type is ITX flag if break point are less than $itxdist bases apart and zero clips on either side
my $chrdist = 125;		# if left and right break points are on same chromosome and type is not a deletion, flag if break points are less than $chrdist bases apart

my $num_dis_pairs = 3;	# number of discordnant pairs required
#my $num_normal_dis_pairs = 0;	# number of normal discordnant pairs tolerated



#bed file of hg19 refseq genes/isoforms

my $bed = "/.mounts/labs/prod/archive/projects/PCSI/filter_crest/hg19_RefSeq_sorted.bed.gz" ; 

 
my %filtered;
my %final;

#read crest file

open (BREAK, "<$crest") or die "No Crest File Found\n";

while (<BREAK>){
	chomp $_;                  
	next if (/^#/);
	next if (/left/);
	next if (/chrUn/);

	my $flag = 0;
	my $cl_flag = 0;
	my ($lchr,$lbr,$lst,$lcl,$rchr,$rbr,$rst,$rcl,$type,$rest)= split(/\t/,$_);


# 4 crest filtering steps:


#1) flag if there are zero soft clips at left or right break point 

#NOTE: currently commented out

#	if (($lcl ==0) || ($rcl ==0)){
#		$cl_flag = 3;
#		}

#### ITX clip 0 rule big enough
       



#2) if type is deletion flag if break point are less than $dist bases apart and fewer than 2 clips on either side

	if (($type eq "DEL") and (($rbr - $lbr) < $smallDelSize)){
		if (($lcl < $smallDelClipMin) or ($rcl < $smallDelClipMin)){
			$cl_flag = 2;
		}
	}



#3) if type is ITX flag if break point are less than $itxdist bases apart and zero clips on either side

	if (($type eq "ITX") && (abs($rbr- $lbr)<$itxdist)){
		if (($lcl ==0) || ($rcl ==0)){
			$cl_flag = 2;
		}
	}
 


#4) if left and right break points are on same chromosome and type is not a deletion, flag if break points are less than $chrdist bases apart

	if ($lchr eq $rchr){
		if ((abs($rbr - $lbr)<$chrdist) && ($type ne "DEL")){
			$flag = 1;
		}
	}

# large deletions should require at least 3 discorant pairs

	


	my $leftbreak = join (":",$lchr,$lbr);
	my $fullleft = join ("",$leftbreak,$lst);

	my $rightbreak = join (":",$rchr,$rbr);
	my $fullright = join ("",$rightbreak,$rst);

#print "$fullleft $fullright\n";
	my $discordant_count = 0;
#	my $normal_discordant_count = 0;



#run findSypporitngPairs script to identify discordant pair read support for break points




	$discordant_count = `$scriptPath/findSupportingPairs.pl $bam $fullleft $fullright |wc -l`; 
#	$normal_discordant_count = `$scriptPath/findSupportingPairs.pl $normal_bam $fullleft $fullright |wc -l`; 

	chomp $discordant_count;
#	chomp $normal_discordant_count;



#final filtering: based on number of discordant read support and no flags:

	next if ((($discordant_count <$num_dis_pairs) or ($flag !=0) or ($cl_flag !=0)) and ($type ne "DEL"));
#	next if ((($discordant_count <$num_dis_pairs) or ($flag !=0) or ($cl_flag !=0) or ($normal_discordant_count > $num_normal_dis_pairs)) and ($type ne "DEL"));
	next if ( ($type eq "DEL") && ($cl_flag !=0));



#	print "$fullleft\t$fullright\n";

####	tabix annotation of genes
	#add repeat track annotation




# 
#### remove reduandant calls ie reversal of left and right breakpoints


	if ((exists $filtered{$leftbreak} and $filtered{$leftbreak} eq $rightbreak) || (exists $filtered{$rightbreak} and $filtered{$rightbreak} eq $leftbreak)) {
#	if (exists $filtered{$rightbreak}) {
#		print "Duplicate Found: $leftbreak\t$rightbreak\n";
		if (exists $filtered{$leftbreak}){	
#			print "found a duplicate left: $filtered{$leftbreak} and $leftbreak\n";	
			}
		if (exists $filtered{$rightbreak}){
#			print "found a duplicate right: $filtered{$rightbreak} and $rightbreak\n";
                         }

		}
	else{
		$filtered{$leftbreak}=$rightbreak;
		$final{$leftbreak}{'lchr'}=$lchr;
		$final{$leftbreak}{'lbr'}=$lbr;
		$final{$leftbreak}{'lst'}=$lst;
		$final{$leftbreak}{'lcl'}=$lcl;
		$final{$leftbreak}{'rchr'}=$rchr;
		$final{$leftbreak}{'rbr'}=$rbr;
		$final{$leftbreak}{'rst'}=$rst;
		$final{$leftbreak}{'rcl'}=$rcl;
		$final{$leftbreak}{'type'}=$type;
		$final{$leftbreak}{'left_break_iso'}="";	
		$final{$leftbreak}{'left_break_gene'}="";
		$final{$leftbreak}{'right_break_iso'}="";
		$final{$leftbreak}{'right_break_gene'}="";

		$final{$leftbreak}{'count'}=$discordant_count;

}
#}

}

#my $tabix = Tabix->new(-data => $bed, -index => "$bed.tbi");
#my $iter_l;
#my $iter_r;
#my $iter_up;
#my $iter_down;


foreach my $ele (keys %final){

	my $left_chr = $final{$ele}{'lchr'};
	my $left_start = $final{$ele}{'lbr'} -1;
	my $left_end = $final{$ele}{'lbr'};
	
	my $up_start= $left_start-1001;
	my $up_end= $left_start-1;


	my $right_chr = $final{$ele}{'rchr'};
	my $right_start = $final{$ele}{'rbr'} -1;
	my $right_end = $final{$ele}{'rbr'};

	my $down_start= $right_end +1;
	my $down_end= $right_end+1001;


# gene anootation for plus or minus 1000 bp of break


#	$iter_l = $tabix->query($left_chr,$left_start,$left_end);
#	$iter_r = $tabix->query($right_chr,$right_start,$right_end);

#	$iter_up = $tabix->query($left_chr,$up_start,$up_end);
#	$iter_down= $tabix->query($right_chr,$down_start,$down_end);

#	                if (defined $iter_l->{"_"})               # happens if the contig isn't in the bed file
 #               {   
  #                      if (my $bedReturn_l = $tabix->read($iter_l))
   #                     {

#			my @return_l = split (/\t/,$bedReturn_l);
#			my $gene_l= $return_l[3];
#			my $iso_l= $return_l[4];
		 
   
 #                        $final{$ele}{'left_break_gene'}=$gene_l; 
  #                       $final{$ele}{'left_break_iso'}=$iso_l; 
   # }
#}


#	                if (defined $iter_r->{"_"})               # happens if the contig isn't in the bed file
 #               {   
#                        if (my $bedReturn_r = $tabix->read($iter_r))
#                        {  

#			my @return_r = split (/\t/,$bedReturn_r);
#			my $gene_r= $return_r[3];
#			my $iso_r= $return_r[4];
 
#                         $final{$ele}{'right_break_gene'}=$gene_r; 
#                         $final{$ele}{'right_break_iso'}=$iso_r; 
 #   }
#}

# if (defined $iter_up->{"_"})               # happens if the contig isn't in the bed file
 #                {       
 #                        if (my $bedReturn_up = $tabix->read($iter_up))
  #                       {
 
#                         my @return_up = split (/\t/,$bedReturn_up);
 #                        my $gene_up= $return_up[3];
  #                       my $iso_up= $return_up[4];
 
#                          $final{$ele}{'upsteam_gene'}=$gene_up;
 #                         $final{$ele}{'upstream_iso'}=$iso_up;
  #   }
# }


# if (defined $iter_down->{"_"})               # happens if the contig isn't in the bed file
 #                {       
  #                       if (my $bedReturn_down = $tabix->read($iter_down))
   #                      {
 
    #                     my @return_down = split (/\t/,$bedReturn_down);
     #                    my $gene_down= $return_down[3];
      #                   my $iso_down= $return_down[4];
 
#                          $final{$ele}{'downsteam_gene'}=$gene_down;
  #                        $final{$ele}{'downstream_iso'}=$iso_down;
 #    }
# }



}

#print Dumper %final;

foreach my $bp (keys %final){

#	print "$final{$bp}{'lchr'}\t$final{$bp}{'lbr'}\t$final{$bp}{'lst'}\t$final{$bp}{'lcl'}\t$final{$bp}{'rchr'}\t$final{$bp}{'rbr'}\t$final{$bp}{'rst'}\t$final{$bp}{'rcl'}\t$final{$bp}{'type'}\t$final{$bp}{'left_break_gene'}\t$final{$bp}{'right_break_gene'}\t$final{$bp}{'upstream_gene'}\t$final{$bp}{'downstream_gene'}\n";
	print "$final{$bp}{'lchr'}\t$final{$bp}{'lbr'}\t$final{$bp}{'lst'}\t$final{$bp}{'lcl'}\t$final{$bp}{'rchr'}\t$final{$bp}{'rbr'}\t$final{$bp}{'rst'}\t$final{$bp}{'rcl'}\t$final{$bp}{'type'}\t$final{$bp}{'count'}\n";




}
