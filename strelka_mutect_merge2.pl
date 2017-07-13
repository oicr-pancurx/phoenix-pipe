#! /usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;




### this script will take the strelka and mutect vcf files and do a merge based on common snv calls
### it will rewrite the header in a manner that maintains information from both sources
### and rebuild the info/format/calls

#my %opts=(
#	translate=>{
#					mutect=>{DP=>"DPM"},
#				}	
#
#
#);






my %opts=(
	canonical=>1,  ### only show canonical
);

GetOptions(
	\%opts,
	'strelka=s'    	=> \$opts{strelka},
	'mutect=s'      => \$opts{mutect},
);

## check for options
die "strelka vcf file not provided or not found" if(  (!defined($opts{strelka})) || ! -e $opts{strelka});
die "mutect vcf file not provided or not found"  if(  (!defined($opts{mutect} )) || ! -e $opts{mutect});

my %gtypes=(ref=>"0/0",het=>"0/1",hom=>"1/1");

### read these both in and create an output that reflects information from both sources
### add in FORMAT fields to the header as necessar
### fix genotype calls in Mutect

my %headers;my %headerline;my %calls;

### load the vcf information
my %vcf=load_vcf($opts{strelka},$opts{mutect});


#### these are headers that need to be injected
####  this should have been handled when the strelka and mutect output was filtered and annotated
our %headers2inject=(
	INFO=>{
		ANNOVAR        =>"##INFO=<ID=ANNOVAR,Number=.,Type=String,Description=\"Annotation from Annovar\">",
		ANNOVAR_EXONIC =>"<ID=ANNOVAR_EXONIC,Number=.,Type=String,Description=\"Annotation from Annovar\">",
		DBSNP_GMAF     =>"<ID=DBSNP_GMAF,Number=1,Type=String,Description=\"dbSNP, global minor allele frequency\">",
		DBSNP_ALLELE   =>"<ID=DBSNP_ALLELE,Number=1,Type=String,Description=\"Allele Identified in dbSNP137\">",
		DBSNP_STRAND   =>"<ID=DBSNP_STRAND,Number=1,Type=String,Description=\"dbSNP strand\">",
		DBSNP_G5A      =>"<ID=DBSNP_G5A,Number=1,Type=String,Description=\"Identified in dbSNP G5A subset\">",
		DBSNP_G5	   =>"<ID=DBSNP_G5,Number=1,Type=String,Description=\"Identified in dbSNP G5 subset\">",
		TRACK          =>"<ID=TRACK,Number=1,Type=String,Description=\"Additional tracks related to variant position\">"
	}
);
for my $key(keys %headers2inject){
	for my $id(keys %{$headers2inject{$key}}){
		$vcf{header}{$key}{$id}=$headers2inject{$key}{$id};
	}
}	



#print Dumper(%vcf);exit;

%vcf=merge_headers(\%vcf,\%opts);


exit;
#print Dumper($vcf{header});<STDIN>;
print Dumper($vcf{headerblock});<STDIN>;

exit;



$vcf{snv}   = merge_snv($vcf{snv});
$vcf{indel} = merge_indel($vcf{indel});
%{$vcf{calls}}=(%{$vcf{snv}},%{$vcf{indel}});










print $vcf{mergedheader};
print $vcf{colheaders}{strelka}."\n";
my @calls=sort_calls(keys %{$vcf{calls}});
for my $call(@calls){
	print "$vcf{calls}{$call}{merged}{line}\n";
}



######################################

#### sort the calls

sub sort_calls{
	my @calls=@_;
	my %callhash;
	my @sortedcalls;   # the list of calls, sorted
	map{/chr(.*):(.*):(.*:.*)/;$callhash{$1}{$2}=$_} @calls;
	for my $chr(1..21,"X","Y"){   ### keep only these chromosomes
		for my $pos(sort {$a<=>$b} keys %{$callhash{$chr}}){
			push(@sortedcalls,$callhash{$chr}{$pos})
		}
	}
	return @sortedcalls;
}



#### this will load the vcf files into a single hash
#### requires both the strelka and mutect files
#### strelka file will be loaded first...all records
#### mutect file will be loaded second...intesecting records
#### non-intersecting strelka records then removed
sub load_vcf{
	my %vcf;
	($vcf{strelka},$vcf{mutect})=@_;
	
	#my %ids;  ### hash to store ID tags found, these cannot be used more than once. If mutect uses the same tag, it will be modified with "M"
	

	for my $source(qw/strelka mutect/){  ### order is important
		(open my $VCF,"<",$vcf{$source}) || die "could not open $source vcf file $vcf{$source}";
		while(<$VCF>){
			chomp;
			if(/^##/){   	### the headers #########################
				
				my($metakey,$metaval)=/^##(.*?)=(.*)/;   ### key value pairs, separated by =
				if($metaval=~/^\<.*\>$/){   #### is the string an ID block?
					my ($id)=$metaval=~/ID=([^,<>]+)/;
					### description needs to be modified to include source
					if($metaval=~/Description/){
						$metaval=~s/(Description=\".*?)\"/$1;source=$source\"/;
					}else{    ### add description key
						$metaval=~s/\>/,Description=\"source=$source\">/;
					}
					$vcf{merge}{header}{$metakey}{$id}{$source}=$metaval;					
				}else{
					$vcf{merge}{header}{$metakey}{noid}{$source}=$metaval;
				}
			}elsif(/^#/){  	### the column headerline ########################
				$vcf{colheaders}{$source}=$_;

			}else{   		#### the data records ##########################
				my @f=split /\t/;
				my $call=join(":",@f[0,1,3,4]);
				my $class=(length($f[3])==1 && length($f[4])==1) ? "snv" : "indel";
				$vcf{merge}{$class}{$call}{$source}{line}=$_;
				$vcf{merge}{$class}{$call}{$source}{info}=$f[7];
				$vcf{merge}{$class}{$call}{$source}{format}=$f[8];
				$vcf{merge}{$class}{$call}{$source}{normal}=$source eq "strelka" ? $f[9] : $f[10];  ### mutect has tumor first, normal second
				$vcf{merge}{$class}{$call}{$source}{tumor} =$source eq "strelka" ? $f[10] : $f[9];   ### 
				$vcf{merge}{$class}{$call}{$source}{main}=join("\t",@f[0..6]);
			}	
		}
	}
	#### return the merge vcf record hash
	return %{$vcf{merge}};
}



sub merge_indel{
	my ($indel)=@_;
	
	
	my @calls=keys %{$indel};
	for my $call(@calls){
		
		$$indel{$call}{merged}{main}=$$indel{$call}{strelka}{main};
		$$indel{$call}{merged}{info}=$$indel{$call}{strelka}{info};
		
		### add genotype information
		$$indel{$call}{merged}{format}="GT:".$$indel{$call}{strelka}{format};
		
	
		
		my ($norm,$tumor)=$$indel{$call}{strelka}{info}=~/SGT=(.*?)->(.*?);/;
		die "cannot extract sgt info" unless($norm && $tumor);
		my $gt_norm	=$norm eq "ref" ? "0/0" : $norm eq "het" ? "0/1" : "1/1";
		my $gt_tumor=$tumor eq "ref" ? "0/0" : $tumor eq "het" ? "0/1" : "1/1";
		$$indel{$call}{merged}{normal}=$gt_norm.":".$$indel{$call}{strelka}{normal};
		$$indel{$call}{merged}{tumor} =$gt_tumor.":".$$indel{$call}{strelka}{tumor};
		
		
		
	
		$$indel{$call}{merged}{line}=join("\t",@{$$indel{$call}{merged}}{qw/main info format normal tumor/});
	
		
	
	}
	return $indel;	
}


sub merge_snv{

	my ($snv)=@_;
	
	my @calls=keys %$snv;
    
    #my %callhash;
    #my @sortedcalls;   # the list of calls, sorted
	#map{/chr(.*):(.*):(.*:.*)/;$callhash{$1}{$2}=$_} @calls;
	#for my $chr(1..21,"X","Y"){   ### keep only these chromosomes
 	#	for my $pos(sort {$a<=>$b} keys %{$callhash{$chr}}){
 	#		push(@sortedcalls,$callhash{$chr}{$pos})
 	#	}
 	#}	
	#for my $call(@sortedcalls){
	
	for my $call(@calls){
		
		#### MAIN : fields before INFO columns
 		if($$snv{$call}{strelka}{main} ne $$snv{$call}{mutect}{main}){
 			warn "snv call at $call differs.  setting to mutect\n";  ### due to one having an rsid, the other not
 			###print Dumper($$snv{$call});<STDIN>;
 			
 			$$snv{$call}{merged}{main}=$$snv{$call}{mutect}{main};   ### this seems more complete AT LEAST FOR THE ONE EXAMPLE (lazy...write a generalized check)
 		}else{
 			$$snv{$call}{merged}{main}=$$snv{$call}{strelka}{main};
 		}
 	

		#### INFO
		my %info;
		for my $source(qw/strelka mutect/){
			map{$info{$_}++} split /;/,$$snv{$call}{$source}{info};
		}
		$$snv{$call}{merged}{info}=join(";",sort keys %info);
	
	
	   	#### FORMAT and genotype calls 
		my %samps;
		for my $source(qw/strelka mutect/){
		 	
			my @fields=split /:/,$$snv{$call}{$source}{format};
				
			### split the sample info and store as values
			@{$samps{normal}{$source}}{@fields}=split /:/,$$snv{$call}{$source}{normal};
			@{$samps{tumor}{$source}}{@fields}=split /:/,$$snv{$call}{$source}{tumor};
			
		}
		
		### merged calls, field identifiers may change (eg mutect DP -> DPM)
		for my $samp(qw/normal tumor/){
	 		@{$samps{$samp}{merged}}{qw/AU CU GU TU DP FDP SDP SUBDP/}	=@{$samps{$samp}{strelka}}{qw/AU CU GU TU DP FDP SDP SUBDP/};
	 		@{$samps{$samp}{merged}}{qw/BQ DPM FA SS AD GT/}			=@{$samps{$samp}{mutect}}{qw/BQ DP FA SS AD GT/};
 			
 		}
	 	$samps{normal}{merged}{GT}="0/0" if($samps{normal}{merged}{GT} eq "0");
 
 		$$snv{$call}{merged}{format}=join(":",qw/GT BQ DP DPM FDP SDP SUBDP FA SS AD AU CU GU TU/);
		$$snv{$call}{merged}{normal}=join(":",@{$samps{normal}{merged}}{qw/GT BQ DP DPM FDP SDP SUBDP FA SS AD AU CU GU TU/});
		$$snv{$call}{merged}{tumor} =join(":",@{$samps{tumor}{merged}}{qw/GT BQ DP DPM FDP SDP SUBDP FA SS AD AU CU GU TU/});
 	
		$$snv{$call}{merged}{line}=join("\t",@{$$snv{$call}{merged}}{qw/main info format normal tumor/});
 	
 		
	}
	return $snv;
}

sub merge_headers{
	my ($vcf,$opts)=@_;
	my @keys=sort keys %{$$vcf{header}};
	my @orderedkeys=qw/fileformat FORMAT INFO FILTER/;
	
	for my $okey(@orderedkeys){
		@keys=grep{$_ ne $okey} @keys;  ### remove the ordered keys from the keys
	}
	my $headerblock="##fileformat=".$$vcf{header}{fileformat}{noid}{strelka}."\n";   ### assumes these values are the same
	### keys to start with
	for my $key(qw/FORMAT INFO FILTER/){
		print "$key";<STDIN>;
		#print Dumper($$vcf{header}{$key});<STDIN>;
		
		for my $id(sort keys %{$$vcf{header}{$key}}) {
			
			### there should be only a single ID value for each key, 
			### if multiple, then one will need modification
			my $sourcecount = scalar keys %{$$vcf{header}{$key}{$id}};
			if($sourcecount>1){
				print "$key $id " . Dumper($$vcf{header}{$key}{$id});<STDIN>;
			    
			}
			
			
			for my $caller(qw/strelka mutect/){
				if(my $headerline=$$vcf{header}{$key}{$id}{$caller}){
				$headerblock.="##".$key."=".$headerline."\n";
				}
			}	
		}
	}

	for my $key(@keys){
		for my $id(sort keys %{$$vcf{header}{$key}}){
			for my $caller(qw/strelka mutect/){
				$headerblock.="##".$key."=".$vcf{header}{$key}{$id}{$caller}."\n" if(defined $$vcf{header}{$key}{$id}{$caller});
			}
		}
	}
	return $headerblock;

}




				








