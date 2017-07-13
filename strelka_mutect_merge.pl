#! /usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;


### this script will take the strelka and mutect vcf files and do a merge based on common snv calls
### it will rewrite the header in a manner that maintains information from both sources
### and rebuild the info/format/calls

my %opts=(
	translate=>{
					mutect=>{DP=>"DPM"},
				}	


);

our %newheaders=(
	INFO=>"##INFO=<ID=ANNOVAR,Number=.,Type=String,Description=\"Annotation from Annovar\">
##INFO=<ID=ANNOVAR_EXONIC,Number=.,Type=String,Description=\"Annotation from Annovar\">
##INFO=<ID=DBSNP_GMAF,Number=1,Type=String,Description=\"dbSNP, global minor allele frequency\">
##INFO=<ID=DBSNP_ALLELE,Number=1,Type=String,Description=\"Allele Identified in dbSNP137\">
##INFO=<ID=DBSNP_STRAND,Number=1,Type=String,Description=\"dbSNP strand\">
##INFO=<ID=DBSNP_G5A,Number=1,Type=String,Description=\"Identified in dbSNP G5A subset\">
##INFO=<ID=DBSNP_G5,Number=1,Type=String,Description=\"Identified in dbSNP G5 subset\">
##INFO=<ID=TRACK,Number=1,Type=String,Description=\"Additional tracks related to variant position\">
"
);





GetOptions(
	\%opts,
	'strelka=s'    	=> \$opts{strelka},
	'mutect=s'      => \$opts{mutect}
);

my %gtypes=(ref=>"0/0",het=>"0/1",hom=>"1/1");

### read these both in and create an output that reflects information from both sources
### add in FORMAT fields to the header as necessar
### fix genotype calls in Mutect

my %headers;my %headerline;my %calls;

### load the vcf information
my %vcf;
load_vcf(\%opts,\%vcf);




$vcf{snv}   =merge_snv($vcf{snv});
$vcf{indel} =merge_indel($vcf{indel});


$vcf{mergedheader}=merge_headers($vcf{header});
%{$vcf{calls}}=(%{$vcf{snv}},%{$vcf{indel}});







print $vcf{mergedheader};
print $vcf{colheaders}{strelka}."\n";
my @calls=sort_calls(keys %{$vcf{calls}});
for my $call(@calls){
	print "$vcf{calls}{$call}{merged}{line}\n";
}





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


sub load_vcf{
	my ($opts,$href)=@_;
	
	for my $source(qw/strelka mutect/){  ### order is important
		my $file=$$opts{$source};
		(open my $VCF,"<",$file) || die "could not open $source vcf file $file";
		while(<$VCF>){
			chomp;
			if(/^##/){   ### the headers
				
				my($heading,$value)=/^##(.*?)=(.*)/;
				#### check for an id field
				my $id=$value=~/<ID=(.*?),/ ? $1 : "-";    ### store on key and id, noid if field does not inlude one
				if(my $newid=$opts{translate}{$source}{$id}){
					$value=~s/ID=${id}/ID=${newid}/;
					$id=$newid
				}
				
				#### add the source to the value description tag for bracketed values
				#print "$value\n";
				if(/\=\<.*\>/){
					
					### is there already a description
					if(my ($desc)=$value=~/Description=\"(.*)\"/){
						$value=~s/$desc/${desc};source=${source}/;
					}else{
						### add a description
						$value=~s/>/,Description=\"source=${source}\">/;
					}
					
				
				}	
				#print "$value\n";<STDIN>;
				$$href{header}{$heading}{$id}{$source}=$value;


			}elsif(/^#/){  ### the column headerline
				$vcf{colheaders}{$source}=$_;

			}else{   #### the data
				my @f=split /\t/;
				my $call=join(":",@f[0,1,3,4]);
				if($source eq "strelka" || ($source eq "mutect" && $$href{snv}{$call})){  #only keep muect calls already in strelka
					#print "$_";<STDIN>;
					my $class=(length($f[3])==1 && length($f[4])==1) ? "snv" : "indel";
				
					$$href{$class}{$call}{$source}{line}=$_;
					$$href{$class}{$call}{$source}{info}=$f[7];
					$$href{$class}{$call}{$source}{format}=$f[8];
					$$href{$class}{$call}{$source}{normal}=$source eq "strelka" ? $f[9] : $f[10];  ### mutect has tumor first, normal second
					$$href{$class}{$call}{$source}{tumor} =$source eq "strelka" ? $f[10] : $f[9];   ### 
					$$href{$class}{$call}{$source}{main}=join("\t",@f[0..6]);
				}
			}	
		}
	}
	
	### remove any snv calls not in both strelka and mutect
	for my $call(sort keys %{$$href{snv}}){
		delete($$href{snv}{$call}) unless($$href{snv}{$call}{strelka} && $$href{snv}{$call}{mutect});
	}
	
	
		
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
	my ($headers)=@_;
	
	my @keys=sort keys %$headers;
	my @orderedkeys=qw/fileformat FORMAT INFO FILTER/;
	for my $okey(@orderedkeys){
		@keys=grep{$_ ne $okey} @keys;  ### remove the ordered keys from the keys
	}

	my $fullheader="##fileformat=".$$headers{fileformat}{"-"}{strelka}."\n";
	### keys to start with
	for my $key(qw/FORMAT INFO FILTER/){
		for my $id(sort keys %{$$headers{$key}}) {
			for my $caller(qw/strelka mutect/){
				if(my $headerline=$$headers{$key}{$id}{$caller}){
				$fullheader.="##".$key."=".$headerline."\n";
				}
			}	
		}
		if(my $new=$newheaders{$key}){
			$fullheader.=$new;
		}
		
	}

	for my $key(@keys){
		for my $id(sort keys %{$$headers{$key}}){
			for my $caller(qw/strelka mutect/){
				$fullheader.="##".$key."=".$$headers{$key}{$id}{$caller}."\n" if(defined $$headers{$key}{$id}{$caller});
			}
		}
	}
	return $fullheader;

}




				








