#!/usr/bin/perl

use strict;
use warnings;

use Tabix;
use Bio::DB::Fasta;

my $pepLength = 9;
my $peptide;

my %fastaHandles = ("path" => "/oicr/data/reference/genomes/homo_sapiens_mc/UCSC/hg19_random/Genomic/references/fasta/");
my $refSeqFile = "/oicr/data/reference/genomes/homo_sapiens_mc/refSeq/refSeq.tsv";
my $exonFile = "/oicr/data/reference/genomes/homo_sapiens_mc/refSeq/refSeq_exons.sorted.bed.gz";

my $exonTabix = Tabix->new(-data => $exonFile, -index => "$exonFile.tbi");

my %aaCode = 
(
	"TTT" => "F",
	"TTC" => "F",
	"TTA" => "L",
	"TTG" => "L",
	"TCT" => "S",
	"TCC" => "S",
	"TCA" => "S",
	"TCG" => "S",
	"TAT" => "Y",
	"TAC" => "Y",
	"TAA" => "X",
	"TAG" => "X",
	"TGT" => "C",
	"TGC" => "C",
	"TGA" => "X",
	"TGG" => "W",
	"CTT" => "L",
	"CTC" => "L",
	"CTA" => "L",
	"CTG" => "L",
	"CCT" => "P",
	"CCC" => "P",
	"CCA" => "P",
	"CCG" => "P",
	"CAT" => "H",
	"CAC" => "H",
	"CAA" => "Q",
	"CAG" => "Q",
	"CGT" => "R",
	"CGC" => "R",
	"CGA" => "R",
	"CGG" => "R",
	"ATT" => "I",
	"ATC" => "I",
	"ATA" => "I",
	"ATG" => "M",
	"ACT" => "T",
	"ACC" => "T",
	"ACA" => "T",
	"ACG" => "T",
	"AAT" => "N",
	"AAC" => "N",
	"AAA" => "K",
	"AAG" => "K",
	"AGT" => "S",
	"AGC" => "S",
	"AGA" => "R",
	"AGG" => "R",
	"GTT" => "V",
	"GTC" => "V",
	"GTA" => "V",
	"GTG" => "V",
	"GCT" => "A",
	"GCC" => "A",
	"GCA" => "A",
	"GCG" => "A",
	"GAT" => "D",
	"GAC" => "D",
	"GAA" => "E",
	"GAG" => "E",
	"GGT" => "G",
	"GGC" => "G",
	"GGA" => "G",
	"GGG" => "G",
	"NNN" => "!",
);

my %refSeq;

my $l;
my @f;

my $name;
my $chr;
my $strand;
my $starts;
my $ends;
my $codingStart;
my $codingEnd;

my $sample;

my $pos;
my $ref;
my $alt;
my $names;
my ($nuc,$aa,$nucPos,$aaPos,$newNuc,$newAA);

my $annovarString;

# read refseq exon start and ends into hash of arrays
open (FILE, $refSeqFile) or die "Couldn't open $refSeqFile\n";

while ($l = <FILE>)
{
	unless ($l =~ /^#/)
	{
		chomp $l;
		@f = split(/\t/, $l);

		$name = $f[1];
		$chr = $f[2];
		$strand = $f[3];
		$codingStart = $f[6];
		$codingEnd = $f[7];
		$starts = $f[9];
		$ends = $f[10];

		$refSeq{$name}{chr} = $chr;
		$refSeq{$name}{strand} = $strand;
		$refSeq{$name}{cds_start} = $codingStart + 1;
		$refSeq{$name}{cds_end} = $codingEnd;
		@{ $refSeq{$name}{starts} } = split(/,/, $starts);
		@{ $refSeq{$name}{ends} } = split(/,/, $ends);
	}
}
close FILE;

while ($l = <STDIN>)
{
	if ($l =~ /^##DCC=<analyzed_sample_id="(.*)">/)
	{
		$sample = $1;
	}
	elsif ($l =~ /^#/)
	{
	}
	else
	{
		chomp $l;
		@f = split(/\t/, $l);

		$chr = $f[0];
		$pos = $f[1];

		for my $altGolden (split(/,/, $f[4]))
		{
			$names = lookUpNM($chr, $pos, $exonTabix);

			for my $name (split(/,/, $names))
			{
				if (($pos >= $refSeq{$name}{cds_start}) and ($pos <= $refSeq{$name}{cds_end}))
				{
					$ref = $f[3];
					$annovarString = "";
					if ($refSeq{$name}{strand} eq "-")
					{
						$ref = reverseCompliment($ref);
						$alt = reverseCompliment($altGolden);
					}
					else
					{
						$alt = $altGolden;
					}
	
					if ($f[7] =~ /$name:exon.*?:(c.*?:p.*?),/)
					{
						$annovarString = $1;
					}
	
	#				print "Working on $chr:$pos $ref>$alt in $name\n";
					# get sequence
					$nuc = getNucForGene(\%{ $refSeq{$name} }, \%fastaHandles);
					$aa = nucToAA($nuc, \%aaCode);
	#				print $nuc . "\n";
	
					# get position in sequence
					$nucPos = getNucPosInGene(\%{ $refSeq{$name} }, $pos);
	#				print $nucPos . "\n";
					$aaPos = int(($nucPos + 2) / 3);
	
	
					# modify sequence
					$newNuc = editNucString($nuc, $nucPos, $ref, $alt);
	
					# get AA
					$newAA = nucToAA($newNuc, \%aaCode);
	
	#				print $nuc . "\n\n";
	#				print $newNuc . "\n\n";
	
	#				print $aa . "\n\n";
	#				print $newAA . "\n\n";
	
					if ($newAA ne $aa)
					{
						unless ($newAA =~ /!/)
						{
	#						print "\n   c.$ref$nucPos$alt:p." . substr($aa, $aaPos-1, 1) . "$aaPos" . substr($newAA,$aaPos-1,1) . "\n";
	#						print "A: " . $annovarString . "\n";
	
							for (my $i = 1; $i <= $pepLength; $i++)
							{
								$peptide = substr($newAA, $aaPos-$i, $pepLength);
								if (length($peptide) == $pepLength)
								{
									print "$chr\t$pos\t$ref\t$alt\t$peptide\n";
								}
							}
						}
					}
	
					# produce and output 17-mer or whatever
				}
			}
		}

	}
}



sub lookUpNM
{
	my $chr = shift;
	my $pos = shift;
	my $tabix = shift;

	my $iter = $tabix->query($chr,$pos,$pos + 1);
	my $l;
	my $NMlist = "";

    if (defined $iter->{"_"})               # happens if the contig isn't in the bed file
    {
    	while ($l = $tabix->read($iter))
        {
			chomp $l;
			$l =~ s/.*,//;
			$NMlist .= "$l,";
		}
	}

	$NMlist =~ s/,$//;
	return $NMlist;

}


sub nucToAA
{
	my $nuc = $_[0];
	my $code = $_[1];

	my $codon;
	my $aa = "";

	while ($nuc =~ s/^(...)//)
	{
		$codon = $1;

		if (exists $code->{$codon})
		{
			$aa .= $code->{$codon};
		}
		else
		{
			$aa .= "!";
		}
	}

	if (length($nuc) > 0)
	{
		warn "Had $nuc left over after conversion to AA\n";
	}

	return $aa;
}


sub getNucForGene
{
	my $refSeq = $_[0];
	my $fastaRef = $_[1];
	my $nuc = "";

	my $start = $refSeq->{cds_start};
	my $end = $refSeq->{cds_end};

	my $exonStart;
	my $exonEnd;

	for (my $i = 0; $i < scalar @{ $refSeq->{starts} }; $i++)
	{
		$exonStart = $refSeq->{starts}[$i] + 1;
		$exonEnd = $refSeq->{ends}[$i];

		if (($exonEnd < $start) or ($exonStart > $end))
		{
			# not coding
		}
		else
		{
			if (($exonStart < $start) and ($exonEnd > $end))
			{
				$exonStart = $start;
				$exonEnd = $end;
			}
			if ($exonStart < $start)
			{
				$exonStart = $start;
			}
			if ($exonEnd > $end)
			{
				$exonEnd = $end;
			}

			$nuc .= getSequence($refSeq->{chr}, $exonStart, $exonEnd, $fastaRef);
		}
	}

	if ($refSeq->{strand} eq "-")
	{
		return reverseCompliment($nuc);
	}
	else
	{
		return $nuc;
	}

}


sub getNucPosInGene
{
	my $refSeq = $_[0];
	my $pos = $_[1];

	my $runningDist = 0;
	my $newPos;

	my $start = $refSeq->{cds_start} - 1;
	my $end = $refSeq->{cds_end};

	my ($exonStart,$exonEnd);

	for (my $i = 0; $i < scalar @{ $refSeq->{starts} }; $i++)
	{
		$exonStart = $refSeq->{starts}[$i];
		$exonEnd = $refSeq->{ends}[$i];

		if (($exonEnd < $start) or ($exonStart > $end))
		{
			# not coding
		}
		else
		{
			if ($exonStart < $start)
			{
				$exonStart = $start;
			}
			if ($exonEnd > $end)
			{
				$exonEnd = $end;
			}

			if (($pos >= $exonStart) and ($pos <= $exonEnd))
			{
				$newPos = ($pos - $exonStart) + $runningDist;
			}
			$runningDist += $exonEnd - $exonStart;
		}
	}

	if ($refSeq->{strand} eq "-")
	{
		$newPos = $runningDist - $newPos + 1;
	}

	return $newPos;
}



sub editNucString
{
	my $nuc = $_[0];
	my $pos = $_[1];
	my $ref = $_[2];
	my $alt = $_[3];

	my $left = substr($nuc,0,$pos - 1);
	my $right = substr($nuc,$pos + length($ref) - 1);

	return $left . $alt . $right;
}

sub getSequence
{
	my $chr = $_[0];
	my $start = $_[1];
	my $end = $_[2];
	my $fastaHandles = $_[3];

	unless (exists ($fastaHandles->{$chr})) # create a handle for the chromosome fasta, if it doesn't exist
	{
		$fastaHandles->{$chr} = Bio::DB::Fasta->new("$fastaHandles{path}/$chr.fa");
	}

	my $seq = uc($fastaHandles{$chr}->seq($chr, $start, $end));
	return $seq;
}


sub reverseCompliment
{
	my $nuc = $_[0];

	my $rc = reverse($nuc);

	$rc =~ s/A/1/g;
	$rc =~ s/C/2/g;
	$rc =~ s/G/3/g;
	$rc =~ s/T/4/g;

	$rc =~ s/1/T/g;
	$rc =~ s/2/G/g;
	$rc =~ s/3/C/g;
	$rc =~ s/4/A/g;

	return $rc;
}


