#!/usr/bin/perl

use strict;
use warnings;
use Bio::DB::Fasta;

my $vcfFile = $ARGV[0];
my $outPath = $ARGV[1];
my $sample = $ARGV[2];

my %fastaHandles = ("path" => "/oicr/data/reference/genomes/homo_sapiens_mc/UCSC/hg19_random/Genomic/references/fasta/");
my %referenceHash = ("chunkSize" => 10000);

my %comp = (
    "ac" => "tg",
    "ag" => "tc",
    "at" => "ta",
    "ga" => "ct",
    "gc" => "cg",
    "gt" => "ca",
);

my %contextComp = (
    "aa" => "tt",
    "ac" => "gt",
    "ag" => "ct",
    "at" => "at",
    "ca" => "tg",
    "cc" => "gg",
    "cg" => "cg",
    "ct" => "ag",
    "ga" => "tc",
    "gc" => "gc",
    "gg" => "cc",
    "gt" => "ac",
    "ta" => "ta",
    "tc" => "ga",
    "tg" => "ca",
    "tt" => "aa",
);

my $change;
my $context;

my %contextCount;
for $change (qw/ca cg ct ta tc tg/)
{
    for $context (qw/aa ac ag at ca cc cg ct ga gc gg gt ta tc tg tt/)
    {
        $contextCount{"${change}_$context"} = 0;
    }
}


my ($l,$alt,$chr,$pos,$id,$ref,$alts,$qual,$filter,$info,$format,$nGT,$tGT);

# collect trinuc contexts for snvs

open (FILE, $vcfFile) or die "Couldn't open $vcfFile\n";
while ($l = <FILE>)
{
	unless ($l =~ /^#/)
	{
        ($chr,$pos,$id,$ref,$alts,$qual,$filter,$info,$format,$nGT,$tGT) = split(/\t/, $l);

        for $alt (split(/,/, $alts))
        {
			if (length($ref) == length($alt))
			{
				$change = lc("$ref$alt");
				$context = lc(getBase($chr, $pos - 1, \%referenceHash, \%fastaHandles) . getBase($chr, $pos + 1, \%referenceHash, \%fastaHandles));

				if (exists $comp{$change})
				{
					$change = $comp{$change};
					unless ($context =~ /n/)
					{
						$context = $contextComp{$context};
					}
				}
				unless ($context =~ /n/)
				{
					$contextCount{"${change}_$context"}++;
				}

				freeMemChunks($chr, $pos, \%referenceHash);
			}
		}

	}
}
close FILE;


# write trinuc and commands to R file
my @order = qw/A[C>A]A A[C>A]C A[C>A]G A[C>A]T A[C>G]A A[C>G]C A[C>G]G A[C>G]T A[C>T]A A[C>T]C A[C>T]G A[C>T]T A[T>A]A A[T>A]C A[T>A]G A[T>A]T A[T>C]A A[T>C]C A[T>C]G A[T>C]T A[T>G]A A[T>G]C A[T>G]G A[T>G]T C[C>A]A C[C>A]C C[C>A]G C[C>A]T C[C>G]A C[C>G]C C[C>G]G C[C>G]T C[C>T]A C[C>T]C C[C>T]G C[C>T]T C[T>A]A C[T>A]C C[T>A]G C[T>A]T C[T>C]A C[T>C]C C[T>C]G C[T>C]T C[T>G]A C[T>G]C C[T>G]G C[T>G]T G[C>A]A G[C>A]C G[C>A]G G[C>A]T G[C>G]A G[C>G]C G[C>G]G G[C>G]T G[C>T]A G[C>T]C G[C>T]G G[C>T]T G[T>A]A G[T>A]C G[T>A]G G[T>A]T G[T>C]A G[T>C]C G[T>C]G G[T>C]T G[T>G]A G[T>G]C G[T>G]G G[T>G]T T[C>A]A T[C>A]C T[C>A]G T[C>A]T T[C>G]A T[C>G]C T[C>G]G T[C>G]T T[C>T]A T[C>T]C T[C>T]G T[C>T]T T[T>A]A T[T>A]C T[T>A]G T[T>A]T T[T>C]A T[T>C]C T[T>C]G T[T>C]T T[T>G]A T[T>G]C T[T>G]G T[T>G]T/;

my $rowNames = "\"A[C>A]A\",\"A[C>A]C\",\"A[C>A]G\",\"A[C>A]T\",\"A[C>G]A\",\"A[C>G]C\",\"A[C>G]G\",\"A[C>G]T\",\"A[C>T]A\",\"A[C>T]C\",\"A[C>T]G\",\"A[C>T]T\",\"A[T>A]A\",\"A[T>A]C\",\"A[T>A]G\",\"A[T>A]T\",\"A[T>C]A\",\"A[T>C]C\",\"A[T>C]G\",\"A[T>C]T\",\"A[T>G]A\",\"A[T>G]C\",\"A[T>G]G\",\"A[T>G]T\",\"C[C>A]A\",\"C[C>A]C\",\"C[C>A]G\",\"C[C>A]T\",\"C[C>G]A\",\"C[C>G]C\",\"C[C>G]G\",\"C[C>G]T\",\"C[C>T]A\",\"C[C>T]C\",\"C[C>T]G\",\"C[C>T]T\",\"C[T>A]A\",\"C[T>A]C\",\"C[T>A]G\",\"C[T>A]T\",\"C[T>C]A\",\"C[T>C]C\",\"C[T>C]G\",\"C[T>C]T\",\"C[T>G]A\",\"C[T>G]C\",\"C[T>G]G\",\"C[T>G]T\",\"G[C>A]A\",\"G[C>A]C\",\"G[C>A]G\",\"G[C>A]T\",\"G[C>G]A\",\"G[C>G]C\",\"G[C>G]G\",\"G[C>G]T\",\"G[C>T]A\",\"G[C>T]C\",\"G[C>T]G\",\"G[C>T]T\",\"G[T>A]A\",\"G[T>A]C\",\"G[T>A]G\",\"G[T>A]T\",\"G[T>C]A\",\"G[T>C]C\",\"G[T>C]G\",\"G[T>C]T\",\"G[T>G]A\",\"G[T>G]C\",\"G[T>G]G\",\"G[T>G]T\",\"T[C>A]A\",\"T[C>A]C\",\"T[C>A]G\",\"T[C>A]T\",\"T[C>G]A\",\"T[C>G]C\",\"T[C>G]G\",\"T[C>G]T\",\"T[C>T]A\",\"T[C>T]C\",\"T[C>T]G\",\"T[C>T]T\",\"T[T>A]A\",\"T[T>A]C\",\"T[T>A]G\",\"T[T>A]T\",\"T[T>C]A\",\"T[T>C]C\",\"T[T>C]G\",\"T[T>C]T\",\"T[T>G]A\",\"T[T>G]C\",\"T[T>G]G\",\"T[T>G]T\"";

my $values = "";

for my $c (@order)
{
	$c =~ /(.)\[(.)>(.)\](.)/;
	$context = lc($2 . $3) . "_" . lc($1 . $4);

	$values .= $contextCount{$context} . ",";
}

$values =~ s/,$//;

my $rFile = "$outPath/$sample.Rcode";

open (R, ">$rFile") or die "Couldn't open $rFile\n";

print R ".libPaths(new=\"/.mounts/labs/PCSI/R/3.3/\")\n";
print R "library(cosmicSignnls)\n\n";

print R "df = data.frame(c($values), row.names = c($rowNames))\n";
print R "colnames(df) = c(\"$sample\")\n";

print R "bootstrapSig(df, samplecolumn=1, B=2000, s30Subset=c(1,2,3,5,6,8,13,17,18,20,26) )\n";

print R "colnames(df) = c(\"${sample}_allSigs\")\n";
print R "bootstrapSig(df, samplecolumn=1, B=2000, s30Subset=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30) )\n";

close R;


# run R script

`Rscript $rFile`;











# getBase returns the base at a specific position in the reference.  It will initialize fasta handles and pull new chunks of reference if necessary
# input is chromosome, position, reference hash and fasta handles
# output is a single base (and the reference hash and fasta handles may be modified)
sub getBase
{
    my $chr = $_[0];
    my $pos = $_[1];
    my $reference = $_[2];
    my $fastaHandles = $_[3];

    my $chunkStart = int(($pos - 1) / $reference->{"chunkSize"}) * $reference->{"chunkSize"} + 1;       # +1 because the first base in the reference is 1, $pos - 1 so that multiples of chunk size resolve to the correct chunk
    my $chunkEnd = $chunkStart + $reference->{"chunkSize"} - 1;

    unless (exists $reference->{$chr}{$chunkStart}{$pos})       # if the position isn't in our hash, we need to get a new chunk from the reference
    {
        unless (exists ($fastaHandles->{$chr}))     # create a handle for the chromosome fasta, if it doesn't exist
        {
#           warn "Creating fasta handle for $chr\n";
            $fastaHandles->{$chr} = Bio::DB::Fasta->new("$fastaHandles->{path}/$chr.fa");
        }

#       warn "Pulling $chr:$chunkStart-$chunkEnd from fasta\n";
        my $newChunk = uc($fastaHandles->{$chr}->seq($chr, $chunkStart, $chunkEnd));
        my $i = $chunkStart;
        for my $base (split("", $newChunk))
        {
            $reference->{$chr}{$chunkStart}{$i} = $base;
            $i++;
        }
    }
#   warn "returning $reference->{$chr}{$chunkStart}{$pos}\n";
    if (exists $reference->{$chr}{$chunkStart}{$pos})
    {
        return $reference->{$chr}{$chunkStart}{$pos};
    }
    else
    {
        return "N";
    }
}

# getRange returns a string of bases from the reference in the specified range by calling getBase
# input is chromosome, start pos, end pos, reference hash and fasta handles
# output is a string of bases
sub getRange
{
    my $chr = $_[0];
    my $start = $_[1];
    my $end = $_[2];
    my $reference = $_[3];
    my $fastaHandles = $_[4];

    my $seq = "";

    for (my $p = $start; $p <= $end; $p++)
    {
        $seq .= getBase($chr, $p, $reference, $fastaHandles);
#       warn "Got base: $chr:$p\t$seq\n";
    }

    return $seq;
}

# freeMemChunks tests if the next indel to be processed is on a different chromosome or more than a chunk away from the reference sequences currently in memory
#   if there exist chunks that we won't need again (assuming the input is sorted) the chunks will be removed from the reference hash
# input is the chromosome and position of the current indel, and the reference hash
# there is no output
sub freeMemChunks
{
    my $chr = $_[0];
    my $pos = $_[1];
    my $reference = $_[2];

    # delete chunks from non-current chromosomes
    for my $refChr (keys %$reference)
    {
        if (($refChr ne $chr) and ($refChr ne "chunkSize"))
        {
#           warn "deleting all chunks for $refChr.\n";
            delete $reference->{$refChr};
        }
    }

    # delete chunks if they are more than 1.5 chunks away from the current indel
    # 1.5 so that we are at least in the middle of the current chunk before dropping the previous one
    for my $chunkPos (keys %{ $reference->{$chr} })
    {
        if ($chunkPos < ($pos - (1.5 * $reference->{"chunkSize"})))
        {
#           warn "deleting $chr:$chunkPos chunk.\n";
            delete $reference->{$chr}{$chunkPos};
        }
    }

    return;
}


