#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Std;

my %opts = (
	f => '/home/zgu1/Z/refData/cds_final.pt',
	s => 2
);
getopts( 'i:o:f:s:', \%opts );
print(
	qq/
Usage:   AnnoINDEL.pl [options]
Description:
	Annotate chr:poz to gene regions.
Options: 
	-i STRING	input poz file name [required]
	-o STRING	output result file name [required]
	-f STRING	reference file directory [$opts{f}]
	-s STRING	splice site length between exon and intron [$opts{s}]
\n/
);

die("Parameters are needed!\n")
  if ( !$opts{i} || !$opts{o} || !$opts{f} );

my ( $fin, $fout, $refFile, $spliceEdge ) =
  ( $opts{i}, $opts{o}, $opts{f}, $opts{s} );
print("Parameters are:\n\t -i $fin -o $fout -f $refFile -s $spliceEdge !\n");

my $refseqFile = "$refFile";
my %refSeqHash = refSeqHash($refseqFile);

open( IN, $fin ) || die("Could not open file $fin!\n");
my $header = <IN>;
chomp $header;
open( OUT, ">$fout" ) || die("Could not create file $fout!\n");
my @titles = qw(Type Dist2Gene RefSeq Name Start End Strand Ref_codon Var_codon
  Ref_aa Var_aa Aa_poz Peptide_len Exon_poz Exon_len);
my $title = join( "\t", @titles );
print OUT "$header\t$title\n";

while (<IN>) {
	chomp;
	my ( $chr, $start, $refBase, $varBase, $end ) = split;
	$chr = "chr" . $chr if ( $chr !~ /chr/ );
	($varBase) = split( ',', $varBase );
	print OUT "$_\t";
	my $MB3        = int( $start / 3e6 );
	my @genesOnChr = ();
	if ( exists $refSeqHash{$chr}{ $MB3 - 1 } ) {
		push( @genesOnChr, @{ $refSeqHash{$chr}{ $MB3 - 1 } } );
	}
	if ( exists $refSeqHash{$chr}{$MB3} ) {
		push( @genesOnChr, @{ $refSeqHash{$chr}{$MB3} } );
	}
	if ( exists $refSeqHash{$chr}{ $MB3 + 1 } ) {
		push( @genesOnChr, @{ $refSeqHash{$chr}{ $MB3 + 1 } } );
	}
	if (@genesOnChr) {
		my $annoGeneRef =
		  indel_Annotation( $start, $refBase, $varBase, $end, $spliceEdge,
			\@genesOnChr );
		print OUT join( "\t", @$annoGeneRef ), "\n";
	}
	else {
		print OUT join( "\t", (0) x @titles ), "\n";
	}
}
close OUT;
close IN;

##########################################################################################################################################################################
sub indel_Annotation {
	my ( $start, $refBase, $varBase, $end, $spliceEdge, $genesOnChr ) = @_;
	my @pozs = $start .. $end;
	my @refs = split( //, $refBase );
	die("inconsistent refBase and $start-$end poz!\n")
	  unless ( @refs == @pozs );
	my ( $typeIdx, $aaLenIdx ) = ( 0, 12 );
	my @allPozAnno = ();

	for my $i ( 0 .. $#pozs ) {
		my $ipoz = $pozs[$i];
		my $iref = $refs[$i];
		my ( $nearestGene, @insideGenes ) = relatedGenes( $ipoz, $genesOnChr );
		my $outGene;
		if ( @insideGenes > 0 ) {
			my @poz2genes;
			for my $geneRef (@insideGenes) {
				my $pozAnno =
				  poz2intragene( $ipoz, $iref, $iref, $spliceEdge, $geneRef );
				push( @poz2genes, $pozAnno );
			}
			$outGene = tag( $typeIdx, $aaLenIdx, @poz2genes );
		}
		elsif ($nearestGene) {
			$outGene = poz2intergene( $ipoz, $nearestGene );
		}
		else {
			die("Error pos location in Annotation!\n");
		}
		push( @allPozAnno, $outGene );
	}
	return mergeAnno( $refBase, $varBase, @allPozAnno );
}

sub mergeAnno {
	my ( $refBase, $varBase, @geneAnno ) = @_;
	my ( $typeIdx, $aaPozIdx, $aaLenIdx ) = ( 0, 11, 12 );
	my $refLen = length($refBase);
	my $varLen = length($varBase);
	my $varType =
	  abs( $refLen - $varLen ) % 3 == 0 ? "in-frame" : "frame_shift";
	my $insDel = $refLen > $varLen ? "del" : "ins";
	$varType = $varType . "_$insDel";
	for my $geneRef (@geneAnno) {
		my ( $type, $aaPoz ) = @$geneRef[ $typeIdx, $aaPozIdx ];
		chop($aaPoz);
		if ( $type eq 'synonymous' ) {
			$$geneRef[$typeIdx] = $varType;
			if ( $varType =~ /frame_shift/ ) {
				$$geneRef[$aaPozIdx] = $aaPoz . 'fs';
			}
			else {
				$$geneRef[$aaPozIdx] = $aaPoz . $insDel;
			}
		}
	}
	return tagIndel( $typeIdx, $aaPozIdx, $aaLenIdx, @geneAnno );
}

sub tagIndel {
	my ( $typeIdx, $aaPozIdx, $aaLenIdx, @items ) = @_;
	my %funHs = (
		'splice'          => 1,
		'frame_shift_del' => 2,
		'frame_shift_ins' => 2,
		'in-frame_del'    => 3,
		'in-frame_ins'    => 3,
		'intron-NS'       => 3.5,
		'utr5'            => 4,
		'utr3'            => 5,
		'exon'            => 6,
		'intron'          => 7,
		'intergenic'      => 8
	);
	my @tmp = ();
	for (@items) {
		my ( $type, $aaPoz, $aaLen ) = @{$_}[ $typeIdx, $aaPozIdx, $aaLenIdx ];
		my $funPri = $funHs{$type};
		my ($aaPozNum)=$aaPoz=~/(\d+)/;
		push( @tmp, [ $funPri, $aaLen, $aaPozNum, $_ ] );
	}
	@tmp = sort mySortIndel @tmp;
	return $tmp[0]->[3];
}

sub mySortIndel{
	my ( $priA, $lenA, $aapozA ) = @$a;
	my ( $priB, $lenB, $aapozB ) = @$b;
	$priA <=> $priB || $lenB <=> $lenA || $aapozA <=> $aapozB;
	
}

sub tag {
	my ( $typeIdx, $aaLenIdx, @items ) = @_;
	my %funHs = (
		'splice'     => 1,
		'synonymous' => 2,
		'intron-NS'  => 2.5,
		'utr5'       => 3,
		'utr3'       => 4,
		'exon'       => 5,
		'intron'     => 6,
	);
	my @tmp = ();
	for (@items) {
		my ( $fun, $protLen ) = @{$_}[ $typeIdx, $aaLenIdx ];
		my $funPri = $funHs{$fun};
		push( @tmp, [ $funPri, $protLen, $_ ] );
	}
	@tmp = sort mySort @tmp;
	return $tmp[0]->[2];
}
########################################################################################################
#Make hash table for refseq data
sub refSeqHash {
	my ($fname) = @_;
	open( REF, $fname ) || die("Can not open file $fname!\n");
	my %ref;
	while (<REF>) {
		chomp;
		my ( $NM, $gene, $chr, $strand, $start, $end, $range, $seq ) = split;
		my $MB3 = int( $start / 3e6 );
		push(
			@{ $ref{$chr}{$MB3} },
			[ $NM, $gene, $strand, $start, $end, $range, $seq ]
		);
	}
	close REF;
	return %ref;
}

#Walk throuth all the genes in the genes array to select related ones,
#if poz in at intergenic postion, return the nearest gene
#else if poz in several genes (or a gene), return all these genes' information
sub relatedGenes {
	my ( $poz, $ref ) = @_;
	my $nearestGene;
	my $distance = 1e9;    # Initialization of the distance
	my @insideGenes;
	for my $NMRef (@$ref) {
		my ( $start, $end ) = @{$NMRef}[ 3, 4 ];
		my $flag = inRange( $poz, $start, $end );
		if ( $flag == 1 ) {
			my $innerDist = $poz - $end;
			if ( $innerDist < $distance ) {
				$nearestGene = $NMRef;
				$distance    = $innerDist;
			}
		}
		elsif ( $flag == 0 ) {
			push( @insideGenes, $NMRef );
		}
		else {
			my $innerDist = $start - $poz;
			if ( $innerDist < $distance ) {
				$nearestGene = $NMRef;
				$distance    = $innerDist;
			}
		}
	}
	return ( $nearestGene, @insideGenes );
}

# Annotate the intragenic mutation information
sub poz2intragene {
	my ( $poz, $refBase, $varBase, $spliceEdge, $geneRef ) = @_;
	my (
		$varType,   $distance, $NM,      $gene,
		$geneStart, $geneEnd,  $strand,  $cd1,
		$cd2,       $aa1,      $aa2,     $aaPoz,
		$aaLen,     $exonPoz,  $exonLen, $adjacentExon,
		$adjacentDist
	) = (0) x 17;

	# [$NM, $gene, $strand, $start, $end, $range, $seq]
	( $NM, $gene, $strand, $geneStart, $geneEnd ) = @{$geneRef}[ 0 .. 4 ];
	my ( $range, $cds ) = @{$geneRef}[ 5, 6 ];

	my @allExon = split( /\:/, $range );
	my @utr5 = grep /\d\|\d/,     @allExon;
	my @cds  = grep /\d\|\|\d/,   @allExon;
	my @utr3 = grep /\d\|\|\|\d/, @allExon;
	my $cdsSeqLen = length($cds);    #total CDS length
	$aaLen = int( $cdsSeqLen / 3 ) - 1;    #sustract 1 stop condon
	$aaLen += 1
	  if ( $cdsSeqLen % 3 != 0 )
	  ;    #if the cds len is not 3*n, take in the last qulified codon

	#if utr5 is connected with cds, $u5c (or $u3c) is 1, or it is 0
	my ( $u5c, $u3c ) = (0) x 2;
	if ( @utr5 && @cds ) {
		$u5c = utrcds( $strand, \@utr5, \@cds );
	}
	if (@utr3) {
		$u3c = utrcds( $strand, \@cds, \@utr3 );
	}
	$exonLen = @allExon - $u5c - $u3c;

	my ( $type, $geneType, $exonNum ) =
	  inWhichType( $poz, \@utr5, \@cds, \@utr3 );
	$varType = 'intron' if ( $geneType eq 'intron' );

	if ( $type eq 'utr5' ) {
		if ( $geneType eq 'intron' ) {
			( $adjacentExon, $adjacentDist ) =
			  adjacent( $poz, $strand, @utr5[ $exonNum - 1, $exonNum ] );
			if ( $adjacentExon == 1 ) {
				$adjacentExon = $exonNum;
				$exonPoz      = "e$exonNum+$adjacentDist";
			}
			elsif ( $adjacentExon == 2 ) {
				$adjacentExon = $exonNum + 1;
				$exonPoz      = "e$exonNum-$adjacentDist";
			}
		}
		elsif ( $geneType eq 'exon' ) {
			$exonPoz = "e$exonNum";
			$varType = $#cds >= 0 ? 'utr5' : 'exon';
		}
	}
	elsif ( $type eq 'cds' ) {
		my $deltaExonNum = @utr5 - $u5c;
		my $trueExonNum  = $exonNum + $deltaExonNum;
		if ( $geneType eq 'intron' ) {
			( $adjacentExon, $adjacentDist ) =
			  adjacent( $poz, $strand, @cds[ $exonNum - 1, $exonNum ] );
			if ( $adjacentExon == 1 ) {
				$adjacentExon = $trueExonNum;
				$exonPoz      = "e$trueExonNum+$adjacentDist";
			}
			elsif ( $adjacentExon == 2 ) {
				$adjacentExon = $trueExonNum + 1;
				$exonPoz      = "e$trueExonNum-$adjacentDist";
			}
			if ( $adjacentDist <= $spliceEdge ) {
				my $spExon = $adjacentExon == 1 ? $exonNum - 1 : $exonNum;
				( $cd1, $aa1, $aaPoz ) = spAnno( $spExon, $cds, @cds );
			}
		}
		elsif ( $geneType eq 'exon' ) {
			$exonPoz = "e$trueExonNum";
			( $cd1, $cd2, $aa1, $aa2, $aaPoz ) =
			  cdsAnno( $strand, $poz, $varBase, $cds, @cds );
			my ($aaPozN) = $aaPoz =~ /(\d+)/;
			$varType =
			    $aa1 eq $aa2 ? 'synonymous'
			  : $aa2 eq '*'  ? 'nonsense'
			  : $aa1 eq '*'  ? 'stoploss'
			  : $aaPozN == 1 ? 'startloss'
			  :                'missense';
		}
	}
	elsif ( $type eq 'utr3' ) {
		my $deltaExonNum = @cds + @utr5 - $u5c - $u3c;
		my $trueExonNum  = $exonNum + $deltaExonNum;
		if ( $geneType eq 'intron' ) {
			( $adjacentExon, $adjacentDist ) =
			  adjacent( $poz, $strand, @utr3[ $exonNum - 1, $exonNum ] );
			if ( $adjacentExon == 1 ) {
				$adjacentExon = $trueExonNum;
				$exonPoz      = "e$trueExonNum+$adjacentDist";
			}
			elsif ( $adjacentExon == 2 ) {
				$adjacentExon = $trueExonNum + 1;
				$exonPoz      = "e$trueExonNum-$adjacentDist";
			}
		}
		elsif ( $geneType eq 'exon' ) {
			$exonPoz = "e$trueExonNum";
			$varType = 'utr3';
		}
	}
	elsif ( $type eq 'utr5-cds' ) {
		$exonNum = $#utr5 + 1;
		( $adjacentExon, $adjacentDist ) =
		  adjacent( $poz, $strand, $utr5[$#utr5], $cds[0] );
		if ( $adjacentExon == 1 ) {
			$adjacentExon = $exonNum;
			$exonPoz      = "e$exonNum+$adjacentDist";
		}
		elsif ( $adjacentExon == 2 ) {
			$adjacentExon = $exonNum + 1;
			$exonPoz      = "e$exonNum-$adjacentDist";
		}
	}
	elsif ( $type eq 'cds-utr3' ) {
		$exonNum = @utr5 + @cds - $u5c;
		( $adjacentExon, $adjacentDist ) =
		  adjacent( $poz, $strand, $cds[$#cds], $utr3[0] );
		if ( $adjacentExon == 1 ) {
			$adjacentExon = $exonNum;
			$exonPoz      = "e$exonNum+$adjacentDist";
		}
		elsif ( $adjacentExon == 2 ) {
			$adjacentExon = $exonNum + 1;
			$exonPoz      = "e$exonNum-$adjacentDist";
		}
	}

	#adjacentDist == 0 means the position in exon
	if ( $adjacentDist > 0 ){
  		if ($adjacentDist <= $spliceEdge ){
  			$varType = 'splice';
  		}elsif($adjacentDist <= 10){
  			$varType = $varType."-NS";
  		}
	}

	return [
		$varType, $distance, $NM,    $gene,    $geneStart,
		$geneEnd, $strand,   $cd1,   $cd2,     $aa1,
		$aa2,     $aaPoz,    $aaLen, $exonPoz, $exonLen
	];
}

#Annotate the intergenic mutation information
sub poz2intergene {
	my ( $poz, $nearestGene ) = @_;
	my (
		$varType, $distance, $NM,    $gene,    $geneStart,
		$geneEnd, $strand,   $cd1,   $cd2,     $aa1,
		$aa2,     $aaPoz,    $aaLen, $exonPoz, $exonLen
	) = (0) x 17;
	$varType = 'intergenic';

	# $NM, $gene, $strand, $start, $end, $range, $seq
	( $NM, $gene, $strand, $geneStart, $geneEnd ) = @{$nearestGene}[ 0 .. 4 ];
	my $flag = inRange( $poz, $geneStart, $geneEnd );
	if ( $flag == 1 ) {
		$distance = $poz - $geneEnd;
	}
	elsif ( $flag == -1 ) {
		$distance = $geneStart - $poz;
	}
	else {
		die("Error for intergenic annotation!\n");
	}
	return [
		$varType, $distance, $NM,    $gene,    $geneStart,
		$geneEnd, $strand,   $cd1,   $cd2,     $aa1,
		$aa2,     $aaPoz,    $aaLen, $exonPoz, $exonLen
	];
}

#calculate the aa length and the last aa for splice mutation
sub spAnno {
	my ( $spExon, $cds, @cds ) = @_;
	my $cdsLen = 0;
	for my $i ( 0 .. $spExon - 1 ) {
		my ( $s, $e ) = split( /\|\|/, $cds[$i] );
		$cdsLen += abs( $s - $e ) + 1;
	}
	my $aaPoz = int( $cdsLen / 3 );
	my $cd1   = substr( $cds, $aaPoz * 3 - 3, 3 );
	my $aa1   = codon2aa($cd1);
	my $aa1_s = codon2aa_s($cd1);
	$aaPoz = "p.$aa1_s${aaPoz}sp";
	return ( $cd1, $aa1, $aaPoz );
}

sub cdsAnno {
	my ( $strand, $poz, $varBase, $cds, @cds ) = @_;
	my ( $cd1, $cd2, $aa1, $aa2, $aaPoz, $aa1_s, $aa2_s );

	my $relativePoz = position( $strand, $poz, @cds );
	my $mod = $relativePoz % 3;
	$varBase =~ tr/ATCG/TAGC/ if ( $strand eq "-" );
	if ( $mod != 0 ) {
		$cd2 = $cd1 = substr( $cds, $relativePoz - $mod, 3 );
		substr( $cd2, $mod - 1, 1, $varBase );
	}
	else {
		$cd2 = $cd1 = substr( $cds, $relativePoz - 3, 3 );
		substr( $cd2, 2, 1, $varBase );
	}
	if ( length($cd1) == 3 and length($cd2) == 3 ) {
		$aa1   = codon2aa($cd1);
		$aa1_s = codon2aa_s($cd1);
		$aa2   = codon2aa($cd2);
		$aa2_s = codon2aa_s($cd2);
	}
	else {    #for the tail of the cds length is not 3*n,
		( $aa1, $aa2, $aa1_s, $aa2_s ) = ( "NA", "_", "NA", "_" );
	}
	$aaPoz = int( ( $relativePoz - 1 ) / 3 ) + 1;
	$aaPoz = "p.$aa1_s$aaPoz$aa2_s";
	return ( $cd1, $cd2, $aa1, $aa2, $aaPoz );
}

sub adjacent {
	my ( $site, $strand, $preExon, $proExon ) = @_;
	my ( $adjacentExon, $adjacentDist ) = ( 0, 0 );
	my ( $preStart, $preEnd ) = split( /\|+/, $preExon );
	my ( $proStart, $proEnd ) = split( /\|+/, $proExon );
	if ( $strand eq '+' ) {
		my $preDist = $site - $preEnd;
		my $proDist = $proStart - $site;
		if ( $preDist <= $proDist ) {
			( $adjacentExon, $adjacentDist ) = ( 1, $preDist );
		}
		else {
			( $adjacentExon, $adjacentDist ) = ( 2, $proDist );
		}
	}
	else {
		my $preDist = $preStart - $site;
		my $proDist = $site - $proEnd;
		if ( $preDist <= $proDist ) {
			( $adjacentExon, $adjacentDist ) = ( 1, $preDist );
		}
		else {
			( $adjacentExon, $adjacentDist ) = ( 2, $proDist );
		}
	}
	return ( $adjacentExon, $adjacentDist );
}

sub utrcds {
	my ( $strand, $pre, $post ) = @_;
	my @pre   = @$pre;
	my @post  = @$post;
	my @last  = split( /\|+/, $pre[$#pre] );
	my @first = split( /\|+/, $post[0] );
	if ( $strand eq '+' ) {
		return ( $first[0] - $last[1] ) == 1 ? 1 : 0;
	}
	else {
		return ( $last[0] - $first[1] ) == 1 ? 1 : 0;
	}
}

sub inWhichType {
	my ( $poz, $utr5Ref, $cdsRef, $utr3Ref ) = @_;
	my ( $isU5,  $u5Num )  = inWhichRange( $poz, @$utr5Ref );
	my ( $isCds, $cdsNum ) = inWhichRange( $poz, @$cdsRef );
	my ( $isU3,  $u3Num )  = inWhichRange( $poz, @$utr3Ref );

	my $geneType = 'intron';
	if ($isU5) {
		my $type = 'utr5';
		if ( $isU5 == 1 ) {
			$geneType = 'exon';
		}
		return ( $type, $geneType, $u5Num );
	}
	elsif ($isCds) {
		my $type = 'cds';
		if ( $isCds == 1 ) {
			$geneType = 'exon';
		}
		return ( $type, $geneType, $cdsNum );
	}
	elsif ($isU3) {
		my $type = 'utr3';
		if ( $isU3 == 1 ) {
			$geneType = 'exon';
		}
		return ( $type, $geneType, $u3Num );
	}

	my $utr5Cds = $u5Num * $cdsNum;
	my $cdsUtr3 = $cdsNum * $u3Num;
	if ( $utr5Cds == -1 ) {
		my $type = 'utr5-cds';
		my $exonNum = return ( $type, $geneType, -1 );
	}
	elsif ( $cdsUtr3 == -1 ) {
		my $type = 'cds-utr3';
		return ( $type, $geneType, -1 );
	}
	else {
		die("Error position for $poz!\n");
	}
}

sub inWhichRange {
	my ( $poz, @ranges ) = @_;
	my $sign = 0;
	for my $i ( 0 .. $#ranges ) {
		my $flag = inRange( $poz, split( /\|+/, $ranges[$i] ) );
		$sign = $flag if ( $i == 0 );
		if ( $flag == 0 ) {

			# find and return the region
			return ( 1, $i + 1 );
		}
		elsif ( $flag * $sign == -1 ) {

			# not find, but between regions; return the 5' exon
			return ( 2, $i );
		}
	}

	# not found, return the sign of comparison for all the regions
	return ( 0, $sign );
}

sub position {
	my ( $strand, $poz, @range ) = @_;
	my $sum = 0;
	for my $i ( 0 .. $#range ) {
		my ( $s, $e ) = split( /\|\|/, $range[$i] );
		my $flag = inRange( $poz, $s, $e );
		if ( $flag == 0 ) {
			if ( $strand eq "+" ) {
				return $sum + $poz - $s + 1;
			}
			else {
				return $sum + $e - $poz + 1;
			}
		}
		else {
			$sum = $sum + $e - $s + 1;
		}
	}
}

sub inRange {
	my ( $p, $s, $e ) = @_;
	return $p < $s ? -1 : $p > $e ? 1 : 0;
}

sub codon2aa {
	my ($codon) = @_;
	$codon = uc $codon;
	my (%codon2aa) = (
		"TCA" => "Ser",    # Serine
		"TCC" => "Ser",    # Serine
		"TCG" => "Ser",    # Serine
		"TCT" => "Ser",    # Serine
		"TTC" => "Phe",    # Phenylalanine
		"TTT" => "Phe",    # Phenylalanine
		"TTA" => "Leu",    # Leucine
		"TTG" => "Leu",    # Leucine
		"TAC" => "Tyr",    # Tyrosine
		"TAT" => "Tyr",    # Tyrosine
		"TAA" => "*",      # Stop
		"TAG" => "*",      # Stop
		"TGC" => "Cys",    # Cysteine
		"TGT" => "Cys",    # Cysteine
		"TGA" => "*",      # Stop
		"TGG" => "Trp",    # Tryptophan
		"CTA" => "Leu",    # Leucine
		"CTC" => "Leu",    # Leucine
		"CTG" => "Leu",    # Leucine
		"CTT" => "Leu",    # Leucine
		"CCA" => "Pro",    # Proline
		"CCC" => "Pro",    # Proline
		"CCG" => "Pro",    # Proline
		"CCT" => "Pro",    # Proline
		"CAC" => "His",    # Histidine
		"CAT" => "His",    # Histidine
		"CAA" => "Gln",    # Glutamine
		"CAG" => "Gln",    # Glutamine
		"CGA" => "Arg",    # Arginine
		"CGC" => "Arg",    # Arginine
		"CGG" => "Arg",    # Arginine
		"CGT" => "Arg",    # Arginine
		"ATA" => "Ile",    # Isoleucine
		"ATC" => "Ile",    # Isoleucine
		"ATT" => "Ile",    # Isoleucine
		"ATG" => "Met",    # Methionine
		"ACA" => "Thr",    # Threonine
		"ACC" => "Thr",    # Threonine
		"ACG" => "Thr",    # Threonine
		"ACT" => "Thr",    # Threonine
		"AAC" => "Asn",    # Asparagine
		"AAT" => "Asn",    # Asparagine
		"AAA" => "Lys",    # Lysine
		"AAG" => "Lys",    # Lysine
		"AGC" => "Ser",    # Serine
		"AGT" => "Ser",    # Serine
		"AGA" => "Arg",    # Arginine
		"AGG" => "Arg",    # Arginine
		"GTA" => "Val",    # Valine
		"GTC" => "Val",    # Valine
		"GTG" => "Val",    # Valine
		"GTT" => "Val",    # Valine
		"GCA" => "Ala",    # Alanine
		"GCC" => "Ala",    # Alanine
		"GCG" => "Ala",    # Alanine
		"GCT" => "Ala",    # Alanine
		"GAC" => "Asp",    # Aspartic Acid
		"GAT" => "Asp",    # Aspartic Acid
		"GAA" => "Glu",    # Glutamic Acid
		"GAG" => "Glu",    # Glutamic Acid
		"GGA" => "Gly",    # Glycine
		"GGC" => "Gly",    # Glycine
		"GGG" => "Gly",    # Glycine
		"GGT" => "Gly",    # Glycine
	);
	if ( exists $codon2aa{$codon} ) {
		return $codon2aa{$codon};
	}
	else {
		return "NULL";
	}
}

sub codon2aa_s {
	my ($codon) = @_;
	$codon = uc $codon;
	my (%genetic_code) = (

		'TCA' => 'S',    # Serine
		'TCC' => 'S',    # Serine
		'TCG' => 'S',    # Serine
		'TCT' => 'S',    # Serine
		'TTC' => 'F',    # Phenylalanine
		'TTT' => 'F',    # Phenylalanine
		'TTA' => 'L',    # Leucine
		'TTG' => 'L',    # Leucine
		'TAC' => 'Y',    # Tyrosine
		'TAT' => 'Y',    # Tyrosine
		'TAA' => '*',    # Stop
		'TAG' => '*',    # Stop
		'TGC' => 'C',    # Cysteine
		'TGT' => 'C',    # Cysteine
		'TGA' => '*',    # Stop
		'TGG' => 'W',    # Tryptophan
		'CTA' => 'L',    # Leucine
		'CTC' => 'L',    # Leucine
		'CTG' => 'L',    # Leucine
		'CTT' => 'L',    # Leucine
		'CCA' => 'P',    # Proline
		'CCC' => 'P',    # Proline
		'CCG' => 'P',    # Proline
		'CCT' => 'P',    # Proline
		'CAC' => 'H',    # Histidine
		'CAT' => 'H',    # Histidine
		'CAA' => 'Q',    # Glutamine
		'CAG' => 'Q',    # Glutamine
		'CGA' => 'R',    # Arginine
		'CGC' => 'R',    # Arginine
		'CGG' => 'R',    # Arginine
		'CGT' => 'R',    # Arginine
		'ATA' => 'I',    # Isoleucine
		'ATC' => 'I',    # Isoleucine
		'ATT' => 'I',    # Isoleucine
		'ATG' => 'M',    # Methionine
		'ACA' => 'T',    # Threonine
		'ACC' => 'T',    # Threonine
		'ACG' => 'T',    # Threonine
		'ACT' => 'T',    # Threonine
		'AAC' => 'N',    # Asparagine
		'AAT' => 'N',    # Asparagine
		'AAA' => 'K',    # Lysine
		'AAG' => 'K',    # Lysine
		'AGC' => 'S',    # Serine
		'AGT' => 'S',    # Serine
		'AGA' => 'R',    # Arginine
		'AGG' => 'R',    # Arginine
		'GTA' => 'V',    # Valine
		'GTC' => 'V',    # Valine
		'GTG' => 'V',    # Valine
		'GTT' => 'V',    # Valine
		'GCA' => 'A',    # Alanine
		'GCC' => 'A',    # Alanine
		'GCG' => 'A',    # Alanine
		'GCT' => 'A',    # Alanine
		'GAC' => 'D',    # Aspartic Acid
		'GAT' => 'D',    # Aspartic Acid
		'GAA' => 'E',    # Glutamic Acid
		'GAG' => 'E',    # Glutamic Acid
		'GGA' => 'G',    # Glycine
		'GGC' => 'G',    # Glycine
		'GGG' => 'G',    # Glycine
		'GGT' => 'G',    # Glycine
	);
	if ( exists $genetic_code{$codon} ) {
		return $genetic_code{$codon};
	}
	else {
		return 'NULL';
	}
}

sub mySort {
	my ( $priA, $lenA ) = @$a;
	my ( $priB, $lenB ) = @$b;
	$priA <=> $priB || $lenB <=> $lenA;
}
