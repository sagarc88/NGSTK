#Sagar Chhangawala
#! usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;
#Allows bundling of options (ex. -rc instead of -r -c)
Getopt::Long::Configure('bundling'); 

&usage if(@ARGV < 1); #Show usage if no input
my %opt=(); #hash to hold all the options. 
my $result = GetOptions (\%opt, "r", "c","f","s","i=s","w=i","h|help");
&usage if ($opt{h});
&main;
exit;

# Reads the fasta file and passes the seq to another subroutine
# for processing.
sub main{
  my $seq="";
  my $maxlen=80;
  my $filename=$opt{i} if (defined $opt{i});
  my $header;
  my $seqout;
  
  $maxlen=$opt{w} if (defined $opt{w});

  $ARGV[0]=$filename if(defined $opt{i}); #change to filehandle if no STDIN. 

  while (<>) {
    next if(/^$/);
    if (/>/) {
	if ($seq) {    
	  &_revcom($header,\$seq,$maxlen);
	}      
	$header=$_;
	$header =~ s/^>//; # remove ">"
	$header =~ s/\s+$//; # remove trailing whitespace

	$seq = ""; # clear out old sequence
    }         
    else {    
	s/\s+//g; # remove whitespace
	$seq .= $_; # add sequence
    }         
  }    

  if ($seq) { # handle last sequence
    &_revcom($header,\$seq,$maxlen);
  }
}

#print the sequence with the specified length
sub _printseq {
  my ($header,$seq_ref,$length) = @_;
  my $seq=$$seq_ref;
  print "$header\n";
  for ( my $pos = 0 ; $pos < length($seq) ; $pos += $length ) {
      print substr($seq, $pos, $length), "\n";
  }
}

#Reverse complement according to the specified inputs. 
sub _revcom{
  my $header=shift;
  my $dna;
  my $dna_ref=shift;
  $dna=$$dna_ref; #dereference
  my $len=shift;
  
  ($dna=~m/[?\.]/)? die(qq/Sequence cannot contain . or ?/):0;
	
  if($dna=~ m/[^ACGTRYWSKMN\-]/i){
    if($opt{s}){$dna=~s/[^ACGTRYSWKMN\-]/\-/gi;}
      elsif(!$opt{f}){
	my @nonchar=($dna =~ m/[^ACGTRYWSKMN-]/gi);
	  my $error="Sequence contains non-standard characters: ".join(",",@nonchar).
		      "\nUse -f or -s to disregard/remove the characters\n";
	  die(qq/$error/);
      }
  }	
  $dna= (reverse $dna) if ($opt{r}); #reverse only if specified 
  $dna =~ tr/ACGTRYSWKMBDHVNacgtryswkmbdhvn/TGCAYRSWMKVHDBNtgcayrswmkvhdbn/ if($opt{c}); #translate if specified

  &_printseq(">".$header, \$dna, $len); #call prinseq to output sequence to STOUT.	
}

#print usage
sub usage {
die(qq/Version: 1.0
Input: Fasta file (or multi-fasta). Sequences can be uppercase,lowercase or mixed.
Output: Fasta sequence to STDOUT.
Accepted IUPAC Standard characters: ACGTRYSWKMNacgtryswkn-

Usage:   perl revcom.pl [<options>] -i filename or -(if piped)\n
Options: -r 	Reverse the input sequences
	 -c 	Complement the input sequences
	 -f 	Force reverse or complement (used if sequences contain non-IUPAC characters). 
		Will leave the unrecognized characters in the sequence. 
	 -s 	Remove the non-IUPAC characters and then get the Reverse or complement of sequence.
	 -i 	Input file name (must be specified to use file as input)
	 -w 	Specify the length of nucleotides per line in output fasta sequence (default 80)
		(By default it will leave the other characters in its place. use -s to replace non-ACGT characters with -)
	 -h|help Displays this help manual

Examples:
 
To get a reverse complement of input file:
perl revcom.pl -rc -i dna.fasta 
  OR 
perl revcom.pl -rci dna.fasta

To get just a complement of a piped fasta sequence:
perl revcom.pl -c -

If you get an error saying the sequence contains an unknown character, use -f or -s 
perl revcom.pl -rcs -i dna.fasta
\n/);
}
