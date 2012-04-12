#Sagar Chhangawala
use warnings;
use strict;
use Getopt::Long;
Getopt::Long::Configure('bundling'); 

my %opt=(); #hash to hold all the options. 
my $result = GetOptions (\%opt, "i=s");
my %seqs1=();
my %seqs2=();

sub main{
&read_fasta(%seq1);

}

sub read_fasta{
  my $seq_ref=shift;
  my %seq=%$seq_ref;
  my $seq="";
#  my $maxlen=80;
  my $filename=$opt{i} if (defined $opt{i});
  my $header;
  my $seqout;
  
  $maxlen=$opt{w} if (defined $opt{w});

  $ARGV[0]=$filename if(defined $opt{i}); #change to filehandle if no STDIN. 

  while (<>) {
    next if(/^$/);
    if (/>/) {
	if ($seq) {
	  $seq{$header}=$seq;
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
    $seq{$header}=$seq;    
  }
}