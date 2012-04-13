#Sagar Chhangawala
use warnings;
use strict;
use Getopt::Long;
Getopt::Long::Configure('bundling'); 

my %opt=(); #hash to hold all the options. 
my @files; #holds the names of files to process
my @genes; #holds the names of genes to process
my $result = GetOptions (\%opt, "i=s", 
			  "in=s"=>\@files,
			  "genes=s"=>\@genes,
			  "b=s");
&main();

sub main{
  my @bed; #bed file which holds the coordinates.
  my %f1=(); #headers and sequences of file1
  my %f2=(); #headers and sequences of file2
  
  #Do every combination of files specified
  for(my $i=0;$i<scalar(@files);$i++){
    for(my $j=$i+1;$j<scalar(@files);$j++){
	my($f1_r)=&read_fasta($files[$i]); #read the fasta file and put seqs in %f1
	%f1=%$f1_r; #dereference
	my ($f2_r)=&read_fasta($files[$j]);#read the fasta file and put seqs in %f2
	%f2=%$f2_r; #dereference

	my($bed_ref)=&makebed(\%f1,\%f2,) if(!defined($opt{b})); #make the bedfile if not specified.
	@bed=@$bed_ref; #dereference
    }
  }
}

#subroutine to read the fata file
#input - filename.
#output -reference to hash w/ seqs n headers.
sub read_fasta{
  my $filename=shift;
  my $header;
  my $seqout;

  my %seq=();
  my $seq="";

  open(IN,$filename);

  while (<IN>) {
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
  return(\%seq);
}

#Subroutine to make the bed file if not already given by user
#input - hashes of f1 and f2
#output - reference to hash of bed array.
sub makebed{
  my $f1_ref=shift;
  my $f2_ref=shift;

  my %f1=%$f1_ref;
  my %f2=%$f2_ref;
  my @bed;
  my $i=0;
  my $j=0;
  
  #if the gene names are inputed by user
  if(@genes){
    foreach my $gene(@genes){ #make bed entry for each gene
	#check to make sure the gene name exists in both files and seqs length are same.
	if(($f1{$gene} && $f2{$gene})&& (length($f1{$gene})==length($f2{$gene}))){
	    $bed[$i][$j]="$gene";
	    $bed[$i][$j+1]="0\t".length($f1{$gene});
	    $i++;
      }
    }
  }
  return(\@bed);
}

sub usage{

die(qq/Version 1.0
use --in to specify file names. 
ex. --in file1 --in file2

use --genes to specify gene names
ex. --genes chr1 --genes chr2

\n/);
}