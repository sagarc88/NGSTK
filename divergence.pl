#Written by Sagar Chhangawala

use warnings;
use strict;
use Getopt::Long;
Getopt::Long::Configure('bundling'); 

&seq2dist();

sub seq2dist{
my $usage='
Version 1.0
Input: Atleast 2 fasta files containing the sequences of interest. Both files must have the same headers for the genes of interest. 
Output: Table of following columns: GeneID,ts,tv,p-dist,jc_dist,and kimura_dist

Options:
  --in: 	To specify file names. 
		ex. --in file1 --in file2

  --genes:	To specify gene names
		ex. --genes chr1 --genes chr2

  --bed:	To input bed files. The tool only takes 1 bed file. 
		ex. --bed bed.txt.

Example commands:
To get the statistics for all the genes in 2 file
perl divergence.pl --in file1 --in file2 

To get the statistics for the specified coordinates in the bed file. 
perl divergence.pl --in file1 --in file2 --bed bed.txt

To get the statistics for the genes of interest. 
perl divergence.pl --in file1 --in file2 --genes cat --genes mouse

To get the statistics for the genes of interest with specified coordinates in bed file
perl divergence.pl --in file1 --in file2 --genes cat --genes mouse --bed bed.txt';

  if(@ARGV < 1){die($usage);} #Show usage if no input

  my %opt=(); #hash to hold all the options. 
  my @files; #holds the names of files to process
  my @genes; #holds the names of genes to process
  my $result = GetOptions (\%opt, "bed=s", "h|help",
			    "in=s"=>\@files,
			    "gene=s"=>\@genes);
  if(!$result){die("exited due to incorrect options");} #exit program if any problems with options.
  if ($opt{h}){die($usage);}
  die("No files specified. Use --file option to specify file") unless (@files >1); ##die if more than 1 file is not specified.

  & _run_seqs2dist(\%opt,\@files,\@genes);

}

sub _run_seqs2dist{
  my $opt_ref=shift;
  my %opt=%$opt_ref;
  my $files_ref=shift;
  my @files=@$files_ref;
  my $genes_ref=shift;
  my @genes=@$genes_ref;

  my @bed; #bed file which holds the coordinates.
  my %f1=(); #headers and sequences of file1
  my %f2=(); #headers and sequences of file2
  print "GeneID\tts\ttv\tp-dist\tjc_dist\tkimura_dist\n";
  #Do every combination of files specified
  for(my $i=0;$i<scalar(@files);$i++){
    for(my $j=$i+1;$j<scalar(@files);$j++){
	my($f1_r)=&_read_fasta($files[$i],\%opt); #read the fasta file and put seqs in %f1
	%f1=%$f1_r; #dereference
	my ($f2_r)=&_read_fasta($files[$j],\%opt);#read the fasta file and put seqs in %f2
	%f2=%$f2_r; #dereference

	my($bed_ref)=&_makebed(\%f1,\%f2,\%opt,\@genes); #make the bedfile if not specified.
	@bed=@$bed_ref; #dereference

	&_run_analysis(\%f1, \%f2, \@bed, $files[$i],$files[$j],\%opt);
    }
  }
}

#subroutine to read the fata file
#input - filename.
#output -reference to hash w/ seqs n headers.
sub _read_fasta{
  my $filename=shift;
  my $opt_ref=shift;
  my %opt=%$opt_ref;
  my $header;
  my $seqout;

  my %seq=();
  my $seq="";

  open(IN,$filename) or die "Can't open $filename: $!";

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

#Subroutine to make the bed array.
#input - hashes of f1 and f2
#output - reference to hash of bed array.
sub _makebed{
  my $f1_ref=shift;
  my $f2_ref=shift;
  my $opt_ref=shift;
  my %opt=%$opt_ref;
  my $genes_ref=shift;
  my @genes=@$genes_ref; 

  my %f1=%$f1_ref;
  my %f2=%$f2_ref;
  my @bed;
  my $i=0;
  my $j=0;
  
  #if the bedfile is provided, read into array
  if(defined($opt{b})){
    open(BED,$opt{b}) or die "Can't open $opt{b}: $!";
    while(<BED>){
      chomp;
      my @t=split(/\t/,$_);
      #check if the genes exist in both files and the length of sequences
      if(($f1{$t[0]} && $f2{$t[0]})&& (length($f1{$t[0]})==length($f2{$t[0]}))){
	  $bed[$i][$j]="$t[0]";
	  $bed[$i][$j+1]="$t[1]\t$t[2]";
	  $i++;
      }
      else{die("$t[0] either has no partner or the sequence length are different. Please correct the fasta file or bed file\n");}
    }
  }
  elsif(@genes){ #if the gene names are inputed by user, make bedfile w/ only those genes
    foreach my $gene(@genes){ #make bed entry for each gene
	#check to make sure the gene name exists in both files and seqs length are same.
	if(($f1{$gene} && $f2{$gene})&& (length($f1{$gene})==length($f2{$gene}))){
	    $bed[$i][$j]="$gene";
	    $bed[$i][$j+1]="0\t".length($f1{$gene});
	    $i++;
	}
	else{die("$gene has no partner. Please correct the header or input options\n");}
    }
  }
  else{
    #make bed files for all the genes given there's one partner for every gene and
    #seq length are the same.
    for my $gene(keys %f1){
	if(($f1{$gene} && $f2{$gene})&& (length($f1{$gene})==length($f2{$gene}))){
	    $bed[$i][$j]="$gene";
	    $bed[$i][$j+1]="0\t".length($f1{$gene});
	    $i++;
	}
	else{die("Either there is a sequence with no partner or the sequence length is different\n");}
    }
  }
  return(\@bed);
}

#subroutine to get the sequences and it through the 
#algorithms.
sub _run_analysis{
  my $f1_ref=shift;
  my $f2_ref=shift;
  my $bed_ref=shift;
  my $file1=shift;
  my $file2=shift;
  my $opt_ref=shift;
  my %opt=%$opt_ref;

  my $p_val=0;

  my %f1=%$f1_ref;
  my %f2=%$f2_ref;
  my @bed=@$bed_ref;
  for(my $i=0; $i<scalar(@bed);$i++){
    my $gene=$bed[$i][0];
    my $co=$bed[$i][1];
    my @c=split('\t',$co);
    my ($s, $e)=($c[0],$c[1]);
    my ($seq1,$seq2);

    if(defined($f1{$gene})){
      $seq1=&_extractsubseq($s, $e, $f1{$gene});
    }
    else{die "Cannot find $gene";}

    if(defined($f2{$gene})){
      $seq2=&_extractsubseq($s, $e, $f2{$gene});
    }
    else{die "Cannot find $gene";}
    
    my ($p,$L)=&_p_value($seq1, $seq2);
    my $jc=&_jukes_cantor($seq1,$seq2,$p);
    my ($k,$ts,$tv)=&_kimura2($seq1,$seq2,$L);
    my ($f,$f1);

    if($jc eq 'NA'){ $f='%s';}else{$f='%.3f';}
    if($k eq 'NA'){ $f1='%s';}else{$f1='%.3f';}
    my $file=$file1."_".$file2;
    printf "%s\t%s\t%d\t%d\t%.3f\t$f\t$f1\n", $file,$gene,$ts,$tv,$p,$jc,$k;
   
  }
}

#Subroutine to extract the sequence according to the cordinates provided in bed file
#Written by Dr. Jonathan Flowers.
sub _extractsubseq {

  my ($s,$e,$seq) = @_;	#$s = 0-based start
			  #$e = 1-based end
  my $len = length ($seq); 
  my $sublen = $e - $s; 	#$s is zero-based, $e is 1-based
  die("error: cant extract sequence of length < 1\n")if($sublen < 1);
  die("error: -e $e is beyond end of sequence\n")if( ($s + $sublen) > $len); #correct
  my $subseq = substr $seq, $s, $sublen; #zero-indexed offset
  return($subseq);	
}

#Calculates P-distance
#number of differences/total
sub _p_value{
  my $seq1=shift;
  my $seq2=shift;
  my $d=0;
  my $L=0;
  my ($n1,$n2);
  my $len=length($seq1);

  for(my $i=0;$i<$len;$i++){
    $n1=uc substr($seq1,$i,1);
    $n2=uc substr($seq2,$i,1);
    if($n1=~ m/[^ACGT]/i || $n2=~ m/[^ACGT]/i){next;}
    else{
      if ($n1 ne $n2){ $d++;}
      $L++;
    }
  }
  return(($d/$L),$L);
}

#uses Jukes Cantor model to calculate the 
sub _jukes_cantor{
  my $seq1=shift;
  my $seq2=shift;
  my $p=shift;
  my $jc=0;

  if($p>=(3/4)){ $jc='NA';}
  else{ $jc=(-3/4)*log(1-((4/3)*($p)));}
  
  return $jc;
}

sub _kimura2{
  my $seq1=shift;
  my $seq2=shift;
  my $L=shift;
  my $p=0;
  my $q=0;
  my $d=0;
  my ($ts, $tv)=&_tstv($seq1, $seq2);

  $p=($ts/$L);
  $q=($tv/$L);
  my $w1=(1-(2*$p)-$q);
  my $w2=(1-(2*$q));
  if($w1 <= 0|| $w2 <= 0){
    $d='NA';
  }
  else{
    $d=(((-1/2)*log($w1))-(1/4*log($w2)));
  }

  return ($d,$ts,$tv);
}  

sub _tstv{
  my $seq1=shift;
  my $seq2=shift;
  my ($n1,$n2);
  my $ts=0;
  my $tv=0;

  my %vals = ('AG','S', 'CT','S', 'AT','V', 'GC','V', 'AC','V', 'GT','V',
	      'GA','S', 'TC','S', 'TA','V', 'CG','V', 'CA','V', 'TG','V');
  
  for(my $i=0;$i<length($seq1);$i++){
    $n1=uc substr($seq1,$i,1);
    $n2=uc substr($seq2,$i,1);
    
    if($n1=~ m/[^ACGT]/i || $n2=~ m/[^ACGT]/i || $n1 eq $n2){next;}
    elsif($vals{"$n1$n2"} eq 'S'){ $ts++;}
    elsif($vals{"$n1$n2"} eq 'V'){$tv++;}
    else{die("nucleotide combination $n1.$n2 is invalid");}
  }

  return ($ts, $tv);
}