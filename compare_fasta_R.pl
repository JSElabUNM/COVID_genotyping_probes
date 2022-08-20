#!/usr/bin/perl -w
# This script reads several sequences and computes the relative content of G+C of each sequence.

use strict; 

# my $infile = "NC_045512.2.fasta";                               # This is the file path
# open INFILE, $infile or die "Can't open $infile: $!";        # This opens file, but if file isn't there it mentions this will not open
# my $outfile = "Lab1_SeqOutput.txt";             # This is the file's output
# open OUTFILE, ">$outfile" or die "Cannot open $outfile: $!"; # This opens the output file, otherwise it mentions this will not open
# 
# my $sequence = ();  # This sequence variable stores the sequences from the .fasta file
# my $GC = 0;         # This variable checks for G + C content
# 
# my $line;                             # This reads the input file one-line-at-a-time
# 
# while ($line = <INFILE>) {
#     chomp $line;                      # This removes "\n" at the end of each line (this is invisible)
# 
#     if($line =~ /^\s*$/) {         # This finds lines with whitespaces from the beginning to the ending of the sequence. Removes blank line.
#         next;
# 
#     } elsif($line =~ /^\s*#/) {        # This finds lines with spaces before the hash character. Removes .fasta comment
#         next; 
#     } elsif($line =~ /^>/) {           # This finds lines with the '>' symbol at beginning of label. Removes .fasta label
#         next;
#     } else {
#         $sequence = $line;
#     }
# 
#     $sequence =~ s/\s//g;               # Whitespace characters are removed
#     print OUTFILE $sequence;
# }

use Bio::SeqIO;
my $seqio = Bio::SeqIO->new(-file => "NC_045512.2.fasta", '-format' => 'Fasta');
my $ref = "";
my $seqio2 = Bio::SeqIO->new(-file => "in.fa", '-format' => 'Fasta');

while(my $seq = $seqio->next_seq) {
  $ref = $seq->seq;
    # do stuff with $string
}


print "Sample Name\tMatch\tLQ Mat\tNoCalls\tHQ Var\tLQ Var\tCoverage           \tHigh Qual Accuracy\tLow Qual Accuracy\n";

while(my $seq = $seqio2->next_seq) {
  my $new = $seq->seq;
  my $name = $seq->id;
  
  my $match = 0;
  my $lqmatch = 0;
  my $nocall = 0;
  my $hqwrong = 0;
  my $lqwrong = 0;
  
  my $shift = 25;

  for (my $i = 0; $i < length($new)-$shift;$i = $i + 1){
	  my $refbase = substr($ref,$i+$shift,1);
	  my $newbase = substr($new,$i,1);
	  
	  if ($newbase eq "N"){
	  	$nocall++;
	  } else {	  
	     if ($refbase eq $newbase){
		     $match++;
	     } elsif (lc($refbase) eq $newbase){
	  	     $lqmatch++;
	     } else {
	        if ($newbase =~ /[A-Z]/) {
               $hqwrong++;
            }else{	     
	     	   $lqwrong++;
	     	}
	     }
	  }	  
  }

  my $cov = (length($new) - $nocall)/(length($new));
  my $hqaccuracy = 1 - ($hqwrong/($match+$hqwrong));
  my $lqaccuracy = 1;
  
  print "$name\t$match\t$lqmatch\t$nocall\t$hqwrong\t$lqwrong\t$cov\t$hqaccuracy\t$lqaccuracy\n";
}

