use strict;
#use warnings;

# To run this program:
# more test_csv_NT_OUT.csv | perl liquid.pl > OUTFILE

# Need to remove ^M from files for running on a MAC.
# BUT, this might be a zsh problem.
# sed -e "s/^M//" filename > newfilename

print STDERR "Reading ...\n";


my $count=0;
my %hash = {};
my %hashB = {};
my %hashC = {};
my %genome = {};
my %positions = {};
my %accessions = {};
my $taskct = 0;

while (my $row2 = <>) {
	$count++;
	if($count % 500000 == 0){
                print STDERR "read $count lines:\n"
        }

	chomp $row2;
	
	my @entries = split(",", $row2);
	
#	print "$entries[0]\n";
#	print "$taskct\n";

	if (!($entries[0] eq "task_id"))
	{
		my $x = $entries[1];
		my $y = 495-$entries[2];
	
		$hash{$x}{$y}{mean} = $entries[4]+.01;
		$hash{$x}{$y}{stdev} = $entries[5]+.001;
		$hash{$x}{$y}{cv} = $entries[6];
		$hash{$x}{$y}{bg} = $entries[7];
		$hash{$x}{$y}{smb} = $hash{$x}{$y}{mean} - $hash{$x}{$y}{bg};
#		print "$x $y $hash{$x}{$y}{mean}\n";
	}
	 else{
 		$taskct++;
 		if ($taskct == 2){
 			print STDERR "HASH SAVED\n";
 			%hashB = %hash;
 			%hash = {};
 		}
 	}
	
} 
	%hashC = %hash;


my @filename = ('core1_annot_mod.tsv','core2_annot_mod.tsv') ;
$taskct = 0;

foreach my $files (@filename) {
$taskct ++;
print STDERR "$files\n";

if ($taskct == 1) {%hash = %hashB};
if ($taskct == 2) {%hash = %hashC};


open my $fh, '<', $files or die "Cannot open: $!";

my %annotation = {};
my $pos = 0;
my $x_A = 0;
my $x_T = 0;
my $x_C = 0;
my $x_G = 0;
my $y_A = 0;
my $y_T = 0;
my $y_C = 0;
my $y_G = 0;
my $ref_seq = "";
my $strain_ref_base = "N";
my $goodpos = 0;
my $counter = 0;

while (my $line = <$fh>) {
	$counter++;

	chomp $line;
	$line =~ s/\s*\z//; #remove trailing white characters (\n and others):
	my @array = split /\t/, $line;
	my $x = $array[0];
	my $y = $array[1];
	my $name = $array[5];
	my @namesplit = split("_", $name);
	my $sense = $namesplit[1];
	

	$annotation{$x}{$y}{species} =  $array[3];
	$annotation{$x}{$y}{accession} =  $array[9];
	$annotation{$x}{$y}{gene} =  $array[4];
	$annotation{$x}{$y}{base} = $array[8];
	$annotation{$x}{$y}{sense} = $sense;
	$annotation{$x}{$y}{strainpos} = $array[10]+12;
	
	
	if ($array[13] eq "TRUE") {
		if ($array[19] > 1) {
 			$goodpos = 1;
 		}else{
 			$goodpos = 0;		
 		}
		
		$pos = (13 - $array[17] + 1) + $array[20] - 1;
		$strain_ref_base = $array[8];
		
	}

	if ($array[8] eq "A") {
		$x_A = $array[0];
		$y_A = $array[1];
		$annotation{$x}{$y}{seq} = $array[6];
		my $name_A = $array[5];
		my @namesplit_A = split("_", $name_A);
		my $sense_A = $namesplit_A[1];		
	}

	if ($array[8] eq "T") {
		$x_T = $array[0];
		$y_T = $array[1];
		$annotation{$x}{$y}{seq} = $array[6];
		my $name_T = $array[5];
		my @namesplit_T = split("_", $name_T);
		my $sense_T = $namesplit_T[1];		
	}

	if ($array[8] eq "C") {
		$x_C = $array[0];
		$y_C = $array[1];
		$annotation{$x}{$y}{seq} = $array[6];
		my $name_C = $array[5];
		my @namesplit_C = split("_", $name_C);
		my $sense_C = $namesplit_C[1];		
	}

	if ($array[8] eq "G") {
		$x_G= $array[0];
		$y_G = $array[1];
		$annotation{$x}{$y}{seq} = $array[6];
		my $name_G = $array[5];
		my @namesplit_G = split("_", $name_G);
		my $sense_G = $namesplit_G[1];		
	}

	if (($annotation{$x}{$y}{species} eq "SARS-COV2") && ($array[8] eq "G")) {	
	
	if ($goodpos == 1){
	
	$positions{$pos}++; # Need to fix this after correctly getting the position calculated.
	$genome{$pos}{$array[9]}{reference} = $strain_ref_base;
	
	$accessions{$annotation{$x_A}{$y_A}{accession}}++;
	
#	print "$pos\t$positions{$pos}\t$annotation{$x_A}{$y_A}{species}\t$annotation{$x_A}{$y_A}{accession}\n";
# 	print "$annotation{$x_A}{$y_A}{seq}\n";
# 	print "$annotation{$x_T}{$y_T}{seq}\n";
# 	print "$annotation{$x_C}{$y_C}{seq}\n";
# 	print "$annotation{$x_G}{$y_G}{seq}\n";
		
	$annotation{$x}{$y}{position} = $pos;
	
	if ($annotation{$x}{$y}{sense} eq "S") { 
		if (exists $genome{$pos}{$array[9]}{A}){
			if (exists $genome{$pos}{$array[9]}{A2}) {print STDERR "OVERWRITING $taskct $annotation{$x_A}{$y_A}{accession} Sense $x_A,$y_A $pos ...\n";}
				$genome{$pos}{$array[9]}{A2} = $hash{$x_A}{$y_A}{mean};
				$genome{$pos}{$array[9]}{A2bg} = $hash{$x_A}{$y_A}{bg};
				$genome{$pos}{$array[9]}{A2smb} = $hash{$x_A}{$y_A}{smb};
				$genome{$pos}{$array[9]}{A2stdev} = $hash{$x_A}{$y_A}{stdev};
				$genome{$pos}{$array[9]}{A2cv} = $hash{$x_A}{$y_A}{cv};
			} else {
				$genome{$pos}{$array[9]}{A} = $hash{$x_A}{$y_A}{mean};
				$genome{$pos}{$array[9]}{Abg} = $hash{$x_A}{$y_A}{bg};
				$genome{$pos}{$array[9]}{Asmb} = $hash{$x_A}{$y_A}{smb};
				$genome{$pos}{$array[9]}{Astdev} = $hash{$x_A}{$y_A}{stdev};
				$genome{$pos}{$array[9]}{Acv} = $hash{$x_A}{$y_A}{cv};
			}		
			
			if (exists $genome{$pos}{$array[9]}{T}){
				$genome{$pos}{$array[9]}{T2} = $hash{$x_T}{$y_T}{mean};
				$genome{$pos}{$array[9]}{T2bg} = $hash{$x_T}{$y_T}{bg};
				$genome{$pos}{$array[9]}{T2smb} = $hash{$x_T}{$y_T}{smb};
				$genome{$pos}{$array[9]}{T2stdev} = $hash{$x_T}{$y_T}{stdev};
				$genome{$pos}{$array[9]}{T2cv} = $hash{$x_T}{$y_T}{cv};
			} else {
				$genome{$pos}{$array[9]}{T} = $hash{$x_T}{$y_T}{mean};
				$genome{$pos}{$array[9]}{Tbg} = $hash{$x_T}{$y_T}{bg};
				$genome{$pos}{$array[9]}{Tsmb} = $hash{$x_T}{$y_T}{smb};
				$genome{$pos}{$array[9]}{Tstdev} = $hash{$x_T}{$y_T}{stdev};
				$genome{$pos}{$array[9]}{Tcv} = $hash{$x_T}{$y_T}{cv};
			}

			if (exists $genome{$pos}{$array[9]}{C}){
				$genome{$pos}{$array[9]}{C2} = $hash{$x_C}{$y_C}{mean};
				$genome{$pos}{$array[9]}{C2bg} = $hash{$x_C}{$y_C}{bg};
				$genome{$pos}{$array[9]}{C2smb} = $hash{$x_C}{$y_C}{smb};
				$genome{$pos}{$array[9]}{C2stdev} = $hash{$x_C}{$y_C}{stdev};
				$genome{$pos}{$array[9]}{C2cv} = $hash{$x_C}{$y_C}{cv};
			} else {
				$genome{$pos}{$array[9]}{C} = $hash{$x_C}{$y_C}{mean};
				$genome{$pos}{$array[9]}{Cbg} = $hash{$x_C}{$y_C}{bg};
				$genome{$pos}{$array[9]}{Csmb} = $hash{$x_C}{$y_C}{smb};
				$genome{$pos}{$array[9]}{Cstdev} = $hash{$x_C}{$y_C}{stdev};
				$genome{$pos}{$array[9]}{Ccv} = $hash{$x_C}{$y_C}{cv};
			}

			if (exists $genome{$pos}{$array[9]}{G}){
				$genome{$pos}{$array[9]}{G2} = $hash{$x_G}{$y_G}{mean};
				$genome{$pos}{$array[9]}{G2bg} = $hash{$x_G}{$y_G}{bg};
				$genome{$pos}{$array[9]}{G2smb} = $hash{$x_G}{$y_G}{smb};
				$genome{$pos}{$array[9]}{G2stdev} = $hash{$x_G}{$y_G}{stdev};
				$genome{$pos}{$array[9]}{G2cv} = $hash{$x_G}{$y_G}{cv};
			} else {
				$genome{$pos}{$array[9]}{G} = $hash{$x_G}{$y_G}{mean};
				$genome{$pos}{$array[9]}{Gbg} = $hash{$x_G}{$y_G}{bg};
				$genome{$pos}{$array[9]}{Gsmb} = $hash{$x_G}{$y_G}{smb};
				$genome{$pos}{$array[9]}{Gstdev} = $hash{$x_G}{$y_G}{stdev};
				$genome{$pos}{$array[9]}{Gcv} = $hash{$x_G}{$y_G}{cv};
			}
	}
	if ($annotation{$x}{$y}{sense} eq "AS") { 
		if (exists $genome{$pos}{$annotation{$x_A}{$y_A}{accession}}{AS_A}){
			if (exists $genome{$pos}{$annotation{$x_A}{$y_A}{accession}}{AS_A2}) {print STDERR "OVERWRITING $taskct $annotation{$x_A}{$y_A}{accession} AntiSense $x_A,$y_A $pos ...\n";}
				$genome{$pos}{$array[9]}{AS_A2} = $hash{$x_A}{$y_A}{mean};
				$genome{$pos}{$array[9]}{AS_A2bg} = $hash{$x_A}{$y_A}{bg};
				$genome{$pos}{$array[9]}{AS_A2smb} = $hash{$x_A}{$y_A}{smb};
				$genome{$pos}{$array[9]}{AS_A2stdev} = $hash{$x_A}{$y_A}{stdev};
				$genome{$pos}{$array[9]}{AS_A2cv} = $hash{$x_A}{$y_A}{cv};
			} else {
				$genome{$pos}{$array[9]}{AS_A} = $hash{$x_A}{$y_A}{mean};
				$genome{$pos}{$array[9]}{AS_Abg} = $hash{$x_A}{$y_A}{bg};
				$genome{$pos}{$array[9]}{AS_Asmb} = $hash{$x_A}{$y_A}{smb};
				$genome{$pos}{$array[9]}{AS_Astdev} = $hash{$x_A}{$y_A}{stdev};
				$genome{$pos}{$array[9]}{AS_Acv} = $hash{$x_A}{$y_A}{cv};
			}		
			
			if (exists $genome{$pos}{$array[9]}{AS_T}){
				$genome{$pos}{$array[9]}{AS_T2} = $hash{$x_T}{$y_T}{mean};
				$genome{$pos}{$array[9]}{AS_T2bg} = $hash{$x_T}{$y_T}{bg};
				$genome{$pos}{$array[9]}{AS_T2smb} = $hash{$x_T}{$y_T}{smb};
				$genome{$pos}{$array[9]}{AS_T2stdev} = $hash{$x_T}{$y_T}{stdev};
				$genome{$pos}{$array[9]}{AS_T2cv} = $hash{$x_T}{$y_T}{cv};
			} else {
				$genome{$pos}{$array[9]}{AS_T} = $hash{$x_T}{$y_T}{mean};
				$genome{$pos}{$array[9]}{AS_Tbg} = $hash{$x_T}{$y_T}{bg};
				$genome{$pos}{$array[9]}{AS_Tsmb} = $hash{$x_T}{$y_T}{smb};
				$genome{$pos}{$array[9]}{AS_Tstdev} = $hash{$x_T}{$y_T}{stdev};
				$genome{$pos}{$array[9]}{AS_Tcv} = $hash{$x_T}{$y_T}{cv};
			}

			if (exists $genome{$pos}{$annotation{$x_C}{$y_C}{accession}}{AS_C}){
				$genome{$pos}{$array[9]}{AS_C2} = $hash{$x_C}{$y_C}{mean};
				$genome{$pos}{$array[9]}{AS_C2bg} = $hash{$x_C}{$y_C}{bg};
				$genome{$pos}{$array[9]}{AS_C2smb} = $hash{$x_C}{$y_C}{smb};
				$genome{$pos}{$array[9]}{AS_C2stdev} = $hash{$x_C}{$y_C}{stdev};
				$genome{$pos}{$array[9]}{AS_C2cv} = $hash{$x_C}{$y_C}{cv};
			} else {
				$genome{$pos}{$array[9]}{AS_C} = $hash{$x_C}{$y_C}{mean};
				$genome{$pos}{$array[9]}{AS_Cbg} = $hash{$x_C}{$y_C}{bg};
				$genome{$pos}{$array[9]}{AS_Csmb} = $hash{$x_C}{$y_C}{smb};
				$genome{$pos}{$array[9]}{AS_Cstdev} = $hash{$x_C}{$y_C}{stdev};
				$genome{$pos}{$array[9]}{AS_Ccv} = $hash{$x_C}{$y_C}{cv};
			}

			if (exists $genome{$pos}{$annotation{$x_G}{$y_G}{accession}}{AS_G}){
				$genome{$pos}{$array[9]}{AS_G2} = $hash{$x_G}{$y_G}{mean};
				$genome{$pos}{$array[9]}{AS_G2bg} = $hash{$x_G}{$y_G}{bg};
				$genome{$pos}{$array[9]}{AS_G2smb} = $hash{$x_G}{$y_G}{smb};
				$genome{$pos}{$array[9]}{AS_G2stdev} = $hash{$x_G}{$y_G}{stdev};
				$genome{$pos}{$array[9]}{AS_G2cv} = $hash{$x_G}{$y_G}{cv};
			} else {
				$genome{$pos}{$array[9]}{AS_G} = $hash{$x_G}{$y_G}{mean};
				$genome{$pos}{$array[9]}{AS_Gbg} = $hash{$x_G}{$y_G}{bg};
				$genome{$pos}{$array[9]}{AS_Gsmb} = $hash{$x_G}{$y_G}{smb};
				$genome{$pos}{$array[9]}{AS_Gstdev} = $hash{$x_G}{$y_G}{stdev};
				$genome{$pos}{$array[9]}{AS_Gcv} = $hash{$x_G}{$y_G}{cv};
			}
	}
#	print "$x,$y $pos $annotation{$x}{$y}{base} $annotation{$x}{$y}{sense} $genome{$pos}{$array[9]}{reference} $annotation{$x}{$y}{isrefbase} *** $genome{$pos}{$array[9]}{A},$genome{$pos}{$array[9]}{Asmb} $genome{$pos}{$array[9]}{T},$genome{$pos}{$array[9]}{Tsmb} $genome{$pos}{$array[9]}{C},$genome{$pos}{$array[9]}{Csmb} $genome{$pos}{$array[9]}{G},$genome{$pos}{$array[9]}{Gsmb} *** $genome{$pos}{$array[9]}{A2},$genome{$pos}{$array[9]}{A2smb} $genome{$pos}{$array[9]}{T2},$genome{$pos}{$array[9]}{T2smb} $genome{$pos}{$array[9]}{C2},$genome{$pos}{$array[9]}{C2smb} $genome{$pos}{$array[9]}{G2},$genome{$pos}{$array[9]}{G2smb} *** $genome{$pos}{$array[9]}{AS_A},$genome{$pos}{$array[9]}{AS_Asmb} $genome{$pos}{$array[9]}{AS_T},$genome{$pos}{$array[9]}{AS_Tsmb} $genome{$pos}{$array[9]}{AS_C},$genome{$pos}{$array[9]}{AS_Csmb} $genome{$pos}{$array[9]}{AS_G},$genome{$pos}{$array[9]}{AS_Gsmb}\n";
	}
	$pos = 0;
	$goodpos = 0;
	}
}
close $fh;
}

 # Print json file
open my $fo, '>', 'core1_outstats.tsv' or die "Cannot open: $!";

 	my $ct=0;
 	my $ctt = 0;
 	
 	my @strains = keys %accessions;
 	my @allpositions = sort keys %positions;
 	my @sortedpositions = sort {$a <=> $b} @allpositions;
 	
 	print "[\n";
	foreach my $printpos (@sortedpositions) {
 	foreach my $printsacc (@strains) {
		$ct++;
		if ($ct == 1){
			print $fo "Accession\tPosition\tRefBase\tSenseCall\tASenseCall\t";
			print $fo "SenseDifference\tSenseDifferetial\tASenseDifference\tASenseDifferetial\tSenseHighestSignal\tASenseHighestSignal\tSenseLowestSignal\tASenseLowestSignal\n";
		}
		
		if (exists $genome{$printpos}{$printsacc}{A}){
		$ctt++;
		
		if ($ct > 1){
			my @bases = ($genome{$printpos}{$printsacc}{A},$genome{$printpos}{$printsacc}{T},$genome{$printpos}{$printsacc}{C},$genome{$printpos}{$printsacc}{G});
			my @bases2 = ($genome{$printpos}{$printsacc}{A2},$genome{$printpos}{$printsacc}{T2},$genome{$printpos}{$printsacc}{C2},$genome{$printpos}{$printsacc}{G2});
			my @asbases = ($genome{$printpos}{$printsacc}{AS_A},$genome{$printpos}{$printsacc}{AS_T},$genome{$printpos}{$printsacc}{AS_C},$genome{$printpos}{$printsacc}{AS_G});
			my @asbases2 = ($genome{$printpos}{$printsacc}{AS_A2},$genome{$printpos}{$printsacc}{AS_T2},$genome{$printpos}{$printsacc}{AS_C2},$genome{$printpos}{$printsacc}{AS_G2});

			my @sortedbases = sort { $a <=> $b } @bases;
			my @sortedbases2 = sort { $a <=> $b } @bases2;
			my @sortedasbases = sort { $a <=> $b } @asbases;
			my @sortedasbases2 = sort { $a <=> $b } @asbases2;

			my $sdifference = $sortedbases[3]-$sortedbases[0];
			my $sdifference2 = $sortedbases2[3]-$sortedbases2[0];
			my $asdifference = $sortedasbases[3]-$sortedasbases[0];
			my $asdifference2 = $sortedasbases2[3]-$sortedasbases2[0];


			my $sdifferetial = 0;
			my $sdifferetial2 = 0;
			my $asdifferetial = 0;
			my $asdifferetial2 = 0;			
			if ($sdifference > 0) { $sdifferetial = ($sortedbases[3]-$sortedbases[2])/$sdifference;}
			if ($sdifference2 > 0) { $sdifferetial2 = ($sortedbases2[3]-$sortedbases2[2])/$sdifference2;}
			if ($asdifference > 0) { $asdifferetial = ($sortedasbases[3]-$sortedasbases[2])/$asdifference;}
			if ($asdifference2 > 0) { $asdifferetial2 = ($sortedasbases2[3]-$sortedasbases2[2])/$asdifference2;}
			$genome{$printpos}{SenseCall} = "N";
			$genome{$printpos}{ASenseCall} = "N";
			if (($genome{$printpos}{$printsacc}{A} > $genome{$printpos}{$printsacc}{T}) && ($genome{$printpos}{$printsacc}{A} > $genome{$printpos}{$printsacc}{C}) && ($genome{$printpos}{$printsacc}{A} > $genome{$printpos}{$printsacc}{G})){
						$genome{$printpos}{$printsacc}{SenseCall} = "A";	
			}
			if (($genome{$printpos}{$printsacc}{T} > $genome{$printpos}{$printsacc}{A}) && ($genome{$printpos}{$printsacc}{T} > $genome{$printpos}{$printsacc}{C}) && ($genome{$printpos}{$printsacc}{T} > $genome{$printpos}{$printsacc}{G})){
						$genome{$printpos}{$printsacc}{SenseCall} = "T";	
			}
			if (($genome{$printpos}{$printsacc}{C} > $genome{$printpos}{$printsacc}{A}) && ($genome{$printpos}{$printsacc}{C} > $genome{$printpos}{$printsacc}{T}) && ($genome{$printpos}{$printsacc}{C} > $genome{$printpos}{$printsacc}{G})){
						$genome{$printpos}{$printsacc}{SenseCall} = "C";	
			}
			if (($genome{$printpos}{$printsacc}{G} > $genome{$printpos}{$printsacc}{A}) && ($genome{$printpos}{$printsacc}{G} > $genome{$printpos}{$printsacc}{T}) && ($genome{$printpos}{$printsacc}{G} > $genome{$printpos}{$printsacc}{C})){
						$genome{$printpos}{$printsacc}{SenseCall} = "G";	
			}
			if (($genome{$printpos}{$printsacc}{AS_A} > $genome{$printpos}{$printsacc}{AS_T}) && ($genome{$printpos}{$printsacc}{AS_A} > $genome{$printpos}{$printsacc}{AS_C}) && ($genome{$printpos}{$printsacc}{AS_A} > $genome{$printpos}{$printsacc}{AS_G})){
						$genome{$printpos}{$printsacc}{ASenseCall} = "A";	
			}
			if (($genome{$printpos}{$printsacc}{AS_T} > $genome{$printpos}{$printsacc}{AS_A}) && ($genome{$printpos}{$printsacc}{AS_T} > $genome{$printpos}{$printsacc}{AS_C}) && ($genome{$printpos}{$printsacc}{AS_T} > $genome{$printpos}{$printsacc}{AS_G})){
						$genome{$printpos}{$printsacc}{ASenseCall} = "T";	
			}
			if (($genome{$printpos}{$printsacc}{AS_C} > $genome{$printpos}{$printsacc}{AS_A}) && ($genome{$printpos}{$printsacc}{AS_C} > $genome{$printpos}{$printsacc}{AS_T}) && ($genome{$printpos}{$printsacc}{AS_C} > $genome{$printpos}{$printsacc}{AS_G})){
						$genome{$printpos}{$printsacc}{ASenseCall} = "C";	
			}
			if (($genome{$printpos}{$printsacc}{AS_G} > $genome{$printpos}{$printsacc}{AS_A}) && ($genome{$printpos}{$printsacc}{AS_G} > $genome{$printpos}{$printsacc}{AS_T}) && ($genome{$printpos}{$printsacc}{AS_G} > $genome{$printpos}{$printsacc}{AS_C})){
						$genome{$printpos}{$printsacc}{ASenseCall} = "G";	
			}
			
			if ($ctt == 1) {
				print "\t{\n";
			} else {
				print "\t,{\n";			
			}			

			print "\t\t\"pos\": $printpos,\n";
			print "\t\t\"number\": $ct,\n";			
			print "\t\t\"number2\": $ctt,\n";			
			print "\t\t\"accession\": \"$printsacc\",\n";			
			print "\t\t\"Ref_base\": \"$genome{$printpos}{$printsacc}{reference}\",\n";
			print "\t\t\"SenseCall_base\": \"$genome{$printpos}{$printsacc}{SenseCall}\",\n";
			print "\t\t\"ASenseCall_base\": \"$genome{$printpos}{$printsacc}{ASenseCall}\",\n";
			if ($genome{$printpos}{$printsacc}{ASenseCall} eq $genome{$printpos}{$printsacc}{SenseCall}){
				print "\t\t\"Call_same\": \"TRUE\",\n";			
			}else{
				print "\t\t\"Call_same\": \"FALSE\",\n";			
			}
			if ($genome{$printpos}{$printsacc}{SenseCall} eq $genome{$printpos}{$printsacc}{reference}){
				print "\t\t\"SCall_ref\": \"TRUE\",\n";			
			}else{
				print "\t\t\"SCall_ref\": \"FALSE\",\n";			
			}
			if ($genome{$printpos}{$printsacc}{ASenseCall} eq $genome{$printpos}{$printsacc}{reference}){
				print "\t\t\"ASCall_ref\": \"TRUE\",\n";			
			}else{
				print "\t\t\"ASCall_ref\": \"FALSE\",\n";			
			}
			
			print $fo "$printsacc\t$printpos\t$genome{$printpos}{$printsacc}{reference}\t$genome{$printpos}{$printsacc}{SenseCall}\t$genome{$printpos}{$printsacc}{ASenseCall}\t";
			print $fo "$sdifference\t$sdifferetial\t$asdifference\t$asdifferetial\t$sortedbases[3]\t$sortedasbases[3]\t$sortedbases[0]\t$sortedasbases[0]\n";

			
			if ($asdifference > 0) {
				print "\t\t\"AS_Difference\": $asdifference,\n";
				print "\t\t\"AS_Differential\": $asdifferetial,\n";
			}
			if ($sdifference > 0) {
				print "\t\t\"S_Difference\": $sdifference,\n";
				print "\t\t\"S_Differential\": $sdifferetial,\n";
			}
			if ($asdifference2 > 0) {
				print "\t\t\"AS_Difference\": $asdifference2,\n";
				print "\t\t\"AS_Differential\": $asdifferetial2,\n";
			}
			if ($sdifference2 > 0) {
				print "\t\t\"S_Difference\": $sdifference2,\n";
				print "\t\t\"S_Differential\": $sdifferetial2,\n";
			}
			
			if ($asdifference > 0) {
				print "\t\t\"antisense\": [\n";
				print "\t\t\t{\n";
				print "\t\t\t\t\"base\": \"A\",\n";
				print "\t\t\t\t\"intensity\": $genome{$printpos}{$printsacc}{AS_A},\n"; 
				print "\t\t\t\t\"smb\": $genome{$printpos}{$printsacc}{AS_Asmb},\n";
				print "\t\t\t\t\"stdev\": $genome{$printpos}{$printsacc}{AS_Astdev}\n";
#				print "\t\t\t\t\"cv\": $genome{$printpos}{$printsacc}{AS_Acv}\n";
				print "\t\t\t},\n";
				
				print "\t\t\t{\n";
				print "\t\t\t\t\"base\": \"T\",\n";
				print "\t\t\t\t\"intensity\": $genome{$printpos}{$printsacc}{AS_T},\n"; 
				print "\t\t\t\t\"smb\": $genome{$printpos}{$printsacc}{AS_Tsmb},\n";
				print "\t\t\t\t\"stdev\": $genome{$printpos}{$printsacc}{AS_Tstdev}\n";
#				print "\t\t\t\t\"cv\": $genome{$printpos}{$printsacc}{AS_Tcv}\n";
				print "\t\t\t},\n";
				
				print "\t\t\t{\n";
				print "\t\t\t\t\"base\": \"C\",\n";
				print "\t\t\t\t\"intensity\": $genome{$printpos}{$printsacc}{AS_C},\n"; 
				print "\t\t\t\t\"smb\": $genome{$printpos}{$printsacc}{AS_Csmb},\n";
				print "\t\t\t\t\"stdev\": $genome{$printpos}{$printsacc}{AS_Cstdev}\n";
#				print "\t\t\t\t\"cv\": $genome{$printpos}{$printsacc}{AS_Ccv}\n";
				print "\t\t\t},\n";
				
				print "\t\t\t{\n";
				print "\t\t\t\t\"base\": \"G\",\n";
				print "\t\t\t\t\"intensity\": $genome{$printpos}{$printsacc}{AS_G},\n"; 
				print "\t\t\t\t\"smb\": $genome{$printpos}{$printsacc}{AS_Gsmb},\n";
				print "\t\t\t\t\"stdev\": $genome{$printpos}{$printsacc}{AS_Gstdev}\n";
#				print "\t\t\t\t\"cv\": $genome{$printpos}{$printsacc}{AS_Gcv}\n";
				print "\t\t\t}";
				
				if ($asdifference2 > 0) {
					print ",\n";

					print "\t\t\t{\n";
					print "\t\t\t\t\"base\": \"A\",\n";
					print "\t\t\t\t\"intensity\": $genome{$printpos}{$printsacc}{AS_A2},\n"; 
					print "\t\t\t\t\"smb\": $genome{$printpos}{$printsacc}{AS_A2smb},\n";
					print "\t\t\t\t\"stdev\": $genome{$printpos}{$printsacc}{AS_A2stdev}\n";
#					print "\t\t\t\t\"cv\": $genome{$printpos}{$printsacc}{AS_A2cv}\n";
					print "\t\t\t},\n";
				
					print "\t\t\t{\n";
					print "\t\t\t\t\"base\": \"T\",\n";
					print "\t\t\t\t\"intensity\": $genome{$printpos}{$printsacc}{AS_T2},\n"; 
					print "\t\t\t\t\"smb\": $genome{$printpos}{$printsacc}{AS_T2smb},\n";
					print "\t\t\t\t\"stdev\": $genome{$printpos}{$printsacc}{AS_T2stdev}\n";
#					print "\t\t\t\t\"cv\": $genome{$printpos}{$printsacc}{AS_T2cv}\n";
					print "\t\t\t},\n";
				
					print "\t\t\t{\n";
					print "\t\t\t\t\"base\": \"C\",\n";
					print "\t\t\t\t\"intensity\": $genome{$printpos}{$printsacc}{AS_C2},\n"; 
					print "\t\t\t\t\"smb\": $genome{$printpos}{$printsacc}{AS_C2smb},\n";
					print "\t\t\t\t\"stdev\": $genome{$printpos}{$printsacc}{AS_C2stdev}\n";
#					print "\t\t\t\t\"cv\": $genome{$printpos}{$printsacc}{AS_C2cv}\n";
					print "\t\t\t},\n";
				
					print "\t\t\t{\n";
					print "\t\t\t\t\"base\": \"G\",\n";
					print "\t\t\t\t\"intensity\": $genome{$printpos}{$printsacc}{AS_G2},\n"; 
					print "\t\t\t\t\"smb\": $genome{$printpos}{$printsacc}{AS_G2smb},\n";
					print "\t\t\t\t\"stdev\": $genome{$printpos}{$printsacc}{AS_G2stdev}\n";
#					print "\t\t\t\t\"cv\": $genome{$printpos}{$printsacc}{AS_G2cv}\n";
					print "\t\t\t}";
				}
				print "\n\t\t],\n";

			}
	
			if ($sdifference > -.1) {
				print "\t\t\"sense\": [\n";
				print "\t\t\t{\n";
				print "\t\t\t\t\"base\": \"A\",\n";
				print "\t\t\t\t\"intensity\": $genome{$printpos}{$printsacc}{A},\n"; 
				print "\t\t\t\t\"smb\": $genome{$printpos}{$printsacc}{Asmb},\n";
				print "\t\t\t\t\"stdev\": $genome{$printpos}{$printsacc}{Astdev}\n";
#				print "\t\t\t\t\"cv\": $genome{$printpos}{$printsacc}{Acv}\n";
				
				print "\t\t\t},\n";
				
				print "\t\t\t{\n";
				print "\t\t\t\t\"base\": \"T\",\n";
				print "\t\t\t\t\"intensity\": $genome{$printpos}{$printsacc}{T},\n"; 
				print "\t\t\t\t\"smb\": $genome{$printpos}{$printsacc}{Tsmb},\n";
				print "\t\t\t\t\"stdev\": $genome{$printpos}{$printsacc}{Tstdev}\n";
#				print "\t\t\t\t\"cv\": $genome{$printpos}{$printsacc}{Tcv}\n";
				print "\t\t\t},\n";
				
				print "\t\t\t{\n";
				print "\t\t\t\t\"base\": \"C\",\n";
				print "\t\t\t\t\"intensity\": $genome{$printpos}{$printsacc}{C},\n"; 
				print "\t\t\t\t\"smb\": $genome{$printpos}{$printsacc}{Csmb},\n";
				print "\t\t\t\t\"stdev\": $genome{$printpos}{$printsacc}{Cstdev}\n";
#				print "\t\t\t\t\"cv\": $genome{$printpos}{$printsacc}{Ccv}\n";
				print "\t\t\t},\n";
				
				print "\t\t\t{\n";
				print "\t\t\t\t\"base\": \"G\",\n";
				print "\t\t\t\t\"intensity\": $genome{$printpos}{$printsacc}{G},\n"; 
				print "\t\t\t\t\"smb\": $genome{$printpos}{$printsacc}{Gsmb},\n";
				print "\t\t\t\t\"stdev\": $genome{$printpos}{$printsacc}{Gstdev}\n";
#				print "\t\t\t\t\"cv\": $genome{$printpos}{$printsacc}{Gcv}\n";
				print "\t\t\t}";
				
				if ($sdifference2 > 0) {
					print ",\n";

					print "\t\t\t{\n";
					print "\t\t\t\t\"base\": \"A\",\n";
					print "\t\t\t\t\"intensity\": $genome{$printpos}{$printsacc}{A2},\n"; 
					print "\t\t\t\t\"smb\": $genome{$printpos}{$printsacc}{A2smb},\n";
					print "\t\t\t\t\"stdev\": $genome{$printpos}{$printsacc}{A2stdev}\n";
#					print "\t\t\t\t\"cv\": $genome{$printpos}{$printsacc}{A2cv}\n";
					print "\t\t\t},\n";
				
					print "\t\t\t{\n";
					print "\t\t\t\t\"base\": \"T\",\n";
					print "\t\t\t\t\"intensity\": $genome{$printpos}{$printsacc}{T2},\n"; 
					print "\t\t\t\t\"smb\": $genome{$printpos}{$printsacc}{T2smb},\n";
					print "\t\t\t\t\"stdev\": $genome{$printpos}{$printsacc}{T2stdev}\n";
#					print "\t\t\t\t\"cv\": $genome{$printpos}{$printsacc}{T2cv}\n";
					print "\t\t\t},\n";
				
					print "\t\t\t{\n";
					print "\t\t\t\t\"base\": \"C\",\n";
					print "\t\t\t\t\"intensity\": $genome{$printpos}{$printsacc}{C2},\n"; 
					print "\t\t\t\t\"smb\": $genome{$printpos}{$printsacc}{C2smb},\n";
					print "\t\t\t\t\"stdev\": $genome{$printpos}{$printsacc}{C2stdev}\n";
#					print "\t\t\t\t\"cv\": $genome{$printpos}{$printsacc}{C2cv}\n";
					print "\t\t\t},\n";
				
					print "\t\t\t{\n";
					print "\t\t\t\t\"base\": \"G\",\n";
					print "\t\t\t\t\"intensity\": $genome{$printpos}{$printsacc}{G2},\n"; 
					print "\t\t\t\t\"smb\": $genome{$printpos}{$printsacc}{G2smb},\n";
					print "\t\t\t\t\"stdev\": $genome{$printpos}{$printsacc}{G2stdev}\n";
#					print "\t\t\t\t\"cv\": $genome{$printpos}{$printsacc}{G2cv}\n";
					print "\t\t\t}";
				}
				print "\n\t\t]\n";
			}
		print "\t}\n";						
		}
		}
	}
	}
	print "]\n";
close $fo;





