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
my %annotation = {};
my %genome = {};
my %positions = {};

while (my $row2 = <>) {
	$count++;
	if($count % 500000 == 0){
                print STDERR "read $count lines:\n"
        }

	chomp $row2;
	
	my @entries = split(",", $row2);
	my $x = $entries[1];
	my $y = 495-$entries[2];
	
	$hash{$x}{$y}{mean} = $entries[4]+.01;
	$hash{$x}{$y}{stdev} = $entries[5]+.001;
	$hash{$x}{$y}{cv} = $entries[6]+.0001;
	$hash{$x}{$y}{bg} = $entries[7]+.0001;
	$hash{$x}{$y}{smb} = $hash{$x}{$y}{mean} - $hash{$x}{$y}{bg};
	#print "$x $y $hash{$x}{$y}{mean}\n";
} 

open my $fh, '<', 'core0_annot_mod.tsv' or die "Cannot open: $!";

while (my $line = <$fh>) {
	chomp $line;
	$line =~ s/\s*\z//; #remove trailing white characters (\n and others):
	my @array = split /\t/, $line;
	my $x = $array[0];
	my $y = $array[1];
	my $start = $array[10];
	my $stop = $array[11];
	my $pos = ($stop - $start)/2 +  $start;
	my $name = $array[5];
	my @namesplit = split("_", $name);
	my $sense = $namesplit[1];
	

	$annotation{$x}{$y}{species} =  $array[3];
	$annotation{$x}{$y}{accession} =  $array[9];
	$annotation{$x}{$y}{gene} =  $array[4];
	$annotation{$x}{$y}{base} = $array[8];
	$annotation{$x}{$y}{sense} = $sense;


	if (($annotation{$x}{$y}{species} eq "SARS-COV2") && ($annotation{$x}{$y}{accession} eq "EPI_ISL_402125")) {
	
	$positions{$pos}++;
		
	if ($array[13] eq "TRUE") {
		$annotation{$x}{$y}{isrefbase} = 1;
		$genome{$pos}{reference} = $annotation{$x}{$y}{base};
	}
	if ($array[13] eq "FALSE") {$annotation{$x}{$y}{isrefbase} = 0;}
	$annotation{$x}{$y}{position} = $pos;
	
	if ($annotation{$x}{$y}{sense} eq "S") { 
		if ($annotation{$x}{$y}{base} eq "A") {
			if (exists $genome{$pos}{A}){
				if (exists $genome{$pos}{A2}) {print STDERR "OVERWRITING Sense $x,$y $pos ($start $stop - $genome{$pos}{A2}) ...\n";}
				$genome{$pos}{A2} = $hash{$x}{$y}{mean};
				$genome{$pos}{A2stdev} = $hash{$x}{$y}{stdev};
				$genome{$pos}{A2cv} = $hash{$x}{$y}{cv};				
				$genome{$pos}{A2bg} = $hash{$x}{$y}{bg};
				$genome{$pos}{A2smb} = $hash{$x}{$y}{smb};
			} else {
				$genome{$pos}{A} = $hash{$x}{$y}{mean};
				$genome{$pos}{Astdev} = $hash{$x}{$y}{stdev};
				$genome{$pos}{Acv} = $hash{$x}{$y}{cv};				
				$genome{$pos}{Abg} = $hash{$x}{$y}{bg};
				$genome{$pos}{Asmb} = $hash{$x}{$y}{smb};
			}
		}
		if ($annotation{$x}{$y}{base} eq "T") {
			if (exists $genome{$pos}{T}){
				$genome{$pos}{T2} = $hash{$x}{$y}{mean};
				$genome{$pos}{T2stdev} = $hash{$x}{$y}{stdev};
				$genome{$pos}{T2cv} = $hash{$x}{$y}{cv};				
				$genome{$pos}{T2bg} = $hash{$x}{$y}{bg};
				$genome{$pos}{T2smb} = $hash{$x}{$y}{smb};
			} else {
				$genome{$pos}{T} = $hash{$x}{$y}{mean};
				$genome{$pos}{Tstdev} = $hash{$x}{$y}{stdev};
				$genome{$pos}{Tcv} = $hash{$x}{$y}{cv};				
				$genome{$pos}{Tbg} = $hash{$x}{$y}{bg};
				$genome{$pos}{Tsmb} = $hash{$x}{$y}{smb};
			}
		}
		if ($annotation{$x}{$y}{base} eq "C") {
			if (exists $genome{$pos}{C}){
				$genome{$pos}{C2} = $hash{$x}{$y}{mean};
				$genome{$pos}{C2stdev} = $hash{$x}{$y}{stdev};
				$genome{$pos}{C2cv} = $hash{$x}{$y}{cv};								
				$genome{$pos}{C2bg} = $hash{$x}{$y}{bg};
				$genome{$pos}{C2smb} = $hash{$x}{$y}{smb};
			} else {
				$genome{$pos}{C} = $hash{$x}{$y}{mean};
				$genome{$pos}{Cstdev} = $hash{$x}{$y}{stdev};
				$genome{$pos}{Ccv} = $hash{$x}{$y}{cv};				
				$genome{$pos}{Cbg} = $hash{$x}{$y}{bg};
				$genome{$pos}{Csmb} = $hash{$x}{$y}{smb};
			}
		}
		if ($annotation{$x}{$y}{base} eq "G") {
			if (exists $genome{$pos}{G}){
				$genome{$pos}{G2} = $hash{$x}{$y}{mean};
				$genome{$pos}{G2stdev} = $hash{$x}{$y}{stdev};
				$genome{$pos}{G2cv} = $hash{$x}{$y}{cv};				
				$genome{$pos}{G2bg} = $hash{$x}{$y}{bg};
				$genome{$pos}{G2smb} = $hash{$x}{$y}{smb};
			} else {
				$genome{$pos}{G} = $hash{$x}{$y}{mean};
				$genome{$pos}{Gstdev} = $hash{$x}{$y}{stdev};
				$genome{$pos}{Gcv} = $hash{$x}{$y}{cv};				
				$genome{$pos}{Gbg} = $hash{$x}{$y}{bg};
				$genome{$pos}{Gsmb} = $hash{$x}{$y}{smb};
			}
		}
	}
	if ($annotation{$x}{$y}{sense} eq "AS") { 
		if ($annotation{$x}{$y}{base} eq "A") {
			if (exists $genome{$pos}{AS_A}) {
				if (exists $genome{$pos}{AS_A2}) {print STDERR "OVERWRITING Anti-Sense $x,$y $pos ($start $stop - $genome{$pos}{A2}) ...\n";}
				$genome{$pos}{AS_A2} = $hash{$x}{$y}{mean};
				$genome{$pos}{AS_A2stdev} = $hash{$x}{$y}{stdev};
				$genome{$pos}{AS_A2cv} = $hash{$x}{$y}{cv};				
				$genome{$pos}{AS_A2bg} = $hash{$x}{$y}{bg};
				$genome{$pos}{AS_A2smb} = $hash{$x}{$y}{smb};
			} else {
				$genome{$pos}{AS_A} = $hash{$x}{$y}{mean};
				$genome{$pos}{AS_Astdev} = $hash{$x}{$y}{stdev};
				$genome{$pos}{AS_Acv} = $hash{$x}{$y}{cv};				
				$genome{$pos}{AS_Abg} = $hash{$x}{$y}{bg};
				$genome{$pos}{AS_Asmb} = $hash{$x}{$y}{smb};
			}
		}
		if ($annotation{$x}{$y}{base} eq "T") {
			if (exists $genome{$pos}{AS_T}) {
				$genome{$pos}{AS_T2} = $hash{$x}{$y}{mean};
				$genome{$pos}{AS_T2stdev} = $hash{$x}{$y}{stdev};
				$genome{$pos}{AS_T2cv} = $hash{$x}{$y}{cv};				
				$genome{$pos}{AS_T2bg} = $hash{$x}{$y}{bg};
				$genome{$pos}{AS_T2smb} = $hash{$x}{$y}{smb};
			} else {
				$genome{$pos}{AS_T} = $hash{$x}{$y}{mean};
				$genome{$pos}{AS_Tstdev} = $hash{$x}{$y}{stdev};
				$genome{$pos}{AS_Tcv} = $hash{$x}{$y}{cv};				
				$genome{$pos}{AS_Tbg} = $hash{$x}{$y}{bg};
				$genome{$pos}{AS_Tsmb} = $hash{$x}{$y}{smb};
			}
		}
		if ($annotation{$x}{$y}{base} eq "C") {
			if (exists $genome{$pos}{AS_C}) {
				$genome{$pos}{AS_C2} = $hash{$x}{$y}{mean};
				$genome{$pos}{AS_C2stdev} = $hash{$x}{$y}{stdev};
				$genome{$pos}{AS_C2cv} = $hash{$x}{$y}{cv};				
				$genome{$pos}{AS_C2bg} = $hash{$x}{$y}{bg};
				$genome{$pos}{AS_C2smb} = $hash{$x}{$y}{smb};
			} else {
				$genome{$pos}{AS_C} = $hash{$x}{$y}{mean};
				$genome{$pos}{AS_Cstdev} = $hash{$x}{$y}{stdev};
				$genome{$pos}{AS_Ccv} = $hash{$x}{$y}{cv};				
				$genome{$pos}{AS_Cbg} = $hash{$x}{$y}{bg};
				$genome{$pos}{AS_Csmb} = $hash{$x}{$y}{smb};
			}
		}
		if ($annotation{$x}{$y}{base} eq "G") {
			if (exists $genome{$pos}{AS_G}) {
				$genome{$pos}{AS_G2} = $hash{$x}{$y}{mean};
				$genome{$pos}{AS_G2stdev} = $hash{$x}{$y}{stdev};
				$genome{$pos}{AS_G2cv} = $hash{$x}{$y}{cv};				
				$genome{$pos}{AS_G2bg} = $hash{$x}{$y}{bg};
				$genome{$pos}{AS_G2smb} = $hash{$x}{$y}{smb};
			} else {
				$genome{$pos}{AS_G} = $hash{$x}{$y}{mean};
				$genome{$pos}{AS_Gstdev} = $hash{$x}{$y}{stdev};
				$genome{$pos}{AS_Gcv} = $hash{$x}{$y}{cv};				
				$genome{$pos}{AS_Gbg} = $hash{$x}{$y}{bg};
				$genome{$pos}{AS_Gsmb} = $hash{$x}{$y}{smb};
			}
		}
	}
	#print "$x,$y $pos $annotation{$x}{$y}{base} $annotation{$x}{$y}{sense} $genome{$pos}{reference} $annotation{$x}{$y}{isrefbase} *** $genome{$pos}{A},$genome{$pos}{Asmb} $genome{$pos}{T},$genome{$pos}{Tsmb} $genome{$pos}{C},$genome{$pos}{Csmb} $genome{$pos}{G},$genome{$pos}{Gsmb} *** $genome{$pos}{A2},$genome{$pos}{A2smb} $genome{$pos}{T2},$genome{$pos}{T2smb} $genome{$pos}{C2},$genome{$pos}{C2smb} $genome{$pos}{G2},$genome{$pos}{G2smb} *** $genome{$pos}{AS_A},$genome{$pos}{AS_Asmb} $genome{$pos}{AS_T},$genome{$pos}{AS_Tsmb} $genome{$pos}{AS_C},$genome{$pos}{AS_Csmb} $genome{$pos}{AS_G},$genome{$pos}{AS_Gsmb}\n";
	}
}
close $fh;

 # Print json file
open my $fo, '>', 'core0_outstats.tsv' or die "Cannot open: $!";

 	my $ct=0;
 	
 	my @allpositions = sort keys %positions;
 	my @sortedpositions = sort {$a <=> $b} @allpositions;
 	
 	print "[\n";
	foreach my $printpos (@sortedpositions) {
		$ct++;
		if ($ct == 1){
			print $fo "Position\tRefBase\tSenseCall\tASenseCall\t";
			print $fo "SenseDifference\tSenseDifferetial\tASenseDifference\tASenseDifferetial\tSenseHighestSignal\tASenseHighestSignal\tSenseLowestSignal\tASenseLowestSignal\n";
		}
		if ($ct > 1){
			my @bases = ($genome{$printpos}{A},$genome{$printpos}{T},$genome{$printpos}{C},$genome{$printpos}{G});
			my @bases2 = ($genome{$printpos}{A2},$genome{$printpos}{T2},$genome{$printpos}{C2},$genome{$printpos}{G2});
			my @asbases = ($genome{$printpos}{AS_A},$genome{$printpos}{AS_T},$genome{$printpos}{AS_C},$genome{$printpos}{AS_G});
			my @asbases2 = ($genome{$printpos}{AS_A2},$genome{$printpos}{AS_T2},$genome{$printpos}{AS_C2},$genome{$printpos}{AS_G2});

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
			if (($genome{$printpos}{A} > $genome{$printpos}{T}) && ($genome{$printpos}{A} > $genome{$printpos}{C}) && ($genome{$printpos}{A} > $genome{$printpos}{G})){
						$genome{$printpos}{SenseCall} = "A";	
			}
			if (($genome{$printpos}{T} > $genome{$printpos}{A}) && ($genome{$printpos}{T} > $genome{$printpos}{C}) && ($genome{$printpos}{T} > $genome{$printpos}{G})){
						$genome{$printpos}{SenseCall} = "T";	
			}
			if (($genome{$printpos}{C} > $genome{$printpos}{A}) && ($genome{$printpos}{C} > $genome{$printpos}{T}) && ($genome{$printpos}{C} > $genome{$printpos}{G})){
						$genome{$printpos}{SenseCall} = "C";	
			}
			if (($genome{$printpos}{G} > $genome{$printpos}{A}) && ($genome{$printpos}{G} > $genome{$printpos}{T}) && ($genome{$printpos}{G} > $genome{$printpos}{C})){
						$genome{$printpos}{SenseCall} = "G";	
			}
			if (($genome{$printpos}{AS_A} > $genome{$printpos}{AS_T}) && ($genome{$printpos}{AS_A} > $genome{$printpos}{AS_C}) && ($genome{$printpos}{AS_A} > $genome{$printpos}{AS_G})){
						$genome{$printpos}{ASenseCall} = "A";	
			}
			if (($genome{$printpos}{AS_T} > $genome{$printpos}{AS_A}) && ($genome{$printpos}{AS_T} > $genome{$printpos}{AS_C}) && ($genome{$printpos}{AS_T} > $genome{$printpos}{AS_G})){
						$genome{$printpos}{ASenseCall} = "T";	
			}
			if (($genome{$printpos}{AS_C} > $genome{$printpos}{AS_A}) && ($genome{$printpos}{AS_C} > $genome{$printpos}{AS_T}) && ($genome{$printpos}{AS_C} > $genome{$printpos}{AS_G})){
						$genome{$printpos}{ASenseCall} = "C";	
			}
			if (($genome{$printpos}{AS_G} > $genome{$printpos}{AS_A}) && ($genome{$printpos}{AS_G} > $genome{$printpos}{AS_T}) && ($genome{$printpos}{AS_G} > $genome{$printpos}{AS_C})){
						$genome{$printpos}{ASenseCall} = "G";	
			}
			
			if ($ct == 2) {
				print "\t{\n";
			} else {
				print "\t,{\n";			
			}			

			print "\t\t\"pos\": $printpos,\n";
			print "\t\t\"Ref_base\": \"$genome{$printpos}{reference}\",\n";
			print "\t\t\"SenseCall_base\": \"$genome{$printpos}{SenseCall}\",\n";
			print "\t\t\"ASenseCall_base\": \"$genome{$printpos}{ASenseCall}\",\n";
			if ($genome{$printpos}{ASenseCall} eq $genome{$printpos}{SenseCall}){
				print "\t\t\"Call_same\": \"TRUE\",\n";			
			}else{
				print "\t\t\"Call_same\": \"FALSE\",\n";			
			}
			if ($genome{$printpos}{SenseCall} eq $genome{$printpos}{reference}){
				print "\t\t\"SCall_ref\": \"TRUE\",\n";			
			}else{
				print "\t\t\"SCall_ref\": \"FALSE\",\n";			
			}
			if ($genome{$printpos}{ASenseCall} eq $genome{$printpos}{reference}){
				print "\t\t\"ASCall_ref\": \"TRUE\",\n";			
			}else{
				print "\t\t\"ASCall_ref\": \"FALSE\",\n";			
			}
			
			print $fo "$printpos\t$genome{$printpos}{reference}\t$genome{$printpos}{SenseCall}\t$genome{$printpos}{ASenseCall}\t";
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
				print "\t\t\t\t\"intensity\": $genome{$printpos}{AS_A},\n"; 
				print "\t\t\t\t\"smb\": $genome{$printpos}{AS_Asmb},\n";
				print "\t\t\t\t\"stdev\": $genome{$printpos}{AS_Astdev}\n";
#				print "\t\t\t\t\"cv\": $genome{$printpos}{AS_Acv}\n";
				print "\t\t\t},\n";
				
				print "\t\t\t{\n";
				print "\t\t\t\t\"base\": \"T\",\n";
				print "\t\t\t\t\"intensity\": $genome{$printpos}{AS_T},\n"; 
				print "\t\t\t\t\"smb\": $genome{$printpos}{AS_Tsmb},\n";
				print "\t\t\t\t\"stdev\": $genome{$printpos}{AS_Tstdev}\n";
#				print "\t\t\t\t\"cv\": $genome{$printpos}{AS_Tcv}\n";
				print "\t\t\t},\n";
				
				print "\t\t\t{\n";
				print "\t\t\t\t\"base\": \"C\",\n";
				print "\t\t\t\t\"intensity\": $genome{$printpos}{AS_C},\n"; 
				print "\t\t\t\t\"smb\": $genome{$printpos}{AS_Csmb},\n";
				print "\t\t\t\t\"stdev\": $genome{$printpos}{AS_Cstdev}\n";
#				print "\t\t\t\t\"cv\": $genome{$printpos}{AS_Ccv}\n";
				print "\t\t\t},\n";
				
				print "\t\t\t{\n";
				print "\t\t\t\t\"base\": \"G\",\n";
				print "\t\t\t\t\"intensity\": $genome{$printpos}{AS_G},\n"; 
				print "\t\t\t\t\"smb\": $genome{$printpos}{AS_Gsmb},\n";
				print "\t\t\t\t\"stdev\": $genome{$printpos}{AS_Gstdev}\n";
#				print "\t\t\t\t\"cv\": $genome{$printpos}{AS_Gcv}\n";
				print "\t\t\t}";
				
				if ($asdifference2 > 0) {
					print ",\n";

					print "\t\t\t{\n";
					print "\t\t\t\t\"base\": \"A\",\n";
					print "\t\t\t\t\"intensity\": $genome{$printpos}{AS_A2},\n"; 
					print "\t\t\t\t\"smb\": $genome{$printpos}{AS_A2smb},\n";
					print "\t\t\t\t\"stdev\": $genome{$printpos}{AS_A2stdev}\n";
#					print "\t\t\t\t\"cv\": $genome{$printpos}{AS_A2cv}\n";
					print "\t\t\t},\n";
				
					print "\t\t\t{\n";
					print "\t\t\t\t\"base\": \"T\",\n";
					print "\t\t\t\t\"intensity\": $genome{$printpos}{AS_T2},\n"; 
					print "\t\t\t\t\"smb\": $genome{$printpos}{AS_T2smb},\n";
					print "\t\t\t\t\"stdev\": $genome{$printpos}{AS_T2stdev}\n";
#					print "\t\t\t\t\"cv\": $genome{$printpos}{AS_T2cv}\n";
					print "\t\t\t},\n";
				
					print "\t\t\t{\n";
					print "\t\t\t\t\"base\": \"C\",\n";
					print "\t\t\t\t\"intensity\": $genome{$printpos}{AS_C2},\n"; 
					print "\t\t\t\t\"smb\": $genome{$printpos}{AS_C2smb},\n";
					print "\t\t\t\t\"stdev\": $genome{$printpos}{AS_C2stdev}\n";
#					print "\t\t\t\t\"cv\": $genome{$printpos}{AS_C2cv}\n";
					print "\t\t\t},\n";
				
					print "\t\t\t{\n";
					print "\t\t\t\t\"base\": \"G\",\n";
					print "\t\t\t\t\"intensity\": $genome{$printpos}{AS_G2},\n"; 
					print "\t\t\t\t\"smb\": $genome{$printpos}{AS_G2smb},\n";
					print "\t\t\t\t\"stdev\": $genome{$printpos}{AS_G2stdev}\n";
#					print "\t\t\t\t\"cv\": $genome{$printpos}{AS_G2cv}\n";
					print "\t\t\t}";
				}
				print "\n\t\t],\n";

			}
	
			if ($sdifference > -.1) {
				print "\t\t\"sense\": [\n";
				print "\t\t\t{\n";
				print "\t\t\t\t\"base\": \"A\",\n";
				print "\t\t\t\t\"intensity\": $genome{$printpos}{A},\n"; 
				print "\t\t\t\t\"smb\": $genome{$printpos}{Asmb},\n";
				print "\t\t\t\t\"stdev\": $genome{$printpos}{Astdev}\n";
#				print "\t\t\t\t\"cv\": $genome{$printpos}{Acv}\n";
				print "\t\t\t},\n";
				
				print "\t\t\t{\n";
				print "\t\t\t\t\"base\": \"T\",\n";
				print "\t\t\t\t\"intensity\": $genome{$printpos}{T},\n"; 
				print "\t\t\t\t\"smb\": $genome{$printpos}{Tsmb},\n";
				print "\t\t\t\t\"stdev\": $genome{$printpos}{Tstdev}\n";
#				print "\t\t\t\t\"cv\": $genome{$printpos}{Tcv}\n";
				print "\t\t\t},\n";
				
				print "\t\t\t{\n";
				print "\t\t\t\t\"base\": \"C\",\n";
				print "\t\t\t\t\"intensity\": $genome{$printpos}{C},\n"; 
				print "\t\t\t\t\"smb\": $genome{$printpos}{Csmb},\n";
				print "\t\t\t\t\"stdev\": $genome{$printpos}{Cstdev}\n";
#				print "\t\t\t\t\"cv\": $genome{$printpos}{Ccv}\n";
				print "\t\t\t},\n";
				
				print "\t\t\t{\n";
				print "\t\t\t\t\"base\": \"G\",\n";
				print "\t\t\t\t\"intensity\": $genome{$printpos}{G},\n"; 
				print "\t\t\t\t\"smb\": $genome{$printpos}{Gsmb},\n";
				print "\t\t\t\t\"stdev\": $genome{$printpos}{Gstdev}\n";
#				print "\t\t\t\t\"cv\": $genome{$printpos}{Gcv}\n";				
				print "\t\t\t}";
				
				if ($sdifference2 > 0) {
					print ",\n";

					print "\t\t\t{\n";
					print "\t\t\t\t\"base\": \"A\",\n";
					print "\t\t\t\t\"intensity\": $genome{$printpos}{A2},\n"; 
					print "\t\t\t\t\"smb\": $genome{$printpos}{A2smb},\n";
					print "\t\t\t\t\"stdev\": $genome{$printpos}{A2stdev}\n";
#					print "\t\t\t\t\"cv\": $genome{$printpos}{A2cv}\n";					
					print "\t\t\t},\n";
				
					print "\t\t\t{\n";
					print "\t\t\t\t\"base\": \"T\",\n";
					print "\t\t\t\t\"intensity\": $genome{$printpos}{T2},\n"; 
					print "\t\t\t\t\"smb\": $genome{$printpos}{T2smb},\n";
					print "\t\t\t\t\"stdev\": $genome{$printpos}{T2stdev}\n";
#					print "\t\t\t\t\"cv\": $genome{$printpos}{T2cv}\n";					
					print "\t\t\t},\n";
				
					print "\t\t\t{\n";
					print "\t\t\t\t\"base\": \"C\",\n";
					print "\t\t\t\t\"intensity\": $genome{$printpos}{C2},\n"; 
					print "\t\t\t\t\"smb\": $genome{$printpos}{C2smb},\n";
					print "\t\t\t\t\"stdev\": $genome{$printpos}{C2stdev}\n";
#					print "\t\t\t\t\"cv\": $genome{$printpos}{C2cv}\n";					
					print "\t\t\t},\n";
				
					print "\t\t\t{\n";
					print "\t\t\t\t\"base\": \"G\",\n";
					print "\t\t\t\t\"intensity\": $genome{$printpos}{G2},\n"; 
					print "\t\t\t\t\"smb\": $genome{$printpos}{G2smb},\n";
					print "\t\t\t\t\"stdev\": $genome{$printpos}{G2stdev}\n";
#					print "\t\t\t\t\"cv\": $genome{$printpos}{G2cv}\n";										
					print "\t\t\t}";
				}
				print "\n\t\t]\n";
			}
		print "\t}\n";						
		}
	}
	print "]";
close $fo;





