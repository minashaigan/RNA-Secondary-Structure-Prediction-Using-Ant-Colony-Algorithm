#!d:/Applications/Strawberry/perl/bin/perl.exe
#Rna Secondary structure with ant colony
use warnings;
use strict;

open (my $in,"D:/.RNA Secondary Structure Prediction Using Ant Colony Algorithm/Data/RNA-data/adenin.txt") or die "Can't read $!";  
my $RNA = <$in>;
$RNA = uc($RNA);
chomp($RNA);
$RNA=~ s/\R//g;
my $Structure = <$in>;
chomp($Structure);
$Structure=~ s/\R//g;
my $main_energy = CompleteEnergy($Structure, $RNA);
open(my $out,'>',"result.txt") or die "Can't open file for writing: $!"; 
close $out or die "Failed to close file: $!";

my @seq = split('',$RNA);
my @stems = Stems();
my %stack_energy = stack_energy();
my @hi = Heuristic_Information();
my @pheremones = Initialization_Pheremone();
my ($IL, $ants_number) = IL();

print_Stems_inFile();
# print "pheremone\n";
print_pheremones_inFile();
# <STDIN>;
main();

## Functions

# main
sub main {
	my @temp_best_ant = ();
	my $best_energy;
	my $best_structure;
	my $ii = 1;
	
	my @first_position_ants = ();
	my @first_position_ants_allowed = ();
	for(my $i=0; $i<$ants_number; $i++ ){
		my @allowed = (0..scalar(@stems)-1);
		my $first_stem = scalar(@stems)-1-$i;
		my @solution = ();
		push @solution, $first_stem;
		@allowed = Update_allowed(\@allowed,\@solution);
		push @first_position_ants, [@solution];
		push @first_position_ants_allowed, [@allowed];
	}
	# print 
	# for(my $j=0; $j<$ants_number; $j++){
		# print join(",",@{$first_position_ants[$j]}),"\n";
		
		# print "****\n";
	# }
	# <STDIN>;
	while( $ii <= 10 )
	{
		open(my $out,'>>',"result.txt") or die "Can't open file for writing: $!"; 
		print $ii,"\n";
		my @ants = @first_position_ants;
		my @ants_allowed = @first_position_ants_allowed;
		
		my $all_allowed_empty = 0;
		while( !$all_allowed_empty ){
			for(my $i=0; $i<$ants_number; $i++ ){
				my ($first, $second) = Select_Next_Stem(\@{$ants_allowed[$i]}, \@{$ants[$i]});
				@{$ants_allowed[$i]} = @$first;
				@{$ants[$i]} = @$second;
			}
			# update just last edge of all ants
			for(my $i=0; $i<$ants_number; $i++ ){
				Local_Update_Pheremone(\@{$ants_allowed[$i]}, \@{$ants[$i]});
				
			}
			my $flag = 1;
			for(my $i=0; $i<$ants_number; $i++ ){
				if( scalar(@{$ants_allowed[$i]}) ){
					$flag = 0;
					last;
				}
			}
			if( $flag == 1 ){
				$all_allowed_empty = 1;
			}
			# for(my $j=0; $j<$ants_number; $j++){
				# print join(", ",@{$ants[$j]}),"\n";
				
				# print "****\n";
			# }
			# <STDIN>;
		}

		# update edges of best ant
		my ($first, $second) = Global_Best_ant(\@ants);
		$best_energy = $first;
		my @best_ant = @$second;
		$best_structure = Structure(\@best_ant);
		Global_Update_Pheremone(\@best_ant, $best_energy);
		
		my $check_sub_solution = Change_Sub_Solution(\@best_ant, \@temp_best_ant);
		@temp_best_ant = @best_ant;

		if( !$check_sub_solution ){
			$ii++;
		}
		else{
			$ii=0;
			print $out $best_structure;
			print $out " ",$best_energy,"\n";
		}
		
		close $out or die "Failed to close file: $!";
		# print "one iteration finished\n";
		print_pheremones_inFile();
		# <STDIN>;
	}
	print_result($best_structure);
}

# Create Dot Plot
sub Dot_Plot {	
	my @matrix = ();
	for(my $i=0;$i<scalar(@seq);$i++){
		for(my $j=0;$j<scalar(@seq);$j++){
			if(($seq[$i] eq "G" and $seq[$j] eq "C") or ($seq[$i] eq "A" and $seq[$j] eq "U") or ($seq[$i] eq "C" and $seq[$j] eq "G") or ($seq[$i] eq "U" and $seq[$j] eq "A")  or ($seq[$i] eq "G" and $seq[$j] eq "U") or ($seq[$i] eq "U" and $seq[$j] eq "G")){
				$matrix[$i][$j] = 1;
			}
			else{
				$matrix[$i][$j] = 0;
			}
		}
	}
	return @matrix;
}

# Declare Stems
sub Stems {
	my @matrix = Dot_Plot();
	my @stems = ();
	#first declare stems with length >=3
	#u1 : initial ribonucleotide position
	#u2 : final ribonucleotide position
	#u3 : length of the stem
	for(my $i=0;$i<scalar(@seq);++$i){
		for(my $j=scalar(@seq)-1;$j>=0;--$j){
			if($j<=$i){
				last;
			}
			if($matrix[$i][$j] == 1){
				my $start=$i;
				my $end=$j;
				my $k;
				for($k=1;$k<scalar(@seq);++$k){
					if(($end-$start-(2*$k)+1)>=3 and $k>=3){
						push @stems, {u1=>$start, u2=>$end, u3=>$k };
					}
					if($matrix[$i+$k][$j-$k] == 0){
						last;
					}
					if($j-$k<=$i+$k){	
						last;
					}
				}

			}
		}
	}
	my @sorted =  sort { $a->{u3} <=> $b->{u3} } @stems;
	@stems = @sorted;
	return @stems;
}

#print stems in file
sub print_Stems_inFile {
	open(my $out,'>',"Stems.txt") or die "Can't open file for writing: $!"; 
	for(my $i=0;$i<scalar(@stems);$i++)
	{
		print $out $i," => ";
		print $out $stems[$i]{u1}," ",$stems[$i]{u2}," ",$stems[$i]{u3}," ";
		for(my $j=0;$j<$stems[$i]{u3};$j++){
			print $out $seq[$stems[$i]{u1}+$j];
		}
		print $out " ";
		for(my $j=$stems[$i]{u3}-1;$j>=0;$j--){
			print $out $seq[$stems[$i]{u2}-$j];
		}
		print $out "\n";
	}
	close $out or die "Failed to close file: $!";
}

#print pheremones in file
sub print_pheremones_inFile {
	open(my $out,'>',"pheremones.txt") or die "Can't open file for writing: $!"; 
	for(my $i=0;$i<scalar(@stems);$i++)
	{
		for(my $j=0;$j<scalar(@stems);$j++)
		{
			print $out $pheremones[$i][$j],"\n";
		}
	}
	close $out or die "Failed to close file: $!";
}

# Check two stems is consistent or not
sub Consistence {
	my $stem_index1 = $_[0];
	my $stem_index2 = $_[1];
	if( $stems[$stem_index1]{u1}<$stems[$stem_index2]{u1} ){
		if( ($stems[$stem_index1]{u1}+$stems[$stem_index1]{u3} <= $stems[$stem_index2]{u1}) and ($stems[$stem_index2]{u2} <= $stems[$stem_index1]{u2}-$stems[$stem_index1]{u3}) ){
			return 1;
		}
		if( $stems[$stem_index1]{u2} < $stems[$stem_index2]{u1} ){
			return 2;
		}
	}
	else{
		if( ($stems[$stem_index2]{u1}+$stems[$stem_index2]{u3} <= $stems[$stem_index1]{u1}) and ($stems[$stem_index1]{u2} <= $stems[$stem_index2]{u2}-$stems[$stem_index2]{u3}) ){
			return 3;
		}
		if( $stems[$stem_index2]{u2} < $stems[$stem_index1]{u1} ){
			return 4;
		}
	}
	return 0;
}

# Compute Heuristic Information
sub Heuristic_Information {
	my @hi = ();
	for(my $i=0;$i<scalar(@stems);$i++){
		for(my $j=0;$j<scalar(@stems);$j++){
			my $consistent = Consistence($i,$j);
			my $energy_current_stem = calculate_Stem_Energy($i);
			my $energy_next_stem = calculate_Stem_Energy($j);
			
			# if( $consistent){
				# $hi[$i][$j] =  ($energy_current_stem*$energy_next_stem) ;
			# }
			if( $consistent == 1){
				my $d1 = $stems[$j]{u1} - ( $stems[$i]{u1} + $stems[$i]{u3} );
				my $d2 = ( $stems[$i]{u2} - $stems[$i]{u3} ) - $stems[$j]{u2};
				if($d1 == 0 and $d2 == 0){
					$hi[$i][$j] =  ($energy_current_stem*$energy_next_stem) ;
				}
				else{
					$hi[$i][$j] = ($energy_current_stem*$energy_next_stem)/( sqrt( ($d1*$d1) + ($d2*$d2) ) );
				}
			}
			if( $consistent == 3){
				my $d1 = $stems[$i]{u1} - ( $stems[$j]{u1} + $stems[$j]{u3} );
				my $d2 = ( $stems[$j]{u2} - $stems[$j]{u3} ) - $stems[$i]{u2};
				if($d1 == 0 and $d2 == 0){
					$hi[$i][$j] = ($energy_current_stem*$energy_next_stem);
				}
				else{
					$hi[$i][$j] = ($energy_current_stem*$energy_next_stem)/( sqrt( ($d1*$d1) + ($d2*$d2) ) );
				}
			}
			if( $consistent == 2 ){
				my $d3 = $stems[$j]{u1} - $stems[$i]{u2} - 1;
				if($d3 == 0){
					$hi[$i][$j] = ($energy_current_stem*$energy_next_stem);
				}
				else{
					$hi[$i][$j] = ($energy_current_stem*$energy_next_stem)/$d3;
				}
			}
			if( $consistent == 4 ){
				my $d3 = $stems[$i]{u1} - $stems[$j]{u2} - 1;
				if($d3 == 0){
					$hi[$i][$j] = ($energy_current_stem*$energy_next_stem);
				}
				else{
					$hi[$i][$j] = ($energy_current_stem*$energy_next_stem)/$d3;
				}
			}
			if( $consistent == 0 ){
				$hi[$i][$j] = 0;
			}
			
		}	
		
	}
	return @hi;
}

# Compute Initialization Pheremone 
sub Initialization_Pheremone {
	my @pheremones = ();
	# for(my $i=0;$i<scalar(@stems);$i++){
		# for(my $j=0;$j<scalar(@stems);$j++){
			# $pheremones[$i][$j] = 1/(scalar(@stems)*scalar(@seq));
			##$pheremones[$i][$j] = 1/(scalar(@seq));
		# }
	# }
	for(my $i=0;$i<scalar(@stems);$i++){
		my $number_consistent = 0;
		for(my $j=0;$j<scalar(@stems);$j++){
			my $consistent = Consistence($i,$j);
			if( $consistent ){
				$number_consistent++;
			}
		}
		for(my $j=0;$j<scalar(@stems);$j++){
			my $consistent = Consistence($i,$j);
			if( $consistent ){
				$pheremones[$i][$j] = 1/$number_consistent;
			}
			else{
				$pheremones[$i][$j] = 0;
			}
		}
	}
	return @pheremones;
}

# Comput IL
sub IL {
	my $landa = 0.13;
	my $number_stems = scalar(@stems);
	my $number_limit_stems = 0;
	my $IL = $stems[-1]{u3};
	for(my $i=scalar(@stems)-1;$i>=0;$i--){
		if( $stems[$i]{u3}>= $IL ){
			$number_limit_stems++;
		}
		else{
			if( $number_limit_stems >= $landa*$number_stems ){
				return ($IL,$number_limit_stems);
			}
			else{
				$IL = $stems[$i]{u3};
				$number_limit_stems++;
			}
		}
	}
	return "IL NOT FOUND\n";
	<STDIN>;
}

# The choice of the first stem
sub First_Stem {
	my @limit_stems = ();
	for(my $i=scalar(@stems)-1;$i>=0;$i--){
		if($stems[$i]{u3}>= $IL){
			push @limit_stems,$i;
		}
		else{
			last;
		}
	}
	my $rnd_stem = int(rand(scalar(@limit_stems)));
	return $limit_stems[$rnd_stem];
}

# Update consistent stems with solution
sub Update_allowed {
	my @allowed = @{$_[0]};
	my @solution = @{$_[1]};
	for(my $i=0;$i<scalar(@allowed);$i++){
		for(my $j=0;$j<scalar(@solution);$j++){	
			my $consistent = Consistence($allowed[$i],$solution[$j]);
			if( !$consistent ){
				splice(@allowed,$i,1);
				$i--;
				last;
			}
		}
	}
	return @allowed;
}

# The choice of the next stem
sub Select_Next_Stem {
	my @allowed = @{$_[0]};
	my @solution = @{$_[1]};
	my $current_stem = $solution[-1];
	my $q = 0.9;
	my $q_rnd = rand(1);
	if( $q_rnd <= $q and scalar(@allowed)){
		my $max = $pheremones[$current_stem][$allowed[0]]*$hi[$current_stem][$allowed[0]]*$hi[$current_stem][$allowed[0]];
		my $max_allowed = $allowed[0];
		for(my $i=1;$i<scalar(@allowed);$i++){
			my $max_pheromone = $pheremones[$current_stem][$allowed[$i]]*$hi[$current_stem][$allowed[$i]]*$hi[$current_stem][$allowed[$i]];
			if( $max < $max_pheromone ){
				$max_allowed = $allowed[$i];
				$max = $max_pheromone;
			}
		}
		push @solution,$max_allowed;
		@allowed = Update_allowed(\@allowed, \@solution);
	}
	else{
		my $sum = 0;
		for(my $i=0;$i<scalar(@allowed);$i++){
			$sum += $pheremones[$current_stem][$allowed[$i]]*$hi[$current_stem][$allowed[$i]]*$hi[$current_stem][$allowed[$i]];
		}
		my $rnd = rand(1);
		
		my $probabilty_current_stem = 0;
		for(my $i=0;$i<scalar(@allowed);$i++){
			$probabilty_current_stem += ($pheremones[$current_stem][$allowed[$i]]*$hi[$current_stem][$allowed[$i]]*$hi[$current_stem][$allowed[$i]])/$sum;
			if( $rnd<$probabilty_current_stem ){
				push @solution,$allowed[$i];
				@allowed = Update_allowed(\@allowed, \@solution);
				last;
			}
		}
	}
	# print "finished\n";
	return (\@allowed,\@solution);
}

# ACS Local Updating Rule LAST EDGE
sub Local_Update_Pheremone {
	my @allowed = @{$_[0]};
	my @solution = @{$_[1]};
	my $current_stem = $solution[-1];
	if($solution[-2]){
		my $previous_stem = $solution[-2];
		my $rou = 0.1;
		$pheremones[$current_stem][$previous_stem] = ((1-$rou)*$pheremones[$current_stem][$previous_stem]) + ($rou/(scalar(@stems)*scalar(@seq)));
		$pheremones[$previous_stem][$current_stem] = ((1-$rou)*$pheremones[$previous_stem][$current_stem]) + ($rou/(scalar(@stems)*scalar(@seq)));
	}
}

# Update Sub Solution Set
sub Global_Best_ant {
	my @ants = @{$_[0]};
	my @best_solution = ();
	my $min_energy = 1000;
	for(my $i=0;$i<scalar(@ants);$i++){
		my @solution = @{$ants[$i]};
		my $structure = Structure(\@solution);
		my $energy = CompleteEnergy($structure, $RNA);
		if( $energy < $min_energy ){
			@best_solution = @solution;
			$min_energy = $energy;
		}
	}
	return ($min_energy,\@best_solution);
}


#find structure of a solution
sub Structure {
	my @solution = @{$_[0]};
	my @sp ;
	for(my $i=0;$i<scalar(@seq);$i++){
		$sp[$i] = ".";
	}
	for(my $i=0;$i<scalar(@solution);$i++){
		for(my $j=0;$j<$stems[$solution[$i]]{u3};$j++){
			$sp[$stems[$solution[$i]]{u1}+$j]="(";
			$sp[$stems[$solution[$i]]{u2}-$j]=")";
		}
	}
	my $Sp = join("",@sp);
	return $Sp;
}

# Calculate Energy Solution
sub CompleteEnergy{
	my($SP,$RNA)=@_;
	my($ES);	
	open(filehandelO3,'>',"antsolution.txt");
	print(filehandelO3	$RNA,"\n",$SP,"\n","@\n");
	close(filehandelO3);
	$ES= qx/RNAeval.exe <antsolution.txt/;
	$ES=~ s/.*\n//;
	$ES=~ s/.* //;
	$ES=~ s/\(//;
	$ES=~ s/\)//;
	return($ES);
}

# Global_Update_Pheremone
sub Global_Update_Pheremone {
	my @solution = @{$_[0]};
	my $energy = $_[1];
	my $alpha = 0.1;
	for(my $m=0;$m<scalar(@solution);$m++){
		for(my $n=$m+1;$n<scalar(@solution);$n++){
			$pheremones[$solution[$n]][$solution[$m]] = ((1-$alpha)*$pheremones[$solution[$n]][$solution[$m]]) - ($alpha/$energy);
			$pheremones[$solution[$m]][$solution[$n]] = ((1-$alpha)*$pheremones[$solution[$m]][$solution[$n]]) - ($alpha/$energy);
		}
	}
}

# check whether the Sub Solution Set Change or not
sub Change_Sub_Solution {
	my @new = @{$_[0]};
	my @old = @{$_[1]};

	if( scalar(@new) == scalar(@old) ){
		my $new_structure = Structure(\@new);
		my $old_structure = Structure(\@old);
		if( $new_structure eq $old_structure ){
			return 0;
		}
	}
	return 1;
}

# Stack Energy

sub stack_energy {
	my %stack_energy;
	$stack_energy{"AA"}{"UU"} = -0.9;
	$stack_energy{"AC"}{"UG"} = -2.2;
	$stack_energy{"AG"}{"UC"} = -2.1;
	$stack_energy{"AG"}{"UU"} = -0.6;
	$stack_energy{"AU"}{"UA"} = -1.1;
	$stack_energy{"AU"}{"UG"} = -1.4;
	$stack_energy{"CA"}{"GU"} = -2.1;
	$stack_energy{"CC"}{"GG"} = -3.3;
	$stack_energy{"CG"}{"GC"} = -2.4;
	$stack_energy{"CG"}{"GU"} = -1.4;
	$stack_energy{"CU"}{"GA"} = -2.1;
	$stack_energy{"CU"}{"GG"} = -2.1;
	$stack_energy{"GA"}{"CU"} = -2.4;
	$stack_energy{"GA"}{"UU"} = -1.3;
	$stack_energy{"GC"}{"CG"} = -3.4;
	$stack_energy{"GC"}{"UG"} = -2.5;
	$stack_energy{"GG"}{"CC"} = -3.3;
	$stack_energy{"GG"}{"CU"} = -1.5;
	$stack_energy{"GG"}{"UC"} = -2.1;
	$stack_energy{"GG"}{"UU"} = -0.5;
	$stack_energy{"GU"}{"CA"} = -2.2;
	$stack_energy{"GU"}{"CG"} = -2.5;
	$stack_energy{"GU"}{"UA"} = -1.4;
	$stack_energy{"GU"}{"UG"} = 1.3;
	$stack_energy{"UA"}{"AU"} = -1.3;
	$stack_energy{"UA"}{"GU"} = -1;
	$stack_energy{"UC"}{"AG"} = -2.4;
	$stack_energy{"UC"}{"GG"} = -1.5;
	$stack_energy{"UG"}{"AC"} = -2.1;
	$stack_energy{"UG"}{"AU"} = -1;
	$stack_energy{"UG"}{"GC"} = -1.4;
	$stack_energy{"UG"}{"GU"} = 0.3;
	$stack_energy{"UU"}{"AA"} = -0.9;
	$stack_energy{"UU"}{"AG"} = -1.3;
	$stack_energy{"UU"}{"GA"} = -0.6;
	$stack_energy{"UU"}{"GG"} = -0.5;
	return %stack_energy;
}

# Calculate stem energy

sub calculate_Stem_Energy {
	my $stem_index = $_[0];
	
	my $sum_energy = 0;
	for(my $j=0;$j<$stems[$stem_index]{u3}-1;++$j){
		foreach my $top (sort { $a cmp $b} keys %stack_energy){
			foreach my $bottom ( keys %{$stack_energy{$top}}){
				my $pair1 = join('',$seq[$stems[$stem_index]{u1}+$j],$seq[$stems[$stem_index]{u1}+$j+1]);
				my $pair2 = join('',$seq[$stems[$stem_index]{u2}-$j],$seq[$stems[$stem_index]{u2}-$j-1]);
				if($pair1 eq $top and $pair2 eq $bottom)
				{
					$sum_energy = $sum_energy + $stack_energy{$top}{$bottom};
				}
			}
		}	
	}
	return $sum_energy;
}


# print
sub print_result {
	my $best_Structure = $_[0];
	my $j;
	open(my $out,'>>',"result.txt") or die "Can't open file for writing: $!"; 
	print $out $best_Structure," : best structure\n";
	print $out $Structure," : main structure\n";
	print $out "Main Energy : ",$main_energy;
	print $out "Len : ",length($RNA),"\n";
	print $out "IL : ",$IL,"\n";
	my @Stack = ();
	my @MainStr = ();
	my $se2 = 0;
	my $ppv2 = 0;
	
	for(my $i=0;$i<length($RNA);++$i){
		if(substr($Structure,$i,1) eq "."){
			$MainStr[$i] = $i;
		}
		elsif(substr($Structure,$i,1) eq "("){
			push(@Stack,$i);
			$se2++;
		}
		else{
			$j=pop(@Stack);
			$MainStr[$i]=$j;
			$MainStr[$j]=$i;
		}
	}
	@Stack="";
	pop(@Stack);
	my @PredictStr = ();
	for(my $i=0;$i<length($RNA);++$i){
		if(substr($best_Structure,$i,1) eq "."){
			$PredictStr[$i]=$i;
		}
		elsif(substr($best_Structure,$i,1) eq "("){
			push(@Stack,$i);
			$ppv2++;
		}
		else{
			$j=pop(@Stack);
			$PredictStr[$i]=$j;
			$PredictStr[$j]=$i;
		}
	}
	my $TP=0;
	my $FP=0;
	my $FN=0;
	for(my $i=0;$i<length($RNA);++$i){
		if($PredictStr[$i] eq $MainStr[$i] ){
			if($PredictStr[$i] != $i){
				++$TP;
			}
		}
		else{
			++$FP;
		}
	}
	$TP = $TP/2;
	$FP = $ppv2-$TP;
	$FN = $se2-$TP;
	print $out "TP : ",$TP,"  FP : ",$FP,"  FN : ",$FN,"\n";
	
	
	print $out "sensitivity : ",$TP/$se2,"\n";
	
	print $out "PPV : ",$TP/$ppv2,"\n";
	
	print $out "f-measure : ",(2*($TP/$se2)*($TP/$ppv2))/(($TP/$se2)+($TP/$ppv2)),"\n";
		
	close $out or die "Failed to close file: $!";
}