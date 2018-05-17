#Rna Secondary structure with ant colony
use warnings;
use strict;

main();

sub main {
	open (my $in,"D:/.RNA Secondary Structure Prediction Using Ant Colony Algorithm/Data/RNA-data/adenin.txt") or die "Can't read $!";  
	my $RNA = <$in>;
	$RNA = uc($RNA);
	chomp($RNA);
	$RNA=~ s/\R//g;
	my $Structure = <$in>;
	chomp($Structure);
	$Structure=~ s/\R//g;
	open(my $out,'>',"result.txt") or die "Can't open file for writing: $!"; 
	close $out or die "Failed to close file: $!";
	
	my @RNAseq = split('',$RNA);
	my @DotPlot = Dot_Plot(@RNAseq);
	my @Stems = Stems(@RNAseq);
	my @HI = Heuristic_Information(@Stems);
	my @Pheremones = Initialization_Pheremone(@Stems);
	my $IL = IL(@Stems);
	my $i = 1;
	
	my $min_energy = 100000;
	my @sub_solution = ();
	my $sub_solution;
	while( $i <= 10 )
	{
		open(my $out,'>>',"result.txt") or die "Can't open file for writing: $!"; 
		print $i,"\n";
		my $k = 1;
		my @ants_solution = ();
		while( $k <= 20 )
		{
			my @allowed = (0..scalar(@Stems)-1);
			
			my $first_stem = First_Stem($IL, \@Stems);
			my @solution = ();
			push @solution, $first_stem;
			@allowed = Update_allowed(\@allowed, \@solution, \@Stems);
			$k++;
			
			while(scalar(@allowed)){
				my ($first, $second) = Select_Next_Stem(\@allowed, \@solution, \@Stems, \@Pheremones, \@HI);
				@allowed = @$first;
				@solution = @$second;
			}
			push @ants_solution,[@solution];
		}

		my @temp_sub_solution = Update_sub_solution(\@sub_solution, \@ants_solution, \@Stems, $RNA, $min_energy);
		my $check_sub_solution = Change_Sub_Solution(\@sub_solution, \@temp_sub_solution, \@Stems, $RNA);
		
		@sub_solution = @temp_sub_solution;
		$min_energy = Update_min_energy(\@ants_solution, \@Stems, $RNA, $min_energy);
		$sub_solution = Structure(\@sub_solution, \@Stems, length($RNA));

		if( !$check_sub_solution ){

			$i++;
		}
		else{
			$i=0;
			print $out $sub_solution;
			print $out " ",$min_energy,"\n";
		}
		
		@Pheremones = Update_Pheremone(\@Pheremones, \@ants_solution, \@Stems, $min_energy, $RNA);
		# $i++;
		close $out or die "Failed to close file: $!";
	}
	
	print_result($Structure, $sub_solution, $RNA, $IL);

}

# Create Dot Plot
sub Dot_Plot {	
	my @seq = @_;
	my @matrix = ();
	for(my $i=0;$i<scalar(@seq);$i++){
		for(my $j=0;$j<scalar(@seq);$j++){
			if(($seq[$i] eq "G" and $seq[$j] eq "C") or ($seq[$i]eq"A" and $seq[$j]eq"U") or ($seq[$i]eq"C" and $seq[$j]eq"G") or ($seq[$i]eq"U" and $seq[$j]eq"A")  or ($seq[$i] eq "G" and $seq[$j] eq "U") or ($seq[$i] eq "U" and $seq[$j] eq "G")){
				$matrix[$i][$j]=1;
			}
			else{
				$matrix[$i][$j]=0;
			}
		}
	}
	return @matrix;
}

# Declare Stems
sub Stems {
	my @seq = @_;
	my @matrix = Dot_Plot(@seq);
	my @stems = ();
	# first declare stems with length >=3
	# u1 : initial ribonucleotide position
	# u2 : final ribonucleotide position
	# u3 : length of the stem
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
					if($matrix[$i+$k][$j-$k] == 0){
						if(abs($i+$k-1-($j-$k+1)+1) < 3){
							$k--;
						}
						last;
					}
					if($j-$k<=$i+$k){
						$k-=2;
						last;
					}
				}
				my $length=$k;
				if($length>=3){
					push @stems, {u1=>$start, u2=>$end, u3=>$length };
				}
			}
		}
	}
	my @sorted =  sort { $a->{u3} <=> $b->{u3} } @stems;
	@stems = @sorted;
	return @stems;
}

# Check two stems is consistent or not
sub Consistence {
	my %stem1 = %{$_[0]};
	my %stem2 = %{$_[1]};
	if( $stem1{u1}<$stem2{u1} ){
		if( ($stem1{u1}+$stem1{u3} <= $stem2{u1}) and ($stem2{u2} <= $stem1{u2}-$stem1{u3}) ){
			return 1;
		}
		if( $stem1{u2} < $stem2{u1} ){
			return 1;
		}
	}
	else{
		if( ($stem2{u1}+$stem2{u3} <= $stem1{u1}) and ($stem1{u2} <= $stem2{u2}-$stem2{u3}) ){
			return 1;
		}
		if( $stem2{u2} < $stem1{u1} ){
			return 1;
		}
	}
	return 0;
}

# Compute Heuristic Information
sub Heuristic_Information {
	my @stems = @_;
	my @hi = ();
	for(my $i=0;$i<scalar(@stems);$i++){
		for(my $j=0;$j<scalar(@stems);$j++){
			my $consistent = Consistence(\%{$stems[$i]},\%{$stems[$j]});
			if( $consistent ){
				$hi[$i][$j] = ($stems[$j]{u3}*$stems[$j]{u3})/$stems[$i]{u3};
			}
			else{
				$hi[$i][$j] = 0;
			}
		}
	}
	return @hi;
}

# Compute Initialization Pheremone 
sub Initialization_Pheremone {
	my @stems = @_;
	my @pheremones = ();
	for(my $i=0;$i<scalar(@stems);$i++){
		my $number_consistent = 0;
		for(my $j=0;$j<scalar(@stems);$j++){
			my $consistent = Consistence(\%{$stems[$i]},\%{$stems[$j]});
			if( $consistent ){
				$number_consistent++;
			}
		}
		for(my $j=0;$j<scalar(@stems);$j++){
			my $consistent = Consistence(\%{$stems[$i]},\%{$stems[$j]});
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

# Compute IL
sub IL {
	my @stems = @_;
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
				return $IL;
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
	my $IL = $_[0];
	my @stems = @{$_[1]};
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
	my @stems = @{$_[2]};
	for(my $i=0;$i<scalar(@allowed);$i++){
		for(my $j=0;$j<scalar(@solution);$j++){	
			my $consistent = Consistence(\%{$stems[$allowed[$i]]},\%{$stems[$solution[$j]]});
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
	my @stems = @{$_[2]};
	my @pheremones = @{$_[3]};
	my @hi = @{$_[4]};
	my $current_stem = $solution[-1];
	my $sum = 0;
	for(my $i=0;$i<scalar(@allowed);$i++){
		my $consistent = Consistence(\%{$stems[$current_stem]},\%{$stems[$allowed[$i]]});
		if( $consistent ){
			$sum += $pheremones[$current_stem][$allowed[$i]]*$pheremones[$current_stem][$allowed[$i]]*$hi[$current_stem][$allowed[$i]];
		}
	}
	my $x0 = rand(1);
	my $rnd =  4*$x0*(1-$x0);
	# my $rnd = $x0;
	
	my $probabilty_current_stem = 0;
	for(my $i=0;$i<scalar(@allowed);$i++){
		$probabilty_current_stem += ($pheremones[$current_stem][$allowed[$i]]*$pheremones[$current_stem][$allowed[$i]]*$pheremones[$current_stem][$allowed[$i]]*$hi[$current_stem][$allowed[$i]])/$sum;
		if( $rnd<$probabilty_current_stem ){
			push @solution,$allowed[$i];
			@allowed = Update_allowed(\@allowed, \@solution, \@stems);
			last;
		}
	}
	return (\@allowed,\@solution);
}

# Update Sub Solution Set
sub Update_sub_solution {
	my @sub_solution = @{$_[0]};
	my @ants = @{$_[1]};
	my @stems = @{$_[2]};	
	my $RNAseq = $_[3];
	my $min_energy = $_[4];
	
	for(my $i=0;$i<scalar(@ants);$i++){
		my @solution = @{$ants[$i]};
		my $structure = Structure(\@solution, \@stems, length($RNAseq));
		my $energy = CompleteEnergy($structure, $RNAseq);
		if( $energy < (-0.1*$min_energy+$min_energy) ){
		# if( $energy < 0.1*$min_energy ){
			@sub_solution = @solution;
		}
	}
	return @sub_solution;
}

# Update min energy
sub Update_min_energy {
	my @ants = @{$_[0]};
	my @stems = @{$_[1]};	
	my $RNAseq = $_[2];
	my $min_energy = $_[3];

	for(my $i=0;$i<scalar(@ants);$i++){
		my @solution = @{$ants[$i]};
		my $structure = Structure(\@solution, \@stems, length($RNAseq));
		my $energy = CompleteEnergy($structure, $RNAseq);
		if( $energy < (-0.1*$min_energy+$min_energy) ){
		# if( $energy < 0.1*$min_energy ){
			$min_energy = $energy;
		}
	}
	return $min_energy;
}

#find structure of a solution
sub Structure {
	my @solution = @{$_[0]};
	my @stems = @{$_[1]};	
	my $seq_length = $_[2];
	my @sp ;
	for(my $i=0;$i<$seq_length;$i++){
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

# Update Pheremone
sub Update_Pheremone{
	my @pheremones = @{$_[0]};
	my @ants = @{$_[1]};
	my @stems = @{$_[2]};
	my $min_energy = $_[3];
	my $RNAseq = $_[4];
	my $rou = 0.2;
	my $gama = 0.6;
	for(my $i=0;$i<scalar(@stems);$i++){
		for(my $j=0;$j<scalar(@stems);$j++){
			$pheremones[$i][$j] = (1-$rou)*$pheremones[$i][$j];
		}
	}
	for(my $i=0;$i<scalar(@ants);$i++){
		my @solution = @{$ants[$i]};
		my $structure = Structure(\@solution, \@stems, length($RNAseq));
		my $energy = CompleteEnergy($structure, $RNAseq);
		for(my $m=0;$m<scalar(@solution);$m++){
			for(my $n=$m+1;$n<scalar(@solution);$n++){
				my $Q = $min_energy/$gama;
				$pheremones[$n][$m] += $energy/$Q;
				$pheremones[$m][$n] += $energy/$Q;
			}
		}
	}
	return @pheremones;
}

# check whether the Sub Solution Set Change or not
sub Change_Sub_Solution {
	my @new = @{$_[0]};
	my @old = @{$_[1]};
	my @stems = @{$_[2]};
	my $RNAseq = $_[3];
	if( scalar(@new) == scalar(@old) ){
		my $new_structure = Structure(\@new, \@stems, length($RNAseq));
		my $old_structure = Structure(\@old, \@stems, length($RNAseq));
		if( $new_structure eq $old_structure ){
			return 0;
		}
	}
	return 1;
}

sub print_result {
	my $main_structure = $_[0];
	my $best_Structure = $_[1];
	my $RNAseq = $_[2];
	my $IL = $_[3];
	my $j;
	open(my $out,'>>',"result.txt") or die "Can't open file for writing: $!"; 
	print $out $best_Structure," : best structure\n";
	print $out $main_structure," : main structure\n";
	print $out "Len : ",length($RNAseq),"\n";
	print $out "IL : ",$IL,"\n";
	my @Stack = ();
	my @MainStr = ();
	my $se2 = 0;
	my $ppv2 = 0;
	
	for(my $i=0;$i<length($RNAseq);++$i){
		if(substr($main_structure,$i,1) eq "."){
			$MainStr[$i] = $i;
		}
		elsif(substr($main_structure,$i,1) eq "("){
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
	for(my $i=0;$i<length($RNAseq);++$i){
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
	for(my $i=0;$i<length($RNAseq);++$i){
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