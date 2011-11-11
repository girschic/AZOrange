#!/usr/bin/perl -w

###############################################################################
# sdf2gspan.pl converts molecular data from SDF to the gSpan input format.
#
# parameter: <SDF file>
###############################################################################


if(@ARGV == 0){
	die "parameter: <SDF file>\n";
}

my $sdf_file = $ARGV[0];
my $gsp_file = $sdf_file;
if($gsp_file =~ /\..+/){
	$gsp_file =~ s/\..+/.gsp/;
}else{
	$gsp_file = $gsp_file.".gsp";
}
my $conversion_file = $sdf_file;
if($conversion_file =~ /\..+/){
	$conversion_file =~ s/\..+/.convert/;
}else{
	$conversion_file = $conversion_file.".convert";
}

open(SDF_FILE,  $sdf_file)   or die "Can't open $sdf_file: $!";
open(GSP_FILE, ">".$gsp_file)   or die "Can't open $gsp_file: $!";


$/ = "\$\$\$\$\n";
my @lines;
my %atom2number;
my @bond2number;
my @true_bond2number;
my @number2id;
my @number2atom;
my @number2bond;
my @true_number2bond;
my $current_atom_no = 0;
my $current_bond_no = 0;
my $number_atoms;
my $number_bonds;
my $start_index;
my $atom;
my $atom_count;
my $id;
my $current_id_no = 0;

while(<SDF_FILE>){
	@lines = split(/\n/, $_);
	if(@lines > 3){
		$id = $lines[0];
		$number2id[$current_id_no] = $id;
		print GSP_FILE "t # ".$current_id_no."\n";
		# #vertices, #edges
		$lines[3] =~ /^ +(\d+) +(\d+) /;
		# substring (start at position 0, length 3)
		$atoms_temp = substr($lines[3],0,3);
		$atoms_temp =~ /(\d+)/;
		# $1: The number variables are the matches from the last successful match or substitution operator applied
		$number_atoms = $1;
		$bonds_temp = substr($lines[3],3,3);
		$bonds_temp =~ /(\d+)/;
		$number_bonds = $1;
		$start_index = 4;
		$atom_count = 0;
		for(my $i=$start_index; $i<$start_index+$number_atoms; $i++){
			$atom = substr($lines[$i],31,3);
			$atom =~ s/\W//g;
			if(!(exists($atom2number{$atom}))){
				$atom2number{$atom} = $current_atom_no;
				$number2atom[$current_atom_no] = $atom;
				$current_atom_no ++;
			}
			print GSP_FILE "v ".$atom_count." ".$atom2number{$atom}."\n";
			$atom_count++;
		}
		$start_index = $start_index + $number_atoms;
		for(my $i=$start_index; $i<$start_index+$number_bonds; $i++){
			$temp = substr($lines[$i],0,3);
			$temp =~ /(\d+)/;
			$atom1 = $1;
			$temp = substr($lines[$i],3,3);
			$temp =~ /(\d+)/;
			$atom2 = $1;
			$temp = substr($lines[$i],6,3);
			$temp =~ /(\d+)/;
			$bond_type = $1;

			# added by Madeleine Seeland, 09.02.10
			# gSpan expects that there is no gap in the bond labels
			# since the bond label for aromating rings is 3, there will 
			# be a problem if there is e.g. no graph with bond label 2
			# thus this preprocessing script is changed such that the 
			# bond labels are in ascending order and there is no gap between 
			# bond labels
			if(!(exists($bond2number[$temp]))){
				$bond2number[$temp] = $current_bond_no;
				$number2bond[$current_bond_no] = $temp;
				$current_bond_no ++;
			}

#			if(!(exists($true_bond2number[$temp]))){
#				$true_bond2number[$temp] = $temp;
#				$true_number2bond[$temp] = $temp;
#			}

			#$bond_type = $1;
			$bond_type = $bond2number[$temp];
			print GSP_FILE "e ".($atom1-1)." ".($atom2-1)." ".($bond_type)."\n";
			#if($id =~ /50-76-0/){
			#	print $id."e ".($atom1-1)." ".($atom2-1)." ".($bond_type-1)."\n";
			#}
		}
		$current_id_no++;
	}
}

close GSP_FILE;
close SDF_FILE;

open(CONV_FILE, ">".$conversion_file)   or die "Can't open $conversion_file: $!";

print CONV_FILE "molecule ids - graph numbers:\n";
for(my $i=0; $i<@number2id; $i++){
	print CONV_FILE $number2id[$i]."\t".$i."\n";
}
print CONV_FILE "atoms - atom numbers:\n";
for(my $i=0; $i<@number2atom; $i++){
#	if(length($number2atom[$i])>1){
#		print CONV_FILE "[".$number2atom[$i]."]"."\t".$i."\n";
#	}
#	else{
	print CONV_FILE $number2atom[$i]."\t".$i."\n";
#	}
}


print CONV_FILE "edge labels:\n";
for(my $i=0; $i<@number2bond; $i++){
#	if(length($number2atom[$i])>1){
#		print CONV_FILE "[".$number2atom[$i]."]"."\t".$i."\n";
#	}
#	else{
	print CONV_FILE $number2bond[$i]."\t".$i."\n";
#	}
}

close CONV_FILE;
