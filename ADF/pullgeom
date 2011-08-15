#!/usr/bin/perl

#use warnings;

# Gets the last geometry from an ADF output file and prints it to the screen in cleaned format.
$step_to_take = 0;
if ("$ARGV[0]" eq "-n")
{
	shift(@ARGV);
	$step_to_take=$ARGV[0];
	shift(@ARGV);
}

if ($ARGV[0] eq "")
{
	print("Usage: pullgeom [-n step] file\n");
	exit; 
}

#Open specified file for input.
open(IN, "<$ARGV[0]");

$out_string = "";
$in_cart_section = 0;
$in_geom_section = 0;
$atom_count = 0;
$geometry_count = 0;

# Loop over every line in the file.
while (<IN>) {
	$this_line_is_separator=0;
	if ( $_ =~ /^ Coordinates \(Cartesian\)/ )
	{
		#print("Found Cartesian Coords\n");
		$in_cart_section = 1;
		$out_string = "";
		$geometry_count = $geometry_count + 1;
	}	
	
	if ( $_ =~ /^ -----([-]+)/ )
	{
		$this_line_is_separator=1;
		
		if ( $in_cart_section == 1 )
		{
			#print("Found separator\n");
			if ( $in_geom_section == 0 )
			{
				$in_geom_section = 1;
			} else {
				$in_geom_section = 0;
				$in_cart_section = 0;

				#If a step count was specified, and this is it, then stop.
				if (($step_to_take != 0) and ($step_to_take == $geometry_count))
				{
					last;
				}	
			}
		}
	}

  # Can just capture all the info, bohr distances and all, but this is not preferable.
	#if ( ( $in_geom_section == 1 ) and ( $this_line_is_separator == 0 ) )
	#{
	#	$out_string=$out_string.$_;
	#}

  # Capture atom position descriptions
	if ( $_ =~ /^(?: +)(\d+) (\w{1,2})(?: +)(?:[0-9.-]+)(?: +)(?:[0-9.-]+)(?: +)(?:[0-9.-]+)(?: +)([0-9.-]+)(?: +)([0-9.-]+)(?: +)([0-9.-]+)/ )
	{
	  # print("$_"); # (Can check whether the regexp is working.)
		if ( $in_geom_section == 1 )
		{
			$atoms[$1][0]=$2;
			$atoms[$1][1]=$3;
			$atoms[$1][2]=$4;
			$atoms[$1][3]=$5;
			if ( $atom_count < $1 ) { $atom_count = $1; }
		}
 	}		

}

$out_string="";
$i=1;
while( $i<( $atom_count + 1 ) )
{
	$j=0;
	while($j<4)
	{
		if ( $atoms[$i][$j] =~ /(^.)/ )
		{
		  #Pad positive numbers
			if ( "$1" =~ /(^\d)/ )
			{
				$out_string=$out_string."  ".$atoms[$i][$j];
			} else {
			  $out_string=$out_string." ".$atoms[$i][$j];
			}

			#Pad elements with only one letter.
		  if ( ( $j==0 ) and ( length($atoms[$i][$j])<2 ) )
			{
				$out_string=$out_string." ";
			}
		}	
		$j=$j+1;
	}
	$out_string=$out_string."\n";
	$i=$i+1;
}

#print("$atom_count\n");
print("$#atoms\n_\n");
print("$out_string");

