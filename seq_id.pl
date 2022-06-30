#!/usr/bin/perl -w

use strict;
use IO::String;
use Bio::Seq; 
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::Align::AlignI;
use File::Find::Object::Rule;
use POSIX;
use Statistics::Descriptive;

my $input = "";
my $full_out = 0;
my $window = 0;

if ($#ARGV == -1) {
	print "\n";
	print "This script calculates average sequence identity statistics for an alignment.\n";
	print "\nSyntax:\n";
	print " -f X		use file X as input ('all' to use all files in the current dir)\n";
	print " -full		print pairwise differences and their positions\n";
	print " -full_window X	print pairwise differences and their positions in X nucleotide windows\n";
	print "\n";
	exit;
}

foreach my $argnum (0 .. $#ARGV) {
	if ($ARGV[$argnum] eq "-f" and $ARGV[$argnum+1]) {
		$input = $ARGV[$argnum+1];
		if (-e $input) {
			print "\nLoading $input...\n";
		}
		elsif ($input ne "all") {
			$input="";
		}
	}
	if ($ARGV[$argnum] eq "-full") {
		$full_out = 1;
	}
	if ($ARGV[$argnum] eq "-full_window") {
		$full_out = 1;
		$window = $ARGV[$argnum+1];
	}	
}

if ($input) {
	print "\n\n";
	if ($input eq "all") {
		open(LOG, ">multifile_id_log.txt");
		my $ffor = File::Find::Object::Rule->file()->name("*.fas");
		$ffor->maxdepth(1);
		my @filelist = $ffor->in(".");
		foreach my $file (@filelist) {
			calc_seq_id($file);
		}
		close LOG;
	}
	elsif (-e $input) {
		open(LOG, ">$input\_id_log.txt");
		calc_seq_id($input);
		close LOG;
	}
}
else {
	print "\nNo file specified or file does not exists. Pass a filename after the \"-f\" parameter (or 'all' to use all files in the current dir)\n";
}

exit();

sub calc_seq_id {
	my $seq_source = $_[0];
	my $seq_source_short = substr($seq_source,0,length($seq_source)-4);

	print "Loading sequences from $seq_source...\n";
	print LOG "Loading sequences from $seq_source...";

	my $seq_in = Bio::AlignIO->new(-format => 'fasta', -file => $seq_source);
	my $seqfile = $seq_in->next_aln;
	my $seq_count = $seqfile->num_sequences;
	#my $total_id=0;
	my $stat = Statistics::Descriptive::Full->new();
	my @full_stat = ();
	
	my $id = 0;
	my $header_line=" ";
	my @pw_data = ();
	for (my $i1 = 1; $i1 <= $seq_count; $i1++) {
		print "Calculating identity values for sequence $i1\/$seq_count...";
		my $seq1 = $seqfile->get_seq_by_pos($i1);
		my $seq1_id = $seq1->id;
		if ($i1 > 1) {
			$header_line = $header_line.("\t$seq1_id");
		}
		push(@full_stat,$seq1_id);
		$full_stat[$i1-1] = $full_stat[$i1-1].(("\t")x($i1-1));

		for (my $i2 = 1; $i2 <= $seqfile->num_sequences; $i2++) {
			if ($i2 < $i1) {
				next;
			}
			if ($i1 == $i2) {
				$id = 100;
				next;
			}
			if ($i2 > $i1) {
				my $seq2 = $seqfile->get_seq_by_pos($i2);
				my $seq2_id = $seq2->id;
				my $pos_list;
				($id, $pos_list) = id_calc($seq1,$seq2,$full_out);
				if ($full_out) {
					print "\n$seq1_id and $seq2_id differences\:\n";				
					print LOG "\n$seq1_id and $seq2_id differences\:\n";					
					if ($window > 0) {
						my $count=0;
						my $total_count=0;
						my $pos_list_size = @$pos_list;
						for (my $i1=0,my $i2=0; $i1 < $pos_list_size; $i1++) {
							if (@$pos_list[$i1] <= ($i2+1)*$window) {
								$count++;
							}
							if (@$pos_list[$i1] > ($i2+1)*$window and $i1+1 == $pos_list_size) {
								$total_count += $count+1;
								print "".(($i2)*$window+1)."-".(($i2+1)*$window).":\t$count\n";
								print LOG "".(($i2)*$window+1)."-".(($i2+1)*$window).":\t$count\n";
								print "".(($i2+1)*$window+1)."-".(($i2+2)*$window).":\t1\n";
								print LOG "".(($i2+1)*$window+1)."-".(($i2+2)*$window).":\t1\n";
								last;
							}
							if (@$pos_list[$i1] > ($i2+1)*$window or $i1+1 == $pos_list_size) {
								$total_count += $count;
								print "".(($i2)*$window+1)."-".(($i2+1)*$window).":\t$count\n";
								print LOG "".(($i2)*$window+1)."-".(($i2+1)*$window).":\t$count\n";
								if ($i1+1 == $pos_list_size) {
									last;
								}
								$count = 0;							
								$i1--;
								$i2++;
							}
						}
						print "Total:\t$total_count\/$pos_list_size\n";
						print LOG "Total:\t$total_count\/$pos_list_size\n";
					}
					else {
						print "@$pos_list\n";
						print LOG "@$pos_list\n";
						$full_stat[$i1-1] = $full_stat[$i1-1].("\t$id");
					}
				}
			}
			if ($full_out == 0) {
				$stat->add_data($id); 
				#$total_id = $total_id+$id;
			}
		}
		print "done\n";
	}
	if ($window == 0) {
		if ($full_out) {
			print "\n";	
			print LOG "\n";	
			foreach my $i1 (@full_stat) {
				print "$i1\n";
				print LOG "$i1\n";
			}
			print "$header_line\n";
			print LOG "$header_line\n";
		}
		else {
			my $mean_id = $stat->mean();
			my $std_id  = $stat->standard_deviation();

			my $final_mean = sprintf("%.2f", $mean_id);
			my $final_std = sprintf("%.2f", $std_id);

	#		my @aaa = $stat->get_data();
	#		print "\n";
	#		foreach my $a1 (@aaa) {
	#			print "$a1\n";
	#		}

			print "\n";
			print "Mean identity value: $final_mean +/- $final_std\n";
			print LOG "mean identity value: $final_mean +/- $final_std\n";
			print "\n";
		}
	}
	print "\n";
#	my $final_mean = (floor(sprintf("%.2f", $mean_id*100)))/100;
#	my $final_mean = (floor(sprintf("%.2f", $std_id*100)))/100;
	
	#TRIM/ROUND VALUES TO 2 DECIMAL PLACES	
#	my $final_id = $total_id/($seq_count*$seq_count);
#	$final_id = (floor(sprintf("%.2f", $final_id*100)))/100;
#	print "Mean identity value: $final_id\n\n";
#	print LOG "Mean identity value: $final_id\n\n";
}


sub id_calc {
	#TRANSFORM SEQUENCES INTO STRINGS
	my $seq1 = $_[0];
	my $seq2 = $_[1];
	my $return_count = $_[2];
	my @pos = ();
	my $str1;
	my $str2;
	my $seq_str1 = IO::String->new(\$str1);
	my $seq_str2 = IO::String->new(\$str2);
	my $seqOut1 = Bio::SeqIO->new(-format => 'fasta', -fh => $seq_str1 );
	my $seqOut2 = Bio::SeqIO->new(-format => 'fasta', -fh => $seq_str2 );
	$seqOut1->write_seq($seq1);
	$seqOut2->write_seq($seq2);

	#REMOVE THE FASTA HEADER AND WHITESPACES, CONVERT TILDES TO HYPHENS
	$str1 = substr($str1,index($str1,"\n")+1);
	$str1 =~ s/[\n\r\s]+//g;
	$str1 =~ s/~/-/g;

	$str2 = substr($str2,index($str2,"\n")+1);
	$str2 =~ s/[\n\r\s]+//g;
	$str2 =~ s/~/-/g;	
			
	my $count = 0;
	my @s1 = split(//,$str1);
	my @s2 = split(//,$str2);
	my $l = length ($str1);
	my $l_true = 0;

	#OLD WAY: 
	#COUNT THE TRUE LENGTH (EXCLUDE GAPS)
#	for (my $i = 0; $i < $l;$i++) {
#		if ($s1[$i] ne $s2[$i]) {
#			$l_true++;
#		}
#		elsif ($s1[$i] ne "-") {
#			$l_true++;
#		}
#	}
	#NEW WAY
	#COUNT THE TRUTH LENGTH AND DIFFERENCES 
	for (my $i = 0; $i < $l;$i++) {
		if ($s1[$i] ne $s2[$i]) {
			$count++;
			push(@pos,$i+1);
			$l_true++;
		}
		elsif ($s1[$i] ne "-") {
			$l_true++;
		}
	}
	#foreach my $pos1 (@pos) {
		#print "$pos1\n";
	#}
	if ($return_count) {
		return ($count, \@pos);
	}
	else {
		my $id = (1-($count/$l_true))*100;

		#TRIM TO TWO DECIMAL PLACES
		my $id2 = (floor(sprintf("%.1f", $id*10)))/10;

		return ($id2, \@pos);
	}
}




