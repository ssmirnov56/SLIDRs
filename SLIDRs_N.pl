#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use Pod::Text;
##use Data::Dumper;


########################################################

my ($help, $verbose, $debug) = (0,0,0);

my ($fasta_file, $sites_file, $Nmutations, $exactNmutations, $context, $pattern, $report_all) = ('', '', 2, 0, 0, '',0);

my $args = {
	"h" => \$help,
	"v" => \$verbose,
	"d" => \$debug,
	"fasta:s"     	=> \$fasta_file,
	"sites:s"      	=> \$sites_file,
   	"Nmut:i"     	=> \$Nmutations, 
   	"exactNmut"   => \$exactNmutations, 
	"context:i"     => \$context,
        "context_pattern:s"     => \$pattern,
        "report_all"    => \$report_all,
};

GetOptions(%$args);


my $msg2usr='';

$msg2usr .= "\nProvide fasta filename; -fasta" if not $fasta_file;
$msg2usr .= "\nProvide sites filename; -sites" if not $sites_file;


if ($msg2usr or $help) {

warn "$msg2usr\n\n" if $msg2usr and not $help;

print "version N. The copyright/years updated.\n";
print "version M. Reports the total count of hits found.\n";
print "version K. Reports sites with user-defined patterns for the context; -context_pattern flag.\n";
print "version J. Very buggy version, it is the time to rewrite the code with proper structuring with subroutines.\n";
print "version I. Bug fixed for sites search under -exactNmut flag.\n";
print "version H. Minor typos in the manual fixed.\n";
print "version G. Software received name SLIDRs: Segmental Labeling for IDRs.\n";
print "version F. Reports adjacent N- and C- terminal polypetdies for sequence context analysis.\n";
print "version E. Can analyze multiple sites, listed in a file. E.g. LPXTA, LPXTFG, LPXTS for sortases.\n";
print "version D. Processes multi-fasta per file with analysis of DisProt in mind (~700 IDRs).\n"; 
print "version C. Reports IDR stats: #AA, #sites, #AA/site.\n";
print "version B. The number of mutations can be exact (with a flag) or up to a limit (default).\n";
print "version A. Processes a single-sequence fasta file only. The number of mutations is exact (min=max).\n";
print "\nPurpose:\n";
print "Reports sites fitting the pattern achievable with a given number of point mutations.\n";
print "\nUsage:\n";

print "\t$0 [-h] [-v] [-d] -fasta {fasta_filename} -site {target_site} [-Nmut {number of point mutations to fit the site}] [-exactNmut] [-context {value}] [-context_pattern {pattern}] [-report_all]\n";
print "Params:\n";
print "\t-fasta:\t\tProtein sequence fasta filename\n";
print "\t-sites:\t\tFilename for listing the target sites to search. \"X\"=any residue [$sites_file]\n";
print "\t-Nmut:\t\tMaximal number of mutations required for the match [$Nmutations]\n";
print "\t-exactNmut:\tIf activated, the number of mutations required for the match is an exact value [$exactNmutations]\n";
print "\t-context:\tIf positive [deafult $context], reports the oligopeptides of specified size on the N- and C- termini of every site found.\n";
print "\t-context_pattern:\tFlanking sequences are searched for specific patterns. E.g., use WYF for aromatics, KR for basic, ED for acidic residues. [$pattern]\n";
print "\t-report_all:\tReport stats on all sequences whether tirget sites/contexts are found or not\n";

print "(C):\n";
print "\tSerge Smirnov, John Antos, Western Washington University, smirnos\@wwu.edu, antosj\@wwu.edu, 2019-2024\n";
exit(0) if not $debug; 
}

report($args) if $verbose;


die "Fhheewww, made it here !" if (0);

#################################################
#
# M A I N    
#
#################################################


###########################################
##     Reading in the input files        ##
###########################################

my (@arr_fasta_seq, @arr_fasta_header, @arr_sites);

my ($site_size) = (0);



## Read in the one-letter aa seqeunce as $fasta_sequence. 

open (FASTA, $fasta_file) or die "could not open $fasta_file: $!";

my $curr_fasta_seq = -1;

while (<FASTA>) {
    chomp;


    if (/^\>/) {

    	push @arr_fasta_header, $_;

    	$curr_fasta_seq++;

    	$arr_fasta_seq[$curr_fasta_seq]='';

    	next
    } else {

        s/\s//g;

    	my $seq = $arr_fasta_seq[$curr_fasta_seq].$_;

	    $arr_fasta_seq[$curr_fasta_seq] = $seq;
    }

}
close (FASTA);


## Convert $site into @arr_site and get its size. 

open (SITES, $sites_file) or die "could not open $sites_file: $!";
while (<SITES>) {
    chomp;

    s/\s//g;

    push @arr_sites, $_;

    my @arr_site = split("", $_);

    if ($site_size == 0) {

        $site_size = scalar (@arr_site);
    } else {

        die "Target sites must be of the same size, otherwise face consequences and schmonsequences...\n" if ($site_size != scalar (@arr_site));

    }
}

close (SITES);


die "Site size: $site_size\n" if ($debug);


## Run through ALL FASTA entries
my $fastaNum = scalar (@arr_fasta_header);
my $siteNum = scalar (@arr_sites);

my $_hits_found = 0; 

for my $seqN ( 0..$fastaNum-1 ) {

## Convert $fasta_sequence into @arr_sequence and get its size. 

	my ($fasta_sequence, $sequence_description) = ($arr_fasta_seq[$seqN], $arr_fasta_header[$seqN]); 

	my @arr_sequence = split("", $fasta_sequence); 

	my $sequence_size = scalar (@arr_sequence);

    my ($sites_total) = (0);

	

## Run through the sequence and check for site pattern matches given the specified number of mutations

for my $i (0..$sequence_size - $site_size) {



## Try all possible sites
for my $siteN ( 0..$siteNum-1) {

    my $site = $arr_sites[$siteN];

    my @arr_site = split("", $site);



	my $mutations = 0;

	for my $j (0..$site_size-1) {

		my $seq_pos = $i+$j;

		next if ($arr_site[$j] eq 'X');

		next if ($arr_site[$j] eq $arr_sequence[$seq_pos]);

		$mutations++;

		last if ($mutations > $Nmutations);
	}

	next if ($mutations > $Nmutations);
	next if ( ($exactNmutations) and ($mutations != $Nmutations) );

	my ($site_found, $context_N, $context_C)= ('','',''); 

# Reading the target site found
	for my $j (0..$site_size-1) {

		$site_found = $site_found.$arr_sequence[$i+$j];
    }

#Reading the N- and C- terminal context polypeptides
    if ($context>0) {

    for my $j (0..$context-1) {

        $context_N = $context_N.$arr_sequence[$i-$context+$j] if ($i-$context+$j>-1);

        $context_C = $context_C.$arr_sequence[$i+$site_size-1+$j+1] if ($i+$site_size-1+$j+1< $sequence_size);
    }

    if ($pattern) {

        my (@a_context_N) = split("",$context_N);
        my (@a_context_C) = split("",$context_C);
        my ($size_N, $size_C, $hits_N, $hits_C) = (scalar (@a_context_N), scalar (@a_context_C), 0, 0);

        for my $k (0..$size_N-1) { $hits_N++ if ($pattern =~ /$a_context_N[$k]/) }
        for my $k (0..$size_C-1) { $hits_C++ if ($pattern =~ /$a_context_C[$k]/) }

        if ( ($report_all) or ($hits_N > $context-2) or ($hits_C > $context-2) ) {	### if report all no matter the pattern or pattern foun
											### report the site, skip testing other sites, move to the next pos in sequence
        	my $pos2report = $i+1;
        	print "pos= $pos2report\tsite= $site_found\tmutations= $mutations\tTarget_site= $site\tN_context= $context_N\tC_context= $context_C\n";
          	$sites_total++;
		last;
	
        }


    } else {		### context requested but no specific context patterns are needed; report all suitable-size context for stite found
			### skip testing other sites, move to the next position in the sequence

    	my $pos2report = $i+1;
        print "pos= $pos2report\tsite= $site_found\tmutations= $mutations\tTarget_site= $site\tN_context= $context_N\tC_context= $context_C\n";
    	$sites_total++;
	last;
    }


    } else {	### no context requested, report the site found, skip testing other sites, move to the next position in the sequence

    	my $pos2report = $i+1;
    	print "pos= $pos2report\tsite= $site_found\tmutations= $mutations\tTarget_site= $site\n";
    	$sites_total++;
    	last;
    }
}


}

if ($sites_total) {

    my $aa_per_site = $sequence_size/$sites_total;


    if ($aa_per_site =~ /^(\d+)\./) { $aa_per_site = $1 }


    print "Sequence_size= $sequence_size\tSites_found= $sites_total\tAA/site= $aa_per_site\tSEQ: $sequence_description\n"; 
	$_hits_found = $_hits_found + $sites_total;
    next;
} else {

    print "Sequence_size= $sequence_size\tNo_sites_found_for \tSEQ: $sequence_description\n" if ($report_all);
}


}

warn "Sites found: $_hits_found\n";


##print "hProt\n" if $verbose;

##print Data::Dumper->new([\%hProt])->Dump,"\n" if $verbose;

##warn Data::Dumper->new([$aliph->{'VAL'}])->Dump,"\n";


##print Data::Dumper->new([$rh])->Dump,"\n" if 1;

##warn Data::Dumper->new([\%best_pos])->Dump,"\n" if $verbose;
##warn Data::Dumper->new([\%best_misses])->Dump,"\n" if $verbose;

##warn Data::Dumper->new([\%BMRB_CS])->Dump,"\n";
##warn Data::Dumper->new([\%BMRB_SD])->Dump,"\n";


##print STDERR "Delta\tMisses\tOne\tTwo\tMany\tAll\n";
##print "$peak_delta\t$misses\t$ones\t$twos\t$many\t$total\n";





#################################################
#
# S U B s
#
#################################################


sub report { 
	my ($rh_args) = @_;
	
	foreach my $key (keys %$rh_args) {

		print "$key => ",${$rh_args->{$key}}, "\n";
	}
}

