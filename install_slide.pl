#!/usr/bin/perl
# ^^^^^^^^^ 
# CUSTOMIZE THIS TO POINT TO YOUR PERL BINARY!!!!
#
# SLIDE installation script   Volker Schnecke   Fri Sep 17 15:43:08 EDT 1999
#
# to run: ./install_slide.pl

if ( $ENV{"SLIDE_DIR"} eq "" )
{
    print "set environment variable SLIDE_DIR to point to this directory\n";
    exit;
}


$slide_dir = $ENV{"SLIDE_DIR"};
$perl = $^X;

mkdir "$slide_dir";
mkdir "$slide_dir/bin";

# unpack the Perl scripts and change the path to the perl binary
open IN, "$slide_dir/src/scripts/perl_scripts";
print "\n*** extracting Perl scripts ***\n";
while ( <IN> )
{
    if ( /^____FILE\_START/ )
    {
	chop;
	# grep the name of packed Perl script
	@line = split;
	$file = "$slide_dir/bin/$line[1]";
	open OUT, ">$file";
	# read the line that includes the path to the Perl binary
	$_ = <IN>;
	print OUT "\#\!$perl\n";
	next;
    }
    if ( /^____FILE\_END/ )
    {
	# end of this script, close filehandle and change permissions
	close OUT;
	chmod 0755, "$file";
	next;
    }
    print OUT $_;
}
close IN;

# unpack the shell scripts
open IN, "$slide_dir/src/scripts/shell_scripts";
print "\n*** extracting Shell scripts ***\n";
while ( <IN> )
{
    if ( /^____FILE\_START/ )
    {
	chop;
	# grep the name of packed Perl script
	@line = split;
	$file = "$slide_dir/bin/$line[1]";
	open OUT, ">$file";
	next;
    }
    if ( /^____FILE\_END/ )
    {
	# end of this script, close filehandle and change permissions
	close OUT;
	chmod 0755, "$file";
	next;
    }
    print OUT $_;
}
close IN;

# compile the libSLIDE_utils first since the other programs will need to
# link to it
chdir "src/utils";
print "\n*** compiling the static slide utilties library ***\n\n";
system "make";

# compile SLIDE
chdir "../slide";
print "\n*** compiling SLIDE ***\n\n";
system "make";
# compile compute_interaction_centers, check_connectivity, and
# generate_rasmol_script
chdir "../interactions";
print "\n*** compiling interactions auxiliaries ***\n\n";
system "make";
# compile average_template and unbiased_template
chdir "../template";
print "\n*** compiling template auxiliaries ***\n\n";
system "make";
print "\n*** SLIDE installation completed ***\n";
print "\n\n*** Please see license.txt for licensing and use information ***\n\n";
