#!/usr/bin/perl
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with
#  permission. For more information, see
#        http://csem.engin.umich.edu/tools/swmf

# Allow in-place editing
$^I = "";

# Add local directory to search
push @INC, ".";

use strict;

our $Component       = "IE";
our $Code            = "SIM";
our $MakefileDefOrig = 'src/Makefile.def';
our @Arguments       = @ARGV;

my $config = "share/Scripts/Config.pl";

# Git clone missing directories as needed
my $GITCLONE = "git clone"; my $GITDIR = "git\@gitlab.umich.edu:swmf_software";
if (not -f $config and not -f "../../$config"){
	print "Cloning share and util directories from herot ...";
    `$GITCLONE $GITDIR/share; $GITCLONE $GITDIR/util`;
}

if(-f $config){
    require $config;
}else{
    require "../../$config";
}

# Variables inherited from share/Scripts/Config.pl
our %Remaining; # Unprocessed arguments
our $ERROR;
our $WARNING;
our $Help;
our $Verbose;
our $Show;
our $ShowGridSize;
our $NewGridSize;

&print_help if $Help;

foreach (@Arguments){
    if(/^-s$/)                {$Show=1;  next};
    warn "WARNING: Unknown flag $_\n" if $Remaining{$_};
}

# Grid size variables
my $NameFile="src/SIM_main.f90";
my $GridSize;
my ($nTheta, $nPhi);

&get_settings;

&set_grid_size if $NewGridSize and $NewGridSize ne $GridSize;

&show_settings if $Show;

print "Config.pl -g=$nTheta,$nPhi\n" if $ShowGridSize and not $Show;


exit 0;

#############################################################################
sub print_help{

print "Set grid size to nTheta=91, nPhi=361 (1 deg by 1 deg):

    Config.pl -g=91,361

Show settings specific to IE/SIMPLE:

    Config.pl -s

";
    exit 0;
}

#############################################################################
sub get_settings{

    # Read size of the grid from $NameFile
    open(FILE,$NameFile) or die "$ERROR could not open $NameFile\n";
    while(<FILE>){
	next if /^\s*!/; # skip commented out lines
	$nTheta=$1 if /\bnTheta\s*=\s*(\d+)/i;
        $nPhi  =$1 if /\bnPhi\s*=\s*(\d+)/i;
	last if $nTheta and $nPhi;
    }
    close FILE;

    die "$ERROR could not read nTheta from $NameFile\n" 
	unless length($nTheta);
    die "$ERROR could not read nPhi from $NameFile\n" 
	unless length($nPhi);

    $GridSize = "$nTheta,$nPhi";
}

#############################################################################
sub set_grid_size{

    $GridSize = $NewGridSize;

    if($GridSize=~/^\d+,\d+$/){
	($nTheta, $nPhi)= split(',', $GridSize);
    }else{
	die "$ERROR -g=$GridSize should be 2 integers separated with a comma\n"
    }

    print "Writing new grid size into $NameFile...\n";

    @ARGV = ($NameFile);

    while(<>){
	if(/^\s*!/){print; next} # Skip commented out lines
	s/\b(nTheta\s*=[^\d]*)(\d+)/$1$nTheta/i;
	s/\b(nPhi\s*=[^\d]*)(\d+)/$1$nPhi/i;
	print;
    }

}

##############################################################################
sub show_settings{
	print "Grid                   : Spherical 2D\n";
    print "Number of colatitudes  : nTheta = $nTheta\n";
    print "Number of longitudes   : nPhi   = $nPhi\n";
}