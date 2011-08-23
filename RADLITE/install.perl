#!/usr/bin/perl
#=======================================================
#           INSTALL ROUTINE FOR "RADICAL"
#=======================================================
$pwd  = `pwd` ; 
chop($pwd) ;
print "Installing RADICAL in the directory:\n  '$pwd'\n";
$home = $ENV{"HOME"};
$bin  = $home . "/.system/bin" ;
if(!(-e $bin)) {
    $bin  = $home . "/bin" ;
}
if(!(-e $bin)) {
    print "You must have a bin/ directory in your home directory\n" ;
    print "-----Shall I make $bin for you?\n" ;
    $input = <STDIN> ;
    print $input ;
    if($input=~/^[yY]/) {
	print "Creating $bin for you...\n" ;
	system("mkdir $bin") ;
    }
}
$path = $ENV{"PATH"};
if(!($path =~ /$bin/)) {
    print "The $bin directory exists, but it not in the PATH environment variable\n" ;
    print "You must put the \n" ;
    print "$bin directory \n" ;
    print "in the path yourself (in the .tcshrc file or so...)\n" ;
    print "If you do it now, don't forget to type 'rehash'.\n" ;
}
#
# Creating the automatic link for the executable
#
print "  Creating a link 'radical' in '$bin/'\n" ;
$radical    = $pwd . "/radical" ;
$radicallnk = $bin . "/radical" ;
if(!(-e $radicallnk)) {
    print "------ Warning: file $radicallnk did not exist previously. You might want to type 'rehash'\n" ;
}
open(FILE,">$radicallnk") || die "Could not open log file\n" ;
#print FILE "#!/bin/sh\n" ;
#print FILE "$radical" ;
print FILE "#!/usr/bin/perl\n" ;
print FILE "system(\"$radical \@ARGV\");" ;
close (FILE) ;
`chmod u+rwx $radicallnk` ;
#
# Creating the backwards compatibility link to schar
#
$radical    = $pwd . "/radical" ;
$radicallnk = $bin . "/schar" ;
if(!(-e $radicallnk)) {
    print "------ Warning: file $radicallnk did not exist previously. You might want to type 'rehash'\n" ;
}
open(FILE,">$radicallnk") || die "Could not open log file\n" ;
#print FILE "#!/bin/sh\n" ;
#print FILE "$radical" ;
print FILE "#!/usr/bin/perl\n" ;
print FILE "system(\"$radical \@ARGV\");" ;
close (FILE) ;
`chmod u+rwx $radicallnk` ;
#
# Creating the local backwards compatibility link to schar
#
$radical    = $pwd . "/radical" ;
$radicallnk = $pwd . "/schar" ;
open(FILE,">$radicallnk") || die "Could not open log file\n" ;
#print FILE "#!/bin/sh\n" ;
#print FILE "$radical" ;
print FILE "#!/usr/bin/perl\n" ;
print FILE "system(\"$radical \@ARGV\");" ;
close (FILE) ;
`chmod u+rwx $radicallnk` ;
#
# Creating the automatic link for the makesigma
#
print "  Creating a link 'makesigma' in '$bin/'\n" ;
$makesigma    = $pwd . "/makesigma" ;
$makesigmalnk = $bin . "/makesigma" ;
if(!(-e $makesigmalnk)) {
    print "------ Warning: file $makesigmalnk did not exist previously. You might want to type 'rehash'\n" ;
}
open(FILE,">$makesigmalnk") || die "Could not open log file\n" ;
#print FILE "#!/bin/sh\n" ;
#print FILE "$makesigma" ;
print FILE "#!/usr/bin/perl\n" ;
print FILE "system(\"$makesigma \@ARGV\");" ;
close (FILE) ;
`chmod u+rwx $makesigmalnk` ;
#
# Creating the automatic link for the matrixrun
#
print "  Creating a link 'matrixrun' in '$bin/'\n" ;
$matrixrun    = $pwd . "/matrixrun.perl" ;
$matrixrunlnk = $bin . "/matrixrun" ;
if(!(-e $matrixrunlnk)) {
    print "------ Warning: file $matrixrunlnk did not exist previously. You might want to type 'rehash'\n" ;
}
open(FILE,">$matrixrunlnk") || die "Could not open log file\n" ;
#print FILE "#!/bin/sh\n" ;
#print FILE "$matrixrun" ;
print FILE "#!/usr/bin/perl\n" ;
print FILE "system(\"$matrixrun \@ARGV\");" ;
close (FILE) ;
`chmod u+rwx $matrixrunlnk` ;

