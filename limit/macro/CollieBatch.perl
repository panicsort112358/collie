#!/usr/bin/perl

use Cwd;

# Usage for CollieBatch.perl
#
# ./CollieBatch.perl <execName> <massPoint> <outputName>
#
# execName   = executable name without ".exe"
# massPoint  = integer mass point to study, or -1 for all
# outputName = base name for "outputName.root" and "outputName.log"
#
# Change this path to your job directory.
$colliePath = "../";
$collieSetup = $colliePath . "setup_Collie.tcsh";

die "Improper input!" if ($#ARGV < 2);

$exe = $ARGV[0];
$mass = $ARGV[1];
$job = $ARGV[2];

$workdir = getcwd;

$output = $exe;
$output =~ s/.exe//;
$output = $output . "_" . $job . "_" . $mass;

$log = $output . ".log";
$output = $output . ".root";

$pfile = "submit_" . $job . "_" . $mass . ".sh";

print "\n***********************************\n";
print "    Submitting Collie Batch Job\n";
print "Executable: $exe\n";
print "Output: $output\n";
print "Logfile: $log\n";
print "Script: $pfile\n";
print "***********************************\n";

open(FILE, ">$pfile");
print FILE <<EOS;
#!/bin/sh
#CLUED0 -l cput=72:00:00
#CLUED0 -l mem=450mb

echo `/bin/hostname`

. /etc/bashrc

pwd

cd $workdir

pwd

source $collieSetup

set WORKDIR="NONE"
if [ -r /usr/local/etc/pbssetup.sh ]; then
    . /usr/local/etc/pbssetup.sh
    echo "Found pbssetup.sh"
    echo "Workdir is now:"
    echo \$WORKDIR
else
    echo "ERROR: Could not find /usr/local/etc/pbssetup.sh"
    echo "ERROR: Maybe the partition is not mounted correctly"
    echo "ERROR: Exiting..."
    exit 1
fi
if [ \$WORKDIR = "NONE" ]; then
    echo "ERROR: WORKDIR is not defined, problem with batch system!"
    echo "ERROR: Exiting..."
    exit 1
fi
if [ \$WORKDIR = \$HOME ]; then
    echo "ERROR: WORKDIR is defined to be your homedir"
    echo "ERROR: This means there is a problem with setting WORKDIR!"
    echo "ERROR: Exiting..."
    exit 1
fi

./$exe $output $mass > $log

EOS
close(FILE);

system("/usr/local/bin/cluesow", $pfile);
system("sleep 1");
