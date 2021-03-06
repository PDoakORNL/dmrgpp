#!/usr/bin/perl

use strict;
use warnings;
use Utils;

my ($omega0,$total,$omegaStep,$observable,$parallel) = @ARGV;
defined($parallel) or die "USAGE: $0 omegaBegin omegaTotal omegaStep observable nobatch | submit | test\n";

my $templateInput = "inputTemplate.inp";
my $templateBatch = "batchTemplate.pbs";

for (my $i = 0; $i < $total; ++$i) {
	my $omega = $omega0 + $omegaStep * $i;
	print STDERR "$0: About to run for omega = $omega\n";
	runThisOmega($i,$omega,$parallel);
	print STDERR "$0: Finished         omega = $omega\n";
}

sub runThisOmega
{
	my ($ind,$omega,$parallel) = @_;
	my $n = Utils::getLabel($templateInput,"TotalNumberOfSites=");
	my $input = createInput($n,$ind,$omega);
	if ($parallel eq "nobatch") {
		system("./dmrg -f $input -o ':$observable.txt' &> out$ind.txt");
		system("echo '#omega=$omega' >> out$ind.txt");
	} else {
		my $batch = createBatch($ind,$omega,$input);
		submitBatch($batch) if ($parallel eq "submit");
	}
}

sub createInput
{
	my ($n,$ind,$omega)=@_;
	my $file="input$ind.inp";
	open(FOUT,">$file") or die "$0: Cannot write to $file\n";
	my $steps = int($n/2) - 1;
	my $data = "data$ind.txt";
	my $nup = int($n/2);
	my $ndown = $nup;
	my $U = Utils::getLabel($templateInput,"##U=");
        my $hubbardU = setVector($n,$U);
        my $V = Utils::getLabel($templateInput,"##V=");
        my $potentialV = setVector(2*$n,$V);

	open(FILE,"$templateInput") or die "$0: Cannot open $templateInput: $!\n";

	while(<FILE>) {
		next if (/^#/);
		if (/\$([a-zA-Z0-9\[\]]+)/) {
				my $name = $1;
				my $str = "\$".$name;
				my $val = eval "$str";
				defined($val) or die "$0: Undefined substitution for $name\n";
				s/\$\Q$name/$val/g;
		}
		print FOUT;
	}

	close(FILE);
	close(FOUT);

	return $file;
}

sub createBatch
{
        my ($ind,$omega,$input) = @_;
        my $file = "Batch$ind.pbs";
        open(FOUT,">$file") or die "$0: Cannot write to $file: $!\n";

        open(FILE,"$templateBatch") or die "$0: Cannot open $templateBatch: $!\n";

        while(<FILE>) {
                while (/\$\$([a-zA-Z0-9\[\]]+)/) {
                        my $line = $_;
                        my $name = $1;
                        my $str = "\$".$name;
                        my $val = eval "$str";
                        defined($val) or die "$0: Undefined substitution for $name\n";
                        $line =~ s/\$\$$name/$val/;
                        $_ = $line;
                }

                print FOUT;
        }

        close(FILE);
        close(FOUT);

        print STDERR "$0: $file written\n";
        return $file;
}

sub submitBatch
{
        my ($batch) = @_;
        sleep(2);
        system("qsub $batch");
        print STDERR "$0: Submitted $batch\n";
}

sub setVector
{
        my ($sites,$U)=@_;
        my $tmp = "$sites ";
        for (my $i=0;$i<$sites;$i++) {
                $tmp .= " $U ";
        }
        return $tmp;
}

