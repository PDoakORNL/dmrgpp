#!/usr/bin/perl

use strict;
use warnings;
use utf8;

# TODO
# (1) Plot

my ($templateInput, $templateTex) = @ARGV;
defined($templateInput) or die "USAGE: $0 templateInput [templateTex]\n";

my %params = readPreParams($templateInput);
my $impurity = impurityInit($templateInput, \%params);

createInput("testImpurity.ain", $impurity, $templateInput);

if (defined($templateTex)) {
	my $plot = createPlot(\%params);
	$impurity->{plot} = $plot;
	createInput("test.tex", $impurity, $templateTex);
} else {
	print STDERR "$0: No templateTex given, no tex created\n";
}

sub plotLine
{
	my ($ix, $jx, $iy, $jy, $style) = @_;
	my $cix = sprintf("%.2f", $ix);
	my $ciy = sprintf("%.2f", $iy);
	my $cjx = sprintf("%.2f", $jx);
	my $cjy = sprintf("%.2f", $jy);
	return "\\draw[$style] ($cix, $ciy) -- ($cjx, $cjy);\n"
}

sub plotNode
{
	my ($ix, $iy, $radius, $label, $withNumbers) = @_;
	my $cix = sprintf("%.2f", $ix);
	my $ciy = sprintf("%.2f", $iy);
	my $spx = 0.;
	my $spy = 0.;
	my $str = "\\draw[mystycir] ($cix, $ciy) circle ($radius);\n";
	if ($withNumbers) {
		$str .= "\\node at ($cix + $spx, $ciy + $spy) {\\tiny $label };\n";
	}

	return $str;
}
sub createPlot
{
	my ($params) = @_;
	my $plot = "";

	# Params
	my $lx = $params->{nsites}/2;

	# Nodes
	my $scale = 0.2;
	my $y = 4*$scale;
	my $radius = 1*$scale;
	my $xoffset = 4*$scale;
	for (my $i = 0; $i < $lx; ++$i) {
		my $label1 = 2*$i;
		my $label2 = 2*$i + 1;
		my $ix = $i*$xoffset;
		$plot .= plotNode($ix,  0, $radius, $label1, 1);
		$plot .= plotNode($ix, -$y, $radius, $label2, 1);
	}

	return $plot;
}

sub verifyParams
{
	my ($h, $file, @keys) = @_;
	foreach my $k (@keys) {
		die "$0: Undefined key $k in $file\n" unless (exists($h->{$k}));
	}
}

sub readPreParams
{
	my ($templateInput) = @_;
	my %params;
	open(my $fh, "<", $templateInput) or die "$0: Cannot open $templateInput : $!\n";
	while (<$fh>) {
		chomp;
		if (/#([^ =]+) *= *(.*$)/) {
			$params{$1} = $2;
			next;
		}
	}

	close($fh);

	return %params;
}

sub impurityInit
{
	my ($templateInput, $preparams) = @_;

	verifyParams($preparams, $templateInput, qw/lx ly bathsPerSite tbath tphysical U/);

	my %hash = %$preparams;

	$hash{n} = $hash{lx}*$hash{ly}*(1+$hash{bathsPerSite});
	$hash{electronsUp} = int($hash{density}*$hash{n});
	$hash{tmatrix} = buildHoppingMatrix($preparams);
	$hash{hubbardU} = buildUvector($preparams);
	$hash{potentialV} = buildVvector($preparams);
	return \%hash;
}

sub buildVvector
{
	my ($h) = @_;
	my $n =	$h->{lx}*$h->{ly}*(1+$h->{bathsPerSite});
	my $str = "[";
	for (my $i = 0; $i < $n; ++$i) {
		$str .= "0";
		$str .= "," if ($i + 1 < $n);
	}

	$str .= "]";
	return $str;
}

sub buildUvector
{
	my ($h) = @_;
	my $block = 1+$h->{bathsPerSite};
	my $n = $h->{lx}*$h->{ly}*$block;
	my $str = "[";
	for (my $i = 0; $i < $n; ++$i) {
		my $isPhysical = (($i == 0) || ($i % $block == 0));
		my $val = ($isPhysical) ? $h->{U} : 0;
		$str .= $val;
		$str .= "," if ($i + 1 < $n);
	}

	$str .= "]";
	return $str;
}

sub buildHoppingMatrix
{
	my ($h) = @_;
	my @matrix = doConnections($h->{lx}, $h->{ly}, $h->{bathsPerSite},
		$h->{tphysical}, $h->{tbath});

	my $rows = $h->{lx} * $h->{ly} * (1 + $h->{bathsPerSite});
	symmetrizeMatrix(\@matrix, $rows);
	my $smatrix = matrixToString(\@matrix, $rows);
}

sub matrixToString
{
	my ($matrix, $rows) = @_;
	my $cols = $rows;
	my $str = "[";
	for (my $i = 0; $i < $rows; ++$i) {
		$str .= "[";
		for (my $j = 0; $j < $cols; ++$j) {
			my $val = $matrix->[$i][$j];
			$val = 0 unless defined($val);
			$str .= $val;
			$str .= "," if ($j + 1 < $cols);
		}

		$str .= "]";
		$str .= ",\n" if ($i + 1 < $rows);
	}

	$str .= "]";
	return $str;
}

sub symmetrizeMatrix
{
	my ($matrix, $rows) = @_;
	my $cols = $rows;
	for (my $i = 0; $i < $rows; ++$i) {
		for (my $j = 0; $j < $cols; ++$j) {
			my $val = $matrix->[$i][$j];
			if ($j < $i) {
				$val = $matrix->[$j][$i];
			}

			defined($val) or  $val = 0;
			$matrix->[$i][$j] = $val;
		}
	}
}

sub doConnections
{
	my ($lx, $ly, $bathsPerSite, $tphys, $tbath) = @_;

	my @matrix;
	my $block = 1 + $bathsPerSite;
	# Physical sites in the x direction (right)
	for (my $i = 0; $i < $lx; ++$i) {
		my @tmpVec;
		# Physical sites in the y direction (down)
		for (my $j = 0; $j < $ly; ++$j) {
			# Index of physical site
			my $index = $block * $j + $i * $ly * $block;

			# Connect physical to bath sites
			connectPhysicalBath(\@matrix, $index, $bathsPerSite, $tbath);

			# Connect physical to physical (y-direction down)
			# Open BC in y direction
			if ($j + 1 < $ly) {
				my $jindex = $index + $block;
				$matrix[$index][$jindex] = $tphys;
			}

			# Connect physical to physical (x-direction right)
			# Open BC in x direction
			if ($i + 1 < $lx) {
				my $jindex = $index + $ly * $block;
				$matrix[$index][$jindex] = $tphys;
			}
		}
	}

	return @matrix;
}

sub connectPhysicalBath
{
	my ($matrix, $index, $bathsPerSite, $tbath) = @_;
	for (my $b = 1; $b < $bathsPerSite +  1; ++$b) {
		$matrix->[$index][$index + $b] = $tbath;
	}
}

sub addVector
{
	my ($v1, $v2) = @_;
	my $n = scalar(@$v1);
	if ($n != scalar(@$v2)) {
		die "$0: addVector sizes different\n";
	}

	for (my $i = 0; $i < $n; ++$i) {
		$v1->[$i] += $v2->[$i];
	}
}

sub staggeredVectorEx
{
	my ($cmp, $nsites, $value) = @_;
	die "$0: nsites $nsites is not even\n" if ($nsites & 1);

	my $lx = $nsites/2;
	my $total = $nsites - 2;
	my $valX = ($cmp) ? 0 : $value;
	my $valY = ($cmp) ? $value : 0;

	my @v;
	$v[0] = $valY;
	my $count = 0;
	for (my $i = 1; $i < $total; ++$i) {
		$v[$i] = ($count & 1) ? $valY : $valX;
		++$count if (($i & 1) == 0);
	}

	return @v;
}

sub toString
{
	my ($v) = @_;
	my $n = scalar(@$v);
	return "[]" if ($n == 0);

	my $str = "[".$v->[0];
	for (my $i = 1; $i < $n; ++$i) {
		$str .= ", ".$v->[$i];
	}

	return $str."]";
}

sub addConstantToVector
{
	my ($v, $t) = @_;
	my $n = scalar(@$v);
	for (my $i = 0; $i < $n; ++$i) {
		$v->[$i] += $t;
	}
}

sub constantVector
{
	my ($total, $t) = @_;
	my @tvectorX = ($t);
	for (my $i = 1; $i < $total; ++$i) {
		push @tvectorX, $t;
	}

	return @tvectorX;
}

sub staggeredVector
{
	my ($cmp, $total, $t) = @_;
	my $torZero = ($cmp == 0) ? 0 : $t;
	my @v = ($torZero);
	for (my $i = 1; $i < $total; ++$i) {
		my $val = (($i & 1) == $cmp) ? 0 : $t;
		push @v, $val;
	}

	return @v;
}

sub createInput
{
	my ($file, $params, $templateInput) = @_;

	open(FOUT, ">", "$file") or die "$0: Cannot write to $file\n";

	open(FILE, "<", "$templateInput") or die "$0: Cannot open $templateInput: $!\n";

	my $plot = $params->{"plot"};

	while(<FILE>) {
		if (/^#/) {
			print FOUT;
			next;
		}

		if (/\$([a-zA-Z0-9]+)/) {
				my $name = $1;
				my $str = "\$"."params->{\"$name\"}";
				my $val = eval "$str";
				defined($val) or die "$0: Undefined substitution for $name\n";
				s/\$\Q$name/$val/g;
		}
		print FOUT;
	}

	close(FILE);
	close(FOUT);
	print STDERR "$0: File $file has been written\n";

	return $file;
}
