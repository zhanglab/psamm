#!/usr/bin/perl
use strict;
use warnings;

##################################################
#Created by Ying Zhang @ URI, Apr. 2014
#This script parse out ModelSEED reaction IDs for an annotated genomes 
#and write combined information to a reaction file
#
#Usage: ./make_rxn_table.pl <ANNOFILE> <SEEDRXNFILE> <SEEDCPDFILE> 
#-<ANNOFILE>: the input annotation file that were prepared based on genome annotation pipelines
#-<SEEDRXNFILE>: the reaction database downloaded from ModelSEED
#-<SEEDCPDFILE>: the compound database downloaded from ModelSEED
##################################################

if($#ARGV+1 < 3){
    print "ERROR: Usage: ./make_rxn_table.pl <ANNOFILE> <SEEDRXNFILE> <SEEDCPDFILE>\n";
    exit;
}

my $annofile = shift;
my $seedrxnfile = shift;
my $seedcpdfile = shift;
my $outrxn = "$annofile\_rxn";
my $outcpd = "$annofile\_cpd";
my $outgene = "$annofile\_gene";

###########Initiation###########################
my @tmplns = ();
#-initiate the parent node info of each taxon id
@tmplns = `cat $seedcpdfile`;
my %cpdinfo = ();
my %cpdmap = (); #map SEED compound names to SEED compound ids
foreach(@tmplns){
    chomp;
    s/\r$//;
    my ($seedcid,$names,$formula,$mass,$keegmaps,$keggcid) = split /\t/;
        for my$name(split /,<br>/,$names){
	    $name =~ s/&#39;//g;
	    $name =~ s/<92>//g;
            $cpdmap{$name} = $seedcid;
        }
        $cpdinfo{$seedcid} = "$formula\t$mass\t$keggcid\t$names";
    }

#-initiate taxon id mapped to scientific names
@tmplns = `cat $seedrxnfile`;
my %rxninfo = ();
foreach(@tmplns){
    chomp;
    s/\r$//;
    my @tmp = split /\t/;
    my ($seedrid,$name,$equation,$keggmaps,$ec,$keggrid) = ($tmp[0],$tmp[1],$tmp[2],$tmp[5],$tmp[6],$tmp[7]);
    $equation =~ s/&#39;//g;
    $equation =~ s/<92>//g;
    my $equation2 = &convert_equation($equation,%cpdmap);
    $rxninfo{$seedrid} = "$name\t$ec\t$equation\t$equation2\t$keggrid\t$keggmaps";
}

#-initiate the summary information
@tmplns = `cut -f11-14,17,19 $annofile`;
shift @tmplns;
my %genomerxn = (); #store the rxn info in a genome
my %geneinfo = (); #store the gene annotation
foreach(@tmplns){
    chomp;
    my ($anno,$subsys,$class1,$class2,$rxnlist,$gis) = split /\t/;
    next if $rxnlist eq '-';
    next if $gis eq '-';
    for my$seedrid(split /,/, $rxnlist){
        unless(defined $genomerxn{$seedrid}){ $genomerxn{$seedrid} = $gis }
 	else{ $genomerxn{$seedrid} .= ",$gis"; }
    }
    for my$gi(split /,/,$gis){
    	$geneinfo{$gi} = "$anno\t$subsys\t$class1\t$class2";
    }
}
############Initiation END###########################

my $titlerxn = "SEED_rid\tRXN_name\tEC\tEquation_cpdid\tEquation_cpdname\tKEGG_rid\tKEGG_maps\tGene_ids\n";
my $titlecpd = "SEED_cid\tFormula\tMass\tKEGG_cid\tCPD_name\n";
my $titlegene = "Gene_id\tAnnotation\tSubsystem\tSubsys_class1\tSubsys_class2\n";

open OUTRXN,">$outrxn"||die;
open OUTCPD,">$outcpd"||die;
open OUTGENE,">$outgene"||die;

print OUTRXN $titlerxn;
print OUTCPD $titlecpd;
print OUTGENE $titlegene;

my %cpds = (); #the compounds visited in the metabolic model
my %genes = (); #the genes visited in the metabolic model
for my $seedrid(keys %genomerxn){
    my $gis = $genomerxn{$seedrid};
    unless(defined $rxninfo{$seedrid}){ print "$seedrid\n"; exit; }
    print OUTRXN "$seedrid\t$rxninfo{$seedrid}\t$gis\n";
    for my $gi(split /,/,$gis){
	$genes{$gi} = 1;
    }
    #-get compound info
    my @tmp = split /\t/,$rxninfo{$seedrid};
    my $equation = $tmp[3];
    my @tmpcpds = &get_cpds($equation);
    for my$c(@tmpcpds){ $cpds{$c} = 1; }
}

for my $gi(sort keys %genes){
    print OUTGENE "$gi\t$geneinfo{$gi}\n"
}

for my $seedcid(sort keys %cpds){
    print OUTCPD "$seedcid\t$cpdinfo{$seedcid}\n";
}

close OUTRXN;
close OUTCPD;
close OUTGENE;

######## SUBS ##################
sub convert_equation{ #convert cpd names in equations into cpd ids
#-$e: the original equation provided in ModelSEED
#-%map: the mapping from cpd names to cpd ids
    my ($e,%map) = @_;
    my $new = ''; #the converted equation
    while($e =~ /\|([^\|]+)\|/){
	my ($before,$name,$after) = ($`,$1,$');
	my $suffix = '';
	#e-extracellular;p-pariplasmic;m-mitocondria;r-endoplasmic;v-vaculor;x-peroxisomal;g-Golgi;n-polymers
	if($name =~ /\[[epmrvxgn]\]$/){ $name = $`; $suffix = $&; }
	unless(defined $map{$name}){ print "ERROR: undefined compound name: $name\n"; exit; }
	$new .= "$before|$map{$name}$suffix|";
	$e = $after;
    }
    $new .= $e;
    return $new;
}

sub get_cpds{ #get cpd ids from an equation
    my ($e) = @_;
    my @c = ();
    while($e =~ /\|([^\|]+)\|/){ 
	my $cpdid = $1; 
	$e = $'; 
	if($cpdid =~ /\[[epmrvxgn]\]$/){ $cpdid = $`; }
	push @c, $cpdid;
    }
    return @c;
}

