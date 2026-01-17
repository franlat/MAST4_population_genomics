#!/usr/bin/perl -w
use strict;
use Getopt::Long;

our ( $IN , $OUT , $L ,$ID, $STAR, $I, $HELP ) = ( "" , "" , 0.8 , 0.95 , 0 , "identity.txt" ,0 );
my $result = &GetOptions ( "i=s"        => \$IN ,
                                "o=s"   => \$OUT ,
                                "l=s"   => \$L ,
                                "d=s"   => \$ID ,
                                                                "star"  => \$STAR,
                                "p=s"     => \$I,
                                "h"     => \$HELP );

my $tmp_file = rand (10000000000000000000000000000000);

usage () if ( ! -f $IN || $OUT eq "" || $HELP );


my $init_count = 0;
my $final_count = 0;
my $TMP_DIR = "tmp_sam"; 
`mkdir $TMP_DIR` if ( ! -d $TMP_DIR );
`samtools view -H $IN > $TMP_DIR/$tmp_file.sam`;

open ( OUT , ">>$TMP_DIR/$tmp_file.sam" );
open ( IN , "samtools view $IN |" ) or die "$!\n";
open (OU , ">$I" ) or die "$!\n";
while ( <IN> ){
        ++$init_count;
        my @data=split(/\t/,$_);
                #$align -> taille de l'alignement
        my $align=0;
                #$match -> nb match et mismatch
                my $match=0;
                my $d=$data[5];
                # pour la taille : compte mismatch/match, insertion et deletion
        while($d=~s/(\d*)[MID]//){
                        $align+=$1;
        }
        while($data[5]=~s/(\d*)M//){
                        $match+=$1;
                }
                #$edit -> nb mismatch
        my $edit=0;
        for my $i (6..$#data){
                                # deux cas pour les mismatchs : sortie de star nM et bwa NM
                if($data[$i]=~/[nN]M\:i\:/){
                        $edit=$data[$i];
                                $edit=~s/[nN]M\:i\://;
                                last;
                                }
        }

                if($STAR) {
                        # Ne garde que les matchs dont le pourcentage d'identite est superieur <C3><A0> $ID
                        my $percent = ( $match - $edit ) / $align;
                        if($percent >= $ID){
                print OU $percent,"\n";
                                ++$final_count;
                                print OUT $_;
                        }
                }else{
                if( ( $match/length($data[9]) >= $L ) && ( ( $align - $edit ) / $align >= $ID ) ){
                print OU ( $align - $edit ) / $align,"\n";
                ++$final_count;
                print OUT $_;
                }
                }
}
close (IN);
close (OUT);
close (OU);
`samtools view -bS -o $OUT $TMP_DIR/$tmp_file.sam; rm $TMP_DIR/$tmp_file.sam`;
print "MinAlignLength: $L\n";
print "MinIdentity: $ID\n";
print "Initial count of reads:\t$init_count\n";
print "Selected reads:\t$final_count\n";


###############################################################################
sub usage {
        my $usage = "
        bamHomoFilter.pl : select reads that mapped over alignment features

        Usage:
        bamHomoFilter.pl -i in.bam -o out.bam [-l MinAlignLength -d MinIdentity]

        Options:

        -l      The minimum alignment size as a ratio of the read length (Default 0.8)
        -d      The minimum identity of the alignment (Default 0.95)
                -star   Only for a star mapping to select with minimum identity

        Warning: samtools needs to be in your \$PATH

";
        warn ( $usage );
        exit (1);
}
