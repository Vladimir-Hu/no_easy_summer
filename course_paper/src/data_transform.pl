#!/usr/bin/perl
$filename = $ARGV[0];
open(FILE,"$filename");
while($lineread = <FILE>){
    if($lineread =~ m/(\d*\.\d*)\ +(\d*\.?\d*)\ +(\d*\.?\d*)?\ +(\d*\.?\d*)?/){
        print("$1 $2 $3 $4\n");
    }
}