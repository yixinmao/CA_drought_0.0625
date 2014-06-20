#!/usr/bin/perl -w
no warnings 'uninitialized';
$sum_data=0;
($infile) = ($ARGV[0]);
@JUNK = `cat $infile`;
$n=0;
for ($i=0;$i<@JUNK;$i++) {
#  ($data) = split(" ",$JUNK[$i]);
  $sum_data+=$JUNK[$i];
  $n++;
 }
printf "%f\n", $sum_data/$n;
