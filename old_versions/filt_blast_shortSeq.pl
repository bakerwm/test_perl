#!/usr/bin/perl -w

$usage = "perl  filt_blast_shortSeq.pl  input.blast  >out.blast\n";

$blast = shift or die $usage;

open (F, "<$blast") or die;
#open (OUT,">$blast\.filt") or die;

while(<F>){
    chomp;
    @sep = split /\t/;
    @{$ha{$sep[0]}->{$.}} = @sep; # $. is the line number.
}

foreach $seed (sort keys %ha){
    ($max, $mismatch, $gap) = (1, 0, 0);
    foreach $i (sort {$a<=>$b} keys %{$ha{$seed}}){
        @lines = @{$ha{$seed}->{$i}};
        ($q_map, $q_mismatch, $q_gap) = ($lines[3], $lines[4], $lines[5]);
        
        $direction = ($lines[6] - $lines[7])*($lines[8] - $lines[9]);

# skip lines match the following criteria:
# 1. more than 2 mismatch, 2. more than 0 gap, 3. no the largest map length, 4. reverse mapped.        
        next if($q_mismatch >= 2 || $q_gap > 0 || $q_map < $max || $direction < 0); 

        if($q_map == $max){
            next unless($q_mismatch <= $mismatch && $q_gap <= $gap);
            ($max, $mismatch, $gap) = ($q_map, $q_mismatch, $q_gap);
            $target = $i;
        }else{
            ($max, $mismatch, $gap) = ($q_map, $q_mismatch, $q_gap);
            $target = $i;
        }
    }
    $out_line = join"\t",(@{$ha{$seed}->{$target}});
    print $out_line,"\n";
}
