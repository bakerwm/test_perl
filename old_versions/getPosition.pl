#!/usr/bin/perl -w
use Data::Dumper; ##

sub help{
    print STDERR <<EOF;
Usage: perl $0  infile gff  strain_id  > outfile

Note:
    -1      : The input file in the following format:
              <ID> <Strand> <Begin:Cov> <Max:Cov> <End:Cov>
    -2      : The GFF file. GFF3 recommend.
    -3      : The Strain name, help to recognize the GFF file.

Example: perl $0  H37Rv_45SE.temp  NC_000962.gff H37Rv  >H37Rv_45SE_lncRNA.txt
EOF
}

$infile = shift or die &help;   # input seed file.
$gff    = shift or die $usage;  # input Gff file.
$strain = shift or die $usage;  # strain id.

open IN,$gff or die;
$i=1;
while(<IN>){
    next unless(/\tgene\t/);
    @ps=split(/\t/,$_,9);
    
   if($strain =~ /B42|B36|BT/i){                
        $ps[8]=~/Name=(\w+)\;/;      $id=$1;    # This is GFF2 format.
    }else{      
        $ps[8]=~/locus_tag=(\w+)/;   $id=$1;    # This is GFF3 format.
    }

    @{$gff{$i}}=($id, $ps[3], $ps[4], $ps[6]);      # Keys: the order,  Value: id, begin, end, strand;
    $i++;
}
close IN;

#print Dumper(%gff);
#exit;

open IN,$infile or die;
while(<IN>){
    next if(/begin\:cov/i);
    chomp;
    ($s_id, $s_str, $p1, $p2, $p3)=split/\t/;
    ($s_b, $t)=split/\:/,$p1;
    ($t, $s_exp)=split/\:/,$p2;
    ($s_e, $t)=split/\:/,$p3;
    $s_len=$s_e-$s_b+1;

    foreach $j (sort{$a<=>$b} keys %gff){
        ($Ga_id, $Ga_b, $Ga_e, $Ga_str)=@{$gff{$j}};
        $k=$j+1;    last unless(exists $gff{$k});
        ($Gb_id, $Gb_b, $Gb_e, $Gb_str)=@{$gff{$k}};
#        $l=$k+1;
#       ($Gc_id, $Gc_b, $Gc_e, $Gc_str)=@{$gff{$l}};
        
        $L1=$s_b-$Ga_b; $L2=$s_b-$Ga_e; $L4=$s_e-$Ga_e;
        $R1=$s_b-$Gb_b; $R3=$s_e-$Gb_b; $R4=$s_e-$Gb_e;

        if($L1>=0 && $L2<=0){
           $g1=$Ga_id;  # the following script is to determin $gap_2 and $dir.
            if($L4<=0){
                $gap_1=$L1;
                if($s_str eq $Ga_str){$des='IM';}else{$des='AS';}
                $g2=$Ga_id; $gap_2= -$L4;  $dir="\/$Ga_str\/$s_str\/$Ga_str\/";
            }elsif($L4>0 && $R3<0){
                $gap_1=$L2;
                if($s_str eq $Ga_str){$des='PM';}else{$des='AS';}
                $g2=$Gb_id; $gap_2= -$R3;  $dir="\/$Ga_str\/$s_str\/$Gb_str\/";
            }elsif($R3>=0 && $R4 <=0){
                $gap_1=$L2;
                $g2=$Gb_id; $gap_2= -$R3; $dir="\/$Ga_str\/$s_str\/$Gb_str\/";
                if($s_str eq $Ga_str ){
                    if($s_str eq $Gb_str){$des='PM2'; }
                    else{$des='PM'; }
                }else{
                    if($s_str eq $Gb_str){$des='PM';  }
                    else{$des='AS2'; }
                }
            }elsif($R4>=0){
                $gap_1=$L2;
                $g2=$Gb_id; $gap_2= -$R4; $dir="\/$Ga_str\/$s_str\/$Gb_str\/";
                if($s_str eq $Ga_str){
                    if($s_str eq $Gb_str){$des='PM3';}
                    else{$des='PM'; }                        
                }else{
                    if($s_str eq $Gb_str){$des='PM';}
                    else{$des='AS3'; }
                }
            }else{
                $gap_1=$L2;
                $des='EX'; $g2=$Gb_id; $gap_2= -$R4; $dir="\/$Ga_str\/$s_str\/$Gb_str\/";
            }
       
       }elsif($L2>0 && $R1<0){
           $g1=$Ga_id; $gap_1=$L2; $g2=$Gb_id; $gap_2= -$R3; $dir="\/$Ga_str\/$s_str\/$Gb_str\/";
           if($R3<0){
                $des='IGR'; 
           }elsif($R3>=0 && $R4<=0){
                if($s_str eq $Gb_str){$des='PM'; }
                else{$des='AS';}
           }else{
                $gap_1=$R1; $gap_2= -$R4;
                if($s_str eq $Gb_str){$des='CM'; }
                else{$des='AS2';}               
           }
        
       }else{
       
       }
    }
    $out=join"\t",($s_id, $s_exp, $s_len, $s_b, $s_e, $s_str, $g1, $gap_1, $g2, $gap_2, $dir, $des);
    print $out,"\n";
}
close IN;
