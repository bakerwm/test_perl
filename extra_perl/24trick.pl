#!/usr/local/bin/perl  -w
use strict;
use warnings;

my $usage1="Input 4 intger numbers between 1 and 13:";
print "Input 4 intger numbers between 1 and 13: ";
while(<STDIN>){
    my @input=split /\s/;
    my @data=&Check(@input);
    my @list=&Order(@data);
    foreach(@list){
        my @lines=split/\:/,$_;
        my @Func=&Sep(@lines);
        foreach my $lib(@Func){
            my $res=eval"$lib";
            if($res==24){
                print "$lib = $res\n";
                exit;
            }
        }
    }
}

sub Check{
    my @num=@_;
    foreach(@num){
        if($_>13 || $_<1 || @_!=4){
            die  $usage1;
        }
        $_=int($_);
    }
    return @num;
}

sub  Order{
    my @val=@_;
    my %ha=();
    my @key=(1,2,3,4);
    @ha{@key}=@val;
    my @arr=();
    foreach my $i(sort {$a<=>$b} keys %ha){
       my %t1=%ha;
       delete($t1{$i});
       foreach my $j(sort {$a<=>$b} keys %t1){
           my %t2=%t1;
           delete($t2{$j});
           foreach my $k(sort {$a<=>$b} keys %t2){
               my %t3=%t2;
               delete($t3{$k});
               foreach my $m(sort {$a<=$b} keys %t3){
                   my $out=join"\:",($ha{$i},$t1{$j},$t2{$k},$t3{$m});
                   push @arr,$out;
               }
           }
       }
    }
    return @arr;
}

sub Sep{
    my @tab=('+','-','*','/');
    my @arrs=@_;
    my @brrs=my @crrs =my @drrs =@arrs;
    my $d1=shift @arrs;
    my @hit1=($d1);
    while(1){
        last if(@arrs < 1);
        my $d2=shift @arrs;
        my @hit1New=();
        foreach my $num (@hit1){
            foreach my $t(@tab){
                my $new="$num $t $d2";
                push @hit1New,$new;
            }
        }
        @hit1=@hit1New;
    }
    
    @brrs=("\($brrs[0]", "$brrs[1]\)", "\($brrs[2]", "$brrs[3]\)");
    my $e1=shift @brrs;
    my @hit2=($e1);
    while(1){
        last if(@brrs < 1);
        my $e2=shift @brrs;
        my @hit2New=();
        foreach my $num (@hit2){
            foreach my $t(@tab){
                my $new="$num $t $e2";
                push @hit2New,$new;
            }
        }
        @hit2=@hit2New;
    }

    @crrs=("\($crrs[0]", "$crrs[1]\)", $crrs[2], $crrs[3]);
    my $f1=shift @crrs;
    my @hit3=($f1);
    while(1){
        last if(@crrs < 1);
        my $f2=shift @crrs;
        my @hit3New=();
        foreach my $num (@hit3){
            foreach my $t (@tab){
                my $new="$num $t $f2";
                push @hit3New,$new;
            }
        }
        @hit3=@hit3New;
    }
    
    @drrs=("\($drrs[0]", $drrs[1], "$drrs[2]\)", $drrs[3]);
    my $g1=shift @drrs;
    my @hit4=($g1);
    while(1){
        last if(@drrs < 1);
        my $g2=shift @drrs;
        my @hit4New=();
        foreach my $num (@hit4){
            foreach my $t(@tab){
                my $new="$num $t $g2";
                push @hit4New,$new;
            }
        }
        @hit4=@hit4New;
    }
    return (@hit1, @hit2, @hit3, @hit4);
}
