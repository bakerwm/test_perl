#!/usr/bin/perl -w
use strict;
use warnings;
use POSIX qw(strftime);
use File::Basename qw(basename dirname);
use File::Spec::Functions qw(catdir catfile);

# Tables: http://www.textfixer.com/html/html-table-generator.php
# Input the list filefie and dir of images (named by id)

# Create multiple html files
# 1. Home.html
# 2. by_IDs.html
# 3. by_details.html
# 4. Html_lists.html
# 5. html_lists/A_to_B.html
#
# 6. css/style.css

my $usage = 'perl $0 infile.txt out.dir
infile.txt : a tab-separated file with at least 6-columns
out.dir    : store the html files of the images
';

my $infile = shift or die "Input a tab file:\n $usage";
my $outdir = shift or die "A directory to store the html files";
$outdir =~ s/(^\s+|\s+$)//g;

## check existence of image/
die "Not found [image/] in folder in $outdir ..." unless -d "$outdir/images";

## parse the current Date
my $cur_date = strftime("%Y-%m-%d %H:%M:%S\n", localtime(time));
#my $cur_date = '2015-03-19';

#my $css = catdir($outdir, "css");
my $nav_bars = catdir($outdir, "nav_bars");
my $html_split = catdir($outdir, "html_split");

## Checking output directory
#mkdir $css unless -d $css;
mkdir $nav_bars unless -d $nav_bars;
mkdir $html_split unless -d $html_split;

# parsing infile
my %seq_ids = ();
open F, "< $infile" or die "$!";
while(<F>){
    chomp;
    my $id = (split /\t/, $_)[0];
    $seq_ids{$id} = $_;
}
close F;

#####################################
# The elements of html files
my $html_header = '<!DOCTYPE html>
<html>
  <head>
    <title>Reads Coverage Maps</title>
    <meta name="R drop plots" content="Reads coverage.">
    <meta name="R coverage maps" content="Reads coverage.">
    
    <style type="text/css">
      .tftable {font-size:15px;color:#333333;width:90%;border-width: 1px;border-color: #9dcc7a;border-collapse: collapse;}
      .tftable th {font-size:14px;background-color:#abd28e;border-width: 1px;padding: 10px;border-style: solid;border-color: #9dcc7a;text-align:left;}
      .tftable tr {background-color:#ffffff;}
      .tftable td {font-size:14px;border-width: 1px;padding: 10px;border-style: solid;border-color: #9dcc7a;}
      .tftable tr:hover {background-color:#ffff99;}
    </style>
    
    <style>
      nav {
        line-height: 20px;
        background-color: #E6E6E6;
        width: 90%;
        float: left;
        padding: 1px;
      }
      
      nav ul {
        list-style-type: none;
        margin: 0;
        padding: 0;
        overflow: hidden;
      }
      
      nav ul li {
        float: left;
      }
      
      nav > ul > li > a:link, nav > ul > li > a:visited {
        display: block;
        width: 100px;
        font-weight: bold;
        color: #FFFFFF;
        background-color: #6E6E6E;
        text-align: center;
        padding: 2px;
        text-decoration: none;
      }
      
      nav > ul > li > a:hover, a:active {
        background-color: #98bf21; # green
      }
            
      header {
        background-color: #E6E6E6;
        color: #424242;
        text-align: left;
        padding: 1px;
        width: 90%;
      }
      
      body {
        width: 80%;
        margin-left: 200px;
        margin-right: 200px;
        font-family: "Helvetica", sans-serif;
      }
      
      section {
        width: 90%;
        float: left;
        padding: 1px; 
      }
      
      footer {
        height: 40px;
        width: 80%;
        #background-color:black;
        color:#6E6E6E;
        #clear:both;
        text-align:center;     
      }
    </style>
    
    </head>
<body>
';

my $nav_bar = '
<nav>
  <ul>
    <li><a href="index.html">Home</a></li>
    <li><a href="nav_bars/1.by_details.html">Seq List</a></li>
    <li><a href="nav_bars/2.by_IDs.html">ID List</a></li>
    <li><a href="nav_bars/3.html_split.html">Tables</a></li>
  </ul>
</nav>
';

my $nav_bar_sub = '
<nav>
  <ul>
    <li><a href="../index.html">Home</a></li>
    <li><a href="../nav_bars/1.by_details.html">Seq List</a></li>
    <li><a href="../nav_bars/2.by_IDs.html">ID List</a></li>
    <li><a href="../nav_bars/3.html_split.html">Tables</a></li>
  </ul>
</nav>
';

my $html_section = 'Content';
my $html_footer = '
<footer>
  <p><small>Generated at: ' . $cur_date . ' by @Ming, IBP</small></p>
</footer>

</body>
</html>
';

########################################
my $image_cont = '';
my @f1_table_new_tabs;

my %by_id_sub = ();
my %tab_12 = ();
my $f1_count = 1;
foreach my $id (sort keys %seq_ids){
    ## 1.by_ids.html
    my $image_dir = catdir($outdir, "images");
    my $image = catfile($image_dir, "$id\.png");
    warn "Image of $id - $image not found... \n" unless -e $image;
    my $by_id_group = sprintf "%d", (($f1_count-1)/50 + 1);
    push @{$by_id_sub{$by_id_group}}, $id;
       
    ## 2.by_details.html
    my @tabs = split /\t/, $seq_ids{$id}, 7;
    my @tabs_links = id_to_link($tabs[0]);
    $tabs[0] = $tabs_links[0];
    my $new_tab = join "\t", ($f1_count, @tabs);
    push @f1_table_new_tabs, $new_tab;
        
    ## 3.html_split.html
    my $tab_group = sprintf "%d", (($f1_count-1)/12 + 1);
    push @{$tab_12{$tab_group}}, $id;
    
    $f1_count ++;
}

####################################
# Generate 0. index.html #
####################################
my $home_header = '
<header>
<h1>Draw reads coverage maps of each sequence</h1>
</header>
';

my $home_section = '
<section>
<h2>Examples:</h2>
</section>

<figure>
  <img src="http://7xi7ro.com1.z0.glb.clouddn.com/Example.png" alt="PNG" width="500" height="375">
  <img src="http://7xi7ro.com1.z0.glb.clouddn.com/Example-freey.png" alt="PNG" width="500" height="375">
  <figcaption>Fig-1. Reads coverage of map of Rvnr01. <i>(left: fixed Y-axis, right: free Y-axis scal)</i></figcaption>
</figure>

<section>
<pre>
<p title="About Coverage map">
This file will present the coverage maps of a sequence.

x-axis: the coordination of the genome,
y-axis: the read counts in each library

Library:
1. 18-40 nt: the library was constructed using RNAs from 18 nt to 40 nt.
2. 40-80 nt: 40 nt to 80 nt RNAs included in this library.
3. 80-140 nt: 80 nt to 140 nt RNAs included in this library.
4. >140 nt: longer than 140 nt RNAs were included in this library.
</p>
</pre>

<pre><code>
Notes:
1. All of the files and images were created by <b>R</b> and <b>PERL</b>.
2. The structure of the folder are:

dir/ 
 |-index.html  # home page
 |
 |- html_split/ 
 |  |-seed0001-seed0012.html ... #pages 
 | 
 |-images/ 
 |  |-seed0001.png... # images 
 |  
 |-nav_bars/
 |  |-1.by_details.html # display details of seq 
 |  |-2.by_IDs.html     # search images by id 
 |  |-3.html_split.html # organize the split_html files  

  
</code></pre>

</section>
';

my $home_html = join "\n", ($html_header, $home_header, $nav_bar, $home_section, $html_footer);

my $home_html_file = catfile($outdir, "index.html");
open OUT, "> $home_html_file" or die "$!";
print OUT $home_html, "\n";
close OUT;

####################################
# Generate 1. by_details.html
####################################
my $det_header = '<header> <h1>List of genome info of sequences</h1> </header>
';

my $table_header = tabs_to_table('header', 'NUM', 'ID', 'Chr', 'Length /nt', 'Start', 'End', 'Strand', 'Note');
my @table_lines;
foreach my $t (@f1_table_new_tabs){
    my @tmp_tabs = split /\t/, $t, 8;
    my $tmp_line = tabs_to_table('line', @tmp_tabs);
    push @table_lines, $tmp_line;
}
my $table_line = join "\n", (@table_lines);

my $f1_table_section = $table_header . "\n" . $table_line;
my $table_block = join "\n", ('<table class="tftable" border="1">', 
                              $f1_table_section, '</table>');
## image to ../image
$table_block =~ s/\=\"images/\=\"\.\.\/images/g;

my $det_section = join "\n", ('<section>', $table_block, '</section>');
my $det_html = join "\t", ($html_header, $det_header, $nav_bar_sub, 
                           $det_section, $html_footer);

my $det_html_file = catfile("$outdir/nav_bars", "1.by_details.html");
open OUT, "> $det_html_file" or die "$!";
print OUT $det_html, "\n";
close OUT;

####################################
# Generate 2. by_IDs.html #
####################################
my $id_header = '
<header>
<h1>List of IDs</h1>
</header>
';

my @id_sections = '';
my $id_count = 1;
foreach my $id (sort {$a<=>$b} keys %by_id_sub){
    my @ids = sort @{$by_id_sub{$id}};
    my $id_start = $ids[0];
    my $id_end = $ids[-1];
    my $id_group_header = join "", ('<h3>', $id_count, '. ', $id_start,
                                   ' to ', $id_end, '</h3>');
    my @id_links = id_to_link(@ids);
    my $id_group_section = join "\n", @id_links;
    
    push @id_sections, ($id_group_header, $id_group_section, "\n");
    $id_count ++;
}
my $id_section = join "\n", ('<section>', @id_sections, '</section>');

## image to ../image
$id_section =~ s/\=\"images/\=\"\.\.\/images/g;

#$id_section = $id_header . "\n" . $id_section;
my $id_html = join "\n", ($html_header, $id_header, $nav_bar_sub, $id_section,
                          $html_footer);

my $id_html_file = catfile("$outdir/nav_bars", "2.by_IDs.html");
open OUT, "> $id_html_file" or die "$!";
print OUT $id_html, "\n";
close OUT;

####################################
# Generate 3. html_split.html / Table
####################################
my $tab_header = '
<header>
<h1>Table of Sequences</h1>
</header>
';

my $tab_table_header = '<tr><th>Num</th><th>File Name</th><th>Content</th></tr>';
my @tab_table_lines;
my $tab_count = 1;

foreach my $tab (sort {$a<=>$b} keys %tab_12){
    my @tab_ids = sort @{$tab_12{$tab}};
    my $tab_ids_start = $tab_ids[0];
    my $tab_ids_end = $tab_ids[-1];
    my $tab_html_file = catfile("../html_split", "$tab_ids_start\-$tab_ids_end\.html");
    my $tab_html_file_name = basename($tab_html_file);
    my $tab_id = join ", ", @tab_ids;
    my $tab_html_link = '<a href="'. $tab_html_file . '">' . $tab_html_file_name . '</a>';
    
    my $tab_table_line = tabs_to_table('line', $tab_count, $tab_html_link, $tab_id);
    push @tab_table_lines, $tab_table_line;
    $tab_count ++;
    
    #######################################
    ## Generate 4. Create separate html files
    #######################################
    my $tab_sub_table_section = ids_to_table_block(@tab_ids);
    my @tab_sub_fig_sections = ids_to_image_block(@tab_ids);
    my $tab_sub_fig_section = join "\n", @tab_sub_fig_sections;
    my $tab_sub_header = '
    <header>
    <h1> Table of ' . $tab_html_file_name . '
    </header>';
    
    my $tab_sub_fig_header = '<h3> Figures. Coverage maps of ' . $tab_html_file_name . '</h3>';
    my $tab_sub_section = join "\n", ('<section>', $tab_sub_table_section, $tab_sub_fig_header, 
                                      $tab_sub_fig_section, '</section>');
    # image to ../image
    $tab_sub_section =~  s/\=\"images/\=\"\.\.\/images/g;
    
    my $tab_sub_html = join "\n", ($html_header, $tab_sub_header, $nav_bar_sub,
                                   $tab_sub_section, $html_footer);
    
    my $tab_sub_html_file = catfile("$outdir\/html_split", $tab_html_file_name);
    open OUT, "> $tab_sub_html_file" or die "$!";
    print OUT $tab_sub_html, "\n";
    close OUT;
}
my $tab_table_section = join "\n", ('<section>', '<table class="tftable" border="1">', 
                                 $tab_table_header,
                                 @tab_table_lines,
                                 '</table>', </section>);
my $tab_html = join "\n", ($html_header, $tab_header, $nav_bar_sub, 
                           $tab_table_section, $html_footer);

my $tab_html_file = catfile("$outdir/nav_bars", "3.html_split.html");
open OUT, "> $tab_html_file" or die "$!";
print OUT $tab_html, "\n";
close OUT;


sub ids_to_image_block{
    my @ids = @_;
    # *.png and *-freey.png
    my @img_blocks;
    foreach my $i (@ids){
        warn ("PNG of $i not found") unless -e "$outdir/images/$i\.png";
        warn ("freey.PNG of $i not fount") unless -e "$outdir/images/$i\-freey\.png";
        my $cap = join "", ('<figcaption>', $i, '</figcaption>');
        my $img1 = join "", ('<img src="images/', $i, '.png" alt="', $i, 
                             '" style="width:500px;height:375px;border:0">');
        my $img2 = join "", ('<img src="images/', $i, '-freey.png" alt="', $i, 
                             '-freey" style="width:500px;height:375px;border:0">');
        my $i_block = join "\n", ('<figure>', $cap, $img1, $img2, '</figure>');
        push @img_blocks, ($i_block, "\n");
    }
    return @img_blocks;
}

sub ids_to_table_block{
    my @ids = @_;
    my $header = '<tr><th>Num</th><th>ID</td><th>Chr</td><th>Length /nt</td><th>Start</td><th>End</td><th>Strand</td><th>Note</td></tr>';
    my @lines;
    my $count = 1;
    foreach my $i (@ids){
        warn "$i not found..." unless(exists $seq_ids{$i});
        my @tabs = split /\t/, $seq_ids{$i}, 7;
        my $i_line = tabs_to_table('line', $count, @tabs);
        push @lines, $i_line;
        $count ++;
    }
    my $table_block = join "\n", ('<table class="tftable" border="1">', $header, @lines, '</table>');
    return $table_block;
}

sub id_to_link{
    my @in_ids = @_;
    #<a href="pics/t1.png">ID1</a>
    my @in_ids_links;
    foreach my $i (@in_ids){
        my $id_to_image = catfile("$outdir\/images", "$i\.png");
        warn("Image of $i not found...\n") unless -e $id_to_image;
        my $id_link = '<a href="'. "images\/$i\.png" . '">' . $i . '</a>';
        push @in_ids_links, $id_link;
    }
    return @in_ids_links;
}

sub tabs_to_table{
    my $table_type = shift(@_); # header or line
    die"Input table type should be [header] or [line].\n" unless($table_type eq 'header' || $table_type eq 'line');
    my @tabs = @_;
    my @tabs_fmt = ();
    foreach my $tab (@tabs){
        my $tab_fmt = '<td>'. $tab . '</td>';
        push @tabs_fmt, $tab_fmt;
    }
    my $line = join "", ('<tr>', @tabs_fmt, '</tr>');
    $line =~ s/<td>/<th>/g if($table_type eq 'header');
    return $line;  
}
