








#### GeneID2UniProt
#sub GeneID2UniProt{
#    my ($from, $to, $input_ids) = @_;
#    my $base = 'http://www.uniprot.org';
#    my $tool = 'mapping';
#    my $params = {
#        from => $from,
#        to => $to,
#        format => 'tab',
#        query => $input_ids
#    };
#    #Please set your email address here to help us debug in case of problems.
#    my $contact = 'wangmcas@gmail.com';
#    my $agent = LWP::UserAgent->new(agent => "libwww-perl $contact");
#    push @{$agent->requests_redirectable}, 'POST';
#
#    my $response = $agent->post("$base/$tool/", $params);
#    while (my $wait = $response->header('Retry-After')) {
#        print STDERR "Waiting ($wait)...\n";
#        sleep $wait;
#        $response = $agent->get($response->base);
#    }
#    my %result;
#    if($response->is_success){
#        my @lines = split (/\n/, $response->content);
#        foreach my $line (@lines){
#            next unless( my ($f, $t) = split /\s+/, $line);
#            $result{$f} = $t if $f ne 'From';
#        }
#    }else{
#        die 'Failed, got ' . $response->status_line . ' for ' .
#            $response->request->uri . "\n";
#    }
#    return %result;
#}
#


use strict;
use warnings;
use LWP::UserAgent;

my $base = 'http://www.uniprot.org';
my $tool = 'mapping';

my $params = {
    from   => 'P_ENTREZGENEID',
    to     => 'ACC',
    format => 'tab',
#    query  => '4441282 4442831 4442568 4441019 4441304 4442231 4441242 4442154 4441662 4442417 4442077'
    query => '14515852 14515879 14515880 14515883 14515893 14515898'
};

my $contact = 'wangmcas@gmail.com'; # Please set your email address here to help us debug in case of problems.
my $agent = LWP::UserAgent->new(agent => "libwww-perl $contact");
push @{$agent->requests_redirectable}, 'POST';

my $response = $agent->post("$base/$tool/", $params);

while (my $wait = $response->header('Retry-After')) {
    print STDERR "Waiting ($wait)...\n";
    sleep $wait;
    $response = $agent->get($response->base);
}

$response->is_success ?
    print $response->content :
    die 'Failed, got ' . $response->status_line .
    ' for ' . $response->request->uri . "\n";





