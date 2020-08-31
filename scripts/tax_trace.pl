#!/usr/bin/env perl  

use strict;

die "Usage:perl $0   <node.dmp>   <name.dmp>  <taxids>  <export>" if( @ARGV != 4 );


my ($node, $name, $taxids, $export) = @ARGV;


open (NODE,         $node) ||die "$!"; 
my %parents = ();
my %rank    = (); 
while(<NODE>){
   my @its = split /\s+\|\s+/ , $_;
   $parents{ $its[0] } = $its[1];
   $rank{ $its[0] }    = $its[2];
}
close NODE;


open (NAME ,       $name) || die "$!"; 
my %names    = (); 
while(<NAME>){
    
    my @its = split /\s+\|\s+/ , $_;    
    if(/scientific name/){
        $names{  $its[0] }  = $its[1];
    }

}
close NAME;


open( TAXID,  $taxids     )   || die "$!";
open( EXPORT, ">$export" )   || die "$!";

while (<TAXID>) {

  $_ =~ s/\s+$//;
  my @its = split /\t/, $_;
  my ($seqId, $taxId) = ($its[0], $its[1]);
  
  if(! exists $names{ $taxId }){
      print "$seqId, no taxonomy found\n";
	  print EXPORT qq{no taxonomy found\n};
      next;
  }

  my ($node_path, $name_path)  = taxon_trace( $taxId );
  #print EXPORT qq{$name_path\n};
  print EXPORT qq{$seqId\t$names{ $taxId }\t$name_path\t$node_path\n};
}
close TAXID;
close EXPORT;


exit;

sub taxon_trace{

    my $node      = shift; 
    my @rank      = ();
    my @name_path = ();
    
    while(1){               
        
        push @rank, $rank{ $node };
        push @name_path, $names{ $node };

        if($node == 1){
           last;
        }
        
        if(exists $parents{ $node }){
            $node = $parents{ $node };
        }else{
           warn "$node \t something may be wrong!\n";
           exit;  
        }
    }
    return  ( join("|", reverse @rank), join("|", reverse @name_path) );
}