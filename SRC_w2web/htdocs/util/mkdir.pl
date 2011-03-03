#!/usr/bin/perl

require "../../libs/w2web.pl";
&GetInput;
&GetSession;

my $absdir = "$dir/$newdir";
$umps = qx( mkdir $absdir );

$DIR = $absdir;

redirectURL("/util/dir.pl?SID=$SID&dir=$absdir&cd=$cd&f=$f");


