#!/usr/bin/perl

require "../../libs/w2web.pl";
require "expert.pm";

&GetInput;
&GetSession;

my $chatbot   = new Chatbot::Eliza;
srand( time ^ ($$ + ($$ << 15)) );    # seed the random number generator

if ( $Comment ) {
  $prompt = $chatbot->transform( $Comment ) ;
} else {
  $prompt = $chatbot->transform('Hello');
}


$OUT = <<__STOP__;

<h3>your online "expert"</h3>

<p><b>Expert:</b> $prompt</p>

<form action="expert.pl" method="POST">
<p>
<textarea name="Comment" wrap=yes rows=3 cols=60></textarea>
<br>
<input type=submit value="pose question">
</p>
</form>

__STOP__
PrintPage("expert", $OUT);
