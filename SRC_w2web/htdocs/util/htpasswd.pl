#!/usr/bin/perl

require "../../libs/w2web.pl";

# password file to manage
$pwdfile="$w2webdir/conf/w2web.users";
$pidfile="$w2webdir/logs/w2web.pid";


&GetInput;

# Get all uids
open(PWD, '<'.$pwdfile) or (&CGIError("Can't open pwd file. Reason: $!")and exit);
flock(PWD, 1) or (&CGIError("Can't get file lock. Reason: $!") and exit);
while(<PWD>) {
  chomp;
  (!/[\W]/ || /\#/) and next;
  ($uid, $pwd) = split /\:/;
  $db{$uid} = $pwd;
}
close PWD;

# Delete uid?
for(keys %db) { (delete $db{$_} and $save = $_ and last)  if($FORM{$_}); }
if($save) {
  open(PWD, '>'.$pwdfile) or (print "Can't open pwd file for writing. Reason: $!" and exit);
  flock(PWD, 2) or (print "Can't get exclusive file lock. Reason: $!" and exit);
  for(sort keys %db) { print PWD $_.":".$db{$_}."\n"; }
  close PWD; 
  $msg = "User $save deleted. ";
}

# Add uid?
if($FORM{'add'} && $FORM{'uid'} && $FORM{'pwd'}) {
  @chars = ('a'..'z', 'A'..'Z', 0..9);
  ($msg .= $FORM{'uid'}." not saved; password is too long!" and last)
    if(length($FORM{'pwd'}) > 13);  # too long?
  ($msg .= $FORM{'uid'}." not saved; password contains invalid chars!" and last)
    if($FORM{'pwd'} =~ /[\W]/);  # invalid chars?
  $msg .= $FORM{'uid'}." saved.";
  $msg .= " (Warning: password isn't very secure!)"
    if($FORM{'pwd'} !~ /[a-zA-Z]/ || $FORM{'pwd'} !~ /[\d]/ ||
      length($FORM{'pwd'}) < 6);  # secure?
  $FORM{'pwd'} = crypt($FORM{'pwd'}, $chars[rand(@chars)].$chars[rand(@chars)]);
  open(PWD, '>>'.$pwdfile) or (print "Can't open pwd file for writing. Reason: $!" and exit);
  flock(PWD, 2) or (print "Can't get exclusive file lock. Reason: $!" and exit);
  print PWD $FORM{'uid'}.":".$FORM{'pwd'}."\n";
  close PWD;
  $db{$FORM{'uid'}} = $FORM{'pwd'};
}

# Output
for(sort keys %db) {
  $output .= qq|<tr><td ALIGN=RIGHT><input type="checkbox" name="$_" value="1"></td><td>$_</td></tr>\n|;
}
$output = qq|<tr><td colspan=2><small>[no Users]</small></td></tr>\n|  if(!$output);

$OUT = <<__STOP__;
<h1>User management</h1>
<div>$msg</div>
<form action="htpasswd.pl" method="post">
<table border>
<tr><th>Delete</th><th>User</th></tr>
$output
</table>
<p>
<input type="checkbox" name="add" value="1">Add User:<BR>
User: 
<tt><input type="text" name="uid" value="" size="9" maxlength=8></tt> 
Password: 
<tt><input type="password" name="pwd" value="" size="9" maxlength=13></tt>
</p>
<input type="submit" VALUE="Execute!">
__STOP__

$pid=qx(cat $pidfile);
$cnt= kill 1, $pid;
$OUT .= "$cnt processes signalled\n";
PrintPage("Password","$OUT");
exit;

