#!/usr/bin/perl
require "libs/libs.pl";

$OUT .=  <<__STOP__;
<font face="Verdana" size="+2" color="Navy">Welcome To <b>HConf</b></font>
<p>
If you are seeing this, that means that you installed correctly,
and have successfully logged in!

<h2>Menu:</h2>
<table cellspacing="2" cellpadding="2" border="0">
<tr>
    <td><a href="/hconf/password.pl?change">Change Your Password</a> (Coming Soon)</td>
    <td><b>If you installed using the default password, it is recommended that you change your password now! If you don't, it is very possible that others may gain access to your system.</b></td>
</tr>
<tr>
    <td><a href="/hconf/processes.pl">View and Edit Processes</a></td>
    <td>You may take a look at your system's processes.</td>
</tr>
<tr>
    <td><a href="/hconf/users.pl">View, Edit and Delete Users</a></td>
    <td>Just what it says!</td>
</tr>
<tr>
    <td></td>
    <td></td>
</tr>
</table>

__STOP__
PrintPage("Welcome to HConf!", $OUT); # PrintPage(Title, Content) does not return anything and has an exit()
