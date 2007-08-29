<html>
<head><title>MCB Retreat</title></head>
<body>
<?php
$db="retreat";
$link = mysql_connect("localhost","root","");
if (! $link)
die("Couldn't connect to MySQL");
mysql_select_db($db , $link)
or die("Couldn't open $db: ".mysql_error());
$result = mysql_query( "SELECT * FROM register" )
or die("SELECT Error: ".mysql_error());
$num_rows = mysql_num_rows($result);
//
print "<table width=829  border=2 align=center cellpadding=1 cellspacing=0>\n";
print " <tr bgcolor=#42C1FF>";
print "   <td>ID</td>";
print "    <td>Name</td>";
print     " <td>Email </td>";
print     " <td>PI Name </td>";
print     " <td>T-Shirt Size </td>";
print     " <td>Roomate1</td>";
print     " <td>Roomate2</td>";
print     " <td>Roomate3</td>";
print     " <td>Roomate4</td>";
print "  </tr>";
//print "<table width=200 border=1>\n";
while ($get_info = mysql_fetch_row($result)){
print "<tr bgcolor=#B9EBF9>\n";
foreach ($get_info as $field)
print "\t<td>$field</td>\n";
print "</tr>\n";
}
print "</table>\n";
mysql_close($link);
?>
</body>
</html>
