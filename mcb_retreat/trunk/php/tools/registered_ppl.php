<html>
<head><title>MCB Retreat</title></head>
<body>
<?php
$link = pg_connect("host=localhost dbname=yhdb user=yh password=123456")
		or die('Could not connect: ' . pg_last_error());

$result = pg_query( "select name, pi, nights, roommate1, roommate2, roommate3, roommate4 from retreat.register order by name,pi" )
or die("SELECT Error: ".pg_last_error());
$num_rows = pg_num_rows($result);
//
print "<table width=829  border=2 align=center cellpadding=1 cellspacing=0>\n";
print " <tr bgcolor=#42C1FF>";
print "    <td>Name</td>";
print     " <td>PI Name </td>";
print     " <td>Nights</td>";
print     " <td>Roomate1</td>";
print     " <td>Roomate2</td>";
print     " <td>Roomate3</td>";
print     " <td>Roomate4</td>";
print "  </tr>";
//print "<table width=200 border=1>\n";
while ($get_info = pg_fetch_row($result)){
print "<tr bgcolor=#B9EBF9>\n";
foreach ($get_info as $field)
print "\t<td>$field</td>\n";
print "</tr>\n";
}
print "</table>\n";
print "";

pg_close($link);
?>
</body>
</html>