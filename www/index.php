
<!-- This is the project specific website template -->
<!-- It can be changed as liked or replaced by other content -->

<?php

$domain=ereg_replace('[^\.]*\.(.*)$','\1',$_SERVER['HTTP_HOST']);
$group_name=ereg_replace('([^\.]*)\..*$','\1',$_SERVER['HTTP_HOST']);
$themeroot='http://r-forge.r-project.org/themes/rforge/';

echo '<?xml version="1.0" encoding="UTF-8"?>';
?>
<!DOCTYPE html
	PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
	"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en   ">

  <head>
	<meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
	<title><?php echo $group_name; ?></title>
	<link href="<?php echo $themeroot; ?>styles/estilo1.css" rel="stylesheet" type="text/css" />
  </head>

<body>

<!-- R-Forge Logo
<table border="0" width="100%" cellspacing="0" cellpadding="0">
<tr><td>
<a href="/"><img src="<?php echo $themeroot; ?>/images/logo.png" border="0" alt="R-Forge Logo" /> </a> </td> </tr>
</table>
 -->

<!-- get project title  -->
<!-- own website starts here, the following may be changed as you like -->

<?php if ($handle=fopen('http://'.$domain.'/export/projtitl.php?group_name='.$group_name,'r')){
$contents = '';
while (!feof($handle)) {
	$contents .= fread($handle, 8192);
}
fclose($handle);
echo $contents; } ?>

<!-- end of project description -->

<h1>Getting started</h1>

<p>You can install the package and run a simulation using the
following commands in R:</p>

<pre>
install.packages("animation","ggplot2") # for dependencies on CRAN
install.packages("nicholsonppp", repos="http://R-Forge.R-project.org")
sim <- sim.drift.selection()
df <- sim2df(sim)
## Plot 6 loci evolving over time:
loci.over.time(interesting.loci(df))
## Make an animation that summarizes the simulation:
evolution.animation(df)
</pre>

<h2>Animations of allele frequency evolution simulation</h2>

<table>
  <tr>
    <td><img src="ani-thumb.png" /></td>
    <td><table>
	<tr>
	  <th>Date</th>
	  <th>Image size</th>
	  <th>Generations</th>
	</tr>
	<tr>
	  <td><a href="2009-08-10/index.htm">2009-08-10</a></td>
	  <td>1200x1000</td>
	  <td>200</td>
	</tr>
	<tr>
	  <td><a href="2009-08-19/index.htm">2009-08-19</a></td>
	  <td>1000x800</td>
	  <td>200</td>
	</tr>
      </table>
    </td>
  </tr>
</table>

<h2>Compiler notes</h2>

<p>The package has been compiled successfully using gfortran and
ifort. We prefer using ifort when possible, since it makes the program
run about 10x faster. To use ifort with R you have to recompile
R. Download the sources then do:</p>

<pre>
FC=ifort ./configure && make
</pre>

<p>Also make sure to add your ifort lib (something like
intel/Compiler/11.0/081/lib/ia32) to the list of libraries (normally
/etc/ld.so.conf), then run ldconfig. If you don't have root you may
want to just set LD_LIBRARY_PATH to the ifort lib.</p>

<p>The <strong>project summary page</strong> you can find <a href="http://<?php echo $domain; ?>/projects/<?php echo $group_name; ?>/"><strong>here</strong></a>. </p>

</body>
</html>
