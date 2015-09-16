
### Create seq coverage plots
ln -fs ../geneCovPlot.R .
Rscript geneCovPlot.R covdir test.txt output TRUE
rm covdir/test0*bin covdir/*.desc

