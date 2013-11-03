#!/bin/bash
#Generate figures for the manuscript using the program's output in ../fig.

outdir="../plosone"
indir="../fig"

if [[ ! -d  $outdir ]];
then
   mkdir $outdir
fi

# Convert pdfs to tifs ImageMagick

in=(hists_nolegend.pdf meandisp_scatter_nolegend.pdf mse.pdf phi_vs_phi_II_nolegend.pdf phi_vs_phi_V_nolegend.pdf auc1.pdf auc2.pdf auc3.pdf auc4.pdf auc5.pdf auc6.pdf)

for (( i=0; i<${#in[@]} ; i++ ))
do
  convert -strip -units PixelsPerInch -density 300 -resample 300 -alpha off -colorspace RGB -depth 8 -trim -bordercolor white -border 1% -resize '2049x2758>' -resize '980x980<' +repage -compress lzw $indir/${in[$i]} $outdir/Figure$[$i+1].tiff
done