import os
from sys import argv
from Bio import SeqIO
#from Bio.SeqUtils import GC
from bokeh.plotting import figure, output_file, show, ColumnDataSource
from bokeh.models import BoxAnnotation, LabelSet, Label
import pandas as pd


#script, map_file = argv
script, ref_file, readfile1 = argv

os.system('minimap2 -ax map-ont %s %s > %s_consensus_mapping.sam'%(ref_file,readfile1,ref_file))


#sort the mapping file
os.system('samtools sort %s_consensus_mapping.sam -O sam -T temp -o  mapping.sorted.sam'%ref_file)

#Extract coverage information from the SORTED mapping file (in SAM format)

os.system('samtools depth -a mapping.sorted.sam > samtools_depth.txt')

#analyzing coverage data with pandas
data = pd.read_csv("samtools_depth.txt",sep='\t',names=['contig','position','depth'])
meancov = data.groupby(['contig'])['depth'].mean()
contig_length = data.groupby(['contig'])['position'].max()
result = pd.concat([contig_length,meancov],axis=1)
result.rename(columns = {'position':'length','depth':'average_coverage'},inplace=True)
result.to_csv("contig_statistics.txt",sep="\t")
result
# calculating average coverage across the genome
len_cov = result['length']*result['average_coverage']
genome_coverage = len_cov.sum()/result['length'].sum()

#cleaning up temp files

os.system('rm smalt_mapping.sam')
os.system('rm temp*')

# plotting with bokeh
output_file("contig_coverage.html")

# create a new plot with a title and axis labels
#creating a new variable to highlight short contigs in plot
result['color'] = ["navy" if int(i) > 1000 else "red" for i in result['length']]

Source = ColumnDataSource(result)
TOOLTIPS = [("Length", "@length"),("Coverage","@average_coverage"),("contig","@contig")]
p = figure(width = 800, height= 800, title="overall average coverage = %s"%genome_coverage, tooltips=TOOLTIPS, x_axis_label='contig size', y_axis_label='coverage',x_range=(0,max(result['length'])),y_range=(0,1.1*max(result['average_coverage'])))

# add a line renderer with legend and line thickness
p.scatter(x='length', y='average_coverage',source = Source, size=10,color="color",alpha=0.5)


low_box = BoxAnnotation(top=float(genome_coverage)/3, fill_alpha=0.1, fill_color='red')
mid_box = BoxAnnotation(bottom=float(genome_coverage)/3, top=genome_coverage*10, fill_alpha=0.1, fill_color='green')
high_box = BoxAnnotation(bottom=float(genome_coverage)*10, fill_alpha=0.1, fill_color='red')

#labels = LabelSet(x='length', y='average_coverage', text='contig', level='glyph',
#              x_offset=5, y_offset=5, source=ColumnDataSource(result), render_mode='canvas')

p.add_layout(low_box)
p.add_layout(mid_box)
p.add_layout(high_box)
#p.add_layout(labels)

# show the results
show(p)
