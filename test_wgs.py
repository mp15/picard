#!/usr/bin/env python

import os

N = 1*10**5

def run(c, return_output=False, return_exit_code=False):
	print("==>  %s" % c)
	if return_output:
		import commands
		exit_code, result = commands.getstatusoutput(c)
	else:
		result = os.system(c)

	if return_output or return_exit_code:
		print("Result: %(result)s \n" % locals())
		return result

outputs = ["metrics.txt", "metrics_fast.txt"]
outputs0, outputs1 = outputs

for output in outputs:
	if os.path.isfile(output):
		os.remove(output)
input_sam = "~/tng/wgs_bams/B084HABXX.3.aligned.duplicates_marked.bam"
input_sam = "~/tng/wgs_bams/B084HABXX.3.aligned.duplicates_marked.subset.sam"

print("\n\n------------\n\n")
run("picard CollectWgsMetrics I= %(input_sam)s O=%(outputs0)s R= ~/hg19.fa STOP_AFTER=%(N)d" % locals())
print("\n\n------------\n\n")
run("picard CollectFastWgsMetrics I= %(input_sam)s O=%(outputs1)s R= ~/hg19.fa STOP_AFTER=%(N)d" % locals())

assert run("diff -I 'Started\|VALIDATION_STRINGENCY\|picard.analysis.Collect' metrics.txt metrics_fast.txt", return_exit_code=True) == 0
assert int(run("wc -l metrics.txt | cut -f 1 -d \\ ", return_output=True)) == 263

print("Success!!")
