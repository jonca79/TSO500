#!/bin/python3.6
import sys
from pysam import VariantFile

vcf_in = VariantFile(sys.argv[1])  #dosen't matter if bgziped or not. Automatically recognizes

new_header = vcf_in.header
# import pdb; pdb.set_trace()
new_header.info.add("DP","1","Integer","Sum of AD fields")
new_header.info.add("AF","1","Float","Alt AD / sum(AD)")
#start new vcf with the new_header
vcf_out = VariantFile(sys.argv[2], 'w', header=new_header)

for record in vcf_in.fetch():
    dp = sum(record.samples[0].get("AD"))
    record.info["DP"]=dp
    af=record.samples[0].get("AD")[1]/dp
    record.info["AF"]=af
    vcf_out.write(record)
