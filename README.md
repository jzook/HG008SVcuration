# Curation of structural variants in HG008-T

The issues in this repository contain curations of all HG008-T somatic SVs, mostly in GRCh38 coordinates except for those that do not have GRCh38 coordinates. The SV curation process generally followed was:
1. Created issues for all minda SVs from multiple technologies and variant callers
2. Each person gets assigned issues to curate in ribbon/splitthreader and IGV with the goal of noting a few things:
    1. Click on the link to Ribbon, or if you already have Ribbon/splitthreader open, it is faster to go to that tab instead. In splitthreader, go to the variant analysis tab and click on the current variant, which should open the illumina coverage plot for the region of interest. Take a screenshot at a zoom level that let's you easily see the entire variant. Zoom out as well to see if any other SVs/CNAs exist in the region that suggest something more complex. Then click "jump to variant in Ribbon" on the right side to see support of reads in Ribbon. Click on one of the reads that cross the breakpoint and contain the SV and take screenshot of the top and bottom view. Note the following things:
        1. Does an SV similar to this one exist within a few kb?
            1. if not, what is the correct SV, if it is not already on our list
        2. If there is true translocation, does it cause a CNA? If so, what is it?
        3. Zoom out in splitthreader to at least 5Mbp to see if there are other changes in coverage and/or SV calls that suggest this variant might be complex. If the variant is complex but has well-defined breakpoints, you can still try to note these, but we may not include this in the initial benchmark.
        4. Is this inside a larger CNA or LOH region? The mean coverage of diploid regions in splitthreader is ~200x
    2. Next, go to IGV
        1. If it's an accurate SV call, are the breakpoint(s) exactly correct?
        2. If it is an insertion, BLAT the inserted sequence in IGV (e.g., by right clicking an insertion in a read in HiFi) to see if it comes from somewhere else in the genome or is a tandem duplication
    3. Does this SV seem to be in all the cells? It's ok if this is hard to tell at this stage, but it's helpful to add tags if it seems likely to be in all the cells or is clearly mosaic
    4. Copy in any screenshots that support the above answers
    5. Look at the list of tags and add any relevant tags to the issue
    6. Look at the  vcf fields in a Google sheet. Most GitHub calls should have a corresponding row in this sheet. If so, make any modifications to the cells and highlight any changed cells in yellow. If it appears exactly correct, then highlight the row in green. If no row exists in the sheet for this variant, then create a new row for it, fill out the fields as much as possible, and highlight the entire row in yellow. If it is complex, then you can try to add lines for each breakend but highlight them in blue. Add a link to the github issue in the last column. Highlight any likely FP rows in red
        1. This table is used to make a draft benchmark VCF
