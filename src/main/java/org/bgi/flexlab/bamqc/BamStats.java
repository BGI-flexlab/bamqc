package org.bgi.flexlab.bamqc;

import htsjdk.samtools.*;
import org.bgi.flexlab.bamqc.util.Pair;
import org.bgi.flexlab.bamqc.util.StatsUtils;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class BamStats {
    final private static String REPORT_HEADER = "## BGI-lowpass bam quality control, version 1.0\n";

    private String bamFile;
    private ReferencePanelSite referencePanelSite;

    private boolean countSecondaryReads;

    // read size
//    int maxReadSize;
//    int minReadSize;
    long referenceLength;
    long n_sites_covered = 0;
    long n_bases_mapped = 0;
    long n_bases_mapped_chrX = 0;
    long n_bases_mapped_chrY = 0;
    double chrY_cov = 0.0;
    double XY_depth_ratio = 0.0;
    long numSecondaryAlignments = 0;
    long totalReads = 0;
    long totalBases = 0;
    long alignedReads = 0;
    long duplicatedReads = 0;

    public BamStats(String bamFile, ReferencePanelSite referencePanelSite, boolean countSecondaryReads) {
        this.bamFile = bamFile;
        this.referencePanelSite = referencePanelSite;
        this.countSecondaryReads = countSecondaryReads;
    }

    public void run(){

        long startTime = System.currentTimeMillis();
        SamReader reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(bamFile));

        SAMFileHeader header = reader.getFileHeader();

        try {
            SAMFileHeader.SortOrder sortOrder = header.getSortOrder();
            if (sortOrder != SAMFileHeader.SortOrder.coordinate) {
                System.out.println("[WARN] According to header the BAM file is not sorted by coordinate!");
            }
        } catch (IllegalArgumentException ex) {
            System.out.println("[WARN] Non-standard header SortOrder value!");
        }

        referenceLength = header.getSequenceDictionary().getReferenceLength();

        int chr_len;
        int chrX_len;
        int chrY_len = 0;
        double chrX_depth = 0;
        double chrY_depth = 0;
        boolean[] coverage = null;
        boolean is_chrX = false;
        boolean is_chrY = false;

        String pre_chr = "";
        for (SAMRecord read : reader) {
            if(read.getContig() != null && !read.getContig().equals(pre_chr)){
                chr_len = header.getSequenceDictionary().getSequence(read.getContig()).getSequenceLength();
                if(!pre_chr.isEmpty()) {
                    referencePanelSite.count_site_covered(pre_chr, coverage);
                    for (boolean b : coverage) if (b) n_sites_covered++;
                    System.out.println("Processing : " + pre_chr);
                    if (pre_chr.endsWith("X")){
                        chrX_len = header.getSequenceDictionary().getSequence(pre_chr).getSequenceLength();
                        chrX_depth = (double)n_bases_mapped_chrX / chrX_len;
                    }else if (pre_chr.endsWith("Y")){
                        chrY_len = header.getSequenceDictionary().getSequence(pre_chr).getSequenceLength();
                        chrY_depth = (double)n_bases_mapped_chrY / chrY_len;
                        chrY_cov = count_coverage_chrY(coverage);
                    }
                }
                coverage = new boolean[chr_len+1];
                System.gc();
                pre_chr = read.getContig();
                is_chrX = pre_chr.endsWith("X");
                is_chrY = pre_chr.endsWith("Y");
            }

            int readSize = read.getReadLength();
            if(is_chrX) n_bases_mapped_chrX +=readSize;
            if(is_chrY) n_bases_mapped_chrY +=readSize;

            if (read.isSecondaryOrSupplementary()) {
                numSecondaryAlignments++;
                if (!countSecondaryReads) continue;
            }

            //compute read size
            totalBases += readSize;
//            if (readSize > maxReadSize) {
//                maxReadSize = readSize;
//            }
//            if (readSize < minReadSize) {
//                minReadSize = readSize;
//            }
            totalReads++;

            if (read.getDuplicateReadFlag()) {
                duplicatedReads++;
                continue;
            }

            // accumulate only mapped reads
            if (read.getReadUnmappedFlag()) {
                continue;
            }
            alignedReads++;

            if(coverage != null){
                for (int i = read.getAlignmentStart(); i <= read.getAlignmentEnd(); i++) {
                    coverage[i]=true;
                    n_bases_mapped+=1;
                }
            }
        }
        if(coverage != null){
            for (boolean b : coverage) if (b) n_sites_covered++;
            referencePanelSite.count_site_covered(pre_chr, coverage);
            if(is_chrY) {
                chrY_cov = count_coverage_chrY(coverage);
            }
        }
        referencePanelSite.count_site_uncover_chrom();
        System.out.println("Processing : " + pre_chr);

        if(chrY_depth != 0){
            XY_depth_ratio = chrX_depth / chrY_depth;
        }else {
            XY_depth_ratio = Double.MAX_VALUE;
        }

        long overallTime = System.currentTimeMillis();
        System.out.println("Overall analysis time: " + (overallTime - startTime) / 1000 + " s");
    }

    private Pair<List<String>, List<String>> getReportResult() {
        final List<String> names = new ArrayList<>();
        final List<String> values = new ArrayList<>();
        names.add("Total Reads");
        values.add(Long.toString(totalReads));
        names.add("Total Bases");
        values.add(Long.toString(totalBases));
        names.add("Bases Deduplicated and Aligned");
        values.add(Long.toString(n_bases_mapped));
        names.add("Known Sites");
        values.add(Long.toString(referencePanelSite.n_known_sites));
        names.add("Known Sites Covered");
        values.add(Long.toString(referencePanelSite.n_known_sites_covered));
        names.add("Effective Coverage");
        values.add(StatsUtils.realFormat(referencePanelSite.getEffectiveCoverage(), 2));
        names.add("Coverage 1X");
        values.add(StatsUtils.divide(n_sites_covered, referenceLength));
        names.add("Mean Depth");
        values.add(StatsUtils.divide(n_bases_mapped, referenceLength));
        names.add("Aligned Reads ratio");
        values.add(StatsUtils.divide(alignedReads, totalReads));
        names.add("Duplicated Reads ratio");
        values.add(StatsUtils.divide(duplicatedReads, totalReads));
        names.add("Sex");
        String sex_info = String.format("(%s,%s)", StatsUtils.realFormat(XY_depth_ratio, 2),
                StatsUtils.realFormat(chrY_cov, 2));
        if(XY_depth_ratio < 2)
            values.add("M "+sex_info);
        else if (XY_depth_ratio > 4.5)
            values.add("F "+sex_info);
        else {
            double cov = (double)n_sites_covered / referenceLength;
            double chrY_cov_ratio = cov/chrY_cov;
            if(chrY_cov_ratio > 10)
                values.add("F "+sex_info);
            else
                values.add("M "+sex_info);
//            values.add("failure ("+StatsUtils.realFormat(XY_depth_ratio, 2) + ")");
        }

        return Pair.create(names, values);
    }

    public String getReport() {
        Pair<List<String>, List<String>> pair = getReportResult();
        List<String> names = pair.getFirst();
        List<String> values = pair.getSecond();
        StringBuilder outString = new StringBuilder();

        for (int i = 0; i < names.size(); i++) {
            outString.append(String.format("%-37s", names.get(i)));
            outString.append(":  ");
            outString.append(values.get(i));
            outString.append("\n");
        }

        return outString.toString();
    }


    public void writeReport(String outfile) throws IOException {
        File report = new File(outfile);
        FileWriter fileWritter = new FileWriter(report, false);
        fileWritter.write(REPORT_HEADER);
        fileWritter.write(getReport());
        fileWritter.close();
    }

    public void count_coverage(boolean[] coverage){
        for (boolean b : coverage) if (b) n_sites_covered++;
    }

    public double count_coverage_chrY(boolean[] coverage){
        int n_sites_covered_chrY = 0;
        int chrY_size = coverage.length;
        int i=2781479;
        while (i < chrY_size - 330000) {
            if (coverage[i]) n_sites_covered_chrY++;
            i++;
        }
        return (double)n_sites_covered_chrY / (chrY_size-2781479-330000);
    }
}
