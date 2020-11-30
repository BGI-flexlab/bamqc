package org.bgi.flexlab.bamqc;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
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
        SamReader reader = SamReaderFactory.makeDefault().open(new File(bamFile));

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
        boolean[] coverage = null;

        String pre_chr = "";
        for (SAMRecord read : reader) {
            if(read.getContig() != null && !read.getContig().equals(pre_chr)){
                chr_len = header.getSequenceDictionary().getSequence(read.getContig()).getSequenceLength();
                if(!pre_chr.isEmpty()) {
                    referencePanelSite.count_site_covered(pre_chr, coverage);
                    count_coverage(coverage);
                    System.out.println("Processing : " + pre_chr);
                }
                coverage = new boolean[chr_len+1];
                pre_chr = read.getContig();
            }

            if (read.isSecondaryOrSupplementary()) {
                numSecondaryAlignments++;
                if (!countSecondaryReads) continue;
            }

            //compute read size
            int readSize = read.getReadLength();
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
                }
            }
        }
        if(coverage != null){
            count_coverage(coverage);
            referencePanelSite.count_site_covered(pre_chr, coverage);
        }
        referencePanelSite.count_site_uncover_chrom();
        System.out.println("Processing : " + pre_chr);

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
        names.add("Total Known Sites");
        values.add(Long.toString(referencePanelSite.n_known_sites));
        names.add("Known Sites Covered");
        values.add(Long.toString(referencePanelSite.n_known_sites_covered));
        names.add("Effective Coverage");
        values.add(StatsUtils.realFormat(referencePanelSite.getEffectiveCoverage(), 2));
        names.add("Coverage");
        values.add(StatsUtils.divide(n_sites_covered, referenceLength));
        names.add("Aligned Reads ratio");
        values.add(StatsUtils.divide(alignedReads, totalReads));
        names.add("Duplicated Reads ratio");
        values.add(StatsUtils.divide(duplicatedReads, totalReads));
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
        for (boolean b : coverage) {
            if (b) n_sites_covered++;
        }
    }
}
