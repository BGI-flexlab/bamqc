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
    int maxReadSize;
    int minReadSize;

    long numSecondaryAlignments;
    long totalReads;
    long totalBases;
    long alignedReads;
    long duplicatedReads;

    public BamStats(String bamFile, ReferencePanelSite referencePanelSite, boolean countSecondaryReads) {
        this.bamFile = bamFile;
        this.referencePanelSite = referencePanelSite;
        this.countSecondaryReads = countSecondaryReads;
    }

    public void run(){

        long startTime = System.currentTimeMillis();
        SamReader reader = SamReaderFactory.makeDefault().open(new File(bamFile));


        // org.bioinfo.ntools.process header
        String lastActionDone = "Loading sam header...";
        System.out.println(lastActionDone);
        SAMFileHeader header = reader.getFileHeader();

        try {
            SAMFileHeader.SortOrder sortOrder = header.getSortOrder();
            if (sortOrder != SAMFileHeader.SortOrder.coordinate) {
                System.out.println("[WARN] According to header the BAM file is not sorted by coordinate!");
            }
        } catch (IllegalArgumentException ex) {
            System.out.println("[WARN] Non-standard header SortOrder value!");
        }

        int chr_len;
        int[] coverage = new int[0];

        String pre_chr = "";
        for (SAMRecord read : reader) {
            if(read.getContig() != null && !read.getContig().equals(pre_chr)){
                chr_len = header.getSequenceDictionary().getSequence(read.getContig()).getSequenceLength();
                if(!pre_chr.isEmpty())
                    referencePanelSite.count_site_covered(read.getContig(), coverage);
                coverage = new int[chr_len];
                pre_chr = read.getContig();
            }

            if (read.isSecondaryOrSupplementary()) {
                numSecondaryAlignments++;
                if (!countSecondaryReads) continue;
            }

            //compute read size
            int readSize = read.getReadLength();
            totalBases += readSize;
            if (readSize > maxReadSize) {
                maxReadSize = readSize;
            }
            if (readSize < minReadSize) {
                minReadSize = readSize;
            }
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

            for (int i = read.getAlignmentStart(); i <= read.getAlignmentEnd(); i++) {
                coverage[i]++;
            }
        }
        referencePanelSite.count_site_covered(pre_chr, coverage);

        long overallTime = System.currentTimeMillis();
        System.out.println("Overall analysis time: " + (overallTime - startTime) / 1000);
    }

    public long getTotalReads() {
        return totalReads;
    }

    public long getTotalBases() {
        return totalBases;
    }

    public float getAlignedReadsRatio() {
        return (float) alignedReads/totalReads;
    }

    public float getDuplicatedReadsRatio() {
        return (float) duplicatedReads/totalReads;
    }

    public String report() {
        StringBuilder buf = new StringBuilder();
        buf.append("Total Reads: ").append(alignedReads).append("\n");
        buf.append("Total Bases: ").append(totalBases).append("\n");
        buf.append("Aligned Reads ratio: ").append(getAlignedReadsRatio()).append("\n");
        buf.append("Duplicated Reads ratio: ").append(getDuplicatedReadsRatio()).append("\n");
        buf.append("Known Sites Covered: ").append(referencePanelSite.n_sites_covered).append("\n");
        buf.append("Effective Coverage: ").append(referencePanelSite.getEffectiveCoverage()).append("\n");
        return buf.toString();
    }

    private Pair<List<String>, List<String>> getReportResult() {
        final List<String> names = new ArrayList<>();
        final List<String> values = new ArrayList<>();
        names.add("Total Reads");
        values.add(Long.toString(totalReads));
        names.add("Total Bases");
        values.add(Long.toString(totalBases));
        names.add("Aligned Reads ratio");
        values.add(StatsUtils.divide(alignedReads, totalReads));
        names.add("Duplicated Reads ratio");
        values.add(StatsUtils.divide(duplicatedReads, totalReads));
        names.add("Known Sites Covered");
        values.add(Long.toString(referencePanelSite.n_sites_covered));
        names.add("Effective Coverage");
        values.add(StatsUtils.divide(referencePanelSite.n_sites_covered, referencePanelSite.n_sites));
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
}
