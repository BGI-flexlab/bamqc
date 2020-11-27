package org.bgi.flexlab.bamqc;

import htsjdk.samtools.*;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

public class BamStats {

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
}
