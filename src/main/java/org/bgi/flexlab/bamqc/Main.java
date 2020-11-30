package org.bgi.flexlab.bamqc;


import java.io.IOException;

public class Main {
    public static void main(String[]args){
        Options options = new Options(args);
        ReferencePanelSite rps = new ReferencePanelSite(options.getSiteVcfList());
        BamStats bamStats = new BamStats(options.getInfile(), rps, options.isCountSecondaryReads());
        bamStats.run();
        try {
            bamStats.writeReport(options.getOutfile());
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}

