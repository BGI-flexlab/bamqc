package org.bgi.flexlab.bamqc;


import java.io.IOException;

public class Main {
    public static void main(String[]args){
        Options options = new Options(args);
        ReferencePanelSite rps = new ReferencePanelSite(options.getSiteVcfList());
        BamStats bamStats = new BamStats(options.getInfile(), rps, options.isCountSecondaryReads());
        System.out.println("Start ...");
        bamStats.run();
        try {
            bamStats.writeReport(options.getOutfile(), options.getAppVersion());
        } catch (IOException e) {
            e.printStackTrace();
        }
        System.out.println("Done");
    }
}

