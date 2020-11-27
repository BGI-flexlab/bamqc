package org.bgi.flexlab.bamqc;


public class Main {
    public static void main(String[]args){
        Options options = new Options(args);
        ReferencePanelSite rps = new ReferencePanelSite(options.getSiteVcfList());
        BamStats bamStats = new BamStats(options.getInfile(), rps, options.isCountSecondaryReads());
        bamStats.run();
        System.out.println(bamStats.report());
    }
}

