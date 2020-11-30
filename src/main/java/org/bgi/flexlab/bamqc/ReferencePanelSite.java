package org.bgi.flexlab.bamqc;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

import java.io.*;
import java.util.HashMap;
import java.util.Map;

public class ReferencePanelSite {
    /*
    we can compute the fraction of those n sites covered by at least one read fcovered
    and compute that sample’s effective coverage λeff = − ln(1 − fcovered)
     */
    Map<String, String> site_vcf_map;
    long n_known_sites = 0;
    long n_known_sites_covered = 0;

    public ReferencePanelSite(String site_vcf_list) {
        this.site_vcf_map = new HashMap<>();
        try {
            read_vcf_list(site_vcf_list);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private void read_vcf_list(String site_vcf_list) throws IOException {
        FileInputStream inputStream = new FileInputStream(new File(site_vcf_list));
        InputStreamReader inputReader = new InputStreamReader(inputStream);
        BufferedReader bufferReader = new BufferedReader(inputReader);

        String line;
        while ((line = bufferReader.readLine()) != null) {
            String[] field = line.trim().split("\t");
            if(field.length != 2)
                continue;
            site_vcf_map.put(field[0], field[1]);
        }
        bufferReader.close();
    }

    public void count_site_covered(String chr, boolean[] coverage) {
        if(!site_vcf_map.containsKey(chr)) return;

        String site_vcf = site_vcf_map.get(chr);
        VCFFileReader reader = new VCFFileReader(new File(site_vcf), false);

        for (VariantContext vc : reader) {
            if (coverage[vc.getStart()])
                n_known_sites_covered++;
            n_known_sites++;
        }

        site_vcf_map.remove(chr);
    }

    public void count_site_uncover_chrom() {
        VCFFileReader reader;

        for(String chr: site_vcf_map.keySet()){
            reader = new VCFFileReader(new File(site_vcf_map.get(chr)), false);
            for (VariantContext vc : reader) {
                n_known_sites++;
            }
        }
    }

    public double getEffectiveCoverage() {
        double fcovered = (double) this.n_known_sites_covered/this.n_known_sites;
        return - Math.log(1-fcovered);
    }
}
