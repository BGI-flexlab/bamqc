package org.bgi.flexlab.bamqc;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

import java.io.*;
import java.util.HashMap;
import java.util.Map;
import java.util.zip.GZIPInputStream;

public class ReferencePanelSite {
    Map<String, String> site_vcf_map;
    int n_sites = 0;
    int n_sites_covered = 0;


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
            if(field.length == 2)
                continue;
            site_vcf_map.put(field[0], field[1]);
        }
        bufferReader.close();
    }

    public void count_site_covered(String chr, int[] coverage) {
        String site_vcf = site_vcf_map.get(chr);
        VCFFileReader reader = new VCFFileReader(new File(site_vcf), false);

        for (VariantContext vc : reader) {
            if (coverage[vc.getStart()] > 0)
                n_sites_covered++;
            n_sites++;
        }
    }


    public float getEffectiveCoverage() {
        return (float) this.n_sites_covered/this.n_sites;
    }
}