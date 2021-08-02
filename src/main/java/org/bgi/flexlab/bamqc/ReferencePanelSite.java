package org.bgi.flexlab.bamqc;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

import java.io.*;
import java.util.*;

public class ReferencePanelSite {
    /*
    we can compute the fraction of those n sites covered by at least one read fcovered
    and compute that sample’s effective coverage λeff = − ln(1 − fcovered)
     */
    Map<String, List<String>> site_vcf_map;
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
            String region = field[0];
            String sites_file = field[1];
            if(field.length != 2){
                System.err.println("knowsites.vcf.list is bad!");
                System.exit(1);
            }
            if(region.contains(":")){
                String chr = region.split(":")[0];
                if(!site_vcf_map.containsKey(chr)){
                    List<String> sites = new ArrayList<>();
                    sites.add(sites_file);
                    site_vcf_map.put(chr, sites);
                }else
                    site_vcf_map.get(chr).add(sites_file);

            }else {
                if(!site_vcf_map.containsKey(region)){
                    List<String> sites = new ArrayList<>();
                    sites.add(sites_file);
                    site_vcf_map.put(region, sites);
                }else
                    site_vcf_map.get(region).add(sites_file);
            }
        }
        bufferReader.close();
    }

    public void count_site_covered(String chr, boolean[] coverage) {
        if(!site_vcf_map.containsKey(chr)) return;

        List<String> site_vcfs = site_vcf_map.get(chr);
        Set<Integer> positions = new HashSet<>();

        for(String site_vcf : site_vcfs){
            VCFFileReader reader = new VCFFileReader(new File(site_vcf), false);
            for (VariantContext vc : reader) {
                if(!chr.equals(vc.getContig())){
                    System.err.println("[ERROR] site_vcf (" + site_vcf + ") has different chrom : " + vc.getContig());
                    System.exit(2);
                }

                if(positions.contains(vc.getStart()))
                    continue;
                else
                    positions.add(vc.getStart());

                if (coverage[vc.getStart()])
                    n_known_sites_covered++;
                n_known_sites++;
            }
        }
        site_vcf_map.remove(chr);
    }

    public void count_site_uncover_chrom() {
        VCFFileReader reader;

        for(String chr: site_vcf_map.keySet()){
            List<String> site_vcfs = site_vcf_map.get(chr);
            for(String site_vcf : site_vcfs) {
                reader = new VCFFileReader(new File(site_vcf), false);
                for (VariantContext vc : reader) {
                    n_known_sites++;
                }
            }
        }
    }

    public double getEffectiveCoverage() {
        double fcovered = (double) this.n_known_sites_covered/this.n_known_sites;
        return - Math.log(1-fcovered);
    }
}
