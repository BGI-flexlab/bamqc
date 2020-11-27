package org.bgi.flexlab.bamqc;

import org.apache.commons.cli.*;

import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.Properties;

public class Options {

    final String SOFTWARE_NAME = "vcfstats";
    private boolean countSecondaryReads = false;
    private String infile;
    private String outdir = null;
    private String dbsnp;
    private Float qual = 0.0f;
    private String siteVcfList;


    SimpleDateFormat df = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");//设置日期格式
    String compile_date = df.format(new Date());

    org.apache.commons.cli.Options options = new org.apache.commons.cli.Options();
    CommandLine cmdLine;
    CommandLineParser parser = new DefaultParser();

    public Options() {
    }

    public Options(String[] args) {
        parse(args);
    }

    private String helpHeader() {
        StringBuilder sb = new StringBuilder();
        sb.append("\nVersion     : ");
        sb.append(getAppVersion());
        sb.append("\nAuthor      : huangzhibo@genomics.cn");
        sb.append("\nCompile Date: ");
        sb.append(compile_date);
        sb.append("\nNote        : BGI-lowpass bam quality control\n");
        sb.append("\nOptions:\n");
        return sb.toString();
    }

    public void parse(String[] args) {
        String header = helpHeader();
        String footer = "\nPlease report issues at https://github.com/BGI-flexlab/bamqc/issues";

        options.addOption(Option.builder("i")
                .longOpt("input")
                .required(true)
                .hasArg()
                .argName("FILE")
                .desc("input bam(BAM). [request]")
                .build());
        options.addOption(Option.builder("o")
                .longOpt("output")
                .hasArg()
                .argName("String")
                .desc("report outdir [request]")
                .build());
        options.addOption(Option.builder("s")
                .longOpt("site")
                .required(true)
                .hasArg()
                .argName("FILE")
                .desc("reference panel site.vcfs list [request]")
                .build());
        options.addOption(Option.builder("c")
                .longOpt("countSecondaryReads")
                .desc("The secondary alignment reads are counted and ignored by default [false]")
                .build());
        options.addOption(Option.builder("h")
                .longOpt("help")
                .desc("Print this help.")
                .build());

        HelpFormatter formatter = new HelpFormatter();
        formatter.setWidth(2 * HelpFormatter.DEFAULT_WIDTH);

        try {
            cmdLine = parser.parse(options, args);
            if (cmdLine.hasOption("h")) {
                formatter.printHelp("java -jar " + SOFTWARE_NAME + ".jar", header, options, footer, true);
                System.exit(0);
            }
        } catch (ParseException e) {
            formatter.printHelp("java -jar " + SOFTWARE_NAME + ".jar", header, options, footer, true);
            System.exit(0);
        }

        infile = cmdLine.getOptionValue("input");

        if (cmdLine.hasOption("output")) {
            outdir = cmdLine.getOptionValue("output");
        }

        if (cmdLine.hasOption("qual")) {
            qual = Float.parseFloat(cmdLine.getOptionValue("qual", "0.0"));
        }

        if (cmdLine.hasOption("site")) {
            siteVcfList = cmdLine.getOptionValue("site");
        }

        if (cmdLine.hasOption("c")) {
            countSecondaryReads = cmdLine.hasOption("countSecondaryReads");
        }
    }

    public String getInfile() {
        return infile;
    }

    public String getOutdir() {
        return outdir;
    }

    public String getDbsnp() {
        return dbsnp;
    }

    public Float getQual() {
        return qual;
    }

    public boolean isCountSecondaryReads() {
        return countSecondaryReads;
    }

    public static String getAppVersion() {
        String appVersion = "NA";
        Properties properties = new Properties();
        try {
            properties.load(Options.class.getClassLoader().getResourceAsStream("app.properties"));
            if (!properties.isEmpty()) {
                appVersion = properties.getProperty("app.version");
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return appVersion;
    }

    public String getSiteVcfList() {
        return siteVcfList;
    }

    /**
     * 测试
     *
     * @param args
     * @throws Exception
     */
    public static void main(String[] args) throws Exception {
        String[] arg = {"-i", "test.txt", "-i", "test2.txt"};
        Options parameter = new Options();
        parameter.parse(arg);
        System.out.println(parameter.getInfile());
    }
}

