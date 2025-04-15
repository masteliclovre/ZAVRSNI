#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <string>
#include <cctype>
#include <algorithm>
#include <cmath>

double phred_to_prob(int phred);
double prob_to_phred(double prob);

struct Variant {
    std::string chrom;
    int pos;
    std::string ref;
    std::string alt;
    int coverage;
    std::string bases;
    std::string qualities;
    double freq;
    int AC;
    int AN;
    std::vector<int> AD;
    double quality;
};

double phred_to_prob(int phred) {
    return std::pow(10.0, -phred / 10.0);
}

double prob_to_phred(double prob) {
    return (prob > 0) ? -10.0 * std::log10(prob) : 99.0;
}

void parse_var(Variant& var, std::ofstream& vcf_output_file, double threshold) {
    std::unordered_map<std::string, int> variant_counts;
    std::unordered_map<std::string, double> variant_qualities;

    bool isDel = false;
    size_t indel_size = 0;
    std::string indel_sequence;
    size_t i = 0;
    size_t index = 0;
    size_t total_bases = 0;

    while(i < var.bases.size()){
        char base = var.bases[i];
        char qual = (index < var.qualities.size()) ? var.qualities[index] : '!';
        int phred = (qual != '!') ? qual - 33 : 30;

        if(base == '^') {
            i += 2;
            continue;
        }
        if(base == '$') {
            i++;
            continue;
        }
        if(base == '*' || base == '<' || base == '>') {
            i++;
            index++;
            continue;
        }

        if(base == '+' || base == '-') {
            isDel = (base == '-');
            i++;

            // Parse indel size
            indel_size = 0;
            while(i < var.bases.size() && isdigit(var.bases[i])) {
                indel_size = indel_size * 10 + (var.bases[i] - '0');
                i++;
            }

            // Parse indel sequence
            indel_sequence.clear();
            for(size_t j = 0; j < indel_size && i < var.bases.size(); j++, i++) {
                indel_sequence.push_back(std::toupper(var.bases[i]));
            }

            std::string key = (isDel ? "-" : "+") + indel_sequence;
            variant_counts[key]++;
            variant_qualities[key] += phred;
            total_bases++;
            index++;
            continue;
        }

        if(base == '.' || base == ',') {
            variant_counts[var.ref]++;
            variant_qualities[var.ref] += phred;
            total_bases++;
            index++;
        }
        else if(isalpha(base)) {
            std::string alt(1, std::toupper(base));
            variant_counts[alt]++;
            variant_qualities[alt] += phred;
            total_bases++;
            index++;
        }
        i++;
    }

    // Create alleles vector
    std::vector<std::pair<std::string, double>> alleles;
    for(const auto& entry : variant_counts) {
        const std::string& base = entry.first;
        int count = entry.second;
        double freq = static_cast<double>(count) / total_bases;
        double avg_qual = variant_qualities[base] / count;
        alleles.emplace_back(base, freq * avg_qual);
    }

    // Sort alleles by score descending
    std::sort(alleles.begin(), alleles.end(),
        [](const std::pair<std::string, double>& a, const std::pair<std::string, double>& b) {
            return a.second > b.second;
        });

    var.alt = var.ref;
    var.freq = 0.0;
    var.AC = 0;
    var.AN = 2;
    var.AD = {0, 0};
    var.quality = 0.0;

    if(!alleles.empty() && alleles[0].first != var.ref){
        var.alt = alleles[0].first;
        var.freq = static_cast<double>(variant_counts[var.alt]) / total_bases;
        var.AD = {variant_counts[var.ref], variant_counts[var.alt]};

        if(var.freq >= 0.9){
            var.AC = 2;
        }else if(var.freq >= threshold){
            var.AC = 1;
        }

        var.quality = variant_qualities[var.alt] / variant_counts[var.alt];
    }

    std::string genotype;
    if(var.AC > 0){
        if(var.AC == 2){
            genotype = "1/1";
        }else{
            genotype = "0/1";
        }
    }

    if(var.alt != var.ref){
        vcf_output_file << var.chrom << "\t" << var.pos << "\t" << '.' << "\t" << var.ref << "\t" << var.alt << "\t" 
                        << var.quality << "\t" << "PASS" << "\t" 
                        << "DP=" << var.coverage << ";AC=" << var.AC << ";AN=" << var.AN << ";AF=" << var.freq << "\t" 
                        << "GT:DP:AD:PL" << "\t";
        
        double p_ref = phred_to_prob(var.quality);
        double p_het = 0.5;
        double p_hom = 1.0 - p_het;

        int pl_ref = static_cast<int>(prob_to_phred(1.0 - p_ref));
        int pl_het = static_cast<int>(prob_to_phred(1.0 - p_het));
        int pl_hom = static_cast<int>(prob_to_phred(1.0 - p_hom));

        vcf_output_file << genotype << ":" << var.coverage << ":" 
                        << var.AD[0] << "," << var.AD[1] << ":" 
                        << pl_ref << "," << pl_het << "," << pl_hom << "\n";
    }
}

int main(int argc, char* argv[]) {
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <pileup_file> <vcf_output_file> <threshold>\n";
        return 1;
    }

    std::string pileup_file = argv[1];
    std::string vcf_output_file = argv[2];
    double threshold = std::stod(argv[3]);

    std::ifstream in(pileup_file);
    if (!in) {
        std::cerr << "Could not open " << pileup_file << "\n";
        return 1;
    }

    std::ofstream vcf_output(vcf_output_file);
    if (!vcf_output) {
        std::cerr << "Could not create " << vcf_output_file << "\n";
        return 1;
    }

    // VCF header
    vcf_output << "##fileformat=VCFv4.2\n";
    vcf_output << "##FILTER=<ID=PASS,Description=\"All filters passed\">\n";
    vcf_output << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n";
    vcf_output << "##INFO=<ID=AC,Number=1,Type=Integer,Description=\"Allele Count\">\n";
    vcf_output << "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total Allele Number\">\n";
    vcf_output << "##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency\">\n";
    vcf_output << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
    vcf_output << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n";
    vcf_output << "##FORMAT=<ID=AD,Number=2,Type=Integer,Description=\"Allelic Depths\">\n";
    vcf_output << "##FORMAT=<ID=PL,Number=3,Type=Integer,Description=\"Phred-scaled Genotype Likelihoods\">\n";
    vcf_output << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n";

    std::string line;
    while (std::getline(in, line)) {
        std::stringstream ss(line);
        std::string chrom;
        int pos;
        std::string ref;
        int coverage;
        std::string bases;
        std::string qualities;

        ss >> chrom >> pos >> ref >> coverage >> bases >> qualities;
        if (coverage > 0) {
            Variant var = {chrom, pos, ref, ref, coverage, bases, qualities, 0.0, 0, 0};
            std::transform(var.ref.begin(), var.ref.end(), var.ref.begin(), [](unsigned char c) { return std::toupper(c); });
            parse_var(var, vcf_output, threshold);
        }
    }

    std::cout << "VCF file generated: " << vcf_output_file << "\n";
    return 0;
}