#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <cctype>

double phred_to_prob(int phred);
double prob_to_phred(double prob);

struct Variant {
    std::string chrom;
    int pos;          // 1-based position in reference
    std::string ref;  // Reference base(s)
    std::string alt;  // Alternate base(s)
    int coverage = 0;
    std::string bases;
    std::string qualities;
    double freq = 0.0;
    int AC = 0;
    int AN = 2;
    std::vector<int> AD = {0, 0};
    double quality = 0.0;
    int count = 0;  // Number of reads supporting this variant
};

double phred_to_prob(int phred) {
    return std::pow(10.0, -phred / 10.0);
}

double prob_to_phred(double prob) {
    return (prob > 0) ? -10.0 * std::log10(prob) : 99.0;
}

class ReferenceGenome {
    std::unordered_map<std::string, std::string> chromosomes;
    
public:
    void load(const std::string& fasta_file) {
        std::ifstream in(fasta_file);
        if (!in) throw std::runtime_error("Cannot open reference file");
        
        std::string line, current_chrom, current_seq;
        
        while (std::getline(in, line)) {
            if (line.empty()) continue;
            
            if (line[0] == '>') {
                if (!current_chrom.empty()) {
                    chromosomes[current_chrom] = current_seq;
                    current_seq.clear();
                }
                current_chrom = line.substr(1);
                size_t space_pos = current_chrom.find(' ');
                if (space_pos != std::string::npos) {
                    current_chrom = current_chrom.substr(0, space_pos);
                }
            } else {
                current_seq += line;
            }
        }
        
        if (!current_chrom.empty()) {
            chromosomes[current_chrom] = current_seq;
        }
    }
    
    std::string get_sequence(const std::string& chrom, int pos, int length) {
        if (chromosomes.find(chrom) == chromosomes.end()) {
            throw std::runtime_error("Chromosome not found: " + chrom);
        }
        
        const std::string& seq = chromosomes.at(chrom);
        if (pos < 1 || pos + length - 1 > seq.length()) {
            throw std::runtime_error("Requested sequence out of bounds");
        }
        
        return seq.substr(pos - 1, length);
    }
};

class Converter {
    ReferenceGenome ref_genome;
    std::unordered_map<std::string, std::vector<Variant>> variants;
    double threshold = 0.9; // Default threshold
    int min_support = 3;    // Minimum number of reads supporting a variant
    
public:
    void process_sam(const std::string& sam_file, const std::string& ref_fasta, double thresh = 0.9) {
        threshold = thresh;
        ref_genome.load(ref_fasta);

        std::ifstream in(sam_file);
        if (!in) throw std::runtime_error("Cannot open SAM file");
         
        std::string line;
        while (std::getline(in, line)) {
            if (line[0] == '@') continue; // Skip headers
              
            std::istringstream ss(line);
            std::vector<std::string> fields;
            std::string field;
              
            while (std::getline(ss, field, '\t')) {
                fields.push_back(field);
            }
              
            if (fields.size() < 11) continue;
              
            const std::string& chrom = fields[2];
            int pos = std::stoi(fields[3]); // 1-based leftmost mapping position
            const std::string& cigar = fields[5];
            const std::string& seq = fields[9];
            const std::string& qual = fields[10];
              
            process_alignment(chrom, pos, cigar, seq, qual);
        }
    }
    
    void write_vcf(const std::string& vcf_file) {
        std::ofstream out(vcf_file);
        if (!out) throw std::runtime_error("Cannot create VCF file");
        
        // Write VCF header
        out << "##fileformat=VCFv4.2\n";
        out << "##FILTER=<ID=PASS,Description=\"All filters passed\">\n";
        out << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n";
        out << "##INFO=<ID=AC,Number=1,Type=Integer,Description=\"Allele Count\">\n";
        out << "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total Allele Number\">\n";
        out << "##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency\">\n";
        out << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
        out << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n";
        out << "##FORMAT=<ID=AD,Number=2,Type=Integer,Description=\"Allelic Depths\">\n";
        out << "##FORMAT=<ID=PL,Number=3,Type=Integer,Description=\"Phred-scaled Genotype Likelihoods\">\n";
        out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n";
        
        // Write variants
        for (const auto& [chrom, chrom_variants] : variants) {
            for (const Variant& var : chrom_variants) {
                // Only output variants with sufficient support and high enough frequency
                if (var.alt != var.ref && var.count >= min_support && var.freq >= threshold) {
                    std::string genotype = "1/1"; // We're only calling homozygous variants like the original code

                    double p_ref = phred_to_prob(var.quality);
                    double p_het = 0.5;
                    double p_hom = 1.0 - p_het;

                    int pl_ref = static_cast<int>(prob_to_phred(1.0 - p_ref));
                    int pl_het = static_cast<int>(prob_to_phred(1.0 - p_het));
                    int pl_hom = static_cast<int>(prob_to_phred(1.0 - p_hom));

                    out << chrom << "\t" << var.pos << "\t.\t" 
                        << var.ref << "\t" << var.alt << "\t" 
                        << var.quality << "\tPASS\t"
                        << "DP=" << var.coverage << ";AC=2;AN=2" 
                        << ";AF=" << var.freq << "\t"
                        << "GT:DP:AD:PL\t"
                        << genotype << ":" << var.coverage << ":" 
                        << var.AD[0] << "," << var.AD[1] << ":" 
                        << pl_ref << "," << pl_het << "," << pl_hom << "\n";
                }
            }
        }
    }
    
private:
    void process_alignment(const std::string& chrom, int pos, const std::string& cigar, 
                          const std::string& seq, const std::string& qual) {
        std::vector<std::pair<int, char>> cigar_ops = parse_cigar(cigar);
        
        int ref_pos = pos;
        int seq_pos = 0;
        int qual_pos = 0;
        
        for (const auto& [len, op] : cigar_ops) {
            switch(op) {
                case 'M': case '=': case 'X': // Match or mismatch
                    for (int i = 0; i < len; ++i) {
                        std::string ref_base = ref_genome.get_sequence(chrom, ref_pos, 1);
                        std::string alt_base(1, seq[seq_pos]);
                        
                        if (ref_base != alt_base) { // SNP
                            Variant var;
                            var.chrom = chrom;
                            var.pos = ref_pos;
                            var.ref = ref_base;
                            var.alt = alt_base;
                            var.bases = alt_base;
                            var.qualities = (qual_pos < qual.size()) ? std::string(1, qual[qual_pos]) : "!";
                            var.coverage = 1;
                            var.count = 1;
                            add_variant(var);
                        }
                        
                        ref_pos++;
                        seq_pos++;
                        qual_pos++;
                    }
                    break;
                    
                case 'I': // Insertion
                {
                    std::string ref_base = ref_genome.get_sequence(chrom, ref_pos, 1);
                    std::string inserted = seq.substr(seq_pos, len);
                    Variant var;
                    var.chrom = chrom;
                    var.pos = ref_pos;
                    var.ref = ref_base;
                    var.alt = ref_base + inserted;
                    var.bases = "+" + std::to_string(len) + inserted;
                    var.qualities = (qual_pos < qual.size()) ? qual.substr(qual_pos, len) : std::string(len, '!');
                    var.coverage = 1;
                    var.count = 1;
                    add_variant(var);
                    seq_pos += len;
                    qual_pos += len;
                    break;
                }
                    
                case 'D': // Deletion
                {
                    std::string ref_bases = ref_genome.get_sequence(chrom, ref_pos, len);
                    Variant var;
                    var.chrom = chrom;
                    var.pos = ref_pos;
                    var.ref = ref_bases;
                    var.alt = ref_bases.substr(0, 1);
                    var.bases = "-" + std::to_string(len) + ref_bases.substr(1);
                    var.qualities = (qual_pos < qual.size()) ? std::string(len, qual[qual_pos]) : std::string(len, '!');
                    var.coverage = 1;
                    var.count = 1;
                    add_variant(var);
                    ref_pos += len;
                    qual_pos++;
                    break;
                }
                    
                case 'S': case 'H': // Soft or hard clip - skip in reference
                    seq_pos += len;
                    qual_pos += len;
                    break;
                    
                case 'N': // Skipped region in reference
                    ref_pos += len;
                    break;
            }
        }
    }
    
    std::vector<std::pair<int, char>> parse_cigar(const std::string& cigar) {
        std::vector<std::pair<int, char>> operations;
        int len = 0;
        
        for (char c : cigar) {
            if (isdigit(c)) {
                len = len * 10 + (c - '0');
            } else {
                operations.emplace_back(len, c);
                len = 0;
            }
        }
        
        return operations;
    }
    
    void add_variant(Variant& var) {
        std::string key = var.chrom + ":" + std::to_string(var.pos) + ":" + var.ref + ">" + var.alt;
        
        // Parse the variant information
        parse_var(var);
        
        // Check if variant already exists
        bool found = false;
        for (auto& existing_var : variants[var.chrom]) {
            if (existing_var.pos == var.pos && existing_var.ref == var.ref && existing_var.alt == var.alt) {
                // Update counts
                existing_var.coverage += var.coverage;
                existing_var.bases += var.bases;
                existing_var.qualities += var.qualities;
                existing_var.count++;
                parse_var(existing_var); // Recalculate frequencies
                found = true;
                break;
            }
        }
        
        if (!found) {
            // Add new variant
            variants[var.chrom].push_back(var);
        }
    }
    
    void parse_var(Variant& var) {
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

        var.freq = 0.0;
        var.AC = 0;
        var.AD = {0, 0};
        var.quality = 0.0;

        if(!alleles.empty() && alleles[0].first != var.ref){
            var.alt = alleles[0].first;
            var.freq = static_cast<double>(variant_counts[var.alt]) / total_bases;
            var.AD = {variant_counts[var.ref], variant_counts[var.alt]};
            var.quality = variant_qualities[var.alt] / variant_counts[var.alt];

            // Only call homozygous variants (like the original code)
            if(var.freq >= threshold) {
                var.AC = 2;
            }
        }
    }
};

int main(int argc, char* argv[]) {
    if (argc != 4 && argc != 5) {
        std::cerr << "Usage: " << argv[0] << " <reference_genome> <sam_sorted_file> <vcf_output_file> [threshold]\n";
        return 1;
    }

    std::string reference_genome = argv[1];
    std::string sam_sorted_file = argv[2];
    std::string vcf_output_file = argv[3];
    double threshold = (argc == 5) ? std::stod(argv[4]) : 0.9; // Default to 0.9 for more stringent calling

    std::ifstream in(sam_sorted_file);
    if (!in) {
        std::cerr << "Could not open " << sam_sorted_file << "\n";
        return 1;
    }

    std::ofstream vcf_output(vcf_output_file);
    if (!vcf_output) {
        std::cerr << "Could not create " << vcf_output_file << "\n";
        return 1;
    }

    try {
        Converter converter;
        converter.process_sam(sam_sorted_file, reference_genome, threshold);
        converter.write_vcf(vcf_output_file);
        std::cout << "VCF file generated: " << vcf_output_file << "\n";
    } catch(const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}