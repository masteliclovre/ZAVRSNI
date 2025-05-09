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
    int pos;
    std::string ref;
    std::string alt;
    int coverage = 0;
    std::string bases;
    std::string qualities;
    double freq = 0.0;
    int AC = 0;
    int AN = 2;
    std::vector<int> AD = {0, 0};
    double quality = 0.0;
    int map_qual = 0;
    int count = 0;
};

/**
 * Converts a Phred quality score to a probability.
 * 
 * The Phred score is a logarithmic representation of the base-calling error
 * probabilities, defined as -10 * log10(probability of error).
 * 
 * @param phred An integer representing the Phred quality score.
 * @return The probability of the base call being incorrect.
 */
double phred_to_prob(int phred) {
    if (phred < 0) phred = 0;  // Ensure phred is not negative
    return std::pow(10.0, -static_cast<double>(phred) / 10.0);
}

/**
 * Converts a probability to a Phred quality score.
 *
 * If the probability is 0, the Phred score is 99.
 *
 * @param prob The probability of the base call being incorrect.
 * @return The Phred quality score.
 */
double prob_to_phred(double prob) {
    return (prob > 0) ? -10.0 * std::log10(prob) : 99.0;
}

class ReferenceGenome {
    std::unordered_map<std::string, std::string> chromosomes;
    
public:
    void load(const std::string& fasta_file) {
        //ucitaj dokument, ako ne uspijes, onda izbaci gresku
        std::ifstream in(fasta_file);
        if (!in) throw std::runtime_error("Cannot open reference file");
        
        std::string line, current_chrom, current_seq;
        
        while (std::getline(in, line)) {
            //prazne linije preskoci
            if (line.empty()) continue;
            
            //provjeri ako si na prvoj liniji fasta dokumenta, ako jesi onda dodaj njegovo ime u mapu the njegovu sekvecu
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
    void process_sam(const std::string& sam_file, const std::string& ref_fasta, double threshold) {
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
              
            const std::string& chrom = fields[2]; //reference name
            int pos = std::stoi(fields[3]); //1-based mapping position
            int map_qual = std::stoi(fields[4]); //mapping quality
            const std::string& cigar = fields[5];
            const std::string& seq = fields[9];
            const std::string& qual = fields[10];
              
            process_alignment(chrom, pos, map_qual, cigar, seq, qual);
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
        out << "##INFO=<ID=MQ,Number=1,Type=Float,Description=\"Average Mapping Quality\">\n";
        out << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
        out << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n";
        out << "##FORMAT=<ID=AD,Number=2,Type=Integer,Description=\"Allelic Depths\">\n";
        out << "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n";
        out << "##FORMAT=<ID=PL,Number=3,Type=Integer,Description=\"Phred-scaled Genotype Likelihoods\">\n";
        out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n";
        
        // Write variants
        for (const auto& [chrom, chrom_variants] : variants) {
            for (const Variant& var : chrom_variants) {
                // Only output variants with sufficient support and high enough frequency
                if (var.alt != var.ref && var.count >= min_support && var.freq >= threshold) {
                    std::string genotype = var.AC == 2 ? "1/1" : "0/1";
                    
                    // Calculate genotype likelihoods (PL)
                    double p_err = phred_to_prob(var.quality);
                    double p_hom_ref = (var.AC == 0) ? 1.0 - p_err : p_err;
                    double p_het = 0.5;
                    double p_hom_alt = (var.AC == 2) ? 1.0 - p_err : p_err;
                    
                    int pl_hom_ref = static_cast<int>(prob_to_phred(1.0 - p_hom_ref));
                    int pl_het = static_cast<int>(prob_to_phred(1.0 - p_het));
                    int pl_hom_alt = static_cast<int>(prob_to_phred(1.0 - p_hom_alt));
                    
                    // Genotype quality is the difference between best and next best genotype
                    int gq = std::min(pl_hom_ref - pl_het, pl_hom_alt - pl_het);
                    gq = std::max(0, std::min(gq, 99)); // Cap at 99
                    
                    out << chrom << "\t" << var.pos << "\t.\t" 
                        << var.ref << "\t" << var.alt << "\t" 
                        << static_cast<int>(var.quality) << "\tPASS\t"
                        << "DP=" << var.coverage << ";AC=" << var.AC << ";AN=2" 
                        << ";AF=" << var.freq << ";MQ=" << static_cast<int>(var.map_qual) << "\t"
                        << "GT:DP:AD:GQ:PL\t"
                        << genotype << ":" << var.coverage << ":" 
                        << var.AD[0] << "," << var.AD[1] << ":" 
                        << gq << ":"
                        << pl_hom_ref << "," << pl_het << "," << pl_hom_alt << "\n";
                }
            }
        }
    }
    
private:
void process_alignment(const std::string& chrom, int pos, int map_qual, const std::string& cigar, 
                       const std::string& seq, const std::string& qual) {
    std::vector<std::pair<int, char>> cigar_ops = parse_cigar(cigar);

    int ref_pos = pos;
    int seq_pos = 0;
    int qual_pos = 0;

    // Check if quality is missing (marked as *)
    bool has_quality = (qual != "*");

    for (const auto& [len, op] : cigar_ops) {
        // Skip if the sequence position is out of the sequence
        if (seq_pos >= seq.size()) break;

        switch (op) {
            case 'M': { // Match or mismatch
                for (int i = 0; i < len && seq_pos < seq.size(); i++) {
                    try {
                        std::string ref_base = ref_genome.get_sequence(chrom, ref_pos, 1);
                        std::string alt_base(1, seq[seq_pos]);

                        if (ref_base != alt_base) { // SNP
                            Variant var;
                            var.chrom = chrom;
                            var.pos = ref_pos;
                            var.ref = ref_base;
                            var.alt = alt_base;
                            var.bases = alt_base;
                            var.map_qual = map_qual;
                            var.qualities = (has_quality && qual_pos < qual.size()) ? 
                                            std::string(1, qual[qual_pos]) : "?";
                            var.coverage = 1;
                            var.count = 1;
                            add_variant(var);
                        }

                        ref_pos++;
                        seq_pos++;
                        if (has_quality && qual_pos < qual.size()) qual_pos++;
                    } catch (const std::exception& e) {
                        std::cerr << "Error processing M operation at " << chrom << ":" << ref_pos
                                  << ": " << e.what() << std::endl;
                        ref_pos++;
                        seq_pos++;
                        if (has_quality && qual_pos < qual.size()) qual_pos++;
                    }
                }
                break;
            }

            case 'I': { // Insertion
                if (seq_pos + len > seq.size()) {
                    std::cerr << "Warning: Insertion operation exceeds sequence length at "
                              << chrom << ":" << ref_pos << std::endl;
                    break;
                }

                try {
                    std::string ref_base = ref_genome.get_sequence(chrom, ref_pos, 1);
                    std::string inserted = seq.substr(seq_pos, len);
                    Variant var;
                    var.chrom = chrom;
                    var.pos = ref_pos;
                    var.ref = ref_base;
                    var.alt = ref_base + inserted;
                    var.bases = "+" + std::to_string(len) + inserted;
                    var.qualities = (has_quality && qual_pos + len <= qual.size()) ?
                                    qual.substr(qual_pos, len) : std::string(len, '?');
                    var.coverage = 1;
                    var.count = 1;
                    add_variant(var);
                    seq_pos += len;
                    if (has_quality && qual_pos + len <= qual.size()) qual_pos += len;
                } catch (const std::exception& e) {
                    std::cerr << "Error processing I operation at " << chrom << ":" << ref_pos
                              << ": " << e.what() << std::endl;
                    seq_pos += len;
                    if (has_quality && qual_pos + len <= qual.size()) qual_pos += len;
                }
                break;
            }

            case 'D': { // Deletion
                try {
                    std::string ref_bases = ref_genome.get_sequence(chrom, ref_pos, len);
                    Variant var;
                    var.chrom = chrom;
                    var.pos = ref_pos;
                    var.ref = ref_bases;
                    var.alt = ref_bases.substr(0, 1);
                    var.bases = "-" + std::to_string(len) + ref_bases.substr(1);
                    var.qualities = (has_quality && qual_pos < qual.size()) ?
                                    std::string(len, qual[qual_pos]) : std::string(len, '?');
                    var.coverage = 1;
                    var.count = 1;
                    add_variant(var);
                    ref_pos += len;
                    if (has_quality && qual_pos < qual.size()) qual_pos++;
                } catch (const std::exception& e) {
                    std::cerr << "Error processing D operation at " << chrom << ":" << ref_pos
                              << ": " << e.what() << std::endl;
                    ref_pos += len;
                    if (has_quality && qual_pos < qual.size()) qual_pos++;
                }
                break;
            }

            case 'S': case 'H': { // Soft or hard clip - skip in reference
                seq_pos += len;
                if (has_quality && qual_pos + len <= qual.size()) qual_pos += len;
                break;
            }

            case 'N': { // Skipped region in reference
                ref_pos += len;
                break;
            }

            default: {
                std::cerr << "Warning: Unknown CIGAR operation " << op << std::endl;
                break;
            }
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
        std::unordered_map<std::string, double> variant_qual_sums;
        std::unordered_map<std::string, double> variant_qual_prods;

        bool isDel = false;
        size_t indel_size = 0;
        std::string indel_sequence;
        size_t i = 0;
        size_t index = 0;
        size_t total_bases = 0;

        while(i < var.bases.size()) {
            char base = var.bases[i];
            
            // Handle quality scores safely
            char qual = (index < var.qualities.size()) ? var.qualities[index] : '?';
            int phred = (qual == '!') ? 0 :        // '!' is Phred 0
                    (qual == '?') ? 30 :        // Default Q30 when missing
                    (qual - 33 < 0) ? 0 :       // Ensure not negative
                    (qual - 33 > 93) ? 93 :     // Cap at Phred 93
                    qual - 33;                  // Standard Phred+33
            
            double prob_error = phred_to_prob(phred);
            double prob_correct = 1.0 - prob_error;

            // Skip quality markers and special characters
            if(base == '^') {  // Start of read segment
                i += 2;  // Skip the base following '^'
                continue;
            }
            if(base == '$') {  // End of read segment
                i++;
                continue;
            }
            if(base == '*' || base == '<' || base == '>') {  // Unmapped or reference skips
                i++;
                if (index < var.qualities.size()) index++;
                continue;
            }

            // Handle indels
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
                variant_qual_sums[key] += phred;
                variant_qual_prods[key] += std::log10(prob_correct);
                total_bases++;
                if (index < var.qualities.size()) index++;
                continue;
            }

            // Handle regular bases
            if(base == '.' || base == ',') {  // Reference matches
                variant_counts[var.ref]++;
                variant_qual_sums[var.ref] += phred;
                variant_qual_prods[var.ref] += std::log10(prob_correct);
                total_bases++;
                if (index < var.qualities.size()) index++;
            }
            else if(isalpha(base)) {  // Mismatches
                std::string alt(1, std::toupper(base));
                variant_counts[alt]++;
                variant_qual_sums[alt] += phred;
                variant_qual_prods[alt] += std::log10(prob_correct);
                total_bases++;
                if (index < var.qualities.size()) index++;
            }
            i++;
        }

        // Create alleles vector with combined metrics
        std::vector<std::pair<std::string, double>> alleles;
        for(const auto& entry : variant_counts) {
            const std::string& base = entry.first;
            int count = entry.second;
            if (count == 0) continue;

            double freq = static_cast<double>(count) / total_bases;
            double avg_qual = variant_qual_sums[base] / count;
            double log_prob = variant_qual_prods[base];
            
            // Combined score considering frequency, quality, and probability
            double score = freq * avg_qual * (1.0 - phred_to_prob(avg_qual));
            alleles.emplace_back(base, score);
        }

        // Sort alleles by score descending
        std::sort(alleles.begin(), alleles.end(),
            [](const std::pair<std::string, double>& a, const std::pair<std::string, double>& b) {
                return a.second > b.second;
            });

        // Reset variant metrics
        var.freq = 0.0;
        var.AC = 0;
        var.AD = {0, 0};
        var.quality = 0.0;

        if(!alleles.empty()) {
            // The reference allele
            var.AD[0] = variant_counts[var.ref];
            
            // Find the best non-reference allele
            for(const auto& allele : alleles) {
                if(allele.first != var.ref) {
                    var.alt = allele.first;
                    var.freq = static_cast<double>(variant_counts[var.alt]) / total_bases;
                    var.AD[1] = variant_counts[var.alt];
                    
                    // Calculate variant quality as the average of supporting read qualities
                    if (variant_counts[var.alt] > 0) {
                        var.quality = variant_qual_sums[var.alt] / var.AD[1];
                        
                        // Calculate genotype quality (GQ) - probability that this is not a sequencing error
                        double error_prob = phred_to_prob(var.quality);
                        var.quality = prob_to_phred(error_prob);
                    }
                    
                    if(var.freq >= threshold) {
                        var.AC = 2; // Homozygous variant
                    }
                    break;
                }
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
    //postavi threshold na 0.9, ako nije zadan na komandnoj liniji
    double threshold = (argc == 5) ? std::stod(argv[4]) : 0.9;

    //provjeri moze li se sam file otvoriti
    std::ifstream in(sam_sorted_file);
    if (!in) {
        std::cerr << "Could not open " << sam_sorted_file << "\n";
        return 1;
    }

    //provjeri moze li se u vcf file pisati
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