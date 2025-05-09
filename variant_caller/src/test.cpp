#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <unordered_map>
#include <algorithm>
#include <cctype>

struct PileupEntry {
    std::string chrom;
    int pos;
    char ref_base;
    int coverage;
    std::string bases;
    std::string qualities;
};

class SamToPileupConverter {
    std::unordered_map<std::string, std::string> ref_sequences;
    
public:
    void load_reference(const std::string& fasta_file) {
        std::ifstream in(fasta_file);
        if (!in) throw std::runtime_error("Cannot open reference file");
        
        std::string line, current_chrom, current_seq;
        
        while (std::getline(in, line)) {
            if (line.empty()) continue;
            
            if (line[0] == '>') {
                if (!current_chrom.empty()) {
                    ref_sequences[current_chrom] = current_seq;
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
            ref_sequences[current_chrom] = current_seq;
        }
    }
    
    char get_reference_base(const std::string& chrom, int pos) {
        if (ref_sequences.find(chrom) == ref_sequences.end()) {
            throw std::runtime_error("Chromosome not found: " + chrom);
        }
        
        const std::string& seq = ref_sequences.at(chrom);
        if (pos < 1 || pos > seq.length()) {
            throw std::runtime_error("Position out of bounds: " + std::to_string(pos));
        }
        
        return seq[pos - 1];
    }
    
    void process_sam(const std::string& sam_file, const std::string& output_file) {
        std::ifstream in(sam_file);
        if (!in) throw std::runtime_error("Cannot open SAM file");
        
        std::unordered_map<std::string, std::unordered_map<int, PileupEntry>> pileup_data;
         
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
            int pos = std::stoi(fields[3]);
            const std::string& cigar = fields[5];
            const std::string& seq = fields[9];
            const std::string& qual = fields[10];
            
            process_read(chrom, pos, cigar, seq, qual, pileup_data);
        }
        
        write_pileup(output_file, pileup_data);
    }
    
private:
    void process_read(const std::string& chrom, int pos, const std::string& cigar, 
                 const std::string& seq, const std::string& qual,
                 std::unordered_map<std::string, std::unordered_map<int, PileupEntry>>& pileup_data) {
        std::vector<std::pair<int, char>> cigar_ops = parse_cigar(cigar);
        int ref_pos = pos;
        int seq_pos = 0;
        int qual_pos = 0;
        bool has_quality = (qual != "*");

        for (const auto& [len, op] : cigar_ops) {
            if (seq_pos >= seq.size()) break;

            switch (op) {
                case 'M': case '=': case 'X': { // Match/mismatch
                    for (int i = 0; i < len && seq_pos < seq.size(); i++) {
                        try {
                            char ref_base = get_reference_base(chrom, ref_pos);
                            char base = seq[seq_pos];
                            char quality = (has_quality && qual_pos < qual.size()) ? qual[qual_pos] : '~';
                            
                            // Add to pileup
                            auto& pos_entry = pileup_data[chrom][ref_pos];
                            if (pos_entry.chrom.empty()) {
                                pos_entry.chrom = chrom;
                                pos_entry.pos = ref_pos;
                                pos_entry.ref_base = ref_base;
                            }
                            
                            // Determine base representation
                            if (base == ref_base) {
                                pos_entry.bases += '.';
                            } else {
                                pos_entry.bases += base;
                            }
                            pos_entry.qualities += quality;
                            pos_entry.coverage++;
                            
                            ref_pos++;
                            seq_pos++;
                            if (has_quality && qual_pos < qual.size()) qual_pos++;
                        } catch (const std::exception& e) {
                            // Skip if reference position is invalid
                            ref_pos++;
                            seq_pos++;
                            if (has_quality && qual_pos < qual.size()) qual_pos++;
                        }
                    }
                    break;
                }
                
                case 'I': { // Insertion
                    if (seq_pos + len > seq.size()) break;
                    
                    try {
                        char ref_base = get_reference_base(chrom, ref_pos - 1); // Base before insertion
                        std::string inserted = seq.substr(seq_pos, len);
                        
                        // Add to pileup at the previous position
                        auto& pos_entry = pileup_data[chrom][ref_pos - 1];
                        if (pos_entry.chrom.empty()) {
                            pos_entry.chrom = chrom;
                            pos_entry.pos = ref_pos - 1;
                            pos_entry.ref_base = ref_base;
                        }
                        
                        pos_entry.bases += "+" + std::to_string(len) + inserted;
                        if (has_quality && qual_pos + len <= qual.size()) {
                            pos_entry.qualities += qual.substr(qual_pos, len);
                        } else {
                            pos_entry.qualities += std::string(len, '~');
                        }
                        pos_entry.coverage++;
                        
                        seq_pos += len;
                        if (has_quality && qual_pos + len <= qual.size()) qual_pos += len;
                    } catch (...) {
                        seq_pos += len;
                        if (has_quality && qual_pos + len <= qual.size()) qual_pos += len;
                    }
                    break;
                }
                
                case 'D': { // Deletion
                    try {
                        std::string deleted_bases;
                        for (int i = 0; i < len; i++) {
                            deleted_bases += get_reference_base(chrom, ref_pos + i);
                        }
                        
                        // Add to pileup at the current position
                        auto& pos_entry = pileup_data[chrom][ref_pos];
                        if (pos_entry.chrom.empty()) {
                            pos_entry.chrom = chrom;
                            pos_entry.pos = ref_pos;
                            pos_entry.ref_base = get_reference_base(chrom, ref_pos);
                        }
                        
                        pos_entry.bases += "-" + std::to_string(len) + deleted_bases;
                        pos_entry.qualities += (has_quality && qual_pos < qual.size()) ? qual[qual_pos] : '~';
                        pos_entry.coverage++;
                        
                        ref_pos += len;
                        if (has_quality && qual_pos < qual.size()) qual_pos++;
                    } catch (...) {
                        ref_pos += len;
                        if (has_quality && qual_pos < qual.size()) qual_pos++;
                    }
                    break;
                }
                
                case 'S': case 'H': { // Soft or hard clip
                    seq_pos += len;
                    if (has_quality && qual_pos + len <= qual.size()) qual_pos += len;
                    break;
                }
                
                case 'N': { // Skipped region
                    ref_pos += len;
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
    
    void write_pileup(const std::string& output_file, 
                     const std::unordered_map<std::string, std::unordered_map<int, PileupEntry>>& pileup_data) {
        std::ofstream out(output_file);
        if (!out) throw std::runtime_error("Cannot create pileup file");
        
        // Collect all chromosomes and positions to sort them
        std::vector<std::string> chromosomes;
        for (const auto& chrom_entry : pileup_data) {
            chromosomes.push_back(chrom_entry.first);
        }
        std::sort(chromosomes.begin(), chromosomes.end());
        
        for (const auto& chrom : chromosomes) {
            // Get all positions for this chromosome and sort them
            std::vector<int> positions;
            for (const auto& pos_entry : pileup_data.at(chrom)) {
                positions.push_back(pos_entry.first);
            }
            std::sort(positions.begin(), positions.end());
            
            // Write entries for this chromosome
            for (int pos : positions) {
                const auto& entry = pileup_data.at(chrom).at(pos);
                out << entry.chrom << "\t" << entry.pos << "\t" << entry.ref_base << "\t"
                    << entry.coverage << "\t" << entry.bases << "\t" << entry.qualities << "\n";
            }
        }
    }
};

int main(int argc, char* argv[]) {
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <reference_genome> <sam_file> <pileup_output>\n";
        return 1;
    }

    std::string reference_file = argv[1];
    std::string sam_file = argv[2];
    std::string pileup_output = argv[3];

    try {
        SamToPileupConverter converter;
        converter.load_reference(reference_file);
        converter.process_sam(sam_file, pileup_output);
        std::cout << "Pileup file generated: " << pileup_output << "\n";
    } catch(const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}