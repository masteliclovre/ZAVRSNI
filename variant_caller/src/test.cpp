#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <stdexcept>

struct Variant {
   std::string chrom;
   int pos;          // 1-based position in reference
   std::string ref;  // Reference base(s)
   std::string alt;  // Alternate base(s)
   std::string type; // SNP, INS, DEL
   int count = 1;
};

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
   
public:
   void process_sam(const std::string& sam_file, const std::string& ref_fasta) {
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
           
         process_alignment(chrom, pos, cigar, seq);
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
       out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n";
       
       // Write variants
       for (const auto& [chrom, chrom_variants] : variants) {
           for (const Variant& var : chrom_variants) {
               out << chrom << "\t" << var.pos << "\t.\t" 
                   << var.ref << "\t" << var.alt << "\t.\tPASS\t"
                   << "DP=" << var.count << ";AC=" << var.count << ";AN=2\t"
                   << "GT:DP:AD\t1/1:" << var.count << ":" << var.count << "\n";
           }
       }
   }
   
private:
   void process_alignment(const std::string& chrom, int pos, const std::string& cigar, const std::string& seq) {
       std::vector<std::pair<int, char>> cigar_ops = parse_cigar(cigar);
       
       int ref_pos = pos;
       int seq_pos = 0;
       
       for (const auto& [len, op] : cigar_ops) {
           switch(op) {
               case 'M': case '=': case 'X': // Match or mismatch
                   for (int i = 0; i < len; ++i) {
                       std::string ref_base = ref_genome.get_sequence(chrom, ref_pos, 1);
                       std::string alt_base(1, seq[seq_pos]);
                       
                       if (ref_base != alt_base) { // SNP
                           add_variant(chrom, ref_pos, ref_base, alt_base, "SNP");
                       }
                       
                       ref_pos++;
                       seq_pos++;
                   }
                   break;
                   
               case 'I': // Insertion
               {
                   std::string ref_base = ref_genome.get_sequence(chrom, ref_pos, 1);
                   std::string inserted = seq.substr(seq_pos, len);
                   add_variant(chrom, ref_pos, ref_base, ref_base + inserted, "INS");
                   seq_pos += len;
                   break;
               }
                   
               case 'D': // Deletion
               {
                   std::string ref_bases = ref_genome.get_sequence(chrom, ref_pos, len);
                   add_variant(chrom, ref_pos, ref_bases, ref_bases.substr(0, 1), "DEL");
                   ref_pos += len;
                   break;
               }
                   
               case 'S': case 'H': // Soft or hard clip - skip in reference
                   seq_pos += len;
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
   
   void add_variant(const std::string& chrom, int pos, const std::string& ref, const std::string& alt, const std::string& type) {
       std::string key = chrom + ":" + std::to_string(pos) + ":" + ref + ">" + alt;
       
       // Check if variant already exists
       for (auto& var : variants[chrom]) {
           if (var.pos == pos && var.ref == ref && var.alt == alt) {
               var.count++;
               return;
           }
       }
       
       // Add new variant
       variants[chrom].push_back({chrom, pos, ref, alt, type});
   }
};

int main(int argc, char* argv[]) {

   if (argc != 4) {
      std::cerr << "Usage: " << argv[0] << " <reference_genome> <sam_sorted_file> <vcf_output_file>\n";
      return 1;
   }

   std::string reference_genome = argv[1];
   std::string sam_sorted_file = argv[2];
   std::string vcf_output_file = argv[3];

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
      converter.process_sam(sam_sorted_file, reference_genome);
      converter.write_vcf(vcf_output_file);
      std::cout << "VCF file generated: " << vcf_output_file << "\n";
   }catch(const std::exception& e) {
      std::cerr << "Error: " << e.what() << "\n";
      return 1;
   }

  return 0;
}
