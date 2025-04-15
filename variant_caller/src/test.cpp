#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>

// Struktura za varijantu
struct Variant {
    std::string type;  // Tip varijante: SNP, I (insertion), D (deletion)
    int pos;           // Pozicija u referenci
    std::string alt;   // Alternativna baza ili "-"
};

// Funkcija za parsiranje CIGAR stringa
std::vector<std::string> parseCigar(const std::string& cigar) {
    std::vector<std::string> operations;
    int length = 0;
    for (char c : cigar) {
        if (isdigit(c)) {
            length = length * 10 + (c - '0');
        } else {
            operations.push_back(std::to_string(length) + c);
            length = 0;
        }
    }
    return operations;
}

// Funkcija za analizu SAM linije i pozivanje varijanti
void processSamLine(const std::string& samLine, const std::string& referenceGenome) {
    std::istringstream ss(samLine);
    std::string field;
    
    // Spliting line by tabs
    std::vector<std::string> fields;
    while (std::getline(ss, field, '\t')) {
        fields.push_back(field);
    }

    int pos = std::stoi(fields[3]);  // Pozicija u referenci (0-based)
    std::string cigar = fields[5];   // CIGAR string
    std::string sequence = fields[9]; // Read sequence

    std::vector<std::string> operations = parseCigar(cigar);

    // Daljnja logika za određivanje varijanti na temelju CIGAR stringa
    // (npr. umetanje, delecija, SNP)
}

// Funkcija za generiranje VCF formata
void generateVCF(const std::vector<Variant>& variants) {
    std::ofstream vcf_output("output.vcf");
    vcf_output << "##fileformat=VCFv4.2\n";
    vcf_output << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

    for (const auto& var : variants) {
        vcf_output << "lambda\t" << var.pos + 1 << "\t.\t" << var.alt << "\n";
    }

    vcfFile.close();
}

int main() {
    // Učitavanje sam datoteke i procesiranje linija
    std::ifstream samFile("sample_sorted.sam");
    std::string line;
    std::vector<Variant> variants;

    while (std::getline(samFile, line)) {
        // Ako linija nije komentar
        if (line[0] != '@') {
            processSamLine(line, "reference_genome.fasta");
        }
    }

    // Generiranje VCF rezultata
    generateVCF(variants);

    return 0;
}
