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

//struktura za spremanje podataka o varijanti
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

//pretvara phred rezultat u vjerojatnost
double phred_to_prob(int phred) {
    return std::pow(10.0, -phred / 10.0);
}


 //pretvara vjerojatnost u phred rezultat
double prob_to_phred(double prob) {
    return (prob > 0) ? -10.0 * std::log10(prob) : 99.0;
}

void parse_var(Variant& var, std::ofstream& vcf_output_file, double threshold) {
    //mapa za spremanje određene varijante i broj njenih ponavljanja
    std::unordered_map<std::string, int> variant_counts;
    //mapa za spremanje određene varijante i njihovih kvaliteta
    std::unordered_map<std::string, double> variant_qualities;

    bool isDel = false;
    size_t indel_size = 0;
    std::string indel_sequence;
    size_t i = 0;
    size_t index = 0;
    size_t total_bases = 0;

    //obradujemo bazu po bazu varijanti na odredenoj referentoj poziciji
    while(i < var.bases.size()){
        char base = var.bases[i]; //spremamo bazu
        char qual = (index < var.qualities.size()) ? var.qualities[index] : '!'; //ako postoji kvaliteta i nju spremamo
        int phred = (qual != '!') ? qual - 33 : 30; //ako nema kvalitete, postavljamo defaultnu vrijednost

        //ignoriramo znakove ^, $, *, <, >, jer ne sadrže informacije relevantne za našu analizu
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

        //ako je trenutna baza + ili - znači da je to indel i moramo ga obraditi
        if(base == '+' || base == '-') {
            isDel = (base == '-');
            i++;

            //zapisujemo duljinu indela tako što provjeravamo je li idući znak broj, ako nije došli smo do zadnje znamenke duljine
            indel_size = 0;
            while(i < var.bases.size() && isdigit(var.bases[i])) {
                indel_size = indel_size * 10 + (var.bases[i] - '0');
                //inkrementiramo poziciju u bazama nakon svake pročitane znamenke
                i++;
            }

            //ispraznimo prijašnji niz zapisan u stringu te sad zapisujemo novu sekvencu, duljine koju smo pronašli korak prije
            indel_sequence.clear();
            for(size_t j = 0; j < indel_size && i < var.bases.size(); j++, i++) {
                indel_sequence.push_back(std::toupper(var.bases[i]));
            }

            //zapisujemo varijantu u mapu i inkrementiramo ukupni broj baza i index za kvalitete
            std::string key = (isDel ? "-" : "+") + indel_sequence;
            variant_counts[key]++;
            variant_qualities[key] += phred;
            total_bases++;
            index++;
            continue;
        }

        //baza je jednaka referentoj bazi 
        if(base == '.' || base == ',') {
            variant_counts[var.ref]++;
            variant_qualities[var.ref] += phred;
            total_bases++;
            index++;
        }
        //baza je SNP 
        else if(isalpha(base)) {
            std::string alt(1, std::toupper(base));
            variant_counts[alt]++;
            variant_qualities[alt] += phred;
            total_bases++;
            index++;
        }
        i++;
    }

    //izracunavamo frekvenciju pojavljivanja varijacije kao omjer pojavljivanja s ukupnim brojem variranih baza
    std::vector<std::pair<std::string, double>> alleles;
    for(const auto& entry : variant_counts) {
        const std::string& base = entry.first;
        int count = entry.second;
        double freq = static_cast<double>(count) / total_bases;
        //kvaliteta varijacije je zbroj kvaliteta svih ponavljanja varijacije podjeljen s ukupnim brojem ponavljanja
        double avg_qual = variant_qualities[base] / count;
        alleles.emplace_back(base, freq * avg_qual);
    }

    //sortiramo novonastalu mapu varijacija po rezultatu freq*avg_qual od veceg prema manjem
    std::sort(alleles.begin(), alleles.end(), [](const std::pair<std::string, double>& a, const std::pair<std::string, double>& b) {
            return a.second > b.second;
        });

    //postavljamo defaultne vrijednosti
    var.alt = var.ref;
    var.freq = 0.0;
    var.AC = 0;
    var.AN = 2;
    var.AD = {0, 0};
    var.quality = 0.0;

    //uzimamo najbolju varijaciju ako postoji i razlicita je od referente, onda spremamo njene vrijednosti
    if(!alleles.empty() && alleles[0].first != var.ref){
        var.alt = alleles[0].first;
        var.freq = static_cast<double>(variant_counts[var.alt]) / total_bases;
        //broj pojavljivanja referentnog nukleotida i alternativnog
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

    //ako smo uspijeli nać varijaciju koja nam je korisna zapisujemo ju na standardni oblik u vcf datoteku
    if(var.alt != var.ref){
        vcf_output_file << var.chrom << "\t" << var.pos << "\t" << '.' << "\t" << var.ref << "\t" << var.alt << "\t" 
                        << var.quality << "\t" << "PASS" << "\t" 
                        << "DP=" << var.coverage << ";AC=" << var.AC << ";AN=" << var.AN << ";AF=" << var.freq << "\t" 
                        << "GT:DP:AD:PL" << "\t";
        
        double p_ref = phred_to_prob(var.quality); //vjerojatnost da je ref ispravan
        double p_het = 0.5; //pretpostavljamo vjerojatnost da je varijacija heterozigotna
        double p_hom = 1.0 - p_het; //pretpostavljamo vjerojatnost da je varijacija homozigotna

        //phred-scaled vjerojatnosti genotipova
        int pl_ref = static_cast<int>(prob_to_phred(1.0 - p_ref));
        int pl_het = static_cast<int>(prob_to_phred(1.0 - p_het));
        int pl_hom = static_cast<int>(prob_to_phred(1.0 - p_hom));

        vcf_output_file << genotype << ":" << var.coverage << ":" 
                        << var.AD[0] << "," << var.AD[1] << ":" 
                        << pl_ref << "," << pl_het << "," << pl_hom << "\n";
    }
}

int main(int argc, char* argv[]) {
    //provjera broja argumenata i ispis poruke u slučaju neispravnosti
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <pileup_file> <vcf_output_file> <threshold>\n";
        return 1;
    }

    std::string pileup_file = argv[1];
    std::string vcf_output_file = argv[2];
    double threshold = std::stod(argv[3]);

    //ako se datoteke ne mogu otvoriti/stvoriti ispisuj poruku o grešci i prekini program
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

    // VCF headeri
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
    //liniju po liniju pileup datoteke parsiramo
    while (std::getline(in, line)) {
        std::stringstream ss(line);
        std::string chrom;
        int pos;
        std::string ref;
        int coverage;
        std::string bases;
        std::string qualities;

        ss >> chrom >> pos >> ref >> coverage >> bases >> qualities;
        pos -= 1;
        //ako je na poziciji referentnog nukleotida postoji varijacija zapisujemo ju u variantu i parsiramo
        if (coverage > 0) {
            Variant var = {chrom, pos, ref, ref, coverage, bases, qualities, 0.0, 0, 0};
            //pazimo da su sva slova velika za lakše procesiranje
            std::transform(var.ref.begin(), var.ref.end(), var.ref.begin(), [](unsigned char c) { return std::toupper(c); });
            parse_var(var, vcf_output, threshold);
        }
    }

    std::cout << "VCF file generated: " << vcf_output_file << "\n";
    return 0;
}