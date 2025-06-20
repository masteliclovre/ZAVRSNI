#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <string>
#include <cctype>
#include <algorithm>
#include <cmath>
#include <chrono>

std::string name = "UNKNOWN";
std::string length_str = "0";

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
    int AC;//broj varijanti
    int AN;//ukupan broj varijanti
    std::vector<int> AD;
    double quality;
};

struct PileupEntry {
    std::string chrom; //naziv kromosoma
    int pos; //pozicija na kojoj se nalazi u genomu (1-indeksiranje)
    char ref_base; //referentna baza na toj poziciji
    int coverage; //broj očitanja na toj poziciji
    std::string bases; //baze koje su pročitane zapisane u obliku . (poklapanje), slovo za SNP, + i - niz slova su indeli
    std::string qualities; // kvalitete pronađenih baza na tim pozicijama
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
    if (!alleles.empty() && alleles[0].first != var.ref) {
        std::string alt_key = alleles[0].first;
    
        if (alt_key[0] == '+') {
            // Insertion
            std::string inserted = alt_key.substr(1);  // izbrisi '+'
            var.alt = var.ref + inserted;
        } else if (alt_key[0] == '-') {
            // Deletion
            std::string deleted = alt_key.substr(1);  // izbrisi '-'
            var.alt = var.ref.substr(0, 1);
            var.ref = var.ref.substr(0, 1) + deleted;
        } else {
            //SNP
            var.alt = alt_key;
        }
    
        var.freq = static_cast<double>(variant_counts[alt_key]) / total_bases;
        var.AD = {variant_counts[var.ref.substr(0, 1)], variant_counts[alt_key]};
    
        if (var.freq >= 0.9) {
            var.AC = 2;
        } else if (var.freq >= threshold) {
            var.AC = 1;
        }
    
        var.quality = variant_qualities[alt_key] / variant_counts[alt_key];
    }
    

    //ako smo uspijeli nać varijaciju koja nam je korisna zapisujemo ju na standardni oblik u vcf datoteku
    if(var.alt != var.ref){
        vcf_output_file << var.chrom << "\t" << var.pos +1 << "\t.\t" << var.ref << "\t" << var.alt 
           << "\t.\t" << var.quality << "\tDP=" << var.coverage << ";AF=" << var.freq << "\n";
    }
}

class SamToPileupConverter {
    std::unordered_map<std::string, std::string> ref_sequences;
    
public:
    //funkcija za učitavanje fasta dokumenta u mapu
    void load_reference(const std::string& fasta_file) {
        //provjera, ako se fasta dokument može otvoriti za daljnje korištenje
        std::ifstream in(fasta_file);
        if (!in) throw std::runtime_error("Cannot open reference file");
        
        std::string line;
        std::string current_chrom;
        std::string current_seq;
        
        //dohvaća iz učitanog dokumenta liniju po liniju dok ne učita kraj dokumenta
        while (std::getline(in, line)) {
            //prazne linije preskačemo
            if (line.empty()) continue;
            
            //ako linija započinje '>', znamo da smo došli na novi kromosom
            if (line[0] == '>') {
                //ako varijabla current_chrom nije prazna, moramo ju mapirati s ucitanom sekvencom te clearati current_seq za novi kromosom
                if (!current_chrom.empty()) {
                    ref_sequences[current_chrom] = current_seq;
                    current_seq.clear();
                }
                //zapisivanje imena novom kromosomu
                current_chrom = line.substr(1);
                size_t space_pos = current_chrom.find(' ');
                if (space_pos != std::string::npos) {
                    current_chrom = current_chrom.substr(0, space_pos);
                }
            //liniju po liniju učitavamo sekvencu kromosoma
            } else {
                current_seq += line;
            }
        }
        //zadnje učitanje kromosoma mapiramo isto kao i prijašnje
        if (!current_chrom.empty()) {
            ref_sequences[current_chrom] = current_seq;
        }
    }
    
    //funkcija dohvaca nukleotid na navedenoj pozicji inace ispisuje error poruku
    char get_reference_base(const std::string& chrom, int pos) {
        //provjeravamo ako imamo kromosom pod navedenim imenom 
        if (ref_sequences.find(chrom) == ref_sequences.end()) {
            throw std::runtime_error("Chromosome not found: " + chrom);
        }
        
        const std::string& seq = ref_sequences.at(chrom);
        
        //gledamo ako je predana pozicija ispravno napisana te ako se nalazi unutar opsega sekvence
        if (pos < 1 || pos > seq.length()) {
            throw std::runtime_error("Position out of bounds: " + std::to_string(pos));
        }
        
        //ako postoji takva pozicija u sekvenci, vraćamo nazad nukleotid na toj poziciji
        return seq[pos - 1];
    }
    
    //funkcija učitava liniju po liniju poravnatih podataka parsira ju i salje na obradu očitanja
    void process_sam(const std::string& sam_file, const std::string& output_file) {
        //provjera ako je moguće otvoriti navedeni sam file
        std::ifstream in(sam_file);
        if (!in) throw std::runtime_error("Cannot open SAM file");
        
        //mapa za pohranu pileup podataka
        std::unordered_map<std::string, std::unordered_map<int, PileupEntry>> pileup_data;
         
        std::string line;
        while (std::getline(in, line)) {
            //preskačemo linije koje pocinju '@', jer su to komentari
            if(line.substr(0, 3) == "@SQ"){
                name = line.substr(7, line.find("LN:") - 8);
                length_str = line.substr(line.find("LN:") + 3);
            }
            if (line[0] == '@') continue;
            
            std::istringstream ss(line);
            std::vector<std::string> fields;
            std::string field;
            
            //rastavljamo liniju po tabovima i zapisujemo u vektor
            while (std::getline(ss, field, '\t')) {
                fields.push_back(field);
            }
            
            //polja mora biti minimalno 11, inače nije ispravno formatirana linija i možemo ju preskočiti
            if (fields.size() < 11) continue;
              
            const std::string& chrom = fields[2]; //kromosom
            int pos = std::stoi(fields[3]); //pozicija
            const std::string& cigar = fields[5]; //niz operacija match, delete, insert (M,D,I,S) s brojem ponavljanja upisanim ispred znaka
            const std::string& seq = fields[9]; //sekvenca
            const std::string& qual = fields[10]; //kvalitete
            
            process_read(chrom, pos, cigar, seq, qual, pileup_data);
        }
        
        write_pileup(output_file, pileup_data);
    }
    
private:
    void process_read(const std::string& chrom, int pos, const std::string& cigar, const std::string& seq, const std::string& qual, std::unordered_map<std::string, std::unordered_map<int, PileupEntry>>& pileup_data) {
        //pretvorba cigar stringa u vektor sastavljen od operacija i broj njihovih ponavljanja
        std::vector<std::pair<int, char>> cigar_ops = parse_cigar(cigar);
        int ref_pos = pos; //pozicija iz linije sam filea da možemo poravnti sekvence
        int seq_pos = 0;
        int qual_pos = 0;
        bool has_quality = (qual != "*");

        for (const auto& [len, op] : cigar_ops) {
            //u slučaju da smo došli do kraja sekvence kromosoma, prekidamo obradu
            if (seq_pos >= seq.size()) break;

            switch (op) {
                case 'M': { // Match/mismatch
                    //Obradujemo sekvencu uzimajuci u obzir broj ponavljanja slova M koji predstavlja match/mismatch te maksimalnu duljinu sekvence koju ne smijemo prekoračiti
                    for (int i = 0; i < len && seq_pos < seq.size(); i++) {
                        try {
                            char ref_base = get_reference_base(chrom, ref_pos); //nukleotid na referentnoj poziciji kormosoma s kojim uspoređujemo
                            char base = seq[seq_pos]; //trenutan nukleotid na sekvenci koja se obrađuje
                            char quality = (has_quality && qual_pos < qual.size()) ? qual[qual_pos] : '~'; //provjera kvalitete
                            
                            //Pretrazivanje zapisa u pileup mapi uz ime i referentnu poziciju
                            auto& pos_entry = pileup_data[chrom][ref_pos];
                            //ako već nije dodan zapis za ovaj kromosom i referentnu poziciju, dodajemo novi
                            if (pos_entry.chrom.empty()) {
                                pos_entry.chrom = chrom; //kromosom
                                pos_entry.pos = ref_pos; //pozicija na referentnom kromosom
                                pos_entry.ref_base = ref_base; //referentni nukleotid
                            }
                            
                            //Usporedba referentnog i trenutnog nukleotida
                            //ako su isti dodajemo znak '.'
                            if (base == ref_base) {
                                pos_entry.bases += '.';
                            //inače dodajemo nukleotid koji ga je zamjenio
                            } else {
                                pos_entry.bases += base;
                            }
                            //zapisujemo isto tako i kvalitetu nukleotida, ako postoji
                            pos_entry.qualities += quality;
                            //inkrementiramo broj očitanja za svaki zapis
                            pos_entry.coverage++;
                            
                            //inkrementiramo pozicije za referetnu sekvencu i sekvencu koju obrađujemo
                            ref_pos++;
                            seq_pos++;
                            //te ako je postoji zapis kvalitete, onda i nju inkrementiramo
                            if (has_quality && qual_pos < qual.size()) qual_pos++;
                        } catch (const std::exception& e) {
                            //U slučaju greške preskačemo poziciju i inkrementiramo referetnu poziciju i poziciju trenutne sekvence
                            ref_pos++;
                            seq_pos++;
                            if (has_quality && qual_pos < qual.size()) qual_pos++;
                        }
                    }
                    break;
                }
                
                case 'I': { // Insertion
                    //ako je pozicija u sekvenci zajedno s duljinom insertion-a veća od maksimalne duljine sekvence, prekidamo obradu
                    if (seq_pos + len > seq.size()) break;
                    
                    try {
                        //tražimo nukleotid na referentnoj poziciji kormosoma s kojim uspoređujemo
                        char ref_base = get_reference_base(chrom, ref_pos - 1);
                        //uzimamo iz trenutne sekvence podniz koji predstavlja insertion
                        std::string inserted = seq.substr(seq_pos, len);
                        
                        //Pretrazivanje zapisa u pileup mapi uz ime i referentnu poziciju koja je za jedan u nazad od referentne pozicije jer insertion nije nova pozicija već dodatak na staru
                        auto& pos_entry = pileup_data[chrom][ref_pos - 1];
                        if (pos_entry.chrom.empty()) {
                            pos_entry.chrom = chrom;
                            pos_entry.pos = ref_pos - 1;
                            pos_entry.ref_base = ref_base;
                        }
                        
                        //zapisujemo bazu s plus koji predstavlja insertion i podniz trenutne sekvence
                        pos_entry.bases += "+" + std::to_string(len) + inserted;
                        if (has_quality && qual_pos + len <= qual.size()) {
                            pos_entry.qualities += qual.substr(qual_pos, len);
                        } else {
                            pos_entry.qualities += std::string(len, '~');
                        }
                        //inkrementiramo broj očitanja na referentnoj poziciji
                        pos_entry.coverage++;
                        //inkrementiramo poziciju na trenutnoj sekvenci za duljinu insertion-a
                        seq_pos += len;
                        if (has_quality && qual_pos + len <= qual.size()) qual_pos += len;
                    } catch (...) {
                        //U slučaju greške preskačemo poziciju i inkrementiramo referetnu poziciju i poziciju trenutne sekvence za duljinu insertion-a
                        seq_pos += len;
                        if (has_quality && qual_pos + len <= qual.size()) qual_pos += len;
                    }
                    break;
                }
                
                case 'D': { // Deletion
                    try {
                        std::string deleted_bases;
                        //dohvaćamo iz kromosoma sve baze koje su izbrisane
                        for (int i = 0; i < len; i++) {
                            deleted_bases += get_reference_base(chrom, ref_pos + i);
                        }
                        
                        //Pretrazivanje zapisa u pileup mapi uz ime i referentnu poziciju
                        auto& pos_entry = pileup_data[chrom][ref_pos];
                        if (pos_entry.chrom.empty()) {
                            pos_entry.chrom = chrom;
                            pos_entry.pos = ref_pos;
                            pos_entry.ref_base = get_reference_base(chrom, ref_pos);
                        }
                        //zapisujemo bazu s minus koji predstavlja deletion popraćenu podnizom koji je izbrisan
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

                //predstavlja dio koji pokazuje što preskačemo na trenutnoj sekvenci, a odnosi se na početak i kraj sekvence koju smo poravnali, može biti dulji i kraći (hard and soft)
                case 'S': case 'H': { // Soft or hard clip
                    seq_pos += len;
                    if (has_quality && qual_pos + len <= qual.size()) qual_pos += len;
                    break;
                }
                
                //mogući dijelovu unutar referentne sekvence koje preskaćemo
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
    
    //uzimamo podatke koje smo prikupili u mapu pileup-a i zapisujemo ih u datoteku
    void write_pileup(const std::string& output_file, const std::unordered_map<std::string, std::unordered_map<int, PileupEntry>>& pileup_data) {
        //provjera ako je moguće stvoriti/otvoriti datoteku 
        std::ofstream out(output_file);
        if (!out) throw std::runtime_error("Cannot create pileup file");
        
        //Uzimamo kromosome iz mape te ih sortiramo po abecedi da održimo redoslijed u ispisu
        std::vector<std::string> chromosomes;
        for (const auto& chrom_entry : pileup_data) {
            chromosomes.push_back(chrom_entry.first);
        }
        std::sort(chromosomes.begin(), chromosomes.end());
        
        for (const auto& chrom : chromosomes) {
            //Uzimamo pozicije iz mape i sortiramo i njih da održimo redoslijed u ispisu, zato što mape nisu sortirane 
            std::vector<int> positions;
            for (const auto& pos_entry : pileup_data.at(chrom)) {
                positions.push_back(pos_entry.first);
            }
            std::sort(positions.begin(), positions.end());
            
            //Zapisujemo podatke u datoteku
            for (int pos : positions) {
                const auto& entry = pileup_data.at(chrom).at(pos);
                out << entry.chrom << "\t" << entry.pos << "\t" << entry.ref_base << "\t"
                    << entry.coverage << "\t" << entry.bases << "\t" << entry.qualities << "\n";
            }
        }
    }
};

int main(int argc, char* argv[]) {

    auto start = std::chrono::high_resolution_clock::now();

     //provjera broja argumenata i ispis poruke u slučaju neispravnosti
     if (argc != 5) {
        std::cerr << "Usage: " << argv[0] << " <reference_genome> <sam_file> <vcf_output_file> <threshold>\n";
        return 1;
    }

    std::string reference_file = argv[1];
    std::string sam_file = argv[2];
    std::string vcf_output_file = argv[3];
    double threshold = std::stod(argv[4]);
    std::string pileup_output = "mympileup.mpileup";

    try {
        SamToPileupConverter converter;
        converter.load_reference(reference_file);
        converter.process_sam(sam_file, pileup_output);
        //poruka koja nam pokazuje da je pileup datoteka uspjesno napravljena
        std::cout << "Pileup file generated: " << pileup_output << "\n";
    } catch(const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    std::string pileup_file = pileup_output;

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

    // VCF header
    vcf_output << "##fileformat=VCFv4.2\n";
    vcf_output << "##FILTER=<ID=PASS,Description=\"All filters passed\">\n";
    vcf_output << "##contig=<ID="<< name << ",length="<< length_str <<">\n";
    vcf_output << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n";
    vcf_output << "##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency\">\n";
    vcf_output << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
    vcf_output << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

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

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Time taken: " << elapsed.count() << " s\n";

    return 0;
}