#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <unordered_map>
#include <algorithm>
#include <cctype>

struct PileupEntry {
    std::string chrom; //naziv kromosoma
    int pos; //pozicija na kojoj se nalazi u genomu (1-indeksiranje)
    char ref_base; //referentna baza na toj poziciji
    int coverage; //broj očitanja na toj poziciji
    std::string bases; //baze koje su pročitane zapisane u obliku . (poklapanje), slovo za SNP, + i - niz slova su indeli
    std::string qualities; // kvalitete pronađenih baza na tim pozicijama
};

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
    
    //funkcija učitava linijju po liniju poravnatih podataka parsira ju i salje na obradu očitanja
    void process_sam(const std::string& sam_file, const std::string& output_file) {
        //provjera ako je moguće otvoriti navedeni sam file
        std::ifstream in(sam_file);
        if (!in) throw std::runtime_error("Cannot open SAM file");
        
        //mapa za pohranu pileup podataka
        std::unordered_map<std::string, std::unordered_map<int, PileupEntry>> pileup_data;
         
        std::string line;
        while (std::getline(in, line)) {
            //preskačemo linije koje pocinju '@', jer su to komentari
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
    //provjera broja argumenata i ispis poruke u slučaju neispravnosti
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
        //poruka koja nam pokazuje da je pileup datoteka uspjesno napravljena
        std::cout << "Pileup file generated: " << pileup_output << "\n";
    } catch(const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}