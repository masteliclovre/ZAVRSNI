#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

int main(int argc, char* argv[]) {

    if(argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <pileup_file>\n";
        return 1;
    }

    std::ifstream in(argv[1]);
    if(!in) {
        std::cerr << "Could not open " << argv[1] << "\n";
        return 1;
    }

    std::string line;
    while(std::getline(in, line)) {
        std::cout << line << "\n";
    }

    return 0;
}