#include "spoa/spoa.hpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>

int main(int argc, char** argv) {
    const auto filePath = argv[5];
    std::vector<std::string> sequences;
    fprintf(stderr, "Reading query input file: %s\n", filePath);
    
    std::ifstream input(filePath);
    if(!input.good()){
      fprintf(stderr, "Bad input file: %s\n", filePath);
      exit(1);
    }

    std::string line, name, content;
    int64_t id = 0;
    input.seekg(0, input.beg);
    while (std::getline(input, line).good()) {
      if (line.empty() || line[0] == '>') { // Identifier marker
	if (!name.empty() ) { // Print out what we read from the last entry
	    sequences.push_back(content);
	    name.clear();
	}
	if (!line.empty()) {
	    name = line.substr(1);
	}
	content.clear();
      } else if (!name.empty()) {
	if (line.find(' ') != std::string::npos) { // Invalid sequence--no spaces allowed
	  name.clear();
	  content.clear();
	} else {
	  content += line;
	}
      } 
    }
    if (!name.empty()) {
      sequences.push_back(content);
    }  

    auto alignment_engine = spoa::createAlignmentEngine(static_cast<spoa::AlignmentType>(atoi(argv[1])),
        atoi(argv[2]), atoi(argv[3]), atoi(argv[4]));

    auto graph = spoa::createGraph();

    for (const auto& it: sequences) {
        auto alignment = alignment_engine->align_sequence_with_graph(it, graph);
        graph->add_alignment(alignment, it);
    }

    std::string consensus = graph->generate_consensus();

    fprintf(stderr, "Consensus (%zu)\n", consensus.size());
    fprintf(stderr, "%s\n", consensus.c_str());

    std::vector<std::string> msa;
    graph->generate_multiple_sequence_alignment(msa);

    fprintf(stderr, "Multiple sequence alignment\n");
    for (const auto& it: msa) {
        fprintf(stderr, "%s\n", it.c_str());
    }

    return 0;
}
