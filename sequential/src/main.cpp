#include "arg_parser.hpp"
#include <iostream>

int main(int argc, char *argv[]) {
  try {
    utils::CmdArgs args = utils::ArgParser::parse(argc, argv);
    std::cout << "Wave simulation parameters loaded successfully." << std::endl;
  } catch (const utils::ParserException &e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }

  return 0;
}