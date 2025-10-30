#include "arg_parser.hpp"
#include <cstring>
#include <gtest/gtest.h>
#include <vector>

char **makeArgv(const std::vector<std::string> &args) {
  char **argv = new char *[args.size()];
  for (size_t i = 0; i < args.size(); ++i) {
    argv[i] = new char[args[i].size() + 1];
    std::strcpy(argv[i], args[i].c_str());
  }
  return argv;
}

void freeArgv(char **argv, int argc) {
  for (int i = 0; i < argc; ++i) {
    delete[] argv[i];
  }
  delete[] argv;
}

TEST(ArgParserTest, ValidInputWithTau) {
  std::vector<std::string> args = {
      "program", "--Lx", "1.0",  "--Ly", "2.0", "--Lz", "3.0",   "--Nx", "10",
      "--Ny",    "20",   "--Nz", "30",   "--T", "10.0", "--tau", "0.1"};
  char **argv = makeArgv(args);

  utils::CmdArgs params = utils::ArgParser::parse(args.size(), argv);

  EXPECT_DOUBLE_EQ(params.Lx, 1.0);
  EXPECT_DOUBLE_EQ(params.Ly, 2.0);
  EXPECT_DOUBLE_EQ(params.Lz, 3.0);
  EXPECT_EQ(params.Nx, 10);
  EXPECT_EQ(params.Ny, 20);
  EXPECT_EQ(params.Nz, 30);
  EXPECT_DOUBLE_EQ(params.T, 10.0);
  EXPECT_DOUBLE_EQ(params.tau, 0.1);
  EXPECT_EQ(params.K, 100);

  freeArgv(argv, args.size());
}

TEST(ArgParserTest, ValidInputWithK) {
  std::vector<std::string> args = {
      "program", "--Lx", "1.0",  "--Ly", "1.0", "--Lz", "1.0", "--Nx", "100",
      "--Ny",    "100",  "--Nz", "100",  "--T", "5.0",  "--K", "50"};
  char **argv = makeArgv(args);

  utils::CmdArgs params = utils::ArgParser::parse(args.size(), argv);

  EXPECT_EQ(params.K, 50);
  EXPECT_DOUBLE_EQ(params.tau, 0.1);

  freeArgv(argv, args.size());
}

TEST(ArgParserTest, ComputeKWithNonDivisibleTau) {
  std::vector<std::string> args = {
      "program", "--Lx", "1.0",  "--Ly", "1.0", "--Lz", "1.0",   "--Nx", "100",
      "--Ny",    "100",  "--Nz", "100",  "--T", "10.0", "--tau", "0.3"};
  char **argv = makeArgv(args);

  utils::CmdArgs params = utils::ArgParser::parse(args.size(), argv);

  EXPECT_DOUBLE_EQ(params.tau, 0.3);
  EXPECT_EQ(params.K, 34);

  freeArgv(argv, args.size());
}

TEST(ArgParserTest, MissingRequiredParam) {
  std::vector<std::string> args = {"program", "--Lx", "1.0", "--Ly", "1.0"};
  char **argv = makeArgv(args);

  EXPECT_THROW(
      { utils::ArgParser::parse(args.size(), argv); }, utils::ParserException);

  freeArgv(argv, args.size());
}

TEST(ArgParserTest, NegativeLx) {
  std::vector<std::string> args = {
      "program", "--Lx", "-1.0", "--Ly", "1.0", "--Lz", "1.0",   "--Nx", "10",
      "--Ny",    "10",   "--Nz", "10",   "--T", "1.0",  "--tau", "0.1"};
  char **argv = makeArgv(args);

  EXPECT_THROW(
      { utils::ArgParser::parse(args.size(), argv); }, utils::ParserException);

  freeArgv(argv, args.size());
}

TEST(ArgParserTest, NegativeNx) {
  std::vector<std::string> args = {
      "program", "--Lx", "1.0",  "--Ly", "1.0", "--Lz", "1.0",   "--Nx", "-10",
      "--Ny",    "10",   "--Nz", "10",   "--T", "1.0",  "--tau", "0.1"};
  char **argv = makeArgv(args);

  EXPECT_THROW(
      { utils::ArgParser::parse(args.size(), argv); }, utils::ParserException);

  freeArgv(argv, args.size());
}

TEST(ArgParserTest, ZeroT) {
  std::vector<std::string> args = {
      "program", "--Lx", "1.0",  "--Ly", "1.0", "--Lz", "1.0",   "--Nx", "10",
      "--Ny",    "10",   "--Nz", "10",   "--T", "0.0",  "--tau", "0.1"};
  char **argv = makeArgv(args);

  EXPECT_THROW(
      { utils::ArgParser::parse(args.size(), argv); }, utils::ParserException);

  freeArgv(argv, args.size());
}

TEST(ArgParserTest, BothTauAndK) {
  std::vector<std::string> args = {"program", "--Lx", "1.0",  "--Ly", "1.0",
                                   "--Lz",    "1.0",  "--Nx", "10",   "--Ny",
                                   "10",      "--Nz", "10",   "--T",  "1.0",
                                   "--tau",   "0.1",  "--K",  "10"};
  char **argv = makeArgv(args);

  EXPECT_THROW(
      { utils::ArgParser::parse(args.size(), argv); }, utils::ParserException);

  freeArgv(argv, args.size());
}

TEST(ArgParserTest, NeitherTauNorK) {
  std::vector<std::string> args = {"program", "--Lx", "1.0",  "--Ly", "1.0",
                                   "--Lz",    "1.0",  "--Nx", "10",   "--Ny",
                                   "10",      "--Nz", "10",   "--T",  "1.0"};
  char **argv = makeArgv(args);

  EXPECT_THROW(
      { utils::ArgParser::parse(args.size(), argv); }, utils::ParserException);

  freeArgv(argv, args.size());
}

TEST(ArgParserTest, UnknownArgument) {
  std::vector<std::string> args = {
      "program", "--Lx",  "1.0",  "--Ly",      "1.0",  "--Lz", "1.0",
      "--Nx",    "10",    "--Ny", "10",        "--Nz", "10",   "--T",
      "1.0",     "--tau", "0.1",  "--unknown", "value"};
  char **argv = makeArgv(args);

  EXPECT_THROW(
      { utils::ArgParser::parse(args.size(), argv); }, utils::ParserException);

  freeArgv(argv, args.size());
}

TEST(ArgParserTest, InvalidDoubleValue) {
  std::vector<std::string> args = {
      "program", "--Lx", "abc",  "--Ly", "1.0", "--Lz", "1.0",   "--Nx", "10",
      "--Ny",    "10",   "--Nz", "10",   "--T", "1.0",  "--tau", "0.1"};
  char **argv = makeArgv(args);

  EXPECT_THROW(
      { utils::ArgParser::parse(args.size(), argv); }, utils::ParserException);

  freeArgv(argv, args.size());
}

TEST(ArgParserTest, InvalidIntValue) {
  std::vector<std::string> args = {
      "program", "--Lx", "1.0",  "--Ly", "1.0", "--Lz", "1.0",   "--Nx", "xyz",
      "--Ny",    "10",   "--Nz", "10",   "--T", "1.0",  "--tau", "0.1"};
  char **argv = makeArgv(args);

  EXPECT_THROW(
      { utils::ArgParser::parse(args.size(), argv); }, utils::ParserException);

  freeArgv(argv, args.size());
}
