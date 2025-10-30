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
      "program", "--Lx", "1.0",  "--Ly", "2.0", "--Lz", "3.0",   "--hx", "0.1",
      "--hy",    "0.2",  "--hz", "0.3",  "--T", "10.0", "--tau", "0.1"};
  char **argv = makeArgv(args);

  utils::CmdArgs params = utils::ArgParser::parse(args.size(), argv);

  EXPECT_DOUBLE_EQ(params.Lx, 1.0);
  EXPECT_DOUBLE_EQ(params.Ly, 2.0);
  EXPECT_DOUBLE_EQ(params.Lz, 3.0);
  EXPECT_DOUBLE_EQ(params.hx, 0.1);
  EXPECT_DOUBLE_EQ(params.hy, 0.2);
  EXPECT_DOUBLE_EQ(params.hz, 0.3);
  EXPECT_FALSE(params.Nx.has_value());
  EXPECT_FALSE(params.Ny.has_value());
  EXPECT_FALSE(params.Nz.has_value());
  EXPECT_DOUBLE_EQ(params.T, 10.0);
  EXPECT_DOUBLE_EQ(params.tau, 0.1);
  EXPECT_EQ(params.K, 100);

  freeArgv(argv, args.size());
}

TEST(ArgParserTest, ValidInputWithK) {
  std::vector<std::string> args = {
      "program", "--Lx", "1.0",  "--Ly", "1.0", "--Lz", "1.0", "--hx", "0.01",
      "--hy",    "0.01", "--hz", "0.01", "--T", "5.0",  "--K", "50"};
  char **argv = makeArgv(args);

  utils::CmdArgs params = utils::ArgParser::parse(args.size(), argv);

  EXPECT_EQ(params.K, 50);
  EXPECT_DOUBLE_EQ(params.tau, 0.1);

  freeArgv(argv, args.size());
}

TEST(ArgParserTest, ComputeKWithNonDivisibleTau) {
  std::vector<std::string> args = {
      "program", "--Lx", "1.0",  "--Ly", "1.0", "--Lz", "1.0",   "--hx", "0.01",
      "--hy",    "0.01", "--hz", "0.01", "--T", "10.0", "--tau", "0.3"};
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
      "program", "--Lx", "-1.0", "--Ly", "1.0", "--Lz", "1.0",   "--hx", "0.01",
      "--hy",    "0.01", "--hz", "0.01", "--T", "1.0",  "--tau", "0.1"};
  char **argv = makeArgv(args);

  EXPECT_THROW(
      { utils::ArgParser::parse(args.size(), argv); }, utils::ParserException);

  freeArgv(argv, args.size());
}

TEST(ArgParserTest, NegativeHx) {
  std::vector<std::string> args = {"program", "--Lx", "1.0",  "--Ly",  "1.0",
                                   "--Lz",    "1.0",  "--hx", "-0.01", "--hy",
                                   "0.01",    "--hz", "0.01", "--T",   "1.0",
                                   "--tau",   "0.1"};
  char **argv = makeArgv(args);

  EXPECT_THROW(
      { utils::ArgParser::parse(args.size(), argv); }, utils::ParserException);

  freeArgv(argv, args.size());
}

TEST(ArgParserTest, ZeroT) {
  std::vector<std::string> args = {
      "program", "--Lx", "1.0",  "--Ly", "1.0", "--Lz", "1.0",   "--hx", "0.01",
      "--hy",    "0.01", "--hz", "0.01", "--T", "0.0",  "--tau", "0.1"};
  char **argv = makeArgv(args);

  EXPECT_THROW(
      { utils::ArgParser::parse(args.size(), argv); }, utils::ParserException);

  freeArgv(argv, args.size());
}

TEST(ArgParserTest, BothTauAndK) {
  std::vector<std::string> args = {"program", "--Lx", "1.0",  "--Ly", "1.0",
                                   "--Lz",    "1.0",  "--hx", "0.01", "--hy",
                                   "0.01",    "--hz", "0.01", "--T",  "1.0",
                                   "--tau",   "0.1",  "--K",  "10"};
  char **argv = makeArgv(args);

  EXPECT_THROW(
      { utils::ArgParser::parse(args.size(), argv); }, utils::ParserException);

  freeArgv(argv, args.size());
}

TEST(ArgParserTest, NeitherTauNorK) {
  std::vector<std::string> args = {"program", "--Lx", "1.0",  "--Ly", "1.0",
                                   "--Lz",    "1.0",  "--hx", "0.01", "--hy",
                                   "0.01",    "--hz", "0.01", "--T",  "1.0"};
  char **argv = makeArgv(args);

  EXPECT_THROW(
      { utils::ArgParser::parse(args.size(), argv); }, utils::ParserException);

  freeArgv(argv, args.size());
}

TEST(ArgParserTest, UnknownArgument) {
  std::vector<std::string> args = {
      "program", "--Lx",  "1.0",  "--Ly",      "1.0",  "--Lz", "1.0",
      "--hx",    "0.01",  "--hy", "0.01",      "--hz", "0.01", "--T",
      "1.0",     "--tau", "0.1",  "--unknown", "value"};
  char **argv = makeArgv(args);

  EXPECT_THROW(
      { utils::ArgParser::parse(args.size(), argv); }, utils::ParserException);

  freeArgv(argv, args.size());
}

TEST(ArgParserTest, InvalidDoubleValue) {
  std::vector<std::string> args = {
      "program", "--Lx", "abc",  "--Ly", "1.0", "--Lz", "1.0",   "--hx", "0.01",
      "--hy",    "0.01", "--hz", "0.01", "--T", "1.0",  "--tau", "0.1"};
  char **argv = makeArgv(args);

  EXPECT_THROW(
      { utils::ArgParser::parse(args.size(), argv); }, utils::ParserException);

  freeArgv(argv, args.size());
}

TEST(ArgParserTest, InvalidDoubleValueForHx) {
  std::vector<std::string> args = {
      "program", "--Lx", "1.0",  "--Ly", "1.0", "--Lz", "1.0",   "--hx", "xyz",
      "--hy",    "0.01", "--hz", "0.01", "--T", "1.0",  "--tau", "0.1"};
  char **argv = makeArgv(args);

  EXPECT_THROW(
      { utils::ArgParser::parse(args.size(), argv); }, utils::ParserException);

  freeArgv(argv, args.size());
}

TEST(ArgParserTest, MPIModeValidInput) {
  std::vector<std::string> args = {
      "program", "--Lx", "1.0",  "--Ly", "1.0",  "--Lz",  "1.0", "--hx",
      "0.01",    "--hy", "0.01", "--hz", "0.01", "--Nx",  "2",   "--Ny",
      "2",       "--Nz", "2",    "--T",  "10.0", "--tau", "0.1"};
  char **argv = makeArgv(args);

  utils::CmdArgs params = utils::ArgParser::parse(args.size(), argv, true);

  EXPECT_DOUBLE_EQ(params.Lx, 1.0);
  EXPECT_DOUBLE_EQ(params.hx, 0.01);
  EXPECT_TRUE(params.Nx.has_value());
  EXPECT_EQ(params.Nx.value(), 2);
  EXPECT_EQ(params.Ny.value(), 2);
  EXPECT_EQ(params.Nz.value(), 2);

  freeArgv(argv, args.size());
}

TEST(ArgParserTest, MPIModeNxNotAllowedInNonMPIMode) {
  std::vector<std::string> args = {"program", "--Lx", "1.0",   "--Ly", "1.0",
                                   "--Lz",    "1.0",  "--hx",  "0.01", "--hy",
                                   "0.01",    "--hz", "0.01",  "--Nx", "2",
                                   "--T",     "10.0", "--tau", "0.1"};
  char **argv = makeArgv(args);

  EXPECT_THROW(
      { utils::ArgParser::parse(args.size(), argv, false); },
      utils::ParserException);

  freeArgv(argv, args.size());
}
